! ******************************************************************************
!
!  README:
!  This module contains extra routines for running MESA star simulations that
!  involve common envelope evolution with a neutron star companion. It includes
!  functions for calculating accretion rates, drag forces, orbital evolution, 
!  and gravitational wave strain.
!
!  At every timestep, the current values (for that timestep) of the relevant
!  variables are already known. The module calculates the next values of these
!  variables based on the chosen prescription, and saves them into the star's
!  xtra array for use in the next timestep. The variables include:
!
!    - orbital separation (a)
!    - neutron star mass (M_ns)
!    - accreted mass (M_acc)
!    - orbital angular velocity (omega)
!    - quadrupole moment (Q)
!    - maximum quadrupole moment (Qmax)
!    - tidal bulge quadrupole moment (Qtb)
!    - eccentricity (e)
!    - x and y positions and velocities of both the primary core and the
!      neutron star (xcore, ycore, vxcore, vycore, xcomp, ycomp, vxcomp, vycomp)
!      (for elliptical orbits)
!    - moment of inertia of the envelope (mom_inert)
!    - strain (strain)
!  
!  The prescription variable determines which orbital evolution prescription
!  to use (Holgado et al. 2018 or Bronner et al. 2024).
!
!  Variable of the type s%xtra(...) refer to the star's xtra array where these
!  variables are stored. MESA remembers these values across timesteps.
!
!  NOTE: All units are in cgs unless otherwise specified.
!  The first function that is called is extras_controls, which sets up the
!  necessary variables and arrays. Then, the evolve_companion_and_inject_energy
!  subroutine is called at every timestep to perform the calculations and
!  update the variables.
!
! ******************************************************************************

module run_star_extras

! ==============================================================================

use star_lib
use star_def
use const_def
use math_lib
   
implicit none

! integration prescriptions
integer, parameter :: ENERGY_PRESCRIPTION = 1, &
                      FORCE_PRESCRIPTION  = 2
! drag prescriptions
integer, parameter :: MR15_DRAG = 1, KIM07_DRAG = 2, KIM10_DRAG = 3

! number of orbit steps per MESA step
! later make this an inlist parameter or choose based on a/|adot|
integer, parameter :: Nsteps = 10

! xtra variables - values describing NS and its orbit
integer, parameter :: a         = 1,  M_ns      = 2,  &
                      M_acc     = 3,  omega_env = 4,  &
                      omega     = 5,  Q         = 6,  &
                      Qmax      = 7,  Qtb       = 8,  &
                      e         = 9,  aa        = 10, &
                      bb        = 11, mom_inert = 12, &
                      strain    = 13, xcore     = 14, &
                      ycore     = 15, xcomp     = 16, &
                      ycomp     = 17, vxcore    = 18, &
                      vycore    = 19, vxcomp    = 20, &
                      vycomp    = 21

! input parameters from inlists
real(dp) :: M_ns_initial, omega_initial, e_initial, D, R0, r1, r2, &
            op_const, eta, efactor, M_crust, n_poly, beta_sec, menc, &
            omega_env_factor

integer  :: prescription  

! other variables
real(dp) :: temp, fd, decay_coeff, Req, Rbar, beta, ebind, eorb_change, &
            mdot, v_rel, Ra_br, mdot_br, fd_br, edot_br, mdot_hl_br, &
            mdot_mr15_ratio_br, mdot_mr15_br, fd_mr15_ratio_br, &
            fd_mr15_br, fd_hl_br, v, edot, csound
integer :: azone

contains

! ==============================================================================

! Evolve the companion's orbit through one MESA timestep and inject the
! appropriate energy into the MESA model.

subroutine evolve_companion_and_inject_energy(id, ierr)

integer, intent(in)       :: id
integer, intent(out)      :: ierr

type (star_info), pointer :: s
integer                   :: i
real(dp)                  :: tlocal, dtlocal

! Get pointer to the MESA model

ierr = 0
call star_ptr(id, s, ierr)
if (ierr /= 0) return

! Initialising - These values will be taken from the inlist_project file
   
s%x_ctrl(14)     = s%model_number
M_ns_initial     = s%x_ctrl(1)*Msun   !mass of the NS (gm)
op_const         = s%x_ctrl(2)        !opacity constant (cm^2/g)
eta              = s%x_ctrl(3)        !efficiency factor
efactor          = s%x_ctrl(4)        !multiplication factor for injected energy
M_crust          = s%x_ctrl(5)*Msun   !mass of the crust (gm)
omega_initial    = 2*pi*s%x_ctrl(6)   !initial spin frequency (Hz)
R0               = s%x_ctrl(7)        !equatorial radius (cm)
e_initial        = s%x_ctrl(8)        !initial ellipticity
n_poly           = s%x_ctrl(9)        !polytropic index
beta_sec         = s%x_ctrl(10)       !beta secular
D                = s%x_ctrl(11)*1d3*pc!distance to the source (cm)
omega_env_factor = s%x_ctrl(12)       !omega_env_factor*initial_orbital_period =
                                      !envelope omega
prescription     = int(s%x_ctrl(13))  !general prescription for the CE evolution

! Initial conditions for companion orbit
! NOTES:
! Presumably this test would not work if we began with a model number /= 1?
! Hard-coded constants here should be given names or made into inlist
! parameters

if (s%model_number == 1) then

  s%xtra(a)         = 290*Rsun !9.8d-1*s%R(1)          !semimajor axis (cm)
  s%xtra(M_ns)      = M_ns_initial                     !NS mass (g)
  s%xtra(M_acc)     = 0                                !gm ! accreted mass?
  s%xtra(omega)     = omega_initial                    !NS spin freq (Hz)
  s%xtra(e)         = e_initial                        !NS ellipticity
  s%xtra(Qmax)      = 4*1d39*(s%xtra(M_acc)/M_crust)   !NS quadrupole moment
                                                       !(g cm^2)
  s%xtra(aa)        = R0*(1+s%xtra(e)/2)               !NS major axis (cm)
  s%xtra(bb)        = R0*(1-s%xtra(e)/2)               !NS minor axis (cm)
  s%xtra(mom_inert) = s%xtra(M_ns)*(s%xtra(aa)**2 + s%xtra(bb)**2)/5
                                                       !NS moment of inertia
                                                       !(g cm^2)
 
  call interpolate(s%xtra(a), menc, s%R, s%m)
  s%xtra(omega_env) = omega_env_factor * &
                      SQRT(standard_cgrav*(s%xtra(M_ns) + menc)/ s%xtra(a)**3)

  if (prescription == FORCE_PRESCRIPTION) then
    s%xtra(xcore)  = 0
    s%xtra(ycore)  = 0
    s%xtra(xcomp)  = s%xtra(a)
    s%xtra(ycomp)  = 0
    s%xtra(vxcore) = 0
    s%xtra(vycore) = 0
    s%xtra(vxcomp) = 0
    s%xtra(vycomp) = SQRT(standard_cgrav*(s%xtra(M_ns) + menc)/s%xtra(a))
  end if

  ! NOTES: what's this?
  s%x_ctrl(15) = 2

end if

! Initialize time integration

tlocal  = s%time 
dtlocal = s%dt_next / Nsteps
s%extra_heat(:) = 0.

! Integrate

do i = 1, Nsteps

  select case (prescription)
    case (ENERGY_PRESCRIPTION)
      stop "Energy prescription not yet supported"
      !call advance_energy_prescription(s, tlocal, dtlocal)
    case (FORCE_PRESCRIPTION)
      call advance_force_prescription(s, tlocal, dtlocal)
  end select

  tlocal = tlocal + dtlocal

enddo

! Print diagnostic info

call print_info

return
end subroutine evolve_companion_and_inject_energy

!NOTES: Not sure what the significance of these controls is
!   !updating the total accreted mass and NS mass
!   if (s%x_ctrl(14) /= s%x_ctrl(15)) then
!      call forward_euler(s%xtra(M_acc_curr), mdot, s%dt_next, M_acc_next)     !gm
!      call forward_euler(s%xtra(M_ns_curr), mdot, s%dt_next, M_ns_next)       !gm
!   else 
!      call forward_euler(s%xtra(M_acc_curr), mdot, 0d0, M_acc_next)     !gm
!      call forward_euler(s%xtra(M_ns_curr), mdot, 0d0, M_ns_next)       !gm
!   end if

! ------------------------------------------------------------------------------

! Integrate the companion orbit and other properties through one timestep using
! the energy prescription.

subroutine advance_energy_prescription(s, t, dt)

type (star_info), pointer, intent(inout) :: s
real(dp), intent(in)                     :: t, dt

!TODO

return
end subroutine advance_energy_prescription

! ------------------------------------------------------------------------------

! Integrate the companion orbit and other properties through one timestep using
! the force prescription.

subroutine advance_force_prescription(s, t, dt)

type (star_info), pointer, intent(inout) :: s
real(dp), intent(in)                     :: t, dt

real(dp) :: x(12), dxdt(12), e_inj, r

! Set up solution vector.

x = (/ s%xtra(xcore), s%xtra(ycore), s%xtra(xcomp), s%xtra(ycomp), &
       s%xtra(vxcore), s%xtra(vycore), s%xtra(vxcomp), s%xtra(vycomp), &
       s%xtra(omega),  s%xtra(M_ns),   s%xtra(M_acc), 0. /)

! Forward Euler step. We can do other integration schemes if we do time
! extrapolation of the MESA model in force_prescrip_rhs.

call force_prescrip_rhs(t, x, dxdt, s)

x = x + dt*dxdt

! Unpack solution vector.

s%xtra(xcore)  = x(1)
s%xtra(ycore)  = x(2)
s%xtra(xcomp)  = x(3)
s%xtra(ycomp)  = x(4)
s%xtra(vxcore) = x(5)
s%xtra(vycore) = x(6)
s%xtra(vxcomp) = x(7)
s%xtra(vycomp) = x(8)
s%xtra(omega)  = x(9)
s%xtra(M_ns)   = x(10)
s%xtra(M_acc)  = x(11)
e_inj          = x(12)

! Add energy to MESA grid.

r = sqrt( (s%xtra(xcomp) - s%xtra(xcore))**2 + &
          (s%xtra(ycomp) - s%xtra(ycore))**2 )
call get_relative_velocity(s, r, s%xtra(xcore), s%xtra(ycore), &
                           s%xtra(xcomp), s%xtra(ycomp), &
                           s%xtra(vxcore), s%xtra(vycore), &
                           s%xtra(vxcomp), s%xtra(vycomp), v_rel)
call add_energy_to_mesa(s, r, s%xtra(M_ns), v_rel, -e_inj)

return
end subroutine advance_force_prescription

! ------------------------------------------------------------------------------

! Return the right-hand-side of the differential equations describing the
! companion's orbit and spin for the energy prescription.

subroutine energy_prescrip_rhs(t, x, dxdt, s)

real(dp), intent(in)                     :: t, x(:)
real(dp), intent(out)                    :: dxdt(:)
type (star_info), pointer, intent(inout) :: s

!TODO

return
end subroutine energy_prescrip_rhs

! ------------------------------------------------------------------------------

! Compute the relative velocity of the NS and the donor envelope.

subroutine get_relative_velocity(s, r, xcore, ycore, xcomp, ycomp, &
                                 vxcore, vycore, vxcomp, vycomp, v_rel)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r
real(dp), intent(out)                 :: v_rel

real(dp)                              :: u

!NOTES: significance of "100"? replace with a named constant
!Should this be 100*Rsun?
if (r > 100) then
  v_rel = sqrt( (vxcomp - vxcore + (s%xtra(omega_env)*(ycomp - ycore)))**2 + &
                (vycomp - vycore - (s%xtra(omega_env)*(xcomp - xcore)))**2 )
else 
  call interpolate(r, u, s%R, s%u)
  v_rel = sqrt( (vxcomp - vxcore + (s%xtra(omega_env)*(ycomp - ycore)) - &
                 u*(xcomp - xcore)/r)**2 + &
                (vycomp - vycore - (s%xtra(omega_env)*(xcomp - xcore)) - &
                 u*(ycomp - ycore)/r)**2 )
end if

return
end subroutine get_relative_velocity

! ------------------------------------------------------------------------------

! Return the right-hand-side of the differential equations describing the
! companion's orbit and spin for the force prescription.

subroutine force_prescrip_rhs(t, x, dxdt, s)

real(dp), intent(in)                     :: t, x(:)
real(dp), intent(out)                    :: dxdt(:)
type (star_info), pointer, intent(inout) :: s

real(dp) :: xcore, ycore, xcomp, ycomp, vxcore, vycore, vxcomp, vycomp, &
            omega, mns, macc, r, v_rel, rho, menc, mdot, fd, edot, &
            Eorb, dErob, Ra, mdot_hl, fd_hl, edot_hl, mdot_edd, mdot_hyper, &
            cs, mach, e_inj

xcore   = x(1)
ycore   = x(2)
xcomp   = x(3)
ycomp   = x(4)
vxcore  = x(5)
vycore  = x(6)
vxcomp  = x(7)
vycomp  = x(8)
omega   = x(9)
mns     = x(10)
macc    = x(11)
e_inj   = x(12)

r = sqrt( (xcore - xcomp)**2 + (ycore - ycomp)**2 )

call get_relative_velocity(s, r, xcore, ycore, xcomp, ycomp, &
                           vxcore, vycore, vxcomp, vycomp, v_rel)

! Get accretion rate
call interpolate(r, rho, s%R, s%rho)
call interpolate(r, menc, s%R, s%m)
call hl_and_energy(s, mns, r, v_rel, rho, menc, mdot_hl, fd_hl, edot_hl, &
                   Eorb, dEorb, Ra)
call mr15(s, r, Ra, mdot_hl, fd_hl, mdot, fd)
call edd_and_hyper(mns, mdot_edd, mdot_hyper)

if (eta*mdot >= mdot_hyper .or. eta*mdot <= mdot_edd) then
  mdot = eta*mdot
else
  mdot = mdot_edd
end if

! Get drag and power
call interpolate(r, cs, s%R, s%csound)
mach = v_rel / cs
call kim2010(s, mach, r, Ra, cs, rho, v_rel, mns, fd)
edot = fd * v_rel

dxdt(1)  = vxcore
dxdt(2)  = vycore
dxdt(3)  = vxcomp
dxdt(4)  = vycomp
dxdt(5)  = standard_cgrav*mns*(xcomp - xcore)/r**3
dxdt(6)  = standard_cgrav*mns*(ycomp - ycore)/r**3
!NOTES: significance of 100*Rsun? see note above
if (r > 100*Rsun) then
  dxdt(7)  = standard_cgrav*menc*(xcomp - xcore)/r**3 + &
            fd*(vxcomp - vxcore + s%xtra(omega_env)*(ycomp - ycore))/(mns*v_rel)
  dxdt(8)  = standard_cgrav*menc*(ycomp - ycore)/r**3 + &
            fd*(vycomp - vycore - s%xtra(omega_env)*(xcomp - xcore))/(mns*v_rel)
else
  dxdt(7)  = standard_cgrav*menc*(xcomp - xcore)/r**3 + &
             fd*(vxcomp - vxcore + s%xtra(omega_env)*(ycomp - ycore) - &
             (s%u(azone)*(xcomp - xcore))/r)/(mns*v_rel)
  dxdt(8)  = standard_cgrav*menc*(ycomp - ycore)/r**3 + &
             fd*(vycomp - vycore - s%xtra(omega_env)*(xcomp - xcore) - &
             (s%u(azone)*(ycomp - ycore)/r)/(mns*v_rel)
endif
!TODO
dxdt(9)  = ...
dxdt(10) = mdot
dxdt(11) = mdot
dxdt(12) = edot

return
end subroutine force_prescrip_rhs

! ------------------------------------------------------------------------------

! Given the NS position and the amount of energy to inject, add thermal energy
! to the MESA model according to the Bronner prescription.

subroutine add_energy_to_mesa(s, r, mns, v_rel, e_inj)

type (star_info), pointer, intent(inout) :: s
real(dp), intent(in)                     :: r, mns, v_rel, e_inj, junk

real(dp) :: rho, menc, Ra, kernel_sum
integer  :: n, j1, j2, j
real(dp), allocatable :: temp_heat(:)

! Get the accretion radius and the indices of zones we will modify

n = size(s%R)
call interpolate(r, rho, s%R, s%rho)
call interpolate(r, menc, s%R, s%m)
call hl_and_energy(s, mns, r, v_rel, rho, menc, junk, junk, junk, junk, junk, &
                   Ra)
call hunt(s%R, n, r-Ra, j1)
call hunt(s%R, n, r+Ra, j2)
j1 = max(j1, 1)
j2 = max(j2, 1)

! Set up the temporary array to keep track of the energy added

allocate(temp_heat(n))
temp_heat(:) = 0.
kernel_sum = 0.

! Compute the heating kernel

do j = j1, j2
  ! NOTES: not sure I have translated this correctly.
  temp_heat(j) = exp(-((r-s%R(j))/Ra)**2)
  kernel_sum = kernel_sum + temp_heat(j)*s%dm(j)
enddo

do j = j1, j2
  s%extra_heat(j)%val = s%extra_heat(j)%val + &
                        efactor * temp_heat(j) * e_inj / kernel_sum
enddo

deallocate(temp_heat)

return
end subroutine add_energy_to_mesa

! ------------------------------------------------------------------------------

! Kim 2007 drag force calculation

subroutine kim2007(s, mach_in, r_in, rho_in, v_in, M_in, fd_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: M_in, mach_in, a_in
real(dp, intent(out)                  :: fd_out

if (mach_in < 1.0) then
  fd_out = 0.7706*LOG((1+mach_in)/(1.0004 - 0.9185*mach_in)) - 1.473*mach_in
else if (mach_in >= 1.0 .and. mach_in < 4.4) then
  fd_out = LOG(330*(r_in*(mach_in-0.71)**5.72)/(R0*mach_in**9.58))
else
  fd_out = LOG(r_in / (R0*(0.11*mach_in + 1.65)))
end if
fd_out = fd_out*(4*pi*rho_in*((standard_cgrav*M_in)**2)/(v_in**2))

return
end subroutine kim2007

! ------------------------------------------------------------------------------

! Kim & Kim 2010 drag force calculation

subroutine kim2010(s, mach_in, a_in, Ra_in, cs_in, rho_in, v_in, M_in, fd_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: M_in, mach_in, a_in, Ra_in, v_in, &
                                         cs_in, rho_in
real(dp), intent(out)                 :: fd_out

integer  :: k
real(dp) :: i_var, beta, eta_b, cd

beta = standard_cgrav*M_in/(a_in*cs_in**2)   !dimensionless
eta_b = beta/(mach_in**2 - 1)                !dimensionless
cd = 0.002

if (mach_in < 1.01) then
  i_var = 0.7706*LOG((1+mach_in)/(1.0004 - 0.9185*mach_in)) - 1.473*mach_in
else if (mach_in >= 1.01 .and. mach_in < 4.4) then
  i_var = LOG(330*(a_in*(mach_in-0.71)**5.72)/(Ra_in*mach_in**9.58))
else
  i_var = LOG(a_in / (Ra_in*(0.11*mach_in + 1.65)))
end if

if (eta_b > 0.1 .and. mach_in > 1.01) then
  fd_out = -cd*0.7*4*pi*rho_in* &
           (1 + 0.46*(beta**1.1)/(mach_in**2 - 1)**0.11)* &
           (standard_cgrav*(M_in**2))/((v_in**2)*(eta_b**0.5))
else
  fd_out = -cd*4*pi*rho_in*(standard_cgrav*(M_in**2))*i_var/(v_in**2)
end if

return
end subroutine kim2010

! ------------------------------------------------------------------------------

! MR15 accretion and drag force calculation

subroutine mr15(s, r_in, Ra_in, mdot_hl, fd_hl, mdotut_out, fd_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r_in, Ra_in, mdot_hl, fd_hl
real(dp), intent(out)                 :: fd_out

real(dp), parameter :: f1 = 1.91791946, f2 = -1.52814698, f3 = 0.75992092, &
                       mu1 = -2.14034214, mu2 = 1.94694764, mu3 = 1.19007536, &
                       mu4 = 1.05762477
real(dp)            :: rsc, eps_rho, fd_mr15_ratio, mdot_mr15_ratio

call interpolate(r_in, rsc, s%R, s%scale_height)

eps_rho         = Ra_in / rsc ! dimensionless
fd_mr15_ratio   = f1 + f2*eps_rho + f3*(eps_rho**2) ! dimensionless
mdot_mr15_ratio = 10**(mu1 + mu2/(1 + mu3*eps_rho + mu4*(eps_rho**2))) !dimless
mdot_out        = mdot_mr15_ratio*mdot_hl !g/s
fd_out          = fd_mr15_ratio*fd_hl !dyne (g cm s^-2)

return
end subroutine mr15

! ------------------------------------------------------------------------------

! Hoyle-Lyttleton accretion and energy calculation  

subroutine hl_and_energy(s, M, r_in, v_in, rho_in, menc_in, mdot_hl, fd_hl, &
                         edot_hl, Eorb, dEorb, Ra)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: M, r_in, v_in, rho_in, menc_in
real(dp), intent(out)                 :: mdot_hl, fd_hl, edot_hl, Eorb, dEorb, &
                                         Ra

! mass is in g, radius is in cm, time is in s

Ra      = 2*standard_cgrav*M/(v_in**2)       !cm
mdot_hl = pi*(Ra**2)*rho_in*v_in             !g/s
fd_hl   = mdot_hl*v_in                       !dyne (g cm s^-2)
edot_hl = fd_hl*v_in                         !erg/s
      
Eorb    = standard_cgrav*M*menc_in/(2*r_in)  !erg
dEorb   = (standard_cgrav*M/(2*r_in))*(menc_in/r_in - 4*pi*(r_in**2)*rho_in)
                                             !erg/cm

return
end subroutine hl_and_energy

! ------------------------------------------------------------------------------

! Eddington and hypercritical accretion rate calculation

subroutine edd_and_hyper(M, mdot_edd, mdot_hyper)

real(dp), intent(in)  :: M
real(dp), intent(out) :: mdot_edd, mdot_hyper

! mass is in g, radius is in cm, time is in s

mdot_edd   = 3.5*(1d-8)*(M/(1.33*Msun))*(0.34/op_const)*Msun/secyer      !gm/s
mdot_hyper = 8.9*(1d-5)*((op_const/0.34)**(-0.73))*Msun/secyer           !gm/s

return
end subroutine edd_and_hyper

! ------------------------------------------------------------------------------

! Find the equatorial radius and beta secular for the NS

!TODO
subroutine Req_and_beta(e1, Req1, Rbar1, beta1)

real(dp), intent(in)  :: e1
real(dp), intent(out) :: Req1, Rbar1, beta1
   
beta1 = 3*(1 - ((e1*sqrt(1-e1**2))/(asin(e1))))/(2*e1**2) - 1  !dimensionless
Rbar1 = R0*((asin(e1) * ((1-e1**2)**(1/6)) * (1-beta1)) / e1)**(-1 * n_poly / &
        (3-n_poly))                                            !cm
Req1 = Rbar1/((1-e1**2)**(1/6))                                !cm

return
end subroutine Req_and_beta

! ------------------------------------------------------------------------------

! Evaluate the omega function and solve for the inverse

!TODO
subroutine omega_function(e1, M, omega1)

real(dp), intent(in)  :: e1, M
real(dp), intent(out) :: omega1
real(dp)              :: rho_bar, qn

rho_bar = 3*M/(4*pi*(R0**3))      !gm/cm^3
qn = (1-n_poly/5)                 !dimensionless
   
omega1 = sqrt(2*pi*standard_cgrav*rho_bar* &
              ((sqrt(1-e1**2)*(3-2*e1**2)*asin(e1)/(e1**3)) - &
               3*(1-e1**2)/(e1**2) )/qn) !Hz

return
end subroutine omega_function

! ------------------------------------------------------------------------------

! Solve for the inverse of the omega function

!TODO
subroutine omega_func_solve_inverse(e_in, omega_val, tol, M, e_out)

real(dp), intent(in)  :: e_in, omega_val, tol, M
real(dp), intent(out) :: e_out

real(dp) :: e1, omega1
integer  :: iterations

iterations = 0

e1 = e_in !dimensionless
call omega_function(e1, M, omega1)

do while (abs(omega1 - omega_val) >= tol .and. iterations < 100)
  e1 = e1*omega_val/omega1
  call omega_function(e1, M, omega1)
  iterations = iterations + 1
end do

e_out = e1 !dimensionless

return
end subroutine omega_func_solve_inverse

! ------------------------------------------------------------------------------

! Evaluate the spin evolution and quadrupole moment evolution

!TODO
subroutine omega_and_q(id, ierr, M, omega_function)

integer, intent(in) :: id
integer, intent(out) :: ierr
type (star_info), pointer :: s

!subroutine variables
real(dp), intent(in) :: M, omega_function
integer :: i

!calling star pointer
ierr = 0
call star_ptr(id, s, ierr)
if (ierr /= 0) return
Qmax_next = 4*1d39*(M_acc_next/M_crust)             !g cm^2
                     
if (s%x_ctrl(14) /= s%x_ctrl(15)) then
  call forward_euler(s%xtra(omega_curr), omega_function, s%dt_next, omega_next)  !Hz
else 
  call forward_euler(s%xtra(omega_curr), omega_function, 0d0, omega_next)          !Hz
end if

call omega_func_solve_inverse(s%xtra(e_curr), omega_next, 1d-9, M, e_next)

aa_next = R0*(1 + e_next/2)                        !cm
bb_next = R0*(1 - e_next/2)                        !cm
mom_inert_next = M*(aa_next**2 + bb_next**2)/5     !gm cm^2

call Req_and_beta(e_next, Req, Rbar, beta)
Qtb_next = sqrt((5*(clight**5)*mdot*sqrt(standard_cgrav*M*Req))/(32*standard_cgrav*(omega_next**5))) !g cm^2
if (Qtb_next > Qmax_next) then
  Q_next = Qmax_next*exp(-1*s%time*mdot/M_crust)   !g cm^2
  decay_coeff = exp(-1*s%time*mdot/M_crust)        !dimensionless
else
  Q_next = Qtb_next !g cm^2
  decay_coeff = exp(-1*s%time*mdot/M_crust) !dimensionless
end if

return
end subroutine omega_and_q

! ------------------------------------------------------------------------------

! Print some useful information to stdout

subroutine print_info()

print *, 'model number                  = ', s% model_number
print *, 'orbital separation            = ', s%xtra(a_curr)/Rsun
print *, 'omega                         = ', s%xtra(omega_curr)
print *, 'Q                             = ', s%xtra(Q_curr)
print *, 'strain                        = ', s%xtra(strain_curr)

select case (prescription)
  case (ENERGY_PRESCRIPTION)
    print *, 'mach                  = ', v/csound
    print *, 'edot                  = ', edot
  case (FORCE_PRESCRIPTION)
    print *, 'mach                  = ', v_rel/csound
    print *, 'edot                  = ', -edot_br
    ! print *, 'fd_br                  = ', fd_br
    print *, 'xcore_curr            = ', s%xtra(xcore_curr)
    print *, 'ycore_curr            = ', s%xtra(ycore_curr)
    print *, 'xcomp_curr            = ', s%xtra(xcomp_curr)
    print *, 'ycomp_curr            = ', s%xtra(ycomp_curr)
    print *, 'vxcore_curr           = ', s%xtra(vxcore_curr)
    print *, 'vycore_curr           = ', s%xtra(vycore_curr)
    print *, 'vxcomp_curr           = ', s%xtra(vxcomp_curr)
    print *, 'vycomp_curr           = ', s%xtra(vycomp_curr)
    print *, 'total injected energy = ', &
             efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
end select

print *, '##############################################################################'
print *, 'ONE TIMESTEP DONE'
print *, '##############################################################################'

return
end subroutine print_info

! ------------------------------------------------------------------------------

!TODO
!DONE: to evaluate the gravitational wave strain
subroutine evaluate_strain(id, ierr)
   !star variables
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s

   !subroutine variables
   real(dp) :: omega_func_curr

   !calling star pointer
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   if (s%model_number == 1) then
      call Req_and_beta(s%xtra(e_curr), Req, Rbar, beta)
      omega_next = s%xtra(omega_curr)                             !Hz
      s%xtra(Qtb_curr) = sqrt((5*(clight**5)*mdot*sqrt(standard_cgrav*s%xtra(M_ns_curr)*Req))/(32*standard_cgrav*((s%xtra(omega_curr))**5))) !g cm^2
      Qtb_next = s%xtra(Qtb_curr)                                 !g cm^2
      Qmax_next = s%xtra(Qmax_curr)                               !g cm^2
      s%xtra(Q_curr) = s%xtra(Qtb_curr)                           !g cm^2
      Q_next = s%xtra(Q_curr)                                     !g cm^2
      e_next = s%xtra(e_curr)                                     !dimensionless
      aa_next = s%xtra(aa_curr)                                   !cm
      bb_next = s%xtra(bb_curr)                                   !cm
      mom_inert_next = s%xtra(mom_inert_curr)                     !gm cm^2

      decay_coeff = 1                                             !dimensionless

      s%xtra(strain_curr) = 2*standard_cgrav*((s%xtra(omega_curr))**2)*(s%xtra(Q_curr))/(D*(clight**4)) !dimensionless
      strain_next = s%xtra(strain_curr)                            !dimensionless

   else 
      if (s%xtra(e_curr) > 0.817) then 
         s%xtra(e_curr) = 0.817                               !dimensionless
      end if

      call Req_and_beta(s%xtra(e_curr), Req, Rbar, beta)
      omega_func_curr = (mdot*sqrt(standard_cgrav*s%xtra(M_ns_curr)*Req) - (32*standard_cgrav*(s%xtra(omega_curr)**5)*(s%xtra(Q_curr)**2))/(5*clight**5))/s%xtra(mom_inert_curr)
      call omega_and_q(id, ierr, M_ns_next, omega_func_curr)
      strain_next = 2*standard_cgrav*(omega_next**2)*Q_next/(D*(clight**4))  !dimensionless

   end if

end subroutine evaluate_strain

! ==============================================================================

! Utility routines

! Find location in a monotonically varying array corresponding to a given
! value, or 0 or N if out of range.

subroutine hunt(xx, n, x, jlo)

  INTEGER  :: jlo,n
  REAL(dp) :: x,xx(n)
  INTEGER  :: inc,jhi,jm
  LOGICAL  :: ascnd
  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
    jlo=0
    jhi=n+1
    goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1   jhi=jlo+inc
    if(jhi.gt.n)then
      jhi=n+1
    else if(x.ge.xx(jhi).eqv.ascnd)then
      jlo=jhi
      inc=inc+inc
      goto 1
    endif
  else
    jhi=jlo
2   jlo=jhi-inc
    if(jlo.lt.1)then
      jlo=0
    else if(x.lt.xx(jlo).eqv.ascnd)then
      jhi=jlo
      inc=inc+inc
      goto 2
    endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.xx(jm).eqv.ascnd)then
    jlo=jm
  else
    jhi=jm
  endif
  goto 3

end subroutine hunt

! ------------------------------------------------------------------------------

! Linear interpolation from an array.

subroutine interpolate(x, y, xtable, ytable)

real(dp), intent(in)  :: x, xtable(:), ytable(:)
real(dp), intent(out) :: y

integer               :: jlo, n
real(dp)              :: s

n = size(xtable)
call hunt(xtable, n, x, jlo)

if ((jlo >= 1) .and. (jlo < n)) then
  s = (ytable(jlo+1)-ytable(jlo)) / (xtable(jlo+1)-xtable(jlo))
  y = ytable(jlo) + s*(x - xtable(jlo))
elseif (jlo < 1) then
  y = ytable(1)
else
  y = ytable(n)
endif

return
end subroutine interpolate

! ==============================================================================

! MESA hooks
   
subroutine extras_controls(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! this is the place to set any procedure pointers you want to change
   ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

   ! the extras functions in this file will not be called
   ! unless you set their function pointers as done below.
   ! otherwise we use a null_ version which does nothing (except warn).

   s% other_energy => evolve_companion_and_inject_energy
   s% extras_startup => extras_startup
   s% extras_start_step => extras_start_step
   s% extras_check_model => extras_check_model
   s% extras_finish_step => extras_finish_step
   s% extras_after_evolve => extras_after_evolve
   s% how_many_extra_history_columns => how_many_extra_history_columns
   s% data_for_extra_history_columns => data_for_extra_history_columns
   s% how_many_extra_profile_columns => how_many_extra_profile_columns
   s% data_for_extra_profile_columns => data_for_extra_profile_columns  

   s% how_many_extra_history_header_items => how_many_extra_history_header_items
   s% data_for_extra_history_header_items => data_for_extra_history_header_items
   s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
   s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

end subroutine extras_controls

! ------------------------------------------------------------------------------

subroutine extras_startup(id, restart, ierr)
   integer, intent(in) :: id
   logical, intent(in) :: restart
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

end subroutine extras_startup

! ------------------------------------------------------------------------------

integer function extras_start_step(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_start_step = 0
end function extras_start_step

! ------------------------------------------------------------------------------

! returns either keep_going, retry, or terminate.
integer function extras_check_model(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_check_model = keep_going         
!TODO
   if (.false. .and. s% star_mass_h1 < 0.35d0) then
      ! stop when star hydrogen mass drops to specified level
      extras_check_model = terminate
      write(*, *) 'have reached desired hydrogen mass'
      return
   end if
   ! if you want to check multiple conditions, it can be useful
   ! to set a different termination code depending on which
   ! condition was triggered.  MESA provides 9 customizeable
   ! termination codes, named t_xtra1 .. t_xtra9.  You can
   ! customize the messages that will be printed upon exit by
   ! setting the corresponding termination_code_str value.
   ! termination_code_str(t_xtra1) = 'my termination condition'

   if (s% xtra(a_curr) .lt. s%R(s%nz)) then
      extras_check_model = terminate
      s% termination_code = t_xtra1
      termination_code_str(t_xtra1) = 'orbital separation less than core size'
      return
   end if

   ! by default, indicate where (in the code) MESA terminated
   if (extras_check_model == terminate) s% termination_code = t_extras_check_model
end function extras_check_model

! ------------------------------------------------------------------------------

integer function how_many_extra_history_columns(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

!TODO
   if (prescription == 1) then
      how_many_extra_history_columns = 44
   else
      how_many_extra_history_columns = 44
   end if
end function how_many_extra_history_columns

! ------------------------------------------------------------------------------

subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len=maxlen_history_column_name) :: names(n)
   real(dp) :: vals(n)
   integer, intent(out) :: ierr
   integer :: k, zone_temp
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! note: do NOT add the extras names to history_columns.list
   ! the history_columns.list is only for the built-in history column options.
   ! it must not include the new column names you are adding here.

!TODO
   names(1) = 'a_curr'
   names(2) = 'Injected_E_per_timestep'
   names(3) = 'mdot_ns'
   names(4) = 'azone'
   names(5) = 'r1'
   names(6) = 'r2'
   names(7) = 'M_ns'
   names(8) = 'M_acc'
   names(9) = '-ebind_curr'
   names(10) = 'delta_Eorb'
   names(11) = 'Qmax'
   names(12) = 'Qtb'
   names(13) = 'Q'
   names(14) = 'omega'
   names(15) = 'Req'
   names(16) = 'Rbar'
   names(17) = 'beta'
   names(18) = 'e'
   names(19) = 'aa'
   names(20) = 'bb'
   names(21) = 'mom_inert'
   names(22) = 'strain'
   names(23) = 'mass_at_a'
   names(24) = 'R_azone'
   names(25) = 'decay_coeff'
   names(26) = 'v_esc'
   names(27) = 'u_k'
   names(28) = 'csound_ns'
   names(29) = 'lnT_ns'
   names(30) = 'rho_ns'
   names(31) = 'lnPgas_ns'
   names(32) = 'gamma1_ns'
   names(33) = 'scale_height_ns'
   names(34) = 'R99'
   names(35) = 'R95'
   names(36) = 'R90'
   names(37) = 'mdot_edd'
   names(38) = 'v_ns'
   names(39) = 'Ra'
   names(40) = 'mdot_hl'
   names(41) = 'mdot_MR15'
   names(42) = 'fd_ns'
   names(43) = 'fd_hl'
   names(44) = 'fd_MR15'

   vals(1) = s%xtra(a_curr)
   vals(2) = efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
   vals(3) = mdot
   vals(4) = azone
   vals(5) = r1
   vals(6) = r2
   vals(7) = s%xtra(M_ns_curr)
   vals(8) = s%xtra(M_acc_curr)
   vals(9) = -ebind
   vals(10) = eorb_change
   vals(11) = s%xtra(Qmax_curr)
   vals(12) = s%xtra(Qtb_curr)
   vals(13) = s%xtra(Q_curr)
   vals(14) = s%xtra(omega_curr)
   vals(15) = Req
   vals(16) = Rbar
   vals(17) = beta
   vals(18) = s%xtra(e_curr)
   vals(19) = s%xtra(aa_curr)
   vals(20) = s%xtra(bb_curr)
   vals(21) = s%xtra(mom_inert_curr)
   vals(22) = s%xtra(strain_curr)
   vals(23) = s%m(azone)
   vals(24) = s%R(azone)
   vals(25) = decay_coeff
   vals(26) = sqrt(2*standard_cgrav*s% m(azone)/s% R(azone))
   vals(27) = s%u(azone)
   vals(28) = s%csound(azone)
   vals(29) = s%lnT(azone)
   vals(30) = s%rho(azone)
   vals(31) = s%lnPgas(azone)
   vals(32) = s%gamma1(azone)
   vals(33) = s%scale_height(azone) 
   call find_zone(s%m(1:s%nz), 0.99*s%m(1), zone_temp)
   vals(34) = s%R(zone_temp)

   call find_zone(s%m(1:s%nz), 0.95*s%m(1), zone_temp)
   vals(35) = s%R(zone_temp)

   call find_zone(s%m(1:s%nz), 0.90*s%m(1), zone_temp)
   vals(36) = s%R(zone_temp)
   vals(37) = mdot_edd(azone)

   if (prescription == 1) then
      vals(38) = v(azone)
      vals(39) = Ra(azone)
      vals(40) = mdot_hl(azone)
      vals(41) = mdot_mr15(azone)
      vals(42) = fd_arr(azone)
      vals(43) = fd_hl(azone)
      vals(44) = fd_mr15(azone)
   else if (prescription == 2) then
      vals(38) = v_rel
      vals(39) = Ra_br
      vals(40) = mdot_hl_br
      vals(41) = mdot_mr15_br
      vals(42) = fd_br
      vals(43) = fd_hl_br
      vals(44) = fd_mr15_br
   end if

end subroutine data_for_extra_history_columns

! ------------------------------------------------------------------------------

integer function how_many_extra_profile_columns(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
!TODO
   if (prescription == 1) then
      how_many_extra_profile_columns = 14
   else if (prescription == 2) then
      how_many_extra_profile_columns = 0
   end if
end function how_many_extra_profile_columns

! ------------------------------------------------------------------------------

subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
   integer, intent(in) :: id, n, nz
   character (len=maxlen_profile_column_name) :: names(n)
   real(dp) :: vals(nz,n)
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   integer :: k
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   
   ! note: do NOT add the extra names to profile_columns.list
   ! the profile_columns.list is only for the built-in profile column options.
   ! it must not include the new column names you are adding here.

   ! here is an example for adding a profile column
   !if (n /= 1) stop 'data_for_extra_profile_columns'
   !names(1) = 'beta'
   !do k = 1, nz
   !   vals(k,1) = s% Pgas(k)/s% P(k)
   !end do
!TODO
   if (prescription == 1) then
      names(1) = 'mdot_hl'
      names(2) = 'mdot_MR15'
      names(3) = 'mdot'
      names(4) = 'fd_hl'
      names(5) = 'fd_MR15'
      names(6) = 'fd'
      names(7) = 'edot_hl'
      names(8) = 'edot'
      names(9) = 'eps_rho'
      names(10) = 'v_ns'
      names(11) = 'Ra'
      names(12) = 'Eorb'
      names(13) = 'dEorb'
      names(14) = 'f'

      do k = 1, s%nz
         vals(k,1) = mdot_hl(k)
         vals(k,2) = mdot_mr15(k)
         vals(k,3) = mdot_arr(k)
         vals(k,4) = fd_hl(k)
         vals(k,5) = fd_mr15(k)
         vals(k,6) = fd_arr(k)
         vals(k,7) = edot_hl(k)
         vals(k,8) = edot(k)
         vals(k,9) = eps_rho(k)
         vals(k,10) = v(k)
         vals(k,11) = Ra(k)
         vals(k,12) = Eorb(k)
         vals(k,13) = dEorb(k)
         vals(k,14) = f(k)
      end do
   end if
end subroutine data_for_extra_profile_columns

! ------------------------------------------------------------------------------

integer function how_many_extra_history_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_history_header_items = 0
end function how_many_extra_history_header_items

! ------------------------------------------------------------------------------

subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len=maxlen_history_column_name) :: names(n)
   real(dp) :: vals(n)
   type(star_info), pointer :: s
   integer, intent(out) :: ierr
   ierr = 0
   call star_ptr(id,s,ierr)
   if(ierr/=0) return

   ! here is an example for adding an extra history header item
   ! also set how_many_extra_history_header_items
   ! names(1) = 'mixing_length_alpha'
   ! vals(1) = s% mixing_length_alpha

end subroutine data_for_extra_history_header_items

! ------------------------------------------------------------------------------

integer function how_many_extra_profile_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_profile_header_items = 0
end function how_many_extra_profile_header_items

! ------------------------------------------------------------------------------

subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
   integer, intent(in) :: id, n
   character (len=maxlen_profile_column_name) :: names(n)
   real(dp) :: vals(n)
   type(star_info), pointer :: s
   integer, intent(out) :: ierr
   ierr = 0
   call star_ptr(id,s,ierr)
   if(ierr/=0) return

   ! here is an example for adding an extra profile header item
   ! also set how_many_extra_profile_header_items
   ! names(1) = 'mixing_length_alpha'
   ! vals(1) = s% mixing_length_alpha

!TODO
   names(1) = 'a_next'
   names(2) = 'mdot_ns'
   names(3) = 'azone'
   names(4) = 'r1'
   names(5) = 'r2'
   names(6) = 'M_acc'
   names(7) = 'Q'
   names(8) = 'omega'
   names(9) = 'v_ns'
   names(10) = 'fd_ns'

   vals(1) = s%xtra(a_curr)
   vals(2) = mdot
   vals(3) = azone
   vals(4) = r1
   vals(5) = r2
   vals(6) = s%xtra(M_acc_curr)
   vals(7) = s%xtra(Q_curr)
   vals(8) = s%xtra(omega_curr)

   if (prescription == 1) then 
      vals(9) = v(azone)
      vals(10) = fd_arr(azone)
   else if (prescription == 2) then
      vals(9) = v_rel
      vals(10) = fd_br
   end if

end subroutine data_for_extra_profile_header_items

! ------------------------------------------------------------------------------


! returns either keep_going or terminate.
! note: cannot request retry; extras_check_model can do that.
integer function extras_finish_step(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   extras_finish_step = keep_going

   ! to save a profile, 
      ! s% need_to_save_profiles_now = .true.
   ! to update the star log,
      ! s% need_to_update_history_now = .true.

   ! see extras_check_model for information about custom termination codes
   ! by default, indicate where (in the code) MESA terminated
   if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step
end function extras_finish_step

! ------------------------------------------------------------------------------

subroutine extras_after_evolve(id, ierr)
   integer, intent(in) :: id
   integer, intent(out) :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
end subroutine extras_after_evolve

! ==============================================================================

end module run_star_extras
