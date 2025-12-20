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
integer, parameter :: ia         = 1,  iM_ns      = 2,  &
                      iM_acc     = 3,  iomega_env = 4,  &
                      iomega     = 5,  iQ         = 6,  &
                      iQmax      = 7,  iQtb       = 8,  &
                      ie         = 9,  iaa        = 10, &
                      ibb        = 11, imom_inert = 12, &
                      istrain    = 13, ixcore     = 14, &
                      iycore     = 15, ixcomp     = 16, &
                      iycomp     = 17, ivxcore    = 18, &
                      ivycore    = 19, ivxcomp    = 20, &
                      ivycomp    = 21

! input parameters from inlists
real(dp) :: M_ns_initial, omega_initial, e_initial, D, R0, r1, r2, &
            op_const, eta, efactor, M_crust, n_poly, beta_sec, menc, &
            omega_env_factor
integer  :: prescription  

! other variables
!real(dp) :: temp, fd, decay_coeff, Req, Rbar, beta, ebind, eorb_change, &
!            mdot, v_rel, Ra_br, mdot_br, fd_br, edot_br, mdot_hl_br, &
!            mdot_mr15_ratio_br, mdot_mr15_br, fd_mr15_ratio_br, &
!            fd_mr15_br, fd_hl_br, v, edot, csound, vx_rel, vy_rel
!integer :: azone

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
! NOTES: Hard-coded constants here should be given names or made into inlist
! parameters

if (s%model_number == 1) then

  s%xtra(ia)         = 290*Rsun !9.8d-1*s%R(1)          !semimajor axis (cm)
  s%xtra(iM_ns)      = M_ns_initial                     !NS mass (g)
  s%xtra(iM_acc)     = 0                                !gm ! accreted mass?
  s%xtra(iomega)     = omega_initial                    !NS spin freq (Hz)
  s%xtra(ie)         = e_initial                        !NS ellipticity
  s%xtra(iQmax)      = 4*1d39*(s%xtra(M_acc)/M_crust)   !NS quadrupole moment
                                                       !(g cm^2)
  s%xtra(iaa)        = R0*(1+s%xtra(e)/2)               !NS major axis (cm)
  s%xtra(ibb)        = R0*(1-s%xtra(e)/2)               !NS minor axis (cm)
  s%xtra(imom_inert) = s%xtra(M_ns)*(s%xtra(aa)**2 + s%xtra(bb)**2)/5
                                                       !NS moment of inertia
                                                       !(g cm^2)
 
  call interpolate(s%xtra(ia), menc, s%R, s%m)
  s%xtra(iomega_env) = omega_env_factor * &
                      SQRT(standard_cgrav*(s%xtra(iM_ns) + menc)/ s%xtra(ia)**3)

  if (prescription == FORCE_PRESCRIPTION) then
    s%xtra(ixcore)  = 0
    s%xtra(iycore)  = 0
    s%xtra(ixcomp)  = s%xtra(ia)
    s%xtra(iycomp)  = 0
    s%xtra(ivxcore) = 0
    s%xtra(ivycore) = 0
    s%xtra(ivxcomp) = 0
    s%xtra(ivycomp) = SQRT(standard_cgrav*(s%xtra(iM_ns) + menc)/s%xtra(ia))
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

!TODO
! Update s%xtra() quantities and other things to be recorded
call update_quantities

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

x = (/ s%xtra(ixcore),  s%xtra(iycore),  s%xtra(ixcomp),  s%xtra(iycomp), &
       s%xtra(ivxcore), s%xtra(ivycore), s%xtra(ivxcomp), s%xtra(ivycomp), &
       s%xtra(iomega),  s%xtra(iM_ns),   s%xtra(iM_acc),  0. /)

! Forward Euler step. We can do other integration schemes if we do time
! extrapolation of the MESA model in force_prescrip_rhs.

call force_prescrip_rhs(t, x, dxdt, s)

x = x + dt*dxdt

! Unpack solution vector.

s%xtra(ixcore)  = x(1)
s%xtra(iycore)  = x(2)
s%xtra(ixcomp)  = x(3)
s%xtra(iycomp)  = x(4)
s%xtra(ivxcore) = x(5)
s%xtra(ivycore) = x(6)
s%xtra(ivxcomp) = x(7)
s%xtra(ivycomp) = x(8)
s%xtra(iomega)  = x(9)
s%xtra(iM_ns)   = x(10)
s%xtra(iM_acc)  = x(11)
e_inj           = x(12)

! Add energy to MESA grid. Note that this is done using updated positions/
! velocities; consider whether this should be done using pre-update values,
! and if higher-order integrator is used, whether we should iterate to
! convergence for this timestep.

r = sqrt( (s%xtra(ixcomp) - s%xtra(ixcore))**2 + &
          (s%xtra(iycomp) - s%xtra(iycore))**2 )
call get_relative_velocity(s, r, s%xtra(ixcore), s%xtra(iycore), &
                           s%xtra(ixcomp), s%xtra(iycomp), &
                           s%xtra(ivxcore), s%xtra(ivycore), &
                           s%xtra(ivxcomp), s%xtra(ivycomp), &
                           vx_rel, vy_rel, v_rel)
call add_energy_to_mesa(s, r, s%xtra(iM_ns), v_rel, e_inj)

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
                                 vxcore, vycore, vxcomp, vycomp, &
                                 vx_rel, vy_rel, v_rel)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r
real(dp), intent(out)                 :: vx_rel, vy_rel, v_rel

real(dp)                              :: u

vx_rel = (vxcomp - vxcore + s%xtra(iomega_env)*(ycomp - ycore))**2
vy_rel = (vycomp - vycore - s%xtra(iomega_env)*(xcomp - xcore))**2

!NOTES: replace 100*Rsun with a named constant
if (r <= 100*Rsun) then
  call interpolate(r, u, s%R, s%u)
  vx_rel = vx_rel - u*(xcomp - xcore)/r
  vy_rel = vy_rel - u*(ycomp - ycore)/r
end if

v_rel = sqrt( vx_rel**2 + vy_rel**2 )

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
            e_inj, omegadot, vx_rel, vy_rel, e, beta, Q

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

! Get separation and relative velocity
r = sqrt( (xcore - xcomp)**2 + (ycore - ycomp)**2 )
call get_relative_velocity(s, r, xcore, ycore, xcomp, ycomp, &
                           vxcore, vycore, vxcomp, vycomp, &
                           vx_rel, vy_rel, v_rel)

! Interpolate needed quantities from MESA model
call interpolate(r, rho, s%R, s%rho)
call interpolate(r, menc, s%R, s%m)

! Get accretion rate using MR15
! Ignore fd and edot from this model
call mr15(s, r, mns, rho, v_rel, mdot, fd, edot)

! Get drag and power using Kim & Kim (2010)
call kim2010(s, r, mns, rho, v_rel, fd, edot)

! Get spinup rate
call get_spinup_rate(s, mns, macc, omega, omegadot, e, beta, Q)

dxdt(1)  = vxcore
dxdt(2)  = vycore
dxdt(3)  = vxcomp
dxdt(4)  = vycomp
dxdt(5)  = standard_cgrav*mns*(xcomp - xcore)/r**3
dxdt(6)  = standard_cgrav*mns*(ycomp - ycore)/r**3
dxdt(7)  = standard_cgrav*menc*(xcore - xcomp)/r**3 - fd*vx_rel/(mns*v_rel)
dxdt(8)  = standard_cgrav*menc*(ycore - ycomp)/r**3 - fd*vy_rel/(mns*v_rel)
dxdt(9)  = omegadot
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
real(dp), intent(in)                     :: r, mns, v_rel, e_inj

real(dp) :: rho, Ra, kernel_sum, junk
integer  :: n, j1, j2, j
real(dp), allocatable :: temp_heat(:)

! Get the accretion radius and the indices of zones we will modify

n = size(s%R)
call interpolate(r, rho, s%R, s%rho)
call get_hl_accretion(mns, rho, v_rel, junk, junk, junk, Ra)
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

subroutine kim2007(s, r, M, rho, v_rel, fd_out, edot_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r, M, rho, v_rel
real(dp, intent(out)                  :: fd_out, edot_out

real(dp) :: cs, mach

call interpolate(r, cs, s%R, s%csound)
mach = v_rel / cs

if (mach < 1.0) then
  fd_out = 0.7706*LOG((1+mach)/(1.0004 - 0.9185*mach)) - 1.473*mach
else if (mach >= 1.0 .and. mach < 4.4) then
  fd_out = LOG(330*(r*(mach-0.71)**5.72)/(R0*mach**9.58))
else
  fd_out = LOG(r / (R0*(0.11*mach + 1.65)))
end if
! NOTES: units on fd_out don't seem right -- comes out as g cm s-2, should be
! g cm-1 s-2
fd_out = fd_out * 4*pi*rho*(standard_cgrav*M/v_rel)**2

edot_out = fd_out * v_rel

return
end subroutine kim2007

! ------------------------------------------------------------------------------

! Kim & Kim 2010 drag force calculation

subroutine kim2010(s, r, M, rho, v_rel, fd_out, edot_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r, M, rho, v_rel
real(dp), intent(out)                 :: fd_out, edot_out

integer  :: k
real(dp) :: i_var, beta, eta_b, cd
real(dp) :: cs, mach, mdot_hl, fd_hl, edot_hl, Ra

call get_hl_accretion(M, rho, v_rel, mdot_hl, fd_hl, edot_hl, Ra)
call interpolate(r, cs, s%R, s%csound)
mach = v_rel / cs

beta  = standard_cgrav*M/(r*cs**2)   !dimensionless
eta_b = beta/(mach**2 - 1)           !dimensionless
cd    = 0.002

if (mach < 1.01) then
  i_var = 0.7706*LOG((1+mach)/(1.0004 - 0.9185*mach)) - 1.473*mach
else if (mach >= 1.01 .and. mach < 4.4) then
  i_var = LOG(330*(r*(mach-0.71)**5.72)/(Ra*mach**9.58))
else
  i_var = LOG(r / (Ra*(0.11*mach + 1.65)))
end if

if (eta_b > 0.1 .and. mach > 1.01) then
  fd_out = -cd*0.7*4*pi*rho* &
           (1 + 0.46*(beta**1.1)/(mach**2 - 1)**0.11)* &
           (standard_cgrav*M**2)/((v_rel**2)*(eta_b**0.5))
else
  fd_out = -cd*4*pi*rho*(standard_cgrav*M**2)*i_var/(v_rel**2)
end if

edot_out = fd_out * v_rel

return
end subroutine kim2010

! ------------------------------------------------------------------------------

! MR15 accretion and drag force calculation

subroutine mr15(s, r, M, rho, v_rel, mdot_out, fd_out, edot_out)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: r, M, rho, v_rel
real(dp), intent(out)                 :: mdot_out, fd_out, edot_out

real(dp), parameter :: f1  = 1.91791946,  f2  = -1.52814698, f3 = 0.75992092, &
                       mu1 = -2.14034214, mu2 = 1.94694764,  mu3 = 1.19007536,&
                       mu4 = 1.05762477
real(dp)            :: rsc, eps_rho, fd_mr15_ratio, mdot_mr15_ratio
real(dp)            :: mdot_hl, fd_hl, edot_hl, Ra, mdot_edd, mdot_hyper

call get_hl_accretion(M, rho, v_rel, mdot_hl, fd_hl, edot_hl, Ra)
call edd_and_hyper(M, mdot_edd, mdot_hyper)
call interpolate(r, rsc, s%R, s%scale_height)

eps_rho         = Ra / rsc
fd_mr15_ratio   = f1 + f2*eps_rho + f3*eps_rho**2
mdot_mr15_ratio = 10**(mu1 + mu2/(1 + mu3*eps_rho + mu4*eps_rho**2))
mdot_out        = eta * mdot_mr15_ratio * mdot_hl
fd_out          = eta * fd_mr15_ratio * fd_hl

if ((mdot_edd < mdot_out) .and. (mdot_out < mdot_hyper)) then
  mdot_out = mdot_edd
  fd_out   = mdot_edd * v_rel 
end if

edot_out = fd_out * v_rel

return
end subroutine mr15

! ------------------------------------------------------------------------------

! Hoyle-Lyttleton accretion radius, accretion rate, drag force, and energy
! loss rate

subroutine get_hl_accretion(M, rho, v_rel, mdot_hl, fd_hl, edot_hl, Ra)

real(dp), intent(in)                  :: M, rho,, v_rel
real(dp), intent(out)                 :: mdot_hl, fd_hl, edot_hl, Ra

Ra      = 2*standard_cgrav*M/v_rel**2       !cm
mdot_hl = pi*(Ra**2)*rho*v_rel              !g/s
fd_hl   = mdot_hl*v_rel                     !dyne (g cm s^-2)
edot_hl = fd_hl*v_rel                       !erg/s

return
end subroutine get_hl_accretion

! ------------------------------------------------------------------------------

! Eddington and hypercritical accretion rate calculation

subroutine edd_and_hyper(M, mdot_edd, mdot_hyper)

real(dp), intent(in)  :: M
real(dp), intent(out) :: mdot_edd, mdot_hyper

! mass is in g, radius is in cm, time is in s

mdot_edd   = 3.5d-8*(M/Msun)*(0.34/op_const)*Msun/secyer      !gm/s
mdot_hyper = 8.9d-5*((op_const/0.34)**(-0.73))*Msun/secyer    !gm/s

return
end subroutine edd_and_hyper

! ------------------------------------------------------------------------------

! Solve equations 34/35 of Holgado et al. for the ellipticity and spin parameter
! of a Maclaurin spheroid with given mass and spin frequency. The spin parameter
! is not permitted to exceed the value corresponding to secular instability.
! NOTES: if omega wants to drive beta to be > beta_sec, shouldn't this limit
! omega as well as e? If so, what's the best way to handle this?

subroutine get_e_beta_given_omega(M, omega, e, beta)

real(dp), intent(in)  :: M, omega
real(dp), intent(out) :: e, beta

real(dp), parameter :: tol = 1.e-6, kappan = 1.
integer, parameter  :: max_iter = 100
real(dp)            :: rhobar, qn, ofctn_val, emin, emax, enxt, onxt
integer             :: n_iter

n_iter    = 0
onxt      = huge(1.)
qn        = (1 - n_poly/5)*kappan               ! Holgado et al. eqn 36
rhobar    = 3*M/(4*pi*R0**3)                    ! Holgado et al. after eqn 35
ofctn_val = omega**2 / (2*pi*standard_cgrav*rhobar/qn)
                                                ! Holgado et al. eqn 35
emin      = 1.e-9
emax      = 0.817   ! NOTES: does this correspond to beta = beta_sec?
                    ! Got the 0.817 number from original code, and it
                    ! corresponds to beta = 0.14, which seems to match
                    ! figure 3 of Holgado et al.

! bisection
do while ((abs(onxt) > tol) .and. (n_iter < max_iter))
  enxt = (emin + emax) / 2
  onxt = sqrt(1-e**2)/e**3*(3-2*e**2)*asin(e) - 3*(1-e**2)/e**2 - ofctn_val
                                                ! Holgado et al. eqn 35
  if (onxt > 0.) then
    emax = enxt
  else
    emin = enxt
  n_iter = n_iter + 1
enddo

e    = enxt
beta = 3./(2*e**2)*(1 - e*(1 - e**2)**0.5/asin(e)) - 1
                                                ! Holgado et al. eqn 34

return
end subroutine get_e_beta_given_omega

! ------------------------------------------------------------------------------

! Get the rate of change of the spin frequency

subroutine get_spinup_rate(s, M, Macc, omega, omegadot, e, beta, Q)

type (star_info), pointer, intent(in) :: s
real(dp), intent(in)                  :: M, Macc
real(dp), intent(out)                 :: omegadot, e, beta, Q

real(dp)            :: Rbar, Req, Rz, Imom, Qtb, Qmax, Nacc, Ngw

! Given M and omega, solve for e and beta via equations 34/35 of Holgado et al.
call get_e_beta_given_omega(M, omega, e, beta)

Rbar = R0 * (asin(e)/e * (1 - e**2)**(1./6.) * (1 - beta))**(-n_poly/(3-n_poly))
                                             ! Holgado et al. eqn 38
Req  = Rbar / (1 - e**2)**(1./6.)            ! Holgado et al. eqn 37
Qtb  = sqrt((5./32.)*clight**5/(standard_cgrav*omega**5)*Mdot * &
            sqrt(standard_cgrav*M*Req))      ! Holgado et al. eqn 19
Qmax = 4.e39*Macc/(0.05*Msun)                ! Holgado et al. eqn 40

if (Qtb > Qmax) then
  Rz       = Req*(1 - e)                     ! from defn of ellipticity
  Imom     = (M/5.) * (Rz**2 + Req**2)       ! Holgado et al. after eqn 16
  Q        = Qmax                            ! Holgado et al. eqn 41
                                             ! NOTES: inconsistency with eqn 16?
  Nacc     = Mdot*sqrt(standard_cgrav*M*Req) ! Holgado et al. eqn 4
  Ngw      = -(32./5.)*standard_cgrav*omega**5*Q**2/clight**5
                                             ! Holgado et al. eqn 14/15
  omegadot = (Nacc + Ngw) / Imom             ! Holgado et al. eqn 17
else
  Q        = Qtb                             ! Holgado et al. eqn 41
  omegadot = 0.
endif

return
end subroutine get_spinup_rate

! ------------------------------------------------------------------------------

! Update all the s%xtra() quantities used to track the model, as well as any
! quantities that will be stored to history files

!TODO
subroutine update_quantities()

return
end subroutine update_quantities

! ------------------------------------------------------------------------------

! Print some useful information to stdout

!TODO
subroutine print_info()

print *, 'model number                  = ', s% model_number
print *, 'orbital separation            = ', s%xtra(ia)/Rsun
print *, 'omega                         = ', s%xtra(iomega)
print *, 'Q                             = ', s%xtra(iQ)
print *, 'strain                        = ', s%xtra(istrain)

select case (prescription)
  case (ENERGY_PRESCRIPTION)
    print *, 'mach                  = ', v/csound
    print *, 'edot                  = ', edot
  case (FORCE_PRESCRIPTION)
    print *, 'mach                  = ', v_rel/csound
    print *, 'edot                  = ', -edot_br
    ! print *, 'fd_br                  = ', fd_br
    print *, 'xcore_curr            = ', s%xtra(ixcore)
    print *, 'ycore_curr            = ', s%xtra(iycore)
    print *, 'xcomp_curr            = ', s%xtra(ixcomp)
    print *, 'ycomp_curr            = ', s%xtra(iycomp)
    print *, 'vxcore_curr           = ', s%xtra(ivxcore)
    print *, 'vycore_curr           = ', s%xtra(ivycore)
    print *, 'vxcomp_curr           = ', s%xtra(ivxcomp)
    print *, 'vycomp_curr           = ', s%xtra(ivycomp)
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
      call Req_and_beta(s%xtra(ie), Req, Rbar, beta)
      omega_next = s%xtra(iomega)                             !Hz
      s%xtra(iQtb) = sqrt((5*(clight**5)*mdot*sqrt(standard_cgrav*s%xtra(iM_ns)*Req))/(32*standard_cgrav*((s%xtra(iomega))**5))) !g cm^2
      Qtb_next = s%xtra(iQtb)                                 !g cm^2
      Qmax_next = s%xtra(iQmax)                               !g cm^2
      s%xtra(Q_curr) = s%xtra(iQtb)                           !g cm^2
      Q_next = s%xtra(iQ)                                     !g cm^2
      e_next = s%xtra(ie)                                     !dimensionless
      aa_next = s%xtra(iaa)                                   !cm
      bb_next = s%xtra(ibb)                                   !cm
      mom_inert_next = s%xtra(imom_inert)                     !gm cm^2

      decay_coeff = 1                                             !dimensionless

      s%xtra(istrain) = 2*standard_cgrav*((s%xtra(iomega))**2)*(s%xtra(iQ))/(D*(clight**4)) !dimensionless
      strain_next = s%xtra(istrain)                            !dimensionless

   else 
      if (s%xtra(ie) > 0.817) then 
         s%xtra(ie) = 0.817                               !dimensionless
      end if

      call Req_and_beta(s%xtra(ie), Req, Rbar, beta)
      omega_func_curr = (mdot*sqrt(standard_cgrav*s%xtra(iM_ns)*Req) - (32*standard_cgrav*(s%xtra(iomega)**5)*(s%xtra(iQ_)**2))/(5*clight**5))/s%xtra(imom_inert)
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

   if (s% xtra(ia) .lt. s%R(s%nz)) then
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

   vals(1) = s%xtra(ia)
   vals(2) = efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
   vals(3) = mdot
   vals(4) = azone
   vals(5) = r1
   vals(6) = r2
   vals(7) = s%xtra(iM_ns)
   vals(8) = s%xtra(iM_acc)
   vals(9) = -ebind
   vals(10) = eorb_change
   vals(11) = s%xtra(iQmax)
   vals(12) = s%xtra(iQtb)
   vals(13) = s%xtra(iQ)
   vals(14) = s%xtra(iomega)
   vals(15) = Req
   vals(16) = Rbar
   vals(17) = beta
   vals(18) = s%xtra(ie)
   vals(19) = s%xtra(iaa)
   vals(20) = s%xtra(ibb)
   vals(21) = s%xtra(imom_inert)
   vals(22) = s%xtra(istrain)
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
