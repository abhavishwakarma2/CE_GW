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
!    - torque balance quadrupole moment (Qtb)
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
                      imom_inert = 9,  istrain    = 10, &
                      ixcore     = 11, iycore     = 12, &
                      ixcomp     = 13, iycomp     = 14, &
                      ivxcore    = 15, ivycore    = 16, &
                      ivxcomp    = 17, ivycomp    = 18  

! input parameters from inlists
real(dp) :: M_ns_initial, omega_initial, e_initial, D, R0, r1, r2, &
            op_const, eta, efactor, M_crust, n_poly, beta_sec, menc_global, &
            omega_env_factor

!constant parameters
real(dp), parameter ::  ia_init_rsun = 261, Qmax_const = 4*1d39, &
                        v_rel_const = 0, aorb_init_ratio = 9.0d-1

integer  :: prescription, drag_prescription

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
   drag_prescription = int(s%x_ctrl(15)) !drag prescription to use for force prescription

   ! Initial conditions for companion orbit

   if (s%model_number == 1) then

      s%xtra(ia)         = ia_init_rsun*Rsun                            !semimajor axis (cm)
      s%xtra(iM_ns)      = M_ns_initial                                 !NS mass (g)
      s%xtra(iM_acc)     = 0.                                           !gm ! accreted mass?
      s%xtra(iomega)     = omega_initial                                !NS spin freq (Hz)
      !s%xtra(ie)         = e_initial                                   !NS ellipticity
      s%xtra(iQmax)      = Qmax_const*(s%xtra(iM_acc)/M_crust)          !NS quadrupole moment
                                                                        !(g cm^2)
      !s%xtra(iaa)        = R0*(1+s%xtra(ie)/2)                         !NS major axis (cm)
      !s%xtra(ibb)        = R0*(1-s%xtra(ie)/2)                         !NS minor axis (cm)
      s%xtra(imom_inert) = s%xtra(iM_ns)*(2*R0**2)/5
                                                                        !NS moment of inertia
                                                                        !(g cm^2)
      
      call interpolate(s%xtra(ia), menc_global, s%R, s%m)
      s%xtra(iomega_env) = omega_env_factor * &
                           SQRT(standard_cgrav*(s%xtra(iM_ns) + menc_global)/ s%xtra(ia)**3) !along z axis, envelope angular velocity (Hz)

      if (prescription == FORCE_PRESCRIPTION) then
         s%xtra(ixcore)  = 0.
         s%xtra(iycore)  = 0.
         s%xtra(ixcomp)  = s%xtra(ia)
         s%xtra(iycomp)  = 0.
         s%xtra(ivxcore) = 0.
         s%xtra(ivycore) = 0.
         s%xtra(ivxcomp) = 0.
         s%xtra(ivycomp) = SQRT(standard_cgrav*(s%xtra(iM_ns) + menc_global)/s%xtra(ia))
      end if

   end if

   ! Initialize time integration

   tlocal  = s%time 
   dtlocal = s%dt_next / Nsteps
   s%extra_heat(:)%val = 0.

   ! Integrate

   do i = 1, Nsteps
      select case (prescription)
         case (ENERGY_PRESCRIPTION)
            call advance_energy_prescription(s, tlocal, dtlocal)
         case (FORCE_PRESCRIPTION)
            call advance_force_prescription(s, tlocal, dtlocal)
      end select

      tlocal = tlocal + dtlocal
   enddo

   ! Print diagnostic info

   call print_info(s)

   return
end subroutine evolve_companion_and_inject_energy

! ------------------------------------------------------------------------------

! Integrate the companion orbit and other properties through one timestep using
! the energy prescription.

subroutine advance_energy_prescription(s, t, dt)

   type (star_info), pointer, intent(inout) :: s
   real(dp), intent(in)                     :: t, dt

   real(dp) :: x(5), dxdt(5), e_inj, a_orb, v_rel, menc

   x = (/ s%xtra(ia), s%xtra(iomega), s%xtra(iM_ns), s%xtra(iM_acc) , 0d0 /)

   call energy_prescrip_rhs(s, t, x, dxdt)

   x = x + dt*dxdt

   s%xtra(ia)     = x(1)
   s%xtra(iomega) = x(2)
   s%xtra(iM_ns)  = x(3)
   s%xtra(iM_acc) = x(4)
   e_inj          = x(5)

   a_orb = s%xtra(ia)

   call interpolate(a_orb, menc, s%R, s%m)
   call get_relative_velocity_energy_rhs(s, a_orb, s%xtra(iM_ns), menc, v_rel)
   call add_energy_to_mesa(s, s%xtra(ia), s%xtra(iM_ns), v_rel, e_inj)

   return
end subroutine advance_energy_prescription

! ------------------------------------------------------------------------------

! Integrate the companion orbit and other properties through one timestep using
! the force prescription.

subroutine advance_force_prescription(s, t, dt)

   type (star_info), pointer, intent(inout) :: s
   real(dp), intent(in)                     :: t, dt

   real(dp) :: x(12), dxdt(12), e_inj, r, v_rel
   real(dp) :: vx_rel, vy_rel 

   ! Set up solution vector.

   x = (/ s%xtra(ixcore),  s%xtra(iycore),  s%xtra(ixcomp),  s%xtra(iycomp), &
         s%xtra(ivxcore), s%xtra(ivycore), s%xtra(ivxcomp), s%xtra(ivycomp), &
         s%xtra(iomega),  s%xtra(iM_ns),   s%xtra(iM_acc),  0d0 /)

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

   s%xtra(ia) = r !A: added this to keep track of separation in the force prescription as well
   call get_relative_velocity_force_rhs(s, r, s%xtra(ixcore), s%xtra(iycore), &
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

subroutine energy_prescrip_rhs(s, t, x, dxdt)

   real(dp), intent(in)                     :: t, x(:)
   real(dp), intent(out)                    :: dxdt(:)
   type (star_info), pointer, intent(inout) :: s

   real(dp) :: a_orb, omega, mns, macc, v_rel, rho, menc, mdot, fd, edot, &
               e_inj, omegadot, e, beta, Q, strain

   !print *, 'energy_prescrip_rhs called'
   a_orb    = x(1)
   omega    = x(2)
   mns      = x(3)
   macc     = x(4)
   e_inj    = x(5)

   ! Interpolate needed quantitites from MESA model
   call interpolate(a_orb, rho, s%R, s%rho)
   call interpolate(a_orb, menc, s%R, s%m)

   ! Relative velocity
   call get_relative_velocity_energy_rhs(s, a_orb, mns, menc, v_rel)
   
   ! Get accretion rate using MR15
   ! Ignore fd and edot from this model
   call mr15(s, a_orb, mns, rho, v_rel, mdot, fd, edot)

   ! Get drag and power using Kim & Kim (2010)
   if (drag_prescription == MR15_DRAG) then
      print *, 'Using MR15 drag prescription'
      call mr15(s, a_orb, mns, rho, v_rel, mdot, fd, edot)
   else if (drag_prescription == KIM07_DRAG) then
      print *, 'Using Kim 2007 drag prescription'
      call kim2007(s, a_orb, mns, rho, v_rel, fd, edot)
   else if (drag_prescription == KIM10_DRAG) then
      print *, 'Using Kim & Kim 2010 drag prescription'
      call kim2010(s, a_orb, mns, rho, v_rel, fd, edot)
   end if

   ! Get spinup rate
   call get_spinup_rate(s, mns, macc, mdot, omega, omegadot, e, beta, Q)

   dxdt(1) = -1 * edot * (((standard_cgrav*mns*menc)/(2 * a_orb**2) - (standard_cgrav*mns/(2*a_orb))*(4 * pi * a_orb**2 * rho))**(-1))
   dxdt(2) = omegadot
   dxdt(3) = mdot
   dxdt(4) = mdot
   dxdt(5) = edot

   ! Get strain
   call get_strain(s, Q, omega, D, strain)

   return
end subroutine energy_prescrip_rhs

! ------------------------------------------------------------------------------

! Return the right-hand-side of the differential equations describing the
! companion's orbit and spin for the force prescription.

subroutine force_prescrip_rhs(t, x, dxdt, s)

   real(dp), intent(in)                     :: t, x(:)
   real(dp), intent(out)                    :: dxdt(:)
   type (star_info), pointer, intent(inout) :: s

   real(dp) :: xcore, ycore, xcomp, ycomp, vxcore, vycore, vxcomp, vycomp, &
               omega, mns, macc, r, v_rel, rho, menc, mdot, fd, edot, &
               e_inj, omegadot, vx_rel, vy_rel, e, beta, Q, strain

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
   call get_relative_velocity_force_rhs(s, r, xcore, ycore, xcomp, ycomp, &
                              vxcore, vycore, vxcomp, vycomp, &
                              vx_rel, vy_rel, v_rel)

   ! Interpolate needed quantities from MESA model
   call interpolate(r, rho, s%R, s%rho)
   call interpolate(r, menc, s%R, s%m)

   ! Get accretion rate using MR15
   ! Ignore fd and edot from this model
   call mr15(s, r, mns, rho, v_rel, mdot, fd, edot)

   ! Get drag and power using Kim & Kim (2010)
   if (drag_prescription == MR15_DRAG) then
      print *, 'Using MR15 drag prescription'
      call mr15(s, r, mns, rho, v_rel, mdot, fd, edot)
   else if (drag_prescription == KIM07_DRAG) then
      print *, 'Using Kim 2007 drag prescription'
      call kim2007(s, r, mns, rho, v_rel, fd, edot)
   else if (drag_prescription == KIM10_DRAG) then
      print *, 'Using Kim & Kim 2010 drag prescription'
      call kim2010(s, r, mns, rho, v_rel, fd, edot)
   end if

   ! Get spinup rate
   call get_spinup_rate(s, mns, macc, mdot, omega, omegadot, e, beta, Q)

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

   ! Get strain
   call get_strain(s, Q, omega, D, strain)

   return
end subroutine force_prescrip_rhs

! ------------------------------------------------------------------------------

! Compute the relative velocity of the NS and the donor envelope.

subroutine get_relative_velocity_energy_rhs(s, r, mns, menc, v_rel)

   type (star_info), pointer, intent(in) :: s
   real(dp), intent(in)                  :: r, mns, menc
   real(dp), intent(out)                 :: v_rel

   v_rel = SQRT(standard_cgrav*(mns + menc)/r) - s%xtra(iomega_env)*r

   return
end subroutine get_relative_velocity_energy_rhs

subroutine get_relative_velocity_force_rhs(s, r, xcore, ycore, xcomp, ycomp, &
                                 vxcore, vycore, vxcomp, vycomp, &
                                 vx_rel, vy_rel, v_rel)

   type (star_info), pointer, intent(in) :: s
   real(dp), intent(in)                  :: r, xcore, ycore, xcomp, ycomp, &
                                            vxcore, vycore, vxcomp, vycomp
   real(dp), intent(out)                 :: vx_rel, vy_rel, v_rel

   real(dp)                              :: u

   vx_rel = vxcomp - vxcore + s%xtra(iomega_env)*(ycomp - ycore) 
   vy_rel = vycomp - vycore - s%xtra(iomega_env)*(xcomp - xcore)

   if (r <= v_rel_const*Rsun) then
      call interpolate(r, u, s%R, s%u)
      vx_rel = vx_rel - u*(xcomp - xcore)/r
      vy_rel = vy_rel - u*(ycomp - ycore)/r
   end if

   v_rel = sqrt( vx_rel**2 + vy_rel**2 )

   return
end subroutine get_relative_velocity_force_rhs


! ------------------------------------------------------------------------------

! Given the NS position and the amount of energy to inject, add thermal energy
! to the MESA model according to the Bronner prescription.

subroutine add_energy_to_mesa(s, r, mns, v_rel, e_inj)

   type (star_info), pointer, intent(inout) :: s
   real(dp), intent(in)                     :: r, mns, v_rel, e_inj

   real(dp) :: rho, Ra, kernel_sum, junk
   integer  :: n, j1, j2, j
   real(dp), allocatable :: temp_heat(:)

   ! Get the accretion radius and the indices of zones we will 
   Ra = 2*standard_cgrav*mns/v_rel**2

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

   do j = j2, j1
   ! NOTES: not sure I have translated this correctly.
   ! I think it is correct
   temp_heat(j) = exp(-((r-s%R(j))/Ra)**2)
   kernel_sum = kernel_sum + temp_heat(j)*s%dm(j)

   enddo

   do j = j2, j1
   s%extra_heat(j)%val = s%extra_heat(j)%val + &
                           efactor * temp_heat(j) * e_inj / (kernel_sum * s%dt_next)
   enddo

   print *, 'Energy injected into MESA: ', e_inj, 'erg'
   print *, 'sum of injected energy in kernels: ', SUM(s%extra_heat(j2:j1)%val * s%dm(j2:j1) * s%dt_next), 'erg'

   deallocate(temp_heat)

   return
end subroutine add_energy_to_mesa

! ------------------------------------------------------------------------------

! Kim 2007 drag force calculation

subroutine kim2007(s, r, M, rho, v_rel, fd_out, edot_out)

   type (star_info), pointer, intent(in) :: s
   real(dp), intent(in)                  :: r, M, rho, v_rel
   real(dp), intent(out)                  :: fd_out, edot_out

   real(dp) :: cs, mach, junk, Ra

   call get_hl_accretion(M, rho, v_rel, junk, junk, junk, Ra)
   call interpolate(r, cs, s%R, s%csound)
   mach = v_rel / cs

   if (mach < 1.0) then
      fd_out = 0.7706*LOG((1+mach)/(1.0004 - 0.9185*mach)) - 1.473*mach
   else if (mach >= 1.0 .and. mach < 4.4) then
      fd_out = LOG(330*(r*(mach-0.71)**5.72)/(Ra*mach**9.58))
   else
      fd_out = LOG(r / (Ra*(0.11*mach + 1.65)))
   end if
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
   real(dp) :: i_var, beta, eta_b, cd, rho_nl
   real(dp) :: cs, mach, junk, Ra

   call get_hl_accretion(M, rho, v_rel, junk, junk, junk, Ra)
   call interpolate(r, cs, s%R, s%csound)
   mach = v_rel / cs

   beta  = standard_cgrav*M/(r*cs**2)   !dimensionless
   eta_b = beta/(mach**2 - 1)           !dimensionless
   cd    = 1

   rho_nl = rho * (1 + 0.46*beta**1.1 / (mach**2 -1)**0.11)

   if (mach < 1.01) then
      i_var = 0.7706*LOG((1+mach)/(1.0004 - 0.9185*mach)) - 1.473*mach
   else if (mach >= 1.01 .and. mach < 4.4) then
      i_var = LOG(330*(r*(mach-0.71)**5.72)/(Ra*mach**9.58))
   else
      i_var = LOG(r / (Ra*(0.11*mach + 1.65)))
   end if

   if (eta_b > 0.1 .and. mach > 1.01) then
      fd_out = cd * (4 * pi * rho_nl * &
            ((standard_cgrav * M / v_rel)**2) ) * (0.7/eta_b**0.5)
   else
      fd_out = cd * (4 * pi * rho * &
            ((standard_cgrav * M / v_rel)**2) ) * i_var
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
                        mu1 = -2.14034214, mu2 = 1.94694764,  mu3 = 1.19007536, &
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

   real(dp), intent(in)                  :: M, rho, v_rel
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

   mdot_edd   = 3.5d-8*(M/(1.33*Msun))*(0.34/op_const)*Msun/secyer      !gm/s
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
   emax      = 0.817    !NOTES: does this correspond to beta = beta_sec?
                        !Got the 0.817 number from original code, and it
                        !corresponds to beta = 0.14, which seems to match
                        !figure 3 of Holgado et al.

   ! bisection
   do while ((abs(onxt) > tol) .and. (n_iter < max_iter))
   enxt = (emin + emax) / 2
   onxt = sqrt(1-e**2)/e**3*(3-2*e**2)*asin(e) - 3*(1-e**2)/e**2 - ofctn_val
                                                   ! Holgado et al. eqn 35
   if (onxt > 0.) then
      emax = enxt
   else
      emin = enxt
   end if

   n_iter = n_iter + 1
   end do

   e    = enxt
   beta = 3./(2*e**2)*(1 - e*(1 - e**2)**0.5/asin(e)) - 1
                                                         ! Holgado et al. eqn 34

   return
end subroutine get_e_beta_given_omega

! ------------------------------------------------------------------------------

! Get the rate of change of the spin frequency

subroutine get_spinup_rate(s, M, Macc, Mdot, omega, omegadot, e, beta, Q)

   type (star_info), pointer, intent(in) :: s
   real(dp), intent(in)                  :: M, Macc, Mdot, omega
   real(dp), intent(out)                 :: omegadot, e, beta, Q

   real(dp)            :: Rbar, Req, Rz, Imom, Qtb, Qmax, Nacc, Ngw

   ! Given M and omega, solve for e and beta via equations 34/35 of Holgado et al.
   call get_e_beta_given_omega(M, omega, e, beta)

   Rbar = R0 * (asin(e)/e * (1 - e**2)**(1./6.) * (1 - beta))**(-n_poly/(3-n_poly))
                                                         ! Holgado et al. eqn 38
   Req  = Rbar / (1 - e**2)**(1./6.)                     ! Holgado et al. eqn 37
   Qtb  = sqrt((5./32.)*clight**5/(standard_cgrav*omega**5)*Mdot * &
               sqrt(standard_cgrav*M*Req))               ! Holgado et al. eqn 19
   Qmax = Qmax_const*Macc/(M_crust)                  ! Holgado et al. eqn 40
   Imom     = (2.*M/5.) * R0**2                          ! Holgado et al. after eqn 16  !(M/5.) * (Rz**2 + Req**2)

   ! I think this should be if (Qtb < Qmax). Eqn 41 of the paper should probably have
   ! max(Qmax, Qtb) instead of min(Qmax, Qtb), if the text following the equation is 
   ! correct.
   if (Qtb > Qmax) then
      !Rz       = Req*(1 - e)                      ! from defn of ellipticity
      Q        = Qmax                              ! Holgado et al. eqn 41
                                                   ! NOTES: inconsistency with eqn 16?
      Nacc     = Mdot*sqrt(standard_cgrav*M*Req)   ! Holgado et al. eqn 4
      Ngw      = -(32./5.)*standard_cgrav*omega**5*Q**2/clight**5
                                                   ! Holgado et al. eqn 14/15
      omegadot = (Nacc + Ngw) / Imom               ! Holgado et al. eqn 17
   else
      Q        = Qtb                               ! Holgado et al. eqn 41
      omegadot = 0.
   endif

   s%xtra(imom_inert) = Imom
   s%xtra(iQmax) = Qmax
   s%xtra(iQtb)  = Qtb
   s%xtra(iQ)    = Q

   return
end subroutine get_spinup_rate

! ------------------------------------------------------------------------------

! Update all the s%xtra() quantities used to track the model, as well as any
! quantities that will be stored to history files

! ------------------------------------------------------------------------------

! Print some useful information to stdout

!TODO
subroutine print_info(s)

   type (star_info), pointer, intent(in) :: s

   print *, 'model number                  = ', s% model_number
   print *, 'orbital separation            = ', s%xtra(ia)/Rsun
   ! print *, 'omega                         = ', s%xtra(iomega)
   ! print *, 'Q                             = ', s%xtra(iQ)
   ! print *, 'strain                        = ', s%xtra(istrain)
   print *, 'total injected energy = ', &
               efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
   print *, 'M_acc                         = ', s%xtra(iM_acc)/Msun
   select case (prescription)
   case (ENERGY_PRESCRIPTION)
      print *, '' 
   case (FORCE_PRESCRIPTION)
      print *, 'xcore_curr            = ', s%xtra(ixcore)
      print *, 'ycore_curr            = ', s%xtra(iycore)
      print *, 'xcomp_curr            = ', s%xtra(ixcomp)
      print *, 'ycomp_curr            = ', s%xtra(iycomp)
      print *, 'vxcore_curr           = ', s%xtra(ivxcore)
      print *, 'vycore_curr           = ', s%xtra(ivycore)
      print *, 'vxcomp_curr           = ', s%xtra(ivxcomp)
      print *, 'vycomp_curr           = ', s%xtra(ivycomp)
   end select

   print *, '##############################################################################'
   print *, 'ONE TIMESTEP DONE'
   print *, '##############################################################################'

   return
end subroutine print_info

! ------------------------------------------------------------------------------

!TODO
!DONE: to evaluate the gravitational wave strain
subroutine get_strain(s, Q, omega, D, strain)

   type (star_info), pointer, intent(in) :: s
   real(dp), intent(in)                  :: Q, omega, D
   real(dp), intent(out)                 :: strain

   !evaluate the strain using the current value of omega and Q
   strain = 2 * standard_cgrav * Q * omega**2 / (D * clight**4)
   
   !update the s%xtra(istrain) quantity
   s%xtra(istrain) = strain

end subroutine get_strain

! ==============================================================================

! Temperature and pressure at the surface of the star, and their derivatives with
! respect to luminosity, radius, mass, and opacity. These are used in the
! surface boundary condition for the MESA model. 

subroutine other_surface_PT(id, &
                  skip_partials, &
                  lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
                  lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap, ierr)

   use const_def, only: dp
   integer, intent(in) :: id
   logical, intent(in) :: skip_partials
   real(dp), intent(out) :: lnT_surf, dlnT_dL, dlnT_dlnR, dlnT_dlnM, dlnT_dlnkap, &
            lnP_surf, dlnP_dL, dlnP_dlnR, dlnP_dlnM, dlnP_dlnkap
   integer, intent(out) :: ierr
   type (star_info), pointer :: s

   real(dp) :: g
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return

   lnT_surf = log(pow(s%L(1) / (4d0 * pi * boltz_sigma * pow2(s%R(1))), 0.25d0))
   dlnT_dL = lnT_surf / s%L(1)
   dlnT_dlnR = -0.5d0
   dlnT_dlnM = 0d0
   dlnT_dlnkap = 0d0

   g = s%M(1) / pow2(s%R(1))
   lnP_surf = log(g / s%opacity(1))
   dlnP_dL = 0d0
   dlnP_dlnR = -2d0
   dlnP_dlnM = 1d0
   dlnP_dlnkap = -1d0

end subroutine other_surface_PT

! ==============================================================================

! Utility routines

! Find location in a monotonically varying array corresponding to a given
! value, or 0 or N if out of range.

subroutine hunt(xx, n, x, jlo)  !why are there numbers in the beginning of some lines?

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
   s% other_surface_pt => other_surface_pt
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

   if (prescription == ENERGY_PRESCRIPTION) then
      how_many_extra_history_columns = 27
   else if (prescription == FORCE_PRESCRIPTION) then
      !we also want to store the current positions and velocities of the NS and companion, so we need 8 more columns for those
      how_many_extra_history_columns = 37
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
   names(1) = 'aorb'
   names(2) = 'Einj_per_dt'
   names(3) = 'M_ns'
   names(4) = 'M_enc'
   names(5) = 'M_acc'
   names(6) = 'rho'
   names(7) = 'Qmax'
   names(8) = 'Qtb'
   names(9) = 'Q'
   names(10) = 'omega'
   names(11) = 'mom_inert'
   names(12) = 'strain'
   names(13) = 'scale_height'
   names(14) = 'csound'
   names(15) = 'lnPgas'
   names(16) = 'lnT'
   names(17) = 'u_k'
   if (prescription == FORCE_PRESCRIPTION) then
      names(18) = 'xcore'
      names(19) = 'ycore'
      names(20) = 'xcomp'
      names(21) = 'ycomp'
      names(22) = 'vxcore'
      names(23) = 'vycore'
      names(24) = 'vxcomp'
      names(25) = 'vycomp'
      names(26) = 'vx_rel'
      names(27) = 'vy_rel'
      names(28) = 'v_rel'
      names(29) = 'mdot_hl'
      names(30) = 'fd_hl'
      names(31) = 'edot_hl'
      names(32) = 'Ra'
      names(33) = 'mdot_mr15'
      names(34) = 'fd_mr15'
      names(35) = 'edot_mr15'
      names(36) = 'fd'
      names(37) = 'edot'
   else if (prescription == ENERGY_PRESCRIPTION) then
      names(18) = 'v_rel'
      names(19) = 'mdot_hl'
      names(20) = 'fd_hl'
      names(21) = 'edot_hl'
      names(22) = 'Ra'
      names(23) = 'mdot_mr15'
      names(24) = 'fd_mr15'
      names(25) = 'edot_mr15'
      names(26) = 'fd'
      names(27) = 'edot'

   end if

   vals(1) = s%xtra(ia)
   vals(2) = efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz)) * s%dt_next
   vals(3) = s%xtra(iM_ns)
   call interpolate(s%xtra(ia), vals(4), s%R, s%m)
   vals(5) = s%xtra(iM_acc)
   call interpolate(s%xtra(ia), vals(6), s%R, s%rho)
   vals(7) = s%xtra(iQmax)
   vals(8) = s%xtra(iQtb)
   vals(9) = s%xtra(iQ)
   vals(10) = s%xtra(iomega)
   vals(11) = s%xtra(imom_inert)
   vals(12) = s%xtra(istrain)
   call interpolate(s%xtra(ia), vals(13), s%R, s%scale_height)
   call interpolate(s%xtra(ia), vals(14), s%R, s%csound)
   call interpolate(s%xtra(ia), vals(15), s%R, s%lnPgas)
   call interpolate(s%xtra(ia), vals(16), s%R, s%lnT)
   call interpolate(s%xtra(ia), vals(17), s%R, s%u)
   if (prescription == FORCE_PRESCRIPTION) then
      vals(18) = s%xtra(ixcore)
      vals(19) = s%xtra(iycore)
      vals(20) = s%xtra(ixcomp)
      vals(21) = s%xtra(iycomp)
      vals(22) = s%xtra(ivxcore)    
      vals(23) = s%xtra(ivycore)
      vals(24) = s%xtra(ivxcomp)
      vals(25) = s%xtra(ivycomp)
      call get_relative_velocity_force_rhs(s, vals(1), s%xtra(ixcore), s%xtra(iycore), s%xtra(ixcomp), s%xtra(iycomp), &
                                         s%xtra(ivxcore), s%xtra(ivycore), s%xtra(ivxcomp), s%xtra(ivycomp), &
                                         vals(26), vals(27), vals(28))
      call get_hl_accretion(vals(3), vals(6), vals(28), vals(29), vals(30), vals(31), vals(32))
      call mr15(s, vals(1), vals(3), vals(6), vals(28), vals(33), vals(34), vals(35))
      if (drag_prescription == MR15_DRAG) then
         vals(36) = vals(34)
         vals(37) = vals(35)
      else if (drag_prescription == KIM07_DRAG) then
         call kim2007(s, vals(1), vals(3), vals(6), vals(28), vals(36), vals(37))
      else if (drag_prescription == KIM10_DRAG) then
         call kim2010(s, vals(1), vals(3), vals(6), vals(28), vals(36), vals(37))
      end if
   else if (prescription == ENERGY_PRESCRIPTION) then
      call get_relative_velocity_energy_rhs(s, vals(1), vals(3), vals(4), vals(18))
      call get_hl_accretion(vals(3), vals(6), vals(18), vals(19), vals(20), vals(21), vals(22))
      call mr15(s, vals(1), vals(3), vals(6), vals(18), vals(23), vals(24), vals(25))
      if (drag_prescription == MR15_DRAG) then
         vals(26) = vals(24)
         vals(27) = vals(25)
      else if (drag_prescription == KIM07_DRAG) then
         call kim2007(s, vals(1), vals(3), vals(6), vals(18), vals(26), vals(27))
      else if (drag_prescription == KIM10_DRAG) then
         call kim2010(s, vals(1), vals(3), vals(6), vals(18), vals(26), vals(27))
      end if
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

   how_many_extra_profile_columns = 0
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
   
end subroutine data_for_extra_profile_columns

! ------------------------------------------------------------------------------

integer function how_many_extra_history_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   how_many_extra_history_header_items = 2
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

   names(1) = 'orbit_prescription'
   vals(1) = prescription

   names(2) = 'drag_prescription'
   vals(2) = drag_prescription

end subroutine data_for_extra_history_header_items

! ------------------------------------------------------------------------------

integer function how_many_extra_profile_header_items(id)
   integer, intent(in) :: id
   integer :: ierr
   type (star_info), pointer :: s
   ierr = 0
   call star_ptr(id, s, ierr)
   if (ierr /= 0) return
   if (prescription == FORCE_PRESCRIPTION) then
      how_many_extra_profile_header_items = 39
   else if (prescription == ENERGY_PRESCRIPTION) then
      how_many_extra_profile_header_items = 29
   end if
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

   names(1) = 'aorb'
   names(2) = 'M_ns'
   names(3) = 'orbit_prescription'
   names(4) = 'drag_prescription'
   names(5) = 'Einj_per_dt'
   names(6) = 'M_enc'
   names(7) = 'M_acc'
   names(8) = 'rho'
   names(9) = 'Qmax'
   names(10) = 'Qtb'
   names(11) = 'Q'
   names(12) = 'omega'
   names(13) = 'mom_inert'
   names(14) = 'strain'
   names(15) = 'scale_height'
   names(16) = 'csound'
   names(17) = 'lnPgas'
   names(18) = 'lnT'
   names(19) = 'u_k'
   if (prescription == FORCE_PRESCRIPTION) then
      names(20) = 'xcore'
      names(21) = 'ycore'
      names(22) = 'xcomp'
      names(23) = 'ycomp'
      names(24) = 'vxcore'
      names(25) = 'vycore'
      names(26) = 'vxcomp'
      names(27) = 'vycomp'
      names(28) = 'vx_rel'
      names(29) = 'vy_rel'
      names(30) = 'v_rel'
      names(31) = 'mdot_hl'
      names(32) = 'fd_hl'
      names(33) = 'edot_hl'
      names(34) = 'Ra'
      names(35) = 'mdot_mr15'
      names(36) = 'fd_mr15'
      names(37) = 'edot_mr15'
      names(38) = 'fd'
      names(39) = 'edot'
   else if (prescription == ENERGY_PRESCRIPTION) then
      names(20) = 'v_rel'
      names(21) = 'mdot_hl'
      names(22) = 'fd_hl'
      names(23) = 'edot_hl'
      names(24) = 'Ra'
      names(25) = 'mdot_mr15'
      names(26) = 'fd_mr15'
      names(27) = 'edot_mr15'
      names(28) = 'fd'
      names(29) = 'edot'
   end if

   vals(1) = s%xtra(ia)
   vals(2) = s%xtra(iM_ns)
   vals(3) = prescription
   vals(4) = drag_prescription
   vals(5) = efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz)) * s%dt_next
   call interpolate(s%xtra(ia), vals(6), s%R, s%m)
   vals(7) = s%xtra(iM_acc)
   call interpolate(s%xtra(ia), vals(8), s%R, s%rho)
   vals(9) = s%xtra(iQmax)
   vals(10) = s%xtra(iQtb)
   vals(11) = s%xtra(iQ)
   vals(12) = s%xtra(iomega)
   vals(13) = s%xtra(imom_inert)
   vals(14) = s%xtra(istrain)
   call interpolate(s%xtra(ia), vals(15), s%R, s%scale_height)
   call interpolate(s%xtra(ia), vals(16), s%R, s%csound)
   call interpolate(s%xtra(ia), vals(17), s%R, s%lnPgas)
   call interpolate(s%xtra(ia), vals(18), s%R, s%lnT)
   call interpolate(s%xtra(ia), vals(19), s%R, s%u)
   if (prescription == FORCE_PRESCRIPTION) then
      vals(20) = s%xtra(ixcore)
      vals(21) = s%xtra(iycore)
      vals(22) = s%xtra(ixcomp)
      vals(23) = s%xtra(iycomp)
      vals(24) = s%xtra(ivxcore)    
      vals(25) = s%xtra(ivycore)
      vals(26) = s%xtra(ivxcomp)
      vals(27) = s%xtra(ivycomp)
      call get_relative_velocity_force_rhs(s, vals(1), s%xtra(ixcore), s%xtra(iycore), s%xtra(ixcomp), s%xtra(iycomp), &
                                         s%xtra(ivxcore), s%xtra(ivycore), s%xtra(ivxcomp), s%xtra(ivycomp), &
                                         vals(28), vals(29), vals(30))
      call get_hl_accretion(vals(3), vals(6), vals(28), vals(31), vals(32), vals(33), vals(34))
      call mr15(s, vals(1), vals(3), vals(6), vals(28), vals(35), vals(36), vals(37))
      if (drag_prescription == MR15_DRAG) then
         vals(38) = vals(36)
         vals(39) = vals(37)
      else if (drag_prescription == KIM07_DRAG) then
         call kim2007(s, vals(1), vals(3), vals(6), vals(28), vals(38), vals(39))
      else if (drag_prescription == KIM10_DRAG) then
         call kim2010(s, vals(1), vals(3), vals(6), vals(28), vals(38), vals(39))
      end if
   else if (prescription == ENERGY_PRESCRIPTION) then
      call get_relative_velocity_energy_rhs(s, vals(1), vals(3), vals(4), vals(20))
      call get_hl_accretion(vals(3), vals(6), vals(18), vals(21), vals(22), vals(23), vals(24))
      call mr15(s, vals(1), vals(3), vals(6), vals(18), vals(25), vals(26), vals(27))
      if (drag_prescription == MR15_DRAG) then
         vals(28) = vals(26)
         vals(29) = vals(27)
      else if (drag_prescription == KIM07_DRAG) then
         call kim2007(s, vals(1), vals(3), vals(6), vals(18), vals(28), vals(29))
      else if (drag_prescription == KIM10_DRAG) then
         call kim2010(s, vals(1), vals(3), vals(6), vals(18), vals(28), vals(29))
      end if
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
