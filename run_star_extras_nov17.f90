! ***********************************************************************
!
!  README:
!  This module contains extra routines for running MESA star
!  simulations that involve common envelope evolution with a neutron star companion.
!  It includes functions for calculating accretion rates, drag forces, orbital evolution, 
!  and gravitational wave strain.
!
!  At every timestep, the current values (for that timestep) of the relevant variables are already known.
!  The module calculates the next values of these variables based on the chosen prescription, and saves 
!  them into the star's xtra array for use in the next timestep. The variables include:
!    - orbital separation (a)
!    - neutron star mass (M_ns)
!    - accreted mass (M_acc)
!    - orbital angular velocity (omega)
!    - quadrupole moment (Q)
!    - maximum quadrupole moment (Qmax)
!    - tidal bulge quadrupole moment (Qtb)
!    - eccentricity (e)
!    - x and y positions and velocities of both the primary core and the neutron star (a1, b1, c1, d1, a2, b2, c2, d2) (for elliptical orbits)
!    - moment of inertia of the envelope (mom_inert)
!    - strain (strain)
!  
!  Variables of the type *_curr refer to the current timestep values, while those with *_next refer to the next timestep values.
!  The prescription variable determines which orbital evolution prescription to use (1 for Holgado et al. 2018, 2 for Bronner et al. 2024).
!
!  Variable of the type s%xtra(...) refer to the star's xtra array where these variables are stored. MESA remembers these values across timesteps.
!
!  NOTE: All units are in cgs unless otherwise specified.
!        The first function that is called is extras_controls, which sets up the necessary variables and arrays.
!        Then, the inject_energy_and_evaluate_strain subroutine is called at every timestep to perform the calculations and update the variables.
!
! ************************************************************************
!
!  LAST MODIFIED : November 17, 2024
!
! ***********************************************************************
 
      module run_star_extras

         use star_lib
         use star_def
         use const_def
         use math_lib
         
         implicit none

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! DECLARATION OF GLOBAL VARIABLES !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !arrays
         real(dp), allocatable :: v(:), Ra(:), mdot_hl(:), fd_hl(:), edot_hl(:), Eorb(:), dEorb(:)                
         real(dp), allocatable :: eps_rho(:), fd_mr15_ratio(:), mdot_mr15_ratio(:), mdot_mr15(:), fd_mr15(:)
         real(dp), allocatable :: mdot_edd(:), mdot_hyper(:), mdot_arr(:), fd_arr(:), edot(:), f(:)

         !xtra variables - values for current step
         integer, parameter :: a_curr = 1, M_ns_curr = 2, M_acc_curr = 3, omega_env = 4                    
         integer, parameter :: omega_curr = 5, Q_curr = 6, Qmax_curr = 7, Qtb_curr = 8       
         integer, parameter :: e_curr = 9, aa_curr = 10, bb_curr = 11, mom_inert_curr = 12   
         integer, parameter :: strain_curr = 13                                              
         integer, parameter :: a1_curr = 14, b1_curr = 15, c1_curr = 16, d1_curr = 17
         integer, parameter :: a2_curr = 18, b2_curr = 19, c2_curr = 20, d2_curr = 21

         !values for next step
         real(dp) :: a_next, M_ns_next, M_acc_next, mdot_next, omega_next, Q_next
         real(dp) :: Qmax_next, Qtb_next, e_next, aa_next, bb_next, strain_next, mom_inert_next
         real(dp) :: a1_next, b1_next, c1_next, d1_next, a2_next, b2_next, c2_next, d2_next

         !input from inlists
         real(dp) :: M_ns_initial, omega_initial, e_initial, D, R0, r1, r2
         real(dp) :: op_const, eta, efactor, M_crust, n_poly, beta_sec                      
         real(dp) :: omega_env_factor, prescription  

         !other variables
         real(dp) :: temp, fd, decay_coeff, Req, Rbar, beta, ebind, eorb_change, mdot       
         real(dp) :: v_rel, Ra_br, mdot_br, fd_br, edot_br, mdot_hl_br
         real(dp) :: mdot_mr15_ratio_br, mdot_mr15_br, fd_mr15_ratio_br, fd_mr15_br, fd_hl_br
         integer :: azone

         contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! FUNCTIONS NOT IN USE !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !DONE: To find the closest lower value to a target in an array
      subroutine find_closest_value_index(arr, tar, n, closest_index)
         implicit none
         real(dp), intent(in) :: arr(:)         ! Array to search
         real(dp), intent(in) :: tar            ! Target value
         integer, intent(in) :: n               ! Size of the array
         integer, intent(out) :: closest_index  ! Index of the closest value
   
         integer :: l
         real(dp) :: min_diff, diff, temp
         
         ! print *, 'find closest value index called'
         ! Initialize closest_index and min_diff
         closest_index = 1
         min_diff = abs(arr(1) - tar)
   
         ! Loop through the array to find the closest value
         do l = 2, n
            !  temp = arr(l) - tar
               diff = abs(arr(l) - tar)
               if (diff < min_diff) then
                  min_diff = diff
                  closest_index = l
               end if
         end do
      end subroutine find_closest_value_index

      !DONE: To evaluate the function F for the differential equation da/dt = F(a)
      subroutine linspace(n, from, to, array)
         real(dp), intent(in) :: from, to
         integer, intent(in) :: n
         real(dp), allocatable, intent(out) :: array(:)
         real(dp) :: range
         integer :: i
         range = to - from
   
         if (n == 0) return
         allocate(array(n))
         if (n == 1) then
               array(1) = from
               return
         end if
   
         do i=1, n
               array(i) = from + range * (i - 1) / (n - 1)
         end do
      end subroutine

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! OPERATIONAL FUNCTIONS !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      !forward euler
      subroutine forward_euler(x_n, dxdt_n, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_n, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + dxdt_n*dt_n
      
      end subroutine forward_euler

      !backward euler
      subroutine backward_euler(x_n, dxdt_m, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_m, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + dxdt_m*dt_n
      
      end subroutine backward_euler

      !trapezoidal rule
      subroutine trapezoidal_rule(x_n, dxdt_n, dxdt_m, dt_n, x_m)
         implicit none
         real(dp), intent(in) :: x_n, dxdt_n, dxdt_m, dt_n
         real(dp), intent(out) :: x_m

         x_m = x_n + (dxdt_n + dxdt_m)*dt_n/2.0_dp
      
      end subroutine trapezoidal_rule

      !find zone
      subroutine find_zone(array, target, index)
         implicit none
         real(dp), intent(in) :: array(:)
         real(dp), intent(in) :: target
         integer, intent(out) :: index
         integer :: i
         real(dp) :: diff

         index = 1
         diff = abs(array(1) - target)
         do i = 2, size(array)
            if (diff > abs(array(i) - target)) then
               diff = abs(array(i) - target)
               index = i
            end if
         end do
      end subroutine find_zone

      ! subroutine find_azone(array, target, index)
      !    implicit none
      !    real(dp), intent(in) :: array(:)
      !    real(dp), intent(in) :: target
      !    integer, intent(out) :: index
      !    integer :: i
      !    real(dp) :: diff

      !    index = 1
      !    do i = 1, size(array)
      !       diff = array(i) - target
      !       if (diff .le. 0) then
      !          index = i
      !          return
      !       end if
      !    end do
      ! end subroutine find_azone

      !DONE: To append elements to an allocatable array
      subroutine AddToList(list, element)
         integer :: i, isize
         integer, intent(in) :: element
         integer, dimension(:), allocatable, intent(inout) :: list
         integer, dimension(:), allocatable :: clist

         if(allocated(list)) then
             isize = size(list)
             allocate(clist(isize+1))
             do i=1,isize          
             clist(i) = list(i)
             end do
             clist(isize+1) = element

             deallocate(list)
             call move_alloc(clist, list)

         else
             allocate(list(1))
             list(1) = element
         end if
      end subroutine AddToList
      
      !DONE: To find array indices where array values lie in a given range
      subroutine find_inbetween_value_index(arr, v1, v2, n, ind_btw)

         real(dp), intent(in) :: arr(:)
         real(dp), intent(in) :: v1, v2
         integer, intent(in) :: n
         integer, allocatable, intent(out) :: ind_btw(:)
         
         integer :: l

         ! print *, 'find inbetween value index called'
         do l = 1, n
            if (arr(l) >= v1 .and. arr(l) <= v2) then
               call AddToList(ind_btw, l)
            end if
         end do
         ! print *, '*16*'
      end subroutine find_inbetween_value_index

      !DONE: fourth order runge kutta solver
      !currently, this considers the orbital_evolution_function_holgado as the function to be solved
      subroutine rk4(id, ierr, a_in, t, dt, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t, dt
         real(dp), intent(out) :: a_out
         real(dp) :: k1, k2, k3, k4
         procedure(orbital_evolution_function_holgado), pointer :: G
         
         !getting pointer to star
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'rk4 called'

         !setting the function pointer
         G => orbital_evolution_function_holgado
         
         !evaluating the k values
         call G(id, ierr, t, a_in, k1)
         call G(id, ierr, t + dt/2, a_in + dt*k1/2, k2)
         call G(id, ierr, t + dt/2, a_in + dt*k2/2, k3)
         call G(id, ierr, t + dt, a_in + dt*k3, k4)
         
         !evaluating the next value of a
         a_out = a_in + dt*(k1 + 2*k2 + 2*k3 + k4)/6
      end subroutine rk4
      
      !DONE: To find the radius range and zones for energy injection
      subroutine find_azone_r1r2(array, a_in, r1_out, r2_out, index_out, rzones_out, nr_out)
         implicit none
         real(dp), intent(in) :: array(:)
         real(dp), intent(in) :: a_in
         real(dp), intent(out) :: r1_out, r2_out
         integer, intent(out) :: index_out, nr_out
         integer, allocatable, intent(out) :: rzones_out(:)

         !zone array - to find the zone closest to s%xtra(a_curr)
         call find_zone(array, a_in, index_out)

         !finding the radius range where energy is to be injected

         if (prescription == 1) then
            r1_out = array(index_out) - Ra(index_out)       !cm
            r2_out = array(index_out) + Ra(index_out)       !cm
         else if (prescription == 2) then
            r1_out = array(index_out) - Ra_br       !cm
            r2_out = array(index_out) + Ra_br       !cm
         end if

         call find_inbetween_value_index(array, r1_out, r2_out, size(array), rzones_out)
         nr_out = size(rzones_out)
      end subroutine find_azone_r1r2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! FUNCTIONS THAT EVALUATE ENERGY INJECTION AND ORBITAL EVOLUTION !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !DONE: Holgado et al. 2018 prescription. Calculated da/dt = F(a) and finds the next a using RK4
      subroutine holgado_prescription(id, ierr, M, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         real(dp), intent(in) :: M
         real(dp), intent(out) :: a_out

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !mass is in g, radius is in cm, time is in s

         do k = 1, s%nz
            v(k) = SQRT(standard_cgrav*(M + s%m(k))/s%R(k))   !cm/s
         end do

         call hl_and_energy(id, ierr, M)            !all cgs
         call mr15(id, ierr)                        !all cgs
         call edd_and_hyper(id, ierr, M)            !all cgs

         do k = 1, s%nz
            if (eta*mdot_mr15(k) >= mdot_hyper(k) .or. eta*mdot_mr15(k) <= mdot_edd(k)) then
               mdot_arr(k) = eta*mdot_mr15(k)       !g/s
            else
               mdot_arr(k) = mdot_edd(k)            !g/s
            end if
            ! fd_arr(k) = eta*fd_mr15(k)              !dyne (g cm s^-2)
            call kim2010(id, ierr, v(k)/s%csound(k), s%R(k), Ra(k), k, v(k), M, fd_arr(k))  !all cgs
            edot(k) = fd_arr(k)*v(k)                !erg/s
            f(k) = -fd_arr(k)*v(k)*(1/dEorb(k))     !s^-1
         end do

         call rk4(id, ierr, s%xtra(a_curr), s%time, s%dt_next, a_out)  !cm

      end subroutine holgado_prescription
      
      !DONE: Bronner et al. 2024 prescription. Calculates da/dt over 10 substeps using elliptical orbit evolution subroutine
      subroutine bronner_prescription(id, ierr, M, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k, azone_temp, x
         real(dp), intent(in) :: M
         real(dp), intent(out) :: a_out
         real(dp) :: a_in, f1, f2, f3, mu1, mu2, mu3, mu4, eps_rho_br, mach_br, dt_solver, a_temp 

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092

         mu1 = -2.14034214
         mu2 = 1.94694764
         mu3 = 1.19007536
         mu4 = 1.05762477

         a_in = SQRT(( s%xtra(a1_curr) - s%xtra(c1_curr) )**2 + (s%xtra(b1_curr) - s%xtra(d1_curr))**2)
         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         if (a_in > 100) then
            v_rel = SQRT( (s%xtra(c2_curr) - s%xtra(a2_curr) + (s%xtra(omega_env)*(s%xtra(d1_curr) - s%xtra(b1_curr))) )**2 + (s%xtra(d2_curr) - s%xtra(b2_curr) - (s%xtra(omega_env)*(s%xtra(c1_curr) - s%xtra(a1_curr))) )**2 )
         else 
            v_rel = SQRT( (s%xtra(c2_curr) - s%xtra(a2_curr) + (s%xtra(omega_env)*(s%xtra(d1_curr) - s%xtra(b1_curr))) - (s%u(azone)*(s%xtra(c1_curr) - s%xtra(a1_curr)))/a_in )**2 + (s%xtra(d2_curr) - s%xtra(b2_curr) - (s%xtra(omega_env)*(s%xtra(c1_curr) - s%xtra(a1_curr))) - (s%u(azone)*(s%xtra(d1_curr) - s%xtra(b1_curr)))/a_in )**2 )
         end if

         Ra_br = 2*standard_cgrav*M/(v_rel**2)                 !cm
         mdot_hl_br = pi*(Ra_br**2)*(s%rho(azone_temp))*v_rel           !g/s
         fd_hl_br = mdot_hl_br*v_rel                           !dyne (g cm s^-2)
         eps_rho_br = Ra_br/s%scale_height(azone_temp)                                                !dimensionless
         mdot_mr15_ratio_br = (10)**(mu1 + mu2/(1 + mu3*eps_rho_br + mu4*(eps_rho_br**2)))   !dimensionless
         mdot_mr15_br = mdot_mr15_ratio_br*mdot_hl_br                                        !g/s
         fd_mr15_ratio_br = f1 + f2*eps_rho_br + f3*(eps_rho_br**2)
         fd_mr15_br = fd_mr15_ratio_br*fd_hl_br                                              !dyne (g cm s^-2)

         call edd_and_hyper(id, ierr, M)            !all cgs

         if (eta*mdot_mr15_br >= mdot_hyper(azone_temp) .or. eta*mdot_mr15_br <= mdot_edd(azone_temp)) then
            mdot_br = eta*mdot_mr15_br       !g/s
         else
            mdot_br = mdot_edd(azone_temp)            !g/s
         end if

         mach_br = v_rel / s% csound(azone_temp)

         call kim2010(id, ierr, mach_br, a_in, Ra_br, azone_temp, v_rel, M, fd_br)  !all cgs

         edot_br = fd_br*v_rel                !erg/s

         ! t_orb = 2*pi*sqrt(a_in**3/(standard_cgrav*(s%m(azone_temp) + s%xtra(M_ns_curr))))
         dt_solver = s%dt_next / 10
         a_temp = a_in
         
         do x = 1, 10
            call orbital_evolution_bronner(id, ierr, s%time, dt_solver, a_temp, azone_temp, fd_br, s%xtra(omega_env), M, v_rel, a_out)  !cm
            a_temp = a_out
         end do
      end subroutine bronner_prescription

      !DONE: Bronner et al. 2024 prescription. Calculates next a for elliptical orbits using forward euler method
      subroutine orbital_evolution_bronner(id, ierr, t, dt, a_in, azone_in, fd_in, omega_in, M_in, v_in, a_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(in) :: t, dt, a_in, fd_in, omega_in, M_in, v_in
         integer, intent(in) :: azone_in
         real(dp), intent(out) :: a_out

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !a and b are x and y positions for the primary core and c and d are x and y positions for the NS
         !all variables of type a1_next are positions, type a2_next are velocities
         a1_next = s%xtra(a1_curr) + s%xtra(a2_curr)*dt
         b1_next = s%xtra(b1_curr) + s%xtra(b2_curr)*dt
         c1_next = s%xtra(c1_curr) + s%xtra(c2_curr)*dt
         d1_next = s%xtra(d1_curr) + s%xtra(d2_curr)*dt

         a2_next = a2_curr + dt*(standard_cgrav*s%xtra(M_ns_curr)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3))
         b2_next = b2_curr + dt*(standard_cgrav*s%xtra(M_ns_curr)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3))
         if (a_in > 100*Rsun) then
            c2_next = c2_curr + dt*(standard_cgrav*s%m(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3) + fd_in*(s%xtra(c2_curr) - s%xtra(a2_curr) + omega_in*(s%xtra(d1_curr) - s%xtra(b1_curr)))/(M_in*v_in))
            d2_next = d2_curr + dt*(standard_cgrav*s%m(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3) + fd_in*(s%xtra(d2_curr) - s%xtra(b2_curr) - omega_in*(s%xtra(c1_curr) - s%xtra(a1_curr)))/(M_in*v_in))
         else
            c2_next = c2_curr + dt*(standard_cgrav*s%m(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr))/(a_in**3) + fd_in*(s%xtra(c2_curr) - s%xtra(a2_curr) + omega_in*(s%xtra(d1_curr) - s%xtra(b1_curr)) - (s%u(azone_in)*(s%xtra(c1_curr) - s%xtra(a1_curr)))/a_in)/(M_in*v_in))
            d2_next = d2_curr + dt*(standard_cgrav*s%m(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr))/(a_in**3) + fd_in*(s%xtra(d2_curr) - s%xtra(b2_curr) - omega_in*(s%xtra(c1_curr) - s%xtra(a1_curr)) - (s%u(azone_in)*(s%xtra(d1_curr) - s%xtra(b1_curr)))/a_in)/(M_in*v_in))
         end if

         a_out = SQRT( (a1_next - c1_next)**2 + (b1_next - d1_next)**2 )

      end subroutine orbital_evolution_bronner

      !DONE: Kim & Kim 2010 drag force calculation
      subroutine kim2010(id, ierr, mach_in, a_in, Ra_in, azone_in, v_in, M_in, fd_out)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         real(dp), intent(in) :: M_in, mach_in, a_in, Ra_in, v_in
         integer, intent(in) :: azone_in
         real(dp), intent(out) :: fd_out
         real(dp) :: i_var, beta, eta_b, cd

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         beta = standard_cgrav*M_in/(a_in*s%csound(azone_in)**2)   !dimensionless
         eta_b = beta/(mach_in**2 - 1)                          !dimensionless
         cd = 0.002

         if (mach_in .lt. 1.01) then
            i_var = 0.7706*LOG((1+mach_in)/(1.0004 - 0.9185*mach_in)) - 1.473*mach_in
         else if (mach_in >= 1.01 .and. mach_in .lt. 4.4) then
            i_var = LOG(330*(a_in*(mach_in-0.71)**5.72)/(Ra_in*mach_in**9.58))
         else
            i_var = LOG(a_in / (Ra_in*(0.11*mach_in + 1.65)))
         end if

         if (eta_b > 0.1 .and. mach_in > 1.01) then
            fd_out = -cd*0.7*4*pi*s%rho(azone_in)*(1 + 0.46*(beta**1.1)/(mach_in**2 - 1)**0.11)*(standard_cgrav*(M_in**2))/((v_in**2)*(eta_b**0.5))
         else
            fd_out = -cd*4*pi*s%rho(azone_in)*(standard_cgrav*(M_in**2))*i_var/(v_in**2)
         end if
      end subroutine kim2010

      !DONE: Hoyle-Lyttleton accretion and energy calculation  
      subroutine hl_and_energy(id, ierr, M)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         real(dp), intent(in) :: M

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !mass is in g, radius is in cm, time is in s
         do k = 1, s% nz
            Ra(k) = 2*standard_cgrav*M/(v(k)**2)                 !cm
            mdot_hl(k) = pi*(Ra(k)**2)*(s%rho(k))*v(k)           !g/s
            fd_hl(k) = mdot_hl(k)*v(k)                           !dyne (g cm s^-2)
            edot_hl(k) = fd_hl(k)*v(k)                           !erg/s
            
            Eorb(k) = standard_cgrav*M*(s%m(k))/(2*s%R(k))                                           !erg
            dEorb(k) = (standard_cgrav*M/(2*s%R(k)))*((s%m(k)/s%R(k)) - 4*pi*(s%R(k)**2)*s%rho(k))   !erg/cm
         end do

      end subroutine hl_and_energy

      !DONE: MR15 accretion and drag force calculation
      subroutine mr15(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer :: k
         real(dp) :: f1, f2, f3, mu1, mu2, mu3, mu4

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         f1 = 1.91791946
         f2 = -1.52814698
         f3 = 0.75992092

         mu1 = -2.14034214
         mu2 = 1.94694764
         mu3 = 1.19007536
         mu4 = 1.05762477

         do k = 1, s% nz
            eps_rho(k) = Ra(k)/s%scale_height(k)                                                !dimensionless
            fd_mr15_ratio(k) = f1 + f2*eps_rho(k) + f3*(eps_rho(k)**2)                          !dimensionless
            mdot_mr15_ratio(k) = (10)**(mu1 + mu2/(1 + mu3*eps_rho(k) + mu4*(eps_rho(k)**2)))   !dimensionless
            mdot_mr15(k) = mdot_mr15_ratio(k)*mdot_hl(k)                                        !g/s
            fd_mr15(k) = fd_mr15_ratio(k)*fd_hl(k)                                              !dyne (g cm s^-2)
         end do
      end subroutine mr15

      !DONE: Eddington and hypercritical accretion rate calculation
      subroutine edd_and_hyper(id, ierr, M)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         real(dp), intent(in) :: M

         !subroutine variables
         integer :: k

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !mass is in g, radius is in cm, time is in s
         do k = 1, s% nz
            mdot_edd(k) = 3.5*(1d-8)*(M/(1.33*Msun))*(0.34/op_const)*Msun/secyer      !gm/s
            mdot_hyper(k) = 8.9*(1d-5)*((op_const/0.34)**(-0.73))*Msun/secyer                         !gm/s
         end do

      end subroutine edd_and_hyper

      !DONE: Kim 2007 drag force calculation
      subroutine kim2007(id, ierr, n, M)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k, n
         real(dp) :: mach(n)
         real(dp), intent(in) :: M

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return 

         
         do k = 1, n
            mach(k) = v(k) / s% csound(k)  
            if (mach(k) .lt. 1.0) then
               fd_arr(k) = 0.7706*LOG((1+mach(k))/(1.0004 - 0.9185*mach(k))) - 1.473*mach(k)
            else if (mach(k) >= 1.0 .and. mach(k) .lt. 4.4) then
               fd_arr(k) = LOG(330*(s%R(k)*(mach(k)-0.71)**5.72)/(R0*mach(k)**9.58))
            else
               fd_arr(k) = LOG(s%R(k) / (R0*(0.11*mach(k) + 1.65)))
            end if
            fd_arr(k) = fd_arr(k)*(4*pi*s%rho(k)*((standard_cgrav*M)**2)/(v(k)**2))
         end do
      end subroutine kim2007

      !DONE: to evaluate the function F for the differential equation da/dt = F(a)
      subroutine orbital_evolution_function_holgado(id, ierr, t, a_in, fa) 
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         
         !subroutine variables
         real(dp), intent(in) :: a_in, t
         real(dp), intent(out) :: fa
         integer :: k, azone_temp
         
         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return 
         
         ! print *, 'orbital evolution function called'

         !finding the zone where the function is to be evaluated
         call find_zone(s%R(1:s%nz), a_in, azone_temp)

         !evaluating the function
         fa = f(azone_temp)

      end subroutine orbital_evolution_function_holgado

      !DONE: Bronner et al. 2024 energy injection prescription
      subroutine energy_injection_bronner(id, ierr, edot_in, nr_in, rzones_in, a_in, r1_in, r2_in)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         integer, intent(in) :: nr_in, rzones_in(:)
         real(dp), intent(in) :: a_in, r1_in, r2_in, edot_in
         real(dp) :: temp_heat
         integer :: i, k

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'energy injection bronner called'

         print *, 'edot_in = ', edot_in
         !evaluating the energy injection
         temp_heat = 0
         do i = 1, nr_in
            k = rzones_in(i)
            if (s%R(k) .ge. s%R(rzones_in(nr_in)) .and. s%R(k) .lt. a_in) then 
               s% extra_heat(k)%val = EXP(-((a_in - s%R(k))/(a_in - s%R(rzones_in(nr_in))))**2)
            else if (s%R(k) .ge. a_in .and. s%R(k) .le. s%R(rzones_in(1))) then
               s% extra_heat(k)%val = EXP(-((a_in - s%R(k))/(a_in - s%R(rzones_in(1))))**2)
            else
               s% extra_heat(k)%val = 0
            end if
            temp_heat = temp_heat + s% extra_heat(k)%val*s%dm(k)
         end do

         s% extra_heat(1:s%nz)%val = efactor * s% extra_heat(1:s%nz)%val * edot_in /temp_heat
         print *, "nr_in = ", nr_in
         print *, 'temp_heat = ', temp_heat
         print *, 'extra_heat = ', SUM(s% extra_heat(1:s%nz)%val)
      end subroutine energy_injection_bronner

      !DONE: To find the binding energy of the star up to a given zone
      subroutine binding_energy(id, ierr, zone_in, ebind_out, eorb_change_out)
         !star variables
         integer, intent(in) :: id, zone_in
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp), intent(out) :: ebind_out, eorb_change_out
         integer :: i

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! print *, 'binding energy called'

         ebind_out = 0
         do i = 1, zone_in
            ebind_out = ebind_out + (s%energy(i) - standard_cgrav*s%m(i)/s%R(i))*s%dm(i)
         end do
         eorb_change_out = Eorb(zone_in) - Eorb(1)

      end subroutine binding_energy

      !DONE: To find the equatorial radius and beta secular
      subroutine Req_and_beta(e1, Req1, Rbar1, beta1)
         real(dp), intent(in) :: e1
         real(dp), intent(out) :: Req1, Rbar1, beta1
         
         beta1 = 3*(1 - ((e1*sqrt(1-e1**2))/(asin(e1))))/(2*e1**2) - 1                             !dimensionless
         Rbar1 = R0*((asin(e1) * ((1-e1**2)**(1/6)) * (1-beta1)) / e1)**(-1 * n_poly / (3-n_poly)) !cm
         Req1 = Rbar1/((1-e1**2)**(1/6))                                                           !cm
      end subroutine Req_and_beta

      !DONE: To evaluate the omega function and solve for the inverse
      subroutine omega_function(e1, M, omega1)
         !subroutine variables
         real(dp), intent(in) :: e1, M
         real(dp), intent(out) :: omega1
         real(dp) :: rho_bar, qn

         rho_bar = 3*M/(4*pi*(R0**3))      !gm/cm^3
         qn = (1-n_poly/5)                 !dimensionless
         
         omega1 = sqrt(2*pi*standard_cgrav*rho_bar*( (sqrt(1-e1**2)*(3-2*e1**2)*asin(e1)/(e1**3)) - 3*(1-e1**2)/(e1**2) )/qn) !Hz
      end subroutine omega_function

      !DONE: To solve for the inverse of the omega function
      subroutine omega_func_solve_inverse(e_in, omega_val, tol, M, e_out)
         !subroutine variables
         real(dp), intent(in) :: e_in, omega_val, tol, M
         real(dp), intent(out) :: e_out

         !local variables
         real(dp) :: e1, omega1
         integer :: iterations

         iterations = 0

         e1 = e_in !dimensionless
         call omega_function(e1, M, omega1)

         do while (abs(omega1 - omega_val) >= tol .and. iterations < 100)
            e1 = e1*omega_val/omega1
            call omega_function(e1, M, omega1)
            iterations = iterations + 1
         end do

         e_out = e1 !dimensionless
      end subroutine omega_func_solve_inverse

      !DONE: to evaluate the spin evolution and quadrupole moment evolution
      subroutine omega_and_q(id, ierr, M, omega_function)
         !star variables
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
      end subroutine omega_and_q

      subroutine inject_energy_and_evaluate_strain(id, ierr)
         !star variables
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         !subroutine variables
         real(dp) :: temp, temp_heat
         integer, allocatable :: rzones(:)
         integer :: i, j, nr, azone_temp

         !calling star pointer
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! INITIALISING AND ALLOCATING
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         s%x_ctrl(14) = s%model_number

         !Initialising - These values will be taken from the inlist_project file
         M_ns_initial = s%x_ctrl(1)*Msun          !mass of the NS (gm)
         op_const = s%x_ctrl(2)                   !opacity constant (cm^2/g)
         eta = s%x_ctrl(3)                        !efficiency factor
         efactor = s%x_ctrl(4)                    !multiplication factor for injected energy
         M_crust = s%x_ctrl(5)*Msun               !mass of the crust (gm)
         omega_initial = 2*pi*s%x_ctrl(6)         !initial spin frequency (Hz)
         R0 = s%x_ctrl(7)                         !equatorial radius (cm)
         e_initial = s%x_ctrl(8)                  !initial ellipticity
         n_poly = s%x_ctrl(9)                     !polytropic index
         beta_sec = s%x_ctrl(10)                  !beta secular
         D = s%x_ctrl(11)*1d3*pc                  !distance to the source (cm)
         omega_env_factor = s%x_ctrl(12)          !omega_env_factor*initial_orbital_period = envelope omega
         prescription = s%x_ctrl(13)              !general prescription for the CE evolution

         !allocating relevant variables
         if (allocated(fd_arr)) then !deallocating if already allocated
            deallocate(fd_arr, mdot_arr, f, edot)
            deallocate(mdot_hl, fd_hl, edot_hl, v, Ra, Eorb, dEorb)
            deallocate(mdot_mr15, fd_mr15, mdot_mr15_ratio, fd_mr15_ratio, eps_rho)
            deallocate(mdot_edd, mdot_hyper)!, omega_arr, e_arr, beta_arr)!, rand)
         end if
         allocate(fd_arr(s%nz), mdot_arr(s%nz), f(s%nz), edot(s%nz))
         allocate(mdot_hl(s% nz), fd_hl(s% nz), edot_hl(s% nz), v(s% nz), Ra(s% nz), Eorb(s%nz), dEorb(s%nz))
         allocate(mdot_mr15(s%nz), fd_mr15(s%nz), mdot_mr15_ratio(s%nz), fd_mr15_ratio(s%nz), eps_rho(s%nz))
         allocate(mdot_edd(s%nz), mdot_hyper(s%nz))

         if (s%model_number == 1) then
            s%xtra(a_curr) = 290*Rsun !9.8d-1*s%R(1)          !cm
            s%xtra(M_ns_curr) = M_ns_initial                  !gm
            s%xtra(M_acc_curr) = 0                            !gm   
            s%xtra(omega_curr) = omega_initial                !spin frequency (Hz)
            s%xtra(e_curr) = e_initial                        !ellipticity
            s%xtra(Qmax_curr) = 4*1d39*(s%xtra(M_acc_curr)/M_crust)   !gm cm^2
            s%xtra(aa_curr) = R0*(1+s%xtra(e_curr)/2)                 !cm
            s%xtra(bb_curr) = R0*(1-s%xtra(e_curr)/2)                 !cm
            s%xtra(mom_inert_curr) = s%xtra(M_ns_curr)*(s%xtra(aa_curr)**2 + s%xtra(bb_curr)**2)/5                 !gm cm^2
            call find_zone(s%R(1:s%nz), s%xtra(a_curr), azone_temp)
 
            s%xtra(omega_env) = omega_env_factor*SQRT(standard_cgrav*(s%xtra(M_ns_curr) + s%m(azone_temp))/(s%R(azone_temp)**3))
            if (prescription == 2) then   ! a and b are x and y positions for the primary core and c and d are x and y positions for the NS
                                          ! all variables of type a1_curr are positions, type a2_curr are velocities
               s%xtra(a1_curr) = 0
               s%xtra(b1_curr) = 0
               s%xtra(c1_curr) = s%R(azone_temp)
               s%xtra(d1_curr) = 0
               s%xtra(a2_curr) = 0
               s%xtra(b2_curr) = 0
               s%xtra(c2_curr) = 0
               s%xtra(d2_curr) = SQRT(standard_cgrav*(s%xtra(M_ns_curr) + s%m(azone_temp))/s%R(azone_temp))
            end if
            s%x_ctrl(15) = 2
         end if

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! EVALUATING CURRENT PARAMETERS
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !doing the orbital evolution
         if (prescription == 1) then
            call holgado_prescription(id, ierr, s%xtra(M_ns_curr), a_next)
         else if (prescription == 2) then
            call bronner_prescription(id, ierr, s%xtra(M_ns_curr), a_next)
         end if

         !injecting energy into the envelope
         call find_azone_r1r2(s%R(1:s%nz), s%xtra(a_curr), r1, r2, azone, rzones, nr)        !finding the relevant zones to inject energy into
         
         if (prescription == 1) then
            call energy_injection_bronner(id, ierr, edot(azone), nr, rzones, s%xtra(a_curr), r1, r2)
            mdot = mdot_arr(azone)  
         else if (prescription == 2) then
            call energy_injection_bronner(id, ierr, -edot_br, nr, rzones, s%xtra(a_curr), r1, r2)
            mdot = mdot_br
         end if

         !calculating the binding energy of the envelope(not used in evaluating anything else)
         call binding_energy(id, ierr, azone, ebind, eorb_change)                                                 !ergs
         
         !updating the total accreted mass and NS mass
         if (s%x_ctrl(14) /= s%x_ctrl(15)) then
            call forward_euler(s%xtra(M_acc_curr), mdot, s%dt_next, M_acc_next)     !gm
            call forward_euler(s%xtra(M_ns_curr), mdot, s%dt_next, M_ns_next)       !gm
         else 
            call forward_euler(s%xtra(M_acc_curr), mdot, 0d0, M_acc_next)     !gm
            call forward_euler(s%xtra(M_ns_curr), mdot, 0d0, M_ns_next)       !gm
         end if

         !calculating the spin and quadrupole moment evolution, and also the strain
         call evaluate_strain(id, ierr)

         !printing relevant information
         print *, 'model number = ', s% model_number
         print *, 'orbital separation = ', s%xtra(a_curr)/Rsun
         print *, 'omega = ', s%xtra(omega_curr)
         print *, 'Q = ', s%xtra(Q_curr)
         print *, 'strain = ', s%xtra(strain_curr)
         if (prescription == 1) then
            print *, 'mach = ', v(azone)/s% csound(azone)
            print *, 'edot = ', edot(azone)
         else if (prescription == 2) then
            print *, 'mach = ', v_rel/s% csound(azone)
            print *, 'edot = ', -edot_br
            ! print *, 'fd_br = ', fd_br
            print *, 'a1_curr = ', s%xtra(a1_curr)
            print *, 'b1_curr = ', s%xtra(b1_curr)
            print *, 'c1_curr = ', s%xtra(c1_curr)
            print *, 'd1_curr = ', s%xtra(d1_curr)
            print *, 'a2_curr = ', s%xtra(a2_curr)
            print *, 'b2_curr = ', s%xtra(b2_curr)
            print *, 'c2_curr = ', s%xtra(c2_curr)
            print *, 'd2_curr = ', s%xtra(d2_curr)
            print *, 'total injected energy = ', efactor*SUM(s% extra_heat(1:s%nz)%val * s% dm(1:s%nz))*s%dt_next
         end if

         print *, '##############################################################################'
         print *, 'ONE TIMESTEP DONE'
         print *, '##############################################################################'


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! SAVING RELEVANT VARIABLES
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! save the updated variables
         s%xtra(a_curr) = a_next
         s%xtra(M_acc_curr) = M_acc_next
         s%xtra(M_ns_curr) = M_ns_next
         s%xtra(omega_curr) = omega_next
         s%xtra(e_curr) = e_next
         s%xtra(Qtb_curr) = Qtb_next
         s%xtra(Qmax_curr) = Qmax_next
         s%xtra(Q_curr) = Q_next
         s%xtra(aa_curr) = aa_next
         s%xtra(bb_curr) = bb_next
         s%xtra(mom_inert_curr) = mom_inert_next
         s%xtra(strain_curr) = strain_next

         s%xtra(a1_curr) = a1_next
         s%xtra(b1_curr) = b1_next
         s%xtra(c1_curr) = c1_next
         s%xtra(d1_curr) = d1_next
         s%xtra(a2_curr) = a2_next
         s%xtra(b2_curr) = b2_next
         s%xtra(c2_curr) = c2_next
         s%xtra(d2_curr) = d2_next

         s%x_ctrl(15) = s%model_number

      end subroutine inject_energy_and_evaluate_strain

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

         s% other_energy => inject_energy_and_evaluate_strain
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
      
      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

      end subroutine extras_startup
      

      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_check_model = keep_going         
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


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         if (prescription == 1) then
            how_many_extra_history_columns = 44
         else
            how_many_extra_history_columns = 44
         end if
      end function how_many_extra_history_columns
      
      
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

      
      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         if (prescription == 1) then
            how_many_extra_profile_columns = 14
         else if (prescription == 2) then
            how_many_extra_profile_columns = 0
         end if
      end function how_many_extra_profile_columns
      
      
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


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


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


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


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
      
      
      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve
      
   
         end module run_star_extras
         
   
