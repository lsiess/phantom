!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles pulsating AGB stars
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres   : *number of boundary spheres (integer)*
!   - n_profile_points    : *number of points in stellar profile calculation (integer)*
!   - n_shells            : *total number of shells (integer, if <0 determined automatically from n_particles)*
!   - n_particles_first   : *particles on first shell (0=disabled, >0 builds shells until <100 particles/shell)*
!   - rho_power           : *density profile exponent: rho ~ r^(-rho_power)*
!   - r_min_on_rstar      : *inner radius as fraction of R_star*
!   - r_max_on_rstar      : *outer radius as fraction of R_star*
!   - pulsation_period    : *pulsation period (days)*
!   - piston_velocity     : *piston velocity amplitude (km/s)*
!   - rho_inner           : *density at inner boundary r_min (cgs)*
!   - wind_type           : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_timestep  : *pulsation timestep as fraction of pulsation period*
!   - phi0                : *initial phase offset (radians)*
!   - wind_shell_spacing  : *fraction of tangential and radial distance between particles*
!   - reinject_enabled    : *enable dynamic reinjection (logical)*
!   - reinject_period_days: *period between reinjections in days*
!   - mass_loss_start     : *start time for mass-loss calculation in years*
!   - mass_loss_end       : *end time for mass-loss calculation in years*
!   - check_radius_au     : *radius within which to count mass (AU)*
!   - meas_int_days       : *interval for mass measurements in days*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units, set_star
!
 use io, only:fatal
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 private

 integer :: iboundary_spheres     = 5
 integer :: n_profile_points      = 10000
 integer :: n_shells              = 15
 integer :: n_particles_first     = 0
 real    :: min_particles_shell   = 100.
 real    :: rho_power             = 4.0
 real    :: r_min_on_rstar        = 0.9
 real    :: r_max_on_rstar        = 1.4
 real    :: dtpulsation           = huge(0.)
 real    :: pulsation_period_days = 300.0
 real    :: piston_velocity_km_s  = 4.0
 real    :: time_puls             = -1.0
 real    :: rho_inner             = 1.0e-12
 integer :: wind_type             = 1
 real    :: pulsation_timestep    = 0.02
 real    :: phi0                  = -3.1415926536d0/2.0
 real    :: wind_shell_spacing    = 1.0
 integer :: var_boundary          = 0
 integer :: save_period           = 0
 integer :: dumps_p_period        = 10

 integer :: reinject_enabled      = 1
 real    :: reinject_period_d     = 0.025
 real    :: meas_int_d            = 0.025
 real    :: mass_loss_start       = 1.0
 real    :: mass_loss_end         = 3.0
 real    :: check_radius_au       = 3.0
 integer :: update_L              = 0
 integer :: verbose               = 1

 integer, parameter :: wind_emitting_sink = 1
 integer, parameter :: max_measurements   = 10000

 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_max
 real :: Mtotal, Msink

 integer, allocatable :: npart_per_shell(:)
 integer, allocatable :: npart_per_boundary_shell(:)

 real, allocatable :: delta_r_gas(:)
 real, allocatable :: delta_r_boundary(:)
 real, allocatable :: shell_radii_gas(:)
 real, allocatable :: shell_radii_bnd(:)

 real :: mass_of_gas_particle      = 0.0
 real :: mass_of_boundary_particle = 0.0

 real, allocatable :: delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 integer :: n_shells_total
 integer :: n_shells_bnd

 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles
 integer              :: active_boundary_spheres

 logical :: reinjection_needed = .false.

 real    :: time_last_reinject = 0.0
 real    :: reinject_period
 integer :: n_reinjections     = 0

 real    :: mass_loss_start_time
 real    :: mass_loss_end_time
 real    :: mass_loss_check_radius
 real    :: measurement_interval
 real    :: time_next_measurement
 real    :: mass_previous_measurement
 real, allocatable :: mass_loss_rates(:)
 integer :: n_measurements              = 0
 real    :: mean_mass_loss_rate         = 0.0
 logical :: mass_loss_rate_calculated   = .false.
 logical :: measurement_active          = .false.
 integer :: particles_to_inject         = 0


 character(len=*), parameter :: label = 'inject_atmosphere'

contains

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres     = 5
 n_profile_points      = 10000
 n_shells              = 15
 n_particles_first     = 0
 min_particles_shell   = 100.
 rho_power             = 4.0
 r_min_on_rstar        = 0.9
 r_max_on_rstar        = 1.4
 dtpulsation           = huge(0.)
 rho_inner             = 1.0e-12
 wind_type             = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s  = 4.0
 time_puls             = -1.0
 pulsation_timestep    = 0.02
 phi0                  = -3.1415926536d0/2.0
 wind_shell_spacing    = 1.0
 var_boundary          = 0
 save_period           = 0
 dumps_p_period        = 10
 reinject_enabled      = 1
 reinject_period_d     = 0.025
 meas_int_d            = 0.025
 mass_loss_start       = 1.0
 mass_loss_end         = 3.0
 check_radius_au       = 3.0
 update_L              = 0
 verbose               = 1

end subroutine set_default_options_inject

!----------------------------------------------------------------
!+
!  Initialize everything
!+
!----------------------------------------------------------------
subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km,years
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,unit_velocity,unit_luminosity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum,npartoftype
 use injectutils,   only:get_neighb_distance
 use wind_pulsating,only:setup_star,calc_stellar_profile,region_mass,interp_stellar_profile
 use dust_formation,only:calc_kappa_max
 use timestep,      only:dtmax

 integer, intent(out) :: ierr
 real    :: Mstar_cgs, Rstar_cgs, Tstar, Lstar_cgs
 real    :: current_radius, dr
 integer :: shell_index, max_shells, n_first, n_shell
 integer :: expected_measurements, i
 integer, parameter  :: max_shells_tmp = 2000
 integer, parameter  :: max_iter_dr    = 100
 real,    parameter  :: tol_dr         = 1.0e-6
 real    :: tmp_dr(max_shells_tmp), tmp_r(max_shells_tmp)
 integer :: tmp_n(max_shells_tmp), int_particles_outer
 logical :: converged, file_exists
 integer :: iunit

 ierr = 0

 if (nptmass < 1) call fatal(label,'need at least one sink particle for central star')

 Mtotal    = xyzmh_ptmass(4, wind_emitting_sink)
 Rstar     = xyzmh_ptmass(iReff, wind_emitting_sink)
 Rstar_cgs = Rstar * au
 Mstar_cgs = Mtotal * solarm
 Tstar     = xyzmh_ptmass(iTeff, wind_emitting_sink)
 Lstar_cgs = xyzmh_ptmass(iLum, wind_emitting_sink) * unit_luminosity

 call calc_kappa_max(Mstar_cgs, Lstar_cgs)

 Msink = Mtotal

 inquire(file='mass_loss_rate.dat', exist=file_exists)

 if ( npartoftype(igas) < 100 .and. file_exists) then
       print *, 'Existing mass loss data file found, but this is a fresh start, so delete'
       open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
       close(iunit, status='delete')
 endif

 if (.not. file_exists) then
    print *,'@1 - dont see the need for redefining xyazhm_ptmass(4,1) ?'
    xyzmh_ptmass(4, wind_emitting_sink) = Msink
 endif

 active_boundary_spheres = iboundary_spheres

 if (wind_type == 2 .and. .not. file_exists) call calculate_period(Mtotal, Rstar, pulsation_period_days)

 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation  = 2.0*pi / pulsation_period
 piston_velocity  = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc       = pulsation_period * piston_velocity / (2.0*pi)

 r_min = r_min_on_rstar * Rstar + deltaR_osc * sin(phi0)
 r_max = r_max_on_rstar * Rstar + deltaR_osc * sin(phi0)
 if (r_min <= 0.) call fatal(label,'r_min must be > 0')

 if (save_period == 1) then
    print *,'@2 - calculation of dtmax should be askked in setup, so remove dumps_p_period'
    dtmax = 1. / (dumps_p_period) * pulsation_period
    print *, 'dtmax: ', dtmax
 endif

 reinject_period        = reinject_period_d * pulsation_period
 measurement_interval   = meas_int_d        * pulsation_period
 mass_loss_start_time   = mass_loss_start   * pulsation_period
 mass_loss_end_time     = mass_loss_end     * pulsation_period
 mass_loss_check_radius = check_radius_au
 time_next_measurement  = mass_loss_start_time
 n_measurements         = 0

 expected_measurements = ceiling((mass_loss_end_time - mass_loss_start_time) / measurement_interval) + 1
 allocate(mass_loss_rates(expected_measurements))
 mass_loss_rates = 0.0

 if (min_particles_shell < 1.0) then
    int_particles_outer = nint(min_particles_shell * n_particles_first)
 else
    int_particles_outer = nint(min_particles_shell)
 endif

 if (n_shells > 0) then
    max_shells = n_shells
 else
    max_shells = 200
 endif
 print *,'@3 - can we get rid of n_shells and n_particles_first ?'

 !+
 ! determine the number of particles per shell and their radial distribution
 !+
 if (n_particles_first > 0) then

    current_radius = r_min
    shell_index    = 0

    do
       shell_index = shell_index + 1
       if (shell_index > max_shells_tmp) then
          print *,'@4 - this cannot be done automatically, needs file edition, can we get rid of it ?'
          call fatal(label,'max_shells_tmp exceeded; increase max_shells_tmp')
       endif

       if (shell_index == 1) then
          n_shell = n_particles_first
       else
          print *,'@5 - fixed an error in the expressionof n_shell'
          n_shell = max(1,nint(real(tmp_n(shell_index-1))*(current_radius/tmp_r(shell_index-1))**(2.*(1.-rho_power/3.)))
       endif

       dr = 0.

       if (shell_index > 1 .and. n_shell < int_particles_outer) then
          shell_index = shell_index - 1
          print *,'@6 - why do you remove dr, it has just been set to zero ?'
          r_max = current_radius - dr
          exit
       endif

       dr = wind_shell_spacing * current_radius * get_neighb_distance(n_shell)

       tmp_dr(shell_index) = dr
       tmp_r(shell_index)  = current_radius
       tmp_n(shell_index)  = n_shell
       current_radius      = current_radius + dr
    enddo

 else
    n_first   = 100
    converged = .false.

    do while (.not. converged)
       n_first        = n_first + 10
       current_radius = r_min
       shell_index    = 0

       do
          shell_index = shell_index + 1
          if (shell_index > max_shells_tmp) &
             call fatal(label,'max_shells_tmp exceeded; increase max_shells_tmp')

          if (shell_index == 1) then
             n_shell = n_first
          else
             n_shell = max(1, nint(real(tmp_n(shell_index-1))* ( current_radius / tmp_r(shell_index-1))**(2.*(3.-rho_power)/3.)))
          endif

          dr = wind_shell_spacing * current_radius * get_neighb_distance(n_shell)

          if (current_radius + dr > r_max) then
             shell_index = shell_index - 1
             exit
          endif

          tmp_dr(shell_index) = dr
          tmp_r(shell_index)  = current_radius
          tmp_n(shell_index)  = n_shell
          current_radius      = current_radius + dr
       enddo

       if (shell_index >= n_shells) converged = .true.
    enddo
 endif

 print *,'@7 - same call ? no need to test n_particles_first ?'
 if (n_particles_first > 0 ) then
    call setup_star(Msink * umass, Tstar, r_max * au, r_min  *au, &
               gmw, gamma, rho_inner, rho_power)
 else
    call setup_star(Msink * umass, Tstar, r_max * au, r_min  *au, &
               gmw, gamma, rho_inner, rho_power)
 endif

 !+
 ! determine hydrostatic profile
 !+
 print *,'@9 - we can easily remove n_profile_points'
 call calc_stellar_profile(n_profile_points)

 mass_of_gas_particle = region_mass(r_min, r_max) / real(sum(tmp_n(1:shell_index)))

 n_shells_total = shell_index
 n_shells_bnd   = min(iboundary_spheres, n_shells_total)

 mass_of_gas_particle      = region_mass(r_min, r_max) / real(sum(tmp_n(1:n_shells_total)))
 mass_of_boundary_particle = mass_of_gas_particle

 allocate(npart_per_boundary_shell(n_shells_bnd))
 allocate(delta_r_boundary(n_shells_bnd))
 allocate(shell_radii_bnd(n_shells_bnd))
 do i = 1, n_shells_bnd
    npart_per_boundary_shell(i) = tmp_n(i)
    delta_r_boundary(i)         = tmp_dr(i)
    shell_radii_bnd(i)          = tmp_r(i)
 enddo

 allocate(npart_per_shell(n_shells_total - n_shells_bnd))
 allocate(delta_r_gas(n_shells_total - n_shells_bnd))
 allocate(shell_radii_gas(n_shells_total - n_shells_bnd))
 do i = 1, n_shells_total - n_shells_bnd
    npart_per_shell(i)  = tmp_n(n_shells_bnd + i)
    delta_r_gas(i)      = tmp_dr(n_shells_bnd + i)
    shell_radii_gas(i)  = tmp_r(n_shells_bnd + i)
 enddo
 n_shells_total = n_shells_total - n_shells_bnd

 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_bnd + n_shells_total))
 if (n_shells_bnd > 0) delta_r_radial(1:n_shells_bnd) = delta_r_boundary
 delta_r_radial(n_shells_bnd+1 : n_shells_bnd+n_shells_total) = delta_r_gas

 massoftype(igas)      = mass_of_gas_particle
 massoftype(iboundary) = mass_of_boundary_particle

 if (file_exists) then
    call read_mass_loss_data()
 endif

 if (verbose == 1) then
    print *, ''
    print *, ' rho_power                        :', rho_power
    print *, ' rho_inner (cgs)                  :', rho_inner
    print *, ' Atmosphere [r_min, r_max] Rstar  :', r_min, r_max
    print *, ' M_atmos / M_total                :', region_mass(r_min, r_max) / Mtotal
    print *, ' M_atmos (Msun)                   :', region_mass(r_min, r_max)
    print *, ' Boundary shells                  :', n_shells_bnd
    print *, ' Gas shells                       :', n_shells_total
    print *, ' Total boundary particles         :', sum(npart_per_boundary_shell)
    print *, ' Total gas particles              :', sum(npart_per_shell)
    print *, ' Innermost boundary N_per_shell   :', npart_per_boundary_shell(1)
    print *, ' Outermost boundary N_per_shell   :', npart_per_boundary_shell(n_shells_bnd)
    print *, ' Innermost gas      N_per_shell   :', npart_per_shell(1)
    print *, ' Outermost gas      N_per_shell   :', npart_per_shell(n_shells_total)
    print *, ' Particle mass (Msun)             :', mass_of_gas_particle
    print *, ''
 endif

 if (verbose == 1) then
    do i = 1, n_shells_bnd
       print *, 'Boundary shell ', i, ': r=', shell_radii_bnd(i)/Rstar, ' Rstar, dr=', delta_r_boundary(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_boundary_shell(i)
    enddo
    do i = 1, n_shells_total
       print *, 'Gas shell      ', i, ': r=', shell_radii_gas(i)/Rstar, ' Rstar, dr=', delta_r_gas(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_shell(i)
    enddo
 endif

end subroutine init_inject

!----------------------------------------------------------------
!+
!  The actual function that is called by phantom
!+
!----------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 use part, only:igas,iboundary,iamtype

 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 dtinject = pulsation_timestep * pulsation_period

 if (npart > 0 .and. .not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
 endif

 if (.not. atmosphere_setup_complete) then
    print *, ''
    print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
    call setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    atmosphere_setup_complete = .true.
    print *, 'Stellar atmosphere setup complete.'
    return
 endif

 if (atmosphere_setup_complete .and. .not. allocated(boundary_particle_ids)) then
    call reconstruct_boundary_info(time,xyzh,npart,xyzmh_ptmass)
    time_last_reinject = time - mod(time, reinject_period)
    call read_mass_loss_data()
 endif

 if (reinject_enabled == 1 .and. .not. mass_loss_rate_calculated) then
    call take_periodic_mass_measurements(time,xyzh,npart,xyzmh_ptmass,npartoftype)
 endif

 if (reinject_enabled == 1 .and. mass_loss_rate_calculated) then
    call check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    if (reinjection_needed) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       time_last_reinject = time
       reinjection_needed = .false.
    endif
 endif

 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!----------------------------------------------------------------
!+
!  Checks how much mass the star has lost, and calculates the mass loss rate
!+
!----------------------------------------------------------------
subroutine take_periodic_mass_measurements(time,xyzh,npart,xyzmh_ptmass,npartoftype)
 use part,   only:igas,iboundary,iphase,iamtype
 use physcon,only:solarm,years,days

 real,    intent(in) :: time
 real,    intent(in) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer, intent(in) :: npartoftype(:)

 real    :: current_mass_within_radius, mass_lost, rate_this_interval
 real    :: sink_mass, x0(3), dx, dy, dz, r
 real    :: sum_rates
 integer :: i

 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active = .true.
    x0        = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4,   wind_emitting_sink)
    mass_previous_measurement = sink_mass
    do i = 1, npart
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r  = sqrt(dx**2 + dy**2 + dz**2)
       if (r <= mass_loss_check_radius) then
          if (iamtype(iphase(i)) == iboundary) then
             mass_previous_measurement = mass_previous_measurement + mass_of_boundary_particle
          else
             mass_previous_measurement = mass_previous_measurement + mass_of_gas_particle
          endif
       endif
    enddo
    time_next_measurement = time + measurement_interval
 endif

 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then
    x0        = xyzmh_ptmass(1:3, wind_emitting_sink)
    sink_mass = xyzmh_ptmass(4,   wind_emitting_sink)
    current_mass_within_radius = sink_mass
    do i = 1, npart
       dx = xyzh(1,i) - x0(1)
       dy = xyzh(2,i) - x0(2)
       dz = xyzh(3,i) - x0(3)
       r  = sqrt(dx**2 + dy**2 + dz**2)
       if (r <= mass_loss_check_radius) then
          if (iamtype(iphase(i)) == iboundary) then
             current_mass_within_radius = current_mass_within_radius + mass_of_boundary_particle
          else
             current_mass_within_radius = current_mass_within_radius + mass_of_gas_particle
          endif
       endif
    enddo
    mass_lost          = mass_previous_measurement - current_mass_within_radius
    rate_this_interval = mass_lost / measurement_interval
    n_measurements     = n_measurements + 1
    mass_loss_rates(n_measurements) = rate_this_interval
    mass_previous_measurement = current_mass_within_radius
    time_next_measurement     = time + measurement_interval
 endif

 if (measurement_active .and. .not. mass_loss_rate_calculated .and. time >= mass_loss_end_time) then
    if (n_measurements > 0) then
       sum_rates = 0.0
       do i = 1, n_measurements
          sum_rates = sum_rates + mass_loss_rates(i)
       enddo
       mean_mass_loss_rate = sum_rates / real(n_measurements)
       particles_to_inject = nint((mean_mass_loss_rate * reinject_period) / mass_of_gas_particle)
       if (particles_to_inject < 1) particles_to_inject = 1
       mass_loss_rate_calculated = .true.
       call write_mass_loss_data()
    endif
 endif

end subroutine take_periodic_mass_measurements

!----------------------------------------------------------------
!+
!  Checks when to reinject
!+
!----------------------------------------------------------------
subroutine check_continuous_reinject(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use units,  only:utime
 use physcon,only:days

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 if ((time - time_last_reinject) < reinject_period .and. time >= mass_loss_start_time) return

 if (verbose == 1) then
    print *, ''
    print *, '-----------------------------------------'
    print *, 'Reinjection triggered at time: ', time * (utime / days), ' days'
    print *, 'Time since last reinject: ', (time - time_last_reinject)*utime/days, ' days'
    print *, '-----------------------------------------'
    print *, ''
 endif

 reinjection_needed = .true.

end subroutine check_continuous_reinject

!----------------------------------------------------------------
!+
!  Inject particles throughout the simulation
!+
!----------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary,iamtype,set_particle_type
 use injectutils,    only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,        only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: old_npart, i
 real    :: r_inject, phase, r_dot, rho, u, T, P
 real    :: x0(3), v0(3)
 real    :: mass_injected

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 r_dot = piston_velocity * cos(phase)

 if (allocated(r_boundary_equilibrium) .and. n_boundary_particles > 0) then
    r_inject = r_boundary_equilibrium(n_boundary_particles) + delta_r_radial(iboundary_spheres + 1)
    r_inject = r_inject + deltaR_osc * sin(phase)
 else
    r_inject = r_min
    do i = 1, iboundary_spheres + 1
       r_inject = r_inject + delta_r_radial(i)
    enddo
    r_inject = r_inject + deltaR_osc * sin(phase)
 endif

 call interp_stellar_profile(r_inject, rho, P, u, T)

 old_npart      = npart
 n_reinjections = n_reinjections + 1

 call inject_geodesic_sphere(n_shells_total + n_reinjections, npart + 1, particles_to_inject, r_inject, r_dot, u, rho, &
                               npart, npartoftype, xyzh, vxyzu, igas, x0, v0, wind_emitting_sink)

 mass_injected = real(npart - old_npart) * mass_of_gas_particle
 xyzmh_ptmass(4, wind_emitting_sink) = xyzmh_ptmass(4, wind_emitting_sink) - mass_injected

 if (verbose == 1) then
    print *, ''
    print *, ' Particles injected         :', (npart - old_npart)
    print *, ' Injection radius           :', r_inject
    print *, ' New total particles        :', npart
    print *, 'Reinjection complete.'
    print *, ''
 endif

end subroutine perform_reinjection

!----------------------------------------------------------------
!+
!  Build the initial atmospheric setup (i.e. build the shells)
!+
!----------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas,iboundary,iphase,iamtype
 use injectutils,    only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,        only:pi,km,au

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i, j, first_particle, nboundary
 real    :: r, r_cur, rho, u, T, P, x0(3), v0(3), v_radial


 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 npart = 0

 r_cur = shell_radii_bnd(1) ! - 0.5*delta_r_boundary(1)
 do i = 1, n_shells_bnd
    r     = shell_radii_bnd(i)
    r_cur = r_cur + delta_r_boundary(i)
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial       = 0.0
    first_particle = npart + 1
    call inject_geodesic_sphere(i, first_particle, npart_per_boundary_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0, wind_emitting_sink)
 enddo

 do i = 1, n_shells_total
    r = shell_radii_gas(i)
    call interp_stellar_profile(r, rho, P, u, T)
    v_radial       = 0.0
    first_particle = npart + 1
    call inject_geodesic_sphere(n_shells_bnd + i, first_particle, npart_per_shell(i), r, v_radial, u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, igas, x0, v0, wind_emitting_sink)
 enddo

 nboundary            = npartoftype(iboundary)
 n_boundary_particles = nboundary

 print *, 'Boundary particles : ', nboundary
 print *, 'Gas particles      : ', npartoftype(igas)

 if (nboundary > 0) then
    allocate(r_boundary_equilibrium(nboundary))
    allocate(boundary_particle_ids(nboundary))

    j = 0
    do i = 1, npart
       if (j >= nboundary) exit
       if (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          r_boundary_equilibrium(j) = sqrt( (xyzh(1,i)-x0(1))**2 + &
                                            (xyzh(2,i)-x0(2))**2 + &
                                            (xyzh(3,i)-x0(3))**2 ) &
                                      - deltaR_osc * sin(phi0)
       endif
    enddo
 endif

end subroutine setup_initial_atmosphere

!----------------------------------------------------------------
!+
!  Reconstructs boundary particle info after resuming from a dump
!+
!----------------------------------------------------------------
subroutine reconstruct_boundary_info(time,xyzh,npart,xyzmh_ptmass)
 use part,   only:iboundary,iphase,iamtype
 use physcon,only:pi

 real,    intent(in) :: time
 real,    intent(inout) :: xyzh(:,:),xyzmh_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: i, j
 real    :: x0(3), r_current, phase

 x0    = xyzmh_ptmass(1:3, wind_emitting_sink)
 phase = omega_pulsation * time + phi0

 n_boundary_particles = 0
 do i = 1, npart
    if (iamtype(iphase(i)) == iboundary) n_boundary_particles = n_boundary_particles + 1
 enddo

 if (n_boundary_particles > 0) then
    allocate(r_boundary_equilibrium(n_boundary_particles))
    allocate(boundary_particle_ids(n_boundary_particles))

    j = 0
    do i = 1, npart
       if (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          r_current = sqrt((xyzh(1,i)-x0(1))**2 + &
                           (xyzh(2,i)-x0(2))**2 + &
                           (xyzh(3,i)-x0(3))**2)
          r_boundary_equilibrium(j) = r_current - deltaR_osc * sin(phase)
       endif
    enddo

    print *, 'Reconstructed boundary particle info:'
    print *, 'Boundary particles: ', n_boundary_particles
 endif

end subroutine reconstruct_boundary_info

!----------------------------------------------------------------
!+
!  Applies the pulsation to the boundary layers
!+
!----------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
  use physcon,        only:pi,solarl
 use wind_pulsating, only:interp_stellar_profile
 use part,           only:iTeff,iLum,iReff

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i, ipart
 real    :: r_eq, r_new, r_current, phase, alpha, deltaR_osc
 real    :: x_hat(3), r_dot, x0(3), v0(3)
 real    :: x, y, z, rho, u, T, P
 real    :: Reff, Teff, Lum

 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 x0 = xyzmh_ptmass(1:3, wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,  wind_emitting_sink)

 phase      = omega_pulsation * time + phi0
 deltaR_osc = pulsation_period * piston_velocity / (2.0 * pi)

 if (time < pulsation_period * time_puls .and. time_puls > 0) then
    alpha = time / (pulsation_period * time_puls)
 else
    alpha = 1.0
 endif

 r_dot = alpha * piston_velocity * cos(phase)

 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    r_eq  = r_boundary_equilibrium(i)

    r_new = r_eq + deltaR_osc * (sin(phi0) * (1.0 - alpha) + alpha * sin(phase))

    x = xyzh(1,ipart) - x0(1)
    y = xyzh(2,ipart) - x0(2)
    z = xyzh(3,ipart) - x0(3)
    r_current = sqrt(x**2 + y**2 + z**2)

    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    xyzh(1,ipart) = r_new * x_hat(1) + x0(1)
    xyzh(2,ipart) = r_new * x_hat(2) + x0(2)
    xyzh(3,ipart) = r_new * x_hat(3) + x0(3)

    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)

    if (var_boundary == 1) then
       call interp_stellar_profile(r_new, rho, P, u, T)
       vxyzu(4,ipart) = u
       xyzh(4,ipart)  = (mass_of_boundary_particle / rho)**(1./3.)
    endif

    if (update_L == 1) then
       Reff = xyzmh_ptmass(iReff,1) + deltaR_osc * (sin(phi0) * (1.0 - alpha) + alpha * sin(phase))
       Teff = xyzmh_ptmass(iTeff,1)
       Lum  = xyzmh_ptmass(iLum,1)
       call get_lum(Lum, Teff, Reff)
       xyzmh_ptmass(iLum,1) = Lum
    endif

 enddo

end subroutine apply_pulsation
!----------------------------------------------------------------
!+
!  Placeholder function
!+
!----------------------------------------------------------------
subroutine update_injected_par

end subroutine update_injected_par

!----------------------------------------------------------------
!+
!  Get luminosity
!+
!----------------------------------------------------------------
subroutine get_lum(Lum,Teff,Reff)
 use physcon, only:au,steboltz,solarl,pi
 use units,   only:unit_luminosity
 real, intent(inout) :: Lum
 real, intent(in)    :: Reff, Teff
 real :: lum_lsun

 lum_lsun = 4.*pi*steboltz*Teff**4*(Reff*au)**2/solarl
 Lum  = lum_lsun*(solarl/unit_luminosity)

end subroutine get_lum

!----------------------------------------------------------------
!+
!  Write mass-loss information to file for resume from dump
!+
!----------------------------------------------------------------
subroutine write_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i

 if (.not. mass_loss_rate_calculated) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'Could not write mass_loss_rate.dat'
    return
 endif

 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) '# Careful, this output is NOT M_sun/yr, but in code units of mass/time!'
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate
 write(iunit,*) Mtotal
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_gas_particle
 write(iunit,*) mass_of_boundary_particle
 do i = 1, n_measurements
    write(iunit,*) mass_loss_rates(i)
 enddo
 close(iunit)

 print *, ''
 write(iprint,*) 'Mass-loss rate data written to mass_loss_rate.dat'
 print *, ' '

end subroutine write_mass_loss_data

!----------------------------------------------------------------
!+
!  Read mass-loss information from file after resuming from dump
!+
!----------------------------------------------------------------
subroutine read_mass_loss_data()
 use io, only:iprint
 integer :: iunit, ierr, i
 logical :: file_exists

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (.not. file_exists) return

 open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)
 read(iunit,*, iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then; close(iunit); return; endif
 read(iunit,*, iostat=ierr) mean_mass_loss_rate
 read(iunit,*, iostat=ierr) Mtotal
 read(iunit,*, iostat=ierr) particles_to_inject
 read(iunit,*, iostat=ierr) n_measurements
 read(iunit,*, iostat=ierr) mass_of_gas_particle
 read(iunit,*, iostat=ierr) mass_of_boundary_particle

 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1, n_measurements
       read(iunit,*, iostat=ierr) mass_loss_rates(i)
       if (ierr /= 0) exit
    enddo
 endif
 close(iunit)

 if (verbose == 1) then
    write(iprint,*) 'Mass-loss rate data read from mass_loss_rate.dat'
    write(iprint,*) ' Mean mass-loss rate          :', mean_mass_loss_rate
    write(iprint,*) ' Gas particle mass            :', mass_of_gas_particle
    write(iprint,*) ' Boundary particle mass       :', mass_of_boundary_particle
    write(iprint,*) ' Particles to inject          :', particles_to_inject
 endif

end subroutine read_mass_loss_data

!----------------------------------------------------------------
!+
!  Use mass-period relation to estimate the pulsation period
!+
!----------------------------------------------------------------
subroutine calculate_period(M, R, pulsation_period_days)
 real, intent(in)  :: M, R
 real, intent(out) :: pulsation_period_days
 real :: logP, logM, logR

 logM = log10(M)
 logR = log10(R * 215.032)
 logP = -1.92 - 0.73*logM + 1.86*logR
 pulsation_period_days = 10.0**logP

 print *, 'Calculated pulsation period (days): ', pulsation_period_days

end subroutine calculate_period

!----------------------------------------------------------------
!+
!  Write options to .in file
!+
!----------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(wind_type,            'wind_type',           'pulsation: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(wind_shell_spacing,   'wind_shell_spacing',  'radial/tangential spacing ratio',iunit)
 call write_inopt(iboundary_spheres,    'iboundary_spheres',   'number of boundary spheres (piston layers)',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period',    'pulsation period (days)',iunit)
 call write_inopt(pulsation_timestep,   'pulsation_timestep',  'pulsation timestep as fraction of period',iunit)
 call write_inopt(piston_velocity_km_s, 'piston_velocity',     'piston velocity amplitude (km/s)',iunit)
 call write_inopt(phi0,                 'phi0',                'initial phase offset (radians)',iunit)
 call write_inopt(rho_inner,            'rho_inner',           'inner boundary density at r_min (cgs)',iunit)
 call write_inopt(rho_power,            'rho_power',           'density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(reinject_enabled,     'reinject_enabled',    'enable dynamic reinjection (0=off, 1=on)',iunit)
 call write_inopt(reinject_period_d,    'reinject_period_d',   'period between reinjections (periods)',iunit)

 call write_inopt(n_profile_points,     'n_profile_points',    'number of points in stellar profile',iunit)
 call write_inopt(n_shells,             'n_shells',            'number of gas shells (if <0 determined from n_particles)',iunit)
 call write_inopt(n_particles_first,    'n_particles_first',   'particles on first shell (0=disabled)',iunit)
 call write_inopt(min_particles_shell,  'min_particles_shell', 'minimum particles per shell when using n_particles_first',iunit)
 call write_inopt(r_min_on_rstar,       'r_min_on_rstar',      'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(r_max_on_rstar,       'r_max_on_rstar',      'gas atmosphere outer radius as fraction of R_star',iunit)
 call write_inopt(time_puls,            'time_puls',           'time for piston to ramp up (in periods, -1=instant)',iunit)
 call write_inopt(var_boundary,         'var_boundary',        'update boundary thermo with pulsation (0=off, 1=on)',iunit)
 call write_inopt(save_period,          'save_period',         'wether to save dumps as fraction of period (0=off, 1=on)',iunit)
 call write_inopt(dumps_p_period,       'dumps_p_period',      'number of dumps per period (if save_period = 1)',iunit)
 call write_inopt(meas_int_d,          'meas_int_d',          'mass measurement interval (periods)',iunit)
 call write_inopt(mass_loss_start,      'mass_loss_start',     'start time for mass-loss calculation (periods)',iunit)
 call write_inopt(mass_loss_end,        'mass_loss_end',       'end time for mass-loss calculation (periods)',iunit)
 call write_inopt(check_radius_au,      'check_radius_au',     'mass-loss counting radius (AU)',iunit)
 call write_inopt(update_L,             'update_L',            'update luminosity with pulsation (0=off, 1=on)',iunit)

end subroutine write_options_inject

!----------------------------------------------------------------
!+
!  Read options from .in file
!+
!----------------------------------------------------------------
subroutine read_options_inject(db,nerr)
 use infile_utils, only:inopts,read_inopt
 use io,           only:warning
 use physcon,      only:pi
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr
 logical, save :: init_opt = .false.
 character(len=*), parameter :: label='read_options'

 if (.not. init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
 call read_inopt(wind_type,'wind_type',db,errcount=nerr,min=1,max=2)
 call read_inopt(wind_shell_spacing,'wind_shell_spacing',db,errcount=nerr,min=epsilon(0.))
 call read_inopt(iboundary_spheres,'iboundary_spheres',db,errcount=nerr,min=0)
! call read_inopt(outer_boundary_au,'outer_boundary',db,errcount=nerr,min=0.)
 !pulsation parameters
 call read_inopt(pulsation_period,'pulsation_period',db,min=0.,errcount=nerr)
 call read_inopt(pulsation_timestep,'pulsation_timestep',db,min=0.,max=1.,errcount=nerr)
 call read_inopt(piston_velocity,'piston_velocity',db,min=0.,errcount=nerr)
 call read_inopt(phi0,'phi0',db,min=-pi,max=pi,errcount=nerr)
 call read_inopt(rho_inner,'rho_inner',db,min=0.,errcount=nerr)
 call read_inopt(rho_power,'rho_power',db,min=0.,errcount=nerr)
 call read_inopt(reinject_enabled,'reinject_enabled',db,min=0,max=1,errcount=nerr)
 call read_inopt(reinject_period_d,'reinject_period_d',db,min=0.,errcount=nerr)
 call read_inopt(n_profile_points,     'n_profile_points',    db,min=11,errcount=nerr)
 call read_inopt(n_shells,             'n_shells',            db,min=0,errcount=nerr)
 call read_inopt(n_particles_first,    'n_particles_first',   db,min=0,errcount=nerr)
 call read_inopt(min_particles_shell,  'min_particles_shell', db,min=0.,errcount=nerr)
 call read_inopt(r_min_on_rstar,       'r_min_on_rstar',      db,min=0.,max=2.,errcount=nerr)
 call read_inopt(r_max_on_rstar,       'r_max_on_rstar',      db,min=0.,errcount=nerr)
 call read_inopt(time_puls,            'time_puls',           db,min=-1.,errcount=nerr)
 call read_inopt(var_boundary,         'var_boundary',        db,min=0,max=1,errcount=nerr)
 call read_inopt(save_period,          'save_period',         db,min=0,max=1,errcount=nerr)
 call read_inopt(dumps_p_period,       'dumps_p_period',      db,min=0,errcount=nerr)
 call read_inopt(meas_int_d,           'meas_int_d',          db,min=0.,errcount=nerr)
 call read_inopt(mass_loss_start,      'mass_loss_start',     db,min=0.,errcount=nerr)
 call read_inopt(mass_loss_end,        'mass_loss_end',       db,min=0.,errcount=nerr)
 call read_inopt(check_radius_au,      'check_radius_au',     db,min=0.,errcount=nerr)
 call read_inopt(update_L,             'update_L',            db,min=0,max=1,errcount=nerr)

end subroutine read_options_inject

end module inject
