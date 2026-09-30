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
!   - iboundary_spheres      : *number of boundary spheres (integer)*
!   - n_profile_points       : *number of points in stellar profile calculation (integer)*
!   - min_mass_fraction      : *criteria determining when to stop building shells, i.e., when M_shell / M_tot < mass_fraction*
!   - rho_power_in           : *density profile exponent: rho ~ r^(-rho_power)*
!   - r_min_on_rstar         : *inner radius as fraction of R_star*
!   - dtpulsation            : *pulsation timestep as fraction of pulsation period*
!   - rho_inner              : *density at inner boundary r_min (cgs)*
!   - iwind                  : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_period_days  : *pulsation period (days) if iwind != 2*
!   - piston_velocity_km_s   : *piston velocity amplitude (km/s)*
!   - phi0                   : *initial phase offset (radians) (best = -pi/2)*
!   - wss                    : *fraction of tangential and radial distance between particles*
!   - save_period            : *wether to save dumps as an multiple of the pulsation period (0=off, 1=on)*
!   - dumps_p_period         : *how many dumps to save every period, if save_period is activated*
!   - reinject_enabled       : *enable reinjection (logical)*
!   - n_inject_period        : *number of reinjections per period*
!   - mass_loss_start        : *start time for mass-loss calculation in pulsation periods*
!   - mass_loss_end          : *end time for mass-loss calculation in pulsation periods*
!   - update_L               : *wether to update the luminosity of the sink particle with the pulsation period (0=off, 1=on)*
!   - use_file_mdot          : *skip measurement phase and use mass_loss_rate.dat directly (0=off, 1=on)*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units
!
 use io,      only:fatal
 use physcon, only:pi
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject, inject_particles, write_options_inject, read_options_inject, &
           set_default_options_inject, update_injected_par
 private

 integer :: iboundary_spheres     = 5
 integer :: n_profile_points      = 10000
 real    :: min_mass_fraction     = 0.01
 real    :: rho_power_in          = 6.0
 real    :: r_min_on_rstar        = 1.0
 real    :: rho_inner             = 1.0e-12
 integer :: iwind                 = 1
 real    :: pulsation_period_days = 300.0
 real    :: piston_velocity_km_s  = 4.0
 real    :: phi0                  = -0.5*pi
 real    :: wss                   = 1.0
 integer :: save_period           = 0
 integer :: dumps_p_period        = 10
 integer :: reinject_enabled      = 1
 integer :: n_inject_period       = 40
 real    :: mass_loss_start       = 6.0
 real    :: mass_loss_end         = 8.0
 integer :: update_L              = 0
 integer :: verbose               = 1
 integer :: use_file_mdot         = 0

 integer,          parameter :: wind_emitting_sink = 1
 real,             parameter :: pulsation_timestep = 0.02
 character(len=*), parameter :: label = 'inject_atmosphere'

 real    :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real    :: r_min, r_max, Mtotal, mass_of_gas_particle
 integer :: n_shells_total, n_shells_bnd
 integer, allocatable :: npart_per_shell(:), npart_per_boundary_shell(:)
 real,    allocatable :: shell_radii_gas(:), shell_radii_bnd(:), delta_r_radial(:)

 logical :: atmosphere_setup_complete = .false.
 logical :: boundary_info_ready       = .false.
 integer :: n_boundary_particles      = 0
 real, allocatable :: r_boundary_equilibrium(:)

 real    :: reinject_period
 real    :: time_last_reinject  = 0.
 integer :: n_reinjections      = 0
 integer :: particles_to_inject = 0

 real    :: mass_loss_start_time, mass_loss_end_time, time_next_measurement
 real    :: mean_mass_loss_rate       = 0.
 integer :: n_measurements            = 0
 integer :: n_escaping_prev           = 0
 logical :: mass_loss_rate_calculated = .false.
 logical :: measurement_active        = .false.
 real, allocatable :: mass_loss_rates(:)

contains

subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres     = 5
 n_profile_points      = 10000
 min_mass_fraction     = 0.01
 rho_power_in          = 6.0
 r_min_on_rstar        = 1.0
 rho_inner             = 1.0e-12
 iwind                 = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s  = 4.0
 phi0                  = -0.5*pi
 wss                   = 1.0
 save_period           = 0
 dumps_p_period        = 10
 reinject_enabled      = 1
 n_inject_period       = 40
 mass_loss_start       = 6.0
 mass_loss_end         = 8.0
 update_L              = 0
 verbose               = 1
 use_file_mdot         = 0

end subroutine set_default_options_inject

!----------------------------------------------------------------
!+
!  Initialize everything
!+
!----------------------------------------------------------------
subroutine init_inject(ierr)
 use physcon,        only:days,au,solarm,km
 use eos,            only:gmw,gamma
 use units,          only:utime,umass,udist,unit_velocity,unit_luminosity
 use part,           only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff,iLum,npartoftype
 use injectutils,    only:get_fibonacci_spacing,find_optimal_rotation
 use wind_pulsating, only:setup_star,calc_stellar_profile,region_mass
 use dust_formation, only:calc_kappa_max
 use timestep,       only:dtmax
 integer, intent(out) :: ierr
 integer, parameter :: max_shells_tmp = 2000
 real    :: tmp_dr(max_shells_tmp), tmp_r(max_shells_tmp)
 integer :: tmp_n(max_shells_tmp)
 real    :: Mstar, Rstar, Tstar, current_radius, dr
 integer :: nshells, n_shell, n_tot, i, iunit
 integer :: n_particles_first              ! derived locally — no longer a module parameter
 logical :: file_exists

 ierr = 0
 if (nptmass < 1) call fatal(label,'need at least one sink particle for central star')

 Mstar = xyzmh_ptmass(4,wind_emitting_sink)
 Rstar = xyzmh_ptmass(iReff,wind_emitting_sink)
 Tstar = xyzmh_ptmass(iTeff,wind_emitting_sink)
 call calc_kappa_max(Mstar*solarm, xyzmh_ptmass(iLum,wind_emitting_sink)*unit_luminosity)

 inquire(file='mass_loss_rate.dat', exist=file_exists)
 if (npartoftype(igas) < 100 .and. file_exists .and. use_file_mdot == 0) then
    print*,'Existing mass loss data file found, but this is a fresh start, so delete'
    open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
    close(iunit, status='delete')
    file_exists = .false.
 endif

 if (iwind == 2 .and. .not. file_exists) call calculate_period(Mstar, Rstar, pulsation_period_days)

 pulsation_period = pulsation_period_days*days/utime
 omega_pulsation  = 2.*pi/pulsation_period
 piston_velocity  = piston_velocity_km_s*km/unit_velocity
 deltaR_osc       = piston_velocity/omega_pulsation

 r_min = r_min_on_rstar*Rstar + deltaR_osc*sin(phi0)
 if (r_min <= 0.) call fatal(label,'r_min must be > 0')

 if (save_period == 1) then
    dtmax = pulsation_period/dumps_p_period
    print*,'dtmax: ',dtmax
 endif

 mass_of_gas_particle  = massoftype(igas)          ! code units (Msun when G=1, dist=au)
 massoftype(iboundary) = mass_of_gas_particle

 call derive_n_particles_first(r_min*udist, rho_inner, rho_power_in, wss, &
                               mass_of_gas_particle*umass, n_particles_first)

 current_radius = r_min
 nshells = 0
 n_tot   = 0
 do
    nshells = nshells + 1
    if (nshells > max_shells_tmp) call fatal(label,'max_shells_tmp exceeded; increase max_shells_tmp')
    if (nshells == 1) then
       n_shell = n_particles_first
    else
       n_shell = max(1, nint(tmp_n(nshells-1)*(current_radius/tmp_r(nshells-1))**(2.*(3.-rho_power_in)/3.)))
    endif
    dr = wss*current_radius*get_fibonacci_spacing(n_shell)
    if (nshells > 2 .and. real(n_shell)/real(n_tot) < min_mass_fraction) then
       nshells = nshells - 1
       r_max   = current_radius - dr
       exit
    endif
    tmp_dr(nshells) = dr
    tmp_r(nshells)  = current_radius
    tmp_n(nshells)  = n_shell
    current_radius  = current_radius + dr
    n_tot           = n_tot + n_shell
 enddo

 call setup_star(Mstar*umass, Tstar, r_max*au, r_min*au, gmw, gamma, rho_inner, rho_power_in)
 call calc_stellar_profile(n_profile_points)

 n_shells_bnd   = min(iboundary_spheres, nshells)
 n_shells_total = nshells - n_shells_bnd
 allocate(npart_per_boundary_shell(n_shells_bnd), shell_radii_bnd(n_shells_bnd))
 allocate(npart_per_shell(n_shells_total), shell_radii_gas(n_shells_total))
 allocate(delta_r_radial(nshells))
 npart_per_boundary_shell = tmp_n(1:n_shells_bnd)
 shell_radii_bnd          = tmp_r(1:n_shells_bnd)
 npart_per_shell          = tmp_n(n_shells_bnd+1:nshells)
 shell_radii_gas          = tmp_r(n_shells_bnd+1:nshells)
 delta_r_radial           = tmp_dr(1:nshells)

 Mtotal = Mstar + sum(tmp_n(1:n_shells_total))*mass_of_gas_particle

 reinject_period      = pulsation_period/n_inject_period
 mass_loss_start_time = mass_loss_start*pulsation_period
 mass_loss_end_time   = mass_loss_end*pulsation_period
 allocate(mass_loss_rates(ceiling((mass_loss_end_time - mass_loss_start_time)/reinject_period) + 1))
 mass_loss_rates = 0.

 if (verbose == 1) then
   print *, ''
   print *, 'Calculated reinject period:', reinject_period
   print *, 'Measurement period:', reinject_period
   print *, 'Mass loss measurement start time:', mass_loss_start_time
   print *, 'Mass loss measurement end time  :', mass_loss_end_time
   print *, 'Rmax                            :', r_max
   print *, 'Expected number of measurements :', size(mass_loss_rates)
   print *, ''
 endif

 ! If use_file_mdot=1, read mass_loss_rate.dat now and mark measurement as done,
 ! so the simulation goes straight to reinjection without any measurement phase.
 ! A hard error is raised if the file is absent, since the user explicitly asked for it.
 if (use_file_mdot == 1 .and. .not. file_exists) &
    call fatal(label,'use_file_mdot=1 but mass_loss_rate.dat not found')
 if (file_exists) then
    call read_mass_loss_data()
    if (use_file_mdot == 1) mass_loss_rate_calculated = .true.
    if (mass_loss_rate_calculated) call find_optimal_rotation(particles_to_inject)
    if (use_file_mdot == 1 .and. verbose == 1) then
       print *, ''
       print *, 'use_file_mdot=1: skipping measurement phase.'
       print *, 'particles_to_inject read from file:', particles_to_inject
       print *, ''
    endif
 endif

 if (verbose == 1) then
    print *, ''
    print *, ' rho_power                        :', rho_power_in
    print *, ' rho_inner (cgs)                  :', rho_inner
    print *, ' Particle mass from setup (Msun)  :', massoftype(igas)
    print *, ' n_particles_first (derived)       :', n_particles_first
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
       print *, 'Boundary shell ', i, ': r=', shell_radii_bnd(i)/Rstar, ' Rstar, dr=', delta_r_radial(i)/Rstar, &
                ' Rstar, N_particles=', npart_per_boundary_shell(i)
    enddo
    do i = 1, n_shells_total
       print *, 'Gas shell      ', i, ': r=', shell_radii_gas(i)/Rstar, ' Rstar, dr=', delta_r_radial(n_shells_bnd+i)/Rstar, &
                ' Rstar, N_particles=', npart_per_shell(i)
    enddo
 endif

end subroutine init_inject

!----------------------------------------------------------------
!+
!  Derive n_particles_first from the target particle mass and the
!  power-law density profile at r_min.
!+
!----------------------------------------------------------------
subroutine derive_n_particles_first(r_min_cgs,rho_inner_cgs,rho_power,wss_in,m_particle_cgs,n_first)
 use injectutils, only:get_fibonacci_spacing
 real,    intent(in)  :: r_min_cgs,rho_inner_cgs,rho_power,wss_in,m_particle_cgs
 integer, intent(out) :: n_first
 integer, parameter :: max_iter = 100
 integer :: iter,n_old
 real    :: C_rho,exponent,dr_cgs

 C_rho    = rho_inner_cgs*r_min_cgs**rho_power
 exponent = 3. - rho_power
 n_old    = 1000
 do iter = 1,max_iter
    dr_cgs  = wss_in*r_min_cgs*get_fibonacci_spacing(n_old)
    n_first = max(1, nint(4.*pi*C_rho*((r_min_cgs + dr_cgs)**exponent - r_min_cgs**exponent) &
                          /exponent/m_particle_cgs))
    if (n_first == n_old) exit
    n_old = (n_old + n_first)/2
 enddo

 if (verbose == 1) then
    print *, ''
    print *, ' derive_n_particles_first: converged in ', iter, ' iterations'
    print *, ' Target particle mass (cgs)    :', m_particle_cgs
    print *, ' First-shell particle count    :', n_first
    print *, ' First-shell dr / r_min        :', dr_cgs / r_min_cgs
    print *, ''
 endif

end subroutine derive_n_particles_first

!----------------------------------------------------------------
!+
!  The actual function that is called by phantom
!+
!----------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 dtinject = pulsation_timestep*pulsation_period

 if (.not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
    if (npart == 0) then
       print *, ''
       print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
       call setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       print *, 'Stellar atmosphere setup complete.'
       return
    endif
 endif

 if (.not. boundary_info_ready) then
    call reconstruct_boundary_info(time,xyzh,xyzmh_ptmass)
    time_last_reinject = time - mod(time,reinject_period)
 endif

 if (reinject_enabled == 1) then
    if (.not. mass_loss_rate_calculated) &
       call take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    if (mass_loss_rate_calculated .and. time - time_last_reinject >= reinject_period &
        .and. time >= mass_loss_start_time) then
       call perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
       time_last_reinject = time
    endif
 endif

 call apply_pulsation(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)

end subroutine inject_particles

!----------------------------------------------------------------
!+
!  Counts particles that are unbound from the AGB star
!+
!----------------------------------------------------------------
function count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass) result(n)
 real,    intent(in) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: n,i
 real    :: x0(3),v0(3),m,r

 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 m  = xyzmh_ptmass(4,wind_emitting_sink)
 n  = 0
 do i = 1,npart
    r = norm2(xyzh(1:3,i) - x0)
    if (r <= 0.) cycle
    if (0.5*sum((vxyzu(1:3,i) - v0)**2) + vxyzu(4,i) - m/r > 0.) n = n + 1
 enddo

end function count_unbound

!----------------------------------------------------------------
!+
!  Checks how much mass the star has lost, and calculates the mass loss rate
!+
!----------------------------------------------------------------
subroutine take_periodic_mass_measurements(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use injectutils, only:find_optimal_rotation
 real,    intent(in) :: time,xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in) :: npart
 integer :: n_escaping,newly_unbound

 if (.not. measurement_active .and. time >= mass_loss_start_time) then
    measurement_active    = .true.
    time_next_measurement = time + reinject_period
    n_escaping_prev       = count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    if (verbose == 1) print*,' Baseline unbound particles at measurement start:',n_escaping_prev
 endif

 if (measurement_active .and. time >= time_next_measurement .and. time < mass_loss_end_time) then
    n_escaping      = count_unbound(xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
    newly_unbound   = max(0, n_escaping - n_escaping_prev)
    n_escaping_prev = n_escaping
    n_measurements  = n_measurements + 1
    mass_loss_rates(n_measurements) = newly_unbound*mass_of_gas_particle/reinject_period
    time_next_measurement = time + reinject_period
    if (verbose == 1) then
       print *, ' Unbound+outflowing particles     :',n_escaping
       print*,' Newly unbound since last snapshot:',newly_unbound
       print*,' Rate this interval               :',mass_loss_rates(n_measurements)
    endif
 endif

 if (measurement_active .and. time >= mass_loss_end_time .and. n_measurements > 0) then
    mean_mass_loss_rate       = sum(mass_loss_rates(1:n_measurements))/n_measurements
    particles_to_inject       = max(1, nint(mean_mass_loss_rate*reinject_period/mass_of_gas_particle))
    mass_loss_rate_calculated = .true.
    call write_mass_loss_data()
    call find_optimal_rotation(particles_to_inject)
 endif

end subroutine take_periodic_mass_measurements

!----------------------------------------------------------------
!+
!  Inject particles throughout the simulation
!+
!----------------------------------------------------------------
subroutine perform_reinjection(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,           only:igas
 use injectutils,    only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: old_npart
 real    :: r_inject,phase,rho,u,T,P,x0(3),v0(3)

 x0    = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0    = vxyz_ptmass(1:3,wind_emitting_sink)
 phase = omega_pulsation*time + phi0

 if (n_boundary_particles > 0) then
    r_inject = r_boundary_equilibrium(n_boundary_particles) + delta_r_radial(iboundary_spheres+1)
 else
    r_inject = r_min + sum(delta_r_radial(1:iboundary_spheres+1))
 endif
 r_inject = r_inject + deltaR_osc*sin(phase)
 call interp_stellar_profile(r_inject,rho,P,u,T)

 old_npart      = npart
 n_reinjections = n_reinjections + 1
 call inject_fibonacci_sphere(n_shells_total + n_reinjections, npart+1, particles_to_inject, r_inject, &
                              piston_velocity*cos(phase), u, rho, npart, npartoftype, xyzh, vxyzu, igas, x0, v0)

 xyzmh_ptmass(4,wind_emitting_sink) = xyzmh_ptmass(4,wind_emitting_sink) - (npart - old_npart)*mass_of_gas_particle

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
 use part,           only:igas,iboundary
 use injectutils,    only:inject_fibonacci_sphere
 use wind_pulsating, only:interp_stellar_profile
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npartoftype(:)
 integer :: i
 real    :: rho,u,T,P,x0(3),v0(3)

 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)

 npart = 0
 do i = 1,n_shells_bnd
    call interp_stellar_profile(shell_radii_bnd(i),rho,P,u,T)
    call inject_fibonacci_sphere(i, npart+1, npart_per_boundary_shell(i), shell_radii_bnd(i), 0., u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, iboundary, x0, v0)
 enddo
 do i = 1,n_shells_total
    call interp_stellar_profile(shell_radii_gas(i),rho,P,u,T)
    call inject_fibonacci_sphere(n_shells_bnd+i, npart+1, npart_per_shell(i), shell_radii_gas(i), 0., u, rho, &
                                 npart, npartoftype, xyzh, vxyzu, igas, x0, v0)
 enddo

 n_boundary_particles = npartoftype(iboundary)
 allocate(r_boundary_equilibrium(n_boundary_particles))
 do i = 1,n_boundary_particles
    r_boundary_equilibrium(i) = norm2(xyzh(1:3,i) - x0) - deltaR_osc*sin(phi0)
 enddo
 boundary_info_ready = .true.

 print*,'Boundary particles : ',n_boundary_particles
 print*,'Gas particles      : ',npartoftype(igas)

end subroutine setup_initial_atmosphere

!----------------------------------------------------------------
!+
!  Reconstructs boundary particle info after resuming from a dump
!+
!----------------------------------------------------------------
subroutine reconstruct_boundary_info(time,xyzh,xyzmh_ptmass)
 use part, only:iboundary,npartoftype
 real, intent(in) :: time,xyzh(:,:),xyzmh_ptmass(:,:)
 integer :: i
 real    :: x0(3),osc

 x0  = xyzmh_ptmass(1:3,wind_emitting_sink)
 osc = deltaR_osc*sin(omega_pulsation*time + phi0)

 n_boundary_particles = npartoftype(iboundary)
 allocate(r_boundary_equilibrium(n_boundary_particles))
 do i = 1,n_boundary_particles
    r_boundary_equilibrium(i) = norm2(xyzh(1:3,i) - x0) - osc
 enddo
 boundary_info_ready = .true.

 if (n_boundary_particles > 0) then
    print *, 'Reconstructed boundary particle info:'
    print *, 'Boundary particles: ', n_boundary_particles
 endif

end subroutine reconstruct_boundary_info

!----------------------------------------------------------------
!+
!  Applies the pulsation to the boundary layers
!+
!----------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass)
 use part, only:iTeff,iLum,iReff
 real, intent(in)    :: time,vxyz_ptmass(:,:)
 real, intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:)
 integer :: i
 real    :: phase,osc,r_dot,x0(3),v0(3),x_hat(3)

 if (n_boundary_particles == 0) return

 x0    = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0    = vxyz_ptmass(1:3,wind_emitting_sink)
 phase = omega_pulsation*time + phi0
 osc   = deltaR_osc*sin(phase)
 r_dot = piston_velocity*cos(phase)

 do i = 1,n_boundary_particles
    ! Radial unit vector from current position
    x_hat        = xyzh(1:3,i) - x0
    x_hat        = x_hat/norm2(x_hat)
    xyzh(1:3,i)  = x0 + (r_boundary_equilibrium(i) + osc)*x_hat
    vxyzu(1:3,i) = v0 + r_dot*x_hat
 enddo

 if (update_L == 1) call get_lum(xyzmh_ptmass(iLum,wind_emitting_sink), xyzmh_ptmass(iTeff,wind_emitting_sink), &
                                 xyzmh_ptmass(iReff,wind_emitting_sink) + osc)

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
 use physcon, only:au,steboltz
 use units,   only:unit_luminosity
 real, intent(out) :: Lum
 real, intent(in)  :: Teff,Reff

 Lum = 4.*pi*steboltz*Teff**4*(Reff*au)**2/unit_luminosity

end subroutine get_lum

!----------------------------------------------------------------
!+
!  Write mass-loss information to file for resume from dump
!+
!----------------------------------------------------------------
subroutine write_mass_loss_data()
 use io,      only:iprint
 use physcon, only:solarm,years
 use units,   only:umass,utime
 integer :: iunit,ierr,i

 open(newunit=iunit, file='mass_loss_rate.dat', status='replace', iostat=ierr)
 if (ierr /= 0) then
    write(iprint,*) 'Could not write mass_loss_rate.dat'
    return
 endif

 write(iunit,*) '# Mass-loss rate data for restart'
 write(iunit,*) r_max
 write(iunit,*) r_max
 write(iunit,*) mass_loss_rate_calculated
 write(iunit,*) mean_mass_loss_rate/(solarm/umass)/(utime/years)
 write(iunit,*) Mtotal
 write(iunit,*) particles_to_inject
 write(iunit,*) n_measurements
 write(iunit,*) mass_of_gas_particle
 do i = 1,n_measurements
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
 integer :: iunit,ierr,i

 open(newunit=iunit, file='mass_loss_rate.dat', status='old', iostat=ierr)
 if (ierr /= 0) return

 read(iunit,*)
 read(iunit,*,iostat=ierr) r_max
 read(iunit,*,iostat=ierr) r_max
 read(iunit,*,iostat=ierr) mass_loss_rate_calculated
 if (ierr /= 0) then
    close(iunit)
    return
 endif
 read(iunit,*,iostat=ierr) mean_mass_loss_rate
 read(iunit,*,iostat=ierr) Mtotal
 read(iunit,*,iostat=ierr) particles_to_inject
 read(iunit,*,iostat=ierr) n_measurements
 read(iunit,*,iostat=ierr) mass_of_gas_particle

 if (n_measurements > 0) then
    if (.not. allocated(mass_loss_rates)) allocate(mass_loss_rates(n_measurements))
    do i = 1,n_measurements
       read(iunit,*,iostat=ierr) mass_loss_rates(i)
       if (ierr /= 0) exit
    enddo
 endif
 close(iunit)

 if (verbose == 1) then
    write(iprint,*) 'Mass-loss rate data read from mass_loss_rate.dat'
    write(iprint,*) ' Mean mass-loss rate          :', mean_mass_loss_rate
    write(iprint,*) ' Gas particle mass            :', mass_of_gas_particle
    write(iprint,*) ' Boundary particle mass       :', mass_of_gas_particle
    write(iprint,*) ' Particles to inject          :', particles_to_inject
 endif

end subroutine read_mass_loss_data

!----------------------------------------------------------------
!+
!  Use mass-period relation to estimate the pulsation period
!+
!----------------------------------------------------------------
subroutine calculate_period(M,R,pulsation_period_days)
 real, intent(in)  :: M,R
 real, intent(out) :: pulsation_period_days

 pulsation_period_days = 10.**(-1.92 - 0.73*log10(M) + 1.86*log10(R*215.032))
 print*,'Calculated pulsation period (days): ',pulsation_period_days

end subroutine calculate_period

!----------------------------------------------------------------
!+
!  Write options to .in file
!+
!----------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,      'n_profile_points',  'number of points in stellar profile',iunit)
 call write_inopt(iboundary_spheres,     'iboundary_spheres', 'number of boundary spheres (piston layers)',iunit)
 call write_inopt(min_mass_fraction,     'min_mass_fraction', 'minimum mass fraction per shell',iunit)
 call write_inopt(rho_power_in,          'rho_power',         'density profile exponent: rho ~ r^(-rho_power)',iunit)
 call write_inopt(r_min_on_rstar,        'r_min_on_rstar',    'gas atmosphere inner radius as fraction of R_star',iunit)
 call write_inopt(rho_inner,             'rho_inner',         'inner boundary density at r_min (cgs)',iunit)
 call write_inopt(iwind,                 'iwind',             'wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days, 'pulsation_period',  'pulsation period (days)',iunit)
 call write_inopt(piston_velocity_km_s,  'piston_velocity',   'piston velocity amplitude (km/s)',iunit)
 call write_inopt(phi0,                  'phi0',              'initial phase offset (radians)',iunit)
 call write_inopt(wss,                   'wss',               'radial/tangential spacing ratio',iunit)
 call write_inopt(save_period,           'save_period',       'wether to save dumps as fraction of period (0=off, 1=on)',iunit)
 call write_inopt(dumps_p_period,        'dumps_p_period',    'number of dumps per period (if save_period = 1)',iunit)
 call write_inopt(reinject_enabled,      'reinject_enabled',  'enable dynamic reinjection (0=off, 1=on)',iunit)
 call write_inopt(n_inject_period,       'n_inject_period',   'period between reinjections (periods)',iunit)
 call write_inopt(mass_loss_start,       'mass_loss_start',   'start time for mass-loss calculation (periods)',iunit)
 call write_inopt(mass_loss_end,         'mass_loss_end',     'end time for mass-loss calculation (periods)',iunit)
 call write_inopt(use_file_mdot,         'use_file_mdot',     'skip measurement phase (0=off, 1=on)',iunit)
 call write_inopt(update_L,              'update_L',          'update luminosity with pulsation (0=off, 1=on)',iunit)
 call write_inopt(verbose,               'verbose',           'enable verbose output (0=off, 1=on)',iunit)

end subroutine write_options_inject

!----------------------------------------------------------------
!+
!  Read options from .in file
!+
!----------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save      :: ngot = 0
 integer, parameter :: noptions = 20

 ierr   = 0
 imatch = .true.
 select case(trim(name))
 case('n_profile_points')
    read(valstring,*,iostat=ierr) n_profile_points
    if (n_profile_points <= 10) call fatal(label,'n_profile_points must be > 10')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    if (iboundary_spheres < 0) call fatal(label,'iboundary_spheres must be >= 0')
 case('min_mass_fraction')
    read(valstring,*,iostat=ierr) min_mass_fraction
    if (min_mass_fraction < 0.) call fatal(label,'min_mass_fraction must be >= 0')
 case('rho_power')
    read(valstring,*,iostat=ierr) rho_power_in
    if (rho_power_in <= 0.) call fatal(label,'rho_power must be > 0')
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 2.) call fatal(label,'r_min_on_rstar must be in (0,2)')
 case('rho_inner')
    read(valstring,*,iostat=ierr) rho_inner
    if (rho_inner <= 0.) call fatal(label,'rho_inner must be > 0')
 case('iwind')
    read(valstring,*,iostat=ierr) iwind
    if (iwind /= 1 .and. iwind /= 2) call fatal(label,'iwind must be 1 or 2')
 case('pulsation_period')
    read(valstring,*,iostat=ierr) pulsation_period_days
    if (pulsation_period_days < 0.) call fatal(label,'pulsation_period must be >= 0')
 case('piston_velocity')
    read(valstring,*,iostat=ierr) piston_velocity_km_s
    if (piston_velocity_km_s < 0.) call fatal(label,'piston_velocity must be >= 0')
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    if (phi0 < -pi .or. phi0 > pi) call fatal(label,'phi0 must be in (-pi,pi)')
 case('wss')
    read(valstring,*,iostat=ierr) wss
    if (wss <= 0. .or. wss > 10.) call fatal(label,'wss must be in (0,10]')
 case('save_period')
    read(valstring,*,iostat=ierr) save_period
    if (save_period /= 0 .and. save_period /= 1) call fatal(label,'save_period must be 0 or 1')
 case('dumps_p_period')
    read(valstring,*,iostat=ierr) dumps_p_period
    if (dumps_p_period < 0) call fatal(label,'dumps_p_period must be > 0')
 case('reinject_enabled')
    read(valstring,*,iostat=ierr) reinject_enabled
    if (reinject_enabled /= 0 .and. reinject_enabled /= 1) call fatal(label,'reinject_enabled must be 0 or 1')
 case('n_inject_period')
    read(valstring,*,iostat=ierr) n_inject_period
    if (n_inject_period <= 0) call fatal(label,'n_inject_period must be > 0')
 case('mass_loss_start')
    read(valstring,*,iostat=ierr) mass_loss_start
    if (mass_loss_start < 0.) call fatal(label,'mass_loss_start must be >= 0')
 case('mass_loss_end')
    read(valstring,*,iostat=ierr) mass_loss_end
    if (mass_loss_end <= 0.) call fatal(label,'mass_loss_end must be > 0')
 case('use_file_mdot')
    read(valstring,*,iostat=ierr) use_file_mdot
    if (use_file_mdot /= 0 .and. use_file_mdot /= 1) call fatal(label,'use_file_mdot must be 0 or 1')
 case('update_L')
    read(valstring,*,iostat=ierr) update_L
    if (update_L /= 0 .and. update_L /= 1) call fatal(label,'update_L must be 0 or 1')
 case('verbose')
    read(valstring,*,iostat=ierr) verbose
    if (verbose /= 0 .and. verbose /= 1) call fatal(label,'verbose must be 0 or 1')
 case default
    imatch = .false.
 end select
 if (imatch) ngot = ngot + 1
 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject