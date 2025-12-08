!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module inject
!
! Handles initial setup of stellar atmosphere with pulsating boundary layers
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   - iboundary_spheres  : *number of boundary spheres (integer)*
!   - n_shells_total     : *total number of atmospheric shells*
!   - iwind_resolution   : *geodesic sphere resolution*
!   - r_min_on_rstar     : *inner radius as fraction of R_star*
!   - pulsation_period   : *pulsation period (days)*
!   - pulsation_amplitude: *fractional pulsation amplitude*
!   - piston_velocity    : *piston velocity amplitude (km/s)*
!   - atmos_mass_fraction: *atmospheric mass as fraction of total stellar mass*
!   - surface_pressure   : *surface pressure (cgs)*
!   - iwind              : *wind type: 1=prescribed, 2=period from mass-radius relation*
!   - pulsation_timestep : *pulsation timestep as fraction of pulsation period*
!   - phi0               : *initial phase offset (radians) (best taken to be -pi/2 to start at minimum radius)*
!   - wss                : *fraction of tangential and radial distance between particles in initial atmosphere setup*
!
! :Dependencies: dim, eos, icosahedron, infile_utils, injectutils, io,
!   part, partinject, physcon, units, set_star
!
 use io,            only:fatal
 implicit none
 character(len=*), parameter, public :: inject_type = 'atmosphere'

 public :: init_inject,inject_particles,write_options_inject,read_options_inject,&
           set_default_options_inject,update_injected_par
 private

!--runtime settings for this module
!
! Read from input file
 integer :: iboundary_spheres = 10
 integer :: n_shells_total = 50
 integer :: n_profile_points = 10000
 integer :: iwind_resolution = 30
 real    :: r_min_on_rstar = 0.9
 real    :: dtpulsation = huge(0.)
 real    :: pulsation_period_days = 300.0  ! Pulsation period in days
 real    :: piston_velocity_km_s = 4.0     ! Piston velocity (in km/s)
 real    :: atmos_mass_fraction = 0.005  ! Atmosphere mass as fraction of total mass
 real    :: surface_pressure = 0.001  ! Surface pressure in cgs units
 integer :: iwind = 1  ! Wind type: 1=prescribed, 2=period from mass-radius relation
 real    :: pulsation_timestep = 0.02
 real    :: phi0 = -3.1415926536d0/2.0  ! Initial phase offset (-pi/2 for starting at minimal radius)
 real    :: wss = 2.0 ! Fraction of the tangential and radial distance between particles in the initial setup

! global variables
 integer, parameter :: wind_emitting_sink = 1
 real :: geodesic_R(0:19,3,3), geodesic_v(0:11,3)
 real :: omega_pulsation, deltaR_osc, pulsation_period, piston_velocity
 real :: Rstar, r_min, r_max, mass_of_particles
 real :: Mtotal, Matmos, Msink  ! Total, atmosphere, and sink masses
 real, allocatable :: delta_r_radial(:)
 integer :: particles_per_sphere, iresolution
 logical :: atmosphere_setup_complete = .false.
 
 ! Store boundary particle information
 real, allocatable    :: r_boundary_equilibrium(:)
 integer, allocatable :: boundary_particle_ids(:)
 integer              :: n_boundary_particles

 character(len=*), parameter :: label = 'inject_atmosphere'

contains

!-----------------------------------------------------------------------
!+
!  Set default options
!+
!-----------------------------------------------------------------------
subroutine set_default_options_inject(flag)
 integer, optional, intent(in) :: flag

 iboundary_spheres = 10
 n_shells_total = 50
 n_profile_points = 10000
 iwind_resolution = 30
 r_min_on_rstar = 0.9
 dtpulsation = huge(0.)
 atmos_mass_fraction = 0.005
 surface_pressure = 0.001
 iwind = 1
 pulsation_period_days = 300.0
 piston_velocity_km_s = 4.0
 pulsation_timestep = 0.02
 phi0 = -3.1415926536d0/2.0
 wss = 2.0

end subroutine set_default_options_inject

!-----------------------------------------------------------------------
!+
!  Initialize atmospheric setup and pulsation parameters
!+
!-----------------------------------------------------------------------
subroutine init_inject(ierr)
 use io,            only:fatal
 use physcon,       only:pi,days,au,solarm,km
 use icosahedron,   only:compute_matrices,compute_corners
 use eos,           only:gmw,gamma
 use units,         only:utime,umass,unit_velocity
 use part,          only:xyzmh_ptmass,massoftype,igas,iboundary,nptmass,iTeff,iReff
 use injectutils,   only:get_parts_per_sphere, get_neighb_distance
 use wind_pulsating,only:setup_star,calc_stellar_profile

 integer, intent(out) :: ierr
 real :: Mstar_cgs, Rstar_cgs, Tstar, delta_r_tangential, current_radius
 integer :: i 

 ierr = 0

 if (nptmass < 1) then
    call fatal(label,'need at least one sink particle for central star')
 endif

 ! Get stellar properties from sink particle
 Mtotal    = xyzmh_ptmass(4,wind_emitting_sink)
 Rstar     = xyzmh_ptmass(iReff,wind_emitting_sink)
 Rstar_cgs = Rstar * au 
 Mstar_cgs = Mtotal * solarm 
 Tstar     = xyzmh_ptmass(iTeff,wind_emitting_sink)

 ! Calculate mass distribution
 Matmos = atmos_mass_fraction * Mtotal
 Msink  = Mtotal - Matmos

 if (iwind == 2) then
    call calculate_period(Mtotal, Rstar, pulsation_period_days)
 endif

 ! Setup pulsation parameters
 pulsation_period = pulsation_period_days * (days / utime)
 omega_pulsation = 2.0*pi / pulsation_period
 piston_velocity = piston_velocity_km_s * (km / unit_velocity)
 deltaR_osc = pulsation_period * piston_velocity / (2.0*pi)

 print *, ''
 print *, 'Initializing pulsating atmosphere injection:'
 print *, 'pulsation period: ', pulsation_period
 print *, 'piston velocity: ', piston_velocity
 print *, 'deltaR_osc: ', deltaR_osc
 print *, ''
 
 ! Setup geodesic sphere parameters
 iresolution = iwind_resolution
 call compute_matrices(geodesic_R)
 call compute_corners(geodesic_v)
 particles_per_sphere = get_parts_per_sphere(iresolution)

 ! Allocate delta_r_radial array
 if (allocated(delta_r_radial)) deallocate(delta_r_radial)
 allocate(delta_r_radial(n_shells_total))

 r_min = r_min_on_rstar * Rstar
 current_radius = r_min
 
 do i = 1, n_shells_total
    delta_r_tangential = current_radius * get_neighb_distance(iresolution)
    delta_r_radial(i) = wss * delta_r_tangential
    current_radius = current_radius + delta_r_radial(i)
 enddo

 r_max = current_radius

 ! Setup stellar structure calculation
 call setup_star(Msink * umass, r_max * au, r_min * au, gmw, gamma, &
                 n_shells_total,surface_pressure, Matmos * umass)
 
 ! Calculate stellar profile for the atmosphere
 call calc_stellar_profile(n_profile_points)

 ! Calculate particle mass from atmospheric mass
 ! Total atmospheric mass distributed over all particles
 mass_of_particles = Matmos / real(n_shells_total * particles_per_sphere)

 print *, ''
 print *, 'Atmospheric particle mass (Msun): ', mass_of_particles
 print *, 'Particles per sphere: ', particles_per_sphere
 print *, 'Amount of particles: ', n_shells_total * particles_per_sphere
 print *, ''

 massoftype(igas) = mass_of_particles
 massoftype(iboundary) = mass_of_particles

 xyzmh_ptmass(4,wind_emitting_sink) = Msink

end subroutine init_inject

!-----------------------------------------------------------------------
!+
!  Main routine: called at the start to setup atmosphere,
!  then called each timestep to handle pulsation
!+
!-----------------------------------------------------------------------
subroutine inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)
 use part,        only:igas,iboundary,iamtype

 real,    intent(in)    :: time,dtlast
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart,npart_old
 integer, intent(inout) :: npartoftype(:)
 real,    intent(out)   :: dtinject

 ! Set timestep constraint for pulsation
 dtinject = pulsation_timestep * pulsation_period

 ! This is neccesary to not re-setup the atmosphere when resuming from a dump
 if (npart > 0 .and. .not. atmosphere_setup_complete) then
    atmosphere_setup_complete = .true.
 endif

 ! Initial setup: create all shells
 if (.not. atmosphere_setup_complete) then
    print *, 'Setting up stellar atmosphere with ', n_shells_total, ' shells.'
    call setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
    atmosphere_setup_complete = .true.
    print *, 'Stellar atmosphere setup complete.'
    return
 endif

 ! Every subsequent call, move the boundary particles
 call apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
!  call apply_pulsation_new(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass,npartoftype)

end subroutine inject_particles

!-----------------------------------------------------------------------
!+
!  Setup initial atmosphere with all shells at t=0
!+
!-----------------------------------------------------------------------
subroutine setup_initial_atmosphere(xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype
 use injectutils, only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon,     only:pi,km, au
 use units,       only:udist, unit_density, unit_ergg, unit_pressure

 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real,    intent(in)    :: xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i,j,first_particle,ipart_type,nboundary
 real    :: r,dr(n_shells_total),rho,u,T,P,x0(3),v0(3),GM,v_radial, r_previous
 logical :: is_boundary

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 ! Shell spacing
!  dr = (r_max - r_min) / real(n_shells_total - 1)
 dr = delta_r_radial
 print *, 'Shell spacing dr:', dr

 r_previous = r_min 

 ! Create shells from inner to outer
 npart = 0
 do i = 1, n_shells_total

    ! Calculate radius for this shell
   !  r = (r_min + (i-1)*dr)
    r = (r_previous + delta_r_radial(i))
    r_previous = r

    ! Determine if this is a boundary or free shell
    is_boundary = (i <= iboundary_spheres)
    
    ! Get stellar properties at this radius from 1D stellar profile
    ! This interpolates on the stellar_1D array calculated by set_star
    call interp_stellar_profile(r, rho, P, u, T)

   !  v_radial = piston_velocity
    v_radial = 0.0
    
    ! Set particle type - this tagging ensures forces are handled correctly
    ! (see partinject.f90 where boundary particles get special treatment)
    if (is_boundary) then
       ipart_type = iboundary
    else
       ipart_type = igas
    endif
    
    ! Inject this shell using geodesic sphere
    first_particle = npart + 1
    call inject_geodesic_sphere(i, first_particle, iresolution, r, v_radial, u, rho, &
                                geodesic_R, geodesic_V, npart, npartoftype, &
                                xyzh, vxyzu, ipart_type, x0, v0)
 enddo

 ! Store information about boundary particles for pulsation
 nboundary = npartoftype(iboundary)
 n_boundary_particles = nboundary

 print *, 'Number of boundary particles: ', nboundary
 print *, 'Number of gas particles: ', npartoftype(igas)
 
 if (nboundary > 0) then
    allocate(r_boundary_equilibrium(nboundary))
    allocate(boundary_particle_ids(nboundary))
    
    ! Store equilibrium radii and IDs of boundary particles
    j = 0
    do i = 1, npart
       if (j >= nboundary) then
          return
       elseif (iamtype(iphase(i)) == iboundary) then
          j = j + 1
          boundary_particle_ids(j) = i
          r_boundary_equilibrium(j) = sqrt( (xyzh(1,i)-x0(1))**2 + &
                                            (xyzh(2,i)-x0(2))**2 + &
                                            (xyzh(3,i)-x0(3))**2 ) &
                                            - deltaR_osc * sin(phi0)  !*omega_pulsation 
       endif
    enddo
 endif

end subroutine setup_initial_atmosphere

!-----------------------------------------------------------------------
!+
!  Apply radial pulsation to boundary particles
!+
!-----------------------------------------------------------------------
subroutine apply_pulsation(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass)
 use physcon, only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(in)    :: npart

 integer :: i,ipart
 real    :: r_eq,r_new,r_current,phase
 real    :: x_hat(3),r_dot
 real    :: x0(3),v0(3),GM
 real    :: x, y, z
 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 
 ! Pulsation amplitude and velocity
 ! R(t) = R0 + U_amp (P/2pi) * sin(omega t)
 ! dR/dt = U_amp * cos(omega t)
 r_dot = piston_velocity * cos(phase)

 ! Update each boundary particle
 do i = 1, n_boundary_particles
    ipart = boundary_particle_ids(i)
    
    ! Equilibrium radius for this particle
    r_eq = r_boundary_equilibrium(i)
    
    ! New radius with pulsation
    r_new = r_eq + deltaR_osc * sin(phase)

    x = xyzh(1,ipart) - x0(1)
    y = xyzh(2,ipart) - x0(2)
    z = xyzh(3,ipart) - x0(3)

    r_current = sqrt(x**2 + y**2 + z**2)
    
    ! Radial unit vector
    x_hat(1) = x / r_current
    x_hat(2) = y / r_current
    x_hat(3) = z / r_current

    xyzh(1,ipart) = r_new * x_hat(1) + x0(1)
    xyzh(2,ipart) = r_new * x_hat(2) + x0(2)
    xyzh(3,ipart) = r_new * x_hat(3) + x0(3)

    ! Update velocity (radial pulsation velocity)
    ! Scale velocity by ratio to equilibrium radius
    ! Add the orbital velocity of the sink
    vxyzu(1,ipart) = r_dot * x_hat(1) + v0(1)
    vxyzu(2,ipart) = r_dot * x_hat(2) + v0(2)
    vxyzu(3,ipart) = r_dot * x_hat(3) + v0(3)
 enddo

end subroutine apply_pulsation

subroutine apply_pulsation_new(time,xyzh,vxyzu,npart,xyzmh_ptmass,vxyz_ptmass,npartoftype)
 use part,        only:igas,iboundary,iphase,iamtype
 use injectutils, only:inject_geodesic_sphere
 use wind_pulsating, only:interp_stellar_profile
 use physcon, only:pi

 real,    intent(in)    :: time
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:),xyzmh_ptmass(:,:),vxyz_ptmass(:,:)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)

 integer :: i,ipart,ipart_type,first_particle
 real    :: r_eq,r_new,r_current,phase
 real    :: x_hat(3),r_dot,rho,u,T,P
 real    :: x0(3),v0(3),GM
 real    :: x, y, z

 if (.not. allocated(boundary_particle_ids)) return
 if (n_boundary_particles == 0) return

 ! Get sink particle position
 x0 = xyzmh_ptmass(1:3,wind_emitting_sink)
 v0 = vxyz_ptmass(1:3,wind_emitting_sink)
 GM = xyzmh_ptmass(4,wind_emitting_sink)

 phase = omega_pulsation * time + phi0
 
 ! Pulsation amplitude and velocity
 ! R(t) = R0 + U_amp (P/2pi) * sin(omega t)
 ! dR/dt = U_amp * cos(omega t)
 r_dot = piston_velocity * cos(phase)

 do i = 1, iboundary_spheres
   first_particle = (i-1) * particles_per_sphere + 1

   r_eq = r_boundary_equilibrium(first_particle)
   r_new = r_eq + deltaR_osc * sin(phase)

   call interp_stellar_profile(r_new, rho, P, u, T)   

   ipart_type = iboundary

   call inject_geodesic_sphere(i, first_particle, iresolution, r_new, r_dot, u, rho, &
                                geodesic_R, geodesic_V, npart, npartoftype, &
                                xyzh, vxyzu, ipart_type, x0, v0)

 enddo

end subroutine apply_pulsation_new


subroutine update_injected_par
 ! -- placeholder function
end subroutine update_injected_par

!-----------------------------------------------------------------------
!+
!  Calculate pulsation period based on stellar mass and radius
!  Using empirical relation from Ostlie & Cox (1986)
!+
!+-----------------------------------------------------------------------

subroutine calculate_period(M, R, pulsation_period_days)
 real, intent(in)  :: M, R
 real              :: logP, logM, logR
 real, intent(out) :: pulsation_period_days

 print *, 'Calculating pulsation period from mass-radius relation:'
 print *, 'Stellar mass (Msun): ', M
 print *, 'Stellar radius (Rsun): ', R

 logM = log10(M)
 logR = log10(R * 215.032)
 logP = -1.92 - 0.73*logM + 1.86*logR
 pulsation_period_days = 10.0**logP

 print *, 'Calculated pulsation period (days): ', pulsation_period_days

end subroutine calculate_period

!-----------------------------------------------------------------------
!+
!  Write options to input file
!+
!-----------------------------------------------------------------------
subroutine write_options_inject(iunit)
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 call write_inopt(n_profile_points,'n_profile_points', 'number of points in stellar profile',iunit)
 call write_inopt(n_shells_total,'n_shells_total', 'total number of atmospheric shells',iunit)
 call write_inopt(iboundary_spheres,'iboundary_spheres', 'number of boundary spheres (inner layers)',iunit)
 call write_inopt(iwind_resolution,'iwind_resolution', 'geodesic sphere resolution (integer)',iunit)
 call write_inopt(r_min_on_rstar,'r_min_on_rstar', 'inner radius as fraction of R_star',iunit)
 call write_inopt(atmos_mass_fraction,'atmos_mass_fraction', 'atmospheric mass as fraction of total stellar mass',iunit)
 call write_inopt(surface_pressure,'surface_pressure', 'surface pressure (cgs)',iunit)
 call write_inopt(iwind,'iwind','wind type: 1=prescribed, 2=period from mass-radius relation',iunit)
 call write_inopt(pulsation_period_days,'pulsation_period','pulsation period (days) (if iwind == 2 this is overwritten)',iunit)
 call write_inopt(piston_velocity_km_s,'piston_velocity','piston velocity amplitude (km/s)',iunit)
 call write_inopt(pulsation_timestep,'pulsation_timestep','pulsation timestep as fraction of pulsation period',iunit)
 call write_inopt(phi0,'phi0','initial phase offset (radians)',iunit)
 call write_inopt(wss,'wss','fraction of radial to tangential distance between particles in initial setup',iunit)

 
end subroutine write_options_inject

!-----------------------------------------------------------------------
!+
!  Read options from input file
!+
!-----------------------------------------------------------------------
subroutine read_options_inject(name,valstring,imatch,igotall,ierr)
 use io, only:fatal
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr

 integer, save :: ngot = 0
 integer, parameter :: noptions = 12
 logical :: init_opt = .false.

 if (.not. init_opt) then
    init_opt = .true.
    call set_default_options_inject()
 endif
 imatch = .true.
 igotall = .false.
 select case(trim(name))
 case('n_profile_points')
    read(valstring,*,iostat=ierr) n_profile_points
    ngot = ngot + 1
    if (n_profile_points <= 10) call fatal(label,'n_profile_points must be > 10')
 case('n_shells_total')
    read(valstring,*,iostat=ierr) n_shells_total
    ngot = ngot + 1
    if (n_shells_total <= 0) call fatal(label,'n_shells_total must be > 0')
 case('iboundary_spheres')
    read(valstring,*,iostat=ierr) iboundary_spheres
    ngot = ngot + 1
    if (iboundary_spheres < 0) call fatal(label,'iboundary_spheres must be >= 0')
    if (iboundary_spheres > n_shells_total) &
       call fatal(label,'iboundary_spheres must be <= n_shells_total')
 case('iwind_resolution')
    read(valstring,*,iostat=ierr) iwind_resolution
    ngot = ngot + 1
    if (iwind_resolution < 1) call fatal(label,'iwind_resolution must be >= 1')
 case('r_min_on_rstar')
    read(valstring,*,iostat=ierr) r_min_on_rstar
    ngot = ngot + 1
    if (r_min_on_rstar <= 0. .or. r_min_on_rstar >= 1.0) &
       call fatal(label,'r_min_on_rstar must be in range (0,1)')
 case('atmos_mass_fraction')
    read(valstring,*,iostat=ierr) atmos_mass_fraction
    ngot = ngot + 1
    if (atmos_mass_fraction <= 0. .or. atmos_mass_fraction >= 1.0) &
       call fatal(label,'atmos_mass_fraction must be in range (0,1)')
case('surface_pressure')
    read(valstring,*,iostat=ierr) surface_pressure
    ngot = ngot + 1
    if (surface_pressure < 0.) call fatal(label,'surface_pressure must be >= 0')
case('iwind')
    read(valstring,*,iostat=ierr) iwind
    ngot = ngot + 1
    if (iwind /= 1 .and. iwind /= 2) call fatal(label,'iwind must be 1 or 2')
case('pulsation_period')
    read(valstring,*,iostat=ierr) pulsation_period_days
    ngot = ngot + 1
    if (pulsation_period_days < 0.) call fatal(label,'pulsation_period must be >= 0')
 case('piston_velocity')
    read(valstring,*,iostat=ierr) piston_velocity_km_s
    ngot = ngot + 1
    if (piston_velocity_km_s < 0.) call fatal(label,'piston_velocity must be >= 0')
 case('pulsation_timestep')
    read(valstring,*,iostat=ierr) pulsation_timestep
    ngot = ngot + 1
    if (pulsation_timestep <= 0. .or. pulsation_timestep > 1.0) call fatal(label,'pulsation_timestep must be in range (0,1]')
 case('phi0')
    read(valstring,*,iostat=ierr) phi0
    ngot = ngot + 1
    if (phi0 < -3.1415926536d0 .or. phi0 > 3.1415926536d0) call fatal(label,'phi0 must be in range (-pi,pi)')
 case('wss')
    read(valstring,*,iostat=ierr) wss
    ngot = ngot + 1
    if (wss <= 0. .or. wss > 10.0) call fatal(label,'wss must be in range (0,10]')
 case default
    imatch = .false.
 end select

 igotall = (ngot >= noptions)

end subroutine read_options_inject

end module inject