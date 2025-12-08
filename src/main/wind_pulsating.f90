!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module wind_pulsating
!
! driver to integrate the hydrostatic equilibrium equations to set the outer layers of an AGB star
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters: None
!
! :Dependencies: eos, io, physcon, table_utils, units
!
 implicit none
 public :: setup_star
 public :: stellar_state,save_stellarprofile,interp_stellar_profile,calc_stellar_profile

 private
 ! Shared variables
 real, parameter :: rho_power = 2.0  ! Density profile exponent (i.e. rho ~ r^(-rho_power))

 ! input parameters
 real :: Mstar_cgs, Rstar_cgs, r_inner, Star_gamma, Star_mu, number_of_steps, P0, rho0, Menv_cgs
 real, dimension(:,:), allocatable, public :: stellar_1D

 ! stellar properties
 type stellar_state
    real :: r, r0, Rstar, rho, P, u, T
    integer :: nsteps
    logical :: error
 end type stellar_state

contains

subroutine setup_star(Mstar_in, Rstar_in, r_min, mu_in, gamma_in, n, surface_pressure, M_env_in)
 use physcon, only:au, solarm
!  use units,   only:umass,udist
!  use eos,     only:gamma, gmw

 real, intent(in)    :: Mstar_in, Rstar_in, r_min, mu_in, gamma_in, surface_pressure, M_env_in
 integer, intent(in) :: n

 Mstar_cgs  = Mstar_in
 Rstar_cgs  = Rstar_in
 r_inner    = r_min  ! Location where the stellar atmosphere is assumed to be inverse quadratic (i.e. inner boundary)
 Star_gamma = gamma_in
 Star_mu    = mu_in
 number_of_steps = n
 P0 = surface_pressure
 Menv_cgs = M_env_in

 print *, "Setting up star with parameters:"
 print *, "Rstar_cgs:", Rstar_cgs
 print *, "Mstar_cgs:", Mstar_cgs
 print *, "Menv_cgs:", Menv_cgs
 print *, "r_inner:", r_inner

end subroutine setup_star

!-----------------------------------------------------------------------
!
!  Initialize variables for stellar profile integration
!
!-----------------------------------------------------------------------
subroutine init_atmosphere(state)
! all quantities in cgs
 use physcon, only:pi, Rg, kboltz, mass_proton_cgs
 type(stellar_state), intent(out) :: state
 real :: C_rho

 ! Initialize at stellar surface
 state%r0     = Rstar_cgs
 state%r      = Rstar_cgs
 state%Rstar  = Rstar_cgs
 state%P      = P0
 C_rho        = Menv_cgs / (4.*pi * (Rstar_cgs - r_inner))
 rho0         = C_rho / Rstar_cgs**rho_power 
 state%rho    = rho0
 state%u      = state%P / (state%rho * (Star_gamma - 1.))
 state%T      = Star_mu * mass_proton_cgs / kboltz * (Star_gamma - 1.) * state%u
 print *, ""
 print *, "Initial stellar surface conditions:"
 print *, " mu:", Star_mu
 print *, " gamma:", Star_gamma
 print *, ""
 state%nsteps = 1
 state%error  = .false.

end subroutine init_atmosphere

!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equilibrium over one radial step
!
!-----------------------------------------------------------------------
subroutine stellar_step(state, r_new)
 use physcon, only:Gg, pi, Rg, kboltz, mass_proton_cgs

 type(stellar_state), intent(inout) :: state
 real, intent(in) :: r_new
 real :: dr, r_mid, rho_mid, mr_mid, dP, C_rho

 dr = r_new - state%r
 r_mid = 0.5 * (state%r + r_new)

 ! Get density and enclosed mass at midpoint
 C_rho = Menv_cgs / (4.*pi * (Rstar_cgs - r_inner))
 rho_mid = C_rho / r_mid**rho_power
 mr_mid = Mstar_cgs + Menv_cgs * ( r_mid - r_inner) / (Rstar_cgs - r_inner)

 ! Integrate hydrostatic equilibrium: dP/dr = -rho * G * M(r) / r^2
 dP = -(Gg * mr_mid * rho_mid / r_mid**2) * dr

 ! Update state
 state%r   = r_new
 state%P   = state%P + dP
 state%rho = C_rho / state%r**rho_power

 ! Calculate thermodynamic quantities
 state%u = state%P / (state%rho * (star_gamma - 1.))
 state%T = star_mu * mass_proton_cgs / kboltz * (star_gamma - 1.) * state%u

 state%nsteps = state%nsteps + 1

end subroutine stellar_step

!-----------------------------------------------------------------------
!
!  Integrate the hydrostatic equilibrium equation
!
!-----------------------------------------------------------------------
subroutine calc_stellar_profile(n)
! all quantities in cgs
 integer, intent(in) :: n
 type(stellar_state) :: state
 integer :: i
 real :: r_new, dr

 ! Initialize stellar structure
 call init_atmosphere(state)

 ! Allocate storage for profile
 if (allocated(stellar_1D)) deallocate(stellar_1D)
 allocate(stellar_1D(5, n))

 ! Store surface values
 stellar_1D(1, n) = state%r
 stellar_1D(2, n) = state%rho
 stellar_1D(3, n) = state%P
 stellar_1D(4, n) = state%u
 stellar_1D(5, n) = state%T

 ! Integrate inward from surface to center in steps of dr
 dr = (Rstar_cgs - r_inner) / real(n-1) 

 do i = n-1, 1, -1
    r_new = Rstar_cgs - real(n-i)*dr 
    call stellar_step(state, r_new)

    ! Store in profile
    stellar_1D(1, i) = state%r
    stellar_1D(2, i) = state%rho
    stellar_1D(3, i) = state%P
    stellar_1D(4, i) = state%u
    stellar_1D(5, i) = state%T
 enddo

 ! Save profile to file
call save_stellarprofile(n, 'stellar_profile1D.dat')

end subroutine calc_stellar_profile

!-----------------------------------------------------------------------
!
!  Interpolate stellar profile at given radius
!
!-----------------------------------------------------------------------
subroutine interp_stellar_profile(r, rho, P, u, T)
 !in/out variables in code units
 use units,       only:udist, unit_density, unit_ergg, unit_pressure
 use table_utils, only:find_nearest_index, interp_1d
 use io,          only:fatal

 real, intent(in)  :: r  ! in code units
 real, intent(out) :: rho, P, u, T
 real :: r_cgs
 integer :: indx, n
 character(len=*), parameter :: label = 'interp_stellar_profile'

 ! Check if profile exists
 if (.not. allocated(stellar_1D)) then
    call fatal(label, 'stellar_1D not allocated. Call setup_star first.')
 endif

 n = size(stellar_1D, 2)
 r_cgs = r * udist

 if (r_cgs <= stellar_1D(1,1)) then
    print *, 'Warning: r_cgs < stellar_1D(1,1). Extrapolating...'
    rho = stellar_1D(2, 1) / unit_density
    P   = stellar_1D(3, 1) / unit_pressure
    u   = stellar_1D(4, 1) / unit_ergg
    T   = stellar_1D(5, 1)
    return
 elseif (r_cgs >= stellar_1D(1, n)) then
    rho = stellar_1D(2, n) / unit_density
    P   = stellar_1D(3, n) / unit_pressure
    u   = stellar_1D(4, n) / unit_ergg
    T   = stellar_1D(5, n)
    return
 endif

 call find_nearest_index(stellar_1D(1,:), r_cgs, indx)

 rho = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1),stellar_1D(2,indx), stellar_1D(2,indx+1)) / unit_density
 P   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1),stellar_1D(3,indx), stellar_1D(3,indx+1)) / unit_pressure
 u   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1),stellar_1D(4,indx), stellar_1D(4,indx+1)) / unit_ergg
 T   = interp_1d(r_cgs, stellar_1D(1,indx), stellar_1D(1,indx+1),stellar_1D(5,indx), stellar_1D(5,indx+1)) 

end subroutine interp_stellar_profile

!-----------------------------------------------------------------------
!
!  Save stellar profile to file
!
!-----------------------------------------------------------------------
subroutine save_stellarprofile(n, filename)
 use physcon, only:au
 use io,      only:iverbose
 integer, intent(in) :: n
 character(*), intent(in) :: filename
 integer :: i, nwrite
 integer, parameter :: iunit = 1338

 if (iverbose >= 1) then
    write(*,'("Saving 1D stellar model to ",A)') trim(filename)
 endif

 open(unit=iunit, file=filename, status='replace')
 call filewrite_stellar_header(iunit, nwrite)

 ! Write profile data
 do i = 1, n
    call filewrite_stellar_state(iunit, nwrite, i)
 enddo
 close(iunit)

end subroutine save_stellarprofile

subroutine filewrite_stellar_header(iunit, nwrite)
 integer, intent(in)  :: iunit
 integer, intent(out) :: nwrite
 character(len=20) :: fmt

 nwrite = 5
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(a15))') 'r','rho','P','u','T'

end subroutine filewrite_stellar_header

subroutine state_to_array(i, array)
 use physcon, only:pi
 integer, intent(in) :: i
 real, intent(out) :: array(:)

 array(1) = stellar_1D(1, i)  ! r
 array(2) = stellar_1D(2, i)  ! rho
 array(3) = stellar_1D(3, i)  ! P
 array(4) = stellar_1D(4, i)  ! u
 array(5) = stellar_1D(5, i)  ! T

end subroutine state_to_array

subroutine filewrite_stellar_state(iunit, nwrite, i)
 integer, intent(in) :: iunit, nwrite, i
 real :: array(nwrite)
 character(len=20) :: fmt

 call state_to_array(i, array)
 write(fmt,*) nwrite
 write(iunit,'('// adjustl(fmt) //'(1x,es14.6E3:))') array(1:nwrite)

end subroutine filewrite_stellar_state

end module wind_pulsating