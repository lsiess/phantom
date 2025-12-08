!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Setup routine for pulsating AGB star with atmospheric wind
!
! :References: None
!
! :Owner: Owen Vermeulen
!
! :Runtime parameters:
!   wind_gamma            - Adiabatic index for wind gas
!   icompanion_star       - Set to 1 for a binary system, 0 for single star
!   semi_major_axis       - Semi-major axis of the binary system
!   eccentricity          - Eccentricity of the binary system
!   primary_Teff          - Primary star effective temperature (K)
!   primary_lum_lsun      - Primary star luminosity (Lsun)
!   primary_mass_msun     - Primary star mass (Msun)
!   primary_Reff_au       - Primary star effective radius (au)
!   primary_racc_au       - Primary star accretion radius (au)
!   secondary_lum_lsun    - Secondary star luminosity (Lsun)
!   secondary_mass_msun   - Secondary star mass (Msun)
!   secondary_Reff_au     - Secondary star effective radius (au)
!   secondary_racc_au     - Secondary star accretion radius (au)
!   mass_of_particles     - Particle mass (Msun, overwritten anyway <>0)
!
! :Dependencies: boundary, dim, infile_utils, io, kernel, options, part,
!   physcon, prompting, ptmass, setup_params, table_utils, timestep, units
!
 implicit none
 public :: setpart
  
 private
 real, public :: wind_gamma
 integer :: icompanion_star
 real :: semi_major_axis, semi_major_axis_au, eccentricity
 real :: primary_Teff, primary_lum_lsun, primary_mass_msun, primary_reff_au, primary_racc_au
 real :: secondary_lum_lsun, secondary_mass_msun, secondary_reff_au, secondary_racc_au
 real :: primary_Reff,primary_lum,primary_mass,primary_racc
 real :: secondary_Reff,secondary_Teff,secondary_lum,secondary_mass,secondary_racc
 real :: default_particle_mass

contains

!-----------------------------------------------------------------------
!+
!  Sets the default parameters for a pulsating AGB star with wind
!+
!+-----------------------------------------------------------------------

subroutine set_default_parameters_wind()
 wind_gamma            = 5./3.
 icompanion_star       = 0
 semi_major_axis       = 4.0
 semi_major_axis_au    = 4.0
 eccentricity          = 0.0
 primary_Teff          = 3000.0
 primary_lum_lsun      = 5315.
 primary_mass_msun     = 1.5
 primary_Reff_au       = 1.
 primary_racc_au       = 0.2
 secondary_lum_lsun    = 0.
 secondary_mass_msun   = 1.0
 secondary_Reff_au     = 0.
 secondary_racc_au     = 0.1
 default_particle_mass = 1.e-10

end subroutine set_default_parameters_wind

!----------------------------------------------------------------
!+
!  Setup routine for pulsating AGB star with wind
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,      only: xyzmh_ptmass, vxyz_ptmass, nptmass, igas, iTeff, iLum, iReff
 use physcon,   only: au, solarm, mass_proton_cgs, kboltz, solarl
 use units,     only: set_units,umass,udist,utime,unit_energ
!  use inject,    only: set_default_options_inject, inject_particles, init_inject
 use inject ,    only: set_default_options_inject
 use wind_pulsating, only: setup_star, calc_stellar_profile, save_stellarprofile
 use setbinary, only: set_binary
 use io,        only: master
 
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=*),  intent(in)    :: fileprefix
 character(len=len(fileprefix)+6) :: filename
 integer :: ierr,k
 logical :: iexist
 
 call set_units(mass=solarm,dist=au,G=1.)
 call set_default_parameters_wind()
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) call set_default_options_inject()

 time = 0.
 filename = trim(fileprefix)//'.setup'
 inquire(file=filename,exist=iexist)
 if (iexist) call read_setupfile(filename,ierr)
 if (.not. iexist .or. ierr /= 0) then
    if (id==master) then
       call write_setupfile(filename)
    endif
 endif

 gamma = wind_gamma

 !
 !--space available for injected gas particles
 !
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.

 if (icompanion_star == 1) then
    call set_binary(primary_mass, &
                    secondary_mass, &
                    semi_major_axis, &
                    eccentricity, &
                    primary_racc, &
                    secondary_racc, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass, ierr)
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
    xyzmh_ptmass(iTeff,2) = secondary_Teff
    xyzmh_ptmass(iReff,2) = secondary_Reff
    xyzmh_ptmass(iLum,2)  = secondary_lum
     print *,'Sink particles summary'
     print *,'  #    mass       racc      lum         Reff'
     do k=1,nptmass
        print '(i4,2(2x,f9.5),2(2x,es10.3))',k,xyzmh_ptmass(4:5,k),xyzmh_ptmass(iLum,k)/(solarl*utime/unit_energ),&
               xyzmh_ptmass(iReff,k)*udist/au
     enddo
     print *,''

 else
    nptmass = 1
    xyzmh_ptmass(4,1)     = primary_mass
    xyzmh_ptmass(5,1)     = primary_racc
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
 endif

 ! This is overwritten anyway, calculated by M_atmosphere/N_particles
 massoftype(igas) = default_particle_mass * (solarm / umass)

!  call init_inject(dummy)
! Setup stellar structure calculation
!  call setup_star(Matmos, Rstar_cgs, r_inner, gmw, gamma, n_shells_total,surface_pressure)
!  call calc_stellar_profile(n_profile_points)
!  call save_stellarprofile(n_profile_points, 'stellar_profile1D.dat')
!  call inject_particles(time,dtlast,xyzh,vxyzu,xyzmh_ptmass,vxyz_ptmass,npart,npart_old,npartoftype,dtinject)

end subroutine setpart

!----------------------------------------------------------------
!+
!  get luminosity and effective radius in code units
!  from various combinations of L, Teff and Reff in physical inuts
!+
!----------------------------------------------------------------
subroutine get_lum_and_Reff(lum_lsun,reff_au,Teff,lum,Reff)
 use physcon, only:au,steboltz,solarl,pi
 use units,   only:udist,unit_luminosity
 real, intent(inout) :: lum_lsun,reff_au,Teff
 real, intent(out)   :: lum,Reff

 if (Teff <= tiny(0.) .and. lum_lsun > 0. .and. Reff_au > 0.) then
    primary_Teff = (lum_lsun*solarl/(4.*pi*steboltz*(Reff_au*au)**2))**0.25
 elseif (Reff_au <= 0. .and. lum_lsun > 0. .and. Teff > 0.) then
    Reff_au = sqrt(lum_lsun*solarl/(4.*pi*steboltz*Teff**4))/au
 elseif (Reff_au > 0. .and. lum_lsun <= 0. .and. Teff > 0.) then
    lum_lsun = 4.*pi*steboltz*Teff**4*(primary_Reff_au*au)**2/solarl
 endif

 lum  = lum_lsun*(solarl/unit_luminosity)
 Reff = Reff_au*(au/udist)

end subroutine get_lum_and_Reff

!----------------------------------------------------------------
!+
!  Write setup parameters to input file
!+
!----------------------------------------------------------------
subroutine write_setupfile(filename)
 use infile_utils, only:write_inopt
 character(len=*), intent(in) :: filename
 integer, parameter           :: iunit = 20


 print "(a)",' writing setup options file '//trim(filename)
 open(unit=iunit,file=filename,status='replace',form='formatted')
 write(iunit,"(a)") '# input file for wind setup routine'

 call get_lum_and_Reff(primary_lum_lsun,primary_Reff_au,primary_Teff,primary_lum,primary_Reff)

 call write_inopt(primary_mass_msun,'primary_mass','primary star mass (Msun)',iunit)
 call write_inopt(primary_racc_au,'primary_racc','primary star accretion radius (au)',iunit)
 call write_inopt(primary_lum_lsun,'primary_lum','primary star luminosity (Lsun)',iunit)
 call write_inopt(primary_Teff,'primary_Teff','primary star effective temperature (K)',iunit)
 call write_inopt(primary_Reff_au,'primary_Reff','primary star effective radius (au)',iunit)
 call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system, 2 for a triple system',iunit)
 if (icompanion_star == 1) then
    call get_lum_and_Reff(secondary_lum_lsun,secondary_Reff_au,secondary_Teff,secondary_lum,secondary_Reff) 
    call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
    call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
    call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
    call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature (K)',iunit)
    call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
    call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
    call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
 endif

 call write_inopt(default_particle_mass,'mass_of_particles','particle mass (Msun, overwritten anyway <>0)',iunit)

 call write_inopt(wind_gamma,'wind_gamma','adiabatic index for wind gas',iunit)

 close(iunit)
 
end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read setup parameters from input file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use physcon,      only:au,steboltz,solarl,solarm,pi
 use units,        only:udist,umass,utime,unit_energ
 use io,           only:error,fatal
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter :: iunit = 21
 type(inopts), allocatable :: db(:)
 integer :: nerr,ichange


 nerr = 0
 ichange = 0
 print "(a)",' reading setup options from '//trim(filename)
 call open_db_from_file(db,filename,iunit,ierr)
 call read_inopt(primary_mass_msun,'primary_mass',db,min=0.,max=1000.,errcount=nerr)
 primary_mass = primary_mass_msun * (solarm / umass)
 call read_inopt(primary_lum_lsun,'primary_lum',db,min=0.,max=1e7,errcount=nerr)
 primary_lum = primary_lum_lsun * (solarl * utime / unit_energ)
 call read_inopt(primary_Teff,'primary_Teff',db,min=0.,errcount=nerr)
 call read_inopt(primary_Reff_au,'primary_Reff',db,min=0.,errcount=nerr)
 primary_Reff = primary_Reff_au * au / udist
 call read_inopt(primary_racc_au,'primary_racc',db,min=0.,errcount=nerr)
 primary_racc = primary_racc_au * au / udist
 if (primary_racc < tiny(0.)) then
    print *,'ERROR: primary accretion radius not defined'
    nerr = nerr+1
 endif

 call read_inopt(icompanion_star,'icompanion_star',db,min=0,errcount=nerr)
 if (icompanion_star == 1) then
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(secondary_lum_lsun,'secondary_lum',db,min=0.,max=1e7,errcount=nerr)
    secondary_lum = secondary_lum_lsun * (solarl * utime / unit_energ)
    call read_inopt(secondary_Teff,'secondary_Teff',db,min=0.,errcount=nerr)
    call read_inopt(secondary_Reff_au,'secondary_Reff',db,min=0.,errcount=nerr)
    secondary_Reff = secondary_Reff_au * au / udist
    call read_inopt(secondary_racc_au,'secondary_racc',db,min=0.,errcount=nerr)
    secondary_racc = secondary_racc_au * au / udist
    if (secondary_racc < tiny(0.)) then
       print *,'ERROR: secondary accretion radius not defined'
       nerr = nerr+1
    endif
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
 endif

 call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)
 call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)
 
 call close_db(db)
 
 call close_db(db)
 ierr = nerr
 call write_setupfile(filename)

 
end subroutine read_setupfile

end module setup