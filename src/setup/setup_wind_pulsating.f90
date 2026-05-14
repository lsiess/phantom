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
 real :: semi_major_axis, semi_major_axis_au, eccentricity,f
 real :: default_particle_mass
 real :: primary_lum_lsun,primary_mass_msun,primary_Reff_au,primary_racc_au
 real :: primary_lum,primary_mass,primary_Reff,primary_racc,primary_Teff
 real :: secondary_lum_lsun,secondary_mass_msun,secondary_Reff_au,secondary_racc_au
 real :: secondary_lum,secondary_mass,secondary_Reff,secondary_racc,secondary_Teff
 real :: primary_veq,primary_veq_km_s,secondary_veq,secondary_veq_km_s,spin(2,3)
 real :: primary_mdot_msun_yr,primary_vwind_km_s,secondary_mdot_msun_yr,secondary_vwind_km_s
 real :: primary_mdot,primary_vwind,primary_wind_temp,secondary_mdot,secondary_vwind,secondary_wind_temp
 real :: primary_alpha,secondary_alpha

contains
!----------------------------------------------------------------
!+
!  Sets the default parameters for a pulsating AGB star with wind
!+
!----------------------------------------------------------------
subroutine set_default_parameters_wind()

 wind_gamma             = 5./3.
 spin                   = 0.
 spin(:,3)              = 1.  !spin along z-axis
 icompanion_star        = 0
 semi_major_axis        = 4.0
 eccentricity           = 0.
 primary_Teff           = 3000.
 ! placeholder default value
 secondary_Teff         = 1000.
 semi_major_axis_au     = 4.0
 f                      = 180.
 default_particle_mass  = 1.e-11
 primary_lum_lsun       = 5315.
 primary_mass_msun      = 1.5
 primary_Reff_au        = 1.2
 primary_racc_au        = 0.2
 primary_mdot_msun_yr   = 0.
 primary_vwind_km_s     = 0.
 primary_wind_temp      = 0.
 primary_alpha          = 0.
 primary_veq            = 0.
 primary_veq_km_s       = 0.
 ! placeholder default value
 secondary_lum_lsun     = 1000.
 secondary_mass_msun    = 1.0
 secondary_Reff_au      = 0.8
 secondary_racc_au      = 0.1
 secondary_mdot_msun_yr = .0
 secondary_vwind_km_s   = 0.
 secondary_wind_temp    = 0.
 secondary_alpha        = 0.
 secondary_veq          = 0.
 secondary_veq_km_s     = 0.

end subroutine set_default_parameters_wind

!----------------------------------------------------------------
!+
!  Setup routine for pulsating AGB star with wind
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use part,            only:xyzmh_ptmass,vxyz_ptmass,nptmass,igas,iTeff,iLum,iReff, &
                           ispinx,ispiny,ispinz,ivwind,imloss,iTwind,iwalpha
 use physcon,         only:au,solarm,mass_proton_cgs,kboltz,solarl,km
 use units,           only:umass,set_units,unit_velocity,utime,unit_energ,udist
 use inject,          only:set_default_options_inject
 use wind_pulsating,  only:setup_star,calc_stellar_profile,save_stellarprofile
 use setbinary,       only:set_binary
 use io,              only:master
 use options,         only:nfulldump
 use timestep,        only:dtmax!
 use eos,             only:gmw
 use infile_utils,    only:get_options
 use kernel,          only:hfact_default
 integer,          intent(in)    :: id
 integer,          intent(inout) :: npart
 integer,          intent(out)   :: npartoftype(:)
 real,             intent(out)   :: xyzh(:,:)
 real,             intent(out)   :: vxyzu(:,:)
 real,             intent(out)   :: massoftype(:)
 real,             intent(out)   :: polyk,gamma,hfact
 real,             intent(inout) :: time
 character(len=*), intent(in)    :: fileprefix
 character(len=len(fileprefix)+6) :: filename
 integer :: ierr,k
 logical :: iexist
 real :: omega_corotate, posang_ascnode, arg_peri, incl
 real :: T_wind

 nfulldump = 1
 dtmax = 0.1
!  ieos = 5
!  gmw = 1.26

 hfact = hfact_default
 call set_units(dist=au,mass=solarm,G=1.)
 call set_default_parameters_wind()
 filename = trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) call set_default_options_inject()

!--general parameters
!
 time = 0.
 if (id==master) print "(/,65('-'),1(/,a),/,65('-'),/)",' Pulsating wind setup'
 call get_options(trim(fileprefix)//'.setup',id==master,ierr,&
      read_setupfile,write_setupfile,setup_interactive)
 if (ierr /= 0) stop 'rerun phantomsetup after editing .setup file'
!
!--space available for injected gas particles
!
 npart = 0
 npartoftype(:) = 0
 xyzh(:,:)  = 0.
 vxyzu(:,:) = 0.
 xyzmh_ptmass(:,:) = 0.
 vxyz_ptmass(:,:) = 0.
 omega_corotate = 0.
 posang_ascnode = 0.
 arg_peri       = 0.
 incl           = 0.

 if (icompanion_star == 1) then
    call set_binary(primary_mass, &
                    secondary_mass, &
                    semi_major_axis, &
                    eccentricity, &
                    primary_racc, &
                    secondary_racc, &
                    xyzmh_ptmass, vxyz_ptmass, nptmass, ierr, &
                    posang_ascnode=posang_ascnode,&
                    arg_peri=arg_peri,&
                    incl=incl,&
                    f=f)
    xyzmh_ptmass(iTeff,1) = primary_Teff
    xyzmh_ptmass(iReff,1) = primary_Reff
    xyzmh_ptmass(iLum,1)  = primary_lum
    xyzmh_ptmass(imloss,1) = primary_mdot
    xyzmh_ptmass(ivwind,1) = primary_vwind
    xyzmh_ptmass(iTwind,1) = primary_wind_temp
    xyzmh_ptmass(iwalpha,1) = primary_alpha
    primary_veq = primary_veq_km_s * (km / unit_velocity)
    xyzmh_ptmass(ispinx,1) = primary_Reff**2*spin(1,1)*primary_veq
    xyzmh_ptmass(ispiny,1) = primary_Reff**2*spin(1,2)*primary_veq
    xyzmh_ptmass(ispinz,1) = primary_Reff**2*spin(1,3)*primary_veq

    xyzmh_ptmass(iTeff,2) = secondary_Teff
    xyzmh_ptmass(iReff,2) = secondary_Reff
    xyzmh_ptmass(iLum,2)  = secondary_lum
    xyzmh_ptmass(imloss,2) = secondary_mdot
    xyzmh_ptmass(ivwind,2) = secondary_vwind
    xyzmh_ptmass(iTwind,2) = secondary_wind_temp
    xyzmh_ptmass(iwalpha,2) = secondary_alpha
    secondary_veq = secondary_veq_km_s * (km / unit_velocity)
    xyzmh_ptmass(ispinx,2) = secondary_Reff**2*spin(2,1)*secondary_veq
    xyzmh_ptmass(ispiny,2) = secondary_Reff**2*spin(2,2)*secondary_veq
    xyzmh_ptmass(ispinz,2) = secondary_Reff**2*spin(2,3)*secondary_veq
    print *,'Sink particles summary'
    print *,'  #    mass       racc      lum         Reff'
    do k=1,nptmass
       print '(i4,2(2x,f9.5),2(2x,es10.3))',k,xyzmh_ptmass(4:5,k),xyzmh_ptmass(iLum,k)/(solarl*utime/unit_energ),&
               xyzmh_ptmass(iReff,k)*udist/au
    enddo
    print *,''

 else
    nptmass = 1
    xyzmh_ptmass(4,1)      = primary_mass
    xyzmh_ptmass(5,1)      = primary_racc
    xyzmh_ptmass(iTeff,1)  = primary_Teff
    xyzmh_ptmass(iReff,1)  = primary_Reff
    xyzmh_ptmass(iLum,1)   = primary_lum
    xyzmh_ptmass(imloss,1) = primary_mdot
    xyzmh_ptmass(ivwind,1) = primary_vwind
    xyzmh_ptmass(iTwind,1) = primary_wind_temp
    xyzmh_ptmass(iwalpha,1) = primary_alpha
    primary_veq = primary_veq_km_s * (km / unit_velocity)
    xyzmh_ptmass(ispinx,1) = primary_Reff**2*spin(1,1)*primary_veq
    xyzmh_ptmass(ispiny,1) = primary_Reff**2*spin(1,2)*primary_veq
    xyzmh_ptmass(ispinz,1) = primary_Reff**2*spin(1,3)*primary_veq
 endif

 ! This is overwritten anyway, calculated by M_atmosphere/N_particles
 massoftype(igas) = default_particle_mass * (solarm / umass)

 T_wind = 0.
 gamma = wind_gamma
 polyk = kboltz*T_wind/(mass_proton_cgs * gmw * unit_velocity**2)

end subroutine setpart

!----------------------------------------------------------------
!+
!  determine which problem to set up interactively
!+
!----------------------------------------------------------------
subroutine setup_interactive()
 use prompting, only:prompt
 use physcon,   only:au,solarm
 use units,     only:umass,udist
 use io,        only:fatal
 integer :: ichoice

 icompanion_star = 0
 call prompt('Add binary?',icompanion_star,0,2)
 if (icompanion_star == 1) then
    print "(a)",'Primary star parameters'
 else
    print "(a)",'Stellar parameters'
 endif
 ! define primary properties, wind and spin characteristics
 call get_sink_properties(primary_mass_msun,primary_mass,primary_racc_au,primary_racc,&
                          primary_Reff_au,primary_Reff,primary_Teff,primary_lum_lsun,primary_lum)
 call get_sink_wind(primary_mdot_msun_yr,primary_vwind_km_s,primary_wind_temp,primary_alpha)
 call get_sink_spin(spin(1,:),primary_veq_km_s)

 if (icompanion_star == 1) then
    print "(/,a)",'Secondary star parameters'
    ! define secondary properties, wind and spin characteristics
    call get_sink_properties(secondary_mass_msun,secondary_mass,secondary_racc_au,&
                             secondary_racc,secondary_Reff_au,secondary_Reff,secondary_Teff,&
                             secondary_lum_lsun,secondary_lum)
    call get_sink_wind(secondary_mdot_msun_yr,secondary_vwind_km_s,secondary_wind_temp,secondary_alpha)
    call get_sink_spin(spin(2,:),secondary_veq_km_s)

    ichoice = 1
    print "(/,a)",'Orbital parameters'
    print "(a)",' 1: semi-axis = 3.7 au, eccentricity = 0, true-anom = 180',' 0: custom'
    call prompt('select semi-major axis, ecccentricity, and true anomaly',ichoice,0,1)
    select case(ichoice)
    case(1)
       semi_major_axis_au = 3.7
       eccentricity       = 0.
    case default
       semi_major_axis_au = 1.
       eccentricity       = 0.
       call prompt('enter semi-major axis in au',semi_major_axis_au,0.,100.)
       call prompt('enter eccentricity',eccentricity,0.)
       call prompt('enter true anomaly',f,0.,360.)
    end select
    semi_major_axis = semi_major_axis_au * au / udist
 endif

end subroutine setup_interactive

!--------------------------------------------------------
!+
! set spin properties of the sink particle :
! rotation axis + spin rate in unit of critical spin rate
!+
!--------------------------------------------------------
subroutine get_sink_spin(xyzspin,wind_rotation_speed_km_s)
 use prompting, only:prompt
 real, intent(inout) :: xyzspin(3)
 real, intent(inout) :: wind_rotation_speed_km_s
 integer :: ichoice
 real :: spin

 ichoice = 2
 print "(a)",' 2: No rotation',&
      ' 1: Equatorial velocity = 10 km/s', &
      ' 0: custom'
 call prompt('select equatorial velocity',ichoice,0,2)
 select case(ichoice)
 case(2)
    wind_rotation_speed_km_s = 0.
 case(1)
    wind_rotation_speed_km_s = 10.
 case default
    ! would be more interesting to input the omega parameter instead
    call prompt('enter wind rotation speed in km/s',wind_rotation_speed_km_s,0.,1000.)
 end select
 if (wind_rotation_speed_km_s /= 0.) then
    ichoice = 1
    print "(a)",'spin orientation (vector)'
    print "(a)",' 1: along z-axis',&
         ' 0: custom'
    call prompt('select spin orientation',ichoice,0,1)
    if (ichoice == 0) then
       call prompt('enter x-component',xyzspin(1),-1.,1.)
       call prompt('enter y-component',xyzspin(2),-1.,1.)
       call prompt('enter z-component',xyzspin(3),-1.,1.)
       ! renormalize vector
       spin = sqrt(xyzspin(1)**2+xyzspin(2)**2+xyzspin(3)**2)
       xyzspin(1) = xyzspin(1)/spin
       xyzspin(2) = xyzspin(2)/spin
       xyzspin(3) = xyzspin(3)/spin
    endif
 endif
end subroutine get_sink_spin

!--------------------------------------------------------
!+
! set wind properties of the sink particle :
! wind mass loss rate, velocity and temperature
!+
!--------------------------------------------------------
subroutine get_sink_wind(wind_mdot_msun_yr,wind_speed_km_s,wind_temp,wind_alpha)
 use prompting, only:prompt
 real, intent(inout) :: wind_mdot_msun_yr,wind_speed_km_s,wind_temp,wind_alpha
 integer :: ichoice

 ichoice = 2
 print "(a)",'Wind properties'
 print "(a)",' 2: No wind',&
      ' 1: Mass loss rate = 1e-7 Msun/yr', &
      ' 0: custom'
 call prompt('select wind mass loss rate',ichoice,0,2)
 select case(ichoice)
 case(2)
    wind_mdot_msun_yr = 0.
 case(1)
    wind_mdot_msun_yr = 1.e-7
 case default
    call prompt('enter wind mass loss rate',wind_mdot_msun_yr,0.,1.)
 end select
 if (wind_mdot_msun_yr /= 0.) then
    ichoice = 1
    print "(a)",' 1: Wind velocity = 10 km/s',&
         ' 0: custom'
    call prompt('select wind velocity',ichoice,0,1)
    select case (ichoice)
    case(1)
       wind_speed_km_s = 10.
    case default
       call prompt('enter wind speed in km/s',wind_speed_km_s,0.,10000.)
    end select
    ichoice = 1
    print "(a)",' 2: Wind temperature = Teff',&
         ' 1: Wind temperature = 3000K',&
         ' 0: custom'
    call prompt('select wind temperature',ichoice,0,2)
    select case (ichoice)
    case(2)
       ! ERROR: primary_wind_temp = -1.000 too small [0.000:0.1000E+09]
       wind_temp = -1.
    case(1)
       wind_temp = 3000.
    case default
       call prompt('enter wind temperature in K',wind_temp,1.,1.e7)
    end select
    ichoice = 1
    print "(a)",' 1: Wind parameter alpha = 0.',&
         ' 0: custom'
    call prompt('select wind velocity',ichoice,0,1)
    select case (ichoice)
    case(1)
       wind_alpha = 0.
    case default
       call prompt('enter wind alpha parameter',wind_alpha,0.,10.)
    end select
 else
    wind_speed_km_s = 0.
    wind_temp = 0.
    wind_alpha = 0.
 endif

end subroutine get_sink_wind

!--------------------------------------------------------
!+
! set properties of the sink particle :
! mass, effective radius, accretion radius, luminosity
! and surface temperature
!+
!--------------------------------------------------------
subroutine get_sink_properties(sink_mass_msun,sink_mass,sink_racc_au,sink_racc,&
                               sink_Reff_au,sink_Reff,sink_Teff,sink_lum_lsun,sink_lum)
 use prompting, only:prompt
 use physcon,   only:au,solarm,solarl
 use units,     only:umass,udist,unit_luminosity
 real, intent(inout) :: sink_mass_msun,sink_racc_au,sink_Reff_au,sink_Teff,sink_lum_lsun
 real, intent(inout) :: sink_mass,sink_racc,sink_Reff,sink_lum
 integer :: ichoice

 ichoice = 1
 print "(a)",' 1: Mass = 1.5 Msun, accretion radius = 0.2568 au',&
      ' 2: Mass = 1 Msun, accretion radius = 0.2568 au',&
      ' 0: custom'
 call prompt('select mass and radius',ichoice,0,3)
 select case(ichoice)
 case(2)
    sink_mass_msun = 1.
    sink_racc_au   = 0.2568
 case(1)
    sink_mass_msun = 1.5
    sink_racc_au   = 0.2568
 case default
    call prompt('enter stellar mass',sink_mass_msun,0.,100.)
    call prompt('enter accretion radius in au ',sink_racc_au,0.)
 end select
 ichoice = 1
 print "(a)",' 1: Effective radius = 1 au',&
    ' 2: Effective radius = 0.1 au',&
    ' 0: custom'
 call prompt('select effective radius',ichoice,0,2)
 select case(ichoice)
 case(2)
    sink_Reff_au = 0.1
 case(1)
    sink_Reff_au = 1.2
 case default
    call prompt('enter effective radius in au ',sink_Reff_au,0.1)
 end select
 ichoice = 1
 print "(a)",' 1: Effective temperature = 3000 K, luminosity = 4770 Lsun',&
    ' 0: custom'
 call prompt('select effective temperature and luminosity',ichoice,0,1)
 select case(ichoice)
 case(1)
    sink_Teff     = 3000.
    sink_lum_lsun = 4770.
 case default
    call prompt('enter effective temperature ',sink_Teff,0.)
    call prompt('enter luminosity ',sink_lum_lsun,0.)
 end select

 sink_mass = sink_mass_msun * (solarm / umass)
 sink_racc = sink_racc_au * (au / udist)
 sink_Reff = sink_Reff_au * (au / udist)
 sink_lum  = sink_lum_lsun * (solarl / unit_luminosity)

end subroutine get_sink_properties

!----------------------------------------------------------------
!+
!  get luminosity and effective radius in code units
!  from various combinations of L, Teff and Reff in physical inuts
!+
!----------------------------------------------------------------
subroutine get_lum_and_Reff(lum_lsun,reff_au,Teff,lum,Reff)
 use physcon, only:au,steboltz,solarl,pi
 use units,   only:udist,unit_luminosity
 use io,      only:warning
 real, intent(inout) :: lum_lsun,reff_au,Teff
 real, intent(out)   :: lum,Reff
 real :: lum_tmp

 if (Teff <= tiny(0.) .and. lum_lsun > 0. .and. reff_au > 0.) then
    Teff = (lum_lsun*solarl/(4.*pi*steboltz*(reff_au*au)**2))**0.25
 elseif (Reff_au <= 0. .and. lum_lsun > 0. .and. Teff > 0.) then
    Reff_au = sqrt(lum_lsun*solarl/(4.*pi*steboltz*Teff**4))/au
 elseif (Reff_au > 0. .and. lum_lsun <= 0. .and. Teff > 0.) then
    lum_lsun = 4.*pi*steboltz*Teff**4*(reff_au*au)**2/solarl
 else
    lum_tmp = 4.*pi*steboltz*Teff**4*(reff_au*au)**2/solarl
    if (abs(lum_tmp-lum_lsun) > 1e-6) then
       print *,'[setup] inconsistent luminosity, should be ',lum_tmp,' instead of ',lum_lsun
       call warning('setup','L /= 4*pi*sigma*R^2*T^4')
    endif
 endif

 lum  = lum_lsun*(solarl/unit_luminosity)
 Reff = Reff_au*(au/udist)

end subroutine get_lum_and_Reff

!----------------------------------------------------------------
!+
!  write parameters to setup file
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
 call write_inopt(primary_mdot_msun_yr,'primary_mdot','primary wind mass loss rate (in Msun/yr)',iunit)
 call write_inopt(primary_vwind_km_s,'primary_vwind','primary wind velocity (in km/s)',iunit)
 call write_inopt(primary_wind_temp,'primary_wind_temp','primary wind temperature (K)',iunit)
 call write_inopt(primary_alpha,'primary_alpha','primary alpha parameter',iunit)
 call write_inopt(primary_veq_km_s,'primary_veq','primary equatorial velocity (in km/s)',iunit)
 if (primary_veq_km_s /= 0) then
    call write_inopt(spin(1,1),'primary_spinx','x-component of spin direction',iunit)
    call write_inopt(spin(1,2),'primary_spiny','y-component of spin direction',iunit)
    call write_inopt(spin(1,3),'primary_spinz','z-component of spin direction',iunit)
 endif
 call write_inopt(icompanion_star,'icompanion_star','set to 1 for a binary system, 2 for a triple system',iunit)
 if (icompanion_star == 1) then
    call get_lum_and_Reff(secondary_lum_lsun,secondary_Reff_au,secondary_Teff,secondary_lum,secondary_Reff)
    call write_inopt(secondary_mass_msun,'secondary_mass','secondary star mass (Msun)',iunit)
    call write_inopt(secondary_racc_au,'secondary_racc','secondary star accretion radius (au)',iunit)
    call write_inopt(secondary_lum_lsun,'secondary_lum','secondary star luminosity (Lsun)',iunit)
    call write_inopt(secondary_Teff,'secondary_Teff','secondary star effective temperature)',iunit)
    call write_inopt(secondary_Reff_au,'secondary_Reff','secondary star effective radius (au)',iunit)
    call write_inopt(secondary_mdot_msun_yr,'secondary_mdot','secondary wind mass loss rate (in Msun/yr)',iunit)
    call write_inopt(secondary_vwind_km_s,'secondary_vwind','secondary wind velocity (in km/s)',iunit)
    call write_inopt(secondary_wind_temp,'secondary_wind_temp','secondary wind temperature (K)',iunit)
    call write_inopt(secondary_alpha,'secondary_alpha','secondary alpha parameter',iunit)
    call write_inopt(secondary_veq_km_s,'secondary_veq','secondary equatorial velocity (in km/s)',iunit)
    if (secondary_veq_km_s /= 0) then
       call write_inopt(spin(2,1),'secondary_spinx','x-component of spin direction',iunit)
       call write_inopt(spin(2,2),'secondary_spiny','y-component of spin direction',iunit)
       call write_inopt(spin(2,3),'secondary_spinz','z-component of spin direction',iunit)
    endif
    call write_inopt(semi_major_axis_au,'semi_major_axis','semi-major axis of the binary system (au)',iunit)
    call write_inopt(eccentricity,'eccentricity','eccentricity of the binary system',iunit)
    call write_inopt(f,'true_anomaly','initial true anomaly of the binary orbit (deg)',iunit)
 endif

 call write_inopt(default_particle_mass,'mass_of_particles','particle mass (Msun, overwritten anyway <>0)',iunit)

 call write_inopt(wind_gamma,'wind_gamma','adiabatic index for wind gas',iunit)
 close(iunit)

end subroutine write_setupfile

!----------------------------------------------------------------
!+
!  Read parameters from setup file
!+
!----------------------------------------------------------------
subroutine read_setupfile(filename,ierr)
 use infile_utils, only:open_db_from_file,inopts,read_inopt,close_db
 use physcon,      only:au,steboltz,solarl,solarm,pi,km,years
 use units,        only:udist,umass,utime,unit_energ,unit_velocity
 use io,           only:error,fatal
 character(len=*), intent(in)  :: filename
 integer,          intent(out) :: ierr
 integer, parameter            :: iunit = 21
 type(inopts), allocatable     :: db(:)
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
 call read_inopt(primary_mdot_msun_yr,'primary_mdot',db,min=0.,max=1.,errcount=nerr)
 primary_mdot = primary_mdot_msun_yr * (solarm / umass) * (utime /years)
 call read_inopt(primary_vwind_km_s,'primary_vwind',db,min=0.,max=1.e4,errcount=nerr)
 primary_vwind = primary_vwind_km_s * (km / unit_velocity)
 call read_inopt(primary_wind_temp,'primary_wind_temp',db,min=0.,max=1.e8,errcount=nerr)
 call read_inopt(primary_alpha,'primary_alpha',db,min=0.,max=10.,errcount=nerr)
 call read_inopt(primary_veq_km_s,'primary_veq',db,min=0.,max=1000.,errcount=nerr)
 if (primary_veq_km_s /= 0) then
    call read_inopt(spin(1,1),'primary_spinx',db,min=-1.,max=1.,errcount=nerr)
    call read_inopt(spin(1,2),'primary_spiny',db,min=-1.,max=1.,errcount=nerr)
    call read_inopt(spin(1,3),'primary_spinz',db,min=-1.,max=1.,errcount=nerr)
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
    call read_inopt(secondary_mdot_msun_yr,'secondary_mdot',db,min=0.,max=1.,errcount=nerr)
    secondary_mdot = secondary_mdot_msun_yr * (solarm / umass) *(utime /years)
    call read_inopt(secondary_vwind_km_s,'secondary_vwind',db,min=0.,max=1.e4,errcount=nerr)
    secondary_vwind = secondary_vwind_km_s * (km / unit_velocity)
    call read_inopt(secondary_wind_temp,'secondary_wind_temp',db,min=0.,max=1.e8,errcount=nerr)
    call read_inopt(secondary_alpha,'secondary_alpha',db,min=0.,max=10.,errcount=nerr)
    call read_inopt(secondary_veq_km_s,'secondary_veq',db,min=0.,max=1000.,errcount=nerr)
    if (secondary_veq_km_s /= 0) then
       call read_inopt(spin(2,1),'secondary_spinx',db,min=-1.,max=1.,errcount=nerr)
       call read_inopt(spin(2,2),'secondary_spiny',db,min=-1.,max=1.,errcount=nerr)
       call read_inopt(spin(2,3),'secondary_spinz',db,min=-1.,max=1.,errcount=nerr)
    endif
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    call read_inopt(secondary_mass_msun,'secondary_mass',db,min=0.,max=1000.,errcount=nerr)
    secondary_mass = secondary_mass_msun * (solarm / umass)
    call read_inopt(semi_major_axis_au,'semi_major_axis',db,min=0.,errcount=nerr)
    semi_major_axis = semi_major_axis_au * au / udist
    call read_inopt(eccentricity,'eccentricity',db,min=0.,errcount=nerr)
    call read_inopt(f,'true_anomaly',db,min=0.,max=360.,errcount=nerr)
 endif

 call read_inopt(default_particle_mass,'mass_of_particles',db,min=0.,errcount=nerr)
 call read_inopt(wind_gamma,'wind_gamma',db,min=1.,max=4.,errcount=nerr)

 call close_db(db)
 ierr = nerr
 call write_setupfile(filename)

end subroutine read_setupfile

end module setup
