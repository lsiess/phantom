!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! various tests of the cooling solver module
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: cooling, cooling_functions, cooling_solver, dim,
!   dust_formation, options, physcon, prompting, units
!

 use cooling
 use cooling_functions
 use cooling_solver
 use physcon,          only:mass_proton_cgs,kboltz,atomic_mass_unit,patm,au,km
 use dust_formation,   only:init_muGamma,set_abundances,kappa_gas,calc_kappa_bowen,&
                              chemical_equilibrium_light,mass_per_H
 use dim,              only:nElements
 use prompting,  only:prompt
 use eos,        only:gmw,gamma,get_temperature_from_u,ieos,init_eos

 implicit none

 character(len=20), parameter, public :: analysistype = 'cooling'
 public :: do_analysis

 private
 integer :: analysis_to_perform
 real    :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
 real    :: eps(nElements) = [1.d0, 1.04d-1, 0.0,  6.d-4, 2.52d-4, 1.17d-4, 3.58d-5, 1.85d-5, 3.24d-5, 8.6d-8]

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:)
 real(kind=8),     intent(in) :: particlemass,time
 integer :: dump_number = 0


 !chose analysis type
 if (dump_number==0) then
    print "(29(a,/))", &
         ' 1) print cooling rates (for winds)', &
         ' 2) generate cooling table (for winds)', &
         ' 3) test integration', &
         ' 4) calc average wind profiles'

    analysis_to_perform = 4

    call prompt('Choose analysis type ',analysis_to_perform,1,4)
 endif

 !analysis
 select case(analysis_to_perform)
 case(1) !test rate
    call get_rate()
 case(2)
    call generate_grid()
 case(3)
    call test_cooling_solvers(dumpfile)
 case(4)
    call get_wind_features(xyzh,vxyzu,npart,time,particlemass,basename(dumpfile),dump_number)
 end select

 !increase dump number counter
 dump_number = dump_number + 1

end subroutine do_analysis


subroutine get_wind_features(xyzh,vxyzu,npart,time,particlemass,dumpfile,dump_number)

 use units, only : utime,udist,unit_velocity,unit_ergg,unit_density
 use part, only:rhoh
 real(kind=8),     intent(in) :: xyzh(:,:),vxyzu(:,:),time,particlemass
 integer,          intent(in) :: npart,dump_number
 character(len=*), intent(in) :: dumpfile

 integer :: i,j,iunit,irmin,irmax,ierr
 real(kind=8) :: rmin,rmax,r(npart),v(npart),vr(npart)
 real(kind=8) :: vrmax,vrmin,dr_bin,data_cols(4),gammai,mui,rhoi,ui
 real(kind=8), dimension(:), allocatable :: rbin,vbin,vrbin,ubin,hbin,Tbin,rhobin
 integer, dimension(:), allocatable :: bin_size
 integer, parameter :: ncols = 7
 integer :: nbin = 100
 character(len=17)   :: columns(ncols),colevol(4)

 !output data
 columns = (/'           r', &
             '           v', &
             '          vr', &
             '           u', &
             '           h', &
             '           T', &
             '         rho'/)
 colevol = (/'      r_vmin', &
             '       vrmin', &
             '      r_vmax', &
             '       vrmax'/)


 if (dump_number==0) then
    call prompt('Number of bins for the average profiles',nbin,1)
    call init_eos(ieos,ierr)
 endif

 allocate(rbin(nbin),vrbin(nbin),vbin(nbin),ubin(nbin),hbin(nbin),Tbin(nbin),rhobin(nbin),bin_size(nbin))

 rmin = 1e99
 rmax = 0.
 do i = 1,npart
    r(i) = sqrt(sum(xyzh(1:3,i)**2))
    v(i) = sqrt(sum(vxyzu(1:3,i)**2))
    vr(i) = sum(xyzh(1:3,i)*vxyzu(1:3,i))/r(i)
    rmin = min(rmin,r(i))
    rmax = max(rmax,r(i))
 enddo

 dr_bin = (rmax-rmin)/(nbin-1)
 rbin  = 0.
 vrbin = 0.
 vbin  = 0.
 ubin  = 0.
 hbin  = 0.
 bin_size = 0
!binning
 do i = 1,npart
    j = int((r(i)-rmin)/dr_bin)+1
    bin_size(j) = bin_size(j)+1
    rbin(j)  = rbin(j)+r(i)
    vbin(j)  = vbin(j)+v(i)
    vrbin(j) = vrbin(j)+vr(i)
    ubin(j)  = ubin(j)+vxyzu(4,i)
    hbin(j)  = hbin(j)+xyzh(4,i)
 enddo
!average binned quantities & make conversion
 do j = 1,nbin
    rhoi  = rhoh(hbin(j)/bin_size(j), particlemass)
    ui    = ubin(j)/bin_size(j)
    rhobin(j) = rhoi*unit_density
    rbin(j)  = rbin(j)/bin_size(j)*udist/au          !in au
    vbin(j)  = vbin(j)/bin_size(j)*unit_velocity/km  !in km/s
    vrbin(j) = vrbin(j)/bin_size(j)*unit_velocity/km !in km/s
    ubin(j)  = ui*unit_ergg         !in cgs
    hbin(j)  = hbin(j)/bin_size(j)*udist/au          !in au
    Tbin(j)  = get_temperature_from_u(ieos,0.,0.,0.,rhoi,ui)
 enddo

 vrmax = 0.
 vrmin = 1.e99
!do analysis on average profile
 do j = 1,nbin
!find max radial velocity & corresponding bin number
    if (vrbin(j)> vrmax) then
       vrmax = vrbin(j)
       irmax = j
    endif
    if (vrbin(j) < vrmin) then
       vrmin = vrbin(j)
       irmin = j
    endif
 enddo


 open(newunit=iunit,file=trim(dumpfile)//'.av',status='replace')
 write(iunit,'("nbin",2x,7(a15,3x))') columns
 do j = 1,nbin
    write(iunit,'(i6,7(es18.11e2,1x))') j,rbin(j),vbin(j),vrbin(j),ubin(j),hbin(j),Tbin(j),rhobin(j)
 enddo
 close(iunit)
 data_cols = (/rbin(irmin),vrmin,rbin(irmax),vrmax/)

!output quantities (1 line per dump)
 call write_time_file('wind_evol',colevol,time,data_cols,4,dump_number)

 deallocate(rbin,vrbin,vbin,ubin,hbin,Tbin,rhobin,bin_size)

end subroutine get_wind_features


function basename(filename) result(base)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=len(filename)) :: base
    integer :: p

    p = index(filename, '/', back=.true.)
    if (p > 0) then
        base = filename(p+1:)
    else
        base = filename
    end if
end function basename

subroutine write_time_file(name_in, cols, time, data_in, ncols, num)
 !outputs a file over a series of dumps
 character(len=*), intent(in) :: name_in
 integer, intent(in)          :: ncols, num
 character(len=*), dimension(ncols), intent(in) :: cols
 character(len=20), dimension(ncols) :: columns
 character(len=40)             :: data_formatter, column_formatter
 character(len(name_in)+9)    :: file_name
 real, intent(in)             :: time
 real, dimension(ncols), intent(in) :: data_in
 integer                      :: i, unitnum

 write(column_formatter, "(a,I2.2,a)") "('#',2x,", ncols+1, "('[',a15,']',3x))"
 write(data_formatter, "(a,I2.2,a)") "(", ncols+1, "(2x,es18.11e2))"
 write(file_name,"(2a,i3.3,a)") name_in, '.ev'

 if (num == 0) then
    unitnum = 1000

    open(unit=unitnum,file=file_name,status='replace')
    do i=1,ncols
       write(columns(i), "(I2,a)") i+1, cols(i)
    enddo

    !set column headings
    write(unitnum, column_formatter) '1         time', columns(:)
    close(unit=unitnum)
 endif

 unitnum=1001+num

 open(unit=unitnum,file=file_name, position='append')

 write(unitnum,data_formatter) time, data_in(:ncols)

 close(unit=unitnum)

end subroutine write_time_file


subroutine test_cooling_solvers(dumpfile)

 use prompting, only:prompt
 use physcon,   only:Rg
 use units,     only:unit_ergg,unit_density
 use options,   only:icooling
 ! For reading the input file:
 use units,            only:utime,umass,udist,set_units
 use infile_utils,     only:open_db_from_file,inopts,read_inopt,close_db
 use initial,          only:initialise
 use readwrite_infile, only:read_infile

 integer, parameter :: ndt = 100
 real :: tstart,tlast,dtstep,dti(ndt),tcool
 real :: rho, T_gas, rho_gas, pH, pH2   !rho in code units
 real :: mu, gamma
 real :: K2, kappa       !cgs
 real :: Q, dlnQ_dlnT
 real :: u,ui,dudt,T_on_u,Tout,dt,tcool0
 real :: T_floor
 ! For reading the input file:
 real :: utime_tmp,umass_tmp,udist_tmp
 character(len=*), intent(in)  ::   dumpfile
 character(len=120)            ::   infile,logfile,evfile,dfile


 integer :: i,imethod,ierr,iunit,ifunct,irate
 character(len=11) :: label

 !default cooling prescriptionHI
 excitation_HI  = 99
 icooling = 1
 icool_method = 1 !0=implicit, 1=explicit, 2=exact solution
 K2    = 0.
 kappa = 0.

 !temperature
 T_gas   = 1.d6
 rho_gas = 1.d-20 !cgs
 rho     = rho_gas/unit_density


 !
 !--store units, otherwise initialise() put them to 1
 !
 utime_tmp = utime
 umass_tmp = umass
 udist_tmp = udist

 !
 ! Read input file
 !
 infile = dumpfile(1:index(dumpfile,'_')-1)//'.in'
 call initialise()
 call set_units(udist_tmp,umass_tmp,utime_tmp)
 call read_infile(infile,logfile,evfile,dfile)


 call init_cooling_solver(ierr)
 call set_abundances
 call init_muGamma(rho_gas, T_gas, mu, gamma, pH, pH2)

 print "(29(a,/))", &
      'Select cooling function',&
      ' 1) Q = -a', &
      ' 2) Q = -b*T', &
      ' 3) Q = -c*T**3', &
      ' 4) Q = -d/T**3', &
      ' 5) HI cooling', &
      ' 6) piecewise law'

 irate = 0
 call prompt('Choose cooling rate ',irate)
 print *,''

 excitation_HI  = 99
 shock_problem = 0
 ifunct = 0
 select case(irate)
 case(2)
    ifunct = 1  !cooling rate depends on T
    label = '_exp.dat'
 case(3)
    ifunct = 3  !cooling rate depends on T**3
    label = '_T3.dat'
 case(4)
    ifunct = -3 !cooling rate depends on T**-3
    label = '_Tm3.dat'
 case(5)
    excitation_HI = 1
    label = '_HI.dat'
 case(6)
    shock_problem = 1
    excitation_HI = 98
    label = '_piece.dat'
 case default
    ifunct = 0  !cst coolint rate
    label = '_linear.dat'
 end select

 if (excitation_HI == 1) then
    Townsend_test = .true.
    T_floor = 10000.
    tstart = 0.01
 else
    K2 = ifunct
    Townsend_test = .false.
    T_floor = 10.
    tstart = 0.01
 endif

 call calc_cooling_rate(Q, dlnQ_dlnT, rho, T_gas, T_gas, mu, gamma, K2, kappa)
 tcool  = abs(kboltz*T_gas/((gamma-1.)*mu*atomic_mass_unit*Q*unit_ergg)) !code unit
 tcool0 = tcool
 T_on_u = (gamma-1.)*mu*unit_ergg/Rg
 ui     = T_gas/T_on_u


!set timesteps
 tlast = 10.
 dtstep = log10(tlast/tstart)/(ndt-1)
 do i = 1,ndt
    dti(i) = log10(tstart)+(i-1)*dtstep
 enddo
 dti = 10.**dti

 !--------------------------------------------------------------
 ! test solver integration over timesteps of different durations
 !--------------------------------------------------------------

 open(newunit=iunit,file='test_cooling_solvers'//trim(label),status='replace')

 do imethod = 0,2
    icool_method = imethod !0=implicit, 1=explicit, 2=exact solution
    print *,'#Tin=',T_gas,', rho_cgs=',rho_gas,', imethod=',icool_method,', cooling funct =',ifunct,excitation_HI,k2
    do i = 1,ndt
       dt = tcool0*dti(i)
       call energ_cooling_solver(ui,dudt,rho,dt,mu,gamma,0.,K2,0.)
       u = ui+dt*dudt
       Tout = max(u*T_on_u,T_floor)
       write(iunit,*) dti(i),dt,Tout,dudt,get_Texact(ifunct,T_gas,dt,tcool0,T_floor)
    enddo
 enddo
 close(iunit)

 !-----------------------------------------------------------
 ! test full integration integration over timestep of different durations
 !-----------------------------------------------------------

 !perform explicit integration
 icool_method = 1
 call integrate_cooling('test_cooling_explicit'//trim(label),ifunct,T_gas,T_floor,tcool0,tstart,ui,rho,mu,gamma)

 !perform implicit integration
 icool_method = 0
 call integrate_cooling('test_cooling_implicit'//trim(label),ifunct,T_gas,T_floor,tcool0,tstart,ui,rho,mu,gamma)

 !perform exact integration
 icool_method = 2
 call integrate_cooling('test_cooling_exact'//trim(label),ifunct,T_gas,T_floor,tcool0,tstart,ui,rho,mu,gamma)

end subroutine test_cooling_solvers

!-----------------------------------------------------------
! time integration of du/dt between t=0 and t=10*tcool
!-----------------------------------------------------------
subroutine integrate_cooling(file_in,ifunct,T_gas,T_floor,tcool0,tstart,ui,rho,mu,gamma)
 use units,   only:unit_ergg
 use physcon, only: Rg

 integer, intent(in) :: ifunct
 real, intent(in) :: tcool0,ui,rho,mu,gamma,T_gas,T_floor,tstart
 character(len=*), intent(in) :: file_in
 real :: time,dudt,dt,Tout,u,tend,dt_fact,T_on_u
 integer :: iunit

 T_on_u  = (gamma-1.)*mu*unit_ergg/Rg
 tend    = 10.*tcool0
 dt_fact = 0.1
 time    = 0.
 Tout    = T_gas
 dt      = tcool0*dt_fact
 u       = ui
 dudt    = 0.

 open(newunit=iunit,file=file_in,status='replace')
 write(iunit,*) tstart,dT,Tout,dudt,get_Texact(99,T_gas,time,tcool0,T_floor)
 do while (time < tend)! .and. Tout > T_floor)
    call energ_cooling_solver(u,dudt,rho,dt,mu,gamma,0.,dble(ifunct),0.)
    u = u+dt*dudt
    Tout = max(u*T_on_u,T_floor)
    time = time+dt
    !update dt based on the new value of tcool (bad for linear cooling!)
    !dt = dt_fact*abs(kboltz*Tout/((gamma-1.)*mu*atomic_mass_unit*dudt*unit_ergg))
    write(iunit,*) time/tcool0,dt,Tout,dudt,get_Texact(ifunct,T_gas,time,tcool0,T_floor)
 enddo
 close(iunit)

end subroutine integrate_cooling

!-------------------------------------------------------
! analytic solution for the evolution of the temperature
!-------------------------------------------------------

real function get_Texact(ifunct,T_gas,time,tcool0,T_floor)

 integer, intent(in) :: ifunct
 real, intent(in) :: T_gas,time,tcool0,T_floor

 select case(ifunct)
 case (0)
    get_Texact = max(T_gas*(1.-time/tcool0),T_floor)
 case (1)
    get_Texact = max(T_gas*exp(-time/tcool0),T_floor)
 case (3)
    get_Texact = max(1./sqrt(1./T_gas**2*(1.+2.*time/tcool0)),T_floor)
 case (-3)
    get_Texact = max(sqrt(sqrt(T_gas**4*(1.-4.*time/tcool0))),T_floor)
 case default
    get_Texact = T_gas
 end select

end function get_Texact

subroutine get_rate

 real :: T_gas, rho_gas, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas
 real :: pH, pH2
 real :: T_dust, v_drift, d2g, a, rho_grain, kappa_dust
 real :: JL
 real :: n_gas

 T_gas      = 1500.
 rho_gas    = 1.d-15

 call set_abundances
 call init_muGamma(rho_gas, T_gas, mu, gamma, pH, pH2)
 nH         = pH  *(patm*MPH(eps, Aw))/(mu*mass_proton_cgs*kboltz*T_gas)
 nH2        = pH2 *(patm*MPH(eps, Aw))/(mu*mass_proton_cgs*kboltz*T_gas)

 n_gas      = rho_gas/(mu*mass_proton_cgs)
 nHe        = 1.d-1*n_gas
 nCO        = 1.d-4*n_gas
 nH2O       = 5.d-5*n_gas
 nOH        = 1.d-7*n_gas

 kappa_gas  = 2.d-4

 T_dust     = 1000.
 v_drift    = 1.d6
 d2g        = 1./200.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_dust = calc_kappa_bowen(T_dust)

 JL         = 2.5d-12     ! Value taken from Barstow et al. 1997

 call print_cooling_rates(T_gas, rho_gas, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_gas, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)

end subroutine get_rate

subroutine generate_grid

 real :: logtmin,logtmax,logT,dlogt,T,crate,nH_tot,rho_cgs
 real :: pC, pC2, pC2H, pC2H2, mu, gamma, T_dust, d2g, v_drift
 real :: nH, nH2, nHe, nCO, nH2O, nOH, a, rho_grain, kappa_g, n_gas, kappa_dust, JL
 integer :: i,iunit
 integer, parameter :: nt = 400, iC=3

 logtmax = log10(1.e8)
 logtmin = log10(1.e0)

 call set_abundances()

 open(newunit=iunit,file='new_cooltable.txt',status='replace')
 write(iunit,"(a)") '#   T   \Lambda_E(T) erg s^{-1} cm^3   \Lambda erg s^{-1} cm^{-3}'
 dlogt = (logtmax - logtmin)/real(nt)
 d2g        = 1.d-2
 v_drift    = 0.
 a          = 1.d-5
 rho_grain  = 2.
 kappa_g    = 2.d-4
 JL         = 2.5d-12     ! Value taken from Barstow et al. 1997
 rho_cgs    = 1.d-15

 do i=1,nt
    logT   = logtmin + (i-1)*dlogt
    T      = 10**logT
    T_dust = T
    kappa_dust = calc_kappa_bowen(T)
    nH_tot     = rho_cgs/mass_per_H
    n_gas      = rho_cgs/(mu*mass_proton_cgs)
    call chemical_equilibrium_light(rho_cgs, T, eps(iC), pC, pC2, pC2H, pC2H2, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH)
    crate = calc_Q(T, rho_cgs, mu, nH, nH2, nHe, nCO, nH2O, nOH, kappa_g, &
                     T_dust, v_drift, d2g, a, rho_grain, kappa_dust, JL)
    !ndens = (rho_cgs/mass_proton_cgs)*5.d0/7.d0
    !print *,rho_cgs, T, mu, gamma, nH, nH2, nHe, nCO, nH2O, nOH, crate
    write(iunit,*) t,crate/nH_tot**2,crate
 enddo
 close(iunit)

end subroutine generate_grid

real function MPH(eps, Aw)

 real, dimension(nElements), intent(inout) :: eps, Aw
 real :: wind_CO_ratio

 wind_CO_ratio = 2.0
 eps(3)        = eps(4) * wind_CO_ratio
 MPH           = atomic_mass_unit*dot_product(Aw,eps)

end function MPH

end module analysis
