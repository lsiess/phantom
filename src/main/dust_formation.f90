!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dust_formation
!
! Dust formation routine : theory of moments
!
! :References: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!
! :Owner: Lionel Siess
!
! :Runtime parameters:
!   - bowen_Tcond   : *dust condensation temperature (K)*
!   - bowen_delta   : *condensation temperature range (K)*
!   - bowen_kmax    : *maximum dust opacity (cm²/g)*
!   - idust_opacity : *compute dust opacity (0=off, 1=bowen)*
!   - kappa_gas     : *constant gas opacity (cm²/g)*
!   - wind_CO_ratio : *wind initial C/O ratio (> 1)*
!
! :Dependencies: chemistry_condensation, chemistry_moments, dim,
!   dump_utils, eos, infile_utils, io, part, physcon
!

 implicit none
 integer, public :: idust_opacity = 0

 public :: set_abundances,calc_kappa_dust,calc_kappa_bowen,&
      read_options_dust_formation,write_options_dust_formation,&
      calc_Eddington_factor,calc_muGamma,init_dust_formation,&
      write_headeropts_dust_formation,read_headeropts_dust_formation
!
!--runtime settings for this module
!
 real, public :: kappa_gas   = 2.d-4
 real, public :: wind_CO_ratio = 2.
 real, parameter, public :: a_init_dust = 1.0d-7 !cgs units 1nm

!number of elements considered in the chemical networks
 integer, public, parameter :: nElements = 11
 real, public :: mass_per_H, eps(nElements)
 real, public :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867, 24.305]
! Indices for initial elemental composition
 integer, parameter, public :: iH=1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10, iMg=11

 private

 character(len=*), parameter :: label = 'dust_formation'
 real :: bowen_kmax  = 2.7991
 real :: bowen_Tcond = 1500.
 real :: bowen_delta = 60.

contains

!----------------------------------------------------------------
subroutine set_abundances
!----------------------------------------------------------------
! initialize abundances
 use physcon, only:atomic_mass_unit
 logical :: ferrarotti = .true. !A&A 2002, 382, 256
 ! These abundances are with respect to H, so take Table 2.1 of Gail and Sedlmayr (Book) and divide by the listed H abundance.
 if (ferrarotti) then
    eps(iH)  = 1.0
    eps(iHe) = 1.04d-1
    eps(iOx) = 6.87d-4
    eps(iN)  = 2.52d-4
    eps(iMg) = 3.85d-5
    eps(iNe) = 1.17d-4
    eps(iSi) = 3.58d-5
    eps(iS)  = 1.85d-5
    eps(iFe) = 3.24d-5
    eps(iTi) = 8.6d-8
 else
    eps(iHe) = 1.009d-1 !1.04d-1
    eps(iOx) = 7.211d-4 !6.64 6.87d-4
    eps(iN)  = 2.106d-4 !2.52d-4
    eps(iMg) = 3.859d-5 !3.85d-5
    eps(iNe) = 1.18d-4  !1.17d-4
    eps(iSi) = 3.561d-5 !3.58d-5
    eps(iS)  = 1.863d-5 !1.85d-5
    eps(iFe) = 3.241d-5 !3.24d-5
    eps(iTi) = 8.621d-8 !8.6d-8
 endif
 eps(iC)  = eps(iOx) * wind_CO_ratio
 mass_per_H = atomic_mass_unit*dot_product(Aw,eps)
 !XH  = atomic_mass_unit*eps(iH)/mass_per_H  ! H mass fraction
 !XHe = atomic_mass_unit*eps(iHe)/mass_per_H ! He mass fraction
end subroutine set_abundances

!----------------------------------------------------------------
subroutine calc_muGamma(rho_cgs,T,mu,gamma,pH_tot,pH,pH2,pressure_cgs)
!----------------------------------------------------------------
! all quantities are in cgs
 use eos, only:ieos
 use physcon, only:patm,kboltz
 real, intent(in)              :: rho_cgs
 real, intent(inout)           :: T, mu, gamma
 real, intent(out)             :: pH, pH_tot, pH2
 real, intent(in), optional    :: pressure_cgs
 real :: KH2
 real :: T_old, mu_old, gamma_old, tol
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma'

 if (present(pressure_cgs)) then
    pH_tot = pressure_cgs/(1.+eps(iHe))/patm
 else
    pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 endif
 if (T > 1.d4) then
    pH  = pH_tot
    pH2 = 0.
    if (ieos == 5) then
       mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
       gamma  = 5./3.
    endif
 elseif (T > 450.) then
! iterate to get consistently pH, T, mu and gamma
    tol       = 1.d-3
    converged = .false.
    isolve    = 0
    pH        = pH_tot ! arbitrary value, overwritten below, to avoid compiler warning
    T_old     = T
    i = 0
    do while (.not. converged .and. i < itermax)
       i = i+1
       KH2 = calc_KH2(T)
       if (present(pressure_cgs)) then
          pH = solve_q(1.,(1.+eps(iHe))/((1.+2.*eps(iHe))*KH2), &
                  -pressure_cgs/patm/((1.+2.*eps(iHe))*KH2))
          pH2    = KH2*pH**2
          pH_tot = pH+2.*pH2
          !rho_cgs = pH_tot*patm*mass_per_H/(T*kboltz)
       else
          pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H)
          pH        = solve_q(2.*KH2, 1., -pH_tot)
          pH2       = KH2*pH**2
       endif
       if (ieos /= 5) return

       mu    = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
       gamma = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
       if (i==1) then
          mu_old    = mu
          gamma_old = gamma
       else
          T_old     = T
          T         = T_old*mu*(gamma-1.)/(mu_old*(gamma_old-1.))
          converged = (abs(T-T_old)/T_old) < tol
          if (i>=itermax .and. .not.converged) then
             if (isolve==0) then
                isolve = isolve+1
                i      = 0
                tol    = 1.d-2
                print *,'[calc_muGamma_cond] cannot converge on T(mu,gamma). Trying with lower tolerance'
             else
                print *,'Told=',T_old,',T=',T,',gamma_old=',gamma_old,',gamma=',gamma,',mu_old=',&
                        mu_old,',mu=',mu,',dT/T=',abs(T-T_old)/T_old
                ! call fatal(label,'cannot converge on T(mu,gamma)')
             endif
          endif
       endif
    enddo
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2 = pH_tot/2.
    pH  = 0.
    if (ieos == 5) then
       mu    = (1.+4.*eps(iHe))/(0.5+eps(iHe))
       gamma = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
    endif
 endif
end subroutine calc_muGamma

pure elemental real function solve_q(a, b, c)
    ! Citardauq form: avoids cancellation. Assumes b > 0 and c < 0,
    ! so the denominator is always >= b > 0.
    real, intent(in) :: a, b, c
    solve_q = max(-2.*c / (b + sqrt(max(b*b - 4.*a*c, 0.))), tiny(solve_q))
end function solve_q

pure real function calc_KH2(T)
! all quantities are in cgs
 real, intent(in) :: T
 real, parameter :: R = 1.987165
 real :: G, d
 real, parameter :: coefs(5) = (/4.25321d+05, -1.07123d+05, 2.69980d+01, 5.48280d-04, -3.81498d-08/)

 G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
 d = min(-G/(R*T),222.)
 calc_KH2 = exp(d)
end function calc_KH2

!----------------------------------------------------------------
subroutine init_dust_formation
!----------------------------------------------------------------
!initialize dust formation arrays
 use dim,  only: do_nucleation,do_condensation,ndust_prop
 use part, only:npart,nucleation,n_nucleation,idmu,idgamma,&
         condensation,n_condensation,icmu,icgamma
 use eos,  only:gamma,gmw
 integer :: i
 real :: tmp_nucleation(n_nucleation),tmp_condensation(n_condensation)

 call set_abundances

 gamma = 5./3.
 !initialize nucleation array
 if (do_nucleation) then
    ndust_prop              = n_nucleation
    tmp_nucleation          = 0.
    tmp_nucleation(idmu)    = gmw
    tmp_nucleation(idgamma) = gamma
    do i=1,npart
       nucleation(:,i) = tmp_nucleation(:)
    enddo
 elseif (do_condensation) then
    ndust_prop                = n_condensation
    tmp_condensation          = 0.
    tmp_condensation(1:6)     = a_init_dust !initial dust radius (1nm)
    tmp_condensation(icmu)    = gmw
    tmp_condensation(icgamma) = gamma
    do i=1,npart
       condensation(:,i) = tmp_condensation(:)
    enddo
 else
    stop '[S-init_dust_formation] something went wrong ??'
 endif

end subroutine init_dust_formation

!------------------------------------
!
!  Bowen dust opacity formula
!
!------------------------------------
pure elemental real function calc_kappa_bowen(Teq)
!all quantities in cgs
 real, intent(in) :: Teq
 real :: dlnT

 dlnT = (Teq-bowen_Tcond)/bowen_delta
 if (dlnT > 50.) then
    calc_kappa_bowen = 0.
 else
    calc_kappa_bowen = bowen_kmax/(1.0 + exp(dlnT)) + kappa_gas
 endif

end function calc_kappa_bowen

!-----------------------------------------------------------------------
!
!  calculate dust opacity
!
!-----------------------------------------------------------------------
pure real function calc_kappa_dust(K3, Tdust, rho_cgs)
!all quantities in cgs
 use physcon, only:atomic_mass_unit
 real, intent(in) :: K3, Tdust, rho_cgs

 real :: kappa_cgs, fac
 real, parameter :: rho_Cdust = 2.62, mc = 12.*atomic_mass_unit
 real, parameter :: Qplanck_abs = 1.6846124267740528e+04
 real, parameter :: Qross_ext = 9473.2722815583311

 fac = max(0.75*K3*mc/(mass_per_H*rho_Cdust),1.e-15)
 !kappa_cgs = Qplanck_abs *fac ! planck
 !kappa_cgs = Qross_ext * fac  ! Rosseland

 ! Gail & Sedlmayr, 1985, A&A, 148, 183, eqs 23,24
 !kappa_cgs = 6.7 * fac * Tdust  ! planck
 kappa_cgs = 5.9 * fac * Tdust  ! Rosseland

 calc_kappa_dust = kappa_cgs + kappa_gas
end function calc_kappa_dust

!-----------------------------------------------------------------------
!
!  calculate alpha, reduced gravity factor
!
!-----------------------------------------------------------------------
pure real function calc_Eddington_factor(Mstar_cgs, Lstar_cgs, kappa_cgs, tau)
!all quantities in cgs
 use physcon, only:c,Gg,pi
 real, intent(in) :: Mstar_cgs,Lstar_cgs,kappa_cgs
 real, intent(in), optional :: tau

 if (present(tau)) then
    calc_Eddington_factor = Lstar_cgs*exp(-tau)/(4.*pi*c*Gg*Mstar_cgs) * kappa_cgs
 else
    calc_Eddington_factor = Lstar_cgs/(4.*pi*c*Gg*Mstar_cgs) * kappa_cgs
 endif
end function calc_Eddington_factor

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_dust_formation(hdr,ierr)
 use dump_utils,        only:dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr

! initial gas composition for dust formation
 call set_abundances
 call add_to_rheader(eps,'epsilon',hdr,ierr) ! array
 call add_to_rheader(Aw,'Amean',hdr,ierr)    ! array
 call add_to_rheader(mass_per_H,'mass_per_H',hdr,ierr) ! real
end subroutine write_headeropts_dust_formation

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_dust_formation(hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 real :: dum(nElements)

 ierr = 0
 call extract('mass_per_H',mass_per_H,hdr,ierr) ! real
 ! it is likely that your dump was generated with an old version of phantom
 ! and the chemical properties not stored. restore and save the default values
 if (mass_per_H < tiny(0.)) then
    print *,'reset dust chemical network properties'
    call set_abundances
    call extract('epsilon',dum(1:nElements),hdr,ierr) ! array
    call extract('Amean',dum(1:nElements),hdr,ierr) ! array
 else
    call extract('epsilon',eps(1:nElements),hdr,ierr) ! array
    call extract('Amean',Aw(1:nElements),hdr,ierr) ! array
 endif
end subroutine read_headeropts_dust_formation

!-----------------------------------------------------------------------
!+
!  Writes input options to the input file
!+
!-----------------------------------------------------------------------
subroutine write_options_dust_formation(iunit)
 use dim,          only:nucleation,condensation
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling dust'
 if (nucleation .or. condensation) then
    call write_inopt(idust_opacity,'idust_opacity','compute dust opacity (0=off, 1=bowen, 2=nucleation, 3=condensation)',iunit)
 else
    call write_inopt(idust_opacity,'idust_opacity','compute dust opacity (0=off, 1=bowen)',iunit)
 endif
 if (idust_opacity == 1) then
    call write_inopt(kappa_gas,'kappa_gas','constant gas opacity (cm²/g)',iunit)
    call write_inopt(bowen_kmax,'bowen_kmax','maximum dust opacity (cm²/g)',iunit)
    call write_inopt(bowen_Tcond,'bowen_Tcond','dust condensation temperature (K)',iunit)
    call write_inopt(bowen_delta,'bowen_delta','condensation temperature range (K)',iunit)
 endif
 if ((nucleation .or. condensation) .and. (idust_opacity == 2 .or. idust_opacity == 3)) then
    call write_inopt(kappa_gas,'kappa_gas','constant gas opacity (cm²/g)',iunit)
    call write_inopt(wind_CO_ratio ,'wind_CO_ratio','wind initial C/O ratio (> 1)',iunit)
 endif

end subroutine write_options_dust_formation

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_dust_formation(db,nerr)
 use io,      only:error
 use dim,     only:nucleation,do_nucleation,inucleation,condensation,&
         do_condensation,icondensation,store_dust_temperature,ndust_prop
 use part,    only:n_nucleation,n_condensation
 use eos,     only:ieos
 use infile_utils, only:inopts,read_inopt
 type(inopts), intent(inout) :: db(:)
 integer,      intent(inout) :: nerr

 call read_inopt(idust_opacity,'idust_opacity',db,errcount=nerr,min=0,max=3)
 ! bowen opacity
 if (idust_opacity == 1) then
    call read_inopt(kappa_gas,'kappa_gas',db,errcount=nerr,min=0.)
    call read_inopt(bowen_kmax,'bowen_kmax',db,errcount=nerr,min=0.)
    call read_inopt(bowen_Tcond,'bowen_Tcond',db,errcount=nerr,min=0.)
    call read_inopt(bowen_delta,'bowen_delta',db,errcount=nerr,min=0.)
 endif
 ! dust nucleation
 if (idust_opacity == 2) then
    do_nucleation = .true.
    do_condensation = .false.
    inucleation = 1
    ndust_prop = n_nucleation
 endif
 ! dust condensation
 if (idust_opacity == 3) then
    do_condensation = .true.
    do_nucleation = .false.
    icondensation = 1
    ndust_prop = n_condensation
 endif
 if ((nucleation .or. condensation) .and. (idust_opacity == 2 .or. idust_opacity == 3)) then
    call read_inopt(kappa_gas,'kappa_gas',db,errcount=nerr,min=0.)
    call read_inopt(wind_CO_ratio,'wind_CO_ratio',db,errcount=nerr,min=0.)
 endif
 if (idust_opacity > 0) store_dust_temperature = .true.
 if ((do_nucleation .or. do_condensation) .and. ieos == 5) call error('read_infile','with nucleation you must use ieos = 2')

end subroutine read_options_dust_formation

end module dust_formation
