!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
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
! :Dependencies: dim, dump_utils, eos, infile_utils, io, part, physcon,
!   units
!

 implicit none
 integer, public :: idust_opacity = 0

 public :: set_abundances,calc_kappa_dust,calc_kappa_bowen,&
      read_options_dust_formation,write_options_dust_formation,&
      calc_Eddington_factor,calc_muGamma,init_muGamma,init_nucleation,&
      write_headeropts_dust_formation,read_headeropts_dust_formation
!
!--runtime settings for this module
!
 real, public :: kappa_gas   = 2.d-4

 private

 character(len=*), parameter :: label = 'dust_formation'
 real :: wind_CO_ratio = 2.
 real :: bowen_kmax  = 2.7991
 real :: bowen_Tcond = 1500.
 real :: bowen_delta = 60.

contains

subroutine set_abundances
  use dim, only: do_nucleation,do_condensation
  use chemistry_moments, only:set_abundances_moments
  use chemistry_condensation, only:set_abundances_condensation

  if (do_nucleation) then
     call set_abundances_moments(wind_CO_ratio)
  elseif (do_condensation) then
     call set_abundances_condensation(wind_CO_ratio)
  endif
end subroutine set_abundances

subroutine init_muGamma(rho_cgs, T, mu, gamma, pH_tot, pH, pH2)
! all quantities are in cgs
 use dim, only: do_nucleation,do_condensation
 use chemistry_moments, only:init_muGamma_moments
 use chemistry_condensation, only:init_muGamma_condensation
 real, intent(in)              :: rho_cgs
 real, intent(inout)           :: T, mu, gamma
 real, intent(out), optional   :: pH_tot, pH,pH2

 if (do_nucleation) then
    call init_muGamma_moments(rho_cgs, T, mu, gamma)
 elseif (do_condensation) then
    call init_muGamma_condensation(rho_cgs, T, mu, gamma, pH_tot, pH, pH2)
 endif
end subroutine init_muGamma

subroutine calc_muGamma(rho_cgs, T, mu, gamma, pH, pH_tot, pH2, pressure)
! all quantities are in cgs
 use dim, only: do_nucleation,do_condensation
 use chemistry_moments, only:calc_muGamma_moments
 use chemistry_condensation, only:calc_muGamma_condensation
 real, intent(in)              :: rho_cgs
 real, intent(inout)           :: T, mu, gamma
 real, intent(out)             :: pH, pH_tot
 real, intent(out), optional   :: pH2, pressure

 if (do_nucleation) then
    call calc_muGamma_moments(rho_cgs, T, mu, gamma, pH, pH_tot)
 elseif (do_condensation) then
    call calc_muGamma_condensation(rho_cgs, T, mu, gamma, pH, pH_tot, pH2, pressure)
 endif
end subroutine calc_muGamma

subroutine init_nucleation
 use part,  only:npart,nucleation,n_nucleation,idmu,idgamma
 use chemistry_moments, only :set_abundances_moments
 use eos,   only:gamma,gmw
 integer :: i
 real :: JKmuS(n_nucleation)

 call set_abundances_moments(wind_CO_ratio)

 !initialize nucleation array
 gamma = 5./3.
 JKmuS = 0.
 JKmuS(idmu)    = gmw
 JKmuS(idgamma) = gamma
 do i=1,npart
    nucleation(:,i) = JKmuS(:)
 enddo

end subroutine init_nucleation

!------------------------------------
!
!  Bowen dust opacity formula
!
!------------------------------------
pure elemental real function calc_kappa_bowen(Teq)
!all quantities in cgs
 real,    intent(in)  :: Teq
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
 use chemistry_moments, only: mass_per_H
 real, intent(in) :: K3, Tdust, rho_cgs

 real :: kappa_cgs, fac
 real, parameter :: rho_Cdust = 2.62, mc = 12.*atomic_mass_unit
 real, parameter :: Qplanck_abs = 1.6846124267740528e+04
 real, parameter :: Qross_ext = 9473.2722815583311

 fac = max(0.75*K3*mc/(mass_per_H*rho_Cdust),1.e-15)
 !kappa_cgs = Qplanck_abs *fac ! planck
 !kappa_cgs = Qross_ext * fac  ! Rosseland

 ! Gail & Sedlmayr, 1985, A&A, 148, 183, eqs 23,24
 !kappa_cgs = 6.7d0 * fac * Tdust  ! planck
 kappa_cgs = 5.9d0 * fac * Tdust  ! Rosseland

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
 use chemistry_moments, only:Aw,mass_per_H,eps
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
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_dust_formation(hdr,ierr)
 use dump_utils, only:dump_h,extract
 use chemistry_moments, only:Aw,mass_per_H,eps,nelements
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
 use dim,          only:nucleation
 use infile_utils, only:write_inopt
 integer, intent(in) :: iunit

 write(iunit,"(/,a)") '# options controlling dust'
 if (nucleation) then
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
 if (nucleation .and. (idust_opacity == 2 .or. idust_opacity == 3)) then
    call write_inopt(kappa_gas,'kappa_gas','constant gas opacity (cm²/g)',iunit)
    call write_inopt(wind_CO_ratio ,'wind_CO_ratio','wind initial C/O ratio (> 1)',iunit)
 endif

end subroutine write_options_dust_formation

!-----------------------------------------------------------------------
!+
!  Reads input options from the input file.
!+
!-----------------------------------------------------------------------
subroutine read_options_dust_formation(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use dim,     only:do_nucleation,inucleation,do_condensation,icondensation,store_dust_temperature
 character(len=*), intent(in)  :: name,valstring
 logical, intent(out) :: imatch,igotall
 integer,intent(out) :: ierr

 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_nucleation'

 imatch  = .true.
 igotall = .false.
 select case(trim(name))
 case('idust_opacity')
    read(valstring,*,iostat=ierr) idust_opacity
    ngot = ngot + 1
    if (idust_opacity == 2) then
       do_nucleation = .true.
       inucleation = 1
    endif
    if (idust_opacity == 3) then
       do_condensation = .true.
       icondensation = 1
    endif
 case('wind_CO_ratio')
    read(valstring,*,iostat=ierr) wind_CO_ratio
    ngot = ngot + 1
    if (wind_CO_ratio < 0.) call fatal(label,'invalid setting for wind_CO_ratio (must be > 0)')
    if (wind_CO_ratio > 0.9 .or. wind_CO_ratio < 1.1) call fatal(label,'wind_CO_ratio must be < 0.9 or > 1.1)')
 case('kappa_gas')
    read(valstring,*,iostat=ierr) kappa_gas
    ngot = ngot + 1
    if (kappa_gas < 0.)    call fatal(label,'invalid setting for kappa_gas (<0)')
    !kgas = kappa_gas / (udist**2/umass)
 case('bowen_kmax')
    read(valstring,*,iostat=ierr) bowen_kmax
    ngot = ngot + 1
    if (bowen_kmax < 0.)    call fatal(label,'invalid setting for bowen_kmax (<0)')
 case('bowen_Tcond')
    read(valstring,*,iostat=ierr) bowen_Tcond
    ngot = ngot + 1
    if (bowen_Tcond < 0.) call fatal(label,'invalid setting for bowen_Tcond (<0)')
 case('bowen_delta')
    read(valstring,*,iostat=ierr) bowen_delta
    ngot = ngot + 1
    if (bowen_delta < 0.) call fatal(label,'invalid setting for bowen_delta (<0)')
 case default
    imatch = .false.
 end select
 igotall = (ngot >= 1)
 if (idust_opacity == 1) igotall = (ngot >= 5)
 if (idust_opacity == 2) igotall = (ngot >= 3)
 if (idust_opacity > 0) store_dust_temperature = .true.

end subroutine read_options_dust_formation

end module dust_formation
