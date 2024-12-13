!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module chemistry_moments
!
! Dust formation using the theory of moments
!
! :References: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!              Siess +, 2202, A&A
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

 implicit none

 public :: chemical_equilibrium_light,psat_C,set_abundances_moments,&
      init_muGamma_moments,calc_muGamma_moments

 !number of elements considered in the nucleation chemical network
 integer, public, parameter :: nElements = 10
 real, public :: mass_per_H, eps(nElements)
 real, public :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867]
 real, public, parameter :: Scrit = 2. ! Critical saturation ratio

 private
 integer, parameter :: nMolecules = 25
 real, parameter :: coefs(5,nMolecules) = reshape([&
       4.25321d+05, -1.07123d+05, 2.69980d+01, 5.48280d-04, -3.81498d-08, & !H2-
       4.15670d+05, -1.05260d+05, 2.54985d+01, 4.78020d-04, -2.82416d-08, & !OH-
       8.66184d+05, -2.27851d+05, 5.61473d+01, 7.62548d-04, -4.95254d-08, & !H2O
       3.30340d+05, -2.59792d+05, 3.23662d+01, 3.33339d-04, -1.69521d-08, & !CO
       5.34072d+05, -3.88536d+05, 6.95497d+01, 0.00000d+00,  0.00000d+00, & !CO2
       1.51591d+06, -4.09711d+05, 1.19952d+02, 5.71281d-04, -4.77554d-08, & !CH4
       6.31402d+05, -2.85395d+05, 5.97427d+01, 1.13748d-04, -2.04586d-08, & !C2H
       8.11613d+05, -3.99041d+05, 9.15250d+01, 2.53889d-05, -2.03949d-08, & !C2H2
       3.80144d+05, -2.28698d+05, 3.09554d+01, 2.71276d-04, -6.21348d-09, & !N2
       1.23626d+06, -2.89588d+05, 8.52686d+01, 5.66638d-04, -3.66862d-08, & !NH3
       4.51527d+05, -1.83350d+05, 2.97771d+01, 1.58267d-04, -1.72510d-08, & !CN
       6.02063d+05, -3.08724d+05, 6.00139d+01, 1.66592d-04, -1.13252d-08, & !HCN
      -1.10204d+05, -7.41179d+04, 2.57470d+01, 1.37542d-04, -2.83506d-09, & !Si2
       6.06066d+04, -1.71746d+05, 5.75688d+01, 1.18547d-04, -6.23584d-08, & !Si3
       1.82780d+05, -1.92839d+05, 3.03804d+01, 4.16079d-04, -1.91035d-08, & !SiO
       2.34924d+05, -2.60546d+05, 6.29047d+01, 2.75378d-04, -3.64630d-08, & !Si2C
       9.83348d+05, -3.16438d+05, 1.13885d+02, 2.90172d-04, -2.37569d-08, & !SiH4
       2.24963d+05, -1.03636d+05, 2.85814d+01, 1.77872d-04, -8.41752d-09, & !S2
       3.65656d+05, -8.77172d-04, 2.41557d+01, 3.40917d-04, -2.09217d-08, & !HS
       7.41911d+05, -1.81063d+05, 5.35114d+01, 4.94594d-04, -3.39231d-08, & !H2S
       1.44507d+05, -1.49916d+05, 2.87381d+01, 3.96186d-04, -2.18002d-08, & !SiS
       2.45408d+05, -7.17752d+04, 2.29449d+01, 4.52999d-04, -2.69915d-08, & !SiH
      -4.38897d+05, -1.58111d+05, 2.49224d+01, 1.08714d-03, -5.62504d-08, & !TiO
      -3.32351d+05, -3.04694d+05, 5.86984d+01, 1.17096d-03, -5.06729d-08, & !TiO2
       2.26786d+05, -1.43775d+05, 2.92429d+01, 1.69434d-04, -1.79867d-08], shape(coefs)) !C2
 integer, parameter, public :: iH = 1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10
 integer, parameter :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, &
      iNH3=10, iCN=11, iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, &
      iHS=19, iH2S=20, iSiS=21, iSiH=22, iTiO=23, iTiO2=24,iC2 = 25, iTiS=26

 character(len=*), parameter :: label = 'dust_formation'

contains


subroutine set_abundances_moments(wind_CO_ratio)
! all quantities in cgs
 use physcon, only:atomic_mass_unit
 real, intent(in) :: wind_CO_ratio
 eps(iH)  = 1.0
 eps(iHe) = 1.04d-1
 eps(iOx) = 6.d-4
 eps(iN)  = 2.52d-4
 eps(iNe) = 1.17d-4
 eps(iSi) = 3.58d-5
 eps(iS)  = 1.85d-5
 eps(iFe) = 3.24d-5
 eps(iTi) = 8.6d-8
 eps(iC)  = eps(iOx) * wind_CO_ratio
 mass_per_H = atomic_mass_unit*dot_product(Aw,eps)
 !XH  = atomic_mass_unit*eps(iH)/mass_per_H  ! H mass fraction
 !XHe = atomic_mass_unit*eps(iHe)/mass_per_H ! He mass fraction
end subroutine set_abundances_moments


!---------------------------------------------------------------
!
!  Compute carbon chemical equilibrum abundance in the gas phase
!
!---------------------------------------------------------------
subroutine chemical_equilibrium_light(rho_cgs, T, epsC, pC, pC2, pC2H, pC2H2, mu, gamma,&
     nH, nH2, nHe, nCO, nH2O, nOH)
! all quantities are in cgs
 use physcon, only:patm,mass_proton_cgs,kboltz
 real, intent(in)    :: rho_cgs, epsC
 real, intent(inout) :: T, mu, gamma
 real, intent(out)   :: pC, pC2, pC2H, pC2H2
 real, intent(out), optional :: nH, nH2, nHe, nCO, nH2O, nOH
 real    :: pH_tot, Kd(nMolecules+1), err, a, b, c, d
 real    :: pH, pCO, pO, pSi, pS, pTi, pN
 real    :: pC_old, pO_old, pSi_old, pS_old, pTi_old, cst
 integer :: i, nit

 call calc_muGamma_moments(rho_cgs, T, mu, gamma, pH, pH_tot)
 if (T > 1.d4) then
    pC    = eps(iC)*pH_tot
    pC2   = 0.
    pC2H  = 0.
    pC2H2 = 0.
    if (present(nH)) then
       nH   = pH  *(patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
       nH2  = 1.d-50
       nHe  = eps(ihe)*pH_tot* (patm*mass_per_H)/(mu*mass_proton_cgs*kboltz*T)
       nCO  = 1.d-50
       nH2O = 1.d-50
       nOH  = 1.d-50
    endif
    return
 endif

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)
 pCO      = epsC*pH_tot

 pC       = 0.
 pC_old   = 0.
 pO       = 0.
 pO_old   = 0.
 pSi      = 0.
 pSi_old  = 0.
 pS       = 0.
 pS_old   = 0.
 pTi      = 0.
 pTi_old  = 0.
 err      = 1.
 nit      = 0

 do while (err > 1.e-5)
! N
    pN  = solve_q(2.*Kd(iN2), &
            1.+pH**3*Kd(iNH3)+pC*(Kd(iCN)+Kd(iHCN)*pH), &
            -eps(iN)*pH_tot)
! C
    pC  = solve_q(2.*(Kd(iC2H)*pH + Kd(iC2H2)*pH**2 + Kd(iC2)), &
            1.+pO**2*Kd(iCO2)+pH**4*Kd(iCH4)+pSi**2*Kd(iSi2C), &
            pCO-epsC*pH_tot)
! O
    pO  = eps(iOx)*pH_tot/(1.+pC*Kd(iCO)+pH*Kd(iOH)+pH**2*Kd(iH2O)+2.*pO*pC*Kd(iCO2)+pSi*Kd(iSiO))
    pCO = Kd(iCO)*pC*pO
! Si
    a   = 3.*Kd(iSi3)
    b   = 2.*Kd(iSi2)+2.*Kd(iSi2C)*pC
    c   = 1.+Kd(iSiH)*pH+Kd(iSiH4)*pH**4+Kd(iSiO)*pO+Kd(iSiS)*pS
    d   = -eps(iSi)*pH_tot
    pSi = -d/(a*pSi**2+b*pSi+c)
! S
    pS  = solve_q(2.*Kd(iS2), &
            1.+Kd(iSiS)*pSi+Kd(iHS)*pH+Kd(iH2S)*pH**2, &
            -eps(iS)*pH_tot)
! Ti
    pTi = eps(iTi)*pH_tot/(1.+Kd(iTiO)*pO+Kd(iTiO2)*pO**2+Kd(iTiS)*pS)

    err = 0.
    if (pC  > 1.e-50) err = err + abs((pC-pC_old)/pC)
    if (pO  > 1.e-50) err = err + abs((pO-pO_old)/pO)
    if (pSi > 1.e-50) err = err + abs((pSi-pSi_old)/pSi)
    if (pS  > 1.e-50) err = err + abs((pS-pS_old)/pS)
    if (pTi > 1.e-50) err = err + abs((pTi-pTi_old)/pTi)

    nit = nit + 1
    if (nit == 200) exit

    pC_old  = pC
    pO_old  = pO
    pSi_old = pSi
    pS_old  = pS
    pTi_old = pTi
 enddo

 pC2   = Kd(iC2)*pC**2
 pC2H  = Kd(iC2H)*pC**2*pH
 pC2H2 = Kd(iC2H2)*pC**2*pH**2

! Convert partial pressures from atm to cgs
 pC    = pC*patm
 pC2   = pC2*patm
 pC2H  = pC2H*patm
 pC2H2 = pC2H2*patm
 if (present(nH)) then
    cst  = mass_per_H/(mu*mass_proton_cgs*kboltz*T)
    if (T < 450.) then
       nH2 = pH_tot/2.      *patm*cst
       nH  = 1.d-99
    else
       nH   = pH            *patm*cst
       nH2  = Kd(iH2)*pH**2 *patm*cst
    endif
    nHe  = eps(ihe)*pH_tot  *patm*cst
    nCO  = Kd(iCO) *pC*pO   *patm*cst
    nH2O = Kd(iH2O)*pH**2*pO*patm*cst
    nOH  = Kd(iOH) *pH*pO   *patm*cst
 endif
end subroutine chemical_equilibrium_light


!----------------------------------------
!
!  Calculate mean molecular weight, gamma
!
!----------------------------------------
subroutine calc_muGamma_moments(rho_cgs, T, mu, gamma, pH, pH_tot)
! all quantities are in cgs
 use io,  only:fatal
 use eos, only:ieos
 use physcon, only:patm,kboltz

 real, intent(in)    :: rho_cgs
 real, intent(inout) :: T, mu, gamma
 real, intent(out)   :: pH, pH_tot
 real :: KH2, pH2, x
 real :: T_old, mu_old, gamma_old, tol
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma_moments'

 pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
 T_old = T
 if (T > 1.d4) then
    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    pH     = pH_tot
    if (ieos /= 17) gamma  = 5./3.
 elseif (T > 450.) then
! iterate to get consistently pH, T, mu and gamma
    tol       = 1.d-3
    converged = .false.
    isolve    = 0
    pH        = pH_tot ! initial value, overwritten below, to avoid compiler warning
    i = 0
    do while (.not. converged .and. i < itermax)
       i = i+1
       pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H)
       KH2       = calc_Kd(coefs(:,iH2), T)
       pH        = solve_q(2.*KH2, 1., -pH_tot)
       pH2       = KH2*pH**2
       mu        = (1.+4.*eps(iHe))/(.5+eps(iHe)+0.5*pH/pH_tot)
       if (ieos == 17) exit !only update mu, keep gamma constant
       x         = 2.*(1.+4.*eps(iHe))/mu
       gamma     = (3.*x+4.+4.*eps(iHe))/(x+4.+4.*eps(iHe))
       converged = (abs(T-T_old)/T_old) < tol
       if (i == 1) then
          mu_old = mu
          gamma_old = gamma
       else
          T = T_old*mu/mu_old/(gamma_old-1.)*2.*x/(x+4.+4.*eps(iHe))
          if (i>=itermax .and. .not.converged) then
             if (isolve==0) then
                isolve = isolve+1
                i      = 0
                tol    = 1.d-2
                print *,'[dust_formation] cannot converge on T(mu,gamma). Trying with lower tolerance'
             else
                print *,'Told=',T_old,',T=',T,',gamma_old=',gamma_old,',gamma=',gamma,',mu_old=',&
                  mu_old,',mu=',mu,',dT/T=',abs(T-T_old)/T_old,', rho=',rho_cgs
                call fatal(label,'cannot converge on T(mu,gamma)')
             endif
          endif
       endif
    enddo
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2    = pH_tot/2.
    pH     = 0.
    mu     = (1.+4.*eps(iHe))/(0.5+eps(iHe))
    if (ieos /= 17)  gamma  = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
 endif

end subroutine calc_muGamma_moments

!--------------------------------------------
!
!  Initialise mean molecular weight and gamma
!
!--------------------------------------------
subroutine init_muGamma_moments(rho_cgs, T, mu, gamma, ppH, ppH2)
! all quantities are in cgs
 use physcon, only:patm,kboltz
 real, intent(in)              :: rho_cgs
 real, intent(inout)           :: T
 real, intent(out)             :: mu, gamma
 real, intent(out), optional   :: ppH, ppH2
 real :: KH2, pH_tot, pH, pH2

 pH_tot = rho_cgs*kboltz*T/(patm*mass_per_H)
 if (T > 1.d5) then
    pH2 = 0.
    pH  = pH_tot
 elseif (T > 450.) then
    KH2 = calc_Kd(coefs(:,iH2), T)
    pH  = solve_q(2.*KH2, 1., -pH_tot)
    pH2 = KH2*pH**2
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules
    pH2 = pH_tot/2.
    pH  = 0.
 endif
 mu    = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
 gamma = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
 call calc_muGamma_moments(rho_cgs, T, mu, gamma, pH, pH_tot)
 if (present(ppH))  ppH = pH
 if (present(ppH2)) ppH2 = pH2

end subroutine init_muGamma_moments


!-------------------------------------------------------------------------
!
!  Compute saturation pressure of carbon clusters C_1, ..., C_5 over graphite
!
!-------------------------------------------------------------------------
pure real function psat_C(T)
! all quantities are in cgs
 real, intent(in) :: T

 real, parameter :: f = 13.8287
 real :: T2,T3,pC1!,pC2,pC3,pC4,pC5

 if (T > 1.d4) then
    Psat_C = huge(Psat_C)
 else
    T2 = T*T
    T3 = T*T*T
    pC1 = exp(-8.61240d4/T + 1.85106d1 + 5.87980d-4*T - 2.51549d-7*T2 + 3.24892d-11*T3 + f)
    !pC2 = exp(-1.01124d5/T + 2.35611d1 + 3.37807d-4*T - 2.94959d-7*T2 + 4.41801D-11*T3 + f)
    !pC3 = exp(-9.90261d4/T + 2.81161d1 - 1.55820d-3*T + 1.60002d-7*T2 - 4.47171D-12*T3 + f)
    !pC4 = exp(-1.17037d5/T + 2.55579d1 - 5.63869d-6*T - 2.13596d-7*T2 + 3.39660D-11*T3 + f)
    !pC5 = exp(-1.18080d5/T + 2.65798d1 + 1.20285d-4*T - 2.68583d-7*T2 + 4.12365D-11*T3 + f)
    psat_C = pC1
 endif
end function psat_C

!-----------------------------
!
!  solve 2nd order polynomial
!
!-----------------------------
pure real function solve_q(a, b, c)
 real, intent(in) :: a, b, c
 real :: delta
 if (-4.*a*c/b**2 > epsilon(0.)) then
    delta = max(b**2-4.*a*c, 0.)
    solve_q = (-b+sqrt(delta))/(2.*a)
 else
    solve_q = -c/b
 endif
 solve_q = max(solve_q,1e-50)
end function solve_q

!------------------------------------
!
!  Compute dissociation coefficients
!
!------------------------------------
pure real function calc_Kd(coefs, T)
! all quantities are in cgs
 real, intent(in) :: coefs(5), T
 real, parameter :: R = 1.987165
 real :: G, d
 G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
 d = min(-G/(R*T),222.)
 calc_Kd = exp(d)
end function calc_Kd

pure real function calc_Kd_TiS(T)
! all quantities are in cgs
 use physcon, only:patm
 real, intent(in) :: T
 real, parameter :: a = 1.3316d1, b = -6.2216, c = 4.5829d-1, d = -6.4903d-2, e = 3.2788d-3
 real :: theta, logKd
 theta = 5040./T
 logKd = a+(b+(c+(d+e*theta)*theta)*theta)*theta
 calc_Kd_TiS = 10.**(-logKd)*patm
end function calc_Kd_TiS

end module chemistry_moments
