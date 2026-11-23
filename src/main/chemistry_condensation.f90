!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module chemistry_condensation
!
! chemistry_condensation
!
! :References: None
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: dump_utils, physcon
!

 use dust_formation, only: iH,iHe,iC,iOx,iN,iNe,iSi,iS,iFe,iTi,iMg,nElements
 implicit none
 public :: network

 !      real :: wind_CO_ratio = 0.34 !Based on Ferrarotti and Gail 2002

 integer, parameter, private :: nMolecules = 70

! Indices for molecules
 integer, parameter, public :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, iNH3=10, iCN=11, &
      iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, iHS=19, iH2S=20, iSiS=21, &
      iSiH=22, iTiO=23, iTiO2=24, iC2=25, iO2=26, iCH=27,iCOH=28, iC2O=29, iCH2=30, iH2CO=31, iCH3=32, &
      iC2H4=33, iNH=34, iNO=35, iNCO=36, iHCNO=37, iC2N=38, iC2N2=39, iHNO=40, iHNO2=41, iHNO3=42, &
      iNH2=43, iNO2=44, INO3=45, iN2O=46, iN2O4=47, iMgH=48, iMgOH=49, iMgO2H2=50, iMgN=51, iMgO=52, &
      iSiC=53, iSiN=54, iSiO2=55, iFeO=56, iFeO2H2=57, iCOS=58, iCS=59, iCS2=60, iFeS=61, iH2SO4=62, &
      iMgS=63, iSN=64, iSO=65, iSO2=66, iSO3=67, iSiH2=68, iSiH3=69, iTiS=70
! Indices for atomic abundances
 integer, parameter, public :: iaH=nMolecules+1, iaHe=nMolecules+2, iaC=nMolecules+3, iaO=nMolecules+4, iaN=nMolecules+5,&
      iaSi=nMolecules+6, iaS=nMolecules+7, iaFe=nMolecules+8, iaMg=nMolecules+9, iaTi=nMolecules+10

 private
 !integer, parameter :: nMolecules = 70
 integer, parameter :: nKd_sharp = 67
 integer, parameter :: nKd_tsuji = 3
 real, parameter :: exp_lim = 700.,log_lim = exp_lim/2.71828182846

 !Coefficients can be found in Sharp, C. M. ; Huebner, W. F. 1990, ApJS, 72, 417
 real, parameter :: coefs(5,nKd_sharp) = reshape([&
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
      3.65656d+05, -8.77172d+04, 2.41557d+01, 3.40917d-04, -2.09217d-08, & !HS
      7.41911d+05, -1.81063d+05, 5.35114d+01, 4.94594d-04, -3.39231d-08, & !H2S
      1.44507d+05, -1.49916d+05, 2.87381d+01, 3.96186d-04, -2.18002d-08, & !SiS
      2.45408d+05, -7.17752d+04, 2.29449d+01, 4.52999d-04, -2.69915d-08, & !SiH
      -4.38897d+05, -1.58111d+05, 2.49224d+01, 1.08714d-03, -5.62504d-08, & !TiO
      -3.32351d+05, -3.04694d+05, 5.86984d+01, 1.17096d-03, -5.06729d-08, & !TiO2
      2.26786d+05, -1.43775d+05, 2.92429d+01, 1.69434d-04, -1.79867d-08, & !C2
      3.17263d+05, -1.21503d+05, 3.11753d+01, 2.53855d-04, -1.73743d-08, & !O2
      4.08366d+05, -8.44491d+04, 2.54152d+01, 1.07036d-04, -1.04264d-08, & !CH
      6.07224d+05, -2.77376d+05, 5.66759d+01, 4.80619d-04, -2.82627d-08, & !COH
      3.56583d+05, -3.36691d+05, 6.30120d+01, -2.15913d-04, 1.12031d-08, & !C2O
      6.99399d+05, -1.88493d+05, 5.33630d+01, 5.15801d-04, -3.28189d-08, & !CH2
      9.91212d+05, -3.70812d+05, 9.03305d+01, 4.39620d-04, -3.01383d-08, & !H2CO
      1.06713d+06, -3.01022d+05, 8.47773d+01, 6.14095d-04, -4.34934d-08, & !CH3
      1.57756d+06, -5.51401d+05, 1.51644d+02, 2.93007d-04, -3.54396d-08, & !C2H4
      4.11181d+05, -7.79746d+04, 2.42428d+01, 4.09611d-04, -2.55866d-08, & !NH
      3.18598d+05, -1.53355d+05, 2.79315d+01, 2.76970d-04, -9.14848d-09, & !NO
      3.22530d+05, -3.08692d+05, 6.20446d+01, 5.46329d-05,  0.00000d+00, & !NCO
      7.64992d+05, -4.26302d+05, 9.19995d+01, 2.86974d-04, -1.52951d-08, & !HCNO
      3.05118d+05, -3.24444d+05, 5.98669d+01, 0.00000d+00,  0.00000d+00, & !C2N
      1.77312d+05, -4.97555d+05, 9.62430d+01, -4.71719d-05, 0.00000d+00, & !C2N2
      6.24780d+05, -2.05700d+05, 5.64594d+01, 4.64490d-04, -2.07625d-08, & !HNO
      6.98937d+05, -3.08668d+05, 8.96184d+01, 2.67463d-04, -1.22202d-08, & !HNO2 (trans)
      8.02373d+05, -3.82689d+05, 1.26355d+02, 1.65050d-05, -2.02960d-10, & !HNO3
      8.79993d+05, -1.77983d+05, 5.28632d+01, 6.37403d-04, -4.63748d-08, & !NH2
      4.70857d+05, -2.28020d+05, 6.18243d+01, 3.24505d-04, -9.54732d-09, & !NO2
      4.42258d+05, -2.78648d+05, 9.82856d+01, 6.50496d-06,  7.35167d-09, & !NO3
      4.31224d+05, -2.69421d+05, 6.44178d+01, 1.28733d-05,  7.77099d-09, & !N2O
      5.92417d+05, -4.67573d+05, 1.63906d+02, -5.25057d-04, 3.57897d-08, & !N2O4
      2.90391d+05, -4.91449d+04, 1.96635d+01, 1.52956d-04, -6.37175d-09, & !MgH
      4.47672d+05, -1.89626d+05, 5.17801d+01, 1.15404d-04,  3.34857d-10, & !MgOH
      7.28023d+05, -4.00948d+05, 1.11212d+02, 5.21735d-05,  1.76309d-11, & !MgO2H2
      2.50267d+05, -8.09495d+04, 2.05587d+01, 4.19383d-05,  1.08947d-08, & !MgN
      -3.68825d+05, -7.91991d+04, 2.24741d+01, -3.33176d-04, 3.46191d-08, & !MgO
      8.92245d+03, -1.06824d+05, 2.65205d+01, 1.67830d-04, -7.50467d-09, & !SiC
      3.30323d+05, -1.34105d+05, 2.82397d+01, -4.05990d-05, 5.40508d-09, & !SiN
      2.61948d+05, -3.02036d+05, 6.59588d+01, 1.07405d-04,  0.00000d+00, & !SiO2
      2.18545d+05, -1.00927d+05, 2.69229d+01, 2.53772d-04,  4.37490d-09, & !FeO
      4.18570d+05, -4.05885d+05, 1.14517d+02, -2.48115d-04, 3.53182d-08, & !FeO2H2
      3.84946d+05, -3.33394d+05, 6.59981d+01, -6.67359d-05, 0.00000d+00, & !COS
      2.54293d+05, -1.72618d+05, 3.05596d+01, 2.67035d-04, -1.59172d-08, & !CS
      2.41613d+05, -2.78285d+05, 6.52867d+01, -1.07214d-04, 0.00000d+00, & !CS2
      4.56757d+04, -7.77456d+04, 2.39279d+01, 3.02718d-04,  6.29763d-09, & !FeS
      9.36433d+05, -5.92570d+05, 1.90572d+02, -2.99786d-04, 1.02593d-08, & !H2SO4
      -3.64897d+05, -6.47278d+04, 2.04150d+01, -1.93425d-04, 2.73835d-08, & !MgS
      2.71848d+05, -1.18293d+05, 2.65189d+01, 1.88122d-04, -5.20791d-09, & !SN
      3.21629d+05, -1.27035d+05, 2.89402d+01, 1.25603d-04, -9.85969d-09, & !SO
      4.16261d+05, -2.59787d+05, 6.30161d+01, 2.30250d-04, -1.21022d-08, & !SO2
      4.54902d+05, -3.43672d+05, 1.00999d+02, 5.09991d-05,  0.00000d+00], shape(coefs)) !SO3

 !Coefficients from Sharp & Huebner, 1990 (ApJS 72, 417), Tsuji (1973) fit formulation
 real, parameter :: coefs_tsuji(6,nKd_tsuji) = reshape([&
      2.46640d+01, -6.88730d+00, 8.37100d-02, -1.00580d-02, 4.92910d-04, 3., & !SiH2
      3.63290d+01, -1.05560d+01, 8.09450d-02, -8.62120d-03, 3.98640d-04, 4., & !SiH3
      1.33160d+01, -6.22160d+00, 4.58290d-01, -6.49030d-02, 3.2788d0-03, 2.], shape(coefs_tsuji)) !TiS
 integer, public, parameter :: ncols = nElements + nMolecules - 1 ! -1 because excluding Neon

contains

!-----------------------------------------------------------------------
!+
!  driver to solve the chemical network
!+
!-----------------------------------------------------------------------
subroutine network(T,rho_cgs,mu,gamma,wind_CO_ratio,pH_tot,fol,fpy,fqu,fir,fsc,fcarb,abundance,pressure_cgs)
! all pressures are in bar
 use physcon, only:patm,kboltz
 use dust_formation, only:calc_muGamma,eps
 real, intent(in)     :: wind_CO_ratio, fol, fpy, fqu, fir, fsc, fcarb
 real, intent(out)    :: mu,gamma,pH_tot
 real, intent(inout)  :: T, rho_cgs
 real, intent(out)    :: abundance(ncols)
 real, intent(in), optional :: pressure_cgs
 integer :: j
 real, dimension(nMolecules)    :: pmol
 real                           :: pelm(nElements),fac

 pelm = 0.
 pmol = 0.

 !compute pH_tot, pH, pH2 and update mu and gamma if ieos = 5
 if (present(pressure_cgs)) then
    call calc_muGamma(rho_cgs,T,mu,gamma,pH_tot,pelm(iH),pmol(iH2),pressure_cgs)
 else
    call calc_muGamma(rho_cgs,T,mu,gamma,pH_tot,pelm(iH),pmol(iH2))
 endif

 pelm(iHe) = eps(iHe)*pH_tot
 call chemical_equilibrium(rho_cgs, T, pmol, pelm, mu, gamma, &
         wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb) !, nMolecules)

 !convert the pressure of the molecules in number density (cgs, hence the factor patm in fac)
 fac = patm/(kboltz*T)
 do j=1,nmolecules
    abundance(j) = max(1.d-50,pmol(j)*fac)
 enddo
 abundance(iaH)  = max(1.d-50,pelm(iH) *fac)
 abundance(iaHe) = max(1.d-50,pelm(iHe)*fac)
 abundance(iaC)  = max(1.d-50,pelm(iC) *fac)
 abundance(iaO)  = max(1.d-50,pelm(iOx)*fac)
 abundance(iaN)  = max(1.d-50,pelm(iN) *fac)
 abundance(iaSi) = max(1.d-50,pelm(iSi)*fac)
 abundance(iaS)  = max(1.d-50,pelm(iS) *fac)
 abundance(iaFe) = max(1.d-50,pelm(iFe)*fac)
 abundance(iaMg) = max(1.d-50,pelm(iMg)*fac)
 abundance(iaTi) = max(1.d-50,pelm(iTi)*fac)

 !abundance(iH2O) = abundance(iH2O) - (3.*fol + 2.*fpy + fqu) * eps(iSi)*fac *pH_tot
 !abundance(iSiO) = abundance(iSiO) - (fol + fpy + fqu) * eps(iSi)*fac * pH_tot
 !abundance(iaMg) = abundance(iaMg) - (2.*fol + fpy) * eps(iSi)*fac * pH_tot
 !if (abundance(iH2O) < 0.0) abundance(iH2O) = 0.0
 !if (abundance(iSiO) < 0.0) abundance(iSiO) = 0.0
 !if (abundance(iaMg) < 0.0) abundance(iaMg) = 0.0

 !!Update the abundance of Si and Oxygen.

end subroutine network

!---------------------------------------------------------------
!
!  Compute partial pressures assuming chemical equilibrum
!
!---------------------------------------------------------------
subroutine chemical_equilibrium(rho_cgs,T,pmol,pelm,mu, gamma, &
      wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb)
! all pressures are in bar
 use dust_formation, only:eps
 real, intent(inout) :: T, pmol(nMolecules), pelm(nElements), rho_cgs, mu, gamma
 real, intent(in)    :: wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb

 integer, parameter :: nit_max  = 100
 integer, parameter :: n_newton = 2     ! try 2 here once validated
 real,    parameter :: tol      = 1.e-6

 real    :: Kd(nKd_sharp+nKd_tsuji)
 real    :: xH, xC, xO, xN, xMg, xSi, xFe, xS, xTi          ! current iterate
 real    :: oC, oO, oN, oMg, oSi, oFe, oS, oTi              ! previous sweep
 real    :: errmax, cC_first, cC_main, eO_first_lt, eO_main, tn
 integer :: i, nit
 logical :: lt1

 if (T > 1.e4) then
    pelm = eps*pH_tot
    pmol = 0.
    return
 endif

! Dissociation constants
 do i=1,nKd_sharp
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 do i=1,nKd_tsuji
    Kd(nKd_sharp+i) = calc_Kd_tsuji(coefs_tsuji(:,i), T)
 enddo

! Loop-invariant constants
 tn          = tiny(1.)
 lt1         = (wind_CO_ratio < 1.0)
 cC_main     = (fcarb-1.0)*eps(iC)*pH_tot
 cC_first    = (eps(iOx)+(fcarb-1.0)*eps(iC))*pH_tot                 ! wind_CO >= 1, first sweep
 eO_first_lt = (eps(iC)+eps(iSi)*(1.+4.*fol+3.*fpy+2.*fqu)-eps(iOx))*pH_tot ! wind_CO < 1, first sweep
 eO_main     = (eps(iSi)*(4.*fol+3.*fpy+2.*fqu)-eps(iOx))*pH_tot

! Scalars (H is never updated)
 xH = pelm(iH);  xC = pelm(iC);  xO = pelm(iOx); xN = pelm(iN)
 xMg = pelm(iMg); xSi = pelm(iSi); xFe = pelm(iFe); xS = pelm(iS); xTi = pelm(iTi)

 do nit = 1, nit_max
    oC = xC; oO = xO; oN = xN; oMg = xMg; oSi = xSi; oFe = xFe; oS = xS; oTi = xTi

    if (lt1) then
       if (nit == 1) then
          call upd_O(.false., eO_first_lt, eps(iOx)*pH_tot, n_newton)
       else
          call upd_O(.true.,  eO_main,     xO,              n_newton)
       endif
       call upd_C(.true., cC_main)
    else
       if (nit == 1) then
          call upd_C(.false., cC_first)
          call upd_O(.true., eO_main, eps(iOx)*pH_tot, n_newton)
       else
          call upd_C(.true., cC_main)
          call upd_O(.true., eO_main, xO, n_newton)
       endif
    endif

    call upd_N()
    call upd_Mg()
    call upd_Si()
    call upd_Fe()
    call upd_S()
    call upd_Ti()

    errmax = max(abs(xC-oC)/(oC+tn),   abs(xO-oO)/(oO+tn),   abs(xN-oN)/(oN+tn),   &
                 abs(xMg-oMg)/(oMg+tn), abs(xSi-oSi)/(oSi+tn), abs(xFe-oFe)/(oFe+tn), &
                 abs(xS-oS)/(oS+tn),    abs(xTi-oTi)/(oTi+tn))
    if (errmax <= tol) exit
 enddo

 pelm(iC) = xC;  pelm(iOx) = xO; pelm(iN) = xN;  pelm(iMg) = xMg
 pelm(iSi) = xSi; pelm(iFe) = xFe; pelm(iS) = xS; pelm(iTi) = xTi

 call upd_pmol()

contains

 !------------------------------------------------------------------
 subroutine upd_O(withCO, e, zz, maxit)
  logical, intent(in) :: withCO
  real,    intent(in) :: e, zz
  integer, intent(in) :: maxit
  real :: xH2, xN2, xC2, a, b, c, d
  xH2 = xH*xH; xN2 = xN*xN; xC2 = xC*xC
  a = 4.*(xN2*Kd(iN2O4) + xH2*xS*Kd(iH2SO4))
  b = 3.*(xH*xN*Kd(iHNO3) + xN*Kd(iNO3) + xS*Kd(iSO3))
  c = 2.*(Kd(iO2) + xC*Kd(iCO2) + xH*xN*Kd(iHNO2) + xN*Kd(iNO2) &
        + xH2*xMg*Kd(iMgO2H2) + xSi*Kd(iSiO2) + xFe*xH2*Kd(iFeO2H2) &
        + xS*Kd(iSO2) + xTi*Kd(iTiO2))
  d = 1. + xH*Kd(iOH) + xH2*Kd(iH2O) + xC*xH*Kd(iCOH) + xC2*Kd(iC2O) &
        + xC*xH2*Kd(iH2CO) + xN*Kd(iNO) + xC*xN*Kd(iNCO) + xC*xH*xN*Kd(iHCNO) &
        + xH*xN*Kd(iHNO) + xN2*Kd(iN2O) + xH*xMg*Kd(iMgOH) + xMg*Kd(iMgO) &
        + xFe*Kd(iFeO) + xC*xS*Kd(iCOS) + xS*Kd(iSO) + xTi*Kd(iTiO)
  if (withCO) d = d + xC*Kd(iCO) + xSi*Kd(iSiO)
  xO = newton_method(a, b, c, d, e, zz, maxit)
 end subroutine upd_O

 !------------------------------------------------------------------
 subroutine upd_C(withCO, c)
  logical, intent(in) :: withCO
  real,    intent(in) :: c
  real :: xH2, xH3, xH4, xN2, xO2, xS2, xSi2, a, b
  xH2 = xH*xH; xH3 = xH2*xH; xH4 = xH2*xH2
  xN2 = xN*xN; xO2 = xO*xO; xS2 = xS*xS; xSi2 = xSi*xSi
  a = 2.*(Kd(iC2) + xH*Kd(iC2H) + xO*Kd(iC2O) + xH2*Kd(iC2H2) + xH4*Kd(iC2H4) &
        + xN*Kd(iC2N) + xN2*Kd(iC2N2))
  b = 1. + xH*Kd(iCH) + xO*xH*Kd(iCOH) + xO2*Kd(iCO2) + xH2*Kd(iCH2) &
        + xO*xH2*Kd(iH2CO) + xH3*Kd(iCH3) + xH4*Kd(iCH4) + xN*Kd(iCN) &
        + xN*xH*Kd(iHCN) + xO*xN*Kd(iNCO) + xO*xH*xN*Kd(iHCNO) &
        + xSi*Kd(iSiC) + xSi2*Kd(iSi2C) + xO*xS*Kd(iCOS) + xS*Kd(iCS) + xS2*Kd(iCS2)
  if (withCO) b = b + xO*Kd(iCO)
  xC = solve_q(a, b, c)
 end subroutine upd_C

 !------------------------------------------------------------------
 subroutine upd_N()
  real :: xH2, xH3, xC2, xO2, xO3, xO4, a, b
  xH2 = xH*xH; xH3 = xH2*xH; xC2 = xC*xC
  xO2 = xO*xO; xO3 = xO2*xO; xO4 = xO2*xO2
  a = 2.*(Kd(iN2) + xC2*Kd(iC2N2) + xO*Kd(iN2O) + xO4*Kd(iN2O4))
  b = 1. + xH*Kd(iNH) + xO*Kd(iNO) + xC*Kd(iCN) + xC*xH*Kd(iHCN) &
        + xO*xC*Kd(iNCO) + xC*xO*xH*Kd(iHCNO) + xC2*Kd(iC2N) &
        + xO*xH*Kd(iHNO) + xO2*xH*Kd(iHNO2) + xO3*xH*Kd(iHNO3) &
        + xH2*Kd(iNH2) + xH3*Kd(iNH3) + xO2*Kd(iNO2) + xO3*Kd(iNO3) &
        + xMg*Kd(iMgN) + xSi*Kd(iSiN) + xS*Kd(iSN)
  xN = solve_q(a, b, -eps(iN)*pH_tot)
 end subroutine upd_N

 !------------------------------------------------------------------
 subroutine upd_Mg()
  real :: xH2, xO2
  xH2 = xH*xH; xO2 = xO*xO
  xMg = (eps(iMg)-(2.*fol+fpy)*eps(iSi))*pH_tot / &
        (1. + xH*Kd(iMgH) + xO*xH*Kd(iMgOH) + xO2*xH2*Kd(iMgO2H2) &
            + xN*Kd(iMgN) + xO*Kd(iMgO) + xS*Kd(iMgS))
 end subroutine upd_Mg

 !------------------------------------------------------------------
 subroutine upd_Si()
  real :: xH2, xH3, xH4, xO2, a, b
  xH2 = xH*xH; xH3 = xH2*xH; xH4 = xH2*xH2; xO2 = xO*xO
  a = 2.*(Kd(iSi2) + xC*Kd(iSi2C))
  b = 1. + xC*Kd(iSiC) + xH*Kd(iSiH) + xH2*Kd(iSiH2) + xH3*Kd(iSiH3) &
        + xH4*Kd(iSiH4) + xN*Kd(iSiN) + xO*Kd(iSiO) + xO2*Kd(iSiO2) + xS*Kd(iSiS)
  xSi = solve_q(a, b, (fol+fpy+fqu+fsc-1.0)*eps(iSi)*pH_tot)
 end subroutine upd_Si

 !------------------------------------------------------------------
 subroutine upd_Fe()
  real :: xH2, xO2
  xH2 = xH*xH; xO2 = xO*xO
  xFe = (1.0-fir)*eps(iFe)*pH_tot / &
        (1. + xO*Kd(iFeO) + xO2*xH2*Kd(iFeO2H2) + xS*Kd(iFeS))
 end subroutine upd_Fe

 !------------------------------------------------------------------
 subroutine upd_S()
  real :: xH2, xO2, xO3, xO4, a, b
  xH2 = xH*xH; xO2 = xO*xO; xO3 = xO2*xO; xO4 = xO2*xO2
  a = 2.*(Kd(iS2) + xC*Kd(iCS2))
  b = 1. + xO*xC*Kd(iCOS) + xC*Kd(iCS) + xFe*Kd(iFeS) + xH*Kd(iHS) &
        + xH2*Kd(iH2S) + xO4*xH2*Kd(iH2SO4) + xMg*Kd(iMgS) + xN*Kd(iSN) &
        + xO*Kd(iSO) + xO2*Kd(iSO2) + xO3*Kd(iSO3) + xSi*Kd(iSiS) + xTi*Kd(iTiS)
  xS = solve_q(a, b, -eps(iS)*pH_tot)
 end subroutine upd_S

 !------------------------------------------------------------------
 subroutine upd_Ti()
  real :: xO2
  xO2 = xO*xO
  xTi = eps(iTi)*pH_tot / (1. + xO*Kd(iTiO) + xO2*Kd(iTiO2) + xS*Kd(iTiS))
 end subroutine upd_Ti

 !------------------------------------------------------------------
 subroutine upd_pmol()
  real :: H2,C2,O2
  H2 = xH*xH
  C2 = xC*xC
  O2 = xO*xO

  pmol(iH2)    = Kd(iH2)*H2
  pmol(iOH)    = Kd(iOH)*pelm(iOx)*pelm(iH)
  pmol(iH2O)   = Kd(iH2O)*pelm(iOx)*H2
  pmol(iCO)    = Kd(iCO)*pelm(iOx)*pelm(iC)
  pmol(iCO2)   = Kd(iCO2)*pelm(iC)*O2
  pmol(iCH4)   = Kd(iCH4)*pelm(iC)*pelm(iH)**4
  pmol(iC2H)   = Kd(iC2H)*pelm(iH)*C2
  pmol(iC2H2)  = Kd(iC2H2)*C2*H2
  pmol(iN2)    = Kd(iN2)*pelm(iN)**2
  pmol(iNH3)   = Kd(iNH3)*pelm(iN)*pelm(iH)**3
  pmol(iCN)    = Kd(iCN)*pelm(iC)*pelm(iN)
  pmol(iHCN)   = Kd(iHCN)*pelm(iH)*pelm(iC)*pelm(iN)
  pmol(iSi2)   = Kd(iSi2)*pelm(iSi)**2
  pmol(iSi3)   = Kd(iSi3)*pelm(iSi)**3
  pmol(iSiO)   = Kd(iSiO)*pelm(iSi)*pelm(iOx)
  pmol(iSi2C)  = Kd(iSi2C)*pelm(iSi)**2*pelm(iC)
  pmol(iSiH4)  = Kd(iSiH4)*pelm(iH)**4*pelm(iSi)
  pmol(iS2)    = Kd(iS2)*pelm(iS)**2
  pmol(iHS)    = Kd(iHS)*pelm(iH)*pelm(iS)
  pmol(iH2S)   = Kd(iH2S)*H2*pelm(iS)
  pmol(iSiS)   = Kd(iSiS)*pelm(iSi)*pelm(iS)
  pmol(iSiH)   = Kd(iSiH)*pelm(iSi)*pelm(iH)
  pmol(iC2)    = Kd(iC2)*C2
  pmol(iO2)    = Kd(iO2)*O2
  pmol(iCH)    = Kd(iCH)*pelm(iH)*pelm(iC)
  pmol(iCOH)   = Kd(iCOH)*pelm(iOx)*pelm(iC)*pelm(iH)
  pmol(iC2O)   = Kd(iC2O)*pelm(iOx)*C2
  pmol(iCH2)   = Kd(iCH2)*H2*pelm(iC)
  pmol(iH2CO)  = Kd(iH2CO)*pelm(iOx)*pelm(iC)*H2
  pmol(iCH3)   = Kd(iCH3)*pelm(iH)**3*pelm(iC)
  pmol(iC2H4)  = Kd(iC2H4)*pelm(iH)**4*C2
  pmol(iNH)    = Kd(iNH)*pelm(iN)*pelm(iH)
  pmol(iNO)    = Kd(iNO)*pelm(iOx)*pelm(iN)
  pmol(iNCO)   = Kd(iNCO)*pelm(iOx)*pelm(iC)*pelm(iN)
  pmol(iHCNO)  = Kd(iHCNO)*pelm(iOx)*pelm(iC)*pelm(iH)*pelm(iN)
  pmol(iC2N)   = Kd(iC2N)*pelm(iN)*C2
  pmol(iC2N2)  = Kd(iC2N2)*pelm(iN)**2*C2
  pmol(iHNO)   = Kd(iHNO)*pelm(iOx)*pelm(iH)*pelm(iN)
  pmol(iHNO2)  = Kd(iHNO2)*O2*pelm(iH)*pelm(iN)
  pmol(iHNO3)  = Kd(iHNO3)*pelm(iOx)**3*pelm(iH)*pelm(iN)
  pmol(iNH2)   = Kd(iNH2)*pelm(iN)*H2
  pmol(iNO2)   = Kd(iNO2)*pelm(iN)*O2
  pmol(iNO3)   = Kd(iNO3)*pelm(iN)*pelm(iOx)**3
  pmol(iN2O)   = Kd(iN2O)*pelm(iN)**2*pelm(iOx)
  pmol(iN2O4)  = Kd(iN2O4)*pelm(iN)**2*pelm(iOx)**4
  pmol(iMgH)   = Kd(iMgH)*pelm(iH)*pelm(iMg)
  pmol(iMgOH)  = Kd(iMgOH)*pelm(iH)*pelm(iMg)*pelm(iOx)
  pmol(iMgO2H2) = Kd(iMgO2H2)*H2*pelm(iMg)*O2
  pmol(iMgN)   = Kd(iMgN)*pelm(iN)*pelm(iMg)
  pmol(iMgO)   = Kd(iMgO)*pelm(iOx)*pelm(iMg)
  pmol(iSiC)   = Kd(iSiC)*pelm(iC)*pelm(iSi)
  pmol(iSiH2)  = Kd(iSiH2)*H2*pelm(iSi)
  pmol(iSiH3)  = Kd(iSiH3)*pelm(iH)**3*pelm(iSi)
  pmol(iSiN)   = Kd(iSiN)*pelm(iN)*pelm(iSi)
  pmol(iSiO2)  = Kd(iSiO2)*O2*pelm(iSi)
  pmol(iFeO)   = Kd(iFeO)*pelm(iOx)*pelm(iFe)
  pmol(iFeO2H2) = Kd(iFeO2H2)*H2*pelm(iFe)*O2
  pmol(iCOS)   = Kd(iCOS)*pelm(iOx)*pelm(iC)*pelm(iS)
  pmol(iCS)    = Kd(iCS)*pelm(iS)*pelm(iC)
  pmol(iCS2)   = Kd(iCS2)*pelm(iC)*pelm(iS)**2
  pmol(iFeS)   = Kd(iFeS)*pelm(iS)*pelm(iFe)
  pmol(iH2SO4) = Kd(iH2SO4)*H2*pelm(iS)*pelm(iOx)**4
  pmol(iMgS)   = Kd(iMgS)*pelm(iS)*pelm(iMg)
  pmol(iSN)    = Kd(iSN)*pelm(iN)*pelm(iS)
  pmol(iSO)    = Kd(iSO)*pelm(iOx)*pelm(iS)
  pmol(iSO2)   = Kd(iSO2)*O2*pelm(iS)
  pmol(iSO3)   = Kd(iSO3)*pelm(iOx)**3*pelm(iS)
  pmol(iTiO)   = Kd(iTiO)*pelm(iOx)*pelm(iTi)
  pmol(iTiO2)  = Kd(iTiO2)*O2*pelm(iTi)
  pmol(iTiS)   = Kd(iTiS)*pelm(iS)*pelm(iTi)
 end subroutine upd_pmol

end subroutine chemical_equilibrium

!---------------------------------------------------------------------
! Positive root of a*x**2 + b*x + c = 0 (Citardauq form).
! b >= 1 in all call sites, so b + sqrt(delta) >= b > 0 and the old
! "denom <= 0" branches were dead code. max(delta,0) is kept because
! c can be > 0 (e.g. fcarb > 1), in which case the root is floored.
!---------------------------------------------------------------------
pure elemental real function solve_q(a, b, c)
    ! Citardauq form: avoids cancellation. Assumes b > 0 and c < 0,
    ! so the denominator is always >= b > 0.
    real, intent(in) :: a, b, c
    solve_q = max(-2.*c / (b + sqrt(max(b*b - 4.*a*c, 0.))), tiny(solve_q))
 end function solve_q

!---------------------------------------------------------------------
! Positive root of a*x**4 + b*x**3 + c*x**2 + d*x + e = 0
! with a,b,c >= 0, d > 0. For e < 0 the polynomial is convex and
! increasing on x > 0, so Newton converges monotonically from the
! right of the root and x_ub = -e/d is a valid upper bound.
! maxit: maximum number of Newton steps (use a small number such as
! 1-3 inside the outer fixed-point loop if you want, 50 = full solve).
!---------------------------------------------------------------------
 real function newton_method(a, b, c, d, e, zz, maxit)
 real,    intent(in) :: a, b, c, d, e, zz
 integer, intent(in) :: maxit
 real, parameter :: tol_step = 1.e-4   ! error after the step ~ tol_step**2
 real    :: x, x_ub, fx, dfx, dx
 integer :: iter

 if (e >= 0.0) then            ! oxygen fully depleted: no positive root
    newton_method = 0.0
    return
 endif

 x_ub = -e / d
 if (c > 0.0) x_ub = min(x_ub, sqrt(-e / c))
 if (a > 0.0) x_ub = min(x_ub, sqrt(sqrt(-e / a)))

 if (zz > 0.0) then
    x = min(zz, x_ub)
 else
    x = x_ub
 endif

 do iter = 1, maxit
    fx  = (((a*x + b)*x + c)*x + d)*x + e
    dfx = ((4.0*a*x + 3.0*b)*x + 2.0*c)*x + d
    dx  = -fx / dfx
    x   = min(x + dx, x_ub)
    if (abs(dx) <= tol_step*x) exit
 enddo
 newton_method = x
end function newton_method


pure real function calc_Kd(coefs, T)
 real, intent(in) :: coefs(5), T
 real, parameter :: R = 1.987165 !in cal/mol/K
 real :: G, d, T_cap
 ! limit T to avoid numerical issues
 T_cap = max(T,400.)
 G = coefs(1)/T_cap + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T_cap)*T_cap)*T_cap
 d = max(-exp_lim,min(-G/(R*T_cap),exp_lim)) !222
 calc_Kd = exp(d)
end function calc_Kd

pure real function calc_Kd_tsuji(coefs, T)
 use physcon, only:patm
 real, parameter :: log10_patm=log10(patm)
 real, intent(in) :: coefs(6), T
 real :: theta, logKd, T_cap, d
 ! limit T to avoid numerical issues
 T_cap = max(T,400.)
 theta = 5040./T_cap
 d = coefs(1)+(coefs(2)+(coefs(3)+(coefs(4)+coefs(5)*theta)*theta)*theta)*theta+(coefs(6)-1.)*log10_patm
 logKd = max(-log_lim,min(d,log_lim))
 calc_Kd_tsuji = 10.**(-logKd)
end function calc_Kd_tsuji

end module chemistry_condensation
