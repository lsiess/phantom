!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module chemistry_condensation

 implicit none
 public :: set_abundances_condensation,init_muGamma_condensation,calc_muGamma_condensation,&
      network,write_headeropts_dust_condensation,read_headeropts_dust_condensation

 !      real :: wind_CO_ratio = 0.34 !Based on Ferrarotti and Gail 2002

! Indices for elements and molecules:
 integer, parameter, public :: nElements = 11
 real, public :: mass_per_H, eps(nElements)
 integer, parameter, public :: iH=1, iHe=2, iC=3, iOx=4, iN=5, iNe=6, iSi=7, iS=8, iFe=9, iTi=10, iMg=11

 integer, parameter, public :: iH2=1, iOH=2, iH2O=3, iCO=4, iCO2=5, iCH4=6, iC2H=7, iC2H2=8, iN2=9, iNH3=10, iCN=11, &
       iHCN=12, iSi2=13, iSi3=14, iSiO=15, iSi2C=16, iSiH4=17, iS2=18, iHS=19, iH2S=20, iSiS=21, &
       iSiH=22, iTiO=23, iTiO2=24, iC2=25, iO2=26, iCH=27,iCOH=28, iC2O=29, iCH2=30, iH2CO=31, iCH3=32, &
       iC2H4=33, iNH=34, iNO=35, iNCO=36, iHCNO=37, iC2N=38, iC2N2=39, iHNO=40, iHNO2=41, iHNO3=42, &
       iNH2=43, iNO2=44, INO3=45, iN2O=46, iN2O4=47, iMgH=48, iMgOH=49, iMgO2H2=50, iMgN=51, iMgO=52, &
       iSiC=53, iSiH2=54, iSiH3=55, iSiN=56, iSiO2=57, iFeO=58, iFeO2H2=59, iCOS=60, iCS=61, &
       iCS2=62, iFeS=63, iH2SO4=64, iMgS=65, iSN=66, iSO=67, iSO2=68, iSO3=69, iTiS=70

 private
 integer, parameter :: nMolecules = 69

 !real :: mass_per_H, eps(nElements)
 real :: Aw(nElements) = [1.0079, 4.0026, 12.011, 15.9994, 14.0067, 20.17, 28.0855, 32.06, 55.847, 47.867, 24.305]

 !Coefficients can be found in the appendix of the book by Gail & Sedlmayr
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
       2.46640d+01, -6.88730d+00, 8.37100d-02, -1.00580d-02, 4.92910d-04, & !SiH2
       3.63290d+01, -1.05560d+01, 8.09450d-02, -8.62120d-03, 3.98640d-04, & !SiH3
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

       integer, public, parameter :: ncols = nElements + nMolecules -2 !Not including Ne and Ti
       !The vector columns include molecules and atoms

contains


subroutine set_abundances_condensation(wind_CO_ratio)
! all quantities in cgs
 use physcon, only:atomic_mass_unit
 real, intent(in) :: wind_CO_ratio
 eps(iH)  = 1.0
 eps(iHe) = 1.009d-1 !1.04d-1
 eps(iOx) = 7.211d-4 !6.87d-4
 eps(iN)  = 2.106d-4 !2.52d-4
 eps(iMg) = 3.859d-5 !3.85d-5
 eps(iNe) = 1.18d-4  !1.17d-4
 eps(iSi) = 3.561d-5 !3.58d-5
 eps(iS)  = 1.863d-5 !1.85d-5
 eps(iFe) = 3.241d-5 !3.24d-5
 eps(iTi) = 8.621d-8 !8.6d-8
 eps(iC)  = eps(iOx) * wind_CO_ratio
 mass_per_H = atomic_mass_unit*dot_product(Aw,eps)
 !XH  = atomic_mass_unit*eps(iH)/mass_per_H  ! H mass fraction
 !XHe = atomic_mass_unit*eps(iHe)/mass_per_H ! He mass fraction
end subroutine set_abundances_condensation

!-----------------------------------------------------------------------
!+
!  write relevant options to the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine write_headeropts_dust_condensation(wind_CO_ratio,hdr,ierr)
 use dump_utils,        only:dump_h,add_to_rheader
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr
 real, intent(in) :: wind_CO_ratio

! initial gas composition for dust formation
 call set_abundances_condensation(wind_CO_ratio)
 call add_to_rheader(eps,'epsilon',hdr,ierr) ! array
 call add_to_rheader(Aw,'Amean',hdr,ierr)    ! array
 call add_to_rheader(mass_per_H,'mass_per_H',hdr,ierr) ! real

end subroutine write_headeropts_dust_condensation

!-----------------------------------------------------------------------
!+
!  read relevant options from the header of the dump file
!+
!-----------------------------------------------------------------------
subroutine read_headeropts_dust_condensation(wind_CO_ratio,hdr,ierr)
 use dump_utils, only:dump_h,extract
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 real :: dum(nElements)
 real, intent(in) :: wind_CO_ratio

 ierr = 0
 call extract('mass_per_H',mass_per_H,hdr,ierr) ! real
 ! it is likely that your dump was generated with an old version of phantom
 ! and the chemical properties not stored. restore and save the default values
 if (mass_per_H < tiny(0.)) then
    print *,'reset dust chemical network properties'
    call set_abundances_condensation(wind_CO_ratio)
    call extract('epsilon',dum(1:nElements),hdr,ierr) ! array
    call extract('Amean',dum(1:nElements),hdr,ierr) ! array
 else
    call extract('epsilon',eps(1:nElements),hdr,ierr) ! array
    call extract('Amean',Aw(1:nElements),hdr,ierr) ! array
 endif
end subroutine read_headeropts_dust_condensation

!-----------------------------------------------------------------------
!+
!  driver to solve the chemical network
!+
!-----------------------------------------------------------------------
subroutine network(T,rho_cgs,mu,gamma,wind_CO_ratio,pH_tot,fol,fpy,fqu,fir,fsc,fcarb,abundance,pressure_cgs)
 use physcon, only:patm,kboltz
 real, intent(in)     :: wind_CO_ratio, fol, fpy, fqu, fir, fsc, fcarb !rho_cgs
 real, intent(out)    :: mu,gamma,pH_tot
 real, intent(out),optional :: abundance(ncols)
 real, intent(inout)  :: T, rho_cgs
 real, intent(in), optional :: pressure_cgs
!       real :: epsC, nH_tot
!       real :: v1
 integer :: j
 real, dimension(nMolecules)    :: pmol  !, nmol=0.
 real                           :: pelm(nElements)

       pelm = 0.
       pmol = 0.

      !  call set_abundances_condensation(wind_CO_ratio)
      !  print *,'SHOULD  BE DONE ONLY ONCE '

      if (present(pressure_cgs)) then
         call init_muGamma_condensation(rho_cgs,T,mu,gamma,pH_tot,pelm(iH),pmol(iH2),pressure_cgs)
      else
         call init_muGamma_condensation(rho_cgs,T,mu,gamma,pH_tot,pelm(iH),pmol(iH2))
      endif

      pelm(iHe) = eps(iHe)*pH_tot
      call chemical_equilibrium(rho_cgs, T, pmol, pelm, mu, gamma, &
           wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb) !, nMolecules)

      if (present(abundance)) then
      !convert the pressure of the molecules in number density
         do j=1,ncols-9
            abundance(j) = max(1.d-50,pmol(j)*patm/(kboltz*T))  !*patm
         enddo

         !convert the pressure of the elements in number density
         abundance(70) = max(1.d-50,pelm(iH)*patm/(kboltz*T)) !*patm
         abundance(71) = max(1.d-50,pelm(iHe)*patm/(kboltz*T)) !*patm
         abundance(72) = max(1.d-50,pelm(iC)*patm/(kboltz*T)) !*patm
         abundance(73) = max(1.d-50,pelm(iOx)*patm/(kboltz*T)) !*patm
         abundance(74) = max(1.d-50,pelm(iN)*patm/(kboltz*T)) !*patm
         abundance(75) = max(1.d-50,pelm(iSi)*patm/(kboltz*T)) !*patm
         abundance(76) = max(1.d-50,pelm(iS)*patm/(kboltz*T)) !*patm
         abundance(77) = max(1.d-50,pelm(iFe)*patm/(kboltz*T)) !*patm
         abundance(78) = max(1.d-50,pelm(iMg)*patm/(kboltz*T)) !*patm

         !abundance(iH2O) = abundance(iH2O) - (3.*fol + 2.*fpy + fqu) * eps(iSi)*patm *pH_tot/(kboltz*T)
         !abundance(iSiO) = abundance(iSiO) - (fol + fpy + fqu) * eps(iSi)*patm * pH_tot/(kboltz*T)
         !abundance(78) = abundance(78) - (2.*fol + fpy) * eps(iSi)*patm * pH_tot/(kboltz*T)
         !if (abundance(iH2O) < 0.0) abundance(iH2O) = 0.0
         !if (abundance(iSiO) < 0.0) abundance(iSiO) = 0.0
         !if (abundance(78) < 0.0)   abundance(78) = 0.0

         !!Update the abundance of Si and Oxygen.

       endif

end subroutine network

!---------------------------------------------------------------
!
!  Compute partial pressures assuming chemical equilibrum
!
!---------------------------------------------------------------
subroutine chemical_equilibrium(rho_cgs,T,pmol,pelm,mu,gamma,&
     wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb)
! all quantities are in cgs
 real, intent(inout) :: T, mu, gamma, pmol(nMolecules), pelm(nElements), rho_cgs
 real, intent(in)    :: wind_CO_ratio, pH_tot, fol, fpy, fqu, fir, fsc, fcarb
 real, dimension(nMolecules)    :: pmol_old
 real    :: Kd(nMolecules+1), err(nElements) !, pH_tot !, a, b, c, d !LUIS err
 real    :: pelm_old(nElements)
 integer :: i, nit, rndnmbr

 pelm_old = 0.
 pmol_old = 0.

 if (T > 1.d4) then
    pelm = eps*pH_tot
    pmol = 0.
    return
 endif

! Dissociation constants
 do i=1,nMolecules
    Kd(i) = calc_Kd(coefs(:,i), T)
 enddo
 Kd(iTiS) = calc_Kd_TiS(T)

 err     = 1.
 nit     = 0
 rndnmbr = 0

 do while (maxval(err) > 1.e-6)

 pelm_old(:) = pelm(:)

 if (wind_CO_ratio < 1.0) then
 if (rndnmbr == 0) then
    pelm(iOx) = newton_method(4.*(pelm(iN)**2*Kd(iN2O4)+pelm(iH)**2*pelm(iS)*Kd(iH2SO4)), &
                         3.*(pelm(iH)*pelm(iN)*Kd(iHNO3)+pelm(iN)*Kd(iNO3)+pelm(iS)*Kd(iSO3)), &
                         2.*(Kd(iO2)+pelm(iC)*Kd(iCO2)+pelm(iH)*pelm(iN)*Kd(iHNO2) &
                             +pelm(iN)*Kd(iNO2)+pelm(iH)**2*pelm(iMg)*Kd(iMgO2H2)+pelm(iSi)*Kd(iSiO2) &
                             +pelm(iFe)*pelm(iH)**2*Kd(iFeO2H2)+pelm(iS)*Kd(iSO2)), &
                         1.+pelm(iH)*Kd(iOH)+pelm(iH)**2*Kd(iH2O)+pelm(iC)*pelm(iH)*Kd(iCOH) &
                         +pelm(iC)**2*Kd(iC2O)+pelm(iC)*pelm(iH)**2*Kd(iH2CO)+pelm(iN)*Kd(iNO) &
                         +pelm(iC)*pelm(iN)*Kd(iNCO)+pelm(iC)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iH)*pelm(iN)*Kd(iHNO) &
                         +pelm(iN)**2*Kd(iN2O)+pelm(iH)*pelm(iMg)*Kd(iMgOH)+pelm(iMg)*Kd(iMgO)+pelm(iFe)*Kd(iFeO) &
                         +pelm(iC)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iSO), &
                         (eps(iC)+eps(iSi)*(1.+4.*fol+3.*fpy+2.*fqu)-eps(iOx))*pH_tot, & !Quitar fol a la molecula, no el elemento
                          eps(iOx)*pH_tot) !pmol(iCO)  = eps(iC)*pH_tot, pmol(iSiO) = eps(iSi)*pH_tot
                 rndnmbr = rndnmbr + 1
 else !rndnmbr
    pelm(iOx) = newton_method(4.*(pelm(iN)**2*Kd(iN2O4)+pelm(iH)**2*pelm(iS)*Kd(iH2SO4)), &
                         3.*(pelm(iH)*pelm(iN)*Kd(iHNO3)+pelm(iN)*Kd(iNO3)+pelm(iS)*Kd(iSO3)), &
                         2.*(Kd(iO2)+pelm(iC)*Kd(iCO2)+pelm(iH)*pelm(iN)*Kd(iHNO2) &
                             +pelm(iN)*Kd(iNO2)+pelm(iH)**2*pelm(iMg)*Kd(iMgO2H2)+pelm(iSi)*Kd(iSiO2) &
                             +pelm(iFe)*pelm(iH)**2*Kd(iFeO2H2)+pelm(iS)*Kd(iSO2)), &
                         1.+pelm(iC)*Kd(iCO)+pelm(iSi)*Kd(iSiO)+pelm(iH)*Kd(iOH)+pelm(iH)**2*Kd(iH2O)+pelm(iC)*pelm(iH)*Kd(iCOH) &
                         +pelm(iC)**2*Kd(iC2O)+pelm(iC)*pelm(iH)**2*Kd(iH2CO)+pelm(iN)*Kd(iNO) &
                         +pelm(iC)*pelm(iN)*Kd(iNCO)+pelm(iC)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iH)*pelm(iN)*Kd(iHNO) &
                         +pelm(iN)**2*Kd(iN2O)+pelm(iH)*pelm(iMg)*Kd(iMgOH)+pelm(iMg)*Kd(iMgO)+pelm(iFe)*Kd(iFeO) &
                         +pelm(iC)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iSO), &
                         (eps(iSi)*(4.*fol+3.*fpy+2.*fqu)-eps(iOx))*pH_tot, &
                         eps(iOx)*pH_tot)
 endif !rndnmbr
 ! pelm(iC) = newton_method(0., &
 !            0., &
 !            2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
 !            +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
 !            1.+pelm(iH)*Kd(iCH)+pelm(iOx)*Kd(iCO)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
 !            +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
 !            +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
 !            +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
 !            +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
 !            (fcarb-1.0)*eps(iC)*pH_tot,eps(iC)*pH_tot)
 pelm(iC) = solve_q(2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
            +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
            1.+pelm(iH)*Kd(iCH)+pelm(iOx)*Kd(iCO)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
            +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
            +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
            +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
            +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
            (fcarb-1.0)*eps(iC)*pH_tot)

else !%% wind_CO is larger than 1.0
  if (rndnmbr == 0) then
   ! pelm(iC) = newton_method(0., &
   !          0., &
   !          2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
   !          +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
   !          1.+pelm(iH)*Kd(iCH)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
   !          +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
   !          +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
   !          +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
   !          +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
   !          (eps(iOx)+(fcarb-1.0)*eps(iC))*pH_tot, & !pmol(iCO)  = eps(iOx)*pH_tot
   !          eps(iC)*pH_tot)
   pelm(iC) = solve_q(2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
            +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
            1.+pelm(iH)*Kd(iCH)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
            +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
            +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
            +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
            +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
            (eps(iOx)+(fcarb-1.0)*eps(iC))*pH_tot)

            rndnmbr = rndnmbr + 1

  else

   ! pelm(iC) = newton_method(0., &
   !          0., &
   !          2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
   !          +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
   !          1.+pelm(iH)*Kd(iCH)+pelm(iOx)*Kd(iCO)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
   !          +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
   !          +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
   !          +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
   !          +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
   !          (fcarb-1.0)*eps(iC)*pH_tot,eps(iC)*pH_tot)
   pelm(iC) = solve_q(2.*(Kd(iC2)+pelm(iH)*Kd(iC2H)+pelm(iOx)*Kd(iC2O)+pelm(iH)**2*Kd(iC2H2)+pelm(iH)**4*Kd(iC2H4) &
            +pelm(iN)*Kd(iC2N)+pelm(iN)**2*Kd(iC2N2)), &
            1.+pelm(iH)*Kd(iCH)+pelm(iOx)*Kd(iCO)+pelm(iOx)*pelm(iH)*Kd(iCOH)+pelm(iOx)**2*Kd(iCO2) &
            +pelm(iH)**2*Kd(iCH2)+pelm(iOx)*pelm(iH)**2*Kd(iH2CO)+pelm(iH)**3*Kd(iCH3)+pelm(iH)**4*Kd(iCH4) &
            +pelm(iN)*Kd(iCN)+pelm(iN)*pelm(iH)*Kd(iHCN)+pelm(iOx)*pelm(iN)*Kd(iNCO) &
            +pelm(iOx)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iSi)*Kd(iSiC)+pelm(iSi)**2*Kd(iSi2C) &
            +pelm(iOx)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iCS)+pelm(iS)**2*Kd(iCS2), &
            (fcarb-1.0)*eps(iC)*pH_tot)
  endif  !rndnmbr

  pelm(iOx) = newton_method(4.*(pelm(iN)**2*Kd(iN2O4)+pelm(iH)**2*pelm(iS)*Kd(iH2SO4)), &
                         3.*(pelm(iH)*pelm(iN)*Kd(iHNO3)+pelm(iN)*Kd(iNO3)+pelm(iS)*Kd(iSO3)), &
                         2.*(Kd(iO2)+pelm(iC)*Kd(iCO2)+pelm(iH)*pelm(iN)*Kd(iHNO2) &
                             +pelm(iN)*Kd(iNO2)+pelm(iH)**2*pelm(iMg)*Kd(iMgO2H2)+pelm(iSi)*Kd(iSiO2) &
                             +pelm(iFe)*pelm(iH)**2*Kd(iFeO2H2)+pelm(iS)*Kd(iSO2)), &
                         1.+pelm(iC)*Kd(iCO)+pelm(iSi)*Kd(iSiO)+pelm(iH)*Kd(iOH)+pelm(iH)**2*Kd(iH2O)+pelm(iC)*pelm(iH)*Kd(iCOH) &
                         +pelm(iC)**2*Kd(iC2O)+pelm(iC)*pelm(iH)**2*Kd(iH2CO)+pelm(iN)*Kd(iNO) &
                         +pelm(iC)*pelm(iN)*Kd(iNCO)+pelm(iC)*pelm(iH)*pelm(iN)*Kd(iHCNO)+pelm(iH)*pelm(iN)*Kd(iHNO) &
                         +pelm(iN)**2*Kd(iN2O)+pelm(iH)*pelm(iMg)*Kd(iMgOH)+pelm(iMg)*Kd(iMgO)+pelm(iFe)*Kd(iFeO) &
                         +pelm(iC)*pelm(iS)*Kd(iCOS)+pelm(iS)*Kd(iSO), &
                         (eps(iSi)*(4.*fol+3.*fpy+2.*fqu)-eps(iOx))*pH_tot, &
                         eps(iOx)*pH_tot)
endif !wind_CO


 ! pelm(iN) = newton_method(0., 0., &
 !            2.*(Kd(iN2)+pelm(iC)**2*Kd(iC2N2) &
 !            +pelm(iOx)*Kd(iN2O) &
 !            +pelm(iOx)**4*Kd(iN2O4)), &
 !            1.+pelm(iH)*Kd(iNH)+pelm(iOx)*Kd(iNO)+pelm(iC)*Kd(iCN)+pelm(iC)*pelm(iH)*Kd(iHCN) &
 !            +pelm(iOx)*pelm(iC)*Kd(iNCO)+pelm(iC)*pelm(iOx)*pelm(iH)*Kd(iHCNO)+pelm(iC)**2*Kd(iC2N) &
 !            +pelm(iOx)*pelm(iH)*Kd(iHNO)+pelm(iOx)**2*pelm(iH)*Kd(iHNO2) &
 !            +pelm(iOx)**3*pelm(iH)*Kd(iHNO3)+pelm(iH)**2*Kd(iNH2)+pelm(iH)**3*Kd(iNH3) &
 !            +pelm(iOx)**2*Kd(iNO2)+pelm(iOx)**3*Kd(iNO3)+pelm(iMg)*Kd(iMgN) &
 !            +pelm(iSi)*Kd(iSiN) &
 !            +pelm(iS)*Kd(iSN), &
 !            -eps(iN)*pH_tot,eps(iN)*pH_tot)
 pelm(iN) = solve_q(2.*(Kd(iN2)+pelm(iC)**2*Kd(iC2N2) &
            +pelm(iOx)*Kd(iN2O) &
            +pelm(iOx)**4*Kd(iN2O4)), &
            1.+pelm(iH)*Kd(iNH)+pelm(iOx)*Kd(iNO)+pelm(iC)*Kd(iCN)+pelm(iC)*pelm(iH)*Kd(iHCN) &
            +pelm(iOx)*pelm(iC)*Kd(iNCO)+pelm(iC)*pelm(iOx)*pelm(iH)*Kd(iHCNO)+pelm(iC)**2*Kd(iC2N) &
            +pelm(iOx)*pelm(iH)*Kd(iHNO)+pelm(iOx)**2*pelm(iH)*Kd(iHNO2) &
            +pelm(iOx)**3*pelm(iH)*Kd(iHNO3)+pelm(iH)**2*Kd(iNH2)+pelm(iH)**3*Kd(iNH3) &
            +pelm(iOx)**2*Kd(iNO2)+pelm(iOx)**3*Kd(iNO3)+pelm(iMg)*Kd(iMgN) &
            +pelm(iSi)*Kd(iSiN) &
            +pelm(iS)*Kd(iSN), &
            -eps(iN)*pH_tot)

 pelm(iMg) = (eps(iMg)-(2.*fol+fpy)*eps(iSi))*pH_tot/(1.+pelm(iH)*Kd(iMgH)+pelm(iOx)*pelm(iH)*Kd(iMgOH) &
             +pelm(iOx)**2*pelm(iH)**2*Kd(iMgO2H2)+pelm(iN)*Kd(iMgN) &
             +pelm(iOx)*Kd(iMgO)+pelm(iS)*Kd(iMgS))


 ! pelm(iSi) = newton_method(0., 0., &
 !             2.*(Kd(iSi2)+pelm(iC)*Kd(iSi2C)), &
 !             1.+pelm(iC)*Kd(iSiC)+pelm(iH)*Kd(iSiH)+pelm(iH)**2*Kd(iSiH2) &
 !             +pelm(iH)**3*Kd(iSiH3)+pelm(iH)**4*Kd(iSiH4)+pelm(iN)*Kd(iSiN)+pelm(iOx)*Kd(iSiO) &
 !             +pelm(iOx)**2*Kd(iSiO2) &
 !             +pelm(iS)*Kd(iSiS), &
 !             (fol+fpy+fqu+fsc-1.0)*eps(iSi)*pH_tot, &
 !             eps(iSi)*pH_tot)
 pelm(iSi) = solve_q(2.*(Kd(iSi2)+pelm(iC)*Kd(iSi2C)), &
             1.+pelm(iC)*Kd(iSiC)+pelm(iH)*Kd(iSiH)+pelm(iH)**2*Kd(iSiH2) &
             +pelm(iH)**3*Kd(iSiH3)+pelm(iH)**4*Kd(iSiH4)+pelm(iN)*Kd(iSiN)+pelm(iOx)*Kd(iSiO) &
             +pelm(iOx)**2*Kd(iSiO2) &
             +pelm(iS)*Kd(iSiS), &
             (fol+fpy+fqu+fsc-1.0)*eps(iSi)*pH_tot)


 pelm(iFe) = (1.0 - fir)*eps(iFe)*pH_tot/(1.+pelm(iOx)*Kd(iFeO)+pelm(iOx)**2*pelm(iH)**2*Kd(iFeO2H2) &
             +pelm(iS)*Kd(iFeS))

 ! pelm(iS) = newton_method(0., 0., 2.*(Kd(iS2)+pelm(iC)*Kd(iCS2)), &
 !            1.+pelm(iOx)*pelm(iC)*Kd(iCOS)+pelm(iC)*Kd(iCS)+pelm(iFe)*Kd(iFeS)+pelm(iH)*Kd(iHS) &
 !            +pelm(iH)**2*Kd(iH2S)+pelm(iOx)**4*pelm(iH)**2*Kd(iH2SO4)+pelm(iMg)*Kd(iMgS)+pelm(iN)*Kd(iSN) &
 !            +pelm(iOx)*Kd(iSO)+pelm(iOx)**2*Kd(iSO2)+pelm(iOx)**3*Kd(iSO3)+pelm(iSi)*Kd(iSiS), &
 !            -eps(iS)*pH_tot,eps(iS)*pH_tot)
 pelm(iS) = solve_q(2.*(Kd(iS2)+pelm(iC)*Kd(iCS2)), &
            1.+pelm(iOx)*pelm(iC)*Kd(iCOS)+pelm(iC)*Kd(iCS)+pelm(iFe)*Kd(iFeS)+pelm(iH)*Kd(iHS) &
            +pelm(iH)**2*Kd(iH2S)+pelm(iOx)**4*pelm(iH)**2*Kd(iH2SO4)+pelm(iMg)*Kd(iMgS)+pelm(iN)*Kd(iSN) &
            +pelm(iOx)*Kd(iSO)+pelm(iOx)**2*Kd(iSO2)+pelm(iOx)**3*Kd(iSO3)+pelm(iSi)*Kd(iSiS), &
            -eps(iS)*pH_tot)

  err = abs((pelm-pelm_old)/(pelm_old+1.0e-60))

  pmol(iH2)  = Kd(iH2)*pelm(iH)**2
  pmol(iOH)  = Kd(iOH)*pelm(iOx)*pelm(iH)
  pmol(iH2O) = Kd(iH2O)*pelm(iOx)*pelm(iH)**2
  pmol(iCO)  = Kd(iCO)*pelm(iOx)*pelm(iC)
  pmol(iCO2)  = Kd(iCO2)*pelm(iC)*pelm(iOx)**2
  pmol(iCH4)  = Kd(iCH4)*pelm(iC)*pelm(iH)**4
  pmol(iC2H)  = Kd(iC2H)*pelm(iH)*pelm(iC)**2
  pmol(iC2H2)  = Kd(iC2H2)*pelm(iC)**2*pelm(iH)**2
  pmol(iN2)  = Kd(iN2)*pelm(iN)**2
  pmol(iNH3)  = Kd(iNH3)*pelm(iN)*pelm(iH)**3
  pmol(iCN)  = Kd(iCN)*pelm(iC)*pelm(iN)
  pmol(iHCN)  = Kd(iHCN)*pelm(iH)*pelm(iC)*pelm(iN)
  pmol(iSi2)  = Kd(iSi2)*pelm(iSi)**2
  pmol(iSi3)  = Kd(iSi3)*pelm(iSi)**3
  pmol(iSiO)  = Kd(iSiO)*pelm(iSi)*pelm(iOx)
  pmol(iSi2C)  = Kd(iSi2C)*pelm(iSi)**2*pelm(iC)
  pmol(iSiH4)  = Kd(iSiH4)*pelm(iH)**4*pelm(iSi)
  pmol(iS2)  = Kd(iS2)*pelm(iS)**2
  pmol(iHS)  = Kd(iHS)*pelm(iH)*pelm(iS)
  pmol(iH2S)  = Kd(iH2S)*pelm(iH)**2*pelm(iS)
  pmol(iSiS)  = Kd(iSiS)*pelm(iSi)*pelm(iS)
  pmol(iSiH)  = Kd(iSiH)*pelm(iSi)*pelm(iH)
  pmol(iC2)  = Kd(iC2)*pelm(iC)**2
  pmol(iO2)  = Kd(iO2)*pelm(iOx)**2
  pmol(iCH)  = Kd(iCH)*pelm(iH)*pelm(iC)
  pmol(iCOH)  = Kd(iCOH)*pelm(iOx)*pelm(iC)*pelm(iH)
  pmol(iC2O)  = Kd(iC2O)*pelm(iOx)*pelm(iC)**2
  pmol(iCH2)  = Kd(iCH2)*pelm(iH)**2*pelm(iC)
  pmol(iH2CO)  = Kd(iH2CO)*pelm(iOx)*pelm(iC)*pelm(iH)**2
  pmol(iCH3)  = Kd(iCH3)*pelm(iH)**3*pelm(iC)
  pmol(iC2H4)  = Kd(iC2H4)*pelm(iH)**4*pelm(iC)**2
  pmol(iNH)  = Kd(iNH)*pelm(iN)*pelm(iH)
  pmol(iNO)  = Kd(iNO)*pelm(iOx)*pelm(iN)
  pmol(iNCO)  = Kd(iNCO)*pelm(iOx)*pelm(iC)*pelm(iN)
  pmol(iHCNO)  = Kd(iHCNO)*pelm(iOx)*pelm(iC)*pelm(iH)*pelm(iN)
  pmol(iC2N)  = Kd(iC2N)*pelm(iN)*pelm(iC)**2
  pmol(iC2N2)  = Kd(iC2N2)*pelm(iN)**2*pelm(iC)**2
  pmol(iHNO)  = Kd(iHNO)*pelm(iOx)*pelm(iH)*pelm(iN)
  pmol(iHNO2)  = Kd(iHNO2)*pelm(iOx)**2*pelm(iH)*pelm(iN)
  pmol(iHNO3)  = Kd(iHNO3)*pelm(iOx)**3*pelm(iH)*pelm(iN)
  pmol(iNH2)  = Kd(iNH2)*pelm(iN)*pelm(iH)**2
  pmol(iNO2)  = Kd(iNO2)*pelm(iN)*pelm(iOx)**2
  pmol(iNO3)  = Kd(iNO3)*pelm(iN)*pelm(iOx)**3
  pmol(iN2O)  = Kd(iN2O)*pelm(iN)**2*pelm(iOx)
  pmol(iN2O4)  = Kd(iN2O4)*pelm(iN)**2*pelm(iOx)**4
  pmol(iMgH)  = Kd(iMgH)*pelm(iH)*pelm(iMg)
  pmol(iMgOH)  = Kd(iMgOH)*pelm(iH)*pelm(iMg)*pelm(iOx)
  pmol(iMgO2H2)  = Kd(iMgO2H2)*pelm(iH)**2*pelm(iMg)*pelm(iOx)**2
  pmol(iMgN)  = Kd(iMgN)*pelm(iN)*pelm(iMg)
  pmol(iMgO)  = Kd(iMgO)*pelm(iOx)*pelm(iMg)
  pmol(iSiC)  = Kd(iSiC)*pelm(iC)*pelm(iSi)
  pmol(iSiH2)  = Kd(iSiH2)*pelm(iH)**2*pelm(iSi)
  pmol(iSiH3)  = Kd(iSiH3)*pelm(iH)**3*pelm(iSi)
  pmol(iSiN)  = Kd(iSiN)*pelm(iN)*pelm(iSi)
  pmol(iSiO2)  = Kd(iSiO2)*pelm(iOx)**2*pelm(iSi)
  pmol(iFeO)  = Kd(iFeO)*pelm(iOx)*pelm(iFe)
  pmol(iFeO2H2)  = Kd(iFeO2H2)*pelm(iH)**2*pelm(iFe)*pelm(iOx)**2
  pmol(iCOS)  = Kd(iCOS)*pelm(iOx)*pelm(iC)*pelm(iS)
  pmol(iCS)  = Kd(iCS)*pelm(iS)*pelm(iC)
  pmol(iCS2)  = Kd(iCS2)*pelm(iC)*pelm(iS)**2
  pmol(iFeS)  = Kd(iFeS)*pelm(iS)*pelm(iFe)
  pmol(iH2SO4)  = Kd(iH2SO4)*pelm(iH)**2*pelm(iS)*pelm(iOx)**4
  pmol(iMgS)  = Kd(iMgS)*pelm(iS)*pelm(iMg)
  pmol(iSN)  = Kd(iSN)*pelm(iN)*pelm(iS)
  pmol(iSO)  = Kd(iSO)*pelm(iOx)*pelm(iS)
  pmol(iSO2)  = Kd(iSO2)*pelm(iOx)**2*pelm(iS)
  pmol(iSO3)  = Kd(iSO3)*pelm(iOx)**3*pelm(iS)

  nit = nit + 1
  if (nit == 200) exit

 enddo
end subroutine chemical_equilibrium


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

pure real function calc_Kd(coefs, T)
! all quantities are in cgs
 real, intent(in) :: coefs(5), T
 real, parameter :: R = 1.987165
 real :: G, d
 G = coefs(1)/T + coefs(2) + (coefs(3)+(coefs(4)+coefs(5)*T)*T)*T
 d = min(-G/(R*T),322.) !222
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

!--------------------------------------------
!
!  Initialise mean molecular weight and gamma
!
!--------------------------------------------
subroutine init_muGamma_condensation(rho_cgs,T,mu,gamma, pH_tot, ppH, ppH2,pressure_cgs)
   ! all quantities are in cgs
    real, intent(in)              :: rho_cgs
    real, intent(inout)           :: T, mu, gamma
    real, intent(out)             :: pH_tot
    real, intent(in),  optional   :: pressure_cgs
    real, intent(out), optional   :: ppH, ppH2
    real :: pH, pH2

    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    gamma  = 5./3.

    !if (present(pressure_cgs)) then
    !  call calc_muGamma_condensation(rho_cgs,T, mu, gamma, pH, pH_tot, pH2, pressure_cgs)
    !else
      call calc_muGamma_condensation(rho_cgs,T, mu, gamma, pH, pH_tot, pH2)
    !endif

    if (present(ppH))  ppH = pH
    if (present(ppH2)) ppH2 = pH2

 end subroutine init_muGamma_condensation

!----------------------------------------
!
!  Calculate mean molecular weight, gamma
!
!----------------------------------------
subroutine calc_muGamma_condensation(rho_cgs,T, mu, gamma, pH, pH_tot, pH2, pressure_cgs)
! all quantities are in cgs
 !use io, only:fatal

 use physcon, only:patm,kboltz
 !real, intent(in)    :: mass_per_H,eps(nElements) !rho_cgs
 real, intent(in) :: rho_cgs
 real, intent(inout) :: mu, gamma, T
 real, intent(out)   :: pH, pH_tot, pH2
 real, intent(in), optional   :: pressure_cgs
 real :: KH2
 real :: T_old, mu_old, gamma_old, tol
 logical :: converged
 integer :: i,isolve
 integer, parameter :: itermax = 100
 character(len=30), parameter :: label = 'calc_muGamma'

 if (present(pressure_cgs)) print *,'WARNING you cannot redefine the density'
 if (T > 1.d4) then
    mu     = (1.+4.*eps(iHe))/(1.+eps(iHe))
    gamma  = 5./3.
    if (present(pressure_cgs)) then
      pH_tot = pressure_cgs/(1.+eps(iHe))/patm
      !rho_cgs = pH_tot*patm*mass_per_H/(T*kboltz)
      !pH = pH_tot
    else
      pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
      !pH     = pH_tot
    endif
    pH = pH_tot
    pH2 = 0.
 elseif (T > 450.) then
! iterate to get consistently pH, T, mu and gamma
    tol       = 1.d-3
    converged = .false.
    isolve    = 0
    !pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H) ! to avoid compiler warning
    !pH        = pH_tot ! arbitrary value, overwritten below, to avoid compiler warning
    i = 0
    do while (.not. converged .and. i < itermax)
       i = i+1
       KH2       = calc_Kd(coefs(:,iH2), T)
       if (present(pressure_cgs)) then
         pH = solve_q(1.,(1.+eps(iHe))/((1.+2.*eps(iHe))*KH2), &
                    -pressure_cgs/patm/((1.+2.*eps(iHe))*KH2))
         pH2 = KH2*pH**2
         pH_tot = pH+2.*pH2
         !rho_cgs = pH_tot*patm*mass_per_H/(T*kboltz)
       else
         pH_tot    = rho_cgs*T*kboltz/(patm*mass_per_H)
         pH        = solve_q(2.*KH2, 1., -pH_tot)
         pH2       = KH2*pH**2
       endif
       mu        = (1.+4.*eps(iHe))*pH_tot/(pH+pH2+eps(iHe)*pH_tot)
       gamma     = (5.*pH+5.*eps(iHe)*pH_tot+7.*pH2)/(3.*pH+3.*eps(iHe)*pH_tot+5.*pH2)
       if (i==1) then
         mu_old    = mu
         gamma_old = gamma
       else
         T_old     = T
         T         = T_old*mu*(gamma-1.)/(mu_old*(gamma_old-1.))
         !T        = T_old    !uncomment this line to cancel iterations
         converged = (abs(T-T_old)/T_old) < tol
         !print *,i,T_old,T,gamma_old,gamma,mu_old,mu,abs(T-T_old)/T_old
         if (i>=itermax .and. .not.converged) then
            if (isolve==0) then
               isolve = isolve+1
               i      = 0
               tol    = 1.d-2
               print *,'[dust_formation] cannot converge on T(mu,gamma). Trying with lower tolerance'
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
   if (present(pressure_cgs)) then
      !print *, "pressure_cgs = ", pressure_cgs
      pH_tot = 2.*pressure_cgs/(1.+2.*eps(iHe))/patm
      !rho_cgs = pH_tot*patm*mass_per_H/(T*kboltz) !*patm
      !print *, "pgas = ", pH_tot*(0.5+eps(iHe))*patm
   else
      pH_tot = rho_cgs*T*kboltz/(patm*mass_per_H)
   endif
    pH2    = pH_tot/2.
    pH     = 0.
    mu     = (1.+4.*eps(iHe))/(0.5+eps(iHe))
    gamma  = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
 endif
end subroutine calc_muGamma_condensation


!Fortran subroutine that implements Newton's method for
!root-finding with the given polynomial function
!f(x) = a*x^4 + b*x^3 + c*x^2 + d*x + e

real function newton_method(aa, bb, cc, dd, ee, zz)
    real, intent(in) :: aa, bb, cc, dd, ee, zz
    real :: a,b,c,d,e,z
    ! Parameters
    integer, parameter :: max_iter = 350 !350
    real, parameter :: tolerance = 1.d-50 !1.d-50
    real, parameter :: tolerance_rel = 1.d-12 !1.d-50

    ! Local variables
    real :: fx, dfx, x
    integer :: iter

    a = aa
    b = bb
    c = cc
    d = dd
    e = ee
    z = zz
    x = z

    ! Newton's method loop
    do iter = 1, max_iter
        ! Compute function value and its derivative
        fx = a*x**4 + b*x**3 + c*x**2 + d*x + e
        dfx = 4*a*x**3 + 3*b*x**2 + 2*c*x + d

        ! Check for convergence
        if (abs(fx) < tolerance .or. abs(fx) < abs(x) * tolerance_rel) exit

        ! Update x using Newton's method
        x = x - fx / dfx

    end do
    newton_method = x
    !print *,iter,max_iter,x,fx,dfx
    !stop

end function newton_method

end module chemistry_condensation
