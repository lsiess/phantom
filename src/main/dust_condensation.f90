!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dust_condensation
!
! Dust formation routine : condensation method
!
! :References: Gail & Sedlmayr textbook Physics and chemistry of Circumstellar dust shells
!              + p[apaers by Ferrarotti & Gail (2000+)
!

use physcon, only:atomic_mass_unit

implicit none

public :: evolve_condensation,dust_growth_condensation

private
 !Coefficients for the calculation of the free entalpy (deltaG), from table A4.1, Gail & Sedlmayr, p651
 real, parameter :: coefficients(5,6) = reshape([&
 3.14987E+05, -3.92876E+06, 1.03657E+03, -1.31876E-02, 0.0d+00, &      !Olivine
 0.00000E+00, -1.86046E+06, 4.54398E+02, -2.75999E-03, 0.0d+00, &      !Quartz
 3.66094E+04, -2.89963E+06, 7.44735E+02, -5.92064E-03, 0.0d+00, &      !Pyroxene
 0.00000E+00, -4.08793E+05, 1.48885E+02, -4.70379E-03, 0.0d+00, &      !Iron
 2.79991E+05, -1.24239E+06, 3.20111E+02, -3.42448E-03, 5.44933E-07, &  !Silicon Carbide
 8.71566E+05, -7.21210E+05, 1.62145E+02, -1.23392E-03, 1.77238E-07],&  !Amorphous Carbon
 shape(coefficients))

!sticking coefficients : Gail & Sedlmayr Book table 12.1, p337
 real, parameter ::   &
     alpha_ol = 0.1,  &
     alpha_qu = 0.05, &
     alpha_py = 0.2,  &
     alpha_ir = 1.0,  &!0.9
     alpha_sc = 0.8,  &
     alpha_carb = 0.3
 !atomic masses, in cgs, from Gail & Sedlmayr Book table 12.1, p337
 real, parameter ::     &
     A_ol = 140.694, &
     A_qu = 60.085,  &
     A_py = 100.389, &
     A_ir = 55.845,  &
     A_sc = 40.10,   &
     A_carb = 12.01
 !density of dust species, in cgs, from Gail & Sedlmayr Book table 12.1, p337
 real, parameter ::     &
     rho_ol = 3.21,  &
     rho_qu = 2.65,  &
     rho_py = 3.19,  &
     rho_ir = 7.87,  &
     rho_sc = 3.21,  &
     rho_carb = 2.20
 !volume of monomers, in cgs, from Gail & Sedlmayr Book table 12.1, p337
 real, parameter :: &
     Vo_ol = A_ol * atomic_mass_unit / rho_ol, &
     Vo_qu = A_qu * atomic_mass_unit / rho_qu, &
     Vo_py = A_py * atomic_mass_unit / rho_py, &
     Vo_ir = A_ir * atomic_mass_unit / rho_ir, &
     Vo_sc = A_sc * atomic_mass_unit / rho_sc, &
     Vo_carb = A_carb * atomic_mass_unit / rho_carb

 real, parameter :: &
     m_SiO  = 44.09*atomic_mass_unit,   &!Si = 28.09  amu, Ox = 16.00 amu
     m_Mg   = 24.31*atomic_mass_unit,   &!Mg = 24.305 amu
     m_H2O  = 18.02*atomic_mass_unit,   &!Ox = 16.00  amu
     m_ir   = 55.845*atomic_mass_unit,  &
     m_Si   = 28.086*atomic_mass_unit,  &
     m_C2H2 = 26.038*atomic_mass_unit, &
     m_Si2C = 68.182*atomic_mass_unit, &
     m_carb = 12.01*atomic_mass_unit


contains

!-----------------------------------------------------------------------
!+
!  set particle dust properties (particle's call)
!+
!-----------------------------------------------------------------------
   subroutine evolve_condensation(dtsph, xyzh, u, dust_prop, rho)
      use dust_formation, only : wind_CO_ratio
      use units,  only:utime,unit_density
      use eos,    only:ieos,get_temperature
      use part,   only:icmu,icgamma,ickappa,ifol,ifqu,ifpy,ifir,ifsc,ifcarb,irol,irqu,irpy,irir,irsc,ircarb

      real,    intent(in) :: dtsph,rho,u,xyzh(4)
      real,    intent(inout) :: dust_prop(:)

      real :: dt_cgs, T, rho_cgs, vxyzui(4)

      dt_cgs    = dtsph* utime
      rho_cgs   = rho*unit_density
      vxyzui(4) = u
      T         = get_temperature(ieos,xyzh,rho,vxyzui,gammai=dust_prop(icgamma),mui=dust_prop(icmu))
      call dust_growth_condensation(T, rho_cgs, dt_cgs,wind_CO_ratio,&
           dust_prop(ifol),dust_prop(ifqu),dust_prop(ifpy),dust_prop(ifir),dust_prop(ifsc),dust_prop(ifcarb),&
           dust_prop(irol),dust_prop(irqu),dust_prop(irpy),dust_prop(irir),dust_prop(irsc),dust_prop(ircarb),&
           dust_prop(ickappa),dust_prop(icmu),dust_prop(icgamma))

   end subroutine evolve_condensation

!-----------------------------------------------------------------------
!+
!  set particle dust properties (particle's call)
!+
!-----------------------------------------------------------------------
   subroutine dust_growth_condensation(T,rho_cgs,dt,wind_CO_ratio,fol,fqu,fpy,fir,fsc,fcarb,&
        r_ol,r_qu,r_py,r_ir,r_sc,r_carb,kappa_dust,mu,gamma,abundance,pH_tot,pressure_cgs)
        use physcon, only:patm,kboltz,pi
        !use dust_formation, only:kappa_gas, mass_per_H, eps
        use chemistry_condensation, only: network, iSiO, iH2O, iH2, eps, &
             iH,iHe,iC,iOx,iN,iSi,iS,iFe,iTi,iMg, iSi2C, mass_per_H
        use dust_formation, only: kappa_gas,a_init_dust
        real, intent(in)          :: dt,wind_CO_ratio
        real, intent(in),optional :: pressure_cgs
        real, intent(inout)       :: T,rho_cgs,fol,fqu,fpy,fir,fsc,fcarb,r_ol,r_qu,r_py,r_ir,r_sc,r_carb
        real, intent(out)         :: kappa_dust
        real, intent(out)         :: mu,gamma
        real, intent(out),optional:: pH_tot,abundance(:)
        real :: Vth_SiO, Vth_Mg, Vth_H2O, Vth_ir, Vth_Si2C, Vth_carb
        real :: pv_ir, pv_SiO, pv_Mg, pv_H2O, pv_carb, pv_Si2C
        real :: Jgr_SiO_ol, Jgr_Mg_ol, Jgr_H2O_ol, Jgr_SiO_qu, Jgr_Mg_qu, Jgr_H2O_qu, Jgr_SiO_py, Jgr_Mg_py, Jgr_H2O_py
        real :: Jgr_ol,  Jgr_qu,  Jgr_py,  Jgr_ir,  Jgr_sc,  Jgr_carb
        real :: Jdec_ol, Jdec_qu, Jdec_py, Jdec_ir, Jdec_sc, Jdec_carb
        real :: kappa_ol,kappa_qu,kappa_py,kappa_ir,kappa_sc,kappa_carb
        real :: fol_dec, fqu_dec ,fpy_dec, fir_dec ,fsc_dec, fcarb_dec
        real :: P_H2, P_Si, G_const(3), root(3) !fol, kappa_dust
        real :: fol_max, fpy_max, fqu_max

        !XXX Aqui debo poner los valores de fol, fpy, fqu para restar las abundancias de Mg, H2O,
        if (present(abundance)) then
           if (present(pressure_cgs)) then
              call network(T,rho_cgs,mu,gamma,wind_CO_ratio,pH_tot, fol, fpy, fqu, fir, fsc, fcarb,abundance, pressure_cgs)
           else
              call network(T,rho_cgs,mu,gamma,wind_CO_ratio,pH_tot, fol, fpy, fqu, fir, fsc, fcarb,abundance)
           endif
        else
            call network(T,rho_cgs,mu,gamma,wind_CO_ratio,pH_tot, fol, fpy, fqu, fir, fsc, fcarb)
        endif

        Vth_SiO = sqrt(kboltz * T / 2.0 / pi / m_SiO) !already in cgs units
        Vth_Mg  = sqrt(kboltz * T / 2.0 / pi / m_Mg)
        Vth_H2O = sqrt(kboltz * T / 2.0 / pi / m_H2O)
        Vth_ir  = sqrt(kboltz * T / 2.0 / pi / m_ir)
        !Vth_Si = sqrt(kboltz * T / 2.0 / pi / m_Si)
        !Vth_C2H2 = sqrt(kboltz * T / 2.0 / pi / m_C2H2)
        Vth_Si2C = sqrt(kboltz * T / 2.0 / pi / m_Si2C)
        Vth_carb = sqrt(kboltz * T / 2.0 / pi / m_carb)

        Jgr_SiO_ol = alpha_ol * abundance(iSiO) * Vth_SiO !abundance comes from chemical equilibrium
        Jgr_Mg_ol  = alpha_ol * abundance(78)   * Vth_Mg !Mg = 78 in subroutine network
        Jgr_H2O_ol = alpha_ol * abundance(iH2O) * Vth_H2O

        Jgr_SiO_qu = alpha_qu * abundance(iSiO) * Vth_SiO !abundance comes from chemical equilibrium
        !Jgr_Mg_qu  = alpha_qu * abundance(78)   * Vth_Mg !Mg = 78 in subroutine network
        Jgr_H2O_qu = alpha_qu * abundance(iH2O) * Vth_H2O

        Jgr_SiO_py = alpha_py * abundance(iSiO) * Vth_SiO !abundance comes from chemical equilibrium
        Jgr_Mg_py  = alpha_py * abundance(78)   * Vth_Mg !Mg = 78 in subroutine network
        Jgr_H2O_py = alpha_py * abundance(iH2O) * Vth_H2O

        !Jgr_Si_sc = alpha_sc * abundance(75) * Vth_Si
        !Jgr_C2H2_sc = alpha_sc * abundance(iC2H2) * Vth_C2H2

        Jgr_ol = min(Jgr_SiO_ol, min(0.5*Jgr_Mg_ol, 0.333*Jgr_H2O_ol))
        Jgr_py = min(Jgr_SiO_py, min(Jgr_Mg_py, 0.5*Jgr_H2O_py))
        Jgr_qu = min(Jgr_SiO_qu, Jgr_H2O_qu)
        Jgr_ir = alpha_ir * abundance(77) * Vth_ir !Fe = 77 in subroutine network
        !Jgr_sc = min(Jgr_Si_sc, Jgr_C2H2_sc)
        Jgr_sc = alpha_sc * abundance(iSi2C) * Vth_Si2C
        Jgr_carb = alpha_carb * abundance(72) *Vth_carb !C = 72 in subroutine network

        !pH_tot = pH_tot !!!!%%%%  *patm !now in cgs
        P_H2 = abundance(iH2)*kboltz*T/patm !H2 pressure from abundance
        P_Si = abundance(75)*kboltz*T/patm !Si = 75 in subroutine network

        !%% Constant term in the law of mass action equation
        G_const(1) = P_H2**3 / (calc_Kp(coefficients(:,1),T) * pH_tot**6) !for olivine
        G_const(2) = P_H2    / (calc_Kp(coefficients(:,2),T) * pH_tot**2) !for Quartz
        G_const(3) = P_H2**2 / (calc_Kp(coefficients(:,3),T) * pH_tot**4) !for Pyroxene

        !%% to determine degree of condensation fol, fqu, fpy
!LS : this test should be removed, check must be made in init since wind_CO_ratio is constant
        if (wind_CO_ratio <= 0.9) then
           ! to prevent overflow
           if (T > 3500.) then
              root = 0.
           else
              call find_root(eps(iMg), eps(iSi), eps(iOx), eps(iC), G_const, root)
           endif
        elseif (wind_CO_ratio >= 1.1) then
            root = 0.
         endif

        fol_dec   = max(0., root(1)) !min(1., max(0., root(1)))
        fqu_dec   = max(0., root(2)) !min(1., max(0., root(2)))
        fpy_dec   = max(0., root(3)) !min(1., max(0., root(3)))
        fir_dec   = min(1., max(0., 1.- 1./(calc_Kp(coefficients(:,4),T) * eps(iFe) * pH_tot)))     !For Iron
        fsc_dec   = min(1., max(0., 1.- P_Si / (calc_Kp(coefficients(:,5),T) * eps(iSi) * pH_tot))) !For SiC
        fcarb_dec = min(1., max(0., 1.-1./wind_CO_ratio-2.*P_H2/(calc_Kp(coefficients(:,6),T)*eps(iC)*pH_tot)))

        !%% to calculate Quartz degree of condensation fqu
        !fqu_dec = solve_q(eps(iSi)**2, &
        !      - eps(iSi)*eps(iOx) + eps(iC)*eps(iSi), &
        !      + eps(iSi)*eps(iOx) - eps(iC)*eps(iSi) - eps(iSi)**2 &
        !      - P_H2 / (calc_Kp(coefficients(:,2), T)*pH_tot**2))

        if (Jgr_ol == Jgr_SiO_ol) then ! JgrSiO < JgrMg and JgrSiO < JgrH2O
            pv_SiO = max(0., (1.-fol_dec) * eps(iSi) * pH_tot) !Check units
            Jdec_ol = alpha_ol * Vth_SiO * pv_SiO *patm / kboltz / T !patm to have cgs units
        elseif (Jgr_ol == 0.5*Jgr_Mg_ol) then
            pv_Mg = max(0., (eps(iMg) - 2.0*fol_dec*eps(iSi))*pH_tot)
            Jdec_ol = alpha_ol * Vth_Mg * pv_Mg *patm / kboltz / T !patm to have cgs units
        elseif (Jgr_ol == 0.333*Jgr_H2O_ol) then
            pv_H2O = max(0.,(eps(iOx)-eps(iC)-(1.0+3.0*fol_dec)*eps(iSi))*pH_tot)
            Jdec_ol = alpha_ol * Vth_H2O * pv_H2O * patm / kboltz / T  !patm to have cgs units
        endif

        if (Jgr_qu == Jgr_SiO_qu) then
            pv_SiO = max(0., (1.-fqu_dec)*eps(iSi)* pH_tot)
            Jdec_qu = alpha_qu * Vth_SiO * pv_SiO *patm / kboltz / T !!patm to have cgs units
        elseif (Jgr_qu == Jgr_H2O_qu) then
            pv_H2O = max(0.,(eps(iOx)-eps(iC)-(1.0+fqu_dec)*eps(iSi))*pH_tot)
            Jdec_qu = alpha_qu * Vth_H2O * pv_H2O *patm / kboltz / T  !patm to have cgs units
        endif

        if (Jgr_py == Jgr_SiO_py) then
            pv_SiO = max(0.,(1.-fpy_dec)*eps(iSi)* pH_tot)
            Jdec_py = alpha_py * Vth_SiO * pv_SiO *patm / kboltz / T  !!patm to have cgs units
        elseif (Jgr_py == Jgr_Mg_py) then
            pv_Mg = max(0., (eps(iMg) - fpy_dec*eps(iSi))*pH_tot)
            Jdec_py = alpha_py * Vth_Mg * pv_Mg *patm / kboltz / T   !patm to have cgs units
        elseif (Jgr_py == 0.5*Jgr_H2O_py) then
            pv_H2O = max(0., (eps(iOx)-eps(iC)-(1.0+2.0*fpy_dec)*eps(iSi))*pH_tot)
            Jdec_py = alpha_py * Vth_H2O * pv_H2O *patm / kboltz / T  !patm to have cgs units
        endif

        pv_ir   = max(0.,  (1.-fir_dec)  * eps(iFe) * pH_tot)
        pv_Si2C = max(0.,  (1.-fsc_dec)  * eps(iSi) * pH_tot)  !0.5 *
        pv_carb = max(0., ((1.-fcarb_dec)* eps(iC)-eps(iOx)) * pH_tot)

        Jdec_ir = alpha_ir * Vth_ir * pv_ir *patm /  kboltz / T   !!patm to have cgs units
        Jdec_sc = 2.0* alpha_sc * Vth_Si2C * pv_Si2C *patm / kboltz / T  !eqn. 21 of Ferrarotti & Gail 2002
        Jdec_carb = alpha_carb * Vth_carb * pv_carb *patm / kboltz / T


        r_ol = max(r_ol + Vo_ol * (Jgr_ol-Jdec_ol) * dt , a_init_dust)
        r_qu = max(r_qu + Vo_qu * (Jgr_qu-Jdec_qu) * dt , a_init_dust)
        r_py = max(r_py + Vo_py * (Jgr_py-Jdec_py) * dt , a_init_dust)
        r_ir = max(r_ir + Vo_ir * (Jgr_ir-Jdec_ir) * dt , a_init_dust)
        r_sc = max(r_sc + Vo_sc * (Jgr_sc-Jdec_sc) * dt , a_init_dust)
        r_carb = max(r_carb + Vo_carb * (Jgr_carb-Jdec_carb) * dt , a_init_dust)


        fol = max(0., 4.*pi*(r_ol**3-a_init_dust**3)*1.d-13/3./Vo_ol/eps(iSi)) !fol
        if (wind_CO_ratio<=0.9) then
           fol_max = min(1.0, min(0.5*eps(iMg)/eps(iSi),(eps(iOx)-eps(iC)-eps(iSi))/3./eps(iSi)))
        elseif (wind_CO_ratio>=1.1) then
            fol_max = 0.0
        endif        !%%fol_max is negative if C/O ratio is larger than 1
        if (fol > fol_max) fol = fol_max
!%%PUEDE SER QUE fol sea negativo porque fol_max es negativo

        fqu = max(0., 4.*pi*(r_qu**3-a_init_dust**3)*1.d-13/3./Vo_qu/eps(iSi)) !fqu
        if (wind_CO_ratio<=0.9) then
           fqu_max = min(1.0, (eps(iOx)-eps(iC)-eps(iSi))/eps(iSi))
        elseif (wind_CO_ratio>=1.1) then
           fqu_max = 0.0
        endif
        if (fqu > fqu_max) fqu = fqu_max

        fpy = max(0., 4.*pi*(r_py**3-a_init_dust**3)*1.d-13/3./Vo_py/eps(iSi)) !fpy
        if (wind_CO_ratio<=0.9) then
           fpy_max = min(1.0, min(eps(iMg)/eps(iSi),(eps(iOx)-eps(iC)-eps(iSi))/2./eps(iSi)))
        elseif (wind_CO_ratio>=1.1) then
            fpy_max = 0.0
        endif
        if (fpy > fpy_max) fpy = fpy_max

        fir = max(0., 4.*pi*(r_ir**3-a_init_dust**3)*1.d-13/3./Vo_ir/eps(iFe)) !fir
        if (fir > 1.0) fir = 1.

        fsc = max(0., 4.*pi*(r_sc**3-a_init_dust**3)*1.d-13/3./Vo_sc/ eps(iSi)) !fsc
        if (fsc > 1.0) fsc = 1.

        fcarb = max(0., 4.*pi*(r_carb**3-a_init_dust**3)*1.d-13/3./Vo_carb/eps(iC)) !fcarb
        if (fcarb > 1.0) fcarb = 1.


        kappa_ol = ((6.147d-07 * T**2.444)**(-2) &
                 + 1.0/((6.957d4 * T**(-2.329))**2  &
                 + sqrt((3.505d-4 * T**0.755)**4 + (1.043d-9 * T**2.523)**4 )))**(-0.5)

        kappa_qu = ( ( 1./(3.898d-6 * T**2.054)**4 + 1./(3.1d5 * T**(-2.622))**4 )**(-0.5) &
                    + ( (2.023d-5 * T**1.074)**4 + (9.394d-11 * T**2.701)**4 )**0.5 )**0.5


        kappa_py = ( (3.773d-5 * T**2.017)**(-2) &
                   + 1.0/( (1.5d6 * T**(-2.679))**2 &
                   + sqrt((8.356d-4 * T**0.7336)**4 + (1.261d-8 * T**2.272 )**4)))**(-0.5)

        kappa_ir = ( (3.341d-5 * T**1.632)**(-4) &
                   + 1.0/ ( (6.405d-4 * T**0.777)**4 + (4.385d-7 * T**1.981)**4 )  )**(-0.25)

        kappa_sc = ( (9.56d-6 * T**1.71)**(-4) &
                   + 1.0 / ( (7.791d-4 * T**0.8808)**4 + (1.171d-6 * T**1.879)**4) )**(-0.25)

        kappa_carb = 5.9 * T * A_carb * eps(iC) * pi / (2.2 * (mass_per_H/atomic_mass_unit))
                    !5.6d-18 * pi * T / rho_cgs !%%% The correct expression for opacity

                    !5.9**(-25.13) * T * pi * A_carb * eps(iC) * (pH_tot*patm/kboltz/T) &
                   !* 2.0*(sqrt(2.5d-5)-sqrt(5.0d-7)) / rho_carb / (mass_per_H/atomic_mass_unit)

        !%% mu = mass_per_H / atomic_mass_unit

        kappa_dust = kappa_gas + fol * kappa_ol + fqu * kappa_qu &
                   + fpy * kappa_py + fir * kappa_ir + fsc * kappa_sc + fcarb * kappa_carb

     end subroutine dust_growth_condensation

subroutine find_root(m, s, o, c, g, root) !, tol, max_iter)
    implicit none
    ! Input parameters
    real, intent(in) :: m, s, o, g(3), c !m = eps(iMg), s = eps(iSi),o = eps(iOx), c = eps(iC)
    real, intent(out) :: root(3)
    real :: tol = 1.0d-6!, intent(in) :: tol
    integer :: max_iter = 1000 !, intent(in) :: max_iter

    ! Local variables
    real :: x(3), fx(3), dfx(3), A(3), B(3), D(3)
    integer :: i

    ! Calculate the interval boundaries
    real :: x_max(3)

    !%only valid when C/O is less than 1, otherwise there is no Silicate dust formation, all Oxygen is locked in CO

    x_max(1) = min(1.0d0, min(5.0d-1*m/s, 1.0d0/3.0d0 * (o-c-s)/s) ) !Max value for Olivine
    x_max(2) = min(1.0d0, (o-c-s)/s) !Max value for Quartz
    x_max(3) = min(1.0d0, min(m/s, 5.0d-1*(o-c-s)/s )) !Max value for Pyroxene

    x(1) = x_max(1)-0.001 ! Initial guess given by the maximum available degree of condensation
    x(2) = x_max(2)-0.001
    x(3) = x_max(3)-0.001

    ! Newton-Raphson iteration
    do i = 1, max_iter
        !%% Calculate A, B, and D for Olivine
        A(1) = m - 2. * x(1) * s
        B(1) = 1.0d0 - x(1)
        D(1) = o - c - (1.0d0 + 3.0d0 * x(1)) * s

        !%% Calculate A, B, and D for Quartz
        A(2) = 1.0d0
        B(2) = 1.0d0 - x(2)
        D(2) = o -c - (1.0d0 + x(2)) * s

        !%% Calculate A, B, and D for Pyroxene
        A(3) = m - x(3) * s
        B(3) = 1.0d0 - x(3)
        D(3) = o - c - (1.0d0 + 2.0d0 * x(3)) * s

        !%% Calculate f(x) and f'(x) for Olivine
        fx(1) = -g(1) + A(1)**2 * B(1) * s * D(1)**3
        dfx(1) = -s * A(1) * D(1)**2 * (4.0d0 * B(1) * D(1) * s + 9.0d0 * s * B(1) * A(1) + A(1) * D(1))

        !%% Calculate f(x) and f'(x) for Quartz
        fx(2) = -g(2) + A(2) * B(2) * s * D(2)
        dfx(2) = -s**2 * B(2) -s * D(2)

        !%% Calculate f(x) and f'(x) for Pyroxene
        fx(3) = -g(3) + A(3) * B(3) * s * D(3)**2
        dfx(3) = -s * D(3) *(4.0d0 * A(3) * B(3) + s * B(3) * D(3) + A(3) * D(3))

        ! Check for zero derivative
        if (ANY(dfx(1:3) == 0.0d0)) then
            print *, "Derivative is zero. No solution found."
            return
        end if

        ! Update x using Newton-Raphson formula
        x(:) = x(:) - fx(:) / dfx(:)

        ! Check for convergence
        if (all(abs(fx(:)) < tol)) then
            root(:) = x(:)
            return
        end if
    end do

    print *, "Maximum iterations reached. No solution found."
    root(1) = x(1)
    root(2) = x(2)
end subroutine find_root

pure real function calc_Kp(coefficients,T)
! all quantities are in cgs
real, intent(in) :: coefficients(5), T
real, parameter :: R = 1.987165
real :: G, d
G = coefficients(1)/T + coefficients(2) + (coefficients(3)+(coefficients(4)+coefficients(5)*T)*T)*T
d = min(-G/(R*T),222.)
calc_Kp = exp(d)
end function calc_kp


end module dust_condensation
