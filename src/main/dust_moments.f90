!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module dust_moments
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

 use physcon, only:kboltz,pi,atomic_mass_unit,mass_proton_cgs,patm

 implicit none

 public :: evolve_moments,evolve_chem,calc_nucleation

 real, public, parameter :: Scrit = 2. ! Critical saturation ratio

 private
 real, parameter :: vfactor = sqrt(kboltz/(2.*pi*atomic_mass_unit*12.01))
 !real, parameter :: vfactor = sqrt(kboltz/(8.*pi*atomic_mass_unit*12.01))

contains

!-----------------------------------------------------------------------
!
!  set particle dust properties (particle's call)
!
!-----------------------------------------------------------------------
subroutine evolve_moments(dtsph, xyzh, u, JKmuS, Tdust, rho)
 use units,  only:utime,unit_density
 use eos,    only:ieos,get_temperature
 use part,   only:idK3,idmu,idgamma,idkappa
 use dust_formation, only:calc_kappa_dust

 real,    intent(in) :: dtsph,Tdust,rho,u,xyzh(4)
 real,    intent(inout) :: JKmuS(:)

 real :: dt_cgs, T, rho_cgs, vxyzui(4)

 dt_cgs    = dtsph* utime
 rho_cgs   = rho*unit_density
 vxyzui(4) = u
 T         = get_temperature(ieos,xyzh,rho,vxyzui,gammai=JKmuS(idgamma),mui=JKmuS(idmu))
 call evolve_chem(dt_cgs, T, rho_cgs, JKmuS)
 JKmuS(idkappa)     = calc_kappa_dust(JKmuS(idK3), Tdust, rho_cgs)

end subroutine evolve_moments

!-----------------------------------------------------------------------
!
!  evolve the chemistry and moments
!
!-----------------------------------------------------------------------
subroutine evolve_chem(dt, T, rho_cgs, JKmuS)
!all quantities in cgs

 use chemistry_moments, only:mass_per_H,iC,iOx,eps,iHe,psat_C,chemical_equilibrium_light
 use part,           only:idJstar,idK0,idK1,idK2,idK3,idmu,idgamma,idsat
 real, intent(in)    :: dt, rho_cgs
 real, intent(inout) :: T, JKmuS(:)

 real :: pC, pC2, pC2H, pC2H2, nH_tot, epsC
 real :: JstarS, taustar, taugr, S
 real :: Jstar_new, K_new(0:3)
 real :: nH, nH2, v1
 real, parameter :: A0 = 20.7d-16
 real, parameter :: alpha2 = 0.34

 nH_tot = rho_cgs/mass_per_H
 epsC   = eps(iC) - JKmuS(idK3)
 if (epsC < 0.) then
    print *,'eps(C) =',eps(iC),', K3=',JKmuS(idK3),', epsC=',epsC,', T=',T,' rho=',rho_cgs
    print *,'JKmuS=',JKmuS
    stop '[S-dust_formation] epsC < 0!'
 endif
 if (T > 450.) then
    call chemical_equilibrium_light(rho_cgs, T, epsC, pC, pC2, pC2H, pC2H2, JKmuS(idmu), JKmuS(idgamma))
    S = pC/psat_C(T)
    if (S > Scrit) then
       !call nucleation(T, pC, pC2, 0., pC2H, pC2H2, S, JstarS, taustar, taugr)
       call calc_nucleation(T, pC, 0., 0., 0., pC2H2, S, JstarS, taustar, taugr)
       JstarS = JstarS/ nH_tot
       call evol_K(JKmuS(idJstar), JKmuS(idK0:idK3), JstarS, taustar, taugr, dt, Jstar_new, K_new)
    else
       Jstar_new  = JKmuS(idJstar)
       K_new(0:3) = JKmuS(idK0:idK3)
    endif
 else
! Simplified low-temperature chemistry: all hydrogen in H2 molecules, all O in CO
    nH  = 0.
    nH2 = nH_tot/2.
    JKmuS(idmu)    = (1.+4.*eps(iHe))*nH_tot/(nH+nH2+eps(iHe)*nH_tot)
    JKmuS(idgamma) = (5.*eps(iHe)+3.5)/(3.*eps(iHe)+2.5)
    pC2H2 = .5*(epsC-eps(iOx))*nH_tot * kboltz * T
    pC2H  = 0.
    S     = 1.d-3
    v1    = vfactor*sqrt(T)
    taugr = kboltz*T/(A0*v1*sqrt(2.)*alpha2*(pC2H+pC2H2))
    call evol_K(0., JKmuS(idK0:idK3), 0., 1., taugr, dt, Jstar_new, K_new)
 endif
 JKmuS(idJstar)   = Jstar_new
 JKmuS(idK0:idK3) = K_new(0:3)
 JKmuS(idsat)     = S

end subroutine evolve_chem

!----------------------------
!
!  Compute nucleation rate
!
!----------------------------
subroutine calc_nucleation(T, pC, pC2, pC3, pC2H, pC2H2, S, JstarS, taustar, taugr)
! all quantities are in cgs
 real, intent(in)  :: T, pC, pC2, pC3, pC2H, pC2H2, S
 real, intent(out) :: JstarS, taustar, taugr
 real, parameter   :: A0 = 20.7d-16
 real, parameter   :: sigma = 1400.
 real, parameter   :: theta_inf = A0*sigma/kboltz
 real, parameter   :: alpha1 = 0.37 !sticking coef for C
 real, parameter   :: alpha2 = 0.34 !sticking coef for C2,C2H,C2H2
 real, parameter   :: alpha3 = 0.08 !sticking coef for C3
 real, parameter   :: Nl_13 = 5.**(1./3.)
 real, parameter   :: mproton = 1.6726485d-24

 real :: ln_S_g, Nstar_inf_13, Nstar_m1_13, Nstar_m1_23, theta_Nstar,Nstar_m1,Nstar_tmp
 real :: dtheta_dNstar, d2lnc_dN2star, Z, A_Nstar, v1, beta, c_star, expon

 ln_S_g = log(S)
 v1     = vfactor*sqrt(T)
 Nstar_inf_13  = 2.*theta_inf/(3.*T*ln_S_g)
 Nstar_tmp     = 1. + sqrt(1. + 2.*Nl_13/Nstar_inf_13) - 2.*Nl_13/Nstar_inf_13
 if (Nstar_tmp > 0.) then
    Nstar_m1      = (Nstar_inf_13**3)/8. * Nstar_tmp**3
    Nstar_m1_13   = 0.5*Nstar_inf_13*Nstar_tmp
    Nstar_m1_23   = Nstar_m1_13**2
    theta_Nstar   = theta_inf/(1. + Nl_13/Nstar_m1_13)
    dtheta_dNstar = Nl_13 * theta_Nstar**2/(3.*theta_inf * Nstar_m1_23**2)
    d2lnc_dN2star = 2.*theta_Nstar/(9.*T*Nstar_m1)*(Nstar_m1_13+2.*Nl_13)/(Nstar_m1_13+Nl_13)**2
    Z             = sqrt(d2lnc_dN2star/(2.*pi))
    A_Nstar       = A0 * (1.+Nstar_m1)**(2./3.)
    beta          = v1/(kboltz*T) * (pC*alpha1 + 4.*alpha2/sqrt(2.)*(pC2 + pC2H + pC2H2) + 9.*alpha3/sqrt(3.)*pC3)
    expon         = Nstar_m1*ln_S_g - theta_Nstar*Nstar_m1_23/T
    if (expon < -100.) then
       c_star = 1.d-70
    else
       c_star = pC/(kboltz*T) * exp(expon)
    endif
    JstarS  = beta * A_Nstar * Z * c_star
    taustar = 1./(d2lnc_dN2star*beta*A_Nstar)
    ! if (isnan(JstarS)) then
    !   print*,i,'(N-1)^1/3=',Nstar_m1_13,'exp=',expon,'T=',T,'theta_N=',theta_Nstar,'d2lnc/dN2=',d2lnc_dN2star,ddd,&
    !        'beta=',beta,'Z=',Z,'c_star=',c_star,'JstarS=',JstarS,'tau*=',taustar
    !   if (isnan(JstarS)) stop
    ! endif
 else
    JstarS  = 0.d0
    taustar = 1.d-30
 endif
 taugr = kboltz*T/(A0*v1*(alpha1*pC*(1.-1./S) + 2.*alpha2/sqrt(2.)*(pC2+pC2H+pC2H2)*(1.-1./S**2)))
end subroutine calc_nucleation

!------------------------------------
!
!  Compute evolution of the moments
!
!------------------------------------
subroutine evol_K(Jstar, K, JstarS, taustar, taugr, dt, Jstar_new, K_new)
! all quantities are in cgs, K and K_new are the *normalized* moments (K/n<H>)
 real, intent(in) :: Jstar, K(0:3), JstarS, taustar, taugr, dt
 real, intent(out) :: Jstar_new, K_new(0:3)

 real, parameter :: Nl_13 = 10. !(lower grain size limit)**1/3
 real :: d, i0, i1, i2, i3, i4, i5, dK0, dK1, dK2, DK3

 d = dt/taustar
 if (d > 500.) then
    i0 = 0.
 else
    i0 = exp(-d)
 endif
 i1 = 1. - i0
 i2 = d - i1
 i3 = d**2/2. - i2
 i4 = d**3/6. - i3
 i5 = d**4/24. - i4
 Jstar_new = Jstar*i0 + JstarS*i1
 dK0 = taustar*(Jstar*i1 + JstarS*i2)
 K_new(0) = K(0) + dK0
 dK1 = taustar**2/(3.*taugr)*(Jstar*i2 + JstarS*i3)
 K_new(1) = K(1) + K(0)*dt/(3.*taugr) + dK1 + Nl_13*dK0
 dK2 = 2.*taustar**3/(3.*taugr)**2 * (Jstar*i3 + JstarS*i4)
 K_new(2) = K(2) + 2.*K(1)*dt/(3.*taugr) + (dt/(3.*taugr))**2*K(0) + dK2 + 2.*Nl_13*dK1 + Nl_13**2*dK0
 dK3 = 3.*dt/(3.*taugr)*K(2) + 3.*(dt/(3.*taugr))**2*K(1) + (dt/(3.*taugr))**3*K(0)  &
     + (6.*taustar**4)/(3.*taugr)**3*(Jstar*i4+JstarS*i5)
 K_new(3) = K(3) + dK3 + Nl_13**3*dK0 + 3.*Nl_13**2*dK1 + 3.*Nl_13*dK2
 !if (any(isnan(K_new))) then
 !  print*,'NaNs in K_new for particle #',i
 !  print *,'dt=',dt,'tau*=',taustar,'taug=',taugr,'d=',d,'i0=',i0,'i1=',i1,'Jstar=',Jstar,'JstarS=',JstarS,&
 !      'k1=',k(1),'dk1=',dk1,'Kn1=',k_new(1),'k2=',k(2),'dk2=',dk2,'Kn2=',k_new(2),'k3=',k(3),'dk3=',dk3,'Kn3=',k_new(3)
 !  stop
 !endif

end subroutine evol_K

end module dust_moments
