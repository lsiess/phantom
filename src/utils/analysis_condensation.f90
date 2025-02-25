!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis

 use chemistry_condensation, only : ncols
 implicit none
 character(len=20), parameter, public :: analysistype = 'condenstation'
 public :: do_analysis

 private

 character(len=*), parameter :: condensation(7) =&
                 (/'fol              ', &
                   'fqu              ', &
                   'fpy              ', &
                   'fir              ', &
                   'fsc              ', &
                   'fcarb            ', &
                   'kappa_dust       '/)

!The vector columns include molecules and atoms
 character(len=10), parameter :: columns(ncols) =&
                 (/'nH2       ', &
                   'nOH       ', &
                   'nH2O      ', &
                   'nCO       ', &
                   'nCO2      ', &
                   'nCH4      ', &
                   'nC2H      ', &
                   'nC2H2     ', &
                   'nN2       ', &
                   'nNH3      ', &
                   'nCN       ', &
                   'nHCN      ', &
                   'nSi2      ', &
                   'nSi3      ', &
                   'nSiO      ', &
                   'nSi2C     ', &
                   'nSiH4     ', &
                   'nS2       ', &
                   'nHS       ', &
                   'nH2S      ', &
                   'nSiS      ', &
                   'nSiH      ', &
                   'nTiO      ', &
                   'nTiO2     ', &
                   'nC2       ', &
                   'nO2       ', &
                   'nCH       ', &
                   'nCOH      ', &
                   'nC20      ', &
                   'nCH2      ', &
                   'nH2CO     ', &
                   'nCH3      ', &
                   'nC2H4     ', &
                   'nNH       ', &
                   'nNO       ', &
                   'nNCO      ', &
                   'nHCNO     ', &
                   'nC2N      ', &
                   'nC2N2     ', &
                   'nHNO      ', &
                   'nHNO2     ', &
                   'nHNO3     ', &
                   'nNH2      ', &
                   'nNO2      ', &
                   'nNO3      ', &
                   'nN20      ', &
                   'nN204     ', &
                   'nMgH      ', &
                   'nMgOH     ', &
                   'nMgO2H2   ', &
                   'nMgN      ', &
                   'nMgO      ', &
                   'nSiC      ', &
                   'nSiH2     ', &
                   'nSiH3     ', &
                   'nSiN      ', &
                   'nSiO2     ', &
                   'nFeO      ', &
                   'nFeO2H2   ', &
                   'nCOS      ', &
                   'nCS       ', &
                   'nCS2      ', &
                   'nFeS      ', &
                   'nH2SO4    ', &
                   'nMgS      ', &
                   'nSN       ', &
                   'nSO       ', &
                   'nSO2      ', &
                   'nSO3      ', &
                   'nH        ', &
                   'nHe       ', &
                   'nC        ', &
                   'nO        ', &
                   'nN        ', &
                   'nSi       ', &
                   'nS        ', &
                   'nFe       ', &
                   'nMg       '/)

contains

  subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
    use dim, only:do_nucleation,do_condensation
    use dust_condensation, only: dust_growth_condensation
    use dust_formation,    only: a_init_dust,wind_CO_ratio,set_abundances
    character(len=*), intent(in) :: dumpfile
    integer,          intent(in) :: num,npart,iunit
    real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
    real,             intent(in) :: particlemass,time

    real :: T, rho_cgs, mu, gamma, abundance(ncols), pH_tot, dt
    real, allocatable :: fol(:), fqu(:), fpy(:), fir(:), fsc(:), fcarb(:), kappa_dust(:)
    real, allocatable :: r_ol(:),r_qu(:),r_py(:),r_ir(:),r_sc(:),r_carb(:)
    integer :: i, k, timestep, nTstep, size_array
    real :: dust_properties(7)
    real :: time_0, time_final, rho_0, Tmin, Tmax, pressure_cgs

    do_condensation = .true.
    do_nucleation = .false.
    ! Set input values
    !We condsider a wind with a density distribution rho(T) = rho_0 * (T/Tmax) **2
    !where rho_0 and Tmax are the density and temperature at the base of the wind
    nTstep = 240  ! number of temperature grid point between Tmin and Tmax
    Tmin  = 600.
    Tmax  = 6000.
    rho_0 = 1.0d-10
    pressure_cgs = -1. !1.d-4  !optional variable for comparison with Gail & Sedlmayr calculations

    ! Set initial parameters
    time_0 = 0.
    time_final = 1.d5 !A larger value creates wiggles in fir and kappa
    dt = 5.0d4 !
    wind_CO_ratio = 2.0 !0.34

    call set_abundances

    size_array = int((time_final - time_0) / dt) + 1
    allocate(fol(size_array),fqu(size_array),fpy(size_array),fir(size_array), &
             fsc(size_array), fcarb(size_array),kappa_dust(size_array))
    allocate(r_ol(size_array),r_qu(size_array),r_py(size_array),r_ir(size_array), &
             r_sc(size_array), r_carb(size_array))

    fol = 0.    !degree of condensation of olivine
    fqu = 0.    !degree of condensation of quartz
    fpy = 0.    !degree of condensation of pyroxene
    fir = 0.    !degree of condensation of solid iron
    fsc = 0.    !degree of condensation of silicon carbide
    fcarb = 0.  !degree of condensation of solid carbon (graphite)
    r_ol = a_init_dust   !radius olivine
    r_qu = a_init_dust   !radius quartz
    r_py = a_init_dust   !radius pyroxene
    r_ir = a_init_dust   !radius solid iron
    r_sc = a_init_dust   !radius silicon carbide
    r_carb = a_init_dust !radius solid carbon (graphite)
    kappa_dust = 0.
    dust_properties = 0.

    do timestep = 1, size_array !int((time_final - time_0) / dt) + 1
      k = 0

      ! Loop over temperature
      do i = 1,nTstep
        T = Tmin + (i-1) * (Tmax-Tmin)/(nTstep-1)
        rho_cgs = rho_0 * (T / Tmax)**2

        if (pressure_cgs > 0.) then
         call dust_growth_condensation(T, rho_cgs, dt,wind_CO_ratio,&
              fol(timestep),fqu(timestep),fpy(timestep),fir(timestep),fsc(timestep),fcarb(timestep),&
              r_ol(timestep),r_qu(timestep),r_py(timestep),r_ir(timestep),r_sc(timestep),r_carb(timestep),&
              kappa_dust(timestep), mu, gamma, abundance, pH_tot, pressure_cgs)
        else
         call dust_growth_condensation(T, rho_cgs, dt,wind_CO_ratio,&
              fol(timestep),fqu(timestep),fpy(timestep),fir(timestep),fsc(timestep),fcarb(timestep),&
              r_ol(timestep),r_qu(timestep),r_py(timestep),r_ir(timestep),r_sc(timestep),r_carb(timestep),&
              kappa_dust(timestep), mu, gamma, abundance, pH_tot)
        endif

        dust_properties(1) = fol(timestep)
        dust_properties(2) = fqu(timestep)
        dust_properties(3) = fpy(timestep)
        dust_properties(4) = fir(timestep)
        dust_properties(5) = fsc(timestep)
        dust_properties(6) = fcarb(timestep)
        dust_properties(7) = kappa_dust(timestep)

        call write_time_file('abundances_' // trim(adjustl(itoa(timestep))), columns, T, abundance, ncols, k)
        call write_time_file('condensation_' // trim(adjustl(itoa(timestep))), condensation, T, dust_properties, 7, k)

        k = k + 1
      end do
    end do

    deallocate(fol,fpy,fqu,fir,fsc, fcarb,r_ol,r_py,r_qu,r_ir,r_sc, r_carb,kappa_dust)

 end subroutine do_analysis

  function itoa(n) result(str)
    integer, intent(in) :: n
    character(len=32) :: str

    write(str, '(I0)') n
    str = adjustl(str)
  end function itoa

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

    open(unit=unitnum, file=file_name, status='replace')
    do i=1,ncols
       write(columns(i), "(I2,a)") i+1, cols(i)
    enddo

    !set column headings
    write(unitnum, column_formatter) '1  Temperature', columns(:)
    close(unit=unitnum)
 endif

 unitnum=1001+num

 open(unit=unitnum, file=file_name, position='append')
 write(unitnum,data_formatter) time, data_in(:ncols)
 close(unit=unitnum)

end subroutine write_time_file

end module analysis
