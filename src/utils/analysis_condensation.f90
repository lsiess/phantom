!-----------------------!
!File Naming with itoa Function:
!A helper function itoa is used to convert the timestep 
!integer to a string, which is then appended to the file 
!names to make them unique for each time step.
!Updating Time and Timestep Count:
!After completing the inner loop over temperatures, 
!time_current is updated by adding dt, and timestep is incremented.
!This structure allows you to generate new output files 
!for each time step, capturing the evolution of the system over time.

!Parallel Execution with OpenMP:

!The outer loop over timestep is parallelized using OpenMP with the !$omp parallel and !$omp do directives.
!The private clause ensures each thread has its own copy of specified variables, preventing race conditions.
!The shared clause lists variables that are shared among all threads.

!Dynamic Scheduling:

!The schedule(dynamic) clause dynamically assigns iterations to threads to help balance the workload.

module testchemistry
  implicit none

contains

  !Daniel-test subroutine test_chemistry(ntests, npass, pressure_cgs)
  subroutine test_chemistry(pressure_cgs)
    use chemistry, only: ncols, write_time_file, columns !network
    use dust_growth, only: O_rich_nucleation, condensation
    !Daniel-test integer, intent(inout) :: ntests, npass
    real :: T, rho_cgs, mu, gamma, abundance(ncols), wind_CO_ratio, pH_tot, dt
    real, allocatable :: fol(:), kappa_dust(:), fqu(:), fpy(:), fir(:), fsc(:), fcarb(:)
    integer :: num, i, timestep, size_array
    real :: dust_properties(7)
    real :: time_0, time_final, time_current, rho_0, T0
    real, optional :: pressure_cgs

    !Daniel-test ntests = ntests + 1
    !Daniel-test npass = npass + 1

    !----XXXX----!
    !We condsider a wind with a density distribution rho(T) = rho_0 * (T/T0) **2
    !where rho_0 and T0 are the density and temperature at the base of the wind


    T0 = 3.0d3
    rho_0 = 1.0d-10
    !pressure_cgs = 1.d-4  !optional variable for comparison with Gail & Sedlmayr calculations

    ! Set initial parameters
    time_0 = 0.0
    time_final = 1.0d3 !1.0d3  ! or whatever final time you want
    dt = 5.0d0 !5.0d2          ! time step size
    wind_CO_ratio = 0.34 !C/O can be larger than or lower than 1.0
    size_array = int((time_final - time_0) / dt) + 1
    allocate(fol(size_array),fpy(size_array),fqu(size_array),fir(size_array), &
             fsc(size_array), fcarb(size_array),kappa_dust(size_array))
    fol = 0.0    !degree of condensation of olivine
    fpy = 0.0    !degree of condensation of pyroxene
    fqu = 0.0    !degree of condensation of quartz
    fir = 0.0    !degree of condensation of solid iron
    fsc = 0.0    !degree of condensation of silicon carbide
    fcarb = 0.0  !degree of condensation of solid carbon (graphite)
    dust_properties = 0.

    !! Parallelize the outer loop (temperature loop) using OpenMP
     !$omp parallel private(i, T, rho_cgs, mu, gamma, abundance, pH_tot, dust_properties, num) &
     !$omp shared(fol, fpy, fqu, fir, kappa_dust, time_0, time_final, dt, wind_CO_ratio, T0, rho_0) 
     !$omp do schedule(dynamic)

    do timestep = 1, int((time_final - time_0) / dt) + 1
      time_current = time_0 + (timestep - 1) * dt
      num = 0

      ! Loop over temperature T from 600 to 3000
      do i = 60, 300 !60
        T = i * 10.0
        rho_cgs = rho_0 * (T / T0)**2

        if (present(pressure_cgs)) then
         call O_rich_nucleation(T,rho_cgs,dt,wind_CO_ratio,fol(timestep),fqu(timestep),fpy(timestep),fir(timestep), &
                  fsc(timestep), fcarb(timestep),kappa_dust(timestep), mu, gamma, pH_tot, abundance, pressure_cgs)
        else
         call O_rich_nucleation(T,rho_cgs,dt,wind_CO_ratio,fol(timestep),fqu(timestep),fpy(timestep),fir(timestep), &
                                fsc(timestep), fcarb(timestep), kappa_dust(timestep), mu, gamma, pH_tot, abundance)
        endif
        
        dust_properties(1) = dust_properties(1) + fol(timestep) != fol(timestep) ONLY 
        dust_properties(2) = dust_properties(2) + fqu(timestep)
        dust_properties(3) = dust_properties(3) + fpy(timestep)
        dust_properties(4) = dust_properties(4) + fir(timestep)
        dust_properties(5) = dust_properties(5) + fsc(timestep)
        dust_properties(6) = dust_properties(6) + fcarb(timestep)
        dust_properties(7) = dust_properties(7) + kappa_dust(timestep)

        ! Ensure file writing is thread-safe using a critical section
        !$omp critical
        call write_time_file('abundances_'//trim(adjustl(itoa(timestep))),columns,T,abundance,ncols,num)
        call write_time_file('condensation_'//trim(adjustl(itoa(timestep))),condensation,T,dust_properties,7,num)
        !$omp end critical

        num = num + 1
      end do
    end do

     !$omp end do
     !$omp end parallel

    deallocate(fol,fpy,fqu,fir,fsc, fcarb,kappa_dust)

  end subroutine test_chemistry

  function itoa(n) result(str)
    integer, intent(in) :: n
    character(len=32) :: str

    write(str, '(I0)') n
    str = adjustl(str)
  end function itoa

end module testchemistry