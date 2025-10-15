!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module raytracer
!
! This module contains all routines required to:
!   - perform radial ray tracing starting from the primary star only
!   - calculate optical depth along the rays given the opacity distribution
!   - interpolate optical depths to all SPH particles
! Applicable both for single and binary star wind simulations
!
! WARNING: This module has only been tested on phantom wind setup
!
! :References: Esseldeurs M., Siess L. et al, 2023, A&A, 674, A122
!
! :Owner: Mats Esseldeurs
!
! :Runtime parameters: None
!
! :Dependencies: dim, healpix, kernel, linklist, part, units
!
 use healpix

 implicit none
 public :: get_all_integrands

 private

contains

 !------------------------------------------------------------------------------------
 !+
 !  MAIN ROUTINE
 !  Returns the optical depth at each particle's location using an outward ray-tracing scheme
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: nptmass:         The number of sink particles
 !  IN: xyzm_ptmass:     The array containing the properties of the sink particle
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa_cgs:       The array containing the opacities of all SPH particles
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !  IN: type:            The type of integrands to use (tau, tau_lucy, column density) in a list as e.g. [.true., .false., .false.]
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !------------------------------------------------------------------------------------
subroutine get_all_integrands(npart, nptmass, xyzmh_ptmass, xyzh, kappa_cgs, order, tau, tau_lucy, column_density, type)
 use part,   only: iReff
 integer, intent(in) :: npart, order, nptmass
 real, intent(in)    :: kappa_cgs(:), xyzh(:,:), xyzmh_ptmass(:,:)
 real, intent(out)   :: tau(:), tau_lucy(:), column_density(:)
 real, dimension(npart)   :: tau_temp, tau_lucy_temp, column_density_temp
 real :: Rinject
 logical, dimension(3), intent(in) :: type


 Rinject = xyzmh_ptmass(iReff,1)
 if (nptmass == 2 ) then
    call get_all_integrands_companion(npart, xyzmh_ptmass(1:3,1), xyzmh_ptmass(iReff,1), xyzh, kappa_cgs, &
            Rinject, xyzmh_ptmass(1:3,2), xyzmh_ptmass(iReff,2), order, tau_temp, tau_lucy_temp, column_density_temp, type)
 else
    call get_all_integrands_single(npart, xyzmh_ptmass(1:3,1), xyzmh_ptmass(iReff,1), xyzh,&
         kappa_cgs, Rinject, order, tau_temp, tau_lucy_temp, column_density_temp, type)
 endif

 if (type(1)) tau(:npart) = tau_temp
 if (type(2)) tau_lucy(:npart) = tau_lucy_temp
 if (type(3)) column_density(:npart) = column_density_temp

end subroutine get_all_integrands

 !---------------------------------------------------------------------------------
 !+
 !  Calculates the optical depth at each particle's location, using the uniform outward
 !  ray-tracing scheme for models containing a single star
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the kappa of all SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !  IN: type:            The type of integrands to use (tau, tau_lucy, tau_lucy_1D) in a list as e.g. [.true., .false., .false.]
 !+
 !  OUT: taus:           The array of optical depths to each SPH particle
 !  OUT: tau_lucy:       The array of optical depths to each SPH particle using the Lucy method
 !  OUT: column_density: The array of column densities to each SPH particle
 !+
 !---------------------------------------------------------------------------------
subroutine get_all_integrands_single(npart, primary, Rstar, xyzh, kappa, Rinject, order, tau, tau_lucy, column_density, type)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart,order
 real, intent(in)    :: primary(3), kappa(:), Rstar, Rinject, xyzh(:,:)
 real, intent(out)   :: tau(:), tau_lucy(:), column_density(:)
 logical, dimension(3), intent(in) :: type


 integer  :: i, nrays, nsides
 real     :: ray_dir(3),part_dir(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau, rays_tau_lucy, rays_column_density
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200 ! maximum number of points along the ray where tau is calculated

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays)) ! distance from the central star of the points on the rays
 allocate(rays_tau(ndim, nrays))  ! value of tau at each point along each ray
 allocate(rays_dim(nrays))        ! effective number of points on the ray (< ndim)
 allocate(rays_tau_lucy(ndim, nrays)) ! value of tau_lucy at each point along each ray
 allocate(rays_column_density(ndim, nrays)) ! value of column density at each point along each ray

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the HEALPix ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

 !$omp parallel default(none) &
 !$omp private(ray_dir) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rinject,rays_dist) &
 !$omp shared(rays_tau,rays_tau_lucy,rays_column_density,rays_dim, type)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying a ray (index i-1 because healpix starts counting from index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !calculate the properties along the ray (tau, distance, number of points)
    call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject, &
               rays_tau(:,i), rays_tau_lucy(:,i), rays_column_density(:,i), rays_dist(:,i), rays_dim(i), type)
 enddo
 !$omp enddo
 !$omp end parallel


 !_----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the HEALPix rays
 !-----------------------------------------------

 !$omp parallel default(none) &
 !$omp private(part_dir) &
 !$omp shared(npart,primary,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_tau_lucy) &
 !$omp shared(rays_column_density,tau,tau_lucy,column_density,type,rays_dim)
 !$omp do
 do i = 1,npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       part_dir = xyzh(1:3,i)-primary
       call interpolate_integrands(nsides, part_dir, rays_tau, rays_tau_lucy, rays_column_density, &
                           rays_dist, rays_dim, tau(i), tau_lucy(i), column_density(i), type)
    else
       tau(i) = -99.
    endif
 enddo
 !$omp enddo
 !$omp end parallel

end subroutine get_all_integrands_single

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth at each particle's location, using the uniform outward
 !  ray-tracing scheme for models containing a primary star and a companion
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: npart:           The number of SPH particles
 !  IN: primary:         The xyz coordinates of the primary star
 !  IN: xyzh:            The array containing the particles position+smooting lenght
 !  IN: kappa:           The array containing the opacity of all the SPH particles
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !  IN: companion:       The xyz coordinates of the companion
 !  IN: Rcomp:           The radius of the companion
 !  IN: order:           The healpix order which is used for the uniform ray sampling
 !+
 !  OUT: tau:            The array of optical depths for each SPH particle
 !+
 !--------------------------------------------------------------------------
subroutine get_all_integrands_companion(npart, primary, Rstar, xyzh, kappa, Rinject, companion, Rcomp, &
                                 order, tau, tau_lucy, column_density, type)
 use part, only : isdead_or_accreted
 integer, intent(in) :: npart, order
 real, intent(in)    :: primary(3), companion(3), kappa(:), Rstar, Rinject, xyzh(:,:), Rcomp
 real, intent(out)   :: tau(:), tau_lucy(:), column_density(:)
 logical, dimension(3), intent(in) :: type

 integer  :: i, nrays, nsides
 real     :: normCompanion,theta0,phi,cosphi,sinphi,theta,sep,root
 real     :: ray_dir(3),part_dir(3),uvecCompanion(3)
 real, dimension(:,:), allocatable  :: rays_dist, rays_tau, rays_tau_lucy, rays_column_density
 integer, dimension(:), allocatable :: rays_dim
 integer, parameter :: ndim = 200 ! maximum number of points along the ray where tau is calculated

 nrays = 12*4**order ! The number of rays traced given the healpix order
 nsides = 2**order   ! The healpix nsides given the healpix order

 allocate(rays_dist(ndim, nrays)) ! distance from the central star of the points on the rays
 allocate(rays_tau(ndim, nrays))  ! value of tau at each point along each ray
 allocate(rays_dim(nrays))        ! effective number of points on the ray (< ndim)
 allocate(rays_tau_lucy(ndim, nrays)) ! value of tau_lucy at each point along each ray
 allocate(rays_column_density(ndim, nrays)) ! value of column density at each point along each ray

 uvecCompanion = companion-primary
 normCompanion = norm2(uvecCompanion)
 uvecCompanion = uvecCompanion/normCompanion
 theta0        = asin(Rcomp/normCompanion)
 phi           = atan2(uvecCompanion(2),uvecCompanion(1))
 cosphi        = cos(phi)
 sinphi        = sin(phi)

 !-------------------------------------------
 ! CONSTRUCT the RAYS given the HEALPix ORDER
 ! and determine the optical depth along them
 !-------------------------------------------

 !$omp parallel default(none) &
 !$omp private(ray_dir,theta,root,sep) &
 !$omp shared(nrays,nsides,primary,kappa,xyzh,Rstar,Rinject,Rcomp,rays_dist) &
 !$omp shared(rays_tau,rays_tau_lucy,rays_column_density,rays_dim,type) &
 !$omp shared(uvecCompanion,normCompanion,cosphi,sinphi,theta0)
 !$omp do
 do i = 1, nrays
    !returns ray_dir, the unit vector identifying a ray (index i-1 because healpix starts counting from index 0)
    call pix2vec_nest(nsides, i-1, ray_dir)
    !rotate ray vectors by an angle = phi so the main axis points to the companion (This is because along the
    !main axis (1,0,0) rays are distributed more uniformally
    ray_dir = (/cosphi*ray_dir(1) - sinphi*ray_dir(2),sinphi*ray_dir(1) + cosphi*ray_dir(2), ray_dir(3)/)
    theta   = acos(dot_product(uvecCompanion, ray_dir))
    !the ray intersects the companion: only calculate tau up to the companion
    !if (theta < theta0) then
    !   root  = sqrt(Rcomp**2-normCompanion**2*sin(theta)**2)
    !   sep   = normCompanion*cos(theta)-root
    !   !print *, "I am intersecting"
    !   call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject, &
    !           rays_tau(:,i), rays_tau_lucy(:,i), rays_column_density(:,i), rays_dist(:,i), rays_dim(i), type, sep)
    !else
    call ray_tracer(primary,ray_dir,xyzh,kappa,Rstar,Rinject, &
               rays_tau(:,i), rays_tau_lucy(:,i), rays_column_density(:,i), rays_dist(:,i), rays_dim(i), type)

 enddo
 !$omp enddo
 !$omp end parallel

 !-----------------------------------------------
 ! DETERMINE the optical depth for each particle
 ! using the values available on the HEALPix rays
 !-----------------------------------------------

 !$omp parallel default(none) &
 !$omp private(part_dir) &
 !$omp shared(npart,primary,cosphi,sinphi,nsides,xyzh,ray_dir,rays_dist,rays_tau,rays_dim) &
 !$omp shared(rays_tau_lucy,rays_column_density,tau,tau_lucy,column_density,type)
 !$omp do
 do i = 1, npart
    if (.not.isdead_or_accreted(xyzh(4,i))) then
       !vector joining the source to the particle
       part_dir = xyzh(1:3,i)-primary
       part_dir = (/cosphi*part_dir(1) + sinphi*part_dir(2),-sinphi*part_dir(1) + cosphi*part_dir(2), part_dir(3)/)
       call interpolate_integrands(nsides, part_dir, rays_tau, rays_tau_lucy, rays_column_density, &
                       rays_dist, rays_dim, tau(i), tau_lucy(i), column_density(i), type)
    else
       tau(i) = -99.
    endif
 enddo
 !$omp enddo
 !$omp end parallel
end subroutine get_all_integrands_companion

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth at the SPH particle's location.
 !  Search for the four closest rays to a particle, perform four-point
 !  interpolation of the optical depth from these rays. Weighted by the
 !  inverse square of the perpendicular distance to the rays.
 !
 !  Relies on healpix, for more information: https://healpix.sourceforge.io/
 !+
 !  IN: nsides:          The healpix nsides of the simulation
 !  IN: vec:             The vector from the primary to the particle
 !  IN: rays_tau:        2-dimensional array containing the cumulative optical
 !                       depth along each ray
 !  IN: rays_tau_lucy:   2-dimensional array containing the cumulative optical
 !                       depth along each ray for tau_lucy
 !  IN: rays_column_density: 2-dimensional array containing the cumulative column
 !                       density along each ray
 !  IN: rays_dist:       2-dimensional array containing the distances from the
 !                       primary along each ray
 !  IN: rays_dim:        The vector containing the number of points defined along each ray
 !+
 !  OUT: tau:            The interpolated optical depth at the particle's location
 !+
 !--------------------------------------------------------------------------
subroutine interpolate_integrands(nsides, vec, rays_tau, rays_tau_lucy, rays_column_density, rays_dist, &
                                  rays_dim, tau, tau_lucy, column_density, type)
 integer, intent(in) :: nsides, rays_dim(:)
 real, intent(in)    :: vec(:), rays_dist(:,:), rays_tau(:,:), rays_tau_lucy(:,:), rays_column_density(:,:)
 real, intent(out)   :: tau, tau_lucy, column_density
 logical, dimension(3), intent(in) :: type

 integer :: rayIndex, neighbours(8), nneigh, i, k
 real    :: tautemp, tauLtemp, column_densitytemp, ray(3), vectemp(3), weight, tempdist(8), distRay_sq, vec_norm2
 logical :: mask(8)

 vec_norm2 = norm2(vec)
 !returns rayIndex, the index of the ray vector of the HEALPix cell that points to the particle (direction vec)
 call vec2pix_nest(nsides, vec, rayIndex)
 !returns ray(3), the unit vector identifying the ray with index number rayIndex
 call pix2vec_nest(nsides, rayIndex, ray)
 !compute optical depth along ray rayIndex(+1)
 call get_integrands_on_ray(vec_norm2, rays_tau(:,rayIndex+1), rays_tau_lucy(:,rayIndex+1), &
               rays_column_density(:,rayIndex+1), rays_dist(:,rayIndex+1), rays_dim(rayIndex+1), &
               tautemp, tauLtemp, column_densitytemp, type)
 !determine distance of the particle to the HEALPix ray
 vectemp       = vec - vec_norm2*ray
 distRay_sq    = dot_product(vectemp,vectemp)
 if (distRay_sq > 0.) then
    tau    = tautemp/distRay_sq
    tau_lucy = tauLtemp/distRay_sq
    column_density = column_densitytemp/distRay_sq
    weight = 1./distRay_sq
 else
    ! the particle sits exactly on the ray, no need to interpolate with the neighbours
    tau    = tautemp
    tau_lucy = tauLtemp
    column_density = column_densitytemp
    return
 endif

 !returns the number nneigh and list of vectors (n) neighbouring the ray number rayIndex
 call neighbours_nest(nsides, rayIndex, neighbours, nneigh)
 !for each neighbouring ray calculate its distance to the particle
 do i=1,nneigh
    call pix2vec_nest(nsides, neighbours(i), ray)
    vectemp     = vec - vec_norm2*ray
    tempdist(i) = dot_product(vectemp,vectemp)
 enddo
 neighbours     = neighbours+1
 mask           = .true.
 if (nneigh <8) mask(nneigh+1:8) = .false.
 !take tau contribution from the 3 closest rays
 do i=1,3
    k       = minloc(tempdist,1,mask)
    mask(k) = .false.
    call get_integrands_on_ray(vec_norm2, rays_tau(:,neighbours(k)), rays_tau_lucy(:,neighbours(k)), &
                  rays_column_density(:,neighbours(k)), rays_dist(:,neighbours(k)), rays_dim(neighbours(k)), &
                  tautemp, tauLtemp, column_densitytemp, type)
    tau    = tau + tautemp/tempdist(k)
    tau_lucy = tau_lucy + tauLtemp/tempdist(k)
    column_density = column_density + column_densitytemp/tempdist(k)
    weight = weight + 1./tempdist(k)
 enddo
 tau = tau / weight
 tau_lucy = tau_lucy / weight
 column_density = column_density / weight

! for some reason it is possible for the interpolation to be greater than 2/3, in that case put it to 2/3
 if (tau_lucy > 2./3.) then
    tau_lucy = 2./3.
 endif
end subroutine interpolate_integrands


 !--------------------------------------------------------------------------
 !+
 !  Interpolation of the optical depth for an arbitrary point on the ray,
 !  at a given distance to the starting point of the ray (primary star).
 !+
 !  IN: distance:        The distance from the staring point of the ray to a
 !                       point on the ray
 !  IN: tau_along_ray:   The vector of cumulative optical depths along the ray
 !  IN: dist_along_ray:  The vector of distances from the primary along the ray
 !  IN: len:             The length of tau_along_ray and dist_along_ray
 !+
 !  OUT: tau:            The optical depth to the given distance along the ray
 !+
 !--------------------------------------------------------------------------
subroutine get_integrands_on_ray(distance, tau_along_ray, tauL_along_ray, column_density_along_ray, &
                           dist_along_ray, len, tau, tauL, column_density, type)
 real, intent(in)    :: distance, tau_along_ray(:), tauL_along_ray(:), column_density_along_ray(:), dist_along_ray(:)
 integer, intent(in) :: len
 real, intent(out)   :: tau, tauL, column_density
 logical, dimension(3), intent(in) :: type

 integer :: L, R, m ! left, right and middle index for binary search


 !otherwise a warning is passed
 tau = 0.0
 tauL = 0.0
 column_density = 0.0

 if (distance  <  dist_along_ray(1)) then
    tau = tau_along_ray(1)
    tauL = tauL_along_ray(1)
    column_density = column_density_along_ray(1)
 elseif (distance  >  dist_along_ray(len)) then
    tau = tau_along_ray(len)
    tauL = tauL_along_ray(len)
    column_density = column_density_along_ray(len)

 else
    L = 2
    R = len
    !bysection search for the index of the closest points on the ray to the specified location
    do while (L < R)
       m = (L + R)/2
       if (dist_along_ray(m) > distance) then
          R = m
       else
          L = m + 1
       endif
    enddo

    !linear interpolation of the optical depth at the the point's location
    if (type(1)) then
       tau = tau_along_ray(L-1)+(tau_along_ray(L)-tau_along_ray(L-1))/ &
                     (dist_along_ray(L)-dist_along_ray(L-1))*(distance-dist_along_ray(L-1))
    endif

    if (type(2)) then
       tauL = tauL_along_ray(L-1)+(tauL_along_ray(L)-tauL_along_ray(L-1))/ &
                     (dist_along_ray(L)-dist_along_ray(L-1))*(distance-dist_along_ray(L-1))
    endif
    if (type(3)) then
       column_density = column_density_along_ray(L-1)+(column_density_along_ray(L)-column_density_along_ray(L-1))/ &
                     (dist_along_ray(L)-dist_along_ray(L-1))*(distance-dist_along_ray(L-1))
    endif

 endif
end subroutine get_integrands_on_ray

 !--------------------------------------------------------------------------
 !+
 !  Calculate the optical depth along a given ray
 !+
 !  IN: primary:         The location of the primary star
 !  IN: ray:             The unit vector of the direction in which the
 !                       optical depth will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing lenght
 !  IN: kappa:           The array containing the particles opacity
 !  IN: Rstar:           The radius of the primary star
 !  IN: Rinject:         The particles injection radius
 !  IN: type:            The type of integrands to use (tau, tau_lucy, tau_lucy_1D) in a list as e.g. [.true., .false., .false.]
 !+
 !  OUT: tau_along_ray:  The vector of cumulative optical depth along the ray
 !  OUT: tauL_along_ray: The vector of cumulative optical depth along the ray for tau_lucy
 !  OUT: column_density_along_ray: The vector of cumulative column density along the ray
 !  OUT: dist_along_ray: The vector of distances from the primary along the ray
 !  OUT: len:            The length of tau_along_ray and dist_along_ray
 !+
 !  OPT: maxDistance:    The maximal distance the ray needs to be traced
 !+
 !--------------------------------------------------------------------------
subroutine ray_tracer(primary, ray, xyzh, kappa, Rstar, Rinject, tau_along_ray, tauL_along_ray, column_density_along_ray, &
                  dist_along_ray, len, type, maxDistance)
 use units, only:unit_opacity, unit_density
 !use part,  only:itauL_alloc
 real, intent(in)     :: primary(3), ray(3), Rstar, Rinject, xyzh(:,:), kappa(:)
 real, optional       :: maxDistance
 real, intent(out)    :: dist_along_ray(:), tau_along_ray(:), tauL_along_ray(:), column_density_along_ray(:)
 integer, intent(out) :: len
 real, parameter :: tau_max = 99.0, tauL_max = 2.0/3.0, column_density_max = 1.0e99
 logical, dimension(3), intent(in) :: type


 real    :: dr, next_dr, h, distance
 real    :: dtaudr, previousdtaudr, nextdtaudr, drhodr, previousdrhodr, nextdrhodr, dtauLdr, nextdtauLdr, previousdtauLdr
 integer :: inext, i, L, R, m ! left, right and middle index for binary search

 h = Rinject/100.
 inext=0
 do while (inext==0)
    h = h*2.
    !find the next point along the ray : index inext
    call find_next(primary+Rinject*ray, h, ray, xyzh, kappa, previousdtaudr, previousdrhodr, dr, inext)
 enddo

 if (.not. type(1)) then
    !if tau is not requested, then return empty arrays
    tau_along_ray = 0.
 endif

 if (.not. type(2)) then
    !if tau_lucy is not requested, then return empty arrays
    tauL_along_ray = 0.
 endif

 if (.not. type(3)) then
    !if column density is not requested, then return empty arrays
    column_density_along_ray = 0.
 endif

 i = 1

 tau_along_ray(i)  = 0.
 tauL_along_ray(i) = 0.
 column_density_along_ray(i) = 0.

 distance          = Rinject
 dist_along_ray(i) = distance

 dtauLdr = 0.
 previousdtauLdr = 0.
 nextdtauLdr = 0.

 drhodr = 0.
 previousdrhodr = 0.
 nextdrhodr = 0.

 do while (hasNext(inext,distance,maxDistance))

    distance = distance+dr
    call find_next(primary + distance*ray, xyzh(4,inext), ray, xyzh, kappa, nextdtaudr, nextdrhodr, next_dr, inext)
    i = i + 1

    if (type(1)) then
       dtaudr            = (nextdtaudr+previousdtaudr)/2.
       previousdtaudr    = nextdtaudr
       tau_along_ray(i) = tau_along_ray(i-1) + dr*dtaudr/unit_opacity
    endif

    if (type(2)) then
       nextdtauLdr = nextdtaudr * (Rstar/distance)**2
       dtauLdr = (nextdtauLdr + previousdtauLdr)/2.
       previousdtauLdr = nextdtauLdr
       tauL_along_ray(i) = tauL_along_ray(i-1) + dr*dtauLdr/unit_opacity
    endif

    if (type(3)) then
       nextdrhodr = drhodr
       drhodr = (nextdrhodr + previousdrhodr)/2.
       previousdrhodr = nextdrhodr
       column_density_along_ray(i) = column_density_along_ray(i-1) + dr*drhodr/unit_density
    endif

    dist_along_ray(i) = distance
    dr                = next_dr

 enddo

 if (type(1) .and. present(maxDistance)) then
    i = i + 1
    tau_along_ray(i)  = tau_max
    dist_along_ray(i) = maxDistance
 endif

 if (type(3) .and. present(maxDistance)) then
    i = i + 1
    column_density_along_ray(i) = column_density_max
    dist_along_ray(i) = maxDistance
 endif
 len = i

 if (type(2)) then
    !reverse integration start from zero inward
    tauL_along_ray(1:len) = tauL_along_ray(len) - tauL_along_ray(1:len)
    !find the first point where tau_lucy < 2/3
    if (tauL_along_ray(1) > 2./3.) then
       !print *, tauL_along_ray
       L = 1
       R = len
       !bysection search for the index of the closest point to tau = 2/3
       do while (L < R)
          m = (L + R)/2
          if (tauL_along_ray(m) < 2./3.) then
             R = m
          else
             L = m + 1
          endif
       enddo
       tauL_along_ray(1:L-1) = 2./3.
       !print *, tauL_along_ray
       !The photosphere is located between ray grid point L and L+1, may be useful information!
    endif
 endif

end subroutine ray_tracer

logical function hasNext(inext, distance, maxDistance)
 integer, intent(in) :: inext
 real, intent(in)    :: distance
 real, optional      :: maxDistance
 if (present(maxDistance)) then
    hasNext = inext /= 0 .and. distance < maxDistance
 else
    hasNext = inext /= 0
 endif
end function hasNext

 !--------------------------------------------------------------------------
 !+
 !  First finds the local optical depth derivative at the starting point, then finds the next
 !                       point on a ray and the distance to this point
 !+
 !  IN: inpoint:         The coordinate of the initial point projected on the ray
 !                       for which the opacity and the next point will be calculated
 !  IN: h:               The smoothing length at the initial point
 !  IN: ray:             The unit vector of the direction in which the next
 !                       point will be calculated
 !  IN: xyzh:            The array containing the particles position+smoothing length
 !  IN: kappa:           The array containing the particles opacity
 !  IN: inext:           The index of the initial point
 !                       (this point will not be considered as possible next point)
 !+
 !  OUT: dtaudr:         The radial optical depth derivative at the given location (inpoint)
 !  OUT: distance:       The distance to the next point
 !  OUT: inext:          The index of the next point on the ray
 !+
 !--------------------------------------------------------------------------
subroutine find_next(inpoint, h, ray, xyzh, kappa, dtaudr, drhodr, distance, inext)
 use linklist, only:getneigh_pos,ifirstincell,listneigh
 use kernel,   only:radkern,cnormk,wkern
 use part,     only:hfact,rhoh,massoftype,igas
 use dim,      only:maxpsph
 real,    intent(in)    :: xyzh(:,:), kappa(:), inpoint(:), ray(:), h
 integer, intent(inout) :: inext
 real,    intent(out)   :: distance, dtaudr, drhodr

 integer, parameter :: nmaxcache = 0
 real  :: xyzcache(0,nmaxcache)

 integer  :: nneigh, i, prev,j
 real     :: dmin, vec(3), dr, raydistance, q, norm_sq

 prev     = inext
 inext    = 0
 distance = 0.

 !for a given point (inpoint), returns the list of neighbouring particles (listneigh) within a radius h*radkern
 call getneigh_pos(inpoint,0.,h*radkern,3,listneigh,nneigh,xyzh,xyzcache,nmaxcache,ifirstincell)

 dtaudr = 0.
 drhodr = 0.
 dmin = huge(0.)
 !loop over all neighbours
 do i=1,nneigh
    j = listneigh(i)
    if (j > maxpsph) cycle
    vec     = xyzh(1:3,j) - inpoint
    norm_sq = dot_product(vec,vec)
    q       = sqrt(norm_sq)/xyzh(4,j)
    !add optical depth contribution from each particle
    dtaudr = dtaudr+wkern(q*q,q)*kappa(j)*rhoh(xyzh(4,j), massoftype(igas))
    drhodr = drhodr+wkern(q*q,q)*rhoh(xyzh(4,j), massoftype(igas))

    ! find the next particle : among the neighbours find the particle located the closest to the ray
    if (j  /=  prev) then
       dr = dot_product(vec,ray) !projected distance along the ray
       if (dr>0.) then
          !distance perpendicular to the ray direction
          raydistance = norm_sq - dr**2
          if (raydistance < dmin) then
             dmin     = raydistance
             inext    = j
             distance = dr
          endif
       endif
    endif
 enddo
 dtaudr = dtaudr*cnormk/hfact**3
 drhodr = drhodr*cnormk/hfact**3
end subroutine find_next

end module raytracer
