!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2026 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module readwrite_starevol
!
! Utility routines for reading stellar profiles from the STAREVOL code
!
! :References: Siess et al. (2000)
!
! :Owner: Lionel Siess
!
! :Runtime parameters: None
!
! :Dependencies: datafiles, fileutils, physcon, table_utils, units
!
 implicit none

 public :: read_starevol,write_starevol

 private

contains

!-----------------------------------------------------------------------
!+
!  Read quantities from STAREVOL structure file
!+
!-----------------------------------------------------------------------
subroutine read_starevol(filepath,rho,r,pres,m,ene,temp,X_in,Z_in,Xfrac,Yfrac,mu,Mstar,ierr,cgsunits)
 use physcon,   only:solarm,solarr
 use fileutils, only:get_nlines,get_ncolumns,string_delete,lcase,read_column_labels
 use table_utils, only:flip_array
 use datafiles, only:find_phantom_datafile
 use units,     only:udist,umass,unit_density,unit_pressure,unit_ergg
 character(len=*),  intent(in)  :: filepath
 integer,           intent(out) :: ierr
 real,              intent(in)  :: X_in,Z_in
 real, allocatable, intent(out) :: rho(:),r(:),pres(:),m(:),ene(:),temp(:),Xfrac(:),Yfrac(:),mu(:)
 real,              intent(out) :: Mstar
 logical,           intent(in), optional :: cgsunits
 integer                                 :: lines,i,ncols,nheaderlines,nlabels
 integer                                 :: iu
 character(len=120)                      :: fullfilepath
 character(len=24), allocatable          :: header(:)
 logical                                 :: iexist,usecgs,isstarevolfile,got_column
 real, allocatable :: dat(:,:)

 usecgs = .false.
 if (present(cgsunits)) usecgs = cgsunits
 !
 !--Get path name
 !
 ierr = 0
 fullfilepath = find_phantom_datafile(filepath,'star_data_files')
 inquire(file=trim(fullfilepath),exist=iexist)
 if (.not.iexist) then
    ierr = 1
    return
 endif
 lines = get_nlines(fullfilepath) ! total number of lines in file

 open(newunit=iu,file=fullfilepath,status='old',iostat=ierr)
 if (ierr /= 0) then
    print "(a,/)",' ERROR opening file '//trim(fullfilepath)
    return
 endif

 call get_ncolumns(iu,ncols,nheaderlines)
 lines = lines - nheaderlines
 if (lines <= 0) then ! file not found
    ierr = 1
    return
 endif

 ! extract column labels from the file header
 allocate(header(ncols),dat(lines,ncols))
 call read_column_labels(iu,nheaderlines,ncols,nlabels,header)
 isstarevolfile = .false.
 do i = 1,ncols
    if (trim(lcase(header(i))) == '#m/msun' .or. &
        trim(lcase(header(i))) == 'm/msun') isstarevolfile = .true.
 enddo
 if (nlabels /= ncols) print*,'WARNING: got ',nlabels,' labels for ',ncols,' columns'

 allocate(m(lines))
 m = -1.
 allocate(r,pres,rho,ene,temp,Xfrac,Yfrac,mu,source=m)

 ! read file forwards, from centre to surface
 do i = 1,lines
    read(iu,*,iostat=ierr) dat(i,1:ncols)
 enddo
 if (ierr /= 0) then
    print "(a,/)",' ERROR reading data from STAREVOL file (new statement)'
    return
 endif

 ! Set mass fractions to fixed inputs if not in file
 Xfrac = X_in
 Yfrac = 1. - X_in - Z_in
 mu = 0.
 do i = 2,ncols
    if (header(i)(1:1) == '#' .and. .not. trim(lcase(header(i)))=='#shell') then
       print '("Detected wrong header entry : ",a," in file ",a)',trim(lcase(header(i))),trim(fullfilepath)
       ierr = 2
       return
    endif
    got_column = .true.
    select case(trim(lcase(header(i))))
    case('Mr')
       m = dat(1:lines,i)
    case('rho')
       rho = dat(1:lines,i)
    case('e')
       ene = dat(1:lines,i)
    case('r')
       r = dat(1:lines,i)
       if (maxval(r) < 1e-5*solarr) r = r * solarr
    case('P')
       pres = dat(1:lines,i)
    case('T')
       temp = dat(1:lines,i)
    case('H1')
       Xfrac = dat(1:lines,i)
    case('HE4')
       Yfrac = dat(1:lines,i)
    case('mu')
       mu = dat(1:lines,i)
    case default
       got_column = .false.
    end select
    if (got_column) print "(1x,i0,': ',a)",i,trim(header(i))
 enddo
 print "(a)"

 close(iu)

 if (min(minval(pres),minval(rho),minval(r),minval(m)) < 0d0) ierr = 4

 if (ierr /= 0) then
    print "(a,/)",' ERROR reading STAREVOL file [missing required columns]'
    return
 endif

 if (.not. usecgs) then
    m = m / umass
    r = r / udist
    pres = pres / unit_pressure
    rho = rho / unit_density
    ene = ene / unit_ergg
 endif

 if (r(1) > r(lines)) then
    call flip_array(r)
    call flip_array(m)
    call flip_array(rho)
    call flip_array(pres)
    call flip_array(ene)
    call flip_array(temp)
    call flip_array(Xfrac)
    call flip_array(Yfrac)
    call flip_array(mu)
 endif

 Mstar = m(lines)

end subroutine read_starevol

!----------------------------------------------------------------
!+
!  Write stellar profile in format readable by read_starevol;
!  used in star setup to write softened stellar profile.
!+
!----------------------------------------------------------------
subroutine write_starevol(outputpath,m,pres,temp,r,rho,ene,Xfrac,Yfrac,csound,mu)
 use physcon, only:solarm
 real,               intent(in) :: m(:),rho(:),pres(:),r(:),ene(:),temp(:)
 character(len=120), intent(in) :: outputpath
 real,               intent(in), optional :: Xfrac(:),Yfrac(:),csound(:),mu(:)
 character(len=200)              :: headers(100)
 integer                         :: i,ncols,noptionalcols,j,iu
 real, allocatable               :: optionalcols(:,:)
 character(len=*), parameter     :: fmtstring = "(5(es24.16e3,2x),es24.16e3)"

 ncols = 6
 headers(1:ncols) = (/'         #M/Msun','       radius_cm','         density', &
                      ' internal_energy','     temperature','        pressure'/)

 ! Add optional columns
 noptionalcols = 0
 allocate(optionalcols(size(r),10))
 if (present(Xfrac)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '        hydrogen'
    optionalcols(:,noptionalcols) = Xfrac
 endif
 if (present(Yfrac)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '          helium'
    optionalcols(:,noptionalcols) = Yfrac
 endif
 if (present(mu)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = 'molecular_weight'
    optionalcols(:,noptionalcols) = mu
 endif
 if (present(csound)) then
    noptionalcols = noptionalcols + 1
    headers(noptionalcols+ncols) = '          csound'
    optionalcols(:,noptionalcols) = csound
 endif

 open(newunit=iu,file=outputpath,status='replace')
 write(iu,'(a)') '# STAREVOL stellar profile written by Phantom'
 write(iu,'(a)') '# Header line'
 write(iu,'(a)') '# Header line'
 write(iu,'(a)') '# Header line'
 do i = 1,noptionalcols+ncols-1
    write(iu,'(a24,2x)',advance="no") trim(adjustl(headers(i)))
 enddo
 write(iu,'(a24)') trim(adjustl(headers(noptionalcols+ncols)))

 do i=1,size(r)
    if (noptionalcols <= 0) then
       write(iu,fmtstring) m(i)/solarm,r(i),rho(i),ene(i),temp(i),pres(i)
    else
       write(iu,fmtstring,advance="no") m(i)/solarm,r(i),rho(i),ene(i),temp(i),pres(i)
       do j=1,noptionalcols
          if (j==noptionalcols) then
             write(iu,'(2x,es24.16e3)') optionalcols(i,j)
          else
             write(iu,'(2x,es24.16e3)',advance="no") optionalcols(i,j)
          endif
       enddo
    endif
 enddo
 close(iu)

end subroutine write_starevol

end module readwrite_starevol
