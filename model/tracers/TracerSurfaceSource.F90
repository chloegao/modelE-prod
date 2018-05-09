#include "rundeck_opts.h"

module TracerSurfaceSource_mod
  use TracerSource_mod
  use timestream_mod, only : timestream
  implicit none
  private

!@var EMstream interface for reading and time-interpolating emissions file. See timestream_mod.
!@var sourceName holds source name from metadata, used to set up diagnostics.
!@var sourceLname holds long source name, used to set up diagnostics.
!@var tracerName tracer name associated with emissions file (often matches trname)
!@var skipReason when gt 0 tags sources to be skipped in trname_## file reading
!@var firstTrip indicates source initialization is needed
!@param itsMegan skip reason is online MEGAN vegetation source
!@param itsCH4MGOL skip reason is CH4 ocean, lake, misc ground source
!@param itsOcean skip reason is online ocean source

  public :: TracerSurfaceSource
  public :: initSurfaceSource
  public :: readSurfaceSource
  public :: itsMegan
  public :: itsOcean
  public :: itsCH4MGOL

  type, extends(TracerSource) :: TracerSurfaceSource
    character(len=30) :: sourceName
    character(len=30) :: sourceLname
    type (timestream) :: EMstream
    character(len=10) :: tracerName
    integer :: skipReason = 0
    logical :: firstTrip = .true.
  end type TracerSurfaceSource

  integer, parameter :: itsMegan=1
  integer, parameter :: itsCH4MGOL=2
  integer, parameter :: itsOcean=3

contains

  subroutine initSurfaceSource(this, tracerName, fileName, sectorNames)
    use SystemTools, only : stLinkStatus,stFileList
    use TracerSource_mod, only: N_MAX_SECT
    use Dictionary_mod, only : sync_param
    USE FILEMANAGER, only: openunit,closeunit
    use pario, only : par_open,par_close,read_attr
    USE DOMAIN_DECOMP_ATM, only: GRID
    use TimeConstants_mod, only: HOURS_PER_DAY
    use SpecialIO_mod, only: write_parallel,read_parallel
    type (TracerSurfaceSource), intent(inout) :: this
    character(len=*), intent(in) :: tracerName
    character(len=*), intent(in) :: fileName
    character(len=300) :: out_line
    character*10, intent(in):: sectorNames(:)
    logical :: diurnalFileExists = .false.

    integer :: nsect, nn, i, j, iu, fid
    integer :: linkstatus, nfiles, ifile, ios, jyr
    integer, parameter :: max_fname_len=128
    character(len=max_fname_len), allocatable :: flist(:)
    character(len=max_fname_len) :: thisline
    character(len=max_fname_len+8) :: fileToRead
    character(len=4) :: c4
    character*32 :: pname
    character*35 :: fname
    character*124 :: tr_sectors_are
    integer :: numTrSectors
    character(len=80) :: name ! sector
    real*8 :: sumDiurnal
    real*8, parameter :: diurnalSumTolerance=1.d-4
    character*80 :: targetVariable

    fileToRead=fileName ! default (e.g. if file is not a directory, or
                        ! the directory search doesn't find a good file)
    call stLinkStatus(trim(fileName), linkstatus)
    if(linkstatus==2) then  ! this is a directory. todo: no hard-coded retcodes
      ! The file is a directory. Do something similar to subroutine check_metadata
      ! in timestream_mod to determine any useable YYYY.nc file in the directory.
      allocate(flist(1000)) ! 1000 files maximum
      call stFileList(trim(fileName),flist,nfiles)
      do ifile=1,nfiles
        thisline = adjustl(flist(ifile))
        if(len_trim(thisline).ne.7) cycle
        if(thisline(5:7).ne.'.nc') cycle
        c4 = thisline(1:4)
        read(c4,*,iostat=ios) jyr
        if(ios.ne.0) cycle
        if(jyr.lt.0) cycle
        ! acceptable file. define it and exit:
        fileToRead=trim(fileName)//'/'//trim(thisline)
        exit
      end do
      deallocate(flist)
    end if ! directory search

    ! continue reading netCDF file:
    this%tracerName = tracerName
    this%sourceName = 'notfound'
    fid = par_open(grid,trim(fileToRead),'read')
    ! First try to read the variable attribute to get source name:
    call read_attr(grid,fid,this%tracerName,'source',i,this%sourceName)
    ! If that fails, look for the source attribute of the variable
    ! that varname attribute points to (like init_stream would):
    if(trim(this%sourceName).eq.'notfound') then
      targetVariable=this%tracerName
      call read_attr(grid,fid,'global',trim(this%tracerName)//'name',&
      & i,targetVariable)
      call read_attr(grid,fid,trim(targetVariable),'source',&
      & i,this%sourceName)
    endif
    ! If that fails, look for a global source attribute:
    if(trim(this%sourceName).eq.'notfound') then
      call read_attr(grid,fid,'global','source',i,this%sourceName)
    endif
    ! If even that fails, stop the model:
    call par_close(grid,fid)
    if(trim(this%sourceName).eq.'notfound') then
      call stop_model('source name not found in file '//trim(fileName),255)
    endif

    ! append ' source' to the long name, and '_src' to the short name
    this%sourceLname = trim(this%sourceName)//' source'
    this%sourceName = trim(this%sourceName)//'_src'

    ! -- begin sector stuff --
    tr_sectors_are = ' '
    pname=trim(trim(fileName)//'_sect')
    call sync_param(pname,tr_sectors_are)
    numTrSectors = 0

    i=1
    do while(i < len(tr_sectors_are))
      j=index(tr_sectors_are(i:len(tr_sectors_are))," ")
      if (j > 1) then
        numTrSectors = numTrSectors + 1
        i=i+j
      else
        i=i+1
      end if
    enddo
    if(numTrSectors > n_max_sect)  &
         &     call stop_model("num_tr_sectors problem",255)
    this%num_tr_sectors = numTrSectors

    if(numTrSectors > 0) then
      read(tr_sectors_are,*) this%tr_sect_name(1:numTrSectors)

      do nsect=1, numTrSectors
        name = trim(this%tr_sect_name(nsect))
        this%tr_sect_index(nsect) = 0
        loop_nn: do nn=1, size(sectorNames)
          if(trim(name) == trim(sectorNames(nn))) then
            this%tr_sect_index(nsect) = nn
            exit loop_nn
          endif
        enddo loop_nn
      enddo
    endif

    ! -- begin diurnal stuff -- 
    fname=trim('diurnal_'//trim(fileName))
    ! governed by file existance:
    inquire(file=trim(fname), exist=diurnalFileExists)
    if(diurnalFileExists)then
       this%applyDiurnalCycle=.true.
       write(out_line,*)'Applying diurnal cycle to file '//fileName
       call write_parallel(trim(out_line))
       call openunit(fname,iu,.false.,.true.)
       call read_parallel(this%diurnalCycle,iu)
       ! check that the diurnal cycle's sum is close to the number of hours
       ! in a day (meaning it's hourly average would be a factor of 1.):
       sumDiurnal=SUM(this%diurnalCycle)
       if(  sumDiurnal > HOURS_PER_DAY + diurnalSumTolerance  &
     & .or. sumDiurnal < HOURS_PER_DAY - diurnalSumTolerance) then
         write(out_line,*) &
     &   trim(fname),' sum is ',sumDiurnal,' not',HOURS_PER_DAY
         call write_parallel(trim(out_line))
         call stop_model('Problem with emissions diurnal cycle.',255)
       end if
       call closeunit(iu)
    end if

  end subroutine initSurfaceSource


  subroutine readSurfaceSource(this, fname, sfc_src, xyear, xday, isChemTracer)
    USE DOMAIN_DECOMP_ATM, only: GRID,  readt_parallel, write_parallel
    use Domain_decomp_atm, only: getDomainBounds
    use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
    use timestream_mod, only : init_stream,read_stream
    use dictionary_mod, only : get_param
    type (TracerSurfaceSource), intent(inout) :: this
    character(*), intent(in) :: fname
    real*8, intent(inout) :: sfc_src(grid%i_strt_halo:,grid%j_strt_halo:)
    integer, intent(in) :: xyear, xday
    logical, intent(in) :: isChemTracer

    integer :: iu,k,ipos,kx,iposDay,kstep=10
    character(len=300) :: out_line
    real*8 :: alpha
    real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO, &
         &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: &
         & sfc_a,sfc_b

    INTEGER :: J_1, J_0, I_0, I_1
    integer :: cyclic_yr,master_yr,nc_emis_use_ppm_interp

    if(this%firstTrip) then
      this%firstTrip = .false.
      call get_param('master_yr',master_yr)
      if (isChemTracer) then
        call get_param('o3_yr',cyclic_yr,default=master_yr)
      else
        call get_param('aer_int_yr',cyclic_yr,default=master_yr)
      end if
      cyclic_yr=ABS(cyclic_yr)
      call get_param('nc_emis_use_ppm_interp',nc_emis_use_ppm_interp,&
        & default=1)
      if (nc_emis_use_ppm_interp==1) then
        call init_stream(grid,this%EMstream,trim(fname), &
           trim(this%tracername),0d0,1d30,'ppm',xyear,xday, &
           cyclic = (cyclic_yr > 0) )
      else
        call init_stream(grid,this%EMstream,trim(fname), &
           trim(this%tracername),0d0,1d30,'linm2m',xyear,xday, &
           cyclic = (cyclic_yr > 0) )
      endif
    endif
    call read_stream(grid,this%EMstream,xyear,xday,sfc_src)

  end subroutine readSurfaceSource

end module TracerSurfaceSource_mod
