!@sum sumfiles is a generic summation program for modelE acc-files.
!@+   See conventions.txt for documentation.
!@auth M. Kelley
      program sumfiles
      implicit none
      include 'netcdf.inc'
      integer :: status,ofid,ivarid,varid
      integer, dimension(:), allocatable :: fids
      character(len=4096) :: ifile,ofile
      integer :: n,nfiles,nlast,iargc,nvars
      integer :: itbeg,itend,itnow,itime0,itime,nday,iyear1,accsize
      integer :: Jyear0,Jmon0,Jday0,Jdate0,Jhour0,
     &           Jyear,Jmon,Jday,Jdate,Jhour
      character(len=4) :: amon0,amon
      real*4 :: days
      character(len=12) :: acc_period
      character(len=30) :: runid,reduction,vname
      character(len=100) :: fromto
      character(len=132) :: xlabel
      integer, dimension(12) :: monacc,monacc1
      real*8, dimension(:), allocatable :: acc,acc_part
      integer :: chunksize
      real*8 :: bynfiles
c
c get the number of input files
c
      nfiles = iargc()
      bynfiles = 1d0/real(nfiles,kind=8)

      if(nfiles.le.1) then
        write(6,*)
     &       'usage: sumfiles files_to_be_summed'
        write(6,*)
     &       '(works on modelE acc files, not on pdE outputs)'
        call exit(-1)
      endif
      allocate(fids(nfiles))

      chunksize = 1024*1024*32

c
c open each input file and find the min/max itime
c
      itbeg=+huge(itbeg)
      itend=-huge(itend)
      monacc(:) = 0
      do n=1,nfiles
        call getarg(n,ifile)
        status = nf_open(trim(ifile),nf_nowrite,fids(n))
c        status = nf__open(trim(ifile),nf_nowrite,chunksize,fids(n))
        if(status.ne.nf_noerr) then
          write(6,*) 'nonexistent/non-netcdf input file ',trim(ifile)
          call exit(-2)
        endif
        call get_var_int(fids(n),'itime0',itnow)
        itbeg = min(itnow,itbeg)
        call get_var_int(fids(n),'itime',itnow)
        if(itnow.gt.itend) then
          itend = itnow
          nlast = n
        endif
        call get_var_int(fids(n),'monacc',monacc1) 
        monacc(:) = monacc(:) + monacc1(:)
      enddo

c
c determine an appropriate name for the averaging period
c
      call get_var_int(fids(nlast),'iyear1',iyear1)
      call get_var_int(fids(nlast),'nday',nday)
      itime  = itend
      itime0 = itbeg
      call getdte(Itime0,Nday,Iyear1,Jyear0,Jmon0,Jday0,Jdate0,Jhour0
     *     ,amon0)
      call getdte(Itime-1,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour
     *     ,amon)
      call aperiod(monacc,jyear0,jyear,acc_period)

c
c copy the structure of the latest input file to the output file
c
      xlabel=''
      status = nf_get_att_text(fids(nlast),nf_global,'xlabel',xlabel)
      runid = xlabel(1:index(xlabel,' ')-1)
      ofile = trim(acc_period)//'.acc'//trim(runid)//'.nc'
      status = nf_create(trim(ofile),nf_clobber,ofid)
      call copy_file_structure(fids(nlast),ofid)

c
c copy the contents of the latest input file to the output file
c and write the appropriate itime0,monacc,fromto to the output file
c
c      call copy_selected_vars(fids(nlast),ofid)
      call copy_shared_vars(fids(nlast),ofid)
      call put_var_int(ofid,'itime0',itime0)
      call put_var_int(ofid,'monacc',monacc)
      days=(itime-itime0)/float(nday)
      write(fromto,902) jyear0,amon0,jdate0,jhour0,
     &     jyear,amon,jdate,jhour,itime,days
      status = nf_put_att_text(ofid,nf_global,'fromto' 
     &     ,len_trim(fromto),fromto)

c
c loop over the fields to be reduced
c
      status = nf_inq_nvars(ofid,nvars)
      do varid=1,nvars
        reduction=''
        status=nf_get_att_text(ofid,varid,'reduction',reduction)
        if(status.ne.nf_noerr) cycle
        status = nf_inq_varname(ofid,varid,vname)
        call get_varsize(ofid,vname,accsize)
        allocate(acc(accsize),acc_part(accsize))
        select case(trim(reduction))
        case ('min')
          acc = +1d30
        case ('max')
          acc = -1d30
        case default
          acc = 0.
        end select
        do n=1,nfiles
          status = nf_inq_varid(fids(n),vname,ivarid)
          status = nf_get_var_double(fids(n),ivarid,acc_part)
          select case(trim(reduction))
          case ('min')
            acc = min(acc,acc_part)
          case ('max')
            acc = max(acc,acc_part)
          case ('avg')
            acc = acc + acc_part*bynfiles
          case default
            acc = acc + acc_part
          end select
        enddo
        status = nf_put_var_double(ofid,varid,acc)
        deallocate(acc,acc_part)
      enddo

c
c close input and output files
c
      status = nf_close(ofid)
      do n=1,nfiles
        status = nf_close(fids(n))
      enddo

      deallocate(fids)

  902 FORMAT ('From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      end program sumfiles

      subroutine getdte(It,Nday,Iyr0,Jyr,Jmn,Jd,Jdate,Jhour,amn)
!@sum  getdte gets julian calendar info from internal timing info
!@auth Gavin Schmidt
      IMPLICIT NONE
      real*8, parameter :: hrday=24.
      integer, parameter :: jmpery=12,jdpery=365
      integer, parameter, dimension(0:jmpery) :: JDendOfM = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      character(len=4), dimension(0:jmpery), parameter :: amonth = (/
     &  'IC  ',
     *  'JAN ','FEB ','MAR ','APR ','MAY ','JUNE',
     *  'JULY','AUG ','SEP ','OCT ','NOV ','DEC '/)
      INTEGER, INTENT(IN) :: It,Nday,Iyr0
      INTEGER, INTENT(OUT) :: Jyr,Jmn,Jd,Jdate,Jhour
      CHARACTER*4, INTENT(OUT) :: amn

      Jyr=Iyr0+It/(Nday*JDperY)
      Jd=1+It/Nday-(Jyr-Iyr0)*JDperY
      Jmn=1
      do while (Jd.GT.JDendOfM(Jmn))
        Jmn=Jmn+1
      end do
      Jdate=Jd-JDendOfM(Jmn-1)
      Jhour=nint(mod(It*hrday/Nday,hrday))
      amn=amonth(Jmn)

      return
      end subroutine getdte

      subroutine aperiod(monacc,yr_start,yr_end,acc_period)
      implicit none
      integer :: monacc(12)
      integer :: yr_start,yr_end
      character(len=12) :: acc_period
      character(len=3), dimension(12), parameter :: amonth = (/
     &  'JAN','FEB','MAR','APR','MAY','JUN',
     &  'JUL','AUG','SEP','OCT','NOV','DEC' /)
      character(len=12) :: mostr
      integer :: mons(12)
      integer :: m,nmo,yr1,ninc,ndec
      logical :: incyr1
     
      if(minval(monacc,mask=monacc>0).ne.maxval(monacc,mask=monacc>0))
     &     stop 'unequal numbers of months'

      mostr=''
      nmo = 0
      do m=1,12
        if(monacc(m).eq.0) cycle
        nmo = nmo + 1
        if(nmo.eq.1) then
          mostr(1:3)=amonth(m)
        else
          mostr(nmo:nmo)=amonth(m)(1:1)
        endif
        mons(nmo) = m
      enddo

      incyr1 = .false.
      if(nmo.eq.12) then
        mostr='ANN'
      elseif(mostr.eq.'JFD' .and. yr_end.gt.yr_start) then
        incyr1 = .true.
        mostr = 'DJF'
      elseif(nmo.eq.2) then
        mostr(3:3) = mostr(2:2)
        mostr(2:2) = '+'
      elseif(nmo.gt.3) then
        if(monacc(1).eq.0 .or. yr_start.eq.yr_end) then
          mostr(3:3) = mostr(nmo:nmo)
          mostr(2:2) = '-'
        else
          incyr1 = .true.
          do m=12,1,-1
            mostr(1:1) = amonth(m)(1:1)
            if(monacc(m-1).eq.0) exit
          enddo
          do m=1,12
            mostr(3:3) = amonth(m)(1:1)
            if(monacc(m+1).eq.0) exit
          enddo
          mostr(2:2) = '-'
        endif
      endif
      yr1 = yr_start
      if(incyr1) yr1 = yr1+1

      acc_period=''
      acc_period(1:3) = mostr(1:3)
      write(acc_period(4:7),'(i4.4)') yr1
      if(yr_end.gt.yr1) write(acc_period(8:12),'(a1,i4.4)') '-',yr_end

c check for gaps
      ninc = count(monacc(2:12).gt.monacc(1:11))
      ndec = count(monacc(2:12).lt.monacc(1:11))
      if(ninc+ndec.gt.2) then
        write(6,*) 'gap'
      endif

      return
      end subroutine aperiod
