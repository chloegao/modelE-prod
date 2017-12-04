!@sum sumfiles is a generic summation program for modelE acc-files.
!@+   See conventions.txt for documentation.
!@auth M. Kelley
      program sumfiles
      implicit none
      include 'netcdf.inc'
      integer :: status,ofid,ivarid,varid
      integer, dimension(:), allocatable :: fids
      character(len=4096) :: ifile,ofile
      integer :: n,nfiles,nfirst,nlast,iargc,nvars,n1Dif,n1To
      integer :: itbeg,itend,itnow,nday,iyear1,accsize
      character(len=6) :: ayear0,ayear
      character(len=3) :: amon0,amon
      real*4 :: days, difs, dif
      character(len=16) :: acc_period
      character(len=30) :: runid,reduction,vname
      character(len=100) :: fromto,fromto0
      character(len=132) :: xlabel
      integer, dimension(12) :: monacc,monacc1
      real*8, dimension(:), allocatable :: acc,acc_part
!csz  integer :: chunksize
      real*8 :: bynfiles
      logical :: do_dif = .true.
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
      end if
      allocate(fids(nfiles))

!csz  chunksize = 1024*1024*32

c
c open each input file and find the min/max itime
c
      itbeg=+huge(itbeg)
      itend=-huge(itend)
      monacc(:) = 0 ; difs = 0 ; dif = 0
      do n=1,nfiles
        call getarg(n,ifile)
        status = nf_open(trim(ifile),nf_nowrite,fids(n))
!csz     status = nf__open(trim(ifile),nf_nowrite,chunksize,fids(n))
        if(status.ne.nf_noerr) then
          write(6,*) 'nonexistent/non-netcdf input file ',trim(ifile)
          call exit(-2)
        end if
        call get_var_int(fids(n),'itime0',itnow)
        if(itnow.lt.itbeg) then
          itbeg = itnow
          nfirst = n
        end if
        call get_var_int(fids(n),'itime',itnow)
        if(itnow.gt.itend) then
          itend = itnow
          nlast = n
        end if
        call get_var_int(fids(n),'monacc',monacc1)
        monacc(:) = monacc(:) + monacc1(:)
        status = nf_get_att_text(fids(n),nf_global,'fromto',fromto)
        n1Dif = index(fromto,'Dif:')+4
        if (fromto(n1Dif:n1Dif+6) == '       ') do_dif = .false.
        if (do_dif) read(fromto(n1Dif:n1Dif+7),*) dif
        difs = difs + dif
      enddo

c
c determine an appropriate name for the averaging period
c
      status = nf_get_att_text(fids(nfirst),nf_global,'fromto',fromto0)
      status = nf_get_att_text(fids(nlast),nf_global,'fromto',fromto)
      read(fromto0(6:16),'(a6,2x,a3)') ayear0,amon0
      n1To = index(fromto,'To:')+3
      read(fromto (n1To:n1To+10),'(a6,2x,a3)') ayear,amon
      call aperiod(monacc,ayear0,ayear,acc_period,amon0,amon)

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
      call put_var_int(ofid,'itime0',itbeg)
      call put_var_int(ofid,'monacc',monacc)
      fromto(1:index(fromto0,'To:')) = fromto0(1:index(fromto0,'To:'))
      n1Dif = index(fromto,'Dif:')+4
      fromto(n1Dif:n1Dif+6) = '       '
      if (do_dif) then
        if (difs < 10000.) then
          write(fromto(n1Dif:n1Dif+6),'(f7.2)') difs
        else if (difs < 1000000.) then
          write(fromto(n1Dif:n1Dif+6),'(f7.0)') difs
        end if
      end if
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

      end program sumfiles


      subroutine aperiod(monacc,ayr0,ayr1,acc_period,amon0,amon)
! Find appropriate name for the accumulation period e.g. MonYear1-Year2
! For NonEarth calendars we assume the months are called A..,B..,C..,..
!      Notes about the output file names:
! restrictions: month_per_year is currently restricted to 12,
!               the months are called JAN FEB ... (all capitals)
! ambivalence:  different periods may have the same name X-Z...
!               e.g. A-M may be Aug-Mar or Aug-May
!               the name of any 12-month period is ANN
      implicit none
      integer :: monacc(12),mon0,mon1,yr0,yr1
      character(len=3) :: amon0,amon
      character(len=6) :: ayr0,ayr1
      character(len=16) :: acc_period
      character(len=26), parameter :: alphabet =
     *  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=3), dimension(12) :: emonth = (/
     &  'JAN','FEB','MAR','APR','MAY','JUN',
     &  'JUL','AUG','SEP','OCT','NOV','DEC' /)
      character(len=12) :: mostr
      integer :: m,mm,nmo,ninc,ndec,n1Dif,n1To
      logical :: incyr0, earth = .true.

      if(minval(monacc,mask=monacc>0).ne.maxval(monacc,mask=monacc>0))
     &     stop 'unequal numbers of months'

      mon0=-1
      do m=1,12
        if(amon0==emonth(m)) mon0 = m
        if(amon0==emonth(m)) exit
      end do
      if(mon0 < 0) then
        mon0=index(alphabet,amon0(1:1))
        earth = .false.
      end if
      if(mon0 .le. 0) stop 'month not named [A-Z]...'

      mostr=''
      nmo = 0 ; incyr0 = .false.
      do mm=mon0,mon0+11 ! mon0+month_in_year - 1
        m = mm ; if(m > 12) m = m-12
        if(monacc(m).eq.0) cycle
        if(mm>12) incyr0 = .true.
        nmo = nmo + 1
        if(nmo.eq.1) then
          mostr(1:3)=amon0
        else
          if (earth) then
            mostr(nmo:nmo)=emonth(m)(1:1)
          else
            mostr(nmo:nmo)=alphabet(m:m)
          end if
        end if
      enddo

      if(nmo.eq.12) then
        mostr='ANN'
      elseif(nmo.eq.2) then
        mostr(3:3) = mostr(2:2)
        mostr(2:2) = '+'
      elseif(nmo.gt.3) then
        mostr(3:3) = mostr(nmo:nmo)
        mostr(2:2) = '-'
        write(*,*) 'labeling may be ambiguous !!'
      end if

      read(ayr0,*) yr0 ; ayr0='      '
      if(incyr0) yr0 = yr0 + 1
      if(yr0.lt.10000) then
        write(ayr0(1:4),'(i4.4)') yr0
      else
        write(ayr0,'(i6)') yr0
      end if

      read(ayr1,*) yr1 ; ayr1='      '
      if(amon=='JAN') yr1 = yr1 - 1
      if(.not.earth.and.amon(1:1)=='A') yr1 = yr1 - 1
      if(yr1.lt.10000) then
        write(ayr1(1:4),'(i4.4)') yr1
      else
        write(ayr1,'(i6)') yr1
      end if

      acc_period=''
      acc_period = mostr(1:3)//trim(adjustl(ayr0))
      if(yr1.gt.yr0)  acc_period = mostr(1:3)//trim(adjustl(ayr0))//
     *   '-'//trim(adjustl(ayr1))

c check for gaps
      ninc = count(monacc(2:12).gt.monacc(1:11))
      ndec = count(monacc(2:12).lt.monacc(1:11))
      if(ninc+ndec.gt.2) then
        write(6,*) 'gap'
      end if

      return
      end subroutine aperiod
