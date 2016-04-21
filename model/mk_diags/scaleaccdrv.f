      program driver
      implicit none
      include 'netcdf.inc'
      integer :: status,ifids(2)
      character(len=200) :: ifile,remap_file
      character(len=80) :: accname
      integer :: i1,i2,nargs
      !integer, external :: iargc
      nargs = iargc()
      if(nargs.ne.2 .and. nargs.ne.3) then
         write(6,*)
     &   'usage: scaleacc acc-file acc_array_name[,name2] [remap_file]'
         stop
      endif
      call getarg(1,ifile)
      call getarg(2,accname)
      i1 = index(ifile,'.acc')
      if(i1.eq.0) i1 = index(ifile,'.subdd')
      i2 = index(ifile,'.nc')
      if(i1.eq.0 .or. i1.gt.i2) then
        stop 'expecting filename of the form *.acc*.nc or *.subdd*.nc'
      endif
      status = nf_open(trim(ifile),nf_nowrite,ifids(1))
      if(status.ne.nf_noerr) then
        write(6,*) 'nonexistent/non-netcdf input file ',trim(ifile)
        stop
      endif
      ifids(2) = -99
      if(nargs.eq.3) then
        call getarg(3,remap_file)
        status = nf_open(trim(remap_file),nf_nowrite,ifids(2))
        if(status.ne.nf_noerr) then
          write(6,*) 'nonexistent/non-netcdf remap file ',
     &          trim(remap_file)
          stop
        endif
      endif
      call scaleacc(ifids,accname,ifile)
      status = nf_close(ifids(1))
      if(ifids(2).ne.-99) status = nf_close(ifids(2))
      end program driver
