#include "rundeck_opts.h"
      subroutine modelE_mainDriver()
      use iso_c_binding

      interface
        subroutine libmodele_refaddr() bind(C)
        end subroutine
      end interface

!@sum Acquire configuration options from the command line and pass to
!@+ the model.
!@auth T. Clune
C**** Command line options
      logical :: qcRestart=.false.
      logical :: coldRestart=.false.
      integer, parameter :: MAX_LEN_IFILE = 32
      character(len=MAX_LEN_IFILE) :: iFile
#if ((! defined(COMPILER_NAG) ) && (! defined(COMPILER_G95) )) || (defined COMPILER_PGI)
      integer, external :: iargc
#endif
      integer :: i
      character(256) :: arg

      ! Print out command line arguments
      do i=1,iargc()
        call getarg(i, arg)
        print *,'ARG ', trim(arg)
      end do

      ! Set up to properly interpret stack traces
      call libmodele_refaddr()

      call read_options(qcRestart, coldRestart, iFile )
      call GISS_modelE(qcRestart, coldRestart, iFile)

      contains

      subroutine read_options(qcRestart, coldRestart, iFile )
!@sum Reads options from the command line
!@auth I. Aleinov
      implicit none
!@var qcRestart true if "-r" is present
!@var iFile is name of the file containing run configuration data
      logical, intent(inout) :: qcRestart
      logical, intent(inout) :: coldRestart
      character(*),intent(out)  :: ifile
      integer, parameter :: MAX_LEN_ARG = 80
      character(len=MAX_LEN_ARG) :: arg, value

      iFile = "";
      do
        call nextarg( arg, 1 )
        if ( arg == "" ) exit          ! end of args
        select case (arg)
        case ("-r")
          qcRestart = .true.
        case ("-cold-restart")
          coldRestart = .true.
        case ("-i")
          call nextarg( value, 0 )
          iFile=value
        ! new options can be included here
        case default
          print *,'Unknown option specified: ', arg
          print *,'Aborting...'
          call stop_model("Unknown option on a command line",255)
        end select
      enddo

      if (iFile == "") then
        print*, 'No configuration file specified on command line: '
        print*, 'Aborting ...'
        call stop_model("No configuration file on command line.",255)
      end if

      return
      end subroutine read_options

      end subroutine modelE_mainDriver
