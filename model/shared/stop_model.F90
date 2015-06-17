#include "rundeck_opts.h"

module stop_model_mod

	implicit none

	interface
		subroutine stop_model_cb(message, retcode)
			implicit none
			!@var message an error message (reason to stop)
			character*(*), intent (in) :: message
			!@var retcode return code to be passed to the calling script
			integer, intent(in) :: retcode
		end subroutine stop_model_cb
	end interface

	procedure(stop_model_cb), pointer :: stop_model_ptr => stop_model_default

CONTAINS

subroutine stop_model_default( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
  use Dictionary_mod
  implicit none
!@var message an error message (reason to stop)
  character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
  integer, intent(in) :: retcode
  integer, parameter :: iu_err = 9
  integer :: rank
#ifdef USE_MPI
  integer :: mpi_err
#  ifdef MPI_DEFS_HACK
#  include "mpi_defs.h"
#  endif
#include "mpif.h"
#endif

#ifdef USE_MPI
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
#else
  rank =0
#endif
  ! skip writing status file for retcode<0
  if ( retcode >= 0 ) call write_run_status( message, retcode )
  if (rank == 0) then
    write (6,'(//2(" ",132("*")/))')
    write (6,*) ' Program terminated due to the following reason:'
    write (6,*) ' >>  ', message, '  <<'
    write (6,'(/2(" ",132("*")/))')
  endif

  call sys_flush(6)

  if ( retcode > 13 ) then
    write (0,*) 'Model crashed due to ',message
#ifdef USE_MPI
    !??? bad: the next line will prevent a job from terminating unless
    !???          all processors reach this point
    !??? bug: without it, jobs don't terminate even if
    !???          all processors reach this point
    call mpi_finalize(mpi_err)
    !??? hopefully, we can get rid of the above line soon
    call mpi_abort(MPI_COMM_WORLD, retcode, iu_err)
#else
    call sys_abort
#endif
  else
#ifdef USE_MPI
    call mpi_finalize(mpi_err)
#endif
    call exit_rc (0)
  endif

    end subroutine stop_model_default


#ifdef USE_FEXCEPTION
	subroutine stop_model_fexception(message, retcode)
		use fexception_mod

		!@var message an error message (reason to stop)
		character*(*), intent (in) :: message
		!@var retcode return code to be passed to the calling script
		integer, intent(in) :: retcode

		call throw(message, retcode)

	end subroutine stop_model_fexception
#endif

end module stop_model_mod

! =====================================================================

! A stub outside a module, to call the pointer.  If everywhere that used stop_model
! were willing to import stop_model_mod, then this would not be needd at all.
subroutine stop_model( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
	use stop_model_mod
implicit none
	!@var message an error message (reason to stop)
	character*(*), intent (in) :: message
	!@var retcode return code to be passed to the calling script
	integer, intent(in) :: retcode

	call stop_model_ptr(message, retcode)
end subroutine stop_model

subroutine throwException(message, retcode)
!@sum Either invokes pFUnit exception for testing or
!@+ stop_model() for run-time testing.
!@auth T. Clune
#ifdef USE_PFUNIT
  use pFUnit_mod, only: throw
#endif
  character(len=*), intent(in) :: message
  integer, intent(in) :: retcode

#ifdef USE_PFUNIT
  call throw(message)
#else
  call stop_model(message, retcode)
#endif
end subroutine throwException

subroutine exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
#if ( defined(COMPILER_NAG) )
  use f90_unix_proc
#endif
  implicit none
  integer, intent(IN) :: code !@var code return code set by user
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC) \
  || ( defined(MACHINE_MAC) && ! defined(COMPILER_XLF) )
       call exit(code) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
  call exit_(code)
#else
  none of supported architectures was specified.
  This will crash the compiling process.
#endif
  return
end subroutine exit_rc

