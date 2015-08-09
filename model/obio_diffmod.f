#include "rundeck_opts.h"

      module obio_diffmod
      implicit none
      private

      public :: obio_listDifferences

      contains

      subroutine obio_listDifferences(operation, phase)
!@sum This routine checks for any changes against the previous state of both the 
!@+   tracer array and the dpinit array.   Each tracer is reported separately.
      
#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo, kdm=>lmo
      use oceanr_dim, only : ogrid
      use obio_com, only: tracer
      use ocn_tracer_com, only : tracerlist
       USE MODEL_COM,  only : nstep=>itime
#else
      use hycom_dim_glob, only : idm, kdm
      use hycom_dim, only: ogrid
      use hycom_scalars, only: nstep
      use hycom_arrays, only: tracer
      use hycom_arrays, only: dpinit
      use hycom_dim_glob, only: numTracers => ntrcr
#endif
      use domain_decomp_1d, only: getDomainBounds, am_i_root

      character(len=*), intent(in) :: operation
      character(len=*), intent(in) :: phase

      logical, save :: init = .false.
      real*8, allocatable, save :: previousTracers(:,:,:,:)
      real*8, allocatable, save :: previousdpinit(:,:,:)

      integer :: j_0, j_1, j_0h, j_1h
      integer :: iTracer
      character(len=50) :: name
#ifdef OBIO_ON_GARYocean
      integer :: numtracers
#endif

#ifdef OBIO_ON_GARYocean
      numtracers=tracerlist%getsize()
#endif
      call getDomainBounds(ogrid, j_strt = j_0, j_stop = j_1,
     &     j_strt_halo=j_0h, j_stop_halo = j_1h)

      if (.not. init) then
        init = .true.

        allocate(previousTracers(idm, j_0h:j_1h, kdm, numTracers))
        previousTracers = tracer

#ifndef OBIO_ON_GARYocean
        allocate(previousdpinit(idm,j_0h:j_1h, kdm))
        previousdpinit = dpinit
#endif
        return ! nothing to compare on the 1st trip
      end if

      select case (trim(phase))
      case ('before')

        previousTracers = tracer
#ifndef OBIO_ON_GARYocean
        previousdpinit = dpinit
#endif

      case default

        if (am_i_root()) 
     &       print*, 'obio_listDifferences for call: ', trim(operation)
        do iTracer = 1, numTracers
          if (am_i_root()) write(name,'(a,1x,a,i10,a,1x,i3.0)')
     .         trim(operation),',nstep = ',nstep,': tracer',iTracer
          call spotDiff3D(name, tracer(:,:,:,iTracer), 
     &         previousTracers(:,:,:,iTracer))
        end do

#ifndef OBIO_ON_GARYocean
        if (am_i_root()) write(name,'(a,1x,a,i10,a,1x)')
     .         trim(operation),',nstep = ',nstep,': dpinit'
        call spotDiff3D(name, dpinit, previousdpinit)
#endif
      end select

      end subroutine obio_listDifferences

      subroutine spotdiff3D(name, array, previous)
!@sum 3D version of similar routine from Rainer, which reports
!@+   locations of min and max differences in an array from the
!@+   previous call.
!@auth T. Clune <Thomas.L.Clune@nasa.gov>
      use domain_decomp_1d, only: am_i_root, getDomainBounds

c
c --- this routine compares 'array' with an earlier version of 'array'
c --- saved during a previous call to this routine
c
      character(len=*), intent(in) :: name
      real*8, intent(in) :: array(:,:,:)
      real*8, intent(out) :: previous(:,:,:)
      real*8               :: valMax, valMin
      integer            :: ijkAtMax(3), ijkAtMin(3)

      integer :: nk

c
      call getLocalMaxMin(array, previous, 
     &     ijkAtMax, valMax, ijkAtMin, valMin)
      call getGlobalExtreme(valMax, ijkAtMax, 'max')
      call getGlobalExtreme(valMin, ijkAtMin, 'min')

      if (am_i_root()) then
        print*,'  ', trim(name), ':'
        print 100,'largest pos change', valMax,' at i,j,k=', ijkAtMax
        print 100,'largest neg change', valMin,' at i,j,k=', ijkAtMin
      end if

 100  format (5x,a,es11.2,a,3i5)
      return

      contains

      subroutine getLocalMaxMin(a, b, 
     &     ijkAtMax, valMax, ijkAtMin, valMin)
#ifdef OBIO_ON_GARYocean
      use oceanres,  only: idm=>imo, kdm=>lmo
      use oceanr_dim, only : ogrid
#else
      use hycom_dim_glob, only : idm, kdm
      use hycom_dim, only: ogrid
#endif
      real*8, intent(in) :: a(:,ogrid%j_strt_halo:,:)
      real*8, intent(inout) :: b(:,ogrid%j_strt_halo:,:)
      integer, intent(out) :: ijkAtMax(3)
      real*8, intent(out) :: valMax
      integer, intent(out) :: ijkAtMin(3)
      real*8, intent(out) :: valMin

      integer :: i, j, k
      real*8 :: diff
      integer :: j_0, j_1

      call getDomainBounds(ogrid, j_strt=j_0, j_stop=j_1)

      valMax=-1.e33
      valMin=+1.e33

      do k = 1, kdm
        do j = j_0, j_1
          do i = 1, idm
            diff = a(i,j,k)-b(i,j,k)
            if      (diff > valMax) then
              valMax = diff
              ijkAtMax = (/ i, j, k /)
            else if (diff < valMin) then
              valMin = diff
              ijkAtMin = (/ i, j, k /)
            end if
            ! save "a" for next use
            b(i,j,k)=a(i,j,k)
          end do
        end do
      end do

      end subroutine getLocalMaxMin

      subroutine getGlobalExtreme(value, ijk, operation)
!@sum This routine uses a hack in MPI to pair real values with an integer 
!@+   (stored as a real).   The MPI_Reduce operation then returns the
!@+   global max (or min) and the associated index.
!@+   We must pass 3 such "pairs" to MPI to get each of i, j, and k.
!@auth T. Clune <Thomas.L.Clune@nasa.gov>

      real*8, intent(inout) :: value
      integer, intent(inout) :: ijk(3)
      character(len=*), intent(in) :: operation
      
#ifdef USE_MPI
      include 'mpif.h'
#endif
      integer :: ierr
      real*8 :: inBuffer(2,3) ! 3 pairs of form {value, index}
      real*8 :: outBuffer(2,3) ! 3 pairs of form {value, index}

      integer :: mpiOperation
#ifdef USE_MPI
      ! use MPI_Reduce to find maxval and associated indices
      inBuffer(1,1:3) = value
      inBuffer(2,1:3) = ijk

      select case (operation)
      case ('max','MAX','Max')
        mpiOperation = MPI_MAXLOC
      case ('min','MIN','Min')
        mpiOperation = MPI_MINLOC
      end select

      call MPI_REDUCE( inBuffer, outBuffer, 30, MPI_2DOUBLE_PRECISION, 
     &     mpiOperation, 0, MPI_COMM_WORLD, ierr ); 

      if (am_i_root()) then
        ! copy results for output
        value = outBuffer(1,1)
        ijk = outBuffer(2,1:3)
      end if
#endif
      end subroutine getGlobalExtreme
      
      end subroutine spotdiff3D

      end module obio_diffmod
