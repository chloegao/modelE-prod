#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      module sstmod

!@sum  Module sstmod contains the arrays/subroutines needed to prescribe
!@+    ocean surface temperature from input files.
!@+    Selection of time period occurs as per the comments concerning sst_yr
!@+    in init_sstmod.
!@auth Original Development Team
!@auth M. Kelley restructuring and netcdf-based input options

      use timestream_mod, only : timestream
      implicit none
      save

!@var sst sea surface temperature (C)
      real*8, dimension(:,:), allocatable :: sst

!@var SSTstream interface for reading and time-interpolating SST files
!@+   See general usage notes in timestream_mod.
!@+   Note regarding comparison of results to runs that use traditional I/O:
!@+   if SST datafiles contain monthly means and the piecewise parabolic
!@+   method is used for monthly->daily interpolation, the presence of OSST_eom
!@+   in the rundeck will prompt read_stream to read end-of-month values from
!@+   OSST_eom rather than computing them on the fly.  OSST and OSST_eom may
!@+   refer to the same file or directory.  The on-the-fly result will differ
!@+   due to roundoff effects.
      type(timestream) :: SSTstream

!@var tocean_4io an array for restart file compatibility with ML ocean
!@+   (see OCNML2.f for definition of its tocean array)
      real*8, dimension(:,:,:), allocatable :: tocean_4io

      logical :: osst_exists=.true.

      contains

      subroutine alloc_sstmod
!@sum alloc_sstmod allocates arrays in module sstmod
      use domain_decomp_atm, only : grid,getDomainBounds
      implicit none
      integer :: i_0h,i_1h,j_0h,j_1h,ier
      call getDomainBounds(grid,j_strt_halo=j_0h,j_stop_halo=j_1h)
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      allocate(sst(i_0h:i_1h,j_0h:j_1h))
      allocate(tocean_4io(3,i_0h:i_1h,j_0h:j_1h))
      sst = 0.
      tocean_4io = 0.
      end subroutine alloc_sstmod

      subroutine init_sstmod(atmocn)
!@sum init_sstmod initializes the SSTstream object
      use domain_decomp_atm, only : grid
      use timestream_mod, only : init_stream
      use model_com, only :  modelEclock, master_yr
      use exchange_types, only : atmocn_xchng_vars
      use dictionary_mod, only : get_param,is_set_param
      use filemanager, only : file_exists
      implicit none
      type(atmocn_xchng_vars) :: atmocn
      integer :: jyear,jday,sst_yr
      logical :: cyclic

      osst_exists = file_exists('OSST')
      if(.not.osst_exists) return

      if(is_set_param('sst_yr')) then
        ! If parameter sst_yr exists, SST data from that year is
        ! selected (only relevant if OSST is a multi-year dataset).
        call get_param( 'sst_yr', sst_yr )
      else
        ! Otherwise, sst_yr is set to ocean_yr or master_yr.
        call get_param( 'ocean_yr', sst_yr, default=master_yr )
      endif
      cyclic = sst_yr /= 0 ! sst_yr==0 implies transient mode.
      sst_yr = abs(sst_yr)
      call modelEclock%get(year=jyear, dayOfYear=jday)
      if(cyclic) jyear = sst_yr
      call init_stream(grid,SSTstream,'OSST','sst',-100d0,100d0,'ppm',
     &       jyear,jday,msk=atmocn%focean,cyclic=cyclic)
      end subroutine init_sstmod

      subroutine set_gtemp_sst(atmocn)
!@sum set_gtemp_sst copies sst into atmocn%gtemp
      use constant, only : tf
      use domain_decomp_atm, only : grid,getDomainBounds
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      use exchange_types, only : atmocn_xchng_vars
      implicit none
      type(atmocn_xchng_vars) :: atmocn
c
      integer :: i,j,j_0,j_1, i_0,i_1
      call getDomainBounds(grid,i_strt=i_0,i_stop=i_1)
      call getDomainBounds(grid,j_strt=j_0,j_stop=j_1)
      do j=j_0,j_1
      do i=i_0,atmocn%imaxj(j)
        if (atmocn%focean(i,j).gt.0) then
          atmocn%gtemp(i,j)=sst(i,j)
          atmocn%gtemp2(i,j)=tocean_4io(2,i,j) ! to preserve identical diagnostics
          atmocn%gtempr(i,j) =sst(i,j)+tf
#ifdef SCM
c         may specify ocean temperature for SCM
          if( SCMopt%Tskin )then
            atmocn%gtemp(I,J) = SCMin%Tskin - TF
            atmocn%gtempr(I,J) = SCMin%Tskin
          endif
#endif
        endif
      enddo
      enddo
      end subroutine set_gtemp_sst

      subroutine def_rsf_sstmod(fid)
!@sum  def_rsf_sstmod defines sstmod array structure in restart files
!@auth M. Kelley
!@ver  beta
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,tocean_4io,'tocean(d3,dist_im,dist_jm)')
      return
      end subroutine def_rsf_sstmod

      subroutine new_io_sstmod(fid,iaction)
!@sum  new_io_sstmod read/write sstmod arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        tocean_4io(1,:,:) = sst(:,:)
        call write_dist_data(grid, fid, 'tocean', tocean_4io, jdim=3)
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'tocean', tocean_4io, jdim=3)
        sst(:,:) = tocean_4io(1,:,:)
      end select
      return
      end subroutine new_io_sstmod

      end module sstmod

      subroutine read_sst(end_of_day,atmocn)
!@sum read_sst invokes procedures to read sea surface temperature from
!@+   input files and perform time interpolation
      use domain_decomp_atm, only : getDomainBounds,grid
      use model_com, only : itime,itimei
      use model_com, only :  modelEclock
      use resolution, only : im,jm
      use seaice, only : tfrez
      use timestream_mod, only : read_stream
      use sstmod, only : SSTstream,SST,osst_exists
      use exchange_types, only : atmocn_xchng_vars
      implicit none
      logical, intent(in) :: end_of_day
      type(atmocn_xchng_vars) :: atmocn
c
      real*8 :: tfo
      integer i,j
      integer :: jyear,jday

      integer :: j_0,j_1, i_0,i_1
      logical :: have_north_pole, have_south_pole

      if(.not.osst_exists) return

      call modelEclock%get(year=jyear, dayOfYear=jday)

      call getDomainBounds(grid,
     &         i_strt=i_0,i_stop=i_1,j_strt=j_0,j_stop=j_1,
     &         have_south_pole=have_south_pole,
     &         have_north_pole=have_north_pole)

      if(.not.(end_of_day.or.itime.eq.itimei)) return

C**** read and time-interpolate
      call read_stream(grid,SSTstream,jyear,jday,SST)

c**** bounds on sst
      do j=j_0,j_1
      do i=i_0,atmocn%imaxj(j)

        if (atmocn%focean(i,j).gt.0) then
          tfo=tfrez(atmocn%sss(i,j))
          if (sst(i,j).lt.tfo) sst(i,j)=tfo
        else
          sst(i,j) = 0.
        endif
      end do
      end do

c**** replicate values at pole
      if(have_north_pole) then
        if (atmocn%focean(1,jm).gt.0) then
          do i=2,im
            sst(i,jm)=sst(1,jm)
          end do
        end if
      end if
      if(have_south_pole) then
        if (atmocn%focean(1,1).gt.0) then
          do i=2,im
            sst(i,1)=sst(1,1)
          end do
        end if
      end if

      return
      end subroutine read_sst


      subroutine setup_ocean
      return
      end subroutine setup_ocean


      subroutine init_ocean(iniocean,istart,atmocn,dynsice)
!@sum init_OCEAN initializes ocean variables
!@auth Original Development Team
!@ver  1.0
      use domain_decomp_atm, only : grid, getDomainBounds
      use model_com, only : kocean,ioread
#ifdef TRACERS_WATER
      use oldtracer_mod, only : trw0
#endif
      use fluxes, only : atmice ! move uisurf,visurf init elsewhere
      use seaice, only : qsfix, osurf_tilt
      use ocnml, only : init_ocnml,set_gtemp_ocnml
      use sstmod, only : init_sstmod,set_gtemp_sst
      use pario, only : par_open, par_close
      use exchange_types, only : atmocn_xchng_vars,iceocn_xchng_vars
      use filemanager, only : file_exists
      implicit none
      logical, intent(in) :: iniocean  ! true if starting from ic.
      integer, intent(in) :: istart
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: dynsice ! not used here
!@var sss0 default sea surface salinity (psu)
      real*8, parameter :: sss0=34.7d0
      integer :: fid
      integer :: i,j
      integer :: i_0,i_1, j_0,j_1

      call getDomainBounds(grid,I_STRT=I_0,I_STOP=I_1)
      call getDomainBounds(grid,j_strt=j_0,j_stop=j_1)

      if (istart.le.0) then
        if(kocean.ge.1) call init_ODEEP(.false.)
        return
      end if

C**** Cold start
      if (istart.le.2 .and. file_exists('GIC')) then
        fid = par_open(grid,'GIC','read')
        call new_io_ocean (fid,ioread)
        call par_close(grid,fid)
      end if

c**** set fluxed arrays for oceans
      do j=j_0,j_1
      do i=i_0,i_1
        if (atmocn%focean(i,j).gt.0) then
          atmocn%sss(i,j) = sss0
#ifdef TRACERS_WATER
          atmocn%gtracer(:,i,j)=trw0()
#endif
        else
          atmocn%sss(i,j) = 0.
        end if
c**** for the time being assume zero surface velocities for drag calc
        atmocn%uosurf(i,j)=0. ; atmocn%vosurf(i,j)=0.
        atmice%uisurf(i,j)=0. ; atmice%visurf(i,j)=0.
c**** also zero out surface height variations
        atmocn%ogeoza(i,j)=0.
      end do
      end do
c**** keep salinity in sea ice constant for fixed-sst and qflux models
      qsfix = .true.
c**** make sure to use geostrophy for ocean tilt term in ice dynamics
c**** (if required). since ocean currents are zero, this implies no sea
c**** surface tilt term.
      osurf_tilt = 0


C**** 
      if (kocean.eq.0) then
        call set_gtemp_sst(atmocn)
        call init_sstmod(atmocn)
      else
        call set_gtemp_ocnml(atmocn)
        ! read ML depths and OHT
        call init_ocnml(iniOCEAN,istart,atmocn)
      endif

      return
      end subroutine init_ocean

      subroutine alloc_ocean
!@sum alloc_ocean calls allocation routines for either the
!@+     (1) prescribed ocean module (kocean=0)
!@+     (2) mixed-layer ocean module (kocean=1)
      use domain_decomp_atm, only : grid
      use sstmod, only  : alloc_sstmod
      use ocnml, only : alloc_ocnml
      use dictionary_mod, only : get_param
      use model_com, only : kocean
      implicit none
      call get_param('kocean',kocean)
      if(kocean.eq.0) then
        call alloc_sstmod
      else
        call alloc_ocnml
        call alloc_odeep(grid)
      endif
      end subroutine alloc_ocean

      subroutine daily_ocean(end_of_day,atmocn)
!@sum daily_ocean calls daily update routines for either the
!@+     (1) prescribed ocean module (kocean=0)
!@+     (2) mixed-layer ocean module (kocean=1)
      use model_com, only : kocean
      use sstmod, only : set_gtemp_sst
      use ocnml, only : daily_ocnml,set_gtemp_ocnml
      use exchange_types, only : atmocn_xchng_vars
      use fluxes, only : atmice
      implicit none
      logical, intent(in) :: end_of_day
      type(atmocn_xchng_vars) :: atmocn

      if (kocean.ge.1) then
        ! update prescribed ml depth, perform associated adjustments
        call daily_ocnml(end_of_day,atmocn,atmice)
        call set_gtemp_ocnml(atmocn)
      else
        ! update prescribed sst
        call read_sst(end_of_day,atmocn)
        call set_gtemp_sst(atmocn)
        call daily_seaice(end_of_day,atmocn,atmice)
      end if

      return
      end subroutine daily_ocean

      subroutine oceans(atmocn,iceocn,dynsice)
!@sum ocean calls routines to apply surface fluxes to either the
!@+     (1) prescribed ocean (kocean=0, no-op)
!@+     (2) mixed-layer ocean (kocean=1)
!@auth original development team
!@ver  1.0
      use model_com, only : kocean
      use diag_com, only : oa
      use domain_decomp_atm, only : grid,getDomainBounds
      use ocnml, only : run_ocnml,set_gtemp_ocnml
      use exchange_types, only : atmocn_xchng_vars,iceocn_xchng_vars
      use fluxes, only : atmice
      implicit none
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      type(iceocn_xchng_vars) :: dynsice ! not used here

      integer i,j, j_0,j_1, i_0,i_1

      call getDomainBounds(grid,i_strt=i_0,i_stop=i_1)
      call getDomainBounds(grid,j_strt=j_0,j_stop=j_1)

      if(kocean.ge.1) then
        ! surface fluxes affect predicted ocean temperature
        call run_ocnml(atmocn,atmice,iceocn)
        call set_gtemp_ocnml(atmocn)
      else
        ! add rvr e to surf. energy budget, set ice formation rate = 0
        do j=j_0,j_1
        do i=i_0,atmocn%imaxj(j)
          if (atmocn%focean(i,j).gt.0) then
            oa(i,j,4)=oa(i,j,4)+
     &         (atmocn%eflowo(i,j)+atmocn%egmelt(i,j))
            iceocn%dmsi(1,i,j)=0.
            iceocn%dmsi(2,i,j)=0.
            iceocn%dhsi(1,i,j)=0.
            iceocn%dhsi(2,i,j)=0.
            iceocn%dssi(1,i,j)=0.
            iceocn%dssi(2,I,J)=0.
#ifdef TRACERS_WATER
            iceocn%dtrsi(:,1,I,J)=0.
            iceocn%dtrsi(:,2,I,J)=0.
#endif
          end if
        end do
        end do
      endif
      return
      end subroutine oceans

      subroutine precip_oc(atmocn,iceocn)
!@sum  precip_oc driver for applying precipitation fluxes to mixed-layer ocean.
!@+    This routine could be folded into oceans, but exists separately for
!@+    historical/diagnostic reasons.
!@auth original development team
!@ver  1.0
      use diag_com, only : oa
      use model_com, only : kocean
      use ocnml, only : precip_ocnml,set_gtemp_ocnml
      use exchange_types, only : atmocn_xchng_vars,iceocn_xchng_vars
      use fluxes, only : atmice
      implicit none
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      where(atmocn%focean.gt.0.) oa(:,:,4) = oa(:,:,4)+atmocn%eprec(:,:)
      if(kocean.ge.1) then
        call precip_ocnml(atmocn,atmice,iceocn)
        call set_gtemp_ocnml(atmocn)
      endif
      return
      end subroutine precip_oc

      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      use model_com, only : kocean
      use sstmod, only : def_rsf_sstmod
      !use ocnml, only : def_rsf_ocnml
      use domain_decomp_atm, only : grid
      implicit none
      integer fid   !@var fid file id
      if(kocean.eq.0) then
        call def_rsf_sstmod(fid)
      else
        call def_rsf_ocnml(fid)
      endif
      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : kocean
      use sstmod, only : new_io_sstmod
      !use ocnml, only : new_io_ocnml
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      if(kocean.eq.0) then
        call new_io_sstmod(fid,iaction)
      else
        call new_io_ocnml(fid,iaction)
      endif
      return
      end subroutine new_io_ocean


      SUBROUTINE DIAGCO (M,atmocn)
!@sum  DIAGCO Keeps track of the ocean conservation properties
!@auth Gary Russell/Gavin Schmidt
      USE MODEL_COM, only : kocean
      USE DIAG_COM, only : icon_OCE
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
!@var M index denoting from where DIAGCO is called
      INTEGER, INTENT(IN) :: M
      type(atmocn_xchng_vars) :: atmocn
C****
C**** THE PARAMETER M INDICATES WHEN DIAGCO IS BEING CALLED
C****     (see DIAGCA)
      REAL*8, EXTERNAL :: conserv_OCE

C**** OCEAN POTENTIAL ENTHALPY
      IF (KOCEAN.ge.1) CALL conserv_DIAG(M,conserv_OCE,icon_OCE)
C****
      RETURN
      END SUBROUTINE DIAGCO

#ifdef TRACERS_WATER
      subroutine tracer_ic_ocean(atmocn)
!@sum tracer_ic_ocean initialise ocean tracer concentration
!@+   called only when tracers turn on
!@auth Gavin Schmidt
      USE MODEL_COM, only : itime
      USE OldTracer_mod, only : trw0,itime_tr0
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      INTEGER i,j,n
      do n=1,atmocn%ntm
        if (itime.eq.itime_tr0(n)) then
          do j=atmocn%j_0,atmocn%j_1
          do i=atmocn%i_0,atmocn%imaxj(j)
            if(atmocn%focean(i,j).gt.0) then
              atmocn%gtracer(n,i,j)=trw0(n)
            end if
          end do
          end do
        end if
      end do
      return
      end subroutine tracer_ic_ocean
#endif
