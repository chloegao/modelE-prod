#include "rundeck_opts.h"
      module flammability_com    
!@sum for routines to calculate flammability potential of surface 
!@+   vegetation. Optionally also altering tracer biomass sources.
!@auth Greg Faluvegi based on direction from Olga Pechony

      implicit none
      save

! main variables:
!@var flammability unitless flammability coefficient
!@var vegetation density read from file for flammability calculation
      real*8, allocatable, dimension(:,:) :: flammability,veg_density
!@param mfcc MODIS fire count calibration (goes with EPFC by veg type
!@+ params. Units=fires/day when multiplied by the unitless flammability
      real*8, parameter :: mfcc=1122.d0
      real*8, parameter :: missing=-1.d30 ! use it from CONST.f

! rest is for the running average:
!@dbparam allowFlammabilityReinit (default 1=YES) allows the
!@+ flammability to initialize to missing when Itime=ItimeI and
!@+ the averaging period to start from the beginning (thus no emissions
!@+ for nday_prec days.) Override this with allowFlammabilityReinit=0
!@+ in the rundeck, and values from the AIC file should be used (if 
!@+ present.) 
!@var maxHR_prec maximum number of sub-daily accumulations
!@param nday_prec number of days in running average for prec
!@var DRAfl daily running average of model prec for flammability
!@var ravg_prec period running average of model prec for flammability
!@var PRSfl period running sum of model prec for flammability
!@var HRAfl hourly running average of model prec for flammability
!@var iHfl "hourly" index for averages of model prec
!@var iDfl "daily"  index for averages of model prec
!@var i0fl ponter to current index in running sum of model prec
!@var first_prec whether in the first model averaging per. for prec
      integer, parameter :: nday_prec=30
      integer :: allowFlammabilityReinit = 1
      real*8, allocatable, dimension(:,:,:):: DRAfl
      real*8, allocatable, dimension(:,:)  :: ravg_prec,PRSfl,iHfl,iDfl,
     &                                        i0fl,first_prec
      real*8, allocatable, dimension(:,:,:):: HRAfl
!@var raP_acc accumulate running avg precip for SUBDD output
      real*8, allocatable, dimension(:,:):: raP_acc
      integer :: maxHR_prec

      end module flammability_com


      subroutine alloc_flammability(grid)
!@SUM  alllocates arrays whose sizes need to be determined
!@+    at run-time
!@auth Greg Faluvegi
      use domain_decomp_atm, only: dist_grid, getDomainBounds
      use dictionary_mod, only : get_param, is_set_param
      use model_com, only: dtsrc
      use TimeConstants_mod, only: SECONDS_PER_HOUR, HOURS_PER_DAY
      use flammability_com, only: flammability,veg_density,
     & first_prec,iHfl,iDfl,i0fl,DRAfl,ravg_prec,PRSfl,HRAfl,
     & nday_prec,maxHR_prec,raP_acc
 
      implicit none
      
      real*8 :: DTsrc_LOCAL
      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

      ! maxHR_prec is to be defined as the number of precip
      ! running-average steps that define one day. (The 1.d0
      ! in the line is left in to represent the number of calls 
      ! per DTsrc timestep. 
      ! I believe the real DTsrc is not available yet so I use
      ! a local copy from the database:

      DTsrc_LOCAL = DTsrc
      if(is_set_param("DTsrc"))call get_param("DTsrc",DTsrc_LOCAL)
      maxHR_prec = NINT(HOURS_PER_DAY*1.d0*SECONDS_PER_HOUR/DTsrc_LOCAL)

      allocate( flammability(I_0H:I_1H,J_0H:J_1H) )
      allocate( veg_density (I_0H:I_1H,J_0H:J_1H) )
      allocate( first_prec  (I_0H:I_1H,J_0H:J_1H) )
      allocate( iHfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( iDfl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( i0fl        (I_0H:I_1H,J_0H:J_1H) )
      allocate( DRAfl       (I_0H:I_1H,J_0H:J_1H,nday_prec) )
      allocate( ravg_prec   (I_0H:I_1H,J_0H:J_1H) )
      allocate( PRSfl       (I_0H:I_1H,J_0H:J_1H) )
      allocate( raP_acc     (I_0H:I_1H,J_0H:J_1H) )
      allocate( HRAfl       (I_0H:I_1H,J_0H:J_1H,maxHR_prec) )

      return
      end subroutine alloc_flammability


      subroutine init_flammability
!@sum initialize flamability, veg density, etc. for fire model
!@auth Greg Faluvegi based on direction from Olga Pechony
      use model_com, only: Itime,ItimeI
      use dictionary_mod, only: sync_param
      use flammability_com, only: flammability, veg_density, first_prec
     & ,missing,allowFlammabilityReinit
      use domain_decomp_atm,only: grid, getDomainBounds, 
     & am_i_root, readt_parallel
      use filemanager, only: openunit, closeunit, nameunit

      implicit none
      character*80 :: title,fname

      integer :: I_1H, I_0H, J_1H, J_0H, iu_data, n

      call getDomainBounds(grid,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      call getDomainBounds(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      call openunit('VEG_DENSE',iu_data,.true.,.true.)
      call readt_parallel(grid,iu_data,nameunit(iu_data),veg_density,1)
      call closeunit(iu_data)

      call sync_param
     &("allowFlammabilityReinit",allowFlammabilityReinit)

      if(Itime==ItimeI .and. allowFlammabilityReinit==1 )then
        first_prec(:,:)=1.d0
        flammability(:,:)=missing
      endif

      return
      end subroutine init_flammability


      subroutine io_flammability(kunit,iaction,ioerr)
!@sum  io_flammabilty reads and writes flammability variables to file
!@auth Greg Faluvegi (based on Jean Lerner io_tracer)
      use resolution, only : im,jm
      use domain_decomp_atm, only : grid
      use model_com, only: ioread,iowrite,irsfic,irsficno,irerun
      use domain_decomp_1d, only: get,am_i_root,
     &     pack_data,unpack_data
      use flammability_com, only: iHfl,iDfl,i0fl,first_prec,PRSfl,
     & DRAfl,HRAfl,maxHR_prec,nday_prec,ravg_prec,flammability,raP_acc

      implicit none

      integer :: kunit   !@var kunit unit number of read/write
      integer :: iaction !@var iaction flag for reading or writing to file
!@var ioerr 1 (or -1) if there is (or is not) an error in i/o
      integer, intent(inout) :: ioerr
!@var header character string label for individual records

      real*8, dimension(:,:), allocatable:: general_glob
      real*8, dimension(:,:,:), allocatable :: DRAfl_glob,HRAfl_glob
      integer :: itm
      character*80 :: header

      INTEGER :: J_0, J_1, J_1H, J_0H

      call getDomainBounds(grid, J_STRT=J_0,     J_STOP=J_1,
     &         J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      if(am_i_root()) allocate(
     &     general_glob(IM,JM),
     &     DRAfl_glob(IM,JM,nday_prec),
     &     HRAfl_glob(IM,JM,maxHR_prec)
     &     )

      SELECT CASE (IACTION)

      CASE (:IOWRITE) ! output to end-of-month restart file

       header='CALCULATE_FLAMMABILITY: DRAfl(i,j,days)'
        do itm=1,nday_prec
         call pack_data(grid,drafl(:,:,itm),drafl_glob(:,:,itm))
        end do
        if(am_i_root())write(kunit,err=10)header,drafl_glob
       header='CALCULATE_FLAMMABILITY: HRAfl(i,j,hours)'
        do itm=1,maxHR_prec
         call pack_data(grid,hrafl(:,:,itm),hrafl_glob(:,:,itm))
        end do
        if(am_i_root())write(kunit,err=10)header,hrafl_glob
       header='CALCULATE_FLAMMABILITY: PRSfl(i,j)'
        call pack_data(grid,prsfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: i0fl(i,j) (real)'
        call pack_data(grid,i0fl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: iDfl(i,j) (real)'
        call pack_data(grid,iDfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: iHfl(i,j) (real)'
        call pack_data(grid,iHfl(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: first_prec(i,j) (real)'
        call pack_data(grid,first_prec(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: ravg_prec(i,j)'
        call pack_data(grid,ravg_prec(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: flammability(i,j)'
        call pack_data(grid,flammability(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob
       header='CALCULATE_FLAMMABILITY: raP_acc(i,j)'
        call pack_data(grid,raP_acc(:,:),general_glob(:,:))
        if(am_i_root())write(kunit,err=10)header,general_glob

      CASE (IOREAD:)          ! input from restart file
        SELECT CASE (IACTION)
        CASE (ioread,irerun,irsfic,irsficno) ! restarts

          if(am_i_root())read(kunit,err=10)header,drafl_glob
          do itm=1,nday_prec
            call unpack_data(grid,drafl_glob(:,:,itm),drafl(:,:,itm))
          end do           
          if(am_i_root())read(kunit,err=10)header,hrafl_glob
          do itm=1,maxHR_prec
            call unpack_data(grid,hrafl_glob(:,:,itm),hrafl(:,:,itm))
          end do
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),prsfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),i0fl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),iDfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),iHfl(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),first_prec(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),ravg_prec(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),flammability(:,:))
          if(am_i_root())read(kunit,err=10)header,general_glob
          call unpack_data(grid,general_glob(:,:),raP_acc(:,:))

        END SELECT
      END SELECT

      call freemem
      return

 10   ioerr=1
      call freemem
      call stop_model('error in io_flammability',255) 
      return

      contains
      subroutine freemem
      if(am_i_root()) deallocate(general_glob,DRAfl_glob,HRAfl_glob)
      end subroutine freemem

      end subroutine io_flammability

#ifdef NEW_IO
      subroutine def_rsf_flammability(fid)
!@sum  def_rsf_flammability defines flammability array structure in 
!@+    restart files
!@auth Greg Faluvegi (directly from M. Kelley's def_rsf_lakes)
!@ver  beta
      use flammability_com
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id

      call defvar(grid,fid,drafl,'drafl(dist_im,dist_jm,nday_prec)')
      call defvar(grid,fid,hrafl,'hrafl(dist_im,dist_jm,maxHR_prec)')
      call defvar(grid,fid,prsfl,'prsfl(dist_im,dist_jm)')
      call defvar(grid,fid,i0fl,'i0fl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iDfl,'iDfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,iHfl,'iHfl(dist_im,dist_jm)') ! real
      call defvar(grid,fid,first_prec,'first_prec(dist_im,dist_jm)')
      call defvar(grid,fid,ravg_prec,'ravg_prec(dist_im,dist_jm)')
      call defvar(grid,fid,flammability,'flammability(dist_im,dist_jm)')
      call defvar(grid,fid,raP_acc,'raP_acc(dist_im,dist_jm)')

      return
      end subroutine def_rsf_flammability

      subroutine new_io_flammability(fid,iaction)
!@sum  new_io_flammability read/write arrays from/to restart files
!@auth Greg Faluvegi (directly from M. Kelley's new_io_lakes)
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use flammability_com
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid, fid, 'drafl', drafl )
        call write_dist_data(grid, fid, 'hrafl', hrafl )
        call write_dist_data(grid, fid, 'prsfl', prsfl )
        call write_dist_data(grid, fid, 'i0fl', i0fl )
        call write_dist_data(grid, fid, 'iDfl', iDfl )
        call write_dist_data(grid, fid, 'iHfl', iHfl )
        call write_dist_data(grid, fid, 'first_prec', first_prec )
        call write_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call write_dist_data(grid, fid, 'flammability', flammability )
        call write_dist_data(grid, fid, 'raP_acc', raP_acc )
      case (ioread)            ! input from restart file
        call read_dist_data(grid, fid, 'drafl', drafl )
        call read_dist_data(grid, fid, 'hrafl', hrafl )
        call read_dist_data(grid, fid, 'prsfl', prsfl )
        call read_dist_data(grid, fid, 'i0fl', i0fl )
        call read_dist_data(grid, fid, 'iDfl', iDfl )
        call read_dist_data(grid, fid, 'iHfl', iHfl )
        call read_dist_data(grid, fid, 'first_prec', first_prec )
        call read_dist_data(grid, fid, 'ravg_prec', ravg_prec )
        call read_dist_data(grid, fid, 'flammability', flammability )
        call read_dist_data(grid, fid, 'raP_acc', raP_acc )
      end select
      return
      end subroutine new_io_flammability
#endif /* NEW_IO */

      subroutine flammability_drv
!@sum driver routine for flammability potential of surface
!@+   vegetation calculation.
!@auth Greg Faluvegi based on direction from Olga Pechony
      use model_com, only: dtsrc
      use resolution, only : jm
      use atm_com, only : pedn
      use domain_decomp_atm,only: grid, getDomainBounds
      use flammability_com, only: flammability,veg_density,ravg_prec,
     & ravg_prec,iHfl,iDfl,i0fl,first_prec,HRAfl,DRAfl,PRSfl,missing,
     & raP_acc

      use fluxes, only: prec,atmsrf
      use constant, only: lhe
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use diag_com, only: ij_flam,aij=>aij_loc

      implicit none

      integer :: J_0S, J_1S, I_0H, I_1H, i, j
      logical :: have_south_pole, have_north_pole     
      real*8 :: qsat ! this is a function in UTILDBL.f
      real*8 :: tsurf,qsurf
                              
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)
      call getDomainBounds(grid,I_STRT_HALO=I_0H,I_STOP_HALO=I_1H)

      if(have_north_pole)flammability(I_0H:I_1H,JM)=missing
      if(have_south_pole)flammability(I_0H:I_1H,1) =missing

      do j=J_0S,J_1S
        do i=I_0H,I_1H

          ! update the precipitation running average, then:
          call prec_running_average(prec(i,j),ravg_prec(i,j), 
     &    iHfl(i,j),iDfl(i,j),i0fl(i,j),first_prec(i,j),HRAfl(i,j,:),
     &    DRAfl(i,j,:),PRSfl(i,j))
    
          ! for sub-daily diag purposes, accumulate the running avg:
          raP_acc(i,j)=raP_acc(i,j)+ravg_prec(i,j)

          ! if the first period has elapsed, calculate the flammability
          if(first_prec(i,j)==0) then
            ! could use tsurf,qsurf over land (from atmlnd%)?
            tsurf = atmsrf%tsavg(i,j)
            qsurf = atmsrf%qsavg(i,j)
            call calc_flammability(tsurf,SECONDS_PER_DAY*ravg_prec(i,j)/dtsrc,
     &           min(1.d0,qsurf/qsat(tsurf,lhe,pedn(1,i,j))),
     &           veg_density(i,j),flammability(i,j)) 
            endif
          ! update diagnostic
          aij(i,j,ij_flam)=aij(i,j,ij_flam)+flammability(i,j)

        enddo
      enddo
  
      return
      end subroutine flammability_drv

