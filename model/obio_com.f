#include "rundeck_opts.h"
      MODULE obio_com
!@sum  obio_com contains the parameters, arrays and definitions
!@+    necessary for the OceanBiology routines
!@auth NR

      USE obio_dim
      use ocalbedo_mod, only: nlt

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo
      use ocean, only : jm
#else
      USE hycom_dim_glob, only: kdm
#endif

      implicit none

c --- dobio       activate Watson Gregg's ocean biology code
      logical dobio
      data dobio/.true./
c

      real, ALLOCATABLE, DIMENSION(:,:)    :: tzoo2d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tfac3d,wshc3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: Fescav3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rmuplsr3d,rikd3d
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: acdom3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: gcmax         !cocco max growth rate
      real, ALLOCATABLE, DIMENSION(:,:)    :: pp2tot_day    !net pp total per day
      real, ALLOCATABLE, DIMENSION(:,:)    :: tot_chlo      !tot chlorophyl at surf. layer
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: rhs_obio      !rhs matrix
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: chng_by       !integr tendency for total C

#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
!      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrmflo_loc      ! riverine nitrate mass flow rate (kg/s)
      real, ALLOCATABLE, DIMENSION(:,:)    :: rnitrconc_loc      ! riverine nitrate concentration (kg/kg)
#endif
#ifdef DIC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdicconc_loc       ! riverine dic concentration (kg/kg)
#endif
#ifdef DOC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rdocconc_loc       ! riverine doc concentration (kg/kg)
#endif 
#ifdef SILI_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rsiliconc_loc      ! riverine silica concentration (kg/kg)
#endif
#ifdef IRON_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rironconc_loc      ! riverine iron concentration (kg/kg)
#endif
#ifdef POC_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: rpocconc_loc       ! riverine poc concentration (kg/kg)
#endif
#ifdef ALK_RUNOFF
      real, ALLOCATABLE, DIMENSION(:,:)    :: ralkconc_loc       ! riverine alkalinity concentration (mol/kg)
#endif

#endif

#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      real, ALLOCATABLE, DIMENSION(:,:) :: pCO2av, pCO2av_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav
      real, ALLOCATABLE, DIMENSION(:,:) :: pp2tot_dayav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: cexpav, cexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: caexpav, caexpav_loc
      real, ALLOCATABLE, DIMENSION(:,:) :: ao_co2fluxav,ao_co2fluxav_loc
#endif
      real, ALLOCATABLE, DIMENSION(:,:,:,:):: tracer

      integer :: nstep0=0

#ifdef OBIO_ON_GARYocean

      !test point
!!    integer, parameter :: itest=16, jtest=45    !equatorial Pacific                  2deg ocean
!!    integer, parameter :: itest=32, jtest=20    !southern ocean; Pacific          
      integer, parameter :: itest=1,  jtest=jm/2     !equator Pacific
#else
!     integer, parameter :: itest=(220,320) equator Atlant; (245,275) 0.6S;274.5E Nino3
!     integer, parameter :: itest=316, jtest=258    !257.5E;-50.7S
      integer, parameter :: itest=243, jtest=1      !equator,dateline
      real :: diag_counter
#endif


      integer, parameter :: EUZ_DEFINED=1


      real, parameter :: rlamz=1.0            !Ivlev constant
      real, parameter :: greff=0.25           !grazing efficiency     
      real, parameter :: drate=0.05/24.0      !phytoplankton death rate/hr
      real, parameter :: dratez1=0.1/24.0     !zooplankton death rate/hr
      real, parameter :: dratez2=0.5/24.0     !zooplankton death rate/hr
      real, parameter :: regen=0.25           !regeneration fraction

      real ::  obio_deltath,obio_deltat       !time steps in hours because all rates are in hrs
      real ::  co2mon(26,12)        !26 years 1979-2004, 12 months

      integer npst,npnd   !starting and ending array index for PAR
      data npst,npnd /3,17/

      real WtoQ(nlt)           !Watts/m2 to quanta/m2/s conversion

! reduced rank arrays for obio_model calculations
      integer ihra_ij

      real cexp, caexp
      real temp1d(kdm),dp1d(kdm),obio_P(kdm,ntyp+n_inert)
     .                 ,det(kdm,ndet),car(kdm,ncar),avgq1d(kdm)
     .                 ,gcmax1d(kdm),saln1d(kdm),p1d(kdm+1)
     .                 ,alk1d(kdm),flimit(kdm,nchl,5)

      real atmFe_ij,covice_ij
      integer inwst,inwnd,jnwst,jnwnd     !starting and ending indices 
                                          !for daylight 
                                          !in i and j directions
      real acdom(kdm,nlt)                 !absorption coefficient of CDOM
      real P_tend(kdm,ntyp+n_inert)       !bio tendency (dP/dt)

#ifdef TRACERS_Alkalinity
      real A_tend(kdm)
#endif
#ifdef OBIO_RUNOFF

#ifdef NITR_RUNOFF
      real rnitrconc_ij
!     .    , rnitrmflo_ij
#endif
#ifdef DIC_RUNOFF
      real rdicconc_ij
#endif
#ifdef DOC_RUNOFF
      real rdocconc_ij
#endif
#ifdef SILI_RUNOFF
      real rsiliconc_ij
#endif
#ifdef IRON_RUNOFF
      real rironconc_ij
#endif
#ifdef POC_RUNOFF
      real rpocconc_ij
#endif
#ifdef ALK_RUNOFF
      real ralkconc_ij
#endif
#endif

      real rmuplsr(kdm,nchl)                  !growth+resp 
      real D_tend(kdm,ndet)                   !detrtial tendency
      real obio_ws(kdm+1,nchl)                !phyto sinking rate
      real tfac(kdm)                          !phyto T-dependence
      real pnoice(kdm)                        !pct ice-free
      real wsdet(kdm+1,ndet)                  !detrital sinking rate
      real rikd(kdm,nchl)                     !photoadaption state
      real tzoo                               !herbivore T-dependence
      real Fescav(kdm)                        !iron scavenging rate

C if NCHL_DEFINED > 3
      real wshc(kdm)                          !cocco sinking rate
C endif

      real :: C_tend(kdm,ncar)                !carbon tendency
      real :: pCO2_ij                         !partial pressure of CO2
      real :: gro(kdm,nchl)                   !realized growth rate
      integer :: day_of_month, hour_of_day

      real :: rhs(kdm,ntrac,17)         !secord arg-refers to tracer definition 
                                        !we are not using n_inert (always ntrac-1)
                                        !second argument refers to process that
                                        !contributes to tendency 

      real :: pp2_1d(kdm,nchl)          !net primary production

      real*8 :: co2flux
      integer kzc
      real*8 :: carb_old,iron_old    !prev timesetep total carbon inventory

#ifdef restoreIRON
!per year change dI/I, + for sink, - for source
      !!!real*8 :: Iron_BC = 0.002
      real*8 :: Iron_BC = -0.005
#endif

      real*8, dimension(:, :, :), allocatable :: ze

      contains

      subroutine build_ze
#ifdef OBIO_ON_GARYocean
      use oceanr_dim, only: ogrid
      use oceanres, only: kdm=>lmo
      use ocean, only : zoe=>ze
#else
      use hycom_dim, only: ogrid, kdm
      USE hycom_arrays, only : dpinit
      USE hycom_scalars, only: onem
#endif
      implicit none
      integer :: k

      if (.not.allocated(ze)) allocate(ze(ogrid%i_strt:ogrid%i_stop,
     &                                ogrid%j_strt:ogrid%j_stop, 0:kdm))
      do k=0, kdm
#ifdef OBIO_ON_GARYocean
        ze(:, :, k)=zoe(k)
#else
        if (k>0) then
          ze(:, :, k)=dpinit(ogrid%i_strt:ogrid%i_stop,
     &                 ogrid%j_strt:ogrid%j_stop, k)/onem+ze(:, :, k-1)
        else
          ze(:, :, k)=0
        endif
#endif
      end do
      end subroutine build_ze

      END MODULE obio_com

!------------------------------------------------------------------------------
      subroutine alloc_obio_com
      USE obio_com
      USE obio_dim

#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only :kdm=>lmo
#else
      USE hycom_dim, only: kdm,ogrid
#endif

      implicit none

c**** Extract domain decomposition info
      INTEGER :: j_0,j_1,i_0,i_1

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP


      ALLOCATE(tracer(i_0:i_1,j_0:j_1,kdm,ntrac))

      call alloc_obio_forc

      ALLOCATE(tzoo2d(i_0:i_1,j_0:j_1))
      ALLOCATE(wshc3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(Fescav3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(rmuplsr3d(i_0:i_1,j_0:j_1,kdm,nchl),
     &            rikd3d(i_0:i_1,j_0:j_1,kdm,nchl))
      ALLOCATE(acdom3d(i_0:i_1,j_0:j_1,kdm,nlt))
      ALLOCATE(tfac3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(gcmax(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(pp2tot_day(i_0:i_1,j_0:j_1))
      ALLOCATE(tot_chlo(i_0:i_1,j_0:j_1))
      ALLOCATE(rhs_obio(i_0:i_1,j_0:j_1,ntrac,17))
      ALLOCATE(chng_by(i_0:i_1,j_0:j_1,14))

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      ALLOCATE(rnitrmflo_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(rnitrconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef DIC_RUNOFF
      ALLOCATE(rdicconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef DOC_RUNOFF
      ALLOCATE(rdocconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef SILI_RUNOFF
      ALLOCATE(rsiliconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef IRON_RUNOFF
      ALLOCATE(rironconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef POC_RUNOFF
      ALLOCATE(rpocconc_loc(i_0:i_1,j_0:j_1))
#endif
#ifdef ALK_RUNOFF
      ALLOCATE(ralkconc_loc(i_0:i_1,j_0:j_1))
#endif
#endif

#ifndef OBIO_ON_GARYocean   /* NOT for Russell ocean */
      ALLOCATE(pCO2av(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pp2tot_dayav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(cexpav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(caexpav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(pCO2av_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(pp2tot_dayav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(cexpav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(caexpav_loc(i_0:i_1,j_0:j_1))
      ALLOCATE(ao_co2fluxav(ogrid%im_world,ogrid%jm_world))
      ALLOCATE(ao_co2fluxav_loc(i_0:i_1,j_0:j_1))
#endif

      end subroutine alloc_obio_com

!------------------------------------------------------------------------------

#ifdef OBIO_ON_GARYocean

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE OCEANR_DIM, only : grid=>ogrid
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,nstep0
     &     ,tracer,pp2tot_day
      use pario, only : defvar
      implicit none
      integer fid   !@var fid file id
      integer :: n

      call defvar(grid,fid,nstep0,'obio_nstep0')
      call defvar(grid,fid,avgq,'avgq(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gcmax,'gcmax(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,tirrq3d,'tirrq3d(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,ihra,'ihra(dist_imo,dist_jmo)')
      call defvar(grid,fid,pp2tot_day,'pp2tot_day(dist_imo,dist_jmo)')
      return
      end subroutine def_rsf_obio

      subroutine new_io_obio(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      USE OCEANR_DIM, only : grid=>ogrid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      use model_com, only : nstep=>itime
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,nstep0
     &     ,tracer,pp2tot_day
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'obio_nstep0',nstep0)
        call write_dist_data(grid,fid,'avgq',avgq)
        call write_dist_data(grid,fid,'gcmax',gcmax)
        call write_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call write_dist_data(grid,fid,'ihra',ihra)
        call write_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'obio_nstep0',nstep0,
     &       bcast_all=.true.)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'ihra',ihra)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      end select
      return
      end subroutine new_io_obio

      subroutine new_io_obio_inicond
      USE obio_com, only: tracer
      use ocn_tracer_com, only : tracerlist, ocn_tracer_entry
      use ocean, only : lmo,lmm,ze,zmid,focean
      use ocean, only : im,jm
      use oceanr_dim, only : grid=>ogrid
      use pario, only : par_open,par_close
     &     ,read_data,read_dist_data,get_dimlens
      use domain_decomp_1d, only : halo_update
      implicit none

      integer i,j,l,lm,lm_in,lmo_in,n,fid,ii,jj
      logical :: need_zregrid
      real*8, allocatable :: z_in(:)
      real*8, dimension(:,:,:), allocatable :: arr_in,arr_tmp
      integer :: dlens(7),ndims
      integer :: i_0,i_1,j_0,j_1
      integer :: i_0h,i_1h,j_0h,j_1h
      type(ocn_tracer_entry), pointer :: entry
      interface
        Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)
        Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM)
        Real*8, optional :: RZ(LM)
        end Subroutine VLKtoLZ 
      end interface

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo


      fid = par_open(grid,'obio_inicond','read')

      call get_dimlens(grid,fid,'Nitr',ndims,dlens)
      lmo_in = dlens(3)

      allocate(z_in(lmo_in),arr_in(i_0h:i_1h,j_0h:j_1h,lmo_in))
      allocate(arr_tmp(0:im+1,j_0h:j_1h,lmo)) ! temporary until 2D decomp
      arr_in = 0.

      if(lmo_in == lmo) then
        z_in = zmid ! default
        call read_data(grid,fid,'z',z_in,bcast_all=.true.)
        need_zregrid = .not. all(abs(z_in-zmid) < 1d0)
      else
        call read_data(grid,fid,'z',z_in,bcast_all=.true.)
        need_zregrid = .true.
      endif

      tracer(:,:,:,:) = -9999.
      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        call read_dist_data(grid,fid,trim(entry%trname),arr_in)
        if(need_zregrid) then
          do j=j_0,j_1
          do i=i_0,i_1
            if(focean(i,j).le.0) cycle
            lm = lmm(i,j)
            call VLKtoLZ(lmo_in,lm,z_in,ze,
     &           arr_in(i,j,:),tracer(i,j,:,n))
          enddo
          enddo
        else
          tracer(:,:,:,n) = arr_in
        endif
        arr_tmp(1:im,:,:) = tracer(:,:,:,n)
        do j=j_0,j_1
          arr_tmp(0,j,:) = arr_tmp(im,j,:)
          arr_tmp(im+1,j,:) = arr_tmp(1,j,:)
        enddo
        call halo_update(grid,arr_tmp)
        do j=j_0,j_1
        do i=i_0,i_1
          lm = lmm(i,j)
          do l=lm+1,lmo
            tracer(i,j,l,n) = 0.
          enddo
          if(focean(i,j).le.0.) cycle
          do lm_in=0,lmo-1
            if(arr_tmp(i,j,1+lm_in).lt.0.) exit
          enddo
          if(lm_in.ge.lm) cycle
          do l=lm_in+1,lm            
c            do jj=j-1,j+1
            do jj=max(1,j-1),min(jm,j+1) ! temporary limits
            do ii=i-1,i+1
              if(jj.eq.j .and. ii.eq.i) cycle
              if(arr_tmp(ii,jj,l).lt.0.) cycle
              tracer(i,j,l,n) = arr_tmp(ii,jj,l)
            enddo
            enddo
          enddo
          do lm_in=0,lmo-1
            if(tracer(i,j,1+lm_in,n).lt.0.) exit
          enddo
          if(lm_in.ge.lm) cycle
          do l=lm_in+1,lm
            tracer(i,j,l,n) = tracer(i,j,lm_in,n)
          enddo
        enddo
        enddo
      enddo

      call par_close(grid,fid)
      deallocate(arr_in,arr_tmp)

      return
      end subroutine new_io_obio_inicond

#else

      subroutine obio_set_data_after_archiv
      USE obio_com, only:
     .     diag_counter
     .    ,ao_co2fluxav_loc, pp2tot_dayav_loc
     .    ,pCO2av_loc, cexpav_loc
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc
#endif
      implicit none
      diag_counter=0
      ao_co2fluxav_loc=0
      pCO2av_loc = 0
      pp2tot_dayav_loc = 0
      cexpav_loc = 0
#ifdef TRACERS_Alkalinity
      caexpav_loc = 0
#endif
      end subroutine obio_set_data_after_archiv

      subroutine obio_gather_before_archive
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      USE obio_com, only:
     .     ao_co2fluxav_loc, ao_co2fluxav
     .    ,pCO2av_loc, pCO2av
     .    ,pp2tot_dayav_loc, pp2tot_dayav
     .    ,cexpav_loc, cexpav
#ifdef TRACERS_Alkalinity
     .    ,caexpav_loc, caexpav
#endif
      implicit none 

      call pack_data(ogrid, ao_co2fluxav_loc, ao_co2fluxav)
      call pack_data(ogrid, pCO2av_loc, pCO2av)
      call pack_data(ogrid, pp2tot_dayav_loc, pp2tot_dayav)
      call pack_data(ogrid, cexpav_loc, cexpav)
#ifdef TRACERS_Alkalinity
      call pack_data(ogrid, caexpav_loc, caexpav)
#endif
      return
      end subroutine obio_gather_before_archive

      subroutine def_rsf_obio(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE HYCOM_DIM, only : grid=>ogrid
      use pario, only : defvar
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,pCO2av=>pCO2av_loc,pp2tot_day,
     &     ao_co2fluxav=>ao_co2fluxav_loc,
     &     pp2tot_dayav=>pp2tot_dayav_loc,
     &     cexpav=>cexpav_loc, diag_counter,nstep0
      implicit none
      integer fid   !@var fid file id
      character(len=14) :: str2d
      character(len=18) :: str3d
      character(len=20) :: str3d2
      str2d ='(idm,dist_jdm)'
      str3d ='(idm,dist_jdm,kdm)'
      str3d2='(idm,dist_jdm,kdmx2)'

      call defvar(grid,fid,nstep0,'obio_nstep0')
      call defvar(grid,fid,diag_counter,'obio_diag_counter')
      call defvar(grid,fid,avgq,'avgq'//str3d)
      call defvar(grid,fid,gcmax,'gcmax'//str3d)
      call defvar(grid,fid,tirrq3d,'tirrq3d'//str3d)
      call defvar(grid,fid,ihra,'ihra'//str2d)
      call defvar(grid,fid,pCO2av,'pCO2av'//str2d)
      call defvar(grid,fid,pp2tot_dayav,'pp2tot_dayav'//str2d)
      call defvar(grid,fid,ao_co2fluxav,'ao_co2fluxav'//str2d)
      call defvar(grid,fid,cexpav,'cexpav'//str2d)
      call defvar(grid,fid,pp2tot_day,'pp2tot_day'//str2d)

      return
      end subroutine def_rsf_obio

      subroutine new_io_obio(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      USE HYCOM_DIM, only : grid=>ogrid
      USE obio_forc, only : avgq,tirrq3d,ihra
      USE obio_com, only : gcmax,pCO2av=>pCO2av_loc,pp2tot_day,
     &     ao_co2fluxav=>ao_co2fluxav_loc,
     &     pp2tot_dayav=>pp2tot_dayav_loc,
     &     cexpav=>cexpav_loc, diag_counter,nstep0
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'obio_nstep0',nstep0)
        call write_data(grid,fid,'obio_diag_counter',diag_counter)
        call write_dist_data(grid,fid,'avgq',avgq)
        call write_dist_data(grid,fid,'gcmax',gcmax)
        call write_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call write_dist_data(grid,fid,'ihra',ihra)
        call write_dist_data(grid,fid,'pCO2av',pCO2av)
        call write_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call write_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call write_dist_data(grid,fid,'cexpav',cexpav)
        call write_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'obio_nstep0',nstep0,
     &       bcast_all=.true.)
        call read_data(grid,fid,'obio_diag_counter',diag_counter)
        call read_dist_data(grid,fid,'avgq',avgq)
        call read_dist_data(grid,fid,'gcmax',gcmax)
        call read_dist_data(grid,fid,'tirrq3d',tirrq3d)
        call read_dist_data(grid,fid,'ihra',ihra)
        call read_dist_data(grid,fid,'pCO2av',pCO2av)
        call read_dist_data(grid,fid,'pp2tot_dayav',pp2tot_dayav)
        call read_dist_data(grid,fid,'ao_co2fluxav',ao_co2fluxav)
        call read_dist_data(grid,fid,'cexpav',cexpav)
        call read_dist_data(grid,fid,'pp2tot_day',pp2tot_day)
      end select
      return
      end subroutine new_io_obio
#endif /* which ocean */
