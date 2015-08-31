#include "rundeck_opts.h"

#if defined(CUBED_SPHERE) || defined(NEW_IO)
#else
#define USE_ATM_GLOBAL_ARRAYS
#endif

      subroutine setup_ocean
      return
      end subroutine setup_ocean

      SUBROUTINE init_OCEAN(iniOCEAN,istart,atmocn,dynsice)
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast
      USE SEAICE, only : osurf_tilt
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : delt1, salmin
     &  , nstep0, nstep, time0, time, itest, jtest
     &  , iocnmx, brntop, brnbot, ocnmx_factor_s, ocnmx_factor_t
     &  , diapyn, diapyc, jerlv0, thkdff
     &  , bolus_biharm_constant, bolus_laplc_constant
     &  , bolus_laplc_exponential

      USE HYCOM_ARRAYS_GLOB, only: scatter_hycom_arrays
      USE HYCOM_CPLER, only : agrid,tempro2a, ssto2a

      USE hycom_arrays_glob_renamer, only : temp_loc,saln_loc
      USE HYCOM_ATM, only : alloc_hycom_atm
      USE Dictionary_mod
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      implicit none

      logical, intent(in) :: iniOCEAN
      integer, intent(in) :: istart
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: dynsice ! not used here

      integer i,j,ia,ja

      atmocn%need_eflow_gl = .true. ! tell atm that eflow_gl is needed

      agrid => atmocn%grid ! grid used for pack_data/unpack_data in cpler

      aJ_0 = atmocn%J_0
      aJ_1 = atmocn%J_1
      aI_0 = atmocn%I_0
      aI_1 = atmocn%I_1

      aJ_0H = atmocn%J_0H
      aJ_1H = atmocn%J_1H
      aI_0H = atmocn%I_0H
      aI_1H = atmocn%I_1H

      call alloc_hycom_atm(atmocn)

#ifdef CUBED_SPHERE /* should be done for latlon atm also */
C**** Make sure to use geostrophy for ocean tilt term in ice dynamics
C**** (hycom ocean dynamics does not feel the weight of sea ice).
      osurf_tilt = 0
#endif

      call sync_param("itest", itest)
      call sync_param("jtest", jtest)
      call sync_param("iocnmx", iocnmx)
      call sync_param("brntop", brntop)
      call sync_param("brnbot", brnbot)
      call sync_param("ocnmx_factor_s", ocnmx_factor_s)
      call sync_param("ocnmx_factor_t", ocnmx_factor_t)
      call sync_param("diapyn", diapyn)
      call sync_param("diapyc", diapyc)
      call sync_param("jerlv0", jerlv0)
      call sync_param("thkdff", thkdff)
      call sync_param("bolus_biharm_constant",  bolus_biharm_constant)
      call sync_param("bolus_laplc_constant",   bolus_laplc_constant)
      call sync_param("bolus_laplc_exponential",bolus_laplc_exponential)

      if (iocnmx.ge.0.and.iocnmx.le.2 .or. iocnmx.eq.5 .or. iocnmx.eq.6) 
     .                                                              then
        call inikpp
      elseif (iocnmx.eq.3 .or. iocnmx.eq.7) then
        call inigis
      else
         stop 'wrong: need to choose one ocean mixing scheme'
      endif
c
      if (AM_I_ROOT()) then ! work on global grids here
        if ((bolus_laplc_constant*bolus_laplc_exponential==1)
     . .or. (bolus_laplc_constant*bolus_biharm_constant==1)
     . .or. (bolus_laplc_exponential*bolus_biharm_constant==1)) then
          print *,' wrong bolus setting: only one can be true'
          stop 'wrong bolus setting: only one can be true'
        elseif (bolus_laplc_constant==0 .and. bolus_laplc_exponential==0
     .    .and. bolus_biharm_constant==0) then
          stop 'wrong bolus setting: one has to be true'
        end if
c
css   if (istart.eq.2 .or. nstep0.eq.0) call geopar
      call inicon
c
c --- increase temp by 2 deg
c     do 21 j=1,jj
c     do 21 l=1,isp(j)
c     do 21 i=ifp(j,l),ilp(j,l)
c       if (latij(i,j,3).lt.-65..and.lonij(i,j,3).le.5.) then !lat[-90:90],lon[0:360]
c         p(i,j,1)=0.
c         do k=1,kk
c         p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c         temp(i,j,   k)=temp(i,j,   k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemm)
c         temp(i,j,kk+k)=temp(i,j,kk+k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemm)
c         enddo
c       endif
c21   continue
c


!!! I guess I had a good reason for commenting this out... 
!!! (probably should done in inicon) IA

!!!! the following is already done in inicon
!!      if (istart.gt.2) then               !istart=2 has done this in inirfn
!!      DO J=1,JM
!!      DO I=1,IM
!!        IF (FOCEAN(I,J).gt.0.) THEN
!!          GTEMP(1,1,I,J)=asst(I,J)
!!! GTEMPR ??
!!#ifdef TRACERS_GASEXCH_ocean
!!        do nt=1,ntm
!!        GTRACER(nt,1,I,J)=atrac(I,J,nt)
!!        enddo
!!#endif
!!        END IF
!!      END DO
 !!     END DO
!!      endif

      endif ! AM_I_ROOT

      call scatter_hycom_arrays

!!! hack needed for serial inicon
      CALL broadcast(ogrid, delt1 )

      CALL broadcast(ogrid, salmin )
      CALL broadcast(ogrid, nstep0 )
      CALL broadcast(ogrid, nstep )
      CALL broadcast(ogrid, time0 )
      CALL broadcast(ogrid, time )

c moved here from inicon:
      if (nstep0.eq.0) then     ! starting from Levitus
        call ssto2a(temp_loc(:,:,1),atmocn%work1)
        call tempro2a(temp_loc(:,:,1),atmocn%work2)
        call ssto2a(saln_loc(:,:,1),atmocn%sss)
c        call ssto2a(omlhc,mlhc)
c
        do ja=aJ_0,aJ_1
          do ia=aI_0,aI_1
            if (atmocn%focean(ia,ja).gt.0.) then
              atmocn%gtemp(ia,ja)=atmocn%work1(ia,ja)
              atmocn%gtempr(ia,ja)=atmocn%work2(ia,ja)
            endif
          enddo
        enddo
c     call findmx(ip,temp,ii,ii,jj,'ini sst')
c     call findmx(ip,saln,ii,ii,jj,'ini sss')
      endif

      do ja=aJ_0,aJ_1
        do ia=aI_0,aI_1
          if (atmocn%focean(ia,ja).gt.0.) then
            if (nstep0.eq.0 .and. atmocn%sss(ia,ja).le.1.) then
              write(*,'(a,2i3,3(a,f6.1))')'chk low saln at agcm ',ia,ja,
     &             ' sss=',atmocn%sss(ia,ja),
     &             ' sst=',atmocn%gtemp(ia,ja),
     &             ' focean=',atmocn%focean(ia,ja)
              stop 'wrong sss in agcm'
            endif
          endif
        enddo
      enddo

c
      END SUBROUTINE init_OCEAN
c
      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
css   ENTRY ODYNAM
      !! fix later: implicit none

      ENTRY ODIFS
      ENTRY io_ocdiag
      ENTRY new_io_ocdiag
      ENTRY def_rsf_ocdiag
      ENTRY def_meta_ocdiag
      ENTRY write_meta_ocdiag
      ENTRY set_ioptrs_ocnacc_default
      ENTRY set_ioptrs_ocnacc_extended
      ENTRY init_ODEEP
      ENTRY reset_ODIAG
      ENTRY diag_OCEAN
      entry OSTRUC(QTCHNG)
      entry OCLIM(end_of_day)
      entry OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFW,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)
css   entry daily_OCEAN(end_of_day)
c      entry PRECIP_OC
      entry GROUND_OC
      entry io_oda(kunit,it,iaction,ioerr)
css   entry io_ocean(iu_GIC,ioread,ioerr)
css   entry CHECKO(SUBR)
c
      ENTRY ADVSI_DIAG
!!      entry alloc_ocean
c --- not calling ice dynamics
css      ENTRY DYNSI
css      ENTRY ADVSI
css      ENTRY io_icedyn
css      ENTRY io_icdiag
css      ENTRY init_icedyn
css      ENTRY reset_icdiag
css      ENTRY diag_ICEDYN
c
      entry diag_OCEAN_prep
      RETURN
      END SUBROUTINE DUMMY_OCN

      subroutine precip_oc(atmocn,iceocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      return
      end subroutine precip_oc

      subroutine diagco(m,atmocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      integer :: m
      type(atmocn_xchng_vars) :: atmocn
      return
      end subroutine diagco

c
!      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!!@sum  io_ocean outputs ocean related fields for restart
!      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, pack_data, unpack_data,
!     &     broadcast, pack_column, unpack_column
!      USE MODEL_COM, only : ioread,iowrite,irsficno,irsfic
!     *     ,irsficnt,irerun,lhead
!!!      USE FLUXES, only : sss,ogeoza,uosurf,vosurf,dmsi,dhsi,dssi
!      USE HYCOM_DIM_GLOB, only : kk,kdm,idm,jdm
!      USE HYCOM_DIM, only : ogrid
!      USE HYCOM_SCALARS, only : nstep,time,oddev,nstep0,time0,baclin
!     &     ,onem,itest,jtest
!#if (defined TRACERS_AGE_OCEAN) \
!     || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!     .  , diag_counter,itest_trac,jtest_trac
!      USE HYCOM_ARRAYS_GLOB_RENAMER, only : plevav_loc,tracav_loc
!#endif
!#ifdef TRACERS_GASEXCH_ocean
!      use domain_decomp_atm, only : agrid=>grid
!      USE TRACER_GASEXCH_COM, only : atrac_loc
!#endif
!
!#ifdef TRACERS_OceanBiology
!      USE obio_forc, only : avgq,tirrq3d,ihra
!      USE obio_com,  only : gcmax
!     .            ,pCO2av,pCO2av_loc,pp2tot_dayav,pp2tot_dayav_loc
!     .            ,ao_co2fluxav,ao_co2fluxav_loc
!     .            ,cexpav,cexpav_loc,diag_counter
!     .            ,pp2tot_day,pp2tot_day_glob
!     .            ,itest_bio=>itest,jtest_bio=>jtest
!      USE obio_com,  only : tracav_loc, plevav_loc, tracav, plevav
!#endif
!      USE HYCOM_ARRAYS_GLOB
!      USE Dictionary_mod
!      IMPLICIT NONE
!c
!      INTEGER, intent(in) :: kunit   !@var kunit unit number of read/write
!      INTEGER, intent(in) :: iaction !@var iaction flag for reading or writing to file
!!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
!      INTEGER, INTENT(INOUT) :: IOERR
!
!c global arrays for i/o
!      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
!     &     SSS,UOSURF,VOSURF,OGEOZA,asst,atempr
!      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: GTEMPR,DMSI,DHSI,DSSI
!      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: GTEMP
!
!!@var HEADER Character string label for individual records
!      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDYN01"
!#if defined(TRACERS_GASEXCH_ocean) || defined(TRACERS_OceanBiology)
!      integer i,j,k
!!@var TRNHEADER Character string label for individual records
!      CHARACTER*80 :: TRNHEADER, TRNMODULE_HEADER = "TRGASEX-OBIOh"
!#ifdef TRACERS_OceanBiology
!      CHARACTER*80 :: TRN2HEADER,
!     .     TRN2MODULE_HEADER = "TRGASXOBIOhdiags"
!#endif 
!#endif
!#if (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!      integer i,j,k
!!@var TRNHEADER Character string label for individual records
!      CHARACTER*80 :: TRNHEADER, TRNMODULE_HEADER = "OCideal trcrs"
!#endif
!#ifdef TRACERS_OceanBiology
!      real, allocatable :: avgq_glob(:,:,:),tirrq3d_glob(:,:,:),
!     &     gcmax_glob(:,:,:),atrac_glob(:,:,:)
!      integer, allocatable :: ihra_glob(:,:)
!#endif
!
!      call sync_param( "itest", itest)
!      call sync_param( "jtest", jtest)
!
!#ifdef TRACERS_OCEAN
!!@var TRHEADER Character string label for individual records
!      CHARACTER*80 :: TRHEADER, TRMODULE_HEADER = "TROCDYN02"
!c
!      write (TRMODULE_HEADER(lhead+1:80),'(a13,i3,a1,i3,a)')
!     *     'R8 dim(im,jm,',LMO,',',NTM,'):TRMO,TX,TY,TZ'
!#endif
!
!#ifdef TRACERS_OceanBiology
!      if (AM_I_ROOT()) then
!        allocate( avgq_glob(idm,jdm,kdm),tirrq3d_glob(idm,jdm,kdm),
!     &       ihra_glob(idm,jdm), gcmax_glob(idm,jdm,kdm) )
!        allocate(atrac_glob(agrid%im_world,agrid%jm_world,
!     &       size(atrac_loc,3)))
!      endif
!      call pack_data(ogrid, avgq,    avgq_glob)
!      call pack_data(ogrid, tirrq3d, tirrq3d_glob)
!      call pack_data(ogrid, ihra,    ihra_glob)
!      call pack_data(ogrid, gcmax,   gcmax_glob)
!      call pack_data(agrid, atrac_loc, atrac_glob)
!#endif
!
!      ! move to global atm grid
!      call alloc_atm_globals
!      call gather_atm_before_checkpoint
!      call gather_hycom_arrays   !mkb Jun  6
!
!#if (defined TRACERS_OceanBiology) || defined (TRACERS_GASEXCH_ocean) \
!     || (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)   
!      call pack_data(ogrid, tracav_loc, tracav)
!      call pack_data(ogrid, plevav_loc, plevav)
!#endif
!#ifdef TRACERS_OceanBiology
!      call pack_data(ogrid, pCO2av_loc, pCO2av)
!      call pack_data(ogrid, pp2tot_dayav_loc, pp2tot_dayav)      !time integrated pp2tot_day
!      call pack_data(ogrid, ao_co2fluxav_loc,ao_co2fluxav)
!      call pack_data(ogrid, cexpav_loc, cexpav)
!      call pack_data(ogrid, pp2tot_day, pp2tot_day_glob)      !instantaneous pp2tot_day
!#endif
!
!      if (AM_I_ROOT()) then ! work on global grids here
!
!c
!css   write (MODULE_HEADER(lhead+1:80),'(a13,i2,a)') 'R8 dim(im,jm,',
!css  *   LMO,'):M,U,V,G0,GX,GY,GZ,S0,SX,SY,SZ, OGZ,OGZSV'
!c
!      write(*,'(a,i9,f9.0)')'chk ocean write at nstep/day=',nstep,time
!      write (MODULE_HEADER(lhead+1:80),'(a,i8,f8.1,a)')
!     . 'u,v,dp,t,s,th,tb,ub,vb,pb,pb,psi,thk,mxl,uf,vf,df,tcr3+o18+a8'
!
!#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_OceanBiology)
!      write(*,'(a,i9,f9.0)')'chk GASEXCH write at nstep/day=',nstep,time
!      write (TRNMODULE_HEADER(lhead+1:80),'(a29)')
!     *'atrac,avgq,gcmax,tirrq,ihra,tracav,pCO2av,co2flxav,cexpav,diag_c'
!      write (TRN2MODULE_HEADER(lhead+1:80),'(a29)')
!     *     'pp2tot_day,pp2tot_dayav'
!#else
!#ifdef TRACERS_GASEXCH_ocean
!      write(*,'(a,i9,f9.0)')'chk GASEXCH write at nstep/day=',nstep,time
!      write (TRNMODULE_HEADER(lhead+1:80),'(a5)')
!     *     'atrac'    
!#endif
!#ifdef TRACERS_OceanBiology
!      write(*,'(a,i9,f9.0)')'chk OCN BIO write at nstep/day=',nstep,time
!      write (TRNMODULE_HEADER(lhead+1:80),'(a63)')
!     *'avgq,gcmax,tirrq3d,ihra,tracav,pCO2av,ao_co2fluxav,cexpav,diag_counter'
!      write (TRN2MODULE_HEADER(lhead+1:80),'(a29)')
!     *     'pp2tot_day,pp2tot_dayav'
!#endif
!#endif
!
!#if (defined TRACERS_AGE_OCEAN) \
!     || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!      write(*,'(a,i9,f9.0)')'chk TRACERS write at nstep/day=',nstep,time
!      write (TRNMODULE_HEADER(lhead+1:80),'(a63)')
!     *'tracav,plevav,diag_counter'
!#endif
!
!
!      SELECT CASE (IACTION)
!c---------------------------------------------------------------------------------
!      CASE (:IOWRITE)            ! output to standard restart file
!css     WRITE (kunit,err=10) MODULE_HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO
!css  *     ,S0M,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
!css#ifdef TRACERS_OCEAN
!css       WRITE (kunit,err=10) TRMODULE_HEADER,tracer
!css#endif
!        WRITE (kunit,err=10) MODULE_HEADER,nstep,time
!     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
!     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
!     . ,dpav,temav,salav,th3av,ubavav,vbavav,pbavav,sfhtav,eminpav
!     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
!     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
!
!#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_OceanBiology)
!      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
!     . ,atrac_glob,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_counter
!      WRITE (kunit,err=10) TRN2MODULE_HEADER,nstep,time
!     . ,pp2tot_day_glob,pp2tot_dayav
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,2(e12.4,1x))') ' tst1a k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j),atrac_glob(10,20,1),pp2tot_day_glob(i,j)
!      enddo
!#else
!#ifdef TRACERS_GASEXCH_ocean
!      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
!     . ,atrac_glob
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,e12.4)') ' tst1a k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j),atrac_glob(10,20,1)
!      enddo
!#endif
!#ifdef TRACERS_OceanBiology
!      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
!     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_countere
!      WRITE (kunit,err=10) TRN2MODULE_HEADER,nstep,time
!     . ,pp2tot_day_glob,pp2tot_dayav
!      i=itest_bio
!      j=jtest_bio
!      print*,'test point at:',itest_bio,jtest_bio
!
!      do k=1,kdm
!      write(*,'(a,i2,8(e12.4,1x),i3)') ' tst1a k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     .    pp2tot_day_glob(i,j),
!     .    ihra_glob(i,j)
!      enddo
!#endif
!#endif
!
!#if (defined TRACERS_AGE_OCEAN) || defined(TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!      WRITE (kunit,err=10) TRNMODULE_HEADER,nstep,time
!     . ,tracav,plevav,diag_counter
!      i=itest_trac
!      j=jtest_trac
!      print*,'test point at:',itest_trac,jtest_trac
!
!      do k=1,kdm
!      write(*,'(a,i2,6(e12.4,1x))') ' tst1a k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),tracer(i,j,k,1),
!     .    plevav(i,j,k),tracav(i,j,k,1),diag_counter
!      enddo
!#endif
!
!c---------------------------------------------------------------------------------
!      CASE (IOREAD:)            ! input from restart file
!        SELECT CASE (IACTION)
!c       --------------------------------------------------------------------------
!          CASE (IRSFICNO)   ! initial conditions (no ocean data)
!            READ (kunit)
!c       --------------------------------------------------------------------------
!          CASE (ioread,irerun,irsfic) ! restarts
!css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
!css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
!c
!            !!call geopar
!            READ (kunit,err=10) HEADER,nstep0,time0
!     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
!     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
!     . ,dpav,temav,salav,th3av,ubavav,vbavav,pbavav,sfhtav,eminpav
!     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
!     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
!
!      nstep0=time0*86400./baclin+.0001
!      write(*,'(a,i9,f9.0)')'chk ocean read at nstep/day=',nstep0,time0
!      nstep=nstep0
!      time=time0
!c
!            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",HEADER
!     *             ,MODULE_HEADER
!              GO TO 10
!            END IF
!
!#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_OceanBiology)
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,atrac_glob,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_counter
!      READ (kunit,err=10) TRN2HEADER,nstep0,time0
!     . ,pp2tot_day_glob,pp2tot_dayav
!      write(*,'(a,i9,f9.0)')'chk GASEXCH read at nstep/day=',nstep0,time0
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,2(e12.4,1x))') ' tst1b k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &    ihra_glob(i,j),atrac_glob(10,20,1),
!     .    pp2tot_day_glob(i,j)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!            IF (TRN2HEADER(1:LHEAD).NE.TRN2MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRN2HEADER
!     .             ,TRN2MODULE_HEADER
!              GO TO 10
!            END IF
!#else
!#ifdef TRACERS_GASEXCH_ocean
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,atrac_glob
!      i=itest_bio
!      j=jtest
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,e12.4)') ' tst1b k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j),atrac_glob(10,20,1)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!#endif
!#ifdef TRACERS_OceanBiology
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_counter
!      READ (kunit,err=10) TRN2HEADER,nstep0,time0
!     . ,pp2tot_day_glob,pp2tot_dayav
!      i=itest_bio
!      j=jtest_bio
!      print*, 'itest, jtest=',itest_bio,jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,8(e12.4,1x),i3)') ' tst1b k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     .    pp2tot_day_glob(i,j),
!     .    ihra_glob(i,j)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!            IF (TRN2HEADER(1:LHEAD).NE.TRN2MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRN2HEADER
!     .             ,TRN2MODULE_HEADER
!              GO TO 10
!            END IF
!#endif
!#endif
!
!#if (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,tracav,plevav,diag_counter
!      i=itest_trac
!      j=jtest_trac
!      print*,'test point at:',itest_trac,jtest_trac
!
!      do k=1,kdm
!      write(*,'(a,i2,6(e12.4,1x))') ' tst1b k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),tracer(i,j,k,1),
!     .    plevav(i,j,k),tracav(i,j,k,1)
!      enddo
!#endif
!
!#ifdef TRACERS_OCEAN
!            READ (kunit,err=10) TRHEADER,TRMO,TXMO,TYMO,TZMO
!            IF (TRHEADER(1:LHEAD).NE.TRMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRHEADER
!     *             ,TRMODULE_HEADER
!              GO TO 10
!            END IF
!#endif
!
!c       --------------------------------------------------------------------------
!          CASE (irsficnt) ! restarts (never any tracer data)
!css         READ (kunit,err=10) HEADER,MO,UO,VO,G0M,GXMO,GYMO,GZMO,S0M
!css  *           ,SXMO,SYMO,SZMO,OGEOZ,OGEOZ_SV
!c
!            print*,'restarts (never any tracer data -irsficnt)'
!            !!call geopar
!            READ (kunit,err=10) HEADER,nstep0,time0
!     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
!     . ,uflxav,vflxav,diaflx
!#if (defined TRACERS_OceanBiology) \
!     || (defined TRACERS_AGE_OCEAN) \
!     || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!     . ,tracer(:,:,:,1)
!#else
!!Shan Sun's rsf files have one dimensional tracer
!     . ,tracer
!#endif
!     . ,dpinit,oddev,uav,vav,dpuav,dpvav
!     . ,dpav,temav,salav,th3av,ubavav,vbavav,pbavav,sfhtav,eminpav
!     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
!     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi         ! agcm grid
!
!      nstep0=time0*86400./baclin+.0001
!      write(*,'(a,i9,f9.0)')'chk ocean read at nstep/day=',nstep0,time0
!
!            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",HEADER
!     *             ,MODULE_HEADER
!              GO TO 10
!            END IF
!
!      if (nstep.eq.0)go to 222   !for a cold start the AIC file does not have this stuff.
!#if defined(TRACERS_GASEXCH_ocean) && defined(TRACERS_OceanBiology)
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,atrac_glob,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_counter
!      READ (kunit,err=10) TRN2HEADER,nstep0,time0
!     . ,pp2tot_day_glob,pp2tot_dayav
!      write(*,'(a,i9,f9.0)')
!     &     'chk GASEXCH read at nstep/day=',nstep0,time0
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,e12.4)') ' tst2 k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j),atrac_glob(10,20,1)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!            IF (TRN2HEADER(1:LHEAD).NE.TRN2MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRN2HEADER
!     .             ,TRN2MODULE_HEADER
!              GO TO 10
!            END IF
!#else
!#ifdef TRACERS_GASEXCH_ocean
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,atrac_glob
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3,1x,e12.4)') ' tst2 k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j),atrac_glob(10,20,1)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!#endif
!#ifdef TRACERS_OceanBiology
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,avgq_glob,gcmax_glob,tirrq3d_glob,ihra_glob
!     . ,tracav,plevav,pCO2av,ao_co2fluxav,cexpav,diag_counter
!      READ (kunit,err=10) TRN2HEADER,nstep0,time0
!     . ,pp2tot_day_glob,pp2tot_dayav
!      i=itest_bio
!      j=jtest_bio
!      do k=1,kdm
!      write(*,'(a,i2,7(e12.4,1x),i3)') ' tst2 k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),avgq_glob(i,j,k),gcmax_glob(i,j,k),
!     .    tracer(i,j,k,1),tracer(i,j,k,15),tirrq3d_glob(i,j,k),
!     &       ihra_glob(i,j)
!      enddo
!            IF (TRNHEADER(1:LHEAD).NE.TRNMODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRNHEADER
!     .             ,TRNMODULE_HEADER
!              GO TO 10
!            END IF
!            IF (TRN2HEADER(1:LHEAD).NE.TRN2MODULE_HEADER(1:LHEAD)) THEN
!              PRINT*,"Discrepancy in module version ",TRN2HEADER
!     .             ,TRN2MODULE_HEADER
!              GO TO 10
!            END IF
!#endif
!#endif
!
!#if (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)
!      READ (kunit,err=10) TRNHEADER,nstep0,time0
!     . ,tracav,plevav,diag_counter
!      i=itest_trac
!      j=jtest_trac
!      print*,'test point at:',itest_trac,jtest_trac
!
!      do k=1,kdm
!      write(*,'(a,i2,6(e12.4,1x))') ' tst2 k=',k,
!     .    dp(i,j,k)/onem,temp(i,j,k),tracer(i,j,k,1),
!     .    plevav(i,j,k),tracav(i,j,k,1)
!      enddo
!#endif
!
! 222  continue
!
!          END SELECT
!      END SELECT
!
!      endif ! AM_I_ROOT
!      call scatter_atm_after_checkpoint
!      call dealloc_atm_globals
!      CALL broadcast(ogrid, nstep0 )
!      CALL broadcast(ogrid, time0 )
!
!#ifdef TRACERS_OceanBiology
!      call unpack_data(ogrid, avgq_glob, avgq)
!      call unpack_data(ogrid, tirrq3d_glob, tirrq3d)
!      call unpack_data(ogrid, ihra_glob, ihra)
!      call unpack_data(ogrid, gcmax_glob, gcmax)
!      call unpack_data(agrid, atrac_glob, atrac_loc)
!      if (AM_I_ROOT()) then
!        deallocate( avgq_glob,tirrq3d_glob,
!     &       ihra_glob, gcmax_glob, atrac_glob )
!      endif
!#endif
!
!#if (defined TRACERS_OceanBiology) || defined (TRACERS_GASEXCH_ocean) \
!      || (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
!     || (defined TRACERS_ZEBRA)   
!      call unpack_data(ogrid, tracav, tracav_loc)
!      call unpack_data(ogrid, plevav, plevav_loc)
!#endif
!#ifdef TRACERS_OceanBiology
!      call unpack_data(ogrid, pCO2av, pCO2av_loc)
!      call unpack_data(ogrid, pp2tot_dayav, pp2tot_dayav_loc)
!      call unpack_data(ogrid, ao_co2fluxav, ao_co2fluxav_loc)
!      call unpack_data(ogrid, cexpav, cexpav_loc)
!      call unpack_data(ogrid, pp2tot_day_glob, pp2tot_day)
!#endif
!
!      RETURN
! 10   IOERR=1
!      call scatter_atm_after_checkpoint
!      call dealloc_atm_globals
!      ! why do we need return after error?
!      call stop_model("error i/o in io_ocean",255)
!      RETURN
!C****
!      contains
!      subroutine alloc_atm_globals
!      USE RESOLUTION, only : im,jm
!      use FLUXES, only: NSTYPE
!      if(am_i_root()) then
!        ALLOCATE( SSS( im, jm ) )
!        ALLOCATE( UOSURF( im, jm ) )
!        ALLOCATE( VOSURF( im, jm ) )
!        ALLOCATE( OGEOZA( im, jm ) )
!        ALLOCATE( GTEMP( 2 , NSTYPE, im, jm ) )
!        ALLOCATE( GTEMPR( NSTYPE, im, jm ) )
!        ALLOCATE( DMSI(  2  , im, jm ) )
!        ALLOCATE( DHSI(  2  , im, jm ) )
!        ALLOCATE( DSSI(  2  , im, jm ) )
!        ALLOCATE( asst( im, jm ) )
!        ALLOCATE( atempr( im, jm ) )
!      endif
!      end subroutine alloc_atm_globals
!      subroutine dealloc_atm_globals
!      if(am_i_root()) then
!        DEALLOCATE(SSS,UOSURF,VOSURF,
!     &       OGEOZA,GTEMP,GTEMPR,DMSI,DHSI,DSSI,asst,atempr)
!      endif
!      end subroutine dealloc_atm_globals
!      subroutine gather_atm_before_checkpoint
!      USE DOMAIN_DECOMP_ATM, ONLY: GRID
!      use hycom_atm
!      call pack_data( grid,  ASST_loc, ASST )
!      call pack_data( grid,  ATEMPR_loc, ATEMPR )
!      call pack_data( grid,  SSS_loc, SSS )
!      call pack_data( grid,  UOSURF_loc, UOSURF )
!      call pack_data( grid,  VOSURF_loc, VOSURF )
!      call pack_data( grid,  OGEOZA_loc, OGEOZA )
!      call pack_column( grid,  DMSI_loc, DMSI )
!      call pack_column( grid,  DHSI_loc, DHSI )
!      call pack_column( grid,  DSSI_loc, DSSI )
!      end subroutine gather_atm_before_checkpoint
!
!      subroutine scatter_atm_after_checkpoint
!      USE DOMAIN_DECOMP_ATM, ONLY: GRID
!      use hycom_atm
!      call unpack_data( grid,  ASST, ASST_loc )
!      call unpack_data( grid,  ATEMPR, ATEMPR_loc )
!      call unpack_data( grid,  SSS, SSS_loc )
!      call unpack_data( grid,  UOSURF, UOSURF_loc )
!      call unpack_data( grid,  VOSURF, VOSURF_loc )
!c UOSURF and VOSURF are also needed on the ice dynamics A-grid.
!c For the moment, HYCOM only runs with modelE configurations having
!c identical atmosphere and ice dynamics grids, so the atmospheric
!c copy of UOSURF,VOSURF can be used.
!      if(grid_icdyn%have_domain) then ! ice dyn may run on subset of PEs
!        call unpack_data( grid_icdyn,  UOSURF, UOSURF_4DYNSI_loc)
!        call unpack_data( grid_icdyn,  VOSURF, VOSURF_4DYNSI_loc) 
!      endif
!      call unpack_data( grid,  OGEOZA, OGEOZA_loc )
!      call unpack_column( grid,  DMSI, DMSI_loc )
!      call unpack_column( grid,  DHSI, DHSI_loc )
!      call unpack_column( grid,  DSSI, DSSI_loc )
!      end subroutine scatter_atm_after_checkpoint
!
!      END SUBROUTINE io_ocean

#ifdef NEW_IO
      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE HYCOM_DIM, only : grid=>ogrid
      use pario, only : defvar
      USE HYCOM_SCALARS, only : nstep,time,oddev
      USE HYCOM_ARRAYS
      implicit none
      integer fid   !@var fid file id
      integer :: n
      character(len=14) :: str2d
      character(len=18) :: str3d
      character(len=20) :: str3d2
      str2d ='(idm,dist_jdm)'
      str3d ='(idm,dist_jdm,kdm)'
      str3d2='(idm,dist_jdm,kdmx2)'

#if defined(TRACERS_OceanBiology)
      call def_rsf_obio(fid)
#endif

      call defvar(grid,fid,nstep,'nstep')
      call defvar(grid,fid,time,'time')

      call defvar(grid,fid,u,'uo'//str3d2)
      call defvar(grid,fid,v,'vo'//str3d2)
      call defvar(grid,fid,dp,'dp'//str3d2)
      call defvar(grid,fid,temp,'temp'//str3d2)
      call defvar(grid,fid,saln,'saln'//str3d2)
      call defvar(grid,fid,th3d,'th3d'//str3d2)
      call defvar(grid,fid,ubavg,'ubavg(idm,dist_jdm,three)')
      call defvar(grid,fid,vbavg,'vbavg(idm,dist_jdm,three)')
      call defvar(grid,fid,pbavg,'pbavg(idm,dist_jdm,three)')
      call defvar(grid,fid,pbot,'pbot'//str2d)
      call defvar(grid,fid,psikk,'psikk'//str2d)
      call defvar(grid,fid,thkk,'thkk'//str2d)
      call defvar(grid,fid,dpmixl,'dpmixl(idm,dist_jdm,two)')
      call defvar(grid,fid,uflxav,'uflxav'//str3d)
      call defvar(grid,fid,vflxav,'vflxav'//str3d)
      call defvar(grid,fid,diaflx,'diaflx'//str3d)
      call defvar(grid,fid,tracer,'tracer(idm,dist_jdm,kdm,ntrcr)')
      call defvar(grid,fid,dpinit,'dpinit'//str3d)
      call defvar(grid,fid,oddev,'oddev')
      call defvar(grid,fid,uav,'uav'//str3d)
      call defvar(grid,fid,vav,'vav'//str3d)
      call defvar(grid,fid,dpuav,'dpuav'//str3d)
      call defvar(grid,fid,dpvav,'dpvav'//str3d)
      call defvar(grid,fid,dpav,'dpav'//str3d)
      call defvar(grid,fid,temav,'temav'//str3d)
      call defvar(grid,fid,salav,'salav'//str3d)
      call defvar(grid,fid,th3av,'th3av'//str3d)
      call defvar(grid,fid,ubavav,'ubavav'//str2d)
      call defvar(grid,fid,vbavav,'vbavav'//str2d)
      call defvar(grid,fid,pbavav,'pbavav'//str2d)
      call defvar(grid,fid,sfhtav,'sfhtav'//str2d)
      call defvar(grid,fid,eminpav,'eminpav'//str2d)
      call defvar(grid,fid,surflav,'surflav'//str2d)
      call defvar(grid,fid,salflav,'salflav'//str2d)
      call defvar(grid,fid,brineav,'brineav'//str2d)
      call defvar(grid,fid,tauxav,'tauxav'//str2d)
      call defvar(grid,fid,tauyav,'tauyav'//str2d)
      call defvar(grid,fid,dpmxav,'dpmxav'//str2d)
      call defvar(grid,fid,oiceav,'oiceav'//str2d)

c write:
c        WRITE (kunit,err=10) nstep,time
c     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
c     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
c     . ,dpav,temav,salav,th3av,ubavav,vbavav,pbavav,sfhtav,eminpav
c     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
c     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi  ! agcm grid

c read: note it reads in nstep0,time0 instead of nstep,time
c            READ (kunit,err=10) HEADER,nstep0,time0
c     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
c     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
c     . ,dpav,temav,salav,th3av,ubavav,vbavav,pbavav,sfhtav,eminpav
c     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
c     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi  ! agcm grid

      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      USE HYCOM_DIM, only : grid=>ogrid
      USE HYCOM_SCALARS, only : nstep,time,nstep0,time0,baclin,oddev
      USE HYCOM_ARRAYS
      USE HYCOM_ARRAYS_GLOB, only : gather_hycom_arrays   ! for now
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
#if defined(TRACERS_OceanBiology)
        call new_io_obio(fid,iaction)
#endif
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'nstep',nstep)
        call write_data(grid,fid,'time',time)
        call write_dist_data(grid,fid,'uo',u)
        call write_dist_data(grid,fid,'vo',v)
        call write_dist_data(grid,fid,'dp',dp)
        call write_dist_data(grid,fid,'temp',temp)
        call write_dist_data(grid,fid,'saln',saln)
        call write_dist_data(grid,fid,'th3d',th3d)
        call write_dist_data(grid,fid,'ubavg',ubavg)
        call write_dist_data(grid,fid,'vbavg',vbavg)
        call write_dist_data(grid,fid,'pbavg',pbavg)
        call write_dist_data(grid,fid,'pbot',pbot)
        call write_dist_data(grid,fid,'psikk',psikk)
        call write_dist_data(grid,fid,'thkk',thkk)
        call write_dist_data(grid,fid,'dpmixl',dpmixl)
        call write_dist_data(grid,fid,'uflxav',uflxav)
        call write_dist_data(grid,fid,'vflxav',vflxav)
        call write_dist_data(grid,fid,'diaflx',diaflx)
        call write_dist_data(grid,fid,'tracer',tracer)
        call write_dist_data(grid,fid,'dpinit',dpinit)
        call write_data(grid,fid,'oddev',oddev)
        call write_dist_data(grid,fid,'uav',uav)
        call write_dist_data(grid,fid,'vav',vav)
        call write_dist_data(grid,fid,'dpuav',dpuav)
        call write_dist_data(grid,fid,'dpvav',dpvav)
        call write_dist_data(grid,fid,'dpav',dpav)
        call write_dist_data(grid,fid,'temav',temav)
        call write_dist_data(grid,fid,'salav',salav)
        call write_dist_data(grid,fid,'th3av',th3av)
        call write_dist_data(grid,fid,'ubavav',ubavav)
        call write_dist_data(grid,fid,'vbavav',vbavav)
        call write_dist_data(grid,fid,'pbavav',pbavav)
        call write_dist_data(grid,fid,'sfhtav',sfhtav)
        call write_dist_data(grid,fid,'eminpav',eminpav)
        call write_dist_data(grid,fid,'surflav',surflav)
        call write_dist_data(grid,fid,'salflav',salflav)
        call write_dist_data(grid,fid,'brineav',brineav)
        call write_dist_data(grid,fid,'tauxav',tauxav)
        call write_dist_data(grid,fid,'tauyav',tauyav)
        call write_dist_data(grid,fid,'dpmxav',dpmxav)
        call write_dist_data(grid,fid,'oiceav',oiceav)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'nstep',nstep0,bcast_all=.true.)
        call read_data(grid,fid,'time',time0,bcast_all=.true.)
        nstep0=time0*SECONDS_PER_DAY/baclin+.0001
        write(*,'(a,i9,f9.0)')
     &       'chk ocean read at nstep/day=',nstep0,time0
        nstep=nstep0
        time=time0
        call read_dist_data(grid,fid,'uo',u)
        call read_dist_data(grid,fid,'vo',v)
        call read_dist_data(grid,fid,'dp',dp)
        call read_dist_data(grid,fid,'temp',temp)
        call read_dist_data(grid,fid,'saln',saln)
        call read_dist_data(grid,fid,'th3d',th3d)
        call read_dist_data(grid,fid,'ubavg',ubavg)
        call read_dist_data(grid,fid,'vbavg',vbavg)
        call read_dist_data(grid,fid,'pbavg',pbavg)
        call read_dist_data(grid,fid,'pbot',pbot)
        call read_dist_data(grid,fid,'psikk',psikk)
        call read_dist_data(grid,fid,'thkk',thkk)
        call read_dist_data(grid,fid,'dpmixl',dpmixl)
        call read_dist_data(grid,fid,'uflxav',uflxav)
        call read_dist_data(grid,fid,'vflxav',vflxav)
        call read_dist_data(grid,fid,'diaflx',diaflx)
        call read_dist_data(grid,fid,'tracer',tracer)
        call read_dist_data(grid,fid,'dpinit',dpinit)
        call read_data(grid,fid,'oddev',oddev,bcast_all=.true.)
        call read_dist_data(grid,fid,'uav',uav)
        call read_dist_data(grid,fid,'vav',vav)
        call read_dist_data(grid,fid,'dpuav',dpuav)
        call read_dist_data(grid,fid,'dpvav',dpvav)
        call read_dist_data(grid,fid,'dpav',dpav)
        call read_dist_data(grid,fid,'temav',temav)
        call read_dist_data(grid,fid,'salav',salav)
        call read_dist_data(grid,fid,'th3av',th3av)
        call read_dist_data(grid,fid,'ubavav',ubavav)
        call read_dist_data(grid,fid,'vbavav',vbavav)
        call read_dist_data(grid,fid,'pbavav',pbavav)
        call read_dist_data(grid,fid,'sfhtav',sfhtav)
        call read_dist_data(grid,fid,'eminpav',eminpav)
        call read_dist_data(grid,fid,'surflav',surflav)
        call read_dist_data(grid,fid,'salflav',salflav)
        call read_dist_data(grid,fid,'brineav',brineav)
        call read_dist_data(grid,fid,'tauxav',tauxav)
        call read_dist_data(grid,fid,'tauyav',tauyav)
        call read_dist_data(grid,fid,'dpmxav',dpmxav)
        call read_dist_data(grid,fid,'oiceav',oiceav)
c certain initialization routines still work with global
c arrays, so we have to gather
        call gather_checkpointed_hycom_arrays
      end select

      return
      end subroutine new_io_ocean

      subroutine gather_checkpointed_hycom_arrays
      ! TODO: see which global-domain arrays are
      ! really needed in non-parallelized init routines
      use hycom_arrays_glob
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      call pack_data( ogrid,  u_loc, u )
      call pack_data( ogrid,  v_loc, v )
      call pack_data( ogrid,  dp_loc, dp )
      call pack_data( ogrid,  temp_loc, temp )
      call pack_data( ogrid,  saln_loc, saln )
      call pack_data( ogrid,  th3d_loc, th3d )
      call pack_data( ogrid,  ubavg_loc, ubavg )
      call pack_data( ogrid,  vbavg_loc, vbavg )
      call pack_data( ogrid,  pbavg_loc, pbavg )
      call pack_data( ogrid,  pbot_loc, pbot )
      call pack_data( ogrid,  psikk_loc, psikk )
      call pack_data( ogrid,  thkk_loc, thkk )
      call pack_data( ogrid,  dpmixl_loc, dpmixl )
      call pack_data( ogrid,  uflxav_loc, uflxav )
      call pack_data( ogrid,  vflxav_loc, vflxav )
      call pack_data( ogrid,  diaflx_loc, diaflx )
      call pack_data( ogrid,  tracer_loc, tracer )
      call pack_data( ogrid,  dpinit_loc, dpinit )
      call pack_data( ogrid,  uav_loc, uav )
      call pack_data( ogrid,  vav_loc, vav )
      call pack_data( ogrid,  dpuav_loc, dpuav )
      call pack_data( ogrid,  dpvav_loc, dpvav )
      call pack_data( ogrid,  dpav_loc, dpav )
      call pack_data( ogrid,  temav_loc, temav )
      call pack_data( ogrid,  salav_loc, salav )
      call pack_data( ogrid,  th3av_loc, th3av )
      call pack_data( ogrid,  ubavav_loc, ubavav )
      call pack_data( ogrid,  vbavav_loc, vbavav )
      call pack_data( ogrid,  pbavav_loc, pbavav )
      call pack_data( ogrid,  sfhtav_loc, sfhtav )
      call pack_data( ogrid,  eminpav_loc, eminpav )
      call pack_data( ogrid,  surflav_loc, surflav )
      call pack_data( ogrid,  salflav_loc, salflav )
      call pack_data( ogrid,  brineav_loc, brineav )
      call pack_data( ogrid,  tauxav_loc, tauxav )
      call pack_data( ogrid,  tauyav_loc, tauyav )
      call pack_data( ogrid,  dpmxav_loc, dpmxav )
      call pack_data( ogrid,  oiceav_loc, oiceav )
      return
      end subroutine gather_checkpointed_hycom_arrays

#endif /* NEW_IO */

c
      SUBROUTINE CHECKO(SUBR)
#ifdef USE_ATM_GLOBAL_ARRAYS
!@sum  CHECKO Checks whether Ocean are reasonable
!!      USE MODEL_COM, only : im,jm
!!      USE FLUXES, only : gtemp
!!      USE MODEL_COM, only : focean
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only: pack_block,AM_I_ROOT
c      USE HYCOM_ATM, only : gtemp,gtemp_loc
      IMPLICIT NONE
      integer i,j

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

c      call pack_block( grid,  GTEMP_loc, GTEMP )
c      if (AM_I_ROOT()) then

c      print *,'SUBR=',SUBR
c      write(*,'(10f7.2)') ((gtemp(1,1,i,j),i=1,10),j=15,20)
c      write(*,'(a)') 'focean'
c      write(*,'(10f7.2)') ((focean(i,j),i=1,10),j=15,20)

c      endif ! AM_I_ROOT
      ! no need to sctter since nothing changed
#endif /* USE_ATM_GLOBAL_ARRAYS */
      END SUBROUTINE CHECKO
c

      SUBROUTINE daily_OCEAN(end_of_day,atmocn)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
C****
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
      type(atmocn_xchng_vars) :: atmocn
      RETURN
      END SUBROUTINE daily_OCEAN
c
css   REAL*8 FUNCTION TFREZS (SIN)
C****
C**** TFREZS calculates the freezing temperature of sea water as a
C**** function of salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: SIN (1) = salinity (kg NaCl/kg sea water), from .004 to .04
C****
C**** Output: TFREZS (C) = freezing temperature of sea water
C****
css   IMPLICIT NONE
css    REAL*8, INTENT(IN) :: SIN
css      REAL*8 :: A01 = -.0575d0, A02 = -2.154996D-4, A03 =1.710523D-3
css      REAL*8 S,S32
C****
css      S   = SIN*1.D3
css      S32 = S*DSQRT(S)
css      TFREZS = (A01 + A02*S)*S + A03*S32
css      RETURN
css      END
c
      subroutine gather_odiags
C     nothing to gather - ocean prescribed
      implicit none
      return
      end subroutine gather_odiags


      subroutine alloc_ocean
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE HYCOM_DIM, only : init_hycom_grid, alloc_hycom_dim
      USE HYCOM_ARRAYS, only : alloc_hycom_arrays

      USE HYCOM_DIM_GLOB, only : alloc_hycom_dim_glob
      USE HYCOM_ARRAYS_GLOB, only : alloc_hycom_arrays_glob

      USE KPRF_ARRAYS, only : alloc_kprf_arrays, alloc_kprf_arrays_local
      USE HYCOM_ATM, only : alloc_hycom_atm

      implicit none
      

      ! seems like this is ok place to create ocean grid since nobody
      ! uses it before this call...
      call init_hycom_grid

      

      !call alloc_hycom_atm

      call alloc_hycom_dim
      call alloc_hycom_arrays

      call alloc_hycom_dim_glob
      call alloc_hycom_arrays_glob

      call alloc_kprf_arrays
      call alloc_kprf_arrays_local

#ifdef TRACERS_OceanBiology
      call alloc_obio_com
#endif

      call geopar(.true.)

      !!call reset_hycom_arrays

      !if (AM_I_ROOT()) then
 !!!       call geopar
      !endif


      end subroutine alloc_ocean

#ifdef THIS_PART_IS_NOT_READY
      subroutine reset_hycom_arrays
      USE HYCOM_DIM_GLOB, only : ii,jj,kk
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS, only : sswflx
      USE HYCOM_SCALARS, only : lp,huge
      implicit none

      integer i,j,k,ja,jb,ia
      real :: zero = 0.
      

      write (*,*) 'laying out arrays in memory ...'
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,:)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,:)=huge
      vbavg(i,j,:)=huge
       utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=zero
      dpmixl(i,j,:)=zero
      oice(i,j)=zero
      taux(i,j)=zero
      tauy(i,j)=zero
      oflxa2o(i,j)=zero
      osalt(i,j)=zero
      oemnp(i,j)=zero
      ustar(i,j)=zero
      sswflx(i,j)=zero
c
      ubavav(i,j)=zero
      vbavav(i,j)=zero
      pbavav(i,j)=zero
      sfhtav(i,j)=zero
      dpmxav(i,j)=zero
      oiceav(i,j)=zero
      eminpav(i,j)=zero
      surflav(i,j)=zero
      tauxav(i,j)=zero
      tauyav(i,j)=zero
      salflav(i,j)=zero
      brineav(i,j)=zero
c
      u  (i,j,:   )=huge
      v  (i,j,:   )=huge
      uflx(i,j,:)=huge
      vflx(i,j,:)=huge
      ufxcum(i,j,:)=huge
      vfxcum(i,j,:)=huge
      dpinit(i,j,:)=huge
      dpold (i,j,:)=huge
      dp (i,j,:   )=huge
      dpu(i,j,k   )=huge
      dpv(i,j,k   )=huge
      p (i,j,:)=huge
      pu(i,j,:)=huge
      pv(i,j,:)=huge
c
      th3d(i,j,:)=huge
      thstar(i,j,:)=huge
!      do nt=1,ntrcr
        tracer(i,j,:,:)=zero
!      end do
      uav(i,j,:)=zero
      vav(i,j,:)=zero
      dpuav(i,j,:)=zero
      dpvav(i,j,:)=zero
      dpav (i,j,:)=zero
      temav(i,j,:)=zero
      salav(i,j,:)=zero
      th3av(i,j,:)=zero
      uflxav(i,j,:)=zero
      vflxav(i,j,:)=zero
      diaflx(i,j,:)=zero
 209  continue
c
      do 210 j=1,jj
      !!ja=mod(j-2+jj,jj)+1
      !!do 210 l=1,isq(j)
      !!do 210 i=ifq(j,l),ilq(j,l)
      do i=1,ii
      pbot(i  ,j  )=0.
      !pbot(i-1,j  )=0.
      !pbot(i  ,ja )=0.
      !pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      !p(i-1,j  ,1)=0.
      !p(i  ,ja ,1)=0.
      !p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,:   )=0.
      dp(i  ,j  ,k+kk)=0.
      !dp(i-1,j  ,k   )=0.
      !dp(i-1,j  ,k+kk)=0.
      !dp(i  ,ja ,k   )=0.
      !dp(i  ,ja ,k+kk)=0.
      !dp(i-1,ja ,k   )=0.
 !210  !dp(i-1,ja ,k+kk)=0.
 210  continue
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
      do 156 j=1,jj
      do 156 i=1,ii !ifu(j,l),ilu(j,l)
      pu(i,j,:)=0.
c
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      dpu(i,j,:   )=0.
c
 156  continue
c
      do 158 j=1,jj
      !do 158 l=1,isp(j)
      do 158 i=1,ii !ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      dpu(:,:,:)=0.
      uflx(:,:,:)=0.
      ufxcum(:,:,:)=0.
      u(:,:,:)=0.
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      pv(:,:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      dpv(:,:,:   )=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
      vbavg(:,:,:)=0.
      vtotm (:,:)=0.
      vflux (:,:)=0.
      vflux2(:,:)=0.
      vflux3(:,:)=0.
      via(:,:)=0.
      vib(:,:)=0.
c
      dpv(:,:,:   )=0.
      vflx(:,:,:)=0.
      vfxcum(:,:,:)=0.
      v(:,:,:   )=0.
      write (*,*) '... array layout completed'
css   endif                    ! end of nstep=0

      end subroutine reset_hycom_arrays
#endif
