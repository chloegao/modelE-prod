#include "rundeck_opts.h"

#if defined(TRACERS_OCEAN_INDEP) && !defined(RUNTIME_NTM_OCEAN)
/* TRACERS_OCEAN_INDEP enables the correct behaviors elsewhere */
/* in the model, but here we want to avoid the hard-coded */
/* specification of tracer names etc. if RUNTIME_NTM_OCEAN is set */
#define TRACERS_OCEAN_INDEP_HARDCODED
#endif

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      MODULE OCN_TRACER_COM
      SAVE
!@sum OCN_TRACER_COM: sets up tracer quantities for ocean tracers
!@+   This information can either be set up directly from the AGCM
!@+   or can be indpendently defined here.
!@+   Thus the ocean can have either more tracers or less than the AGCM 
!@+   depending on application.
!@+   The default behaviour will be to have the same number of 
!@+   TRACERS_WATER, but TRACERS_OCEAN will be independent.
!@+   Use "TRACERS_OCEAN" to allow ocean tracers.
!@+   Use "TRACERS_WATER" for freshwater tracers from ATM (this can be
!@+         used indpendently from TRACERS_OCEAN if a surface boundary 
!@+         condition is all that is required)
!@+   Use "TRACERS_OCEAN_INDEP" for independently defined ocn tracers
!@+        "TRACERS_AGE_OCEAN" is one partciularly case
!@param conc_from_fw definition for defining surface ocean conc
!@dbparam to_per_mil For printout of tracer concentration in permil


#ifdef TRACERS_OCEAN_INDEP_HARDCODED
C**** this defines tracer parameters that are local to ocean code


#ifdef TRACERS_AGE_OCEAN
      INTEGER, PARAMETER :: ntm=1
      CHARACTER*10 :: trname(ntm) = (/ 'Age       '/)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0
      INTEGER, DIMENSION(ntm) :: ntrocn=0
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .false.
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
#endif  /*TRACERS_AGE_OCEAN */

#ifdef TRACERS_ZEBRA
      INTEGER, PARAMETER ::  ntm=20
      CHARACTER*10 :: trname(ntm) =   (/ 'zebraL06   ', 'zebraL07   ',
     .     'zebraL08   ', 'zebraL09   ', 'zebraL10   ', 'zebraL11   ',
     .     'zebraL12   ', 'zebraL13   ', 'zebraL14   ', 'zebraL15   ',
     .     'zebraL16   ', 'zebraL17   ', 'zebraL18   ', 'zebraL19   ',
     .     'zebraL20   ', 'zebraL21   ', 'zebraL22   ', 'zebraL23   ',
     .     'zebraL24   ', 'zebraL26   '/)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0
      INTEGER, DIMENSION(ntm) :: ntrocn=0
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .false.
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
#endif    /*zebra*/

#ifdef TRACERS_OceanBiology
#ifdef TRACERS_Alkalinity
      INTEGER, PARAMETER :: ntm=16
#else
      INTEGER, PARAMETER :: ntm=15
#endif
      CHARACTER*10 :: trname(ntm) = (/ 'Nitr      ', 'Ammo      ', 
     .     'Sili      ', 'Iron      ', 'Diat      ', 'Chlo      ', 
     .     'Cyan      ', 'Cocc      ', 'Herb      ', 'Inert     ',
     .     'N_det     ', 'S_det     ', 'I_det     ', 'DOC       ',
     .     'DIC       '
#ifdef TRACERS_Alkalinity
     .     ,'Alk       '
#endif
     .     /)
      REAL*8, DIMENSION(ntm) :: trw0=0, trdecay=0
!@var ntrocn scaling exponent for tracers
      INTEGER, DIMENSION(ntm) :: ntrocn = (/ -4,-6,-4,-8,-8,-8,-8,-8,-8,
     .     -4,-6,-6,-10,-6,-3
#ifdef TRACERS_Alkalinity
     .     ,-6
#endif
     .     /)
      INTEGER, DIMENSION(NTM) :: to_per_mil = 0
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .false.
#endif  /* TRACERS_OceanBiology */

C**** a default tracer if nothing else is defined
#if ! (defined TRACERS_OceanBiology || defined TRACERS_AGE_OCEAN || defined TRACERS_ZEBRA)
      INTEGER, PARAMETER :: ntm=1
      CHARACTER*10 :: trname(ntm) = (/ 'Water     '/)
      REAL*8, DIMENSION(ntm) :: trw0=1d0, trdecay=0
      INTEGER, DIMENSION(ntm) :: ntrocn=2
      LOGICAL, DIMENSION(NTM) :: conc_from_fw = .true.
#endif

      LOGICAL, DIMENSION(ntm) :: t_qlimit=.true. 
      LOGICAL, DIMENSION(ntm) :: need_ic=.false.
      INTEGER, DIMENSION(ntm) :: itime_tr0 = 0
      INTEGER :: n_water = 0

#else   /* not TRACERS_OCEAN_INDEP_HARDCODED */

C**** These arrays are allocated/intialized in one of two ways:
C**** (1) by the AGCM, which copies its data into them
C**** (2) if RUNTIME_NTM_OCEAN is defined, the allocation/initialization
C****     of these arrays happens in alloc_ocn_tracer_com() using
C****     information from the rundeck and/or other config files.
C****     Set ocean_trname='name1 name2 ...' to instantiate
C****     a given number of tracers, whose ICs will be read
C****     from the file OCN_TRACER_CONFIG.  This approach is
C****     still being tailored to handle all cases for which
C****     it will prove useful.
      INTEGER :: ntm, n_water
      INTEGER, DIMENSION(:), ALLOCATABLE ::
     &     itime_tr0, ntrocn, to_per_mil
      LOGICAL, DIMENSION(:), ALLOCATABLE ::
     &     t_qlimit, conc_from_fw, need_ic
      REAL*8, DIMENSION(:), ALLOCATABLE ::
     &     trdecay, trw0
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE ::
     &     trname

#endif
      INTEGER :: n_age=0, n_obio=0, n_vent=0, n_wms1=0, n_wms2=0
     .          ,n_wms3=0,n_dets,n_cfc,n_dic


      REAL*8, allocatable, DIMENSION(:) :: expDecayRate

#ifdef TRACERS_SPECIAL_O18
      ! could be made a local var in tracer_ic_ocean?
      integer :: water_tracer_ic=1 ! Read water tracers ic from H2O18ic (=1) or set all to SMOW (=0)
#endif

      END MODULE OCN_TRACER_COM

      SUBROUTINE initOcnTracerCom(DTS)
      use ocn_tracer_com, only: ntm, expDecayRate, trdecay
      real*8, intent(in) :: dts 
      integer n
      allocate(expDecayRate(ntm))
      do n=1,ntm
        if (trdecay(n).gt.0.0) expDecayRate(n)=exp(-trdecay(n)*dts)
      end do
      END SUBROUTINE initOcnTracerCom

      subroutine alloc_ocn_tracer_com
      USE DOMAIN_DECOMP_1D, only : am_i_root
      use ocn_tracer_com
      use Dictionary_mod, only : sync_param,is_set_param,get_param
      ! TODO: it would be better long-term if tracer arrays were
      ! owned by ocn_tracer_com.
      USE OCEAN, only : TRMO
     *       ,oc_tracer_mean
      USE OCEAN, only : use_qus,
     &     TXMO,TYMO,TZMO, TXXMO,TYYMO,TZZMO, TXYMO,TYZMO,TZXMO
      USE OCEANRES, only : IM=>IMO, JM=>JMO, LMO 
      USE OCEANR_DIM, only : J_0H,J_1H
      USE STRAITS, only : NMST,TRMST,TXMST,TZMST,TRME,TXME,TYME,TZME
#ifdef TRACERS_WATER
      ! The ocean model should not have to know how many layers
      ! the sea ice model uses - this dependence is unfriendly
      ! to componentization and will be eliminated at some point,
      ! if the array TRSIST is not eliminated first (as of 11/2201
      ! TRSIST is inactive but still needs to be allocated).  -M.K.
      USE STRAITS, only : TRSIST
      USE SEAICE, only : LMI
#endif
      implicit none
      character(len=128) :: trname_list
      character(len=10) :: trname_(50)
      integer :: i,ier
      integer :: img, jmg, lmg, n
      if (am_i_root()) then
        img = im
        jmg = jm
        lmg = lmo
      else
        img = 1
        jmg = 1
        lmg = 1
      end if
#ifdef RUNTIME_NTM_OCEAN
      if(is_set_param("ocean_trname")) then
        trname_list=''
        call get_param("ocean_trname",trname_list)
        trname_list=adjustl(trname_list)
        ntm=0
        do while(len_trim(trname_list).gt.0)
          i=index(trname_list,' ')
          ntm = ntm + 1
          trname_(ntm) = trname_list(1:i-1)
          trname_list = adjustl(trname_list(i:128))
        enddo
        allocate(trname(ntm)); trname(:) = trname_(1:ntm)
        ! RUNTIME_NTM_OCEAN is currently only being used for
        ! a few simple tracers, so just set defaults for those.
        allocate(trw0(ntm)); trw0(:) = 0.
        allocate(trdecay(ntm)); trdecay(:) = 0.
        allocate(to_per_mil(ntm)); to_per_mil = 0
        allocate(t_qlimit(ntm)); t_qlimit=.true.
        allocate(conc_from_fw(ntm)); conc_from_fw=.false.
        allocate(need_ic(ntm)); need_ic=.false.
      else
        call stop_model('RUNTIME_NTM_OCEAN needs ocean_trname',255)
      endif
#endif
#ifndef TRACERS_OCEAN_INDEP
      allocate(
     &     itime_tr0(ntm), ntrocn(ntm), to_per_mil(ntm),
     &     t_qlimit(ntm), conc_from_fw(ntm),
     &     trdecay(ntm), trw0(ntm)
     &     )
#endif
      allocate(oc_tracer_mean(ntm)) 
      oc_tracer_mean(:) = -999.

      ALLOCATE(TRMST(LMO,NMST,NTM),
     &         TXMST(LMO,NMST,NTM),
     &         TZMST(LMO,NMST,NTM),
     &         TRME(2,NMST,LMO,NTM),
     &         TXME(2,NMST,LMO,NTM),
     &         TYME(2,NMST,LMO,NTM),
     &         TZME(2,NMST,LMO,NTM)
     &        )
      trmst = 0.
      txmst = 0.
      tzmst = 0.
      trme = 0.
      txme = 0.
      tyme = 0.
      tzme = 0.
#ifdef TRACERS_WATER
      ALLOCATE(TRSIST(NTM,LMI,NMST))
      trsist = 0.
#endif

      ALLOCATE( TRMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TXMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TYMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TZMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      trmo = 0.
      txmo = 0.
      tymo = 0.
      tzmo = 0.

      if(use_qus==1) then
      ALLOCATE( TXXMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TYYMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TZZMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TXYMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TYZMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      ALLOCATE( TZXMO(IM,J_0H:J_1H,LMO,NTM), STAT = IER)
      txxmo=0.; tyymo=0.; tzzmo=0.; txymo=0.; tyzmo=0.; tzxmo=0.
      endif

      do n=1,ntm
       if (trname(n).eq.'OceanAge') n_age = n
       if (trname(n).eq.'Ventilatn') n_vent = n
       if (trname(n).eq.'WatrMass1') n_wms1 = n
       if (trname(n).eq.'WatrMass2') n_wms2 = n
       if (trname(n).eq.'WatrMass3') n_wms3 = n
       if (trname(n).eq.'DetSet') n_dets = n
       if (trname(n).eq.'CFC') n_cfc = n
       if (trname(n).eq.'DIConly') then
                  n_dic = n
                 need_ic(n_dic)=.true.
       endif
      enddo

      return
      end subroutine alloc_ocn_tracer_com

#endif
