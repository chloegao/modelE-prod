#include "rundeck_opts.h"

      module ocn_tracer_entry_mod
      implicit none
      type ocn_tracer_entry
       character*10 :: trname
       integer :: trw0=0, trdecay=0, ntrocn=0, to_per_mil=0
       integer :: itime_tr0=0, ntrocn_delta=0
       logical :: conc_from_fw=.false., t_qlimit=.true., need_ic=.false.
       logical :: from_file=.false.
       integer, dimension(:), allocatable :: con_point_idx
       character(len=10), dimension(:), allocatable :: con_point_str
      end type ocn_tracer_entry
      end module ocn_tracer_entry_mod

      module ocn_tracer_vector_mod
      use ocn_tracer_entry_mod
#define _entry type(ocn_tracer_entry)
#include "containers/vector.fh"
      end module ocn_tracer_vector_mod

#if (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      MODULE OCN_TRACER_COM

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


      use ocn_tracer_entry_mod
      use ocn_tracer_vector_mod, only: vector_ocn_tracer_entry=>vector
      SAVE
      type(vector_ocn_tracer_entry) :: tracerlist

      integer :: n_water
      INTEGER :: n_age=0, n_obio=0, n_vent=0, n_wms1=0, n_wms2=0
     .          ,n_wms3=0,n_dets,n_cfc,n_dic,n_gasx


      REAL*8, allocatable, DIMENSION(:) :: expDecayRate

#ifdef TRACERS_SPECIAL_O18
      ! could be made a local var in tracer_ic_ocean?
      integer :: water_tracer_ic=1 ! Read water tracers ic from H2O18ic (=1) or set all to SMOW (=0)
#endif

      contains

      subroutine add_ocn_tracer(i_trname, i_trw0, i_ntrocn,
     &              i_ntrocn_delta, i_conc,
     &              i_from_file, i_con_point_idx, i_con_point_str)
      use ocn_tracer_entry_mod
      implicit none
      character(len=*), intent(in) :: i_trname
      integer, intent(in), optional :: i_trw0, i_ntrocn, i_ntrocn_delta
      logical, intent(in), optional :: i_conc, i_from_file
      integer, dimension(:), intent(in), optional :: i_con_point_idx
      character*10, dimension(:), intent(in), optional:: i_con_point_str
      type(ocn_tracer_entry) :: entry
      integer :: numpts

      entry%trname=i_trname
      if (present(i_trw0)) entry%trw0=i_trw0
      if (present(i_ntrocn)) entry%ntrocn=i_ntrocn
      if (present(i_ntrocn_delta)) then
        entry%ntrocn_delta=i_ntrocn_delta
      else
        entry%ntrocn_delta=entry%ntrocn-10 ! default is 10 orders of magnitude less
      endif
      if (present(i_conc)) entry%conc_from_fw=i_conc
      if (present(i_from_file)) entry%from_file=i_from_file
      numpts=0
      if (present(i_con_point_idx)) numpts=size(i_con_point_idx)
      allocate(entry%con_point_idx(numpts))
      allocate(entry%con_point_str(numpts))
      if (numpts>0) then
        entry%con_point_idx=i_con_point_idx
        entry%con_point_str=i_con_point_str
      endif
      call tracerlist%push_back(entry)
      return
      end subroutine add_ocn_tracer


      END MODULE OCN_TRACER_COM


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
      integer :: i,ier
      integer :: img, jmg, lmg, n
      integer :: numtracers
      type(ocn_tracer_entry), pointer :: entry

      if (am_i_root()) then
        img = im
        jmg = jm
        lmg = lmo
      else
        img = 1
        jmg = 1
        lmg = 1
      end if
      if(is_set_param("ocean_trname")) then
        trname_list=''
        call get_param("ocean_trname",trname_list)
        trname_list=adjustl(trname_list)
        do while(len_trim(trname_list).gt.0)
          i=index(trname_list,' ')
          call add_ocn_tracer(trname_list(1:i-1), i_from_file=.true.)
          trname_list = adjustl(trname_list(i:128))
        enddo
      endif
      numtracers=tracerlist%getsize()
      allocate(oc_tracer_mean(numtracers)) 
      oc_tracer_mean(:) = -999.

      ALLOCATE(TRMST(LMO,NMST,numtracers),
     &         TXMST(LMO,NMST,numtracers),
     &         TZMST(LMO,NMST,numtracers),
     &         TRME(2,NMST,LMO,numtracers),
     &         TXME(2,NMST,LMO,numtracers),
     &         TYME(2,NMST,LMO,numtracers),
     &         TZME(2,NMST,LMO,numtracers)
     &        )
      trmst = 0.
      txmst = 0.
      tzmst = 0.
      trme = 0.
      txme = 0.
      tyme = 0.
      tzme = 0.
#ifdef TRACERS_WATER
      ALLOCATE(TRSIST(numtracers,LMI,NMST))
      trsist = 0.
#endif

      ALLOCATE( TRMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TXMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TYMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TZMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      trmo = 0.
      txmo = 0.
      tymo = 0.
      tzmo = 0.

      if(use_qus==1) then
      ALLOCATE( TXXMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TYYMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TZZMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TXYMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TYZMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      ALLOCATE( TZXMO(IM,J_0H:J_1H,LMO,numtracers), STAT = IER)
      txxmo=0.; tyymo=0.; tzzmo=0.; txymo=0.; tyzmo=0.; tzxmo=0.
      endif

      do n=1,numtracers
       entry=>tracerlist%at(n)
       if (entry%trname.eq.'OceanAge') n_age = n
       if (entry%trname.eq.'Ventilatn') n_vent = n
       if (entry%trname.eq.'GASX') n_gasx = n
       if (entry%trname.eq.'WatrMass1') n_wms1 = n
       if (entry%trname.eq.'WatrMass2') n_wms2 = n
       if (entry%trname.eq.'WatrMass3') n_wms3 = n
       if (entry%trname.eq.'DetSet') n_dets = n
       if (entry%trname.eq.'aoCFC') n_cfc = n
!       if (entry%trname.eq.'DIConly') then
!         n_dic = n
!         entry%need_ic=.true.
!       endif
      enddo

      return
      end subroutine alloc_ocn_tracer_com
#endif
