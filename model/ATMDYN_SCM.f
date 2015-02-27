#include "rundeck_opts.h"

      module ATMDYN
      implicit none

      contains

      SUBROUTINE init_ATMDYN
      return
      end SUBROUTINE init_ATMDYN

      SUBROUTINE DYNAM
      USE RESOLUTION, only: im,lm,ls1
      USE SOMTQ_COM,  only: tmom,mz
      USE ATM_COM,    only: t,p,q,PMID,PEDN,MUs,MVs,MWs
      USE DOMAIN_DECOMP_ATM, only : grid

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     TZ,PIJL

      INTEGER L

      do L=1,LM
         MUs(:,:,L) = 0.
         MVs(:,:,L) = 0.
         MWs(:,:,L) = 0.
      ENDDO

      call update_SCM_inputs

      CALL CALC_PIJL(LM,P,PIJL)
      CALL CALC_AMPK(LM)

      call SCM_FORCN
  
      CALL tq_zmom_init(T,Q,PMID,PEDN)

      DO L=1,LM
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO

      CALL PGF_SCM(T,TZ,PIJL)

      return
      END SUBROUTINE DYNAM

      SUBROUTINE SCM_FORCN
c     apply large-scale forcings to T, Q, U, V

      USE MODEL_COM,  only: DTSRC
      USE ATM_COM,    only: P,T,Q,PK,U,V,PMID
      USE RESOLUTION, only: LM,PTOP
      USE DYNAMICS,   only: SIG
      USE CONSTANT,   only: KAPA, OMEGA, GRAV, RGAS
      USE GEOM, only : sinlat2d
      USE SCM_COM,    only: SCMopt,SCMin
#ifdef CACHED_SUBDD
      USE DOMAIN_DECOMP_ATM, only : grid,getDomainBounds
      use subdd_mod, only : subdd_groups,subdd_ngroups,subdd_type
     &     ,inc_subdd,find_groups
#endif

      IMPLICIT NONE

      real*8, dimension(LM) :: Tabs(LM),
     &                         SCM_ver_u_adv,SCM_ver_v_adv,
     &                         SCM_nudge_T,SCM_nudge_Q,
     &                         SCM_force_T,SCM_force_Q
      real*8 f_cor

      INTEGER L
#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      integer :: i,i_0,i_1,j,j_0,j_1,k
      type(subdd_type), pointer :: subdd
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &                  sddarr3d
      CALL getDomainBounds(grid,
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
#endif

c     operate on absolute temperature
      do L = 1,LM
        Tabs(L) = T(1,1,L)*PK(L,1,1)
      enddo

      SCM_force_T = 0.
      SCM_force_Q = 0.
      SCM_nudge_T = 0.
      SCM_nudge_Q = 0.

      do L = 1,LM

c       fix winds if specified and no Coriolis acceleration

        if( SCMopt%wind .and. .not. SCMopt%geo )then
          U(1,1,L) = SCMin%U(L)
          V(1,1,L) = SCMin%V(L)
        endif

c       large-scale forcings

        if( SCMopt%omega .or. SCMopt%w )then
c       *** apply omega defined at layer bottom to upwind gradient

          if ( L < LM ) then ! omega assumed zero at top of layer LM
            if ( SCMin%Omega(L+1) > 0. ) then ! upwind gradient above
              SCMin%SadvV(L) = -SCMin%Omega(L)*
     &           (T(1,1,L+1)-T(1,1,L))*PK(L,1,1)/
     &           (PMID(L+1,1,1)-PMID(L,1,1))
              SCMin%QadvV(L) = -SCMin%Omega(L)*
     &           (Q(1,1,L+1)-Q(1,1,L))/
     &           (PMID(L+1,1,1)-PMID(L,1,1))
            endif
          endif

          if ( L > 1 ) then ! no atmospheric gradient through surface
            if ( SCMin%Omega(L) < 0. ) then ! upwind gradient below
              SCMin%SadvV(L) = -SCMin%Omega(L)*
     &           (T(1,1,L)-T(1,1,L-1))*PK(L,1,1)/
     &           (PMID(L,1,1)-PMID(L-1,1,1))
              SCMin%QadvV(L) = -SCMin%Omega(L)*
     &           (Q(1,1,L)-Q(1,1,L-1))/
     &           (PMID(L,1,1)-PMID(L-1,1,1))
            endif
          endif

        else if( .not. SCMopt%ls_v )then
c       *** otherwise no LS vertical flux divergence if not specified
          SCMin%SadvV(L) = 0.
          SCMin%QadvV(L) = 0.
        endif

        if( .not. SCMopt%ls_h )then
c       *** no vertical forcings
          SCMin%TadvH(L) = 0.
          SCMin%QadvH(L) = 0.
        endif

        SCM_force_T(L) = (SCMin%TadvH(L)+SCMin%SadvV(L))*DTSRC
        SCM_force_Q(L) = (SCMin%QadvH(L)+SCMin%QadvV(L))*DTSRC

        Tabs(L) = Tabs(L) + SCM_force_T(L)
        Q(1,1,L) = Q(1,1,L) + SCM_force_Q(L)

        if( Q(1,1,L) < 0. )then
          SCM_force_Q(L) = -Q(1,1,L)
          Q(1,1,L) = 0.0
        endif

c       apply nudging terms

        if( SCMopt%nudge )then
c       *** calculate nudging toward observed profile

          SCM_nudge_T(L) = (SCMin%T(L)-Tabs(L))/SCMopt%tau*DTSRC
          SCM_nudge_Q(L) = (SCMin%Q(L)-Q(1,1,L))/SCMopt%tau*DTSRC

          if( SCMopt%Fnudge )then
            SCM_nudge_T(L) = SCM_nudge_T(L)*SCMin%Fnudge(L)
            SCM_nudge_Q(L) = SCM_nudge_Q(L)*SCMin%Fnudge(L)
          endif

          Tabs(L) = Tabs(L) + SCM_nudge_T(L)
          Q(1,1,L) = Q(1,1,L) + SCM_nudge_Q(L)

        endif

      enddo

c     *** apply changes to actual prognostic variable (potential temperature)
      do L = 1,LM
        T(1,1,L) = Tabs(L)/PK(L,1,1)
      enddo

#ifdef CACHED_SUBDD
C****
C**** Collect some high-frequency outputs
C****
      call find_groups('fijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('dq_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_Q(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_T(l)/PK(l,1,1)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dq_nudge')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_nudge_Q(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_nudge')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_nudge_T(l)/PK(l,1,1)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      end select
      enddo
      enddo
#endif

c     when Coriolis forcing used (computed from geostrophic winds),
c     also possibly apply vertical advection to horizontal winds

      if ( SCMopt%geo ) then

        f_cor = 2.*omega*sinlat2d(1,1)

        SCM_ver_u_adv = 0.
        SCM_ver_v_adv = 0.

        do L = 1,LM

          if( SCMopt%VadvHwind )then
c         *** apply omega defined at layer bottom to upwind gradient

            if ( L < LM ) then ! omega assumed zero at top of layer LM
              if ( SCMin%Omega(L+1) > 0. ) then ! upwind gradient above
                 SCM_ver_u_adv(L) = -SCMin%Omega(L)*
     &              (U(1,1,L+1)-U(1,1,L))/
     &              (PMID(L+1,1,1)-PMID(L,1,1))
                 SCM_ver_v_adv(L) = -SCMin%Omega(L)*
     &              (V(1,1,L+1)-V(1,1,L))/
     &              (PMID(L+1,1,1)-PMID(L,1,1))
              endif
            endif

            if ( L > 1 ) then ! no atmospheric gradient through surface
              if ( SCMin%Omega(L) < 0. ) then ! upwind gradient below
                 SCM_ver_u_adv(L) = -SCMin%Omega(L)*
     &              (U(1,1,L)-U(1,1,L-1))/
     &              (PMID(L,1,1)-PMID(L-1,1,1))
                 SCM_ver_v_adv(L) = -SCMin%Omega(L)*
     &              (V(1,1,L)-V(1,1,L-1))/
     &              (PMID(L,1,1)-PMID(L-1,1,1))
              endif
            endif
          endif

c         apply combined forcings to horizontal winds
          U(1,1,L) = U(1,1,L) +
     &      ( SCM_ver_u_adv(L) +
     &        f_cor*(V(1,1,L)-SCMin%Vg(L)) )*dtsrc
          V(1,1,L) = V(1,1,L) +
     &      ( SCM_ver_v_adv(L) -
     &        f_cor*(U(1,1,L)-SCMin%Ug(L)) )*dtsrc

        enddo      ! L = 1,LM
      endif        ! use geostrophic winds for Coriolis forcing

      return
      END SUBROUTINE SCM_FORCN 

  
      SUBROUTINE SDRAG(DT1)
      REAL*8, INTENT(IN) :: DT1 
      return
      END SUBROUTINE SDRAG 


      SUBROUTINE PGF_SCM (T,SZ,P)
!@SCM-version    For SCM need to calculate geopotential height. 
!                Remove other calculations.
!@sum  PGF Adds pressure gradient forces to momentum
!@auth Original development team
      USE CONSTANT,   only: grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE RESOLUTION, only: im,jm,lm,ls1,psfmpt,ptop
      USE ATM_COM,    only: zatmo, gz, phi
      USE DYNAMICS,   only: sig,bydsig,do_polefix,
     *     dsig,sige,pu,spa
      IMPLICIT NONE

      REAL*8, DIMENSION(1,1,LM):: T
      REAL*8, DIMENSION(1,1,LM) :: P, SZ

      REAL*8 PKE(LS1:LM+1)
      REAL*8 PIJ,PDN,PKDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     *     ,BYDP
      REAL*8 TZBYDP,FLUX,FDNP,FDSP,RFDU,PHIDN,FACTOR
      INTEGER I,J,L,IM1,IP1,IPOLE  !@var I,J,IP1,IM1,L,IPOLE loop variab.

C****
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
C****
C**** VERTICAL DIFFERENCING
C****
      DO L=LS1,LM
      SPA(:,:,L)=0.
      END DO

      DO J=1,1
      DO I=1,1
        PIJ=P(I,J,1)
        PDN=PIJ+PTOP
        PKDN=PDN**KAPA
        PHIDN=ZATMO(I,J)
C**** LOOP OVER THE LAYERS
        DO L=1,LM
          PKPDN=PKDN*PDN
          PKPPDN=PKPDN*PDN
          IF(L.GE.LS1) THEN
            DP=DSIG(L)*PSFMPT
            BYDP=1./DP
            P0=SIG(L)*PSFMPT+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PSFMPT+PTOP
            PKUP=PKE(L+1)
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
          ELSE
            DP=DSIG(L)*PIJ
            BYDP=1./DP
            P0=SIG(L)*PIJ+PTOP
            TZBYDP=2.*SZ(I,J,L)*BYDP
            X=T(I,J,L)+TZBYDP*P0
            PUP=SIGE(L+1)*PIJ+PTOP
            PKUP=PUP**KAPA
            PKPUP=PKUP*PUP
            PKPPUP=PKPUP*PUP
C****   CALCULATE SPA, MASS WEIGHTED THROUGHOUT THE LAYER
            SPA(I,J,L)=RGAS*((X+TZBYDP*PTOP)*(PKPDN-PKPUP)*BYKAPAP1
     *      -X*PTOP*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP
          END IF
C**** CALCULATE PHI, MASS WEIGHTED THROUGHOUT THE LAYER
          PHI(I,J,L)=PHIDN+RGAS*(X*PKDN*BYKAPA-TZBYDP*PKPDN*BYKAPAP1
     *      -(X*(PKPDN-PKPUP)*BYKAPA-TZBYDP*(PKPPDN-PKPPUP)*BYKAPAP2)
     *      *BYDP*BYKAPAP1)
C**** CALULATE PHI AT LAYER TOP (EQUAL TO BOTTOM OF NEXT LAYER)
          PHIDN=PHIDN+RGAS*(X*(PKDN-PKUP)*BYKAPA-TZBYDP*(PKPDN-PKPUP)
     *     *BYKAPAP1)
          PDN=PUP
          PKDN=PKUP
        END DO
      END DO
      END DO

      DO L=1,LM
        GZ(:,:,L)=PHI(:,:,L)
      END DO

      RETURN
      END SUBROUTINE PGF_SCM


c     SUBROUTINE AFLUX (U,V,PIJL)
c     END SUBROUTINE AFLUX


C**** Dummy routines

      SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS( MUs,MVs,dt)
!@sum COMPUTE_DYNAM_AIJ_DIAGNOSTICS Dummy
      use DOMAIN_DECOMP_ATM, only: grid

      real*8, intent(in) :: MUs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: MVs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      return
      END SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS

      end module ATMDYN

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates A-grid column-sum atmospheric kinetic energy,
!@sum  multiplied by cell area
!@auth Gary Russell/Gavin Schmidt
      IMPLICIT NONE

      REAL*8, DIMENSION(1,1) :: RKE

      RKE = 0.
      !call stop_model('calculate a-grid value instead',255)

      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE calc_kea_3d(kea)
!@sum  calc_kea_3d calculates square of wind speed on the A grid
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      IMPLICIT NONE
      REAL*8, DIMENSION(1,1,LM) :: KEA

      RETURN

      END SUBROUTINE calc_kea_3d

      subroutine recalc_agrid_uv
      USE ATM_COM,    only: u,v
      USE ATM_COM,    only: ua=>ualij,va=>valij
      implicit none

      ua(:,1,1)=u(1,1,:)
      va(:,1,1)=v(1,1,:)

      return
      end subroutine recalc_agrid_uv

      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,1,1) :: UR,VR
      real*8, dimension(1,lm) :: ursp,vrsp,urnp,vrnp ! not used
      integer :: l
      if(k.ne.1)
     &     call stop_model('incorrect k in replicate_uv_to_agrid',255)
      do l=1,lm
        ur(1,l,1,1) = u(1,1,l)
        vr(1,l,1,1) = v(1,1,l)
      enddo ! l
      return
      end subroutine replicate_uv_to_agrid

      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,1,1) :: DU,DV
      real*8, dimension(1,lm) :: dusp,dvsp,dunp,dvnp ! not used
      integer :: l

      if(k.ne.1) call stop_model(
     &     'incorrect k in avg_replicated_duv_to_vgrid',255)

      do l=1,lm
        u(1,1,l)=u(1,1,l)+du(1,l,1,1)
        v(1,1,l)=v(1,1,l)+dv(1,l,1,1)
      enddo ! l

      return
      end subroutine avg_replicated_duv_to_vgrid

      SUBROUTINE QDYNAM
      return
      END SUBROUTINE QDYNAM

#ifdef CACHED_SUBDD
      subroutine fijlh_defs(arr,nmax,decl_count)
c
c 3D outputs
c
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use model_com, only: dtsrc
      use constant, only: kapa
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c
      decl_count = 0
c
      arr(next()) = info_type_(
     &  sname = 'dq_ls',
     &  lname = 'moisture tendency from large-scale forcings',
     &  units = 'kg/kg/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_ls',
     &  lname = 'theta tendency from large-scale forcings',
     &  units = 'K/day',
     &  scale = 1000.**kapa/dtsrc*SECONDS_PER_DAY
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dq_nudge',
     &  lname = 'moisture tendency from nudging',
     &  units = 'kg/kg/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_nudge',
     &  lname = 'theta tendency from nudging',
     &  units = 'K/day',
     &  scale = 1000.**kapa/dtsrc*SECONDS_PER_DAY
     &     )
c
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine fijlh_defs
#endif
