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

      call pass_SCMDATA

      CALL CALC_PIJL(LM,P,PIJL)
      CALL CALC_AMPK(LM)

      call SCM_FORCN

      CALL tq_zmom_init(T,Q,PMID,PEDN)

      DO L=1,LM
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO

      CALL PGF_SCM(T,TZ,PIJL)

      call FCONV

      return
      END SUBROUTINE DYNAM

      SUBROUTINE SCM_FORCN
c     apply advective forcings from ARM Variational analysis to T and Q

      USE MODEL_COM,  only: DTSRC     
      USE ATM_COM,    only: P,T,Q,PK
      USE RESOLUTION, only: LM
      USE DYNAMICS,   only: SIG
      USE CONSTANT,   only: KAPA 
      USE CLOUDS,     only: SCM_DEL_T, SCM_DEL_Q
      USE SCMCOM,     only: SG_HOR_TMP_ADV, SG_VER_S_ADV, SG_HOR_Q_ADV,
     &              SG_VER_Q_ADV,iu_scm_prt,NSTEPSCM

      IMPLICIT NONE



      INTEGER L

cccccc is there some other variable they keep or function for doing this
      do L = 1,LM
         T(1,1,L) = T(1,1,L)*PK(L,1,1) 
c        write(iu_scm_prt,*) 'FORCN -old tq  ',L,T(I_TARG,J_TARG,L),
c    &               Q(I_TARG,J_TARG,L)*1000.0 
      enddo
     

      do L = 1,LM
c        write(iu_scm_prt,*) 'tadvs ',L,SG_HOR_TMP_ADV(L),
c    *                       SG_VER_S_ADV(L)
         SCM_DEL_T(L) = SG_HOR_TMP_ADV(L)*DTSRC + SG_VER_S_ADV(L)*DTSRC   
         T(1,1,L) = T(1,1,L) + SCM_DEL_T(L)
c        write(iu_scm_prt,*) 'add tadv delT T ',L,SCM_DEL_T(L),
c    &               T(I_TARG,J_TARG,L)    
      enddo   
      do L = 1,LM
c        write(iu_scm_prt,*) 'qadvs ',L,SG_HOR_Q_ADV(L),SG_VER_Q_ADV(L) 
         SCM_DEL_Q(L) = SG_HOR_Q_ADV(L)*DTSRC + SG_VER_Q_ADV(L)*DTSRC    
         Q(1,1,L) = Q(1,1,L) + SCM_DEL_Q(L)
         if (Q(1,1,L).lt.0.0) then
            write(99,51) NSTEPSCM,L,Q(1,1,L)
  51        format(1x,'SCM_FORCN NSTEP  L Q ',
     &               2(i5),f10.7) 
            SCM_DEL_Q(L) = -Q(1,1,L)
            Q(1,1,L) = 0.0
         endif
      enddo
    
      do L = 1,LM
c        write(iu_scm_prt,*) 'FORCN - new tq  ',L,T(I_TARG,J_TARG,L),
c    &               Q(I_TARG,J_TARG,L)*1000.0 
         T(1,1,L) = T(1,1,L)/PK(L,1,1)    
      enddo

      RETURN

      END SUBROUTINE SCM_FORCN 

  
      SUBROUTINE FCONV
C*****
C     for single column model
C     compute CONV=Horizontal Mass Convergence
C     as filled in subroutine AFLUX in the GCM for use in
C     the CONDSE and MSTCNV Subroutines
C     Use the Wind Divergence from the ARM data
C     CONV = Wind Divergence*dSigma*P*DelArea
C
C     NOTE:    Wind Divergence is calculated for the area of the
C              ARM site. Therefore we need to take into account the
C              difference between the GCM grid box area and the ARM
C              Site.   Oklahoma site (SGP)  300 x 365 KM = 109500KM**2
C                      GCM 2 x 2.5 degrees (smaller for SGP)
C                          ~ 222.63 * 2223.42 = 49739.01 KM**2
C                     SGP/GCM = 2.2
C                     
c              Note: for NSA  domain for the variational analysis
c                    is  230KM (longitudinal) x 100KM (latitudinal)
c                          230x100 = 23000
c                    GCM 2x2.5 degree grid box ~ 21266
c               area/box = (sin(q1)-sin(q2))*2(pi)R**2/144
c                        72 degrees-71degrees
c                    ARMFAC = NSA/GCM = 23000/21266 ~ 1.08
c   
c              What about for TWP site ? ? ?
c
c
c
  
      USE RESOLUTION, only: LM
      USE DYNAMICS,   only: DSIG
      USE ATM_COM,    only: P
      USe GEOM,       only: AXYP   
      USE SCMCOM,     only: SG_WINDIV, SG_CONV
   
      IMPLICIT NONE


      real*4 ARMFAC 
      integer L

      DATA ARMFAC/1.0/
c     DATA ARMFAC/2.2/
c     DATA ARMFAC/1.08/
      

c     want to fill SD (IDUM,JDUM)  check out 

      DO L=1,LM
         SG_CONV(L) = SG_WINDIV(L)*DSIG(L)*P(1,1)
     &                 *AXYP(1,1)*ARMFAC
      ENDDO

      return

      end SUBROUTINE FCONV  

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
      USE SCMCOM,     only: iu_scm_prt
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
c     do L=1,LM
c        write(iu_scm_prt,*) 'PGF_SCM  L GZ ',L,GZ(1,1,L)
c     enddo
C****
C
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
