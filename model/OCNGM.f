#include "rundeck_opts.h"

      MODULE GM_COM
!@sum  GM_COM variables related to GM isopycnal and Redi fluxes
!@auth Gavin Schmidt/Dan Collins
!@ver  2009/02/13
      USE CONSTANT, only : radius,omega,grav
      USE OCEAN, only : im,jm,lmo,lmm,lmu,lmv,dts,cospo,sinpo,ze,dxypo
     *     ,mo,dypo,dyvo,dxpo,dzo
      USE KPP_COM, only : kpl
      USE OCEAN_DYN, only  : dh  ! ,vbar
      USE ODIAG, only : oij=>oij_loc,ij_gmsc 

      USE DOMAIN_DECOMP_1D, ONLY : GETDomainBounds, HALO_UPDATE, NORTH,
     *                          SOUTH, PACK_DATA, AM_I_ROOT,
     *                          UNPACK_DATA, GLOBALSUM
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE

      SAVE

!@var AI0,AI1,AI2,AI3 Cmponents of GM mixing coeff = F(isopycnal slopes)
!@var SIX0-SIX3,SIY0-SIY3: Slopes calculated from 4 triads of density.
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     ASX0,ASX1,ASX2,ASX3,AIX0,AIX1,AIX2,AIX3,
     *     ASY0,ASY1,ASY2,ASY3,AIY0,AIY1,AIY2,AIY3,
     *     S2X0,S2X1,S2X2,S2X3,S2Y0,S2Y1,S2Y2,S2Y3

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     DXZ, BXZ, CXZ,CDXZ, EXZ, DYZ,BYZ,
     *     CYZ,CDYZ,EYZ,CEXZ,CEYZ

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     BXX, BYY, BZZ, AZX, BZX, CZX,AEZX,
     *     EZX,CEZX, AZY, BZY,  CZY,AEZY, EZY,CEZY

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     BYDZV,BYDH,DZV

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     RHOX,RHOY,RHOMZ,BYRHOZ

!@var AINV Calculated Isopycnal thickness diffusion (m^2/s)
!@var ARIV Calculated Redi diffusion (m^2/s)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     AINV,ARIV

!@var ARAI Scaling for Redi diffusion terms to GM diffusion term (1)
      !REAL*8, PARAMETER :: ARAI = 1d0
!@var RGMI ratio of GM (thickness) diffusivity to isoneutral (Redi) diffusivity
      REAL*8 :: RGMI
!@var QCROSS true if cross terms should be calculated
      !LOGICAL, PARAMETER :: QCROSS = .NOT. (ARAI.eq.1d0) ! i.e..FALSE.
      LOGICAL :: QCROSS
!@var AMU = Visbeck scheme scaling parameter (1)
      REAL*8, PARAMETER :: AMU = 0.13d0
!@var SLIM = Upper limit of isopycnal slopes (stability parameter)
      REAL*8, PARAMETER :: SLIM=2d-3, BYSLIM=1./SLIM

      REAL*8, PARAMETER :: eps=TINY(1.D0)

#if (defined CONSTANT_MESO_DIFFUSIVITY)
!@dbparam meso_diffusivity_const constant diffusivity (m2/s)
!@+       for sensitivity studies
      real*8 :: meso_diffusivity_const
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
!@dbparam meso_diffusivity_z0
!@dbparam meso_diffusivity_zscale
!@dbparam meso_diffusivity_deep
      real*8 ::
     & meso_diffusivity_z0,meso_diffusivity_zscale,meso_diffusivity_deep
#elif (defined CONSTANT_MESO_LENSCALE)
!@dbparam meso_lenscale_const a fixed length scale (meters) to use in
!@+       lieu of Rossby radius RD in AINV = AMU * RD**2 * BYTEADY
      real*8 :: meso_lenscale_const
#endif
      REAL*8, ALLOCATABLE, DIMENSION(:) ::
     *     BYDYP,BYDXP,BYDYV

!@var VBAR specific volume (ref to mid point pressure)
!@var dVBARdZ specific volume vertical difference (ref to lower point pressure)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: VBAR, dVBARdZ,G3D,S3D,P3D
!@var RHOZ1K density gradient over top 1km
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: RHOZ1K

!@var LUP level corresponding to 1km depth
      INTEGER :: LUP

      contains

      SUBROUTINE ALLOC_GM_COM
        implicit none
c**** allocate arrays
        allocate( ASX0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASX3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIX3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ASY3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AIY3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2X3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y0(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y1(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y2(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S2Y3(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( DXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CDXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CDYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEXZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEYZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( BXX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZZ (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AEZX(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EZX (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEZX(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( AEZY(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( EZY (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( CEZY(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( BYDZV(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYDH (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DZV  (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( VBAR   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( DVBARDZ(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( G3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( S3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( P3D   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOZ1K (IM,grid%j_strt_halo:grid%j_stop_halo) )

        allocate( RHOX   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOY   (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( RHOMZ  (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( BYRHOZ (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( AINV (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )
        allocate( ARIV (IM,grid%j_strt_halo:grid%j_stop_halo,LMO) )

        allocate( BYDYP (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDXP (grid%j_strt_halo:grid%j_stop_halo) )
        allocate( BYDYV (grid%j_strt_halo:grid%j_stop_halo) )

      END SUBROUTINE ALLOC_GM_COM

      END MODULE GM_COM

      SUBROUTINE GMKDIF(RGMI_in)
!@sum GMKDIF calculates density gradients and tracer operators
!@+   for Redi and GM isopycnal skew fluxes
!@auth Dan Collins/Gavin Schmidt
      USE GM_COM
      USE Dictionary_mod
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: RGMI_in

      INTEGER I,J,L,IM1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      INTEGER, SAVE :: IFIRST = 1

      if ( IFIRST.ne.0 ) then
        IFIRST = 0
        call ALLOC_GM_COM
#if (defined CONSTANT_MESO_DIFFUSIVITY)
       call get_param('meso_diffusivity_const',meso_diffusivity_const)
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
       call get_param('meso_diffusivity_z0',meso_diffusivity_z0)
       call get_param('meso_diffusivity_zscale',meso_diffusivity_zscale)
       call get_param('meso_diffusivity_deep',meso_diffusivity_deep)
#elif (defined CONSTANT_MESO_LENSCALE)
       call get_param( 'meso_lenscale_const', meso_lenscale_const )
#endif
      endif

      RGMI = RGMI_in
      QCROSS = .NOT. (RGMI.eq.1d0)

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H)

C**** Initialize SLOPES common block of coefficients
      DO L=1,LMO
        ASX0(:,:,L)=0. ;ASX1(:,:,L)=0. ; ASX2(:,:,L)=0. ; ASX3(:,:,L)=0.
        AIX0(:,:,L)=0. ;AIX1(:,:,L)=0. ; AIX2(:,:,L)=0. ; AIX3(:,:,L)=0.
        ASY0(:,:,L)=0. ;ASY1(:,:,L)=0. ; ASY2(:,:,L)=0. ; ASY3(:,:,L)=0.
        AIY0(:,:,L)=0. ;AIY1(:,:,L)=0. ; AIY2(:,:,L)=0. ; AIY3(:,:,L)=0.
        S2X0(:,:,L)=0. ;S2X1(:,:,L)=0. ; S2X2(:,:,L)=0. ; S2X3(:,:,L)=0.
        S2Y0(:,:,L)=0. ;S2Y1(:,:,L)=0. ; S2Y2(:,:,L)=0. ; S2Y3(:,:,L)=0.
         BXX(:,:,L)=0. ; BYY(:,:,L)=0. ;  BZZ(:,:,L)=0.
      END DO

C**** Calculate all diffusivities
      CALL ISOSLOPE4

      CALL GET_PSI_DIAG

      CALL HALO_UPDATE(grid,AIY0 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,AIY1 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=SOUTH)

C**** Surface layer is calculated directly with appropriate AI,S = 0
      DO L=1,LMO       !Calculating all layers!
      DO J=J_0STG,J_1STG
      IM1 = IM
      DO I=1,IM
      IF(LMM(I,J).lt.L) GO TO 110
C**** FXZPR  Tracer flux at center point in x-direction from tracer
C**** z-gradients by right-side slopes (F=A*S* -gradT)
C**** FXZPR(I,J,L) = AI0(I,J,L)*SIX0(I,J,L)*(TRM(I,J,LP1)/MO(I,J,LP1) -
C****                                        TRM(I,J,L)  /MO(I,J,L))
C****              + AI1(I,J,L)*SIX1(I,J,L)*(TRM(I,J,L)  /MO(I,J,L)   -
C****                                        TRM(I,J,LM1)/MO(I,J,LM1))
C**** FXZPL(I,J,L) = AI2(I,J,L)*SIX2(I,J,L)*
C****                              (TRM(I,J,LP1)/MO(I,J,LP1) -
C****                               TRM(I,J,L)  /MO(I,J,L))
C****                + AI3(I,J,L)*SIX3(I,J,L)*
C****                              (TRM(I,J,L)  /MO(I,J,L)   -
C****                               TRM(I,J,LM1)/MO(I,J,LM1))
C**** FXZPL X-direction tracer flux at center point from z-gradients
C**** by left side slopes.
C**** FXZ   X-direction tracer flux from tracer z-gradients, at center
C**** of right grid wall.
C**** FXZ(I,J,L) = (FXZPR(I,J,L) + FXZPL(IP1,J,L))*1/4 (Four slopes)
C**** TRMXZ+(I,J,L) = TRMXZ-(I,J,L) - (FXZ(I,J,L)-FXZ(IM1,J,L))*DTS/MO
C****
C**** Calculate coefficients of Tracer Operator
C**** ASX0,ASX1,ASX2,ASX3 (ASY0..) A*S for 4 triads in x (y) direction
C**** FXZ   coefficients
      IF (QCROSS) THEN
        DXZ(I,J,L)     = -ASX1(I,J,L) * BYDH(I,J,L)
        CDXZ(IM1,J,L)  = -ASX3(I,J,L) * BYDH(I,J,L)
        BXZ(I,J,L)     = (ASX1(I,J,L) - ASX0(I,J,L))*BYDH(I,J,L)
        CXZ(IM1,J,L)   = (ASX3(I,J,L) - ASX2(I,J,L))*BYDH(I,J,L)
        EXZ(I,J,L)     =  ASX0(I,J,L) * BYDH(I,J,L)
        CEXZ(IM1,J,L)  =  ASX2(I,J,L) * BYDH(I,J,L)
C**** Flux in Y-direction
C**** ASY1(I,JM-1,L), ASY0(I,JM-1,LM1) and ASY3(I,2,L), ASY2(I,2,LM1)??
        DYZ(I,J,L)     = -ASY1(I,J,L) * BYDH(I,J,L)
        BYZ(I,J,L)     = (ASY1(I,J,L) - ASY0(I,J,L))*BYDH(I,J,L)
        CDYZ(I,J-1,L)  = -ASY3(I,J,L) * BYDH(I,J,L)
        CYZ(I,J-1,L)   = (ASY3(I,J,L) - ASY2(I,J,L))*BYDH(I,J,L)
        CEYZ(I,J-1,L)  =  ASY2(I,J,L) * BYDH(I,J,L)
        EYZ(I,J,L)     =  ASY0(I,J,L) * BYDH(I,J,L)
      END IF
C**** Diagonal terms of Horizontal fluxes are decoupled
C**** from downgradients (-ve gradT)  !!!Sign!!!
      BXX(IM1,J,L) = (AIX2(I,J,L) + AIX0(IM1,J  ,L) +
     *                AIX3(I,J,L) + AIX1(IM1,J  ,L))
      BYY(I,J-1,L) = (AIY2(I,J,L) + AIY0(I  ,J-1,L) +
     *                AIY3(I,J,L) + AIY1(I  ,J-1,L))
C**** Z-direction fluxes
C**** Diagonal (Includes AI and DYV/DYP)
      IF (L.gt.KPL(I,J)) THEN
      IF (L.gt.1) BZZ(I,J,L-1) = (S2X1(I,J,L  ) + S2X3(I,J,L  ) +
     *                            S2X0(I,J,L-1) + S2X2(I,J,L-1) +
     *                            S2Y1(I,J,L  ) + S2Y3(I,J,L  ) +
     *                            S2Y0(I,J,L-1) + S2Y2(I,J,L-1))
C**** Off-diagonal
C**** FZXbot(I,J,L) = ASX0(I,J,L)*(TR(I  ,J,L) - TR(IP1,J,L))
C****               + ASX2(I,J,L)*(TR(IM1,J,L) - TR(I  ,J,L))
C**** FZXtop(I,J,L) = ASX1(I,J,L)*(TR(I  ,J,L) - TR(IP1,J,L))
C****               + ASX3(I,J,L)*(TR(IM1,J,L) - TR(I  ,J,L))
C**** Coeficients for (I,J,L) are multiplied by T(IP1,J,L) and T(IM1)
      AZX(I,J,L)   =  ASX2(I,J,L)
      BZX(I,J,L)   =  ASX0(I,J,L) - ASX2(I,J,L)
      CZX(I,J,L)   = -ASX0(I,J,L)
      AEZX(I,J,L)  =  ASX3(I,J,L)
      EZX(I,J,L)   =  ASX1(I,J,L) - ASX3(I,J,L)
      CEZX(I,J,L)  = -ASX1(I,J,L)
C**** y-coefficients *DYV(J-1)/DYV(J-1) or *DYV(J)/DYV(J)
      AZY(I,J,L)   =  ASY2(I,J,L)
      BZY(I,J,L)   =  ASY0(I,J,L) - ASY2(I,J,L)
      CZY(I,J,L)   = -ASY0(I,J,L)
      AEZY(I,J,L)  =  ASY3(I,J,L)
      EZY(I,J,L)   =  ASY1(I,J,L) - ASY3(I,J,L)
      CEZY(I,J,L)  = -ASY1(I,J,L)
      ELSE
      IF (L.gt.1) BZZ(I,J,L-1) = 0.
      AZX(I,J,L)   = 0.
      BZX(I,J,L)   = 0.
      CZX(I,J,L)   = 0.
      AEZX(I,J,L)  = 0.
      EZX(I,J,L)   = 0.
      CEZX(I,J,L)  = 0.
      AZY(I,J,L)   = 0.
      BZY(I,J,L)   = 0.
      CZY(I,J,L)   = 0.
      AEZY(I,J,L)  = 0.
      EZY(I,J,L)   = 0.
      CEZY(I,J,L)  = 0.
      END IF
  110 IM1 = I
      END DO
      END DO
      END DO

!     calculate CDYZ,CYZ,CEYZ,BYY at J_1 from J_1 + 1
      IF (QCROSS) THEN
        CALL HALO_UPDATE(grid,ASY3 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
        CALL HALO_UPDATE(grid,BYDH (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
        CALL HALO_UPDATE(grid,ASY2 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      END IF
      CALL HALO_UPDATE(grid,AIY2 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      CALL HALO_UPDATE(grid,AIY3 (:,grid%j_strt_halo:grid%j_stop_halo,
     *           :) , FROM=NORTH)
      J=J_1STG+1
      if(J.lt.JM) then
      DO L=1,LMO       !Calculating for all layers & all I's
        IM1 = IM
        DO I=1,IM
          IF(LMM(I,J).lt.L) GO TO 111
          IF (QCROSS) THEN
            CDYZ(I,J-1,L)  = -ASY3(I,J,L) * BYDH(I,J,L)
            CYZ(I,J-1,L)   = (ASY3(I,J,L) - ASY2(I,J,L))*BYDH(I,J,L)
            CEYZ(I,J-1,L)  =  ASY2(I,J,L) * BYDH(I,J,L)
          END IF
          BYY(I,J-1,L) = (AIY2(I,J,L) + AIY0(I  ,J-1,L) +
     *                    AIY3(I,J,L) + AIY1(I  ,J-1,L))
  111     IM1 = I
        END DO
      END DO
      endif

C**** End of Main Loop of GMKDIF
      RETURN
      END SUBROUTINE GMKDIF
C****
      SUBROUTINE GMFEXP (TRM, TXM, TYM, TZM, QLIMIT, GIJL)
!@sum  GMFEXP apply GM fluxes to tracer quantities
!@auth Gavin Schmidt/Dan Collins
      USE GM_COM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     *        INTENT(INOUT) :: TRM,TXM,TYM,TZM
      LOGICAL, INTENT(IN) :: QLIMIT
!@var GIJL Tracer flux
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     *         INTENT(INOUT) :: GIJL
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) :: TR
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     *         FXX,FXZ,FYY,FYZ,FZZ,FZX,FZY
      REAL*8 MOFX, MOFY, MOFZ, RFXT, RFYT, RFZT, DT4, DT4DY, DT4DX,
     *     STRNP
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     *         flux_x, flux_y, flux_z
      INTEGER I,J,L,IM1,IP1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H)

      DT4 = 0.25*DTS
C**** Explicit calculation of "fluxes" (R(T))
C**** Calculate tracer concentration
      DO L = 1,LMO
        DO J = J_0,J_1
          DO I = 1,IM
            IF(L.le.LMM(I,J)) THEN
              TR(I,J,L)=TRM(I,J,L)/(DXYPO(J)*MO(I,J,L))
            ELSE
              TR(I,J,L)=0.0
            END IF
            FXX(I,J,L) = 0. ; FXZ(I,J,L) = 0.
            FYY(I,J,L) = 0. ; FYZ(I,J,L) = 0.
            FZZ(I,J,L) = 0. ; FZX(I,J,L) = 0. ; FZY(I,J,L) = 0.
          END DO
        END DO
      END DO

      CALL HALO_UPDATE(grid,TR(:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH+SOUTH)

      DO L=1,LMO
      DO J=J_0S,J_1S
C**** 1/DYPO(DXP) from Y(X) gradient of T for FZY(FZX)
      DT4DY = 0.25*DTS*BYDYP(J)
      DT4DX = 0.25*DTS*BYDXP(J)
      IM1 = IM - 1
      I   = IM
      DO IP1=1,IM
      IF(LMM(I,J).le.0) GO TO 530
C**** Loop for Fluxes in X-direction
      IF(LMU(IM1,J).lt.L) GO TO 510
C**** Calculate fluxes
C**** Diagonal (Includes 1/DXP for gradTX, Not divFX)
      FXX(IM1,J,L) = (DT4DX * BXX(IM1,J,L)) * (TR(IM1,J,L) - TR(I,J,L))
C**** Off-diagonal
      IF (QCROSS) THEN
        IF(L.gt.1) FXZ(IM1,J,L) = FXZ(IM1,J,L) +
     *       DT4 * (DXZ(IM1,J,L) * TR(IM1,J,L-1) +
     *             CDXZ(IM1,J,L) * TR(I,J,L-1))
        FXZ(IM1,J,L) = FXZ(IM1,J,L) +
     *       DT4 * (BXZ(IM1,J,L) * TR(IM1,J,L) +
     *              CXZ(IM1,J,L)  * TR(I,J,L))
C**** Skip for L+1 greater than LMM(IM1,J)
        IF(LMM(IM1,J).gt.L) FXZ(IM1,J,L) =  FXZ(IM1,J,L) +
     *       DT4 * EXZ(IM1,J,L) * TR(IM1,J,L+1)
C**** Skip for L+1 greater than LMM(I,J)
        IF(LMM(I,J).gt.L) FXZ(IM1,J,L) =  FXZ(IM1,J,L) +
     *       DT4 * CEXZ(IM1,J,L) * TR(I,J,L+1)
        !FXZ(IM1,J,L) =  FXZ(IM1,J,L)*(1d0-ARAI)
        FXZ(IM1,J,L) =  -FXZ(IM1,J,L)*(RGMI-1d0)
      END IF
  510 CONTINUE
C**** END of FX
C**** Loop for Fluxes in Y-direction
      IF(LMV(I,J).lt.L) GO TO 520
C**** Calculate fluxes in Y-direction
C**** Diagonal:  Must be divided by DYPO(J) for divF!
      FYY(I,J,L) = DT4*BYY(I,J,L)*(TR(I,J,L) - TR(I,J+1,L))*BYDYV(J)
C**** Off-diagonal:   Divide by DYPO(J) for divF!
C**** YZ terms are divided by appropriate DZV for Z-gradient of T
      IF (QCROSS) THEN
        IF(L.gt.1) FYZ(I,J,L) = FYZ(I,J,L) +
     *       DT4 * (DYZ(I,J,L) * TR(I,J  ,L-1) +
     *             CDYZ(I,J,L) * TR(I,J+1,L-1))
        FYZ(I,J,L) = FYZ(I,J,L) +
     *       DT4 * (BYZ(I,J,L) * TR(I,J  ,L) +
     *              CYZ(I,J,L) * TR(I,J+1,L))
C**** Skip for L+1 greater than LMM(I,J)
        IF(LMM(I,J).gt.L) FYZ(I,J,L) =  FYZ(I,J,L) +
     *       DT4 * EYZ(I,J,L) * TR(I,J,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
        IF(LMM(I,J+1).gt.L) FYZ(I,J,L) =  FYZ(I,J,L) +
     *       DT4 * CEYZ(I,J,L) * TR(I,J+1,L+1)
        !FYZ(I,J,L) =  FYZ(I,J,L) *(1d0-ARAI)
        FYZ(I,J,L) =  -FYZ(I,J,L) *(RGMI-1d0)
      END IF
  520 CONTINUE
C**** END of FY
C**** Loop for Fluxes in Z-direction
      IF(LMM(I,J).le.L) GO TO 530
C**** Calculate fluxes in Z-direction
C**** Diagonal      :  Need BYDH for divF!
      IF(KPL(I,J).le. L) 
     * FZZ(I,J,L) = DT4*BZZ(I,J,L)*(TR(I,J,L+1)-TR(I,J,L))*BYDZV(I,J,L)
C**** Off-diagonal X:  May need to use IM2,IM1,I and F(IM1)
      FZX(I,J,L) = DT4DX * (BZX(I,J,L) * TR(I,J,L) +
     *     AZX(I,J,L) * TR(IM1,J,L) + CZX(I,J,L) * TR(IP1,J,L) +
     *     EZX(I,J,L+1) * TR(I,J,L+1)) !EZX=0 for LMU(IorIM1)=0
C**** Skip for L+1 greater than LMM(I,J)
      IF(LMM(IM1,J).gt.L) FZX(I,J,L) =  FZX(I,J,L) +
     *     DT4DX * AEZX(I,J,L+1) * TR(IM1,J,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
      IF(LMM(IP1,J).gt.L) FZX(I,J,L) =  FZX(I,J,L) +
     *     DT4DX * CEZX(I,J,L+1) * TR(IP1,J,L+1)
C**** Off-diagonal Y:  May need to use JM2,J-1,J and F(J-1)
      FZY(I,J,L) = DT4DY * (BZY(I,J,L) * TR(I,J,L) +
     *     AZY(I,J,L) * TR(I,J-1,L) + CZY(I,J,L) * TR(I,J+1,L) +
     *     EZY(I,J,L+1) * TR(I,J,L+1)) !EZY=0 for LMV(JorJ-1)=0
C**** Skip for L+1 greater than LMM(I,J)
      IF(LMM(I,J-1).gt.L) FZY(I,J,L) =  FZY(I,J,L) +
     *     DT4DY * AEZY(I,J,L+1) * TR(I,J-1,L+1)
C**** Skip for L+1 greater than LMM(I,J+1)
      IF(LMM(I,J+1).gt.L) FZY(I,J,L) =  FZY(I,J,L) +
     *     DT4DY * CEZY(I,J,L+1) * TR(I,J+1,L+1)
  530 IM1 = I
      I = IP1
      END DO
      END DO
      END DO

!     HALO data used in computeFluxes 
      CALL HALO_UPDATE(grid,MO(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid,BYY(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)

      CALL HALO_UPDATE(grid,FYY(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,FYZ(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)

      call computeFluxes(DT4, flux_x, flux_y, flux_z, TXM, TYM,
     &                   TZM, FXX, FXZ, FYY, FYZ, FZZ, FZX, FZY)


      IF (QLIMIT) THEN
        !serial version to preserve original algorithm for TRM update 
        call wrapAdjustFluxes(TRM, flux_x, flux_y, flux_z, GIJL)
      ELSE
        CALL HALO_UPDATE(grid,TRM(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=SOUTH)
      CALL HALO_UPDATE(grid,flux_x(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid,flux_y(:,grid%j_strt_halo:grid%j_stop_halo
     *           ,:) , FROM=NORTH+SOUTH)
        call addFluxes (TRM, flux_x, flux_y, flux_z, GIJL)
      END IF

      RETURN
      END SUBROUTINE GMFEXP
C****
      SUBROUTINE computeFluxes(DT4, flux_x, flux_y, flux_z, TXM, TYM,
     &                         TZM, FXX, FXZ, FYY, FYZ, FZZ, FZX, FZY)
!@sum  computes GM fluxes for tracer quantities
      USE GM_COM, ONLY: grid,GETDomainBounds,IM,JM,LMO,LMM,LMU,LMV,
     &                  DXYPO,MO, kpl,
     &                  BXX, BYY, BZZ, BYDH, RGMI, BYDXP, BYDYP
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT4
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(INOUT) :: flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(INOUT) :: TXM,TYM,TZM
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &             INTENT(IN) :: FXX,FXZ,FYY,FYZ,FZZ,FZX,FZY
      REAL*8 MOFX, MOFY, MOFZ, RFXT, RFYT, RFZT
      INTEGER I,J,L,IM1

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Summation of explicit fluxes to TR, (R(T))
      flux_x = 0.0; flux_y = 0.0; flux_z = 0.0;

      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      IM1 = IM
      DO I=1,IM
      IF(LMM(I,J).le.0) GO TO 610
C**** Loop for Fluxes in X-direction
      IF(LMU(IM1,J).ge.L) THEN
        MOFX = (MO(IM1,J,L) +  MO(I  ,J,L)) * DXYPO(J) *BYDXP(J) *0.5
        RFXT =(FXX(IM1,J,L) + FXZ(IM1,J,L))*MOFX

        flux_x(I,J,L) =RFXT
      endif
C**** Gradient fluxes in X direction affected by diagonal terms
      IF (L.le.LMM(I,J)) TXM(I,J,L) = (TXM(I,J,L) -
     *     3.*(FXX(IM1,J,L)+FXX(I,J,L))*MO(I,J,L)*DXYPO(J)*BYDXP(J))/
     *     (1.+6.*DT4*(BXX(IM1,J,L)+BXX(I,J,L))*BYDXP(J)**2)

C**** Loop for Fluxes in Y-direction
      IF(LMV(I,J-1).ge.L) THEN
        MOFY =((MO(I,J-1,L)*BYDYP(J-1)*DXYPO(J-1)) +
     *         (MO(I,J  ,L)*BYDYP(J  )*DXYPO(J  )))*0.5
        RFYT =(FYY(I,J-1,L) + FYZ(I,J-1,L))*MOFY

        flux_y(I,J,L) =RFYT
      endif
C**** Gradient fluxes in Y direction affected by diagonal terms
      IF (L.le.LMM(I,J)) TYM(I,J,L) = (TYM(I,J,L) -
     *     3.*(FYY(I,J-1,L)+FYY(I,J,L))*MO(I,J,L)*DXYPO(J)*BYDYP(J))/
     *     (1.+6.*DT4*(BYY(I,J-1,L)+BYY(I,J,L))*BYDYP(J)**2)

C**** Summation of explicit fluxes to TR, (R(T)) --- Z-direction
C**** Loop for Fluxes in Z-direction
      IF(LMM(I,J).lt.L) GO TO 610
      IF(LMM(I,J).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
        MOFZ =((MO(I,J,L+1)*BYDH(I,J,L+1)) +
     *         (MO(I,J,L  )*BYDH(I,J,L  ))) * DXYPO(J) *0.5
        !RFZT =(FZZ(I,J,L) +(FZX(I,J,L)+FZY(I,J,L))*(1.d0+ARAI))*MOFZ
        RFZT =(FZZ(I,J,L) +(FZX(I,J,L)+FZY(I,J,L))*(1.d0+RGMI))*MOFZ

        flux_z(I,J,L) =RFZT
      endif
C**** Gradient fluxes in Z direction affected by diagonal terms
c      IF (L.gt.1) THEN
c        TZM(I,J,L) = (TZM(I,J,L) - 3.*(FZZ(I,J,L)+FZZ(I,J,L-1))*
c     *       MO(I,J,L)*DXYPO(J)*BYDH(I,J,L))/(1.+6.*DT4*(BZZ(I,J,L)
c     *       +BZZ(I,J,L-1))*BYDH(I,J,L)**2)
c      ELSE
c        TZM(I,J,L) = (TZM(I,J,L) - 3.*FZZ(I,J,L)*MO(I,J,L)*DXYPO(J)
c     *       *BYDH(I,J,L))/(1.+12.*DT4*BZZ(I,J,L)*BYDH(I,J,L)**2)
c      END IF

C**** END of I and J loops
  610 IM1 = I
      END DO
      END DO

      IF(HAVE_NORTH_POLE) THEN

C****   North Polar box
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,JM-1).ge.L) THEN
            MOFY =((MO(I,JM-1,L)*BYDYP(JM-1)*DXYPO(JM-1)) +
     *             (MO(1,JM  ,L)*BYDYP(JM  )*DXYPO(JM  )))*0.5
            RFYT =(FYY(I,JM-1,L) + FYZ(I,JM-1,L))*MOFY
            flux_y(I,JM,L) = RFYT
          END IF
        END DO

C****   Loop for Fluxes in Z-direction
        IF(LMM(1,JM).ge.L) THEN 
        IF(LMM(1,JM).gt.L) THEN
C****     Calculate new tracer/salinity/enthalpy
          MOFZ =((MO(1,JM,L+1)*BYDH(1,JM,L+1)) +
     *           (MO(1,JM,L  )*BYDH(1,JM,L  ))) * DXYPO(JM) *0.5
          !RFZT =(FZZ(1,JM,L)+FZY(1,JM,L)*(1d0+ARAI))*MOFZ
          RFZT =(FZZ(1,JM,L)+FZY(1,JM,L)*(1d0+RGMI))*MOFZ
          flux_z(1,JM,L) = RFZT
        END IF
C****   Gradient fluxes in Z direction affected by diagonal terms
c        IF (L.gt.1) THEN
c          TZM(1,JM,L) = (TZM(1,JM,L) - 3.*(FZZ(1,JM,L)+FZZ(1,JM,L-1))*
c     *      MO(1,JM,L)*DXYPO(JM)*BYDH(1,JM,L))/(1.+6.*DT4*(BZZ(1,JM,L)
c     *      +BZZ(1,JM,L-1))*BYDH(1,JM,L)**2)
c        ELSE
c          TZM(1,JM,L) = (TZM(1,JM,L) - 3.*FZZ(1,JM,L)*MO(1,JM,L)*
c     &      DXYPO(JM) *BYDH(1,JM,L))/(1.+12.*DT4*BZZ(1,JM,L)*
c     &      BYDH(1,JM,L)**2)
c        END IF
        END IF
      END IF   ! HAVE_NORTH_POLE

 615  IF(HAVE_SOUTH_POLE) THEN

C****   South Polar box
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,2).ge.L) THEN
            MOFY =((MO(I,2,L)*BYDYP(2)*DXYPO(2)) +
     *             (MO(1,1,L)*BYDYP(1)*DXYPO(1)))*0.5
            RFYT =(FYY(I,2,L) + FYZ(I,2,L))*MOFY
            flux_y(I,1,L) = RFYT
          END IF
        END DO

C****   Loop for Fluxes in Z-direction
        IF(LMM(1,1).ge.L) THEN
        IF(LMM(1,1).gt.L) THEN
C****     Calculate new tracer/salinity/enthalpy
          MOFZ =((MO(1,1,L+1)*BYDH(1,1,L+1)) +
     *           (MO(1,1,L  )*BYDH(1,1,L  ))) * DXYPO(1) *0.5
          !RFZT =(FZZ(1,1,L)+FZY(1,1,L)*(1d0+ARAI))*MOFZ
          RFZT =(FZZ(1,1,L)+FZY(1,1,L)*(1d0+RGMI))*MOFZ
          flux_z(1,1,L) = RFZT
        END IF
C****   Gradient fluxes in Z direction affected by diagonal terms
c        IF (L.gt.1) THEN
c          TZM(1,1,L) = (TZM(1,1,L) - 3.*(FZZ(1,1,L)+FZZ(1,1,L-1))*
c     *      MO(1,1,L)*DXYPO(1)*BYDH(1,1,L))/(1.+6.*DT4*(BZZ(1,1,L)
c     *      +BZZ(1,1,L-1))*BYDH(1,1,L)**2)
c        ELSE
c          TZM(1,1,L) = (TZM(1,1,L) - 3.*FZZ(1,1,L)*MO(1,1,L)*
c     &      DXYPO(1) *BYDH(1,1,L))/(1.+12.*DT4*BZZ(1,1,L)*
c     &      BYDH(1,1,L)**2)
c        END IF
        END IF
      END IF   ! HAVE_SOUTH_POLE

 620  END DO   ! L loop

      RETURN
      END SUBROUTINE computeFluxes

      SUBROUTINE wrapAdjustFluxes(TRM, flux_x, flux_y, flux_z, GIJL)
!@sum Applies limiter to GM flux divergence to prevent tracer quantities
!@+   from becoming negative.
      USE GM_COM
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &        INTENT(INOUT) :: TRM, flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     &        INTENT(INOUT) :: GIJL

      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO)
     &        :: conv
      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) ::
     &     convadj_j,convpos_j
      real*8 :: sumadj,sumpos,posadj
      integer :: i,j,l

      INTEGER :: J_0S, J_1S, J_0, J_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid,
     &     J_STRT=J_0, J_STOP=J_1,
     &     J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      call halo_update(grid,flux_y,from=north)

c
c If any flux convergences are adjusted, these flux diagnostics
c will be wrong.  Since the adjustments are typically rare in
c time and space, the diagnostics error will be small.
c Shift x/y indices by 1 for consistency with ocean dynamics
c conventions.
c
      do l=1,lmo
        do j=j_0,j_1
          do i=1,im-1
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(i+1,j,l)
          enddo
          i=im
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(  1,j,l)
          do i=1,im
            gijl(i,j,l,3) = gijl(i,j,l,3) - flux_z(i,j,l) ! gijl is positive down
          enddo
        enddo
        do j=j_0,min(j_1,jm-1)
          do i=1,im
            gijl(i,j,l,2) = gijl(i,j,l,2) + flux_y(i,j+1,l)
          enddo
        enddo
      enddo

c
c Compute flux convergence
c
      do l=lmo,1,-1
        do j=j_0s,j_1s
          do i=1,im-1
            conv(i,j,l) = flux_x(i,j,l)-flux_x(i+1,j,l)
     &                   +flux_y(i,j,l)-flux_y(i,j+1,l)
          enddo
          i=im
            conv(i,j,l) = flux_x(i,j,l)-flux_x(1,j,l)
     &                   +flux_y(i,j,l)-flux_y(i,j+1,l)
        enddo
        if(have_north_pole) conv(1,jm,l) = sum(flux_y(:,jm,l))/im
        if(have_south_pole) conv(1,1,l) = sum(flux_y(:,1,l))/im
        if(l.lt.lmo) then
          do j=j_0s,j_1s
          do i=1,im
            conv(i,j,l  ) = conv(i,j,l  ) + flux_z(i,j,l)
            conv(i,j,l+1) = conv(i,j,l+1) - flux_z(i,j,l)
          enddo
          enddo
          if(have_north_pole) then
            conv(1,jm,l  ) = conv(1,jm,l  ) + flux_z(1,jm,l)
            conv(1,jm,l+1) = conv(1,jm,l+1) - flux_z(1,jm,l)
          endif
          if(have_south_pole) then
            conv(1,1,l  ) = conv(1,1,l  ) + flux_z(1,1,l)
            conv(1,1,l+1) = conv(1,1,l+1) - flux_z(1,1,l)
          endif
        endif
      enddo
      if(have_north_pole) then
        do l=1,lmo
          trm(2:im,jm,l) = trm(1,jm,l)
          conv(2:im,jm,l) = conv(1,jm,l)
        enddo
      endif
      if(have_south_pole) then
        do l=1,lmo
          trm(2:im,1,l) = trm(1,1,l)
          conv(2:im,1,l) = conv(1,1,l)
        enddo
      endif

c
c Reduce the magnitude of flux convergences where they cause
c TRM to become negative.  Keep a tally of the convergence
c adjustments in convadj_j(:)
c
      convadj_j(:) = 0d0
      convpos_j(:) = 0d0
      do l=1,lmo
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).lt.l) cycle
        if(conv(i,j,l).gt.0.) then
          convpos_j(j) = convpos_j(j) +conv(i,j,l)
        elseif(conv(i,j,l).lt.-trm(i,j,l)) then
          convadj_j(j) = convadj_j(j) -trm(i,j,l)-conv(i,j,l)
          conv(i,j,l) = -trm(i,j,l)
        endif
      enddo
      enddo
      enddo

c
c To conserve the global sum of TRM, subtract the adjustments to
c negative convergences from the positive convergences.
c Use a multiplicative factor rather than straight subtraction.
c
      call globalsum(grid,convpos_j,sumpos,all=.true.)
      call globalsum(grid,convadj_j,sumadj,all=.true.)
      if(sumpos.gt.0.) then
        posadj = 1.-sumadj/sumpos
      else
        posadj = 1.
      endif

      if(am_i_root() .and. posadj.ne.1.)
     &     write(6,*) 'posadj ',posadj

c
c Apply the corrected flux convergence
c
      do l=1,lmo
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).lt.l) cycle
        if(conv(i,j,l).gt.0.) conv(i,j,l) = conv(i,j,l)*posadj
        trm(i,j,l) = max(0d0,trm(i,j,l)+conv(i,j,l))
      enddo
      enddo
      enddo

      RETURN
      END SUBROUTINE wrapAdjustFluxes


      SUBROUTINE addFluxes (TRM, flux_x, flux_y, flux_z, GIJL)
!@sum applies GM fluxes to tracer quantities
      USE GM_COM, only: grid, GETDomainBounds, IM, JM, LMO,DXYPO,MO,
     &                        LMM, LMU, LMV, KPL
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &         INTENT(INOUT) :: TRM
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &          INTENT(IN) :: flux_x, flux_y, flux_z
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,3),
     &          INTENT(INOUT) :: GIJL
      REAL*8 STRNP, RFZT, STRSP
      INTEGER :: I, J, L, IM1

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

!updates TRM and GIJL using flux values
      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      IM1 = IM
      DO I=1,IM

      IF(LMM(I,J).le.0) GO TO 613

        IF(LMU(IM1,J).ge.L) THEN
          TRM(I  ,J,L) = TRM(I  ,J,L) + flux_x(I,J,L)
          TRM(IM1,J,L) = TRM(IM1,J,L) - flux_x(I,J,L)
        END IF

        IF(LMV(I,J-1).ge.L) THEN
          TRM(I,J  ,L) = TRM(I,J  ,L) + flux_y(I,J,L)
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
        END IF

C**** END of I and J loops
  613 IM1 = I
      END DO
      END DO

      IF(HAVE_NORTH_POLE) THEN
C****   North Polar box
        STRNP=0.
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,JM-1).ge.L) THEN
C****       Add and Subtract horizontal Y fluxes
            STRNP= STRNP + flux_y(I,JM,L)
            TRM(I,JM-1,L) = TRM(I,JM-1,L) - flux_y(I,JM,L)
          END IF
        END DO
C**** adjust polar box
        TRM(1,JM,L)=TRM(1,JM,L) + STRNP/IM
      END IF

      IF(HAVE_SOUTH_POLE) THEN
C****   South Polar box
        STRSP=0.
C****   Fluxes in Y-direction
        DO I=1,IM
          IF(LMV(I,2).ge.L) THEN
C****       Add and Subtract horizontal Y fluxes
            STRSP= STRSP + flux_y(I,1,L)
            TRM(I,2,L) = TRM(I,2,L) - flux_y(I,1,L)
          END IF
        END DO
C**** adjust polar box
        TRM(1,1,L)=TRM(1,1,L) + STRSP/IM
      END IF
 622  END DO    !L loop


      J = J_1S + 1
      if (J.lt.JM) then

      DO L=1,LMO
C**** Non-Polar boxes
      IM1 = IM
      DO I=1,IM

      IF(LMM(I,J).le.0) GO TO 614

        IF(LMV(I,J-1).ge.L) THEN
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
        END IF

C**** END of L and I loops
  614 IM1 = I
      END DO
      END DO

      endif

C**** Summation of explicit fluxes to TR, (R(T)) --- Z-direction
C**** Extracted from preceding L loop to allow parallelization of that
C**** loop; this loop can't be
      DO L=1,LMO
C**** Non-Polar boxes
      DO J=J_0S,J_1S
      DO I=1,IM
      IF(LMM(I,J).le.0) GO TO 710
C**** Loop for Fluxes in Z-direction
      IF(LMM(I,J).lt.L) GO TO 710
      IF(LMM(I,J).gt.L) THEN
C****   tracer/salinity/enthalpy
        RFZT = flux_z(I,J,L)

C**** Add and Subtract vertical flux. Note +ve upward flux
        TRM(I,J,L  ) = TRM(I,J,L  ) + RFZT
        TRM(I,J,L+1) = TRM(I,J,L+1) - RFZT

      END IF

C**** END of I and J loops
 710  CONTINUE
      END DO
      END DO


      IF(HAVE_NORTH_POLE) THEN
C**** North Polar box
C**** Loop for Fluxes in Z-direction
        IF(LMM(1,JM).lt.L) GO TO 720
        IF(LMM(1,JM).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
          RFZT = flux_z(1,JM,L)
C**** Add and Subtract vertical flux. Note +ve upward flux
          TRM(1,JM,L  ) = TRM(1,JM,L  ) + RFZT
          TRM(1,JM,L+1) = TRM(1,JM,L+1) - RFZT
        END IF
 720    CONTINUE
      END IF     ! HAVE_NORTH_POLE


      IF(HAVE_SOUTH_POLE) THEN
C**** South Polar box
C**** Loop for Fluxes in Z-direction
        IF(LMM(1,1).lt.L) GO TO 730
        IF(LMM(1,1).gt.L) THEN
C**** Calculate new tracer/salinity/enthalpy
          RFZT = flux_z(1,1,L)
C**** Add and Subtract vertical flux. Note +ve upward flux
          TRM(1,1,L  ) = TRM(1,1,L  ) + RFZT
          TRM(1,1,L+1) = TRM(1,1,L+1) - RFZT
        END IF
 730    CONTINUE
      END IF     ! HAVE_SOUTH_POLE
      END DO

c
c Save flux diagnostics.
c Shift x/y indices by 1 for consistency with ocean dynamics
c conventions.
c
      do l=1,lmo
        do j=j_0,j_1
          do i=1,im-1
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(i+1,j,l)
          enddo
          i=im
            gijl(i,j,l,1) = gijl(i,j,l,1) + flux_x(  1,j,l)
          do i=1,im
            gijl(i,j,l,3) = gijl(i,j,l,3) - flux_z(i,j,l) ! gijl is positive down
          enddo
        enddo
        do j=j_0,min(j_1,jm-1)
          do i=1,im
            gijl(i,j,l,2) = gijl(i,j,l,2) + flux_y(i,j+1,l)
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE addFluxes
C****
      SUBROUTINE ISOSLOPE4
!@sum  ISOSLOPE4 calculates the isopycnal slopes from density triads
!@auth Gavin Schmidt/Dan Collins
!@ver  2009/02/13
      USE GM_COM
      IMPLICIT NONE
      INTEGER I,J,L,IM1
      REAL*8 :: AIX0ST,AIX2ST,AIY0ST,AIY2ST,SIX0,SIX2,SIY0,SIY2,
     *          AIX1ST,AIX3ST,AIY1ST,AIY3ST,SIX1,SIX3,SIY1,SIY3
      Real*8 :: byAIDT, DSX0sq,DSX1sq,DSX2sq,DSX3sq,
     *                  DSY0sq,DSY1sq,DSY2sq,DSY3sq
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      LOGICAL :: dothis

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H)

C**** Calculate horizontal and vertical density gradients.
      CALL DENSGRAD
C****
C**** Three grid boxes (triads) are used for each slope. The GM skew
C**** diffusion coefficient is calculated for each triad as well.
C**** The diffusion coefficient is taken from Visbeck et al (1997)
C****
C**** Main Loop over I,J and L
      DO L=1,LMO
      DO J=J_0STG,J_1STG
      IM1=IM
      DO I=1,IM
      IF(LMM(I,J).lt.L) GO TO 800
C**** SIX0, SIY0, SIX2, SIY2: four slopes that use RHOMZ(L)
      IF(L.EQ.LMM(I,J) .or. RHOMZ(I,J,L).eq.0.) THEN
        AIX0ST = 0.
        AIX2ST = 0.
        AIY0ST = 0.
        AIY2ST = 0.
        SIX0 = 0.
        SIX2 = 0.
        SIY0 = 0.
        SIY2 = 0.
        AIX0(I,J,L) = 0.
        AIX2(I,J,L) = 0.
        AIY0(I,J,L) = 0.
        AIY2(I,J,L) = 0.
      ELSE
        AIX0ST = ARIV(I,J,L)
        AIX2ST = ARIV(I,J,L)
        AIY0ST = ARIV(I,J,L)
        AIY2ST = ARIV(I,J,L)
        SIX0 = RHOX(I  ,J,L) * BYRHOZ(I,J,L)
        SIX2 = RHOX(IM1,J,L) * BYRHOZ(I,J,L)
        SIY2 = RHOY(I,J-1,L) * BYRHOZ(I,J,L)
        SIY0 = RHOY(I,J  ,L) * BYRHOZ(I,J,L)
        IF (ARIV(I,J,L).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J,L)+ARIV(I,J,L)))
          DSX0sq = DZV(I,J,L)**2 * byAIDT
          DSX2sq = DZV(I,J,L)**2 * byAIDT
          DSY0sq = DZV(I,J,L)**2 * byAIDT
          DSY2sq = DZV(I,J,L)**2 * byAIDT
          If (SIX0**2 > DSX0sq)  AIX0ST = AIX0ST * DSX0sq / SIX0**2
          If (SIX2**2 > DSX2sq)  AIX2ST = AIX2ST * DSX2sq / SIX2**2
          If (SIY0**2 > DSY0sq)  AIY0ST = AIY0ST * DSY0sq / SIY0**2
          If (SIY2**2 > DSY2sq)  AIY2ST = AIY2ST * DSY2sq / SIY2**2
        END IF
C**** AI are always * layer thickness for vertical gradient in FXX, FYY
        AIX0(I,J,L) = AIX0ST * DZV(I,J,L) * BYDH(I,J,L)
        AIX2(I,J,L) = AIX2ST * DZV(I,J,L) * BYDH(I,J,L)
        AIY0(I,J,L) = AIY0ST * DZV(I,J,L) * BYDH(I,J,L)
        AIY2(I,J,L) = AIY2ST * DZV(I,J,L) * BYDH(I,J,L)
      ENDIF
C**** SIX1, SIY1, SIX3, SIY3: four slopes that use RHOMZ(L-1)
      IF(L.EQ.1.) THEN
        AIX1ST = 0. ; AIX3ST = 0. ; AIY1ST = 0. ;  AIY3ST = 0.
        SIX1 = 0.   ; SIX3 = 0.   ; SIY1 = 0.   ;  SIY3 = 0.
        AIX1(I,J,L) = 0. ; AIX3(I,J,L) = 0.
        AIY1(I,J,L) = 0. ; AIY3(I,J,L) = 0.
      ELSEIF (RHOMZ(I,J,L-1).eq.0.) THEN
        AIX1ST = 0. ; AIX3ST = 0. ; AIY1ST = 0. ;  AIY3ST = 0.
        SIX1 = 0.   ; SIX3 = 0.   ; SIY1 = 0.   ;  SIY3 = 0.
        AIX1(I,J,L) = 0. ; AIX3(I,J,L) = 0.
        AIY1(I,J,L) = 0. ; AIY3(I,J,L) = 0.
      ELSE
        AIX1ST = ARIV(I,J,L)
        AIX3ST = ARIV(I,J,L)
        AIY1ST = ARIV(I,J,L)
        AIY3ST = ARIV(I,J,L)
        SIX1 = RHOX(I  ,J,L) * BYRHOZ(I,J,L-1)
        SIX3 = RHOX(IM1,J,L) * BYRHOZ(I,J,L-1)
        SIY1 = RHOY(I,J  ,L) * BYRHOZ(I,J,L-1)
        SIY3 = RHOY(I,J-1,L) * BYRHOZ(I,J,L-1)
        IF (ARIV(I,J,L).gt.0.) THEN ! limit slopes <ML
          byAIDT = 1 / (4*DTS*(AINV(I,J,L)+ARIV(I,J,L)))
          DSX1sq = DZV(I,J,L-1)**2 * byAIDT
          DSX3sq = DZV(I,J,L-1)**2 * byAIDT
          DSY1sq = DZV(I,J,L-1)**2 * byAIDT
          DSY3sq = DZV(I,J,L-1)**2 * byAIDT
          If (SIX1**2 > DSX1sq)  AIX1ST = AIX1ST * DSX1sq / SIX1**2
          If (SIX3**2 > DSX3sq)  AIX3ST = AIX3ST * DSX3sq / SIX3**2
          If (SIY1**2 > DSY1sq)  AIY1ST = AIY1ST * DSY1sq / SIY1**2
          If (SIY3**2 > DSY3sq)  AIY3ST = AIY3ST * DSY3sq / SIY3**2
        END IF
C**** AI are always * layer thickness for vertical gradient in FXX, FYY
        AIX1(I,J,L) = AIX1ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIX3(I,J,L) = AIX3ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIY1(I,J,L) = AIY1ST * DZV(I,J,L-1) * BYDH(I,J,L)
        AIY3(I,J,L) = AIY3ST * DZV(I,J,L-1) * BYDH(I,J,L)
      ENDIF
C**** AIX0...AIX3, AIY0...AIY3
      ASX0(I,J,L) = AIX0ST * SIX0
      ASX1(I,J,L) = AIX1ST * SIX1
      ASX2(I,J,L) = AIX2ST * SIX2
      ASX3(I,J,L) = AIX3ST * SIX3
      ASY0(I,J,L) = AIY0ST * SIY0
      ASY1(I,J,L) = AIY1ST * SIY1
      ASY2(I,J,L) = AIY2ST * SIY2
      ASY3(I,J,L) = AIY3ST * SIY3
C**** S2X0...S2X3, S2Y0...S2Y3
      S2X0(I,J,L) = AIX0ST * SIX0 * SIX0
      S2X1(I,J,L) = AIX1ST * SIX1 * SIX1
      S2X2(I,J,L) = AIX2ST * SIX2 * SIX2
      S2X3(I,J,L) = AIX3ST * SIX3 * SIX3
      S2Y0(I,J,L) = AIY0ST * SIY0 * SIY0 * BYDYP(J) * DYVO(J)
      S2Y1(I,J,L) = AIY1ST * SIY1 * SIY1 * BYDYP(J) * DYVO(J)
      S2Y2(I,J,L) = AIY2ST * SIY2 * SIY2 * BYDYP(J) * DYVO(J-1)
      S2Y3(I,J,L) = AIY3ST * SIY3 * SIY3 * BYDYP(J) * DYVO(J-1)
#ifdef OCN_GISS_MESO
      AIX0(I,J,L) = ARIV(I,J,L)
      AIX2(I,J,L) = ARIV(I,J,L)
      AIY0(I,J,L) = ARIV(I,J,L)
      AIY2(I,J,L) = ARIV(I,J,L)
      AIX1(I,J,L) = ARIV(I,J,L)
      AIX3(I,J,L) = ARIV(I,J,L)
      AIY1(I,J,L) = ARIV(I,J,L)
      AIY3(I,J,L) = ARIV(I,J,L)
#endif
  800 IM1 = I
      END DO
      END DO
      END DO
      RETURN
      END SUBROUTINE ISOSLOPE4
C****
      SUBROUTINE DENSGRAD
!@sum  DENSGRAD calculates all horizontal and vertical density gradients
!@auth Gavin Schmidt/Dan Collins
      USE GM_COM
      USE FILEMANAGER

      use ocean, only : dzo

      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     &     RHO
      REAL*8  BYRHO,DZVLM1,CORI,BETA,ARHO,ARHOX,ARHOY,ARHOZ,AN,RD
     *     ,BYTEADY,DH0,DZSUMX,DZSUMY,R1,R2,P12
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO) ::
     *   kappam3d
      REAL*8, SAVE :: HUP
      INTEGER I,J,L,IM1,LAV,iu_ODIFF
      INTEGER, SAVE :: IFIRST = 1
      CHARACTER TITLE*80
      Real*8,External   :: VOLGSP
      REAL*8, DIMENSION(IM,JM) ::  AINV_glob

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H, J_1H
      INTEGER :: J_1HR
      LOGICAL :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE 
      logical :: dothis
      real*8 :: k1d(lmo)

#if (defined CONSTANT_MESO_DIFFUSIVITY)
#define USE_1D_DIFFUSIVITY
      k1d(:) = meso_diffusivity_const
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
#define USE_1D_DIFFUSIVITY
      do l=1,lmo
        k1d(l) = meso_diffusivity_z0*
     &       exp(-(.5*(ze(l-1)+ze(l)))/meso_diffusivity_zscale)
     &       + meso_diffusivity_deep
      enddo

#endif

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** set up geometry needed
      IF (IFIRST.eq.1) THEN
        DO J=J_0,J_1
          BYDYP(J)=1d0/DYPO(J)
          BYDXP(J)=1d0/DXPO(J)
        END DO
        !DO J=1,JM-1
        DO J=J_0,J_1S
          BYDYV(J)=1d0/DYVO(J)
        END DO
C**** Calculate level at 1km depth
        LUP=0
   10   LUP=LUP + 1
        IF (ZE(LUP+1).lt.1d3) GOTO 10
        HUP = ZE(LUP)
        IFIRST = 0. 

        CALL HALO_UPDATE(grid,BYDYP (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)
        CALL HALO_UPDATE(grid,BYDYV (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)
        CALL HALO_UPDATE(grid,BYDXP (grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH+NORTH)

      END IF
C****
      call vbar_gm0

      !initialize
      RHO = -0.0; BYDH = -0.0;
      DZV = -0.0; BYDZV = -0.0;

      DO L=1,LMO
      DO J=J_0,J_1
      DO I=1,IM
C**** Skip non-ocean grid points
        IF (L.le.LMM(I,J)) THEN
C**** RHO(I,J,L)  Density=1/specific volume
          BYRHO = VBAR(I,J,L)
            DH0 =   DH(I,J,L)
          IF (J.eq.1)  THEN         ! South pole
            BYRHO = VBAR(1,1,L)
              DH0 =   DH(1,1,L)
          END IF
          IF (J.eq.JM) THEN         ! North pole
            BYRHO = VBAR(1,JM,L)
              DH0 =   DH(1,JM,L)
          END IF
          RHO(I,J,L)  = 1d0/BYRHO
          BYDH(I,J,L) = 1d0/DH0
          IF(L.gt.1) THEN
            DZVLM1     = 0.5* (DH(I,J,L) + DH(I,J,L-1))
            DZV(I,J,L-1) = DZVLM1
            BYDZV(I,J,L-1) = 1d0/DZVLM1
          END IF
        END IF
 911  END DO
      END DO
      END DO

      CALL HALO_UPDATE(grid,RHO (:,grid%j_strt_halo:grid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)

      CALL HALO_UPDATE(grid,G3D)
      CALL HALO_UPDATE(grid,S3D)
      CALL HALO_UPDATE(grid,P3D)

C**** Calculate density gradients

      !initialize
      RHOMZ = -0.0; BYRHOZ = -0.0;
      RHOY = -0.0; RHOX = -0.0;

      J_1HR = min(J_1STG+1, JM)
      DO L=1,LMO
        DO J=J_0STG,J_1HR
          IM1 = IM
          DO I=1,IM
C**** Skip non-ocean grid points
            IF(LMM(I,J).lt.L) GO TO 931
C**** minus vertical gradient
            IF(L.gt.1) THEN
              RHOMZ(I,J,L-1)=MAX(0d0,
     *           -dVBARdZ(I,J,L)*BYDZV(I,J,L-1)/VBAR(I,J,L-1)**2)
              IF(RHOMZ(I,J,L-1).ne.0.)
     *             BYRHOZ(I,J,L-1)=1./RHOMZ(I,J,L-1)
            END IF
C**** Calculate horizontal gradients
            IF(LMV(I,J-1).ge.L) THEN
              p12 = .5d0*(p3d(i,j-1,l)+p3d(i,j,l))
              r1 = 1d0/volgsp(g3d(i,j-1,l),s3d(i,j-1,l),p12)
              r2 = 1d0/volgsp(g3d(i,j  ,l),s3d(i,j  ,l),p12)
              RHOY(I,J-1,L) =
C     *           (RHO(I,J,L) - RHO(I,J-1,L))*BYDYV(J-1)
     *           (r2 - r1)*BYDYV(J-1)
            ENDIF
            IF(LMU(IM1,J).ge.L) THEN
              p12 = .5d0*(p3d(im1,j,l)+p3d(i,j,l))
              r1 = 1d0/volgsp(g3d(im1,j,l),s3d(im1,j,l),p12)
              r2 = 1d0/volgsp(g3d(i  ,j,l),s3d(i  ,j,l),p12)
              RHOX(IM1,J,L) =
C     *           (RHO(I,J,L) - RHO(IM1,J,L))*BYDXP(J)
     *           (r2 - r1)*BYDXP(J)
            ENDIF
  931       IM1 = I
          END DO
        END DO
      END DO
#ifdef OCN_GISS_MESO 
      CALL OCN_mesosc(kappam3d)
#endif
C**** Calculate VMHS diffusion = amu* min(NH/f,equ.rad)^2 /Teady
      AINV = 0.
      ARIV = 0.
      DO J=J_0S,J_1S
        CORI = ABS(2d0*OMEGA*SINPO(J))
        BETA = ABS(2d0*OMEGA*COSPO(J)/RADIUS)
        IM1=IM
        DO I=1,IM
          IF (LMM(I,J).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
            ARHO  = 0.
            ARHOX = 0.
            ARHOY = 0.
            DZSUMX = 0.
            DZSUMY = 0.
            LAV = MIN(LUP,LMM(I,J))
            DO L=1,LAV
              ARHO  = ARHO  + RHO(I,J,L)
              IF(LMU(IM1,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(IM1,J,L)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMU(I  ,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(I  ,J,L)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMV(I,J-1).ge.L) THEN
                ARHOY = ARHOY + RHOY(I,J-1,L)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
              IF(LMV(I,J  ).ge.L) THEN
                ARHOY = ARHOY + RHOY(I,J  ,L)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
            ARHO  = ARHO / REAL(LAV,KIND=8)
            IF (DZSUMX.gt.0.) ARHOX = ARHOX / DZSUMX
            IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
            IF (LAV.gt.1) THEN
              ARHOZ=2*RHOZ1K(I,J)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
            ELSE
              ARHOZ = 0.
            END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
            IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
              ARIV(I,J,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
              ARIV(I,J,:) = kappam3d(I,J,:)
#else
              AN = SQRT(GRAV * ARHOZ / ARHO)
#ifdef CONSTANT_MESO_LENSCALE
              RD = meso_lenscale_const
#else
              RD = AN * HUP / CORI
              IF (RD.gt.ABS(J-.5*(JM+1))*DYPO(J)) RD=SQRT(AN*HUP/BETA)
#endif
              BYTEADY = GRAV * SQRT(ARHOX*ARHOX + ARHOY*ARHOY) / (AN
     *             *ARHO)
              ARIV(I,J,:) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
            END IF
            !ARIV(I,J,:) = ARAI * AINV(I,J,:) ! was = ARI
            AINV(I,J,:) = RGMI * ARIV(I,J,:)
          END IF
          IM1=I
C**** Set diagnostics
          OIJ(I,J,IJ_GMSC) = OIJ(I,J,IJ_GMSC) + AINV(I,J,1) ! GM-scaling
        END DO
      END DO
C**** North pole
      if ( HAVE_NORTH_POLE ) then
        IF (LMM(1,JM).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,JM))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(1,JM,L)
            DO I=1,IM
              IF(LMV(I,JM-1).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(I,JM-1,L))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
            ARHOZ=2*RHOZ1K(1,JM)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
            ARIV(1,JM,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
            ARIV(1,JM,:) = kappam3d(1,JM,:)
#else
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            ARIV(1,JM,:) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
          END IF
          !ARIV(1,JM,:) = ARAI * AINV(1,JM,:) ! was = ARI
          AINV(1,JM,:) = RGMI * ARIV(1,JM,:)
        END IF
c       AINV(2:IM,JM)=AINV(1,JM)
c       ARIV(2:IM,JM)=ARIV(1,JM)
        DO L=1,LMO
          AINV(2:IM,JM,L)=AINV(1,JM,L)
          ARIV(2:IM,JM,L)=ARIV(1,JM,L)
        ENDDO
C**** Set diagnostics
        OIJ(1,JM,IJ_GMSC) = OIJ(1,JM,IJ_GMSC) + AINV(1,JM,1) ! GM-scaling
      endif
C**** South pole
      if ( HAVE_SOUTH_POLE ) then
        IF (LMM(1,1).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,1))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(1,1,L)
            DO I=1,IM
              IF(LMV(I,2).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(I,2,L))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
            ARHOZ=2*RHOZ1K(1,1)/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
#ifdef USE_1D_DIFFUSIVITY
            ARIV(1,1,:) = k1d(:)
#else
#ifdef OCN_GISS_MESO
            ARIV(1,1,:) = kappam3d(1,1,:)
#else
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            ARIV(1,1,:) = AMU * RD**2 * BYTEADY ! was = AIN
#endif
#endif
          END IF
          !ARIV(1,1,:) = ARAI * AINV(1,1,:) ! was = ARI
          AINV(1,1,:) = RGMI * ARIV(1,1,:)
        END IF
c       AINV(2:IM,1)=AINV(1,1)
c       ARIV(2:IM,1)=ARIV(1,1)
        DO L=1,LMO
          AINV(2:IM,1,L)=AINV(1,1,L)
          ARIV(2:IM,1,L)=ARIV(1,1,L)
        ENDDO
C**** Set diagnostics
        OIJ(1,1,IJ_GMSC) = OIJ(1,1,IJ_GMSC) + AINV(1,1,1) ! GM-scaling
      endif
C****
c      IF (IFIRST.eq.1) THEN  !output GM diffusion coefficient
c        CALL PACK_DATA(grid,    AINV  ,    AINV_glob)
c        if( AM_I_ROOT() ) then
c          call openunit('ODIFF',iu_ODIFF,.true.,.false.)
c          TITLE = "Visbeck scaling for GM coefficient m^2/s"
c          WRITE(iu_ODIFF) TITLE,((REAL(AINV_glob(I,J),KIND=4),i=1,im),
c     *                            j=1,jm)
c          call closeunit(iu_ODIFF)
c        endif
c        IFIRST = 0
c      END IF

      RETURN
C****
      END SUBROUTINE DENSGRAD

      subroutine get_psi_diag
c
c Calculate bolus velocity diagnostics and store them in units of
c accumulated mass flux (kg) for consistency with the diagnostics
c for resolved velocities. For convenience, a constant reference
c density is used to convert the eddy-induced streamfunction to
c mass-flux units.  Local densities can be used if the O(1%) discrepancy
c proves to be a confounding factor in analyses of this diagnostic.
c
c Since the bolus velocity does not appear explicitly in the skew-flux
c form of GM, the attempt has not been made here to exactly reproduce
c the effective interpolation stencil for the diffusivities and density
c gradients used in ISOSLOPE4.  However, a future commit will attempt
c to insert code into ISOSLOPE4 to store the diffusivities in effect
c after modification by CFL criteria and ML exclusion (rather than
c mimicking that logic here).

      use gm_com, only : ainv
      use odiag, only : oijl=>oijl_loc,ijl_mfub,ijl_mfvb
      use ocean, only : dzo,dxvo,dypo
      use ocean, only : im,jm,lmo,dts
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use constant, only : rhows
      use domain_decomp_1d, only : getdomainbounds, halo_update, north
      use oceanr_dim, only : grid=>ogrid
      use kpp_com, only : kpl
      use gm_com, only : rhox,rhoy,rhomz,byrhoz
      implicit none
      integer :: i,j,l,n,ip1
      integer :: j_0,j_1,j_0s,j_1s
      logical :: have_north_pole
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     kbyrhoz
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     psi
      real*8 :: wtup,wtdn,rhomid,kmax,rhox_,slx,rhoy_,sly,krat

c**** Extract domain decomposition info
      call getdomainbounds(grid, j_strt = j_0, j_stop = j_1,
     &     j_strt_skp  = j_0s,   j_stop_skp  = j_1s,
     &               have_north_pole = have_north_pole)

      kbyrhoz = 0.
      do l=1,lmo-1
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        do j=j_0,j_1
        do n=1,nbyzm(j,l+1)
        do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
          if(l.gt.kpl(i,j)) then ! GM excludes ML currently
            ! interpolate K/(drho/dz) to layer edges
            kbyrhoz(i,j,l) = (wtup*ainv(i,j,l)+wtdn*ainv(i,j,l+1))*
     &           byrhoz(i,j,l)
          endif
        enddo
        enddo
        enddo
      enddo
      call halo_update(grid,kbyrhoz,from=north)
      call halo_update(grid,byrhoz,from=north) ! needed?

      rhomid = rhows ! for now

      psi = 0.
      do l=1,lmo-1
        kmax = ( .25*(dzo(l)+dzo(l+1))**2 )/dts
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        do j=j_0s,j_1s
        do n=1,nbyzu(j,l+1)
        do i=i1yzu(n,j,l+1),i2yzu(n,j,l+1)
          ip1 = i+1; if(i.eq.im) ip1 = 1
          ! vertically interpolate rhox to layer edges
          rhox_ = (wtup*rhox(i,j,l)+wtdn*rhox(i,j,l+1))
          ! calculate x-slope
          slx = rhox_*(byrhoz(i,j,l)+byrhoz(ip1,j,l))*.5
          ! streamfunc = K
          psi(i,j,l) = rhox_ *(kbyrhoz(i,j,l)+kbyrhoz(ip1,j,l))*.5
          ! mimic the "CFL tapering" in ISOSLOPE4: locally reduce
          ! K so that the vertical component of Redi diffusion does
          ! not violate CFL conditions for numerical stability
          krat = psi(i,j,l)*slx/kmax
          if(krat .gt. 1d0) psi(i,j,l) = psi(i,j,l)/krat
          psi(i,j,l) = psi(i,j,l)*rhomid
        enddo
        enddo
        enddo
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do n=1,nbyzu(j,l)
        do i=i1yzu(n,j,l),i2yzu(n,j,l)
          oijl(i,j,l,ijl_mfub) = oijl(i,j,l,ijl_mfub) +
     &         (psi(i,j,l)-psi(i,j,l-1))*dypo(j)*dts
        enddo
        enddo
        enddo
      enddo

      psi = 0.
      do l=1,lmo-1
        kmax = ( .25*(dzo(l)+dzo(l+1))**2 )/dts
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        do j=j_0s,j_1s
        do n=1,nbyzv(j,l+1)
        do i=i1yzv(n,j,l+1),i2yzv(n,j,l+1)
          ! vertically interpolate rhoy to layer edges
          rhoy_ = (wtup*rhoy(i,j,l)+wtdn*rhoy(i,j,l+1))
          ! calculate y-slope
          sly = rhoy_*(byrhoz(i,j,l)+byrhoz(i,j+1,l))*.5
          ! streamfunc = K
          psi(i,j,l) = rhoy_ *(kbyrhoz(i,j,l)+kbyrhoz(i,j+1,l))*.5
          ! mimic the "CFL tapering" in ISOSLOPE4: locally reduce
          ! K so that the vertical component of Redi diffusion does
          ! not violate CFL conditions for numerical stability
          krat = psi(i,j,l)*sly/kmax
          if(krat .gt. 1d0) psi(i,j,l) = psi(i,j,l)/krat
          psi(i,j,l) = psi(i,j,l)*rhomid
        enddo
        enddo
        enddo
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do n=1,nbyzv(j,l)
        do i=i1yzv(n,j,l),i2yzv(n,j,l)
          oijl(i,j,l,ijl_mfvb) = oijl(i,j,l,ijl_mfvb) +
     &         (psi(i,j,l)-psi(i,j,l-1))*dxvo(j)*dts
        enddo
        enddo
        enddo
      enddo

      return
      end subroutine get_psi_diag

      Subroutine VBAR_GM0
!@sum VBAR_GM0 calculates specific volume and vertical gradients for GM
!@auth Gary Russell/Gavin Schmidt
      Use CONSTANT, only: grav
      Use GM_COM, only: vbar,dvbardz,rhoz1k,lup,g3d,s3d,p3d
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,DXYPO, MO, ZE,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO, OPRESS, FOCEAN
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L,IMAX, J1,JN,JNH, LAV,LAVM
      Real*8, dimension(lmo) :: gup,gdn,sup,sdn,pm,gmd,smd
      Real*8 :: vup,vdn,vupu,vdnu,pe,bym
      Logical*4 QSP,QNP
      Real*8,External   :: VOLGSP,temgs

C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID,OPRESS,FROM=NORTH)
      Call HALO_UPDATE (GRID,  MO, FROM=NORTH)
      Call HALO_UPDATE (GRID, G0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, GZM, FROM=NORTH)
      Call HALO_UPDATE (GRID, S0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, SZM, FROM=NORTH)

      Do J=J1,JNH
        IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
        Do I=1,IMAX
          If (FOCEAN(I,J) == 0)  CYCLE

          Do L=1,LMOM(I,J)
            BYM=1d0/(MO(I,J,L)*DXYPO(J))
            GMD(L)= G0M(I,J,L)*BYM
            GUP(L)=(G0M(I,J,L)-2*z12eH*GZM(I,J,L))*BYM
            GDN(L)=(G0M(I,J,L)+2*z12eH*GZM(I,J,L))*BYM
            SMD(L)= S0M(I,J,L)*BYM
            SUP(L)=(S0M(I,J,L)-2*z12eH*SZM(I,J,L))*BYM
            SDN(L)=(S0M(I,J,L)+2*z12eH*SZM(I,J,L))*BYM
          EndDo
          
C**** Calculate pressure by integrating from the top down
          PE = OPRESS(I,J)
          Do L=1,LMOM(I,J)
            PM(L) = PE
            PE    = PE + MO(I,J,L)*GRAV
            PM(L) = .5d0*(PE + PM(L))
            G3D(I,J,L) = GMD(L)
            S3D(I,J,L) = SMD(L)
            P3D(I,J,L) = PM(L)
          EndDo
C**** Calculate potential specific volume (ref to mid-point pr)
          Do L=LMOM(I,J),1,-1
            VUP = VOLGSP (GUP(L),SUP(L),PM(L))
            VDN = VOLGSP (GDN(L),SDN(L),PM(L))
            VBAR(I,J,L) = (VUP + VDN)*.5
C**** Vertical gradient calculated using lower box mid-point pr
            IF (L.gt.1) then 
              VUPU = VOLGSP (GUP(L-1),SUP(L-1),PM(L))
              VDNU = VOLGSP (GDN(L-1),SDN(L-1),PM(L))
              dVBARdZ(I,J,L) = .5* (VUP + VDN - VUPU - VDNU)
            end if
          EndDo
C**** Vertical potential gradient in top 1km
          LAV = MIN(LUP,LMOM(I,J))
          LAVM = MAX(LAV/2,1)   ! mid depth
          VUP = VOLGSP (GMD(1),SMD(1),PM(LAVM))
          VDN = VOLGSP (GMD(LAV),SMD(LAV),PM(LAVM))
          RHOZ1K(I,J) = (VUP - VDN)/VBAR(I,J,LAVM)**2
        EndDo
      EndDo
C**** Copy VBAR to all longitudes at poles
      If (QNP) Then
        Do L=1,LMOM(1,JM)
          VBAR(2:IM,JM,L) = VBAR(1,JM,L)
          G3D(2:IM,JM,L) = G3D(1,JM,L)
          S3D(2:IM,JM,L) = S3D(1,JM,L)
          P3D(2:IM,JM,L) = P3D(1,JM,L)
          dVBARdZ(2:IM,JM,L) = dVBARdZ(1,JM,L)
        EndDo
      EndIf
      If (QSP) Then
        Do L=1,LMOM(1,1)
          VBAR(2:IM,1,L) = VBAR(1,1,L)
          G3D(2:IM,1,L) = G3D(1,1,L)
          S3D(2:IM,1,L) = S3D(1,1,L)
          P3D(2:IM,1,L) = P3D(1,1,L)
          dVBARdZ(2:IM,1,L) = dVBARdZ(1,1,L)
        EndDo
      EndIf

      Return
      End Subroutine VBAR_GM0
