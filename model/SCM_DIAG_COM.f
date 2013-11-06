!     SCM_DIAG_COM.f    
!     save diagnostics for run of MODELE SCM

      MODULE SCMDIAG  

      USE RESOLUTION , only : LM 
      USE DIAG_COM , only : npres,ntau
      USE CLOUDS_COM, only : ncol

      IMPLICIT NONE


      real*8 LHPSAV(LM),LHPMC(LM),PRESAV(LM),PREMC(LM)

      real*8   CLCVSS(LM),CLCVMC(LM),CLTHCK(LM),CSIZE(LM,2),   
     *         EFFRAD(LM),TAUSSC(LM),TAUMCC(LM),CUMFLX(LM,2,LM),
     *         CUMHET(LM),CUMOST(LM),PRCSS,PRCMC,EVPFLX,SHFLX,
     *         SOILMS,CLDFLG(LM),DWNFLX(LM,2,LM),RHC(LM)
      real*8   clsav(LM)
      real*8 SRDFLBTOP,SRNFLBTOP,SRUFLBTOP,SRUFLBBOT,
     *       TRUFLBTOP,TRDFLBTOP,SRDFLBBOT,SRNFLBBOT,
     *       TRUFLBBOT,TRDFLBBOT,CSSRNTOP,CSTRUTOP,
     *       CSSRNBOT,CSTRNBOT,CSSRDBOT,TRNFLBBOT
      real*8 SCM_PBL_HGT
      real*8, DIMENSION(LM) :: SRFHRLCOL,TRFCRLCOL
c     
c     isccp diagnostics
c
      real*8 isccp_sunlit,isccp_ctp,isccp_tauopt,isccp_lowcld,
     &       isccp_midcld,isccp_highcld,isccp_fq(ntau,npres),
     &       isccp_totcldarea,isccp_boxtau(ncol),isccp_boxptop(ncol)
c
c     cumulus updraft speed
      real*8   WCUSCM(LM,2),WCUALL(LM,2,LM)
      real*8   WCUDEEP(LM,2)
C--- Added by J.W. starting ---C
      real*8   MPLUMESCM(LM,2),MPLUMEALL(LM,2,LM)
      real*8   MPLUMEDEEP(LM,2)
      real*8   ENTSCM(LM,2),ENTALL(LM,2,LM)
      real*8   ENTDEEP(LM,2)
      real*8   DETRAINDEEP(LM,2,LM)
C--- Added by J.W. ending ---C
c
c     precipitating and non-precipitating convective condensate for
c     deep convection 
      real*8 PRCCDEEP(LM,2,LM),NPRCCDEEP(LM,2,LM),TPALL(LM,2,LM)   
      real*8 PRCCGRP(LM,2,LM),PRCCICE(LM,2,LM)
     
c     condensate for all convection
      real*8 mccond(LM,2,LM)
c
      real*8 SCM_LWP_MC,SCM_LWP_SS,SCM_IWP_MC,SCM_IWP_SS,
     &       SCM_WM_MC(LM),SCM_SVWMXL(LM)
c
c     dt and dq over time step from different processes
      real*8 dTtot(LM),dqtot(LM),dTfrc(LM),dqfrc(LM),dTrad(LM)
      real*8 dTHmc(LM),dqmc(LM),dTHbl(LM),dqbl(LM),dTHss(LM),dqss(LM)
      real*8 dTradsw(LM),dTradlw(LM)

      END MODULE SCMDIAG  
