c
c     save diagnostics for run of MODELE SCM

      SUBROUTINE  SCM_DIAG  


      USE RESOLUTION, only: LM
      USE ATM_COM,    only: p,u,v,t,q,qcl,qci,gz,pk
      USE MODEL_COM , only: dtsrc
      USE DYNAMICS,   only: sige,sig
      USE CLOUDS_COM, only: SVLHX,SVLAT,RHSAV,CLDSAV,tauss,taumc,
     &                cldss,cldmc,csizmc,csizss,ncol
      use DIAG_COM, only : npres,ntau,isccp_press,isccp_tau
      USE SCMCOM
      USE SCMDIAG
      USE RAD_COM, only : srhr,trhr
      USE FLUXES, only : atmlnd, atmsrf
      USE CONSTANT, only : SHA, GRAV, kapa, tf 
      USE GEOM, only : axyp 
      USE FILEMANAGER, only : openunit,closeunit
      

      IMPLICIT NONE


C     ALERT!!! LX=57 is for a 53-layer model (i.e., LX=LM+4)
C              LX=40          23-layer model 
C     If LX is changed, need to change data statements in 
C     BLOCK DATA RADPAR that initialize PLB and HLB
CCC   PARAMETER (LX=57,LMOX=12*(1998-1881),MO3X=12*(2050-1950) )
c     PARAMETER (LX=40,LMOX=12*(1998-1881),MO3X=12*(2050-1850+1) )
c     REAL*4 TROAER,VDBCSU,TDUST,EPLMHC,UVLEAN   ! ,FVEG11,FOLGIZ     
c     REAL*4 ATAU98,SIZE98,HTF498,O3CLIM         !  only for offline use
c     CHARACTER*80 ATITLE,VTITLE,DTITLE,TITLE
C
c     COMMON/RADCM1/V6ALB(11,4,7),TAUWC0,TAUIC0,EPSCON,RO3COL,FULGAS(13)
c    A             ,FRAYLE,FALGAE,FCLDTR,FCLDSR,PTLISO,TLGRAD,FGOLDH(13)
c    B             ,WETTRA,WETSRA,FSXAER(5),FTXAER(5),FZSRA(6),FEMTRA(6)
c    C             ,KWVCON,KEEPAL,KEEPRH,KEEP10,KCNORM,KCLDEP,ICE012,NV
c    D             ,MADVES,MRELAY,MOZONE,KO3LON,NO3COL,NORMS0,KSOLAR
c    E             ,KTREND,NTRACE,ITR(8),NL,NLP,MLAT46,MLON72,LASTVC
c    X             ,HLB(LX),RHL(LX),TRACER(LX,8),AERTAU(LX),ZOICE
c    X             ,S00WM2,RATLS0,S0,XXJDAY,JYEARR,JDAYR,ISPARE(98)
c    X             ,DMOICE,DMLICE,TRSLCR,TRDFSL,TRUFSL
C     INPUT DATA
c     COMMON/RADCM2/
c    F              PLB(LX),        TLB(LX),TLT(LX),TLM(LX),U0GAS(LX,12)
c    G             ,ULGAS(LX,12),SHL(LX)
c    H             ,TAUWC(LX),TAUIC(LX),SIZEWC(LX),SIZEIC(LX),CLDEPS(LX)
c    I          ,POCEAN,PEARTH,POICE,PLICE,AGESN(3),SNOWE,SNOWOI,SNOWLI
c    J             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH,      FSPARE(998)
c    K                              ,COSZ,PVT(11),BXA(153),SRBXAL(15,2)
c    L                              ,JLAT,ILON
C     OUTPUT DATA
c    M             ,TRDFLB(LX),TRUFLB(LX),TRNFLB(LX),TRFCRL(LX)
c    N             ,SRDFLB(LX),SRUFLB(LX),SRNFLB(LX),SRFHRL(LX),SRSLHR
c    O             ,SRIVIS,SROVIS,PLAVIS,SRINIR,SRONIR,PLANIR,SRXATM(4)
c    P             ,SRDVIS,SRUVIS,ALBVIS,SRDNIR,SRUNIR,ALBNIR,FSRNFG(4)
c    Q             ,SRTVIS,SRRVIS,SRAVIS,SRTNIR,SRRNIR,SRANIR,FTRUFG(4)
c    R             ,TRDFGW,TRUFGW,TRUFTW,BTEMPW,              DTRUFG(4)
c    S             ,TRSLTS,TRSLTG,TRSLWV,TRSLBS,TTRUFG,LBOTCL,LTOPCL
c    X             ,Z0(12)
c
c-------------------------------------------------------------------------------

      real*8 ARMT(LM),ARMQ(LM)
C
C 
C             Record Layout  
C 
C             NSTEPSCM            Time Stamp - note --- add date/time stamp also    
C             P                   Column Pressure (mb)
C             T          (LM)     Temperature (K) 
C             Q          (LM)     Specific Humidity (Kg/Kg)
C             TSAVG               Ts  -  Surface Air T (K) 
C             GTEMP(1,4,itarg,jtarg) Tskin  -  Sking temperature (C)
C             CLCVSS     (LM)     Cloud Cover SS (by area) 
C             CLCVMC     (LM)     Cloud Cover MC 
C             CLTHCK     (LM)     Cloud Thickness 
C             SCMQCL     (LM)     Cloud Liquid Water Content (Kg/Kg) 
C             SCMQCI     (LM)     Cloud Ice Water Content (Kg/Kg)
C             SVLHXCOL   (LM)     Liquid/Ice Flag (SS) save Latent Heats (j/Kg) 
C             SVLATCOL   (LM)     Liquid/Ice Flag (MC) save Latent Heats (j/Kg)
C             CSIZE      (LM,2)   Particle Size (10**-06m)     1=mc,2=ss 
C             EFFRAD     (LM)     Effective Radius (10**-06m)
C             CUMFLXCOL  (LM)     Cumulus Mass Flux (kg/m**2 /s) 
C             DWNFLXCOL  (LM)     Downdraft Mass Flux (kg/m**2 /s)
C             CUMHET     (LM)     Cumulus Heating  (10**14 W)
C             CUMOST     (LM)     Cumulus Moistening (10**14 W)
C             SRDFLBTOP           INC SW on Top of Atmos (W/m**2) 
C             SRNFLBTOP           SW ABSORBTION BELOW P0 (W/m**2)
C             TRUFLBTOP           UPward LW at P0  (W/m**2)
C             SRDFLBBOT           SW INC on Z0  (W/m**2)
C             SRNFLBBOT           SW ABS on Z0  (W/m**2)
C             TRUFLBBOT           Upward LW at Z0  (W/m**2)
C             TRDFLBBOT           LW INC on Z0  (W/m**2)
C             PRCSS               Precip - Large Scale SS (mm/hour) 
C             PRCMC               Precip - Convective (mm/hour)
C             SCM_PBL_HGT         height of the top of the boundary layer (m)
C             EVPFLX              Evaporation Flux 
C             SHFLX               Sensible Heat Flux 
C             SOILMS              Soil Moisture 
C             SRFHRLCOL(1-LM)*COSZ1  SW Heating 
C             TRFCRLCOL (1-LM)    LW Heating 
C             TAUSSC     (LM)     Cloud Optical Thickness - SS 
C             TAUMCC     (LM)     Cloud Optical Thickness - MC 
C             SG_P       (LM)     Pressure at sigma layers (mb)
C             SG_T       (LM)     Arm T interpolated to Sig levels (K)
C             SG_Q       (LM)     ARM Q interpolated to Sig levels (g/KG)
C             SG_OMEGA   (LM)     ARM Omega (mb/hour)
C             APREC               ARM precipitation (mm/hour)
C             ALH                 ARM Latent heat upwards (W/m2)
C             ASH                 ARM Sensible heat upwards (W/m2)
C             AMEANPS             ARM Area mean Suface Pressure (mb)
C             ATSAIR              ARM Surface Temperature Air (C)
C             ATSKIN              ARM Surface Skin Temperature (C)
C             ARHSAIR             ARM Surface Relative Humidity (%)
C             SG_U       (LM)     ARM u-wind (m/s)
C             SG_V       (LM)     ARM v-wind (m/s)
C             SG_HOR_TMP_ADV(LM)  ARM Horizontal Temperature Advection(K/s) 
C             SG_VER_S_ADV(LM)    ARM Vertical S Advection (K/s)
C             SG_HOR_Q_ADV(LM)    ARM Horizontal Q Advection (kg/kg/s) 
C             SG_VER_Q_ADV(LM)    ARM Vertiacl Q Advection (kg/Kg/s)
C             CLSAV(LM)           SCM cloud fraction SS (by volume)
C             CLDFLG(LM)          SCM Cloud flag from Radia-outcome of Rand(0,1)
C             DWNFLX(LM)          SCM downdraft cloud mass flux (kg/m**2 /s)
C             RHC(LM)             SCM Relative Humidity saved after Cloud 
c                                     routines (Qs(over water))
c             ALWP                ARM MWR cloud liquid water path (cm)
c             ADWDT               ARM d(Column_H2O)/dt (mm/hr)
c             ADWADV              ARM Column_H2O_Advection)_(mm/hr) 
c             ATLWUP              ARM TOA LW Up (W/m**2)
c             ATSWDN              ARM TOA SW Dn (W/m**2)
c             ATSWIN              ARM TOA SW IN (W/m**2)
c             SG_ARSCL(LM)        ARM  ARSCL CLD AMOUNT (%)
c             PRESAV(LM)        SCM Large Scale Precip by layer (kg/kg)
c             PREMC(LM)         SCM MSTCNV Precip by layer (kg/kg)
c             LHPSAV(LM)        SCM Phase of Precip for SS
c             LHPMC(LM)         SCM Phase of Precip for MC
C--- Added by J.W. starting ---C
c             ENTSCM(LM,2)        SCM Entrainment rate for 2 two plume types
c             ENTDEEP(LM,2)       SCM Entrainment rate for deep convection - 2 plumes
c             MPLUMESCM(LM,2)        SCM Mass flux for 2 two plume types
c             MPLUMEDEEP(LM,2)       SCM Mass flux for deep convection - 2 plumes
c             DETRAINDEEP(LM,2,LM)  SCM Detrained convective condensate for Deep conv
C--- Added by J.W. ending ---C
c             WCUSCM(LM,2)        SCM Cumulus updraft speed for 2 two plume types
c             WCUDEEP(LM,2)       SCM  Cumulus updraft speed for deep convection - 2 plumes
c             PRCCDEEP(LM,2,LM)   SCM Precipitating convective condensate for Deep conv
c             NPRCCDEEP(LM,2,LM)  SCM Non-Precipitating conv condensate for Deep conv
c             TPALL(LM,2,LM)      SCM Plume Temperature for Deep Conv
c             MCCOND(LM,2,LM)     SCM convective condensate for deep and shallow  
c             PRCCGRP(LM,2,LM)    SCM deep convective condensate graupel
c             PRCCICE(LM,2,LM)    SCM deep convective condensate ICE
c             SCM_LWP_MC          SCM MC liquid water path (kg/m2)
c             SCM_IWP_MC          SCM MC ice water path (kg/m2)
c             SCM_LWP_SS          SCM SS liquid water path (kg/m2)
c             SCM_IWP_SS          SCM SS ice water path (kg/m2)
c             SCM_WM_MC(LM)       SCM Cloud water for moist convective clouds  kg/kg
c             SRUFLBBOT           Short Wave radiation up at z0 (W/m2)
c             SRUFLBTOP           Short Wave radiation up at p0 (W/m2)
c             TRDFLBTOP           Long Wave radiation down at p0 (W/m2)
c             dTtot(LM)           dT(modelPT)/dt over time step (K/day)
c             dqtot(LM)           dq/dt over time step (kg/kg /day)
c             dTfrc(LM)           dT(modelPT)/dt over time step from FORCN (K/day)
c             dqfrc(LM)           dq/dt over time step from FORCN (kg/kg /day)
c             dTrad(LM)           dT/dti(modelPT) over time step from radiation (K/day)
c             dTHmc(LM)           dTH/dt over time step from mstcnv (K/day)
c             dqmc(LM)            dq/dt over time step from mstcnv (kg/kg/day)
c             dTHbl(LM)           dTH/dt over time step from boundary layer (srf flxs + aturb) (K/day)
c             dqbl(LM)            dq/dt over time step from boundary layer (srf flxs + aturb) (kg/kg/day)
c             dTHss(LM)           dTH/dt over time step from large scale clouds (K/day) 
c             dqss(LM)            dq/dt over time step from large scale clouds (kg/kg/day) 
c

c             isccp record layout 
c
c             isccp_sunlit        ISCCP   1-day 0-night
c             isccp_ctp           ISCCP   cloud top pressure
c             isccp_tauopt        ISCCP   optical thickness
c             isccp_lowcld        ISCCP   low cloud fraction
c             isccp_midcld        ISCCP   mid cloud fraction
c             isccp_highcld       ISCCP   high cloud fraction
c             isccp_fq(ntau,npres)ISCCP  fraction of model grid box
c                                     covered by each of the 49 ISCCP D level cloud 
c                                     types
c             isccp_totcldarea    ISCCP  fraction of model grid box
c                                     columns with cloud somewhere in them
c                                     (sum over all fqI)
c             isccp_boxtau(ncol)  ISCCP optical thickness in each column
c             isccp_boxptop(ncol) ISCCP cloud top pressure (mb) in each column
c            
C              
C--- Added by J.W. starting ---C
      real*8 GZPRT(LM)
C--- Added by J.W. ending ---C
      real*8 TPRT(LM),QPRT(LM),TSURF,TSKIN,SCMQCL(LM),SCMQCI(LM)    
      real*8 TDIFF,QDIFF
      real*8 PCOL, SVLHXCOL(LM),SVLATCOL(LM)    
      real*8 CUMFLXCOL(LM),DWNFLXCOL(LM)
      real*8 daysec,pk1000
      real*8 tt,tr,tmc,tss,tbl
      INTEGER L,LMIN,IC,IU    
      INTEGER IPLUM,IPL,IPLUMSV      
      INTEGER IDEBUG
 
      DATA  daysec/86400./

      if (NSTEPSCM.eq.0) then
          call openunit('scm.save.sige',iu,.true.,.false.)
          WRITE(iu) SIGE
          call closeunit(iu)
          call openunit('scm.output',iu_scm_diag,.true.,.false.)
          if (SCM_ATURB_FLAG.eq.0) then
              write(iu_scm_prt,*) 'RUN with DRYCNV routine '
          elseif (SCM_ATURB_FLAG.eq.1) then
              write(iu_scm_prt,*) 'RUN with ATURB routine '
          endif
      endif

      pk1000 = 1000.**kapa
      PCOL = P(1,1)
      TSURF = atmsrf%TSAVG(1,1)
      TSKIN = atmlnd%GTEMP(1,1)
ccc   SCM_PBL_HGT = atmsrf%dblavg(1,1)   ??????
      SCM_PBL_HGT = 0.0
      do L = 1,LM
C--- Added by J.W. starting ---C
         GZPRT(L) = GZ(1,1,L)
C--- Added by J.W. ending ---C
         TPRT(L) = T(1,1,L)*PK(L,1,1) 
         QPRT(L) = Q(1,1,L)
         SCMQCL(L) = QCL(1,1,L)
         SCMQCI(L) = QCI(1,1,L)
         SVLHXCOL(L) = SVLHX(L,1,1)
         SVLATCOL(L) = SVLAT(L,1,1)
         CLCVSS(L) = CLDSS(L,1,1)
         CLCVMC(L) = CLDMC(L,1,1)
         CLSAV(L) = CLDSAV(L,1,1)   
         TAUSSC(L) = TAUSS(L,1,1)
         TAUMCC(L) = TAUMC(L,1,1)
ccc Now use potential temp  in K/day and q still in kg/kg 
         dTtot(L) = PK1000*dTtot(L)*daysec/dtsrc
         dqtot(L) = dqtot(L)*daysec/dtsrc
         dTfrc(L) = PK1000*dTfrc(L)*daysec/dtsrc
         dqfrc(L) = dqfrc(L)*daysec/dtsrc
         dTrad(L) = PK1000*dTrad(L)*daysec/dtsrc
ccc here change to potential temperature (factor of 1000**kapa)
         dTHmc(L) = PK1000*dTHmc(L)*daysec/dtsrc
         dqmc(L) = dqmc(L)*daysec/dtsrc
         dTHbl(L) = PK1000*dTHbl(L)*daysec/dtsrc
         dqbl(L) = dqbl(L)*daysec/dtsrc
         dTHss(L) = PK1000*dTHss(L)*daysec/dtsrc
         dqss(L) = dqss(L)*daysec/dtsrc
c        write(iu_scm_prt,18) NSTEPSCM,L,
c    *      dTtot(L),dTfrc(L),dTrad(L),
c    *      dTHmc(L),dTHbl(L),dTHss(L) 
c  18    format(1x,'N L dT tot frc rad mc bl ss ',i5,i5,6(f12.3))
      enddo      

      CUMFLXCOL = 0.d0
      DWNFLXCOL = 0.d0
      do LMIN = 1,LM
         do IC = 1,2
            do L=1,LM
               CUMFLXCOL(L) = CUMFLXCOL(L) + CUMFLX(L,IC,LMIN)
               DWNFLXCOL(L) = DWNFLXCOL(L) + DWNFLX(L,IC,LMIN)
            enddo
         enddo
      enddo
c     do L = 1,LM
c        if (CUMFLXCOL(L).gt.0.d0.or.DWNFLXCOL(L).gt.0.d0)
c    &       write(iu_scm_prt,82) L,CUMFLXCOL(L),DWNFLXCOL(L)
c82      format(1x,'DIAG  L CUMFLX DWNFLX ',i5,f12.6,f12.6)
c     enddo


      do L=1,LM
C--- Added by J.W. starting ---C
         ENTSCM(L,1) = 0.0
         ENTSCM(L,2) = 0.0
         MPLUMESCM(L,1) = 0.0
         MPLUMESCM(L,2) = 0.0
C--- Added by J.W. ending ---C
         WCUSCM(L,1) = 0.0
         WCUSCM(L,2) = 0.0
      enddo

      do ic=1,2
         IPLUM = 0
         do LMIN = 1,LM
            do L=1,LM
               if (WCUALL(L,ic,LMIN).ne.0.0d0) then
                   IPLUM = LMIN 
c                  write(iu_scm_prt,24) ic,lmin,iplum
c24                format(1x,
c    *               'WCUPLUM found   ic lmin iplum ',i5,i5,i5)
                   go to 25
               endif
            enddo
         enddo
25       continue
         if (IPLUM.gt.0) then
             do L=1,LM
                WCUSCM(L,ic) = WCUALL(L,ic,IPLUM) 
C--- Added by J.W. starting ---C
                MPLUMESCM(L,ic) = MPLUMEALL(L,ic,IPLUM)
                ENTSCM(L,ic) = ENTALL(L,ic,IPLUM)
C--- Added by J.W. ending ---C
             enddo
         endif
      enddo
c     do ic=1,2
c        do l=1,LM 
c           write(iu_scm_prt,26) ic,l,wcudeep(l,ic)
c26         format(1x,'ic l wcudeep ',i5,i5,f10.4)
c        enddo
c     enddo

C     before writing out diagnostics convert cumulus diagnostics
c     do L = 1,LM
c        write(iu_scm_prt,80) L,CUMHET(L),CUMOST(L)
c 80     format(1x,'L  het mst ',i5,2(f12.3))    
c        CUMHET(L) = CUMHET(L)*10.E-13*SHA*AXYP(1,1)/(GRAV*DTSRC)
c        CUMOST(L) = CUMOST(L)*10.E-13*SHA*AXYP(1,1)/(GRAV*DTSRC)
c     enddo

c     Use the hourly version of the ARM data to save as a diagnostic
      do L = 1,LM
         ARMT(L) = THR(L,NSTEPSCM+IKT)
         ARMQ(L) = QHR(L,NSTEPSCM+IKT)
      enddo

      write(iu_scm_prt,111) isccp_sunlit,isccp_ctp,isccp_tauopt,
     &       isccp_lowcld,isccp_midcld,isccp_highcld,isccp_totcldarea
111   format(1x,'ISCCP diags   sunlit cldptop tau  ',f4.0,f8.2,f8.2,
     &           '    low mid high tot ',4(f8.3))

c     do ic = 1,ncol
c        write(iu_scm_prt,112) ic,isccp_boxptop(ic),isccp_boxtau(ic)
c112     format(1x,' ic   ptop   tau  ',i5,f9.2,f10.4)
c     enddo

      write(iu_scm_prt,113) (isccp_tau(ic),ic=1,7)
 113  format(1x,10x,7(f10.4))
      do L = 1,npres
         write(iu_scm_prt,114) isccp_press(L),
     &                          (isccp_fq(ic,L),ic=1,ntau)
 114     format(1x,I10,7(f10.4))
      enddo

      WRITE(iu_scm_diag) NSTEPSCM,PCOL,TPRT,QPRT,TSURF,TSKIN,CLCVSS,
     *           CLCVMC,SCMQCL,SCMQCI,SVLHXCOL,SVLATCOL,CSIZE,EFFRAD,
     *           CUMFLXCOL,DWNFLXCOL,SRDFLBTOP,SRNFLBTOP,TRUFLBTOP,
     *           SRDFLBBOT,SRNFLBBOT,TRUFLBBOT,TRDFLBBOT,PRCSS,PRCMC,
     *           SCM_PBL_HGT,EVPFLX,SHFLX,SRFHRLCOL,
     *           TRFCRLCOL,TAUSSC,TAUMCC,SG_P,ARMT,ARMQ,
     *           SG_OMEGA,APREC,ALH,ASH,AMEANPS,ATSAIR,ATSKIN,
     *           SG_U,SG_V,SG_HOR_TMP_ADV,
     *           SG_VER_S_ADV,SG_HOR_Q_ADV,SG_VER_Q_ADV,CLSAV,
     *           ATSWDN,ATSWIN,SG_ARSCL,
     *           SRUFLBTOP,SRUFLBBOT,TRDFLBTOP, 
     *           dTtot,dqtot,dTfrc,dqfrc,dTrad,dTHmc,dqmc,
     *           dTHbl,dqbl,dTHss,dqss
     *            ,isccp_sunlit,isccp_ctp,
     *           isccp_tauopt,isccp_lowcld,isccp_midcld,isccp_highcld,
     *           isccp_fq,isccp_totcldarea,isccp_boxtau,isccp_boxptop


C 
      write(iu_scm_prt,99) NSTEPSCM
 99   format(//1x,'END OF TIME STEP      NSTEPSCM  ',i6)
      write(iu_scm_prt,100) NSTEPSCM,PCOL,TSURF,TSKIN+tf,ATSAIR+tf,
     *                      ATSKIN+tf
100   format(1x,'NSTEP P tsurfair tskin  ',i5,5(f10.4))
      write(iu_scm_prt,110) PRCSS,PRCMC
110   format(1x,'PRCSS MC ',2(f10.4))
      write(iu_scm_prt,120) SRNFLBBOT,SRNFLBTOP,SRDFLBBOT,SRDFLBTOP, 
     *              TRUFLBTOP,TRUFLBBOT,TRDFLBBOT
120   format(1x,'rad  sw nbot top dbot top ',4(f10.2),
     *       '    lw utop bot dbot ',3(f10.2)) 
      write(iu_scm_prt,125) EVPFLX,SHFLX
125   format(1x,'EVPFLX SHFLX ',2(f12.6))

c
      IDEBUG = 0
      do L=1,LM
         write(iu_scm_prt,140) L,SG_P(L),TPRT(L),
     +         QPRT(L)*1000.0,SCMQCL(L)*1000.0,SCMQCI(L)*1000.,
     +         SCM_SVWMXL(L)*1000.0,TAUSSC(L),TAUMCC(L),
     +                 CLCVSS(L)*100.,CLCVMC(L)*100.,SG_ARSCL(L)
 140     format(1x,i2,f8.2,' T ',f7.2,' Q',f7.3,' qcl qci det',
     +          3(f7.3),' tauss mc',2(f7.2),' cfss mc',
     +          2(f7.2),' cld',f5.1)
      enddo 

c     do L=1,LM
c        write(iu_scm_prt,140) L,SG_P(L),TPRT(L),ARMT(L),
c    +         QPRT(L)*1000.0,ARMQ(L)*1000.0,SCMQCL(L)*1000.0,
c    +         SCMQCI(L)*1000.0,TAUSSC(L),TAUMCC(L),
c    +         CLCVSS(L)*100.,CLCVMC(L)*100.,SG_ARSCL(L)
c140     format(1x,i3,f8.2,' T ',2(f7.2),'  Q ',2(f7.3),'  qcl qci ',
c    +          2(f7.4),'  tauss mc ',2(f8.3),'  camtss mc ',
c    +          2(f7.3),' arscl ',f6.2)
c     enddo 

      RETURN 

      END SUBROUTINE SCM_DIAG 
