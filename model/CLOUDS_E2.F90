#include "rundeck_opts.h"

#define ALT_CDNC_INPUTS

module CLOUDS

!@sum  CLOUDS column physics of moist conv. and large-scale condensation
!@auth M.S.Yao/A. Del Genio (modifications by Gavin Schmidt)
!@cont MSTCNV,LSCOND,ANVIL_OPTICAL_THICKNESS,MC_CLOUD_FRACTION,
!@+    CONVECTIVE_MICROPHYSICS,MC_PRECIP_PHASE,MASS_FLUX,PRECIP_MP
  use CONSTANT, only : rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6 &
       ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi &
       ,mb2kg
  use RESOLUTION, only : lm
  use MODEL_COM, only : dtsrc,itime
  use TimeConstants_mod, only: SECONDS_PER_HOUR
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
  use CONSTANT, only : kapa,mair,gasc
  use RESOLUTION, only : ptop,psf,ls1
  use DYNAMICS, only : sig,sige
#endif
#ifdef SCM
  use SCMCOM, only: SCM_SAVE_T,SCM_SAVE_Q,SCM_DEL_T, &
       SCM_DEL_Q,SCM_ATURB_FLAG,iu_scm_prt,NRINIT,I_TARG,J_TARG
  use SCMDIAG, only : WCUSCM,WCUALL,WCUDEEP,PRCCDEEP,NPRCCDEEP, &
       MPLUMESCM,MPLUMEALL,MPLUMEDEEP, &
       ENTSCM,ENTALL,ENTDEEP,DETRAINDEEP, &
       TPALL,PRCCGRP,PRCCICE,MCCOND, &
       PRESAV,LHPSAV,PREMC,LHPMC, &
       CUMFLX,DWNFLX,SCM_LWP_MC,SCM_LWP_SS, &
       SCM_IWP_MC,SCM_IWP_SS,SCM_WM_MC
#endif

  use CLOUDS_COM, only : ncol

#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
  use mo_bulk2m_driver_gcm, only: execute_bulk2m_driver
#endif
  use QUSDEF, only : nmom,xymoms,zmoms,zdir
#ifdef TRACERS_ON
  use TRACER_COM, only: NTM, trname,t_qlimit,ntm_soa,ntm_ococean
#ifdef TRACERS_AEROSOLS_OCEAN
  use TRACER_COM, only: n_ococean,n_seasalt1,trm,trpdens
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_WATER
  use TRACER_COM, only: nWATER, tr_wd_TYPE, &
       tr_evap_fact, gases_list,gases_count
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  use TRACER_COM, only: aqchem_list,aqchem_count
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
  use TRACER_COM, only: Ntm_dust
#endif
#endif
#endif

#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
#ifdef TRACERS_AMP
  use CLOUDS_COM, only: NACTC,NAERC
  use AERO_CONFIG, only: NMODES
#endif
#endif
  implicit none
  save
  !**** parameters and constants
  real*8, parameter :: TI=233.16d0   !@param TI pure ice limit
  real*8, parameter :: CLDMIN=.10d0 !@param CLDMIN min MC/LSC region
!@param WMU critical cloud water content for rapid conversion (g m**-3)
  real*8, parameter :: WMU=.25
  real*8, parameter :: WMUL=.5       !@param WMUL WMU over land
  !     REAL*8, PARAMETER :: WMUI=.1d0     !@param WMUI WMU for ice clouds
  real*8 WMUI                          !@param WMUI WMU for ice clouds
  real*8 WMUSI        !@param WMUSI WMU for liquid clouds over sea-ice
  real*8, parameter :: BRCLD=.2d0    !@param BRCLD for cal. BYBR
  real*8, parameter :: FDDET=.25d0 !@param FDDET remainder of downdraft
  real*8, parameter :: DTMIN1=1.d0 !@param DTMIN1 min DT to stop downdraft drop
  real*8, parameter :: SLHE=LHE*BYSHA
  real*8, parameter :: SLHS=LHS*BYSHA
!@param CCMUL multiplier for convective cloud cover
!@param CCMUL1 multiplier for deep anvil cloud cover
!@param CCMUL2 multiplier for shallow anvil cloud cover
!@param COETAU multiplier for convective cloud optical thickness
  real*8, parameter :: CCMUL=2.,CCMUL1=5.,CCMUL2=3.,COETAU=.08d0

  real*8 :: RTEMP,CMX,RCLDX,WMUIX,CONTCE1,CONTCE2,TNX,QNX
  real*8 :: BYBR,BYDTsrc,XMASS,PLAND
!@var BYBR factor for converting cloud particle radius to effect. radius
!@var XMASS dummy variable
!@var PLAND land fraction

  !**** Set-able variables
!@dbparam LMCM max level for originating MC plumes
  integer :: LMCM = -1 ! defaults to LS1-1 if not set in rundeck
!@dbparam ISC integer to turn on computation of stratocumulus clouds
  integer :: ISC = 0  ! set ISC=1 to compute stratocumulus clouds
  !     REAL*8 :: U00MAX = .99d0      ! maximum U00 for water clouds
  !****
  !**** WARNING: U00wtrX AND U00ice ARE NO LONGER USED BY THE GCM. USE U00a and U00b INSTEAD
  !****
!@dbparam U00wtrX multiplies U00ice for critical humidity for water clds
  real*8 :: U00wtrX = 1.0d0     ! default, needed for AR4 runs
!@dbparam U00ice critical humidity for ice cloud condensation
  real*8 :: U00ice = .7d0       ! default, needed for AR4 runs
!@dbparam U00a tuning knob for U00 above 850 mb without moist convection
!@dbparam U00b tuning knob for U00 below 850 mb and in convective regions
  real*8 :: U00a = 0.55d0       ! default
  real*8 :: U00b = 1.00d0       ! default
!@dbparam funio_denominator funio denominator
  real*8 :: funio_denominator=22.d0  ! default
!@dbparam autoconv_multiplier autoconversion rate multiplier
  real*8 :: autoconv_multiplier=1.d0 ! default
!@dbparam radius_multiplier cloud particle radius multiplier
  real*8 :: radius_multiplier=1.d0   ! default
!@dbparam wmui_multiplier critical ice cloud water multiplier
  real*8 :: wmui_multiplier=1.d0     ! default
!@dbparam entrainment_cont1 constant for entrainment rate, plume 1
  real*8 :: entrainment_cont1=.3d0   ! default
!@dbparam entrainment_cont2 constant for entrainment rate, plume 2
  real*8 :: entrainment_cont2=.6d0   ! default
!@dbparam HRMAX maximum distance an air parcel rises from surface
  real*8 :: HRMAX = 1000.d0     ! default (m)
!@dbparam RIMAX maximum ice cloud size
!@dbparam RWMAX maximum water cloud size
  real*8 :: RIMAX = 100.d0, RWMAX = 20.d0      ! microns
!@dbparam RWCldOX multiplies part.size of water clouds over ocean
  real*8 :: RWCldOX=1.d0
!@dbparam RICldX multiplies part.size of ice clouds at 1000mb
!@+       RICldX changes linearly to 1 as p->0mb
  real*8 :: RICldX=1.d0 , xRICld
!@dbparam do_blU00 =1 if boundary layer U00 is treated differently
  integer :: do_blU00=0     ! default is to disable this

#ifdef TRACERS_ON
!@var ntx,NTIX: Number and Indices of active tracers used in convection
  integer, allocatable, dimension(:) :: ntix
  integer ntx
#endif
  !**** ISCCP diag related variables
!!! @parameter ncol used to be set here  20 for gcm runs and 100 for scm runs
!!!     now moved to CLOUDS_COM.f  for portability
!@var tautab look-up table to convert count value to optical thickness
!@var invtau look-up table to convert optical thickness to count value
  real*8 :: tautab(0:255)
  integer :: invtau(-20:45000)

  !**** input variables
  logical DEBUG
!@var RA ratio of primary grid box to secondary gridbox
  real*8, dimension(:), allocatable :: RA !(KMAX)
!@var UM,VM,UM1,VM1,U_0,V_0 velocity related variables(UM,VM)=(U,V)*AIRM
  real*8, dimension(:,:), allocatable :: UM,VM,UM1,VM1 !(KMAX,LM)
  real*8, dimension(:,:), allocatable :: U_0,V_0       !(KMAX,LM)

!@var Miscellaneous vertical arrays set in driver
!@var PLE pressure at layer edge
!@var LHP array of precip phase ! may differ from LHX
  real*8, dimension(LM+1) :: PLE,LHP
  real*8, dimension(LM) :: PL,PLK,AIRM,BYAM,ETAL,TL,QL,TH,RH,WMX &
       ,VSUBL,MCFLX,DGDSM,DPHASE,DTOTW,DQCOND,DGDQM,AQ,DPDT,RH1 &
       ,FSSL,VLAT,DDMFLX,WTURB,TVL,W2L,GZL &
       ,SAVWL,SAVWL1,SAVE1L,SAVE2L,DPHASHLW,DPHADEEP,DGSHLW,DGDEEP &
       ,QDNL,TDNL,U00L
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
  real*8, dimension(LM) :: WMXICE
#endif
  real*8, dimension(LM) :: DQMTOTAL,DQMSHLW,DQMDEEP &
       ,DQCTOTAL,DQCSHLW,DQCDEEP,DQLSC
!@var PL layer pressure (mb)
!@var PLK PL**KAPA
!@var AIRM the layer's pressure depth (mb)
!@var BYAM 1./AIRM
!@var ETAL fractional entrainment rate
!@var TL, QL temperature, specific humidity of the layer
!@var TH potential temperature (K)
!@var RH relative humidity
!@var RH1 relative humidity to compare with the threshold humidity
!@var WMX cloud water mixing ratio (kg/kg)
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
!@var WMXICE ice water mixing ratio (kg/kg)
#endif
!@var VSUBL downward vertical velocity due to cumulus subsidence (cm/s)
!@var MCFLX, DGDSM, DPHASE, DQCOND, DGDQM dummy variables
!@var DDMFLX accumulated downdraft mass flux (mb)
!@var AQ time change rate of specific humidity (s**-1)
!@var DPDT time change rate of pressure (mb/s)
!@var FSSL grid fraction for large-scale clouds
!@var VLAT dummy variable
  real*8, dimension(LM+1) :: PRECNVL
!@var WTURB turbulent vertical velocity (m)
!@var PRECNVL convective precip entering the layer top
  !**** new arrays must be set to model arrays in driver (before MSTCNV)
  real*8, dimension(LM) :: SDL,WML
!@var SDL vertical velocity in sigma coordinate
!@var WML cloud water mixing ratio (kg/kg)
  !**** new arrays must be set to model arrays in driver (after MSTCNV)
  real*8, dimension(LM) :: TAUMCL,SVLATL,CLDMCL,SVLHXL,SVWMXL,SVLAT1
!@var TAUMCL convective cloud optical thickness
!@var SVLATL saved LHX for convective cloud
!@var CLDMCL convective cloud cover
!@var SVLHXL saved LHX for large-scale cloud
!@var SVWMXL saved detrained convective cloud water
  real*8, dimension(LM) :: CSIZEL
!@var CSIZEL cloud particle radius (micron)
#ifdef CLD_AER_CDNC
  real*8, dimension(LM) :: ACDNWM,ACDNIM
!@var ACDNWM,ACDNIM -CDNC - warm and cold moist cnv clouds (cm^-3)
  real*8, dimension(LM) :: ACDNWS,ACDNIS
!@var ACDNWS,ACDNIS -CDNC - warm and cold large scale clouds (cm^-3)
  real*8, dimension(LM) :: CDNC_NENES,CDNC_TOMAS
!@var CDNC_TOMAS, CDNC_NENS -CDNC from Nenes and Seinfel parameterization- warm large scale clouds (cm^-3)
  real*8, dimension(LM) :: AREWS,AREIS,AREWM,AREIM  ! for diag
!@var AREWS and AREWM are moist cnv, and large scale Reff arrays (um)
  real*8, dimension(LM) :: ALWWS,ALWIS,ALWWM,ALWIM  ! for diag
  real*8, dimension(LM) :: CDN3DL,CRE3DL
!@var ALWWM and ALWIM  etc are liquid water contents
!@var SMLWP is LWP
  real*8 SMLWP
!@var SME is the TKE in 1 D from e(l) = egcm(l,i,j)  (m^2/s^2)
  real*8, dimension(LM)::SME
  integer NLSW,NLSI,NMCW,NMCI
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
!@var CTTEM,CD3DL,CL3DL,CI3DL are cld temp, cld thickness,cld water
  real*8, dimension(LM) ::CTEML,CD3DL,CL3DL,CI3DL
#endif
  !**** new arrays must be set to model arrays in driver (before LSCOND)
  real*8, dimension(LM) :: TTOLDL,CLDSAVL,CLDSV1
!@var TTOLDL previous potential temperature
!@var CLDSAVL saved large-scale cloud cover
#ifdef CLD_AER_CDNC
  real*8, dimension(LM)::OLDCDL,OLDCDI
!@var OLDCDL is saved CDNC
!@var OLDCDI is saved ice crystal numbe
#endif
  !**** new arrays must be set to model arrays in driver (after LSCOND)
  real*8, dimension(LM) :: SSHR,DCTEI,TAUSSL,CLDSSL
!@var SSHR,DCTEI height diagnostics of dry and latent heating by MC
!@var TAUSSL large-scale cloud optical thickness
!@var CLDSSL large-scale cloud cover

!@var SM,QM Vertical profiles of (T/p**kappa)*AIRM, q*AIRM
  real*8, dimension(LM) :: SM,QM
  real*8, dimension(NMOM,LM) :: SMOM,QMOM,SMOMMC,QMOMMC, &
       SMOMLS,QMOMLS

#ifdef TRACERS_ON
!@var TM Vertical profiles of tracers
  real*8, allocatable, dimension(:,:) :: TM
  real*8, allocatable, dimension(:,:,:) :: TMOM
!@var TRDNL tracer concentration in lowest downdraft (kg/kg)
  real*8, allocatable, dimension(:,:) :: TRDNL
#ifdef TRACERS_WATER
!@var TRWML Vertical profile of liquid water tracers (kg)
!@var TRSVWML New liquid water tracers from m.c. (kg)
  real*8, allocatable, dimension(:,:) :: TRWML, TRSVWML
!@var TRPRSS super-saturated tracer precip (kg)
!@var TRPRMC moist convective tracer precip (kg)
  real*8, allocatable, dimension(:)    :: TRPRSS,TRPRMC
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
  ! for diagnostics
  real*8, allocatable, dimension(:,:) :: DT_SULF_MC,DT_SULF_SS
#endif
#ifdef TRDIAG_WETDEPO
!@dbparam diag_wetdep switches on/off special diags for wet deposition
  integer :: diag_wetdep=0 ! =off (default) (on: 1)
!@var trcond_mc saves tracer condensation in MC clouds [kg]
!@var trdvap_mc saves tracers evaporated in downdraft of MC clouds [kg]
!@var trflcw_mc saves tracers condensed in cloud water of MC clouds [kg]
!@var trprcp_mc saves tracer precipitated from MC clouds [kg]
!@var trnvap_mc saves reevaporated tracer of MC clouds precip [kg]
!@var trwash_mc saves tracers washed out by collision for MC clouds [kg]
  real*8,allocatable, dimension(:,:) :: trcond_mc,trdvap_mc,trflcw_mc, &
       trprcp_mc,trnvap_mc,trwash_mc
!@var trwash_ls saves tracers washed out by collision for LS clouds [kg]
!@var trprcp_ls saves tracer precipitation from LS clouds [kg]
!@var trclwc_ls saves tracers condensed in cloud water of LS clouds [kg]
!@var trevap_ls saves reevaporated tracers of LS cloud precip [kg]
!@var trclwe_ls saves tracers evaporated from cloud water of LS clouds [kg]
!@var trcond_ls saves tracer condensation in LS clouds [kg]
  real*8,allocatable,dimension(:,:) :: trwash_ls,trevap_ls,trclwc_ls, &
       trprcp_ls,trclwe_ls,trcond_ls
#endif
#else
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
!@var tm_dust vertical profile of dust/mineral tracers [kg]
  real*8,dimension(Lm,Ntm_dust) :: tm_dust
!@var tmom_dust vertical profiles of dust/mineral tracer moments [kg]
  real*8,dimension(nmom,Lm,Ntm_dust) :: tmom_dust
!@var trprc_dust dust/mineral tracer precip [kg]
  real*8,dimension(Lm,Ntm_dust) :: trprc_dust
#endif
#endif
#endif

!@var KMAX index for surrounding velocity
!@var LP50 50mb level
  integer ::  KMAX,LP50
!@var PEARTH fraction of land in grid box
!@var TS average surface temperture (C)
!@var RIS, RI1, RI2 Richardson numbers
  real*8 :: PEARTH,TS,QS,US,VS,RIS,RI1,RI2,DXYPIJ,ROICE
!@var DCL max level of planetary boundary layer
  integer :: DCL

  !**** output variables
  real*8 :: PRCPMC,PRCPSS,HCNDSS,WMSUM
!@var PRCPMC precip due to moist convection
!@var PRCPSS precip due to large-scale condensation
!@var HCNDSS heating due to large-scale condensation
!@var WMSUM cloud liquid water path
#ifdef CLD_AER_CDNC
  real*8 :: WMCLWP,WMCTWP
!@var WMCLWP , WMCTWP moist convective LWP and total water path
#endif
  real*8 :: CLDSLWIJ,CLDDEPIJ
!@var CLDSLWIJ shallow convective cloud cover
!@var CLDDEPIJ deep convective cloud cover
  integer :: LMCMAX,LMCMIN
!@var LMCMAX upper-most convective layer
!@var LMCMIN lowerest convective layer
!@var AIRXL is convective mass flux (mb)
  real*8 AIRXL,PRHEAT
!@var RNDSSL stored random number sequences
  real*8  RNDSSL(3,LM)
!@var prebar1 copy of variable prebar
  real*8 prebar1(Lm+1)

#ifdef TRACERS_ON
  ! The following tracer arrays are workspace for MSTCNV.  They are
  ! declared as permanent arrays here to avoid the expense of initializing
  ! temporary LM,NTM arrays to zero each time MSTCNV is called.  After
  ! completion of MC calculations, MSTCNV resets these arrays to zero in
  ! the layers in which they were used.
!@var DTM,DTMR: Vertical profiles of Tracers changes
  real*8, allocatable, dimension(:,:)      :: DTM, DTMR, TMDNL
  real*8, allocatable, dimension(:,:,:) :: DTMOM, DTMOMR, TMOMDNL
!@var TPOLD saved plume temperature after condensation for tracers
!@+   (this is slightly different from TPSAV)
  real*8, dimension(LM)       :: TPOLD=0
#ifdef TRACERS_WATER
!@var TRCOND tracer mass in condensate
!@var TRCONDV tracer mass in lofted condensate
  real*8, allocatable, dimension(:,:)   :: TRCOND,TRCONDV
#endif
#endif

contains

  subroutine MSTCNV(IERR,LERR,i_debug,j_debug)

!@sum  MSTCNV moist convective processes (precip, convective clouds,...)
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@calls adv1d,QSAT,DQSATDT,THBAR, MASS_FLUX, CONVECTIVE_MICROPHYSICS, PRECIP_MP,
!@+     MC_PRECIP_PHASE, MC_CLOUD_FRACTION, ANVIL_OPTICAL_THICKNESS

    !**** FREE PARAMETERS THAT USERS MIGHT WANT TO VARY INCLUDE
    !****    (1) ADJUSTMENT TIME FOR STABILIZATION OF CLOUD BASE BY CUMULUS MASS
    !****        FLUX (TADJ): DEFAULT VALUE = 1.0 (HOUR)
    !****    (2) SCALING FACTOR FOR ENTRAINMENT STRENGTH (CONTCE): DEFAULT VALUES
    !****        CONTCE1 = 0.3, CONTCE2 = 0.6
    !****    (3) SCALING FACTOR FOR EFFECT OF PRESSURE GRADIENT ON CONVECTIVE
    !****        HORIZONTAL MOMENTUM TRANSPORT (PGRAD): DEFAULT VALUE = 0.7
    !****    (4) FRACTION OF UPDRAFT MASS ASSUMED FOR INITIAL DOWNDRAFT MASS (ETADN):
    !****        DEFAULT VALUE = 1./3.
    !****    (5) FRACTION OF PRECIPITATING CONDENSATE AVAILABLE FOR RE-EVAPORATION
    !****        RATHER THAN INCORPORATION INTO DOWNDRAFT (FDDRT): DEFAULT VALUE = 1.
    !****
    !**** NOTE: THESE PARAMETERS ARE INTENDED FOR RESEARCH PURPOSES ONLY AND
    !**** SHOULD NOT BE USED TO ADJUST THE MODEL TO RADIATION BALANCE
    !****

    implicit none
    !
    !****
    !              *******************************
    !              ***   DECLARATION SECTION   ***
    !              *******************************
    !                        FUNCTIONS
    !                        PARAMETERS
    !                        VARIABLES
    !
    !              *******************************
    !
    !                  ***               ***
    !                  ***   FUNCTIONS   ***
    !                  ***               ***
    !
    real*8 ::   DQSATDT, PRECIP_MP, QSAT, THBAR
    !
!@var DQSATDT     dQSAT/dT
!@var PRECIP_MP   mass density of precipitating condensate (kg/m^3)
!@var QSAT        saturation specific humidity
!@var THBAR       virtual temperature at layer edge
    !
    !              *******************************
    !
    !                  ***                ***
    !                  ***   PARAMETERS   ***
    !                  ***                ***
    !
    real*8,   parameter :: AIRM0=100.d0
    real*8,   parameter :: CK1 = 1.,  CN0=8.d6,  CN0I=8.d6,  CN0G=8.d6
    integer,  parameter :: ITMAX=50
    real*8,   parameter :: FITMAX=1d0/ITMAX
    real*8,   parameter :: PN=1.d0,   RHOG=400., RHOIP=100.
    !
!@param AIRM0 air mass used to compute convective cloud cover
!@param CK1 a tunning const.
!@param CN0, CN0I, CN0G intercepts of Marshall-Palmer particle size dist.
!@param ITMAX max iteration indices
!@param FITMAX set to 1/ITMAX
!@param PN tuning exponential for computing WV
!@param RHOG,RHOIP density of graupel and ice particles
    !
#ifdef CLD_AER_CDNC
    integer, parameter :: SNTM=31  !for tracers for CDNC
#endif
    !
    !              *******************************
    !
    !                  ***               ***
    !                  ***   VARIABLES   ***
    !                  ***               ***
    !
    !              Sections:
    !                        ARRAYS   - DECLARATIONS
    !                                 - DEFINITIONS
    !                        SCALARS  - DECLARATIONS
    !                                 - DEFINITIONS
    !                        VARIABLES FOR CONDITIONAL COMPILIATIONS
    !
    !
    !        *********   ARRAY DECLARATIONS   **********
    !
    real*8, dimension(0:LM) :: CM
    !
    real*8, dimension(LM) :: &
         BUOY, &
         CCM,  CDHEAT, CMNEG, COND, CONDP, CONDP1, CONDV, CONDGP, &
         CONDIP, &
         DQM,  DSM,    DQMR,  DSMR, DET,   DM,     DMR,   DDM,    DDR, &
         !cc  *   ENT,  F,      HEAT1, LHP,  ML,    QM1,    QMT,   QMDNL,  QMOLD,
         ENT,  F,      HEAT1,       ML,    QM1,    QMT,   QMDNL,  QMOLD, &
         SM1,  SMT,    SMDNL, SMOLD,TPSAV, TAUMC1, WCU,   WCU2
    !
    real*8, dimension(KMAX) :: &
         SUMU,SUMV,SUMU1,SUMV1,UMP,VMP,UMDN,VMDN
    !
    real*8, dimension(KMAX,LM) :: DUM,DVM,UMDNL,VMDNL
    !
    real*8, dimension(NMOM) :: &
         SMOMP,QMOMP, SMOMPMAX,QMOMPMAX, SMOMDN,QMOMDN
    !
    real*8, dimension(NMOM,LM) :: &
         DQMOM, DQMOMR, DSMOM, DSMOMR, FMOM, QMOMDNL, QMOMOLD, &
         SMOMDNL, SMOMOLD
    !
    !        *********   ARRAY DEFINITIONS   **********
    !
    !        (0:LM)
!@var CM            air mass of subsidence
    !        (LM)
!@var BUOY
!@var CCM           convective plume mass (mb)
!@var CDHEAT        heating due to condensation
!@var CMNEG
!@var COND,CONDP, CONDP1, CONDV, CONDGP,CONDIP condensate mass density (kg/m^3)
!@var DQM,DSM,DQMR,DSMR Vertical profiles of T/Q and changes
!@var DDM           downdraft mass (mb)
!@var DM, DMR       change in air mass
!@var ENT
!@var F
!@var HEAT1         heating needed for phase change
!@var ML            layer air mass (mb)
!@var SM1, QM1, SMT, QMT dummy variables
!@var SMOLD,QMOLD   profiles prior to any moist convection
!@var SMDNL
!@var TPSAV         array to save plume temperature (set once) (K)
!@var TAUMC1
!@var WCU,  WCU2
    !        (KMAX)
!@var SUMU,SUMV,SUMU1,SUMV1
!@var UMP, VMP      momentum carried by convective plumes
!@var UMDN,VMDN     dummy variables
    !        (KMAX,LM)
!@var DUM, DVM      changes of UM,VM
!@var UMDNL,VMDN
    !        (NMOM)
!@var SMOMP,SMOMPMAX,SMOMDN
!@var QMOMP,QMOMPMAX,QMOMDN
!@var
    !        (NMOM,LM)
!@var DQMOM, DQMOMR, DSMOM, DSMOMR
!@var FMOM
!@var QMOMDNL, QMOMOLD, SMOMDNL, SMOMOLD
    !
    !
    !          *********   SCALAR DECLARATIONS   *********
    !
    logical  BELOW_CLOUD, MC1
    integer, intent(IN)  :: I_DEBUG, J_DEBUG
    integer, intent(OUT) :: IERR,LERR
    integer  IC, ITER, ITYPE, IERRT, K, KSUB, &
         L, LLMIN, LFRZ,  LERRT, LDRAFT, LMIN, LMAX, LDMIN, LM1, &
         MCCONT, MAXLVL, MINLVL, N, NPPL, NSUB

    real*8 &
         ALPHA,ALPHAU, BETA,BETAU,BYKSUB,BYPBLM, &
         CDHM,CDHSUM,CLDM,CLDREF,CDHSUM1,CONTCE,CDHDRT,CONDMU, &
         !
         DCW,DCG,DCI,DWCU,DMSE,DFP,DQSUM,DMMIX,DQ,DMSE1,DQSUM1, &
         DDRAFT,DELTA,DQEVP,DDRUP,DDROLD, &
         EPLUME,ETADN,ETAL1,EVPSUM,EDRAFT, &
         !
         FLAMW,FLAMG,FLAMI,FG,FI,FMC1,FCLW, &
         FMP0,FPLUME,FMP2,FRAT1,FRAT2,FCTYPE,FQCOND1,FPOLD, &
         FDDRT,FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FENTR,FENTRA, &
         FEVAP,FLEFT,FQCOND,FQCONDV,FQEVP,FPRCP,FSEVP,FSSUM, &
         FCONV_tmp,FSUB_tmp,FSSL_tmp, &
         !
         GAMA, HDEP,HPBL, LHX,LHX1, &
         MPLUME,MPLUM1,MCLOUD,MPMAX,MPOLD, &
         MNdO,MNdL,MNdI,MCDNCW,MCDNCI,      & ! Menon
         PBLM,PRCP,PGRAD, &
         !
         QENV,QMO1,QMO2,QDN,QUP,QEDGE,QMN1,QMN2,QMP,QMDN,QMIX, &
         QMPMAX,QMPT,QSATC,QSATMP, RCLD,RCLDE,RHO, &
         !
         SENV,SMO1,SMO2,SDN,SUP,SEDGE,SVDN,SVUP,SVEDG,SMN1,SMN2, &
         SMP,SLH,SMDN,SMIX,SMPMAX,SMPT,SUMAJ,SVMIX,SVM1,SUMDP, &
         !
         TADJ,TEMWM,TEM,TIG,TNX1,TP,TVP,TOLD,TOLD1,TTURB,TRATIO, &
         UMTEMP,VMTEMP,VT, &
         WMDN,WMUP,WMEDG,WMIX,W2TEM,WTEM,WCONST,WCUFRZ,WORK,WMAX,WV
    !
    !
    !            *********   SCALAR DEFINITIONS   *********
    !
    !
!@var     Logical Variables
!@var BELOW_CLOUD      is the current level below cloud?
!@var MC1              true for the first convective event
    !
!@var     Integer Variables
!@var IC               integer for cloud types
!@var ITER             number for iteration
!@var ITYPE            convective cloud types
!@var IERR,LERR        error reports from advection
!@var IERRT,LERRT      error reports from advection
!@var K,L,N            loop variables
!@var KSUB, BYKSUB     number, 1/number of subsidence iterations
!@var LDRAFT           the layer at which the downdraft orginates
!@var LMIN, LMAX       the base, top layers of a convective event
!@var LDMIN            the lowest layer to which the downdraft descends
!@var LFRZ             freezing level
!@var LM1
!@var MCCONT           integer to count convective events
!@var MAXLVL, MINLVL   the lowest, the highest layer of convective events
!@var NPPL             iteration counter for computing MC area partition
!@var NSUB             LMAX - LMIN + 1
    !
!@var     REAL*8  Dummy Variables
!@var       ALPHA,ALPHAU, BETA,BETAU,BYPBLM, CDHDRT,CDHM,CDHSUM,CLDREF,
!@var       DQSUM,DQ,DDRUP,DDROLD, EVPSUM,
!@var       FCDH,FCDH1,FCLD,FCLOUD,FDDL,FDDP,FMP2,FRAT1,FRAT2,  GAMA,
!@var       QMN1,QMN2,QMDN,QMIX,QMPT,QENV,QMO1,QMO2,QDN,QUP,QEDGE,QMPT,QNX,
!@var       SMN1,SMN2,SMDN,SMIX,SMPT,SENV,SMO1,SMO2,SDN,SUP,SEDGE,SMPT,
!@var       SUMAJ,SUMDP,SVDN,SVUP,SVEDG,
!@var       TEMWM,TEM,TIG,TNX,TTURB,TRATIO,  WTEM,WCONST,WMDN,WMUP,WMEDG
    !
!@var     REAL*8 Variables
!@var CONDMU        convective condensate in Kg/m^3
!@var CLDM          subsidence due to convection (mb)
!@var CONTCE        scaling factor for entrainment strength
!@var DCW,DCG,DCI   critical cloud particle sizes for onset of precip
!@var DMSE, DMSE1   difference in moist static energy
!@var DFP           an iterative increment
!@var DMMIX
!@var DDRAFT        downdraft mass (mb)
!@var DELTA         fraction of plume that stays in the layer
!@var DQEVP         amount of condensate that evaporates in downdrafts
!@var EPLUME        mass of entrained air (mb)
!@var ETADN         initial downdraft mass / updraft mass
!@var ETAL1         fractional entrainment rate
!@var EDRAFT        mass of air entrained into downdrafts (mb)
!@var FLAMW,FLAMG,FLAMI Marshall-Palmer lambda for water, graupel and ice
!@var FG, FI        fractions of glaciated mass in graupel and ice
!@var FMC1          fraction of grdibox occupied by moist convection + subsidence
!@var FCLW          fraction of condensate in plume that remains as CLW
!@var FMP0          less entraining convective mass (mb)
!@var FPLUME        convective plume mass / layer mass
!@var FCTYPE        fraction of total convective mass for each convective cloud type
!@var FDDRT         fraction of precipitating condensate avaliable for re-evaporation
!@var FENTR         entrained mass / convective mass
!@var FENTRA        entrained mass / layer mass
!@var FEVAP         fraction of layer mass available for precip evaporation
!@var FLEFT         fraction of plume after removing downdraft mass
!@var FQCOND        fraction of water vapor that condenses in plume
!@var FQCONDV       fraction of condensate that is lofted
!@var FQEVP         fraction of water vapor that evaporates in downdraft
!@var FPRCP         fraction of evaporated precipitation
!@var FSEVP         fraction of energy lost to evaporate
!@var FSSUM         fraction of energy lost to evaporate
!@var HDEP
!@var HPBL, PBLM    PBL height (m) and air mass in PBL (mb)
!@var LHX           latent heat of evaporation or sublimation (J/Kg)
!@var MPLUME,MPLUM1 mass of convective plume (mb)
!@var MCLOUD        air mass available for re-evaporation of precip (mb)
!@var MPMAX         mass of convective plume at the detrainment level (mb)
!@var MPOLD
!@var MNdO,MNdL,MNdI,MCDNCW,MCDNCI      Menon Stuff
!@var QMPMAX,SMPMAX values of QMP, SMP in detrained air
!@var QSATC         saturation vapor mixing ratio
!@var QSATMP        plume's saturation vapor mixing ratio
!@var RCLD,RCLDE    cloud particle radius, effective radius (microns)
!@var RHO           air density
!@var SMP, QMP      plume's SM, QM
!@var SLH           LHX/SHA
!@var TADJ          adjustment time for stablization of cloud base by cumulus mass flux
!@var TP            plume's temperature (K)
!@var TVP
!@var TOLD,TOLD1    old temperatures
!@var VT            precip terminal velocity (m/s)
!@var WMIX
!@var W2TEM
!@var WORK          work done on convective plume
!@var WMAX specified maximum convective updraft speed (m/s)
!@var WV convective updraft speed (m/s)
    !
    !          *******    VARIABLES DEFINED FOR     *******
    !          *******  CONDITIONAL COMPILIATIONS   *******
    !
#ifdef TRACERS_ON
!@var TMOLD: old TM (tracer mass)
    real*8, dimension(LM,NTM)      :: TMOLD, TM1
    real*8, dimension(NMOM,LM,NTM) :: TMOMOLD
    real*8, dimension(NTM) :: TMP, TMPMAX, TENV, TMDN, TM_dum, DTR
    real*8, dimension(NMOM,NTM) :: TMOMP, TMOMPMAX, TMOMDN
    real*8 :: vsum
    integer :: lborrow1
#ifdef TRACERS_WATER
!@var TRPCRP tracer mass in precip
    real*8, dimension(NTM)      :: TRPRCP
!@var FQCONDT fraction of tracer that condenses
!@var FQEVPT  fraction of tracer that evaporates (in downdrafts)
!@var FPRCPT fraction of tracer that evaporates (in net re-evaporation)
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
    real*8 :: FQCONDT(NTM), FWASHT(NTM), FPRCPT(NTM), FQEVPT(NTM)
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var precip_mm precipitation (mm) from the grid box above for washout
    real*8 WMXTR, precip_mm
    ! for tracers in general, added by Koch
    real*8, dimension(NTM) :: THLAW,THWASH,TR_LEF,TMFAC,TR_LEFT
    real*8 CLDSAVT
    integer :: IGAS
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var TMFAC used to adjust tracer moments
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
    real*8 TMP_SUL(LM,NTM)
    ! for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
    real*8 WA_VOL
    real*8, dimension(NTM) ::SULFOUT,SULFIN,SULFINC
    integer :: IAQCH
#endif
    real*8 HEFF
#endif
#endif
    !
#ifdef CLD_AER_CDNC
    real*8, dimension(LM) ::  CONDPC
#endif
    !
#ifdef CLD_AER_CDNC
    real*8 &
         MCDNO1,MCDNL1,CDNCB,fcnv,ATEMP,VVEL &
         ,DSGL(LM,SNTM),DSS(SNTM) &
         ,Repsis,Repsi,Rbeta,RCLD_C &
         ,AIRM_CDNC,TM_CDNC(NTM) &
         ,MNdO_max(LM),MNdL_max(LM) &
         ,MNdO_min(LM),MNdL_min(LM)
#ifdef TRACERS_AMP
    real*8                    :: ncaero (nmodes)
    integer                   ::nm
#endif
#endif
#ifdef CLD_AER_CDNC
!@var MCDNCW,MCDNCI cloud droplet # for warm,cold moist conv clouds (cm^-3)
#endif
#ifdef CLD_AER_CDNC
    !     REAL*8 RHO   ! air density
    !CN0 is the No parameter in the Marshall-Palmer distribution
#endif
    !
    !          *******                                          *******
    !          *******         END DECLARATION SECTION          *******
    !          *******                                          *******
    !
    !
    !****
    !**** MOIST CONVECTION
    !****
    !**** CONVECTION USES A MASS FLUX APPROACH WITH THE MASS PARTITIONED
    !**** INTO TWO PLUMES WITH DIFFERENT ENTRAINMENT RATES.  CONVECTION
    !**** IS INITIATED AT THE LOWEST MOIST CONVECTIVELY UNSTABLE LEVEL
    !**** AND TERMINATES AT THE LEVEL AT WHICH THE UPDRAFT SPEED IS NO
    !**** LONGER POSITIVE.  FURTHER CONVECTION MAY THEN ARISE FROM
    !**** SUCCESSIVELY HIGHER BASE LEVELS.
    !****
    !**** THE PARAMETERIZATION CONTAINS 7 PRIMARY PHYSICS LOOPS:
    !**** (1) MAJOR OUTER LOOP OVER SUCCESSIVE POTENTIAL CLOUD BASE LAYERS IN
    !****     WHICH THE CUMULUS MASS FLUX AND ITS PARTITIONING ARE CALCULATED.
    !**** (2) MAJOR INNER LOOP OVER THE TWO CLOUD TYPES CONTAINING MOST OF THE
    !****     OTHER PHYSICS LOOPS.
    !**** (3) PARCEL ASCENT LOOP IN WHICH CONVECTION IS TRIGGERED, WATER IS
    !****     CONDENSED, AIR IS ENTRAINED AND DETRAINED, DOWNDRAFT INITIATION
    !****     LEVELS ARE DEFINED, UPDRAFT SPEED IS CALCULATED, CONDENSATE IS
    !****     PARTITIONED INTO CLOUD AND PRECIPITATION, AND ASCENT IS TERMINATED.
    !**** (4) DOWNDRAFT DESCENT LOOP.
    !**** (5) COMPENSATING ENVIRONMENTAL SUBSIDENCE LOOP.
    !**** (6) PRECIPITATION AND RE-EVAPORATION LOOP, INCLUDING CONVECTIVE CLOUD
    !****     FRACTION CALCULATION.
    !**** (7) CONVECTIVE CLOUD OPTICAL THICKNESS LOOP, THE ONLY ONE OUTSIDE
    !****     THE MAJOR OUTER LOOP.
    !****
    ierr=0 ; lerr=0
    LMCMIN=0
    LMCMAX=0
    MCCONT=0
    FMC1=0.
    FSSL=1.
    RCLDX=radius_multiplier
    !**** SET USER-CONTROLLED FREE PARAMETERS
    TADJ = 1.d0
    !     PGRAD = 0.7d0
    PGRAD = 0.7                  ! to maintain bit compatibility
    CONTCE1=entrainment_cont1
    CONTCE2=entrainment_cont2
    FDDRT = 1.d0
    !**** initiallise arrays of computed output
    TAUMCL=0
    SVWMXL=0
    SVLATL=0
    SVLAT1=0
    VSUBL=0
    PRECNVL=0
    CLDMCL=0
    CLDSLWIJ=0
    CLDDEPIJ=0
    PRCPMC=0.
    TPSAV=0
    CSIZEL=RWCLDOX*10.*(1.-PEARTH)+10.*PEARTH ! droplet rad in stem
    VLAT=LHE
    SAVWL=0.
    SAVWL1=0.
    LHP=0
#ifdef SCM
    if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
      WCUDEEP=0. ; MPLUMEDEEP=0. ; ENTDEEP=0.
      DETRAINDEEP = 0. ; PRCCDEEP = 0. ; NPRCCDEEP = 0.
      TPALL = 0. ; MCCOND = 0. ; PRCCGRP = 0. ; PRCCICE = 0.
    endif
#endif
    SAVE1L=0.
    SAVE2L=0.
#ifdef TRACERS_WATER
    trsvwml = 0.
    TRPRCP = 0.
    TRPRMC = 0.
#ifdef TRDIAG_WETDEPO
    if (diag_wetdep == 1) then
      !**** initialize diagnostic arrays
      trcond_mc=0.D0
      trdvap_mc=0.D0
      trflcw_mc=0.D0
      trprcp_mc=0.D0
      trnvap_mc=0.D0
      trwash_mc=0.D0
    end if
#endif
#endif
    !**** zero out diagnostics
    MCFLX =0.
    DGDSM=0.
    DGDEEP=0.
    DGSHLW=0.
    DPHASE=0.
    DPHADEEP=0.
    DPHASHLW=0.
    DTOTW=0.
    DQCOND=0.
    DGDQM=0.
    DQMTOTAL=0.
    DQMSHLW=0.
    DQMDEEP=0.
    DQCTOTAL=0.
    DQCSHLW=0.
    DQCDEEP=0.
    DQLSC=0.
    DDMFLX=0.
    TDNL=0.
    QDNL=0.
#ifdef SCM
    if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
      CUMFLX=0.
      DWNFLX=0.
    endif
#endif
    !**** save initial values (which will be updated after subsid)
    SM1=SM
    QM1=QM
#ifdef TRACERS_ON
    TM1(:,1:NTX) = TM(:,1:NTX)
    TRDNL = 0.
#ifdef TRACERS_WATER
    CLDSAVT=0.
#endif
#endif
    !**** SAVE ORIG PROFILES
    SMOLD(:) = SM(:)
    SMOMOLD(:,:) = SMOM(:,:)
    QMOLD(:) = QM(:)
    QMOMOLD(:,:) = QMOM(:,:)
#ifdef TRACERS_ON
    TMOLD(:,1:NTX) = TM(:,1:NTX)
    TMOMOLD(:,:,1:NTX) = TMOM(:,:,1:NTX)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
    DT_SULF_MC(1:NTM,:)=0.
#endif
#ifdef TRACERS_WATER
    ! TR_LEF is an input to get_cond_factor not currently used for MC clouds
    TR_LEF(:)=1.D0
    thlaw(:) = 0.  ! nonzero only for gas tracers
    thwash(:) = 0. ! nonzero only for gas tracers
    tmfac(:) = 0.  ! nonzero only for gas tracers
    fwasht(:) = 0. ! nonzero only for aerosols
    fqcondt(:) = 0.
#endif
#endif
    !**** CALULATE PBL HEIGHT AND MASS
    PBLM=0.
    HPBL=0.
    do L=1,DCL
      if(L.lt.DCL) then
        PBLM=PBLM+AIRM(L)
        HPBL=HPBL+AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
      else
        PBLM=PBLM+.5d0*AIRM(L)
        HPBL=HPBL+.5d0*AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
      end if
    end do
    BYPBLM=1.d0/PBLM
    !**** CALCULATE THRESHOLD RH FOR LARGE-SCALE CLOUD FORMATION IN PBL
    !**** BASED ON SIEBESMA ET AL. (2003, JAS)
    do L=1,LM
      U00L(L)=0.
      if(PL(L).ge.850.d0) then
        LHX=LHE
        !         IF(TL(L).LT.TF) LHX=LHS ! use 10C
        U00L(L)=1.d0-2.*(U00b*.001*.050*(HPBL/500.)* &
             !    *         (.001*SQRT(DXYPIJ))**.33)/QSAT(TL(L),LHX,PL(L))
             (.001*sqrt(DXYPIJ))**.33)/QSAT(283.16d0,LHX,PL(L))
      end if
    end do
    !**** CALCULATE DEL WCU TO TRAVEL HALF LAYER THICKNESS IN ONE
    !**** TIMESTEP; USED LATER TO DETERMINE FRACTION OF CLOUD WATER
    !**** THAT DETRAINS AT CURRENT LAYER
    DWCU=0.
    do L=1,LMCM
      DWCU=DWCU+AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
    end do
    DWCU=0.5*DWCU*BYDTsrc/real(LMCM)

#ifdef CLD_AER_CDNC
    MNdO_max(:)=teeny
    MNdL_max(:)=teeny
    MNdO_min(:)=teeny
    MNdL_min(:)=teeny
#endif

    !****
    !**** BEGIN OUTER LOOP (1) OVER BASE LAYERS
    !****

    CLOUD_BASE: do LMIN=1,LMCM-1
      MAXLVL=0
      MINLVL=LM

      !****
      !**** THE MASS FLUX CLOSURE CALCULATES THE MASS REQUIRED TO RESTORE THE
      !**** CLOUD BASE LEVEL TO NEUTRAL BUOYANCY AFTER SUBSIDENCE (YAO AND
      !**** DEL GENIO 1989, J. CLIM.).  THE RESULTING CUMULUS MASS FLUX IS
      !**** ASSUMED TO OCCUR OVER A CONVECTIVE ADJUSTMENT TIME OF 1 HOUR
      !****

      !**** COMPUTE THE CONVECTIVE MASS OF THE LESS ENTRAINING PART BASED ON THE
      !**** LARGE-SCALE VERTICAL VELOCITY AT CLOUD BASE
      FMP0=-10.*CK1*SDL(LMIN+1)*BYGRAV*XMASS
      if(FMP0.le.0.) FMP0=0.

      !**** CREATE A PLUME IN THE BOTTOM LAYER
      !****
      !**** DETERMINE IF CLOUD BASE LEVEL IS UNSTABLE
      !**** ITERATION TO FIND FPLUME WHICH RESTORES CLOUD BASE TO NEUTRAL STATE
      !****
      SMO1=SM(LMIN)
      QMO1=QM(LMIN)
      SMO2=SM(LMIN+1)
      QMO2=QM(LMIN+1)
      SDN=SMO1*BYAM(LMIN)
      SUP=SMO2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMO1*BYAM(LMIN)
      QUP=QMO2*BYAM(LMIN+1)
      WMDN=WML(LMIN)
      WMUP=WML(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      WMEDG=.5*(WMUP+WMDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      LHX=LHE
      SLH=LHX*BYSHA
      DMSE=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+ &
           SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      if(DMSE.gt.-1d-10) cycle  ! try next level

      !**** MASS_FLUX PERFORMS THE ITERATIONS
      !**** COMPUTES FPLUME, FMP2
      call MASS_FLUX (FPLUME, FMP2, DQSUM, LMIN, &
           LHX, QMO1, QMO2, SLH, SMO1, SMO2, WMDN, WMUP, WMEDG)

      if(FPLUME.le..001) cycle ! try next level

      !****
      !**** BEGIN LOOP (2) THROUGH CLOUD TYPES OF DIFFERENT ENTRAINMENT RATE
      !****

      ITYPE=2           ! always 2 types of clouds: less and more entraining
      FCTYPE=1.

      !**** STABILIZATION IS ASSUMED TO OCCUROVER 1 HOUR, SO ONLY APPLY A FRACTION
      !**** OF THE REQUIRED MASS FLUX IN ONE PHYSICS TIMESTEP
      FMP2=FMP2*min(1d0,DTsrc/(TADJ*SECONDS_PER_HOUR))
      WMAX=50.
      !****
      CLOUD_TYPES:  do IC=1,ITYPE
        !****
        !**** Initialise plume characteristics
        MC1=.false.    ! flag for first convection event
        LHX=LHE
        MPLUME=min(AIRM(LMIN),AIRM(LMIN+1))
        if(MPLUME.gt.FMP2) MPLUME=FMP2

        if(ITYPE.eq.2) then     ! cal. MPLUME for 1st plume and 2nd plume
          FCTYPE=1.
          if(MPLUME.gt.FMP0) FCTYPE=FMP0/MPLUME
          if(IC.eq.2) FCTYPE=1.-FCTYPE
          if(FCTYPE.lt.0.001) cycle CLOUD_TYPES
        end if
        MPLUM1=MPLUME

        !****
        !**** Calculate gridbox area partition between (moist convection + subsidence)
        !**** and (large-scale cloud + clear sky).  A first calculation of the subsidence
        !**** area is performed.  If the convection rises 2 layers or more (MCCONT > 1)
        !**** subsidence is assumed to occur over a greater area relative to the updraft
        !****
        !**** To test code with only one estimate of area partition
        !**** remove NPPL Loop, or set NPPL=1,1
        !****
        AREA_PARTITION: do NPPL=1,2

          if(NPPL.eq.2)   then                  ! Perform 2nd pass if MC
            if(.not.MC1 .or. MCCONT.lt.2)  then ! went more than 2 layers
              cycle AREA_PARTITION              ! else - skip 2nd pass
            else
#ifdef TRACERS_ON
              call reset_tracer_work_arrays(lmin,lmax)
#endif
            end if
          end if

          MPLUME=MPLUM1*FCTYPE

          !**** Calculate fraction for plume either for the first time, or if
          !**** we are redoing the calculation for deeper convection
          if (MCCONT.eq.0  .or. MC1) then

            if (MCCONT.eq.0) then   ! first time
              FSUB_tmp=1.d0+(AIRM(LMIN+1)-100.d0)/200.d0
            else                    ! redoing plume calc
              FSUB_tmp=1.d0+(PL(LMIN)-PL(LMAX)-100.d0)/200.d0
            end if

            FCONV_tmp=min(MPLUM1*BYAM(LMIN+1),1d0)
            if(FSUB_tmp.gt.1.d0/(FCONV_tmp+1.d-20)-1.d0) &
                 FSUB_tmp=1.d0/(FCONV_tmp+1.d-20)-1.d0
            FSUB_tmp=max(1.d0,min(FSUB_tmp,5.d0))

            FSSL_tmp=1.d0-(1.d0+FSUB_tmp)*FCONV_tmp
            FSSL_tmp=max(CLDMIN,min(FSSL_tmp,1d0-CLDMIN))

            FMC1=(1.d0-FSSL_tmp)+teeny
          end if

          !**** guard against possibility of too big a plume
          if (MC1 .or. MCCONT.gt.0) then
            MPLUME=min(0.95d0*AIRM(LMIN)*FMC1,MPLUME)
          end if

          do L=1,LM
            COND(L)=0.    ;    CDHEAT(L)=0. !;  VLAT(L)=LHE
            CONDP(L)=0.   ;  CONDP1(L)=0. ;    CONDGP(L)=0.
            CONDIP(L)=0.  ;  CONDV(L)=0.  ;    HEAT1(L)=0.
            DM(L)=0.      ;  DMR(L)=0.    ;    DDR(L)=0.
            CCM(L)=0.     ;  DDM(L)=0.    ;    TAUMC1(L)=0.
            ENT(L)=0.     ;  DET(L)=0.    ;    BUOY(L)=0.
            WCU(L)=0.     ; SMDNL(L)=0.   ;  QMDNL(L)=0.
#ifdef SCM
            if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
              WCUALL(L,IC,LMIN) = 0.
              MPLUMEALL(L,IC,LMIN) = 0.
              ENTALL(L,IC,LMIN) = 0.
            endif
#endif
          end do
          SMOMDNL(:,:)=0.   ;  QMOMDNL(:,:)=0.
          UMDN(1:KMAX)=0.   ;  VMDN(1:KMAX)=0.
          UMDNL(1:KMAX,:)=0.  ;  VMDNL(1:KMAX,:)=0.
          DUM(1:KMAX,:)=0.    ;  DVM(1:KMAX,:)=0.
          DSM(:) = 0. ; DSMOM(:,:) = 0. ; DSMR(:) = 0. ; DSMOMR(:,:) = 0.
          DQM(:) = 0. ; DQMOM(:,:) = 0. ; DQMR(:) = 0. ; DQMOMR(:,:) = 0.
          !#ifdef TRACERS_ON
          ! (re)zeroing now done post-calculation with calls to
          ! reset_tracer_work_arrays
          !      DTM(:,1:NTX) = 0.   ;   DTMOM(:,:,1:NTX) = 0.
          !      DTMR(:,1:NTX) = 0.  ;  DTMOMR(:,:,1:NTX) = 0.
          !      TMDNL(:,1:NTX) = 0.  ;  TMOMDNL(:,:,1:NTX) = 0.
          !      TPOLD = 0.
          !#endif
          !#ifdef TRACERS_WATER
          !      TRCOND = 0.  ; TRCONDV = 0.
          !#endif

          !**** adjust MPLUME to take account of restricted area of subsidence
          !**** (i.e. MPLUME is now a greater fraction of the relevant airmass.
          MPLUME=min( MPLUME/FMC1, &
               AIRM(LMIN)*0.95d0*QM(LMIN)/(QMOLD(LMIN) + teeny) )
          if(MPLUME.le..001*AIRM(LMIN)) cycle CLOUD_TYPES
          FPLUME=MPLUME*BYAM(LMIN)
          SMP  =  SMOLD(LMIN)*FPLUME
          SMOMP(xymoms)=SMOMOLD(xymoms,LMIN)*FPLUME
          QMP  =  QMOLD(LMIN)*FPLUME
          QMOMP(xymoms)=QMOMOLD(xymoms,LMIN)*FPLUME
          if (TPSAV(LMIN).eq.0) TPSAV(LMIN)=SMP*PLK(LMIN)/MPLUME
          DMR(LMIN)=-MPLUME
          DSMR(LMIN)=-SMP
          DSMOMR(xymoms,LMIN)=-SMOMP(xymoms)
          DSMOMR(zmoms,LMIN)=-SMOMOLD(zmoms,LMIN)*FPLUME
          DQMR(LMIN)=-QMP
          DQMOMR(xymoms,LMIN)=-QMOMP(xymoms)
          DQMOMR(zmoms,LMIN)=-QMOMOLD(zmoms,LMIN)*FPLUME
#ifdef TRACERS_ON
          !**** This is a fix to prevent very occasional plumes that take out
          !**** too much tracer mass. This can impact tracers with very sharp
          !**** vertical gradients
          do n=1,ntx
            TMP(n) = TMOLD(LMIN,n)*FPLUME
            if(t_qlimit(n)) TMP(n) = min(TMP(n),0.95d0*TM(LMIN,n))
          enddo
!TOMAS DEBUG
            DO N=1,NTM
              if(TMP(n).lt.0.) print*,'TMP<0 1',TMP(n),trname(n)
            ENDDO
!TOMAS DEBUG

          TMOMP(xymoms,1:NTX)=TMOMOLD(xymoms,LMIN,1:NTX)*FPLUME
          DTMR(LMIN,1:NTX)=-TMP(1:NTX)
          DTMOMR(xymoms,LMIN,1:NTX)=-TMOMP(xymoms,1:NTX)
          DTMOMR( zmoms,LMIN,1:NTX)=-TMOMOLD(zmoms,LMIN,1:NTX)*FPLUME
          TPOLD(LMIN)=TPSAV(LMIN)  ! initial plume temperature
#endif
          do K=1,KMAX
            UMP(K)=UM(K,LMIN)*FPLUME
            DUM(K,LMIN)=-UMP(K)
            VMP(K)=VM(K,LMIN)*FPLUME
            DVM(K,LMIN)=-VMP(K)
          end do
          !****
          !**** RAISE THE PLUME TO THE TOP OF CONVECTION AND CALCULATE
          !**** ENTRAINMENT, CONDENSATION, AND SECONDARY MIXING
          !****
          CDHSUM=0.
          CDHSUM1=0.
          CDHDRT=0.
          ETADN=0.
          LDRAFT=LM
          EVPSUM=0.
          DDRAFT=0.
          LFRZ=0
          LMAX=LMIN
          CONTCE=CONTCE1
          if(IC.eq.2) CONTCE=CONTCE2
          !**** INITIAL CUMULUS UPDRAFT SPEED DETERMINED FROM SQRT(TKE), BUT
          !**** THE LARGER OF 0.5 M/S OR THE TURBULENT VERTICAL VELOCITY FOR
          !**** THE MORE ENTRAINING PLUME; THE LARGER OF 0.5 M/S OR TWICE THE
          !**** TURBULENT VERTICAL VELOCITY FOR THE LESS ENTRAINING PLUME
          WCU(LMIN)=max(.5D0,WTURB(LMIN))
          if(IC.eq.1) WCU(LMIN)=max(.5D0,2.D0*WTURB(LMIN))
          WCU2(LMIN)=WCU(LMIN)*WCU(LMIN)

          !****
          !**** BEGIN LOOP (3) OVER POSSIBLE CLOUD TOP LEVELS
          !****

          CLOUD_TOP:  do L=LMIN+1,LM

            !****
            !**** TRIGGERING CONDITIONS FOR MOIST CONVECTION
            !****

            !**** (1) TEST WHETHER MASS OF AIR REQUIRED TO STABILIZE CLOUD BASE
            !**** LARGE ENOUGH TO WARRANT PERFORMING CALCULATIONS
            if(MPLUME.le..001*AIRM(L)) exit CLOUD_TOP

            !**** (2) TEST WHETHER VIRTUAL MOIST STATIC ENERGY OF PARCEL LIFTED
            !**** TO NEXT LEVEL EXCEEDS THAT OF ENVIRONMENT AT THAT LEVEL, I.E.,
            !**** WHETHER PARCEL IS BUOYANT
            SDN=SMP/MPLUME
            SUP=SM1(L)*BYAM(L)
            QDN=QMP/MPLUME
            QUP=QM1(L)*BYAM(L)
            WMDN=0.
            WMUP=WML(L)
            SVDN=SDN*(1.+DELTX*QDN-WMDN)
            SVUP=SUP*(1.+DELTX*QUP-WMUP)
            if(PLK(L-1)*(SVUP-SVDN)+SLHE*(QUP-QDN).ge.0.) exit CLOUD_TOP
            !**** Only set TPSAV if it hasn't already been set
            if (TPSAV(L).eq.0) TPSAV(L)=SMP*PLK(L)/MPLUME
            TP=TPSAV(L)
            if(TPSAV(L-1).ge.TF.and.TPSAV(L).lt.TF) LFRZ=L-1

            !**** (3)TEST TO SEE WHETHER LIFTED PARCEL IS ABOVE LIFTING CONDENSATION
            !**** W.R.T. LIQUID WATER, OR W.R.T. ICE FOR HOMOGENEOUS NUCLEATION (T<-40)
            LHX=LHE
            if(TP.lt.TI) LHX=LHS
            QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
            if(QMP.lt.QSATMP) exit CLOUD_TOP    !  no plume

            if(TP.lt.TF.and.LHX.eq.LHE) then
              LHX=LHS
              QSATMP=MPLUME*QSAT(TP,LHX,PL(L))
            end if

            !**** DEFINE DUMMY LATENT HEAT VARIABLE TO AVOID PHASE DISCREPANCY BETWEEN PLUMES
            if (VLAT(L).eq.LHS) LHX=LHS
            VLAT(L)=LHX
            SLH=LHX*BYSHA

            !**** THRESHOLD RH FOR PBL STRATIFORM CLOUDS
            if(TL(L).ge.TF .and. U00L(L).ne.U00a) U00L(L)= &
                 1.d0-2.*(U00b*2.d-4*MPLUME/QSATMP)

            MCCONT=MCCONT+1
            if(MCCONT.eq.1) MC1=.true.
            !****
            !**** IF PLUME MASS IS TOO LARGE FOR UPPER LAYER, LEAVE PART BEHIND IN LOWER LAYER
            !**** AND ADJUST TEMPERATURE, HUMIDITY, AND MOMENTUM THERE
            !****
            if(MPLUME.gt..95*AIRM(L)) then
              DELTA=(MPLUME-.95*AIRM(L))/MPLUME
              DM(L-1)=DM(L-1)+DELTA*MPLUME
              MPLUME=.95*AIRM(L)

              DSM(L-1)=  DSM(L-1)+DELTA*SMP
              SMP = SMP  *(1.-DELTA)
              DSMOM(xymoms,L-1)=DSMOM(xymoms,L-1)+DELTA*SMOMP(xymoms)
              SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)

              DQM(L-1)=  DQM(L-1)+DELTA*QMP
              QMP = QMP  *(1.-DELTA)
              DQMOM(xymoms,L-1)=DQMOM(xymoms,L-1)+DELTA*QMOMP(xymoms)
              QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)

              do K=1,KMAX
                DUM(K,L-1)=DUM(K,L-1)+UMP(K)*DELTA
                DVM(K,L-1)=DVM(K,L-1)+VMP(K)*DELTA
                UMP(K)=UMP(K)-UMP(K)*DELTA
                VMP(K)=VMP(K)-VMP(K)*DELTA
              end do

#ifdef TRACERS_ON
              DTM(L-1,1:NTX) = DTM(L-1,1:NTX)+DELTA*TMP(1:NTX)
              DTMOM(xymoms,L-1,1:NTX)=DTMOM(xymoms,L-1,1:NTX)+DELTA &
                   *TMOMP(xymoms,1:NTX)
              TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
              TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
!TOMAS DEBUG
            DO N=1,NTM
              if(TMP(n).lt.0.) print*,'TMP<0 2',TMP(n),trname(n),DELTA
            ENDDO
!TOMAS DEBUG

#endif
            end if

            !**** WORK DONE BY CONVECTION IN UPPER LAYER REMOVES ENERGY FROM THE PLUME
            WORK=MPLUME*(SUP-SDN)*(PLK(L-1)-PLK(L))/PLK(L-1)
            DSM(L-1)=DSM(L-1)-WORK
            CCM(L-1)=MPLUME

            !**** CONDENSE VAPOR IN THE PLUME AND ADD LATENT HEAT
            call get_dq_cond(smp,qmp,plk(l),mplume,lhx,pl(l),dqsum,fqcond)

            if(DQSUM.gt.0. .and. QMP.gt.teeny) then
              QMOMP(xymoms) =  QMOMP(xymoms)*(1.-FQCOND)
              SMP=SMP+SLH*DQSUM/PLK(L)
              QMP=QMP-DQSUM
            end if

#ifdef TRACERS_ON
            !**** save plume temperature after possible condensation
            TPOLD(L)=SMP*PLK(L)/MPLUME
#endif

            !**** TOTAL CONDENSATE IN LAYER = NEWLY CONDENSED + ADVECTED FROM LOWER LAYER;
            !**** THE LATTER IS CALCULATED IN THE MICROPHYSICS SECTION
            COND(L)=DQSUM
            CDHEAT(L)=SLH*COND(L)          ! calculate CDHEAT before add CONDV
            CDHSUM=CDHSUM+CDHEAT(L)
            COND(L)=COND(L)+CONDV(L-1)     ! add in the vertical transported COND

            if (VLAT(L-1).ne.VLAT(L)) then
              SMP=SMP-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA/PLK(L) ! phase change of uplifted condensate
              CDHEAT(L)=CDHEAT(L)-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA
              CDHSUM=CDHSUM-(VLAT(L-1)-VLAT(L))*CONDV(L-1)*BYSHA
            end if

            !**** CALCULATE SLOPE PARAMETER FOR MARSHALL-PALMER SIZE DISTRIBUTION FOR
            !**** LIQUID, GRAUPEL, ICE (DEL GENIO ET AL. 2005, J. CLIM.)
            CONDMU=100.*COND(L)*PL(L)/(CCM(L-1)*TL(L)*RGAS)
            FLAMW=(1000.d0*PI*CN0/(CONDMU+teeny))**.25
            FLAMG=(400.d0*PI*CN0G/(CONDMU+teeny))**.25
            FLAMI=(100.d0*PI*CN0I/(CONDMU+teeny))**.25

#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
!@auth Menon  saving aerosols mass for CDNC prediction
            do N=1,SNTM
              DSS(N)=1.d-10
              DSGL(L,N)=1.d-10
            enddo
#endif

#ifdef CLD_AER_CDNC
#ifdef ALT_CDNC_INPUTS
            ! aerosols in the updraft
            tm_cdnc(:) = tmp(:)
            airm_cdnc = mplume
#else
            ! ambient aerosols at this level
            tm_cdnc(:) = tm(l,:)
            airm_cdnc = airm(l)
#endif
#endif
            !**** Here we change convective precip due to aerosols
#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
            do N=1,NTX
              select case (trname(ntix(n)))
              case('SO4')
                DSGL(L,1)=tm_cdnc(n)     !n=19
                DSS(1) = DSGL(L,1)
              case('seasalt1')
                DSGL(L,2)=tm_cdnc(n)     !n=21
                DSS(2) = DSGL(L,2)
              case('seasalt2')
                DSGL(L,3)=tm_cdnc(n)     !n=22
                DSS(3) = DSGL(L,3)
              case('OCIA')
                DSGL(L,4)=tm_cdnc(n)     !n=27
                DSS(4) = DSGL(L,4)
              case('OCB')
                DSGL(L,5)=tm_cdnc(n)     !n=28
                DSS(5) = DSGL(L,5)
              case('BCIA')
                DSGL(L,6)=tm_cdnc(n)     !n=24
                DSS(6) = DSGL(L,6)
              case('BCB')
                DSGL(L,7)=tm_cdnc(n)     !n=25
                DSS(7) = DSGL(L,7)
              case('OCII')
                DSGL(L,8)=tm_cdnc(n)     !n=26
                DSS(8) = DSGL(L,8)
              case('BCII')
                DSGL(L,9)=tm_cdnc(n)     !n=23
                DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
              case('Clay')
                DSGL(L,10)=tm_cdnc(n)    !n=23
                DSS(10) = DSGL(L,10)
              case('Silt1')
                DSGL(L,11)=tm_cdnc(n)    !n=23
                DSS(11) = DSGL(L,11)
              case('Silt2')
                DSGL(L,12)=tm_cdnc(n)    !n=23
                DSS(12) = DSGL(L,12)
              case('Silt3')
                DSGL(L,13)=tm_cdnc(n)    !n=23
                DSS(13) = DSGL(L,13)
#endif
#ifdef TRACERS_NITRATE
              case('NO3p')
                DSGL(L,14)=tm_cdnc(n)    !n=23
                DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
                !**** Here are dust particles coated with sulfate
              case('SO4_d1')
                DSGL(L,15)=tm_cdnc(n)    !n=20
                DSS(15) = DSGL(L,15)
              case('SO4_d2')
                DSGL(L,16)=tm_cdnc(n)    !n=21
                DSS(16) = DSGL(L,16)
              case('SO4_d3')
                DSGL(L,17)=tm_cdnc(n)    !n=22
                DSS(17) = DSGL(L,17)
#endif
#ifdef TRACERS_AEROSOLS_SOA
              case('isopp1a')
                DSGL(L,18)=tm_cdnc(n)
                DSS(18) = DSGL(L,18)
              case('isopp2a')
                DSGL(L,19)=tm_cdnc(n)
                DSS(19) = DSGL(L,19)
#ifdef TRACERS_TERP
              case('apinp1a')
                DSGL(L,20)=tm_cdnc(n)
                DSS(20) = DSGL(L,20)
              case('apinp2a')
                DSGL(L,21)=tm_cdnc(n)
                DSS(21) = DSGL(L,21)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
              case('OCocean')
                DSGL(L,22)=tm_cdnc(n)
                DSS(22) = DSGL(L,22)
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_AEROSOLS_VBS
              case('vbsAm2')
                DSGL(L,23)=tm_cdnc(n)
                DSS(23) = DSGL(L,23)
              case('vbsAm1')
                DSGL(L,24)=tm_cdnc(n)
                DSS(24) = DSGL(L,24)
              case('vbsAz')
                DSGL(L,25)=tm_cdnc(n)
                DSS(25) = DSGL(L,25)
              case('vbsAp1')
                DSGL(L,26)=tm_cdnc(n)
                DSS(26) = DSGL(L,26)
              case('vbsAp2')
                DSGL(L,27)=tm_cdnc(n)
                DSS(27) = DSGL(L,27)
              case('vbsAp3')
                DSGL(L,28)=tm_cdnc(n)
                DSS(28) = DSGL(L,28)
              case('vbsAp4')
                DSGL(L,29)=tm_cdnc(n)
                DSS(29) = DSGL(L,29)
              case('vbsAp5')
                DSGL(L,30)=tm_cdnc(n)
                DSS(30) = DSGL(L,30)
              case('vbsAp6')
                DSGL(L,31)=tm_cdnc(n)
                DSS(31) = DSGL(L,31)
#endif  /* TRACERS_AEROSOLS_VBS */
              end select
            end do      !end of n loop for tracers
#endif  /* (TRACERS_AEROSOLS_Koch) and (CLD_AER_CDNC) */
            !** Use MATRIX AMP_actv to decide what the aerosol number conc. is
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
#ifndef TRACERS_TOMAS
#ifdef TRACERS_AMP
            do nm=1,nmodes
              ncaero(nm)=naerc(l,nm)*1.d-6
              !         if(naerc(l,nm).gt.1.d-30) write(6,*)"mat",ncaero(nm),nm
            enddo
            call GET_CC_CDNC_MX(L,nmodes,ncaero,MCDNL1,MCDNO1)
#else
            !** This is for the old mass to number calculations nc. is
            call GET_CC_CDNC(L,AIRM_CDNC,DXYPIJ,PL(L),TL(L),DSS, &
                 MCDNL1,MCDNO1)

#endif  /* (TRACERS_AMP) */
#endif
#ifdef TRACERS_TOMAS
       CALL GET_CC_CDNC_TOMAS(L,I_debug,J_debug,AIRM_CDNC, &
          DXYPIJ,PL(L),TL(L),MCDNL1,MCDNO1)
#endif
            MNdO=MCDNO1
            MNdL=MCDNL1
            MNdO_max(L)=max(MNdO_max(L),MCDNO1)
            MNdL_max(L)=max(MNdL_max(L),MCDNL1)
            if(MNdO_min(L)==teeny) then
              MNdO_min(L)=MCDNO1
              MNdL_min(L)=MCDNL1
            else
              MNdO_min(L)=min(MNdO_min(L),MCDNO1)
              MNdL_min(L)=min(MNdL_min(L),MCDNL1)
            endif
            MNdI = 0.06417127d0
            MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
            MCDNCI=MNdI
            !      if(MCDNCW.gt.0.)write(6,*)"Mst CDNC,a",MCDNCW,MNdO,MNdL,L
            if (MCDNCW.le.20.d0) MCDNCW=20.d0     !set min CDNC, sensitivity test
            if (MCDNCW.ge.2000.d0) MCDNCW=2000.d0     !set max CDNC, sensitivity test
            !** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
            !** for spectral dispersion effects on droplet size distribution
            !** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
            Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
            Repsis=Repsi*Repsi
            Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
            RCLD_C=14.d0/Rbeta       !set Reff to threshold size =14 um (Rosenfeld)
#endif /* (CLD_AER_CDNC) */
            TAUMC1(L)=TAUMC1(L)+COND(L)*FMC1

#ifdef TRACERS_WATER
            !**** CONDENSING TRACERS
            WMXTR=DQSUM*BYAM(L)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
            WA_VOL=COND(L)*1.d2*BYGRAV*DXYPIJ

            if (FPLUME.gt.teeny) then
              TMP_SUL(L,aqchem_list)=TMP(aqchem_list)/FPLUME
            else
              TMP_SUL(L,aqchem_list)=0.
            end if

            call GET_SULFATE(L,TPOLD(L),FPLUME,WA_VOL,WMXTR,SULFIN, &
                 SULFINC,SULFOUT,TR_LEFT,TMP_SUL,TRCOND(1,L), &
                 AIRM,LHX,DT_SULF_MC(1,L),CLDSAVT)

            do iaqch=1,aqchem_count
              n = aqchem_list(iaqch)
              ! first apply chemistry
              ! removal of precursers
              TMP(N)=TMP(N)*(1.+SULFIN(N))
              TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.+SULFIN(N))
              ! formation of sulfate
              TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)
            enddo

#endif
            DO N=1,NTM
              if(TMP(n).lt.0.) print*,'TMP<0 3',TMP(n),trname(n)
            ENDDO

            TM_dum(:) = TMP(:)
            !     CALL GET_COND_FACTOR(
            !    &      NTX,WMXTR,TPOLD(L),TPOLD(L-1),LHX,FPLUME
            !    &     ,FQCOND,FQCONDT,.true.,TRCOND(:,L),TM_dum,THLAW,TR_LEF,PL(L)
            !    &     ,ntix,CLDSAVT)
!TOMAS DEBUG
            DO N=1,NTM
              if(TM_dum(n).lt.0.) print*,'TM_dum<0 1',TM_dum(n),trname(n)
            ENDDO
!TOMAS DEBUG
            call GET_COND_FACTOR( &
                 NTX,WMXTR,TPOLD(L),TPOLD(L-1),LHX,FPLUME &
                 ,FQCOND,FQCONDT,.true.,TRCOND(:,L),TM_dum,THLAW,TR_LEF,PL(L) &
                 ,ntix,FPLUME)
#ifdef TRACERS_AEROSOLS_OCEAN
            if (trm(i_debug,j_debug,L,n_ococean) .gt. 0.d0) then
              do n=1,ntx
                select case (trname(ntix(n)))
                case ('seasalt1', 'OCocean')
                  fqcondt(n)=fqcondt(n)* &
                       trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1) &
                       /(trm(i_debug,j_debug,L,n_ococean)/trpdens(n_ococean)+ &
                       trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1))
                end select
              enddo
            endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
            dtr(1:ntx) = fqcondt(1:ntx)*tmp(1:ntx)
#ifdef TRDIAG_WETDEPO
            if (diag_wetdep == 1) trcond_mc(l,1:ntx)=trcond_mc(l,1:ntx) &
                 +dtr(1:ntx)
#endif
            do N=1,NTX
              TRCOND(N,L) = DTR(N) + TRCOND(N,L) + TRCONDV(N,L-1)
              TMP(N)         = TMP(N)         *(1.-FQCONDT(N))
              TMOMP(xymoms,N)= TMOMP(xymoms,N)*(1.-FQCONDT(N))
            end do
#endif

            !****
            !**** ENTRAINMENT
            !****

            !**** ENTRAINMENT USES THE PARAMETERIZATION OF GREGORY (2001, QJRMS)
            !**** WITH CONSTANT C = 0.3 (0.6) FOR THE LESS (MORE) ENTRAINING PLUME
            !****
            TVP=(SMP/MPLUME)*PLK(L)*(1.+DELTX*QMP/MPLUME)
            BUOY(L)=(TVP-TVL(L))/TVL(L)-COND(L)/MPLUME
            ENT(L)=.16667D0*CONTCE*GRAV*BUOY(L)/(WCU(L-1)*WCU(L-1)+teeny)

            !**** WHEN PARCEL LOSES BUOYANCY TERMINATE ENTRAINMENT AND DETRAIN AN
            !**** AMOUNT OF AIR EQUAL TO THAT WHICH WOULD HAVE BEEN ENTRAINED HAD
            !**** THE BUOYANCY BEEN POSITIVE
            if (ENT(L).lt.0.D0) then
              DET(L)=-ENT(L)
              ENT(L)=0.D0
            end if

            if(ENT(L).gt.0.D0) then    ! non-zero entrainment
              FENTR=1000.D0*ENT(L)*GZL(L)*FPLUME
              if(FENTR+FPLUME.gt.1.) then
                FENTR=1.-FPLUME
                ENT(L)=0.001d0*FENTR/(GZL(L)*FPLUME)
              end if
              if(FENTR.ge.teeny) then    !  Big Enough, Proceed
                MPOLD=MPLUME
                FPOLD=FPLUME
                ETAL1=FENTR/(FPLUME+teeny)
                EPLUME=MPLUME*ETAL1
                !**** Reduce EPLUME so that mass flux is less than mass in box
                if (EPLUME.gt.AIRM(L)*0.975d0-MPLUME) then
                  EPLUME=AIRM(L)*0.975d0-MPLUME
                end if
                MPLUME=MPLUME+EPLUME
                ETAL1=EPLUME/MPOLD
                FENTR=ETAL1*FPOLD
                ENT(L)=0.001d0*FENTR/(GZL(L)*FPOLD)
                FPLUME=FPLUME+FENTR      ! to increase mass flux, remove this formula
                !     FPLUME = MPLUME*BYAM(L) ! and use this instead
                FENTRA = EPLUME*BYAM(L)
                DSMR(L)=DSMR(L)-EPLUME*SUP        ! = DSM(L)-SM(L)*FENTRA
                DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
                DQMR(L)=DQMR(L)-EPLUME*QUP        ! = DQM(L)-QM(L)*FENTRA
                DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
                DMR(L)=DMR(L)-EPLUME
                SMP=SMP+EPLUME*SUP
                SMOMP(xymoms)= SMOMP(xymoms)+ SMOM(xymoms,L)*FENTRA
                QMP=QMP+EPLUME*QUP
                QMOMP(xymoms)= QMOMP(xymoms)+ QMOM(xymoms,L)*FENTRA
#ifdef TRACERS_ON
                DTMR(L,1:NTX) = DTMR(L,1:NTX)-TM(L,1:NTX)*FENTRA
                DTMOMR(:,L,1:NTX) = DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
                TMP(1:NTX) = TMP(1:NTX)+TM(L,1:NTX)*FENTRA
                TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)+TMOM(xymoms,L,1:NTX) &
                     *FENTRA
#endif

                !****
                !**** CONVECTIVE MOMENTUM TRANSPORT IS BASED ON GREGORY ET AL. (1997, QJRMS).
                !**** THE MOMENTUM OF THE RISING PARCEL IS CHANGED BY ENTRAINMENT AND BY THE
                !**** EFFECT OF THE CONVECTIVE SCALE HORIZONTAL PRESSURE GRADIENT, THE
                !**** LATTER PARAMETERIZED AS PROPORTIONAL TO THE CUMULUS MASS FLUX TIMES THE
                !**** VERTICAL SHEAR OF THE HORIZONTAL WIND
                !****
                do K=1,KMAX
                  UMTEMP=PGRAD*MPLUME**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
                  VMTEMP=PGRAD*MPLUME**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
                  UMP(K)=UMP(K)+U_0(K,L)*EPLUME+UMTEMP
                  DUM(K,L)=DUM(K,L)-U_0(K,L)*EPLUME-UMTEMP
                  VMP(K)=VMP(K)+V_0(K,L)*EPLUME+VMTEMP
                  DVM(K,L)=DVM(K,L)-V_0(K,L)*EPLUME-VMTEMP
                end do

              end if         ! big enough check
            end if         ! non-zero entrainment check

            !****
            !**** DETRAINMENT ONLY OCCURS WHEN PARCEL LOSES BUOYANCY; MASS
            !**** CALCULATED ABOVE
            !****

            !**** CHANGE LAYER PROPERTIES DUE TO DETRAINMENT
            if(DET(L).gt.0.D0) then
              DELTA=1000.D0*DET(L)*GZL(L)
              if(DELTA.gt..95D0) then
                DELTA=.95d0
                DET(L)=.001d0*DELTA/GZL(L)
              end if
              DM(L)=DM(L)+DELTA*MPLUME
              MPLUME=MPLUME*(1.D0-DELTA)
              DSM(L)=  DSM(L)+DELTA*SMP
              SMP = SMP  *(1.-DELTA)
              DSMOM(xymoms,L)=DSMOM(xymoms,L)+DELTA*SMOMP(xymoms)
              SMOMP(xymoms) = SMOMP(xymoms)*(1.-DELTA)
              DQM(L)=  DQM(L)+DELTA*QMP
              QMP = QMP  *(1.-DELTA)
              DQMOM(xymoms,L)=DQMOM(xymoms,L)+DELTA*QMOMP(xymoms)
              QMOMP(xymoms) = QMOMP(xymoms)*(1.-DELTA)
              do K=1,KMAX
                UMTEMP=PGRAD*MPLUME**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
                VMTEMP=PGRAD*MPLUME**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
                DUM(K,L)=DUM(K,L)+UMP(K)*DELTA-UMTEMP
                DVM(K,L)=DVM(K,L)+VMP(K)*DELTA-VMTEMP
                UMP(K)=UMP(K)-UMP(K)*DELTA+UMTEMP
                VMP(K)=VMP(K)-VMP(K)*DELTA+VMTEMP
              end do

#ifdef TRACERS_ON
              DTM(L,1:NTX) = DTM(L,1:NTX)+DELTA*TMP(1:NTX)
              DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)+DELTA &
                   *TMOMP(xymoms,1:NTX)
              TMP(1:NTX) = TMP(1:NTX)*(1.-DELTA)
              TMOMP(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)*(1.-DELTA)
#endif
            end if

            !****
            !**** CHECK THE POSSIBILITY OF DOWNDRAFT.  A DOWNDRAFT FORMS WHEN AN
            !**** EQUAL MIXTURE OF CLOUDY AND CLEAR AIR IS NEGATIVELY BUOYANT
            !****
            !**** Define downdraft properties
            !****
            !**** DOWNDRAFT BEGINS 2 LEVELS ABOVE CLOUD BASE
            !****
            if(L-LMIN.gt.1) then

              SMIX=.5*(SUP+SMP/MPLUME)
              QMIX=.5*(QUP+QMP/MPLUME)
              SVMIX=SMIX             ! *(1.+DELTX*QMIX)
              SVUP=SUP               ! *(1.+DELTX*QUP)
              DMMIX=(SVUP-SVMIX)*PLK(L)
              !     +SLHE*(QSAT(SUP*PLK(L),LHX,PL(L))-QMIX)
              if(DMMIX.lt.1d-10) CDHDRT=CDHDRT+CDHEAT(L)

              !**** NO DOWNDRAFT IF BUOYANT
              if(DMMIX.ge.1d-10) then       ! the mixture is negatively buoyant

                LDRAFT=L                      ! the highest downdraft level
                !**** INITIATE DOWNDRAFT WITH SPECIFIED FRACTION OF UPDRAFT MASS
                ETADN=1.d0/3.d0    ! .20d0 ! reduce ETADN to improve computational stability

                !**** To test with code with no downdrafts, set etadn=0. here
                !**** etadn=0.  ! test

                FLEFT=1.-.5*ETADN
                DDRAFT=ETADN*MPLUME
                DDR(L)=DDRAFT
                CDHSUM1=CDHSUM1+CDHDRT*.5*ETADN      ! calculate before CDHDRT
                CDHDRT=CDHDRT-CDHDRT*.5*ETADN+CDHEAT(L)    ! SLH*COND(L)
                FDDP = .5*DDRAFT/MPLUME
                FDDL = .5*DDRAFT*BYAM(L)
                MPLUME=FLEFT*MPLUME
                SMDNL(L)=DDRAFT*SMIX
                SMOMDNL(xymoms,L)=SMOM(xymoms,L)*FDDL +  SMOMP(xymoms)*FDDP
                SMP=FLEFT*SMP
                SMOMP(xymoms)= SMOMP(xymoms)*FLEFT
                QMDNL(L)=DDRAFT*QMIX
                QMOMDNL(xymoms,L)=QMOM(xymoms,L)*FDDL +  QMOMP(xymoms)*FDDP
                QMP=FLEFT*QMP
                QMOMP(xymoms)= QMOMP(xymoms)*FLEFT
                DMR(L) = DMR(L)-.5*DDRAFT
                DSMR(L)=DSMR(L)-.5*DDRAFT*SUP        ! = DSM(L)-SM(L)*FDDL
                DSMOMR(:,L)=DSMOMR(:,L) - SMOM(:,L)*FDDL
                DQMR(L)=DQMR(L)-.5*DDRAFT*QUP        ! = DQM(L)-QM(L)*FDDL
                DQMOMR(:,L)=DQMOMR(:,L) - QMOM(:,L)*FDDL
#ifdef TRACERS_ON
                Tmdnl(l,1:NTX) = tm(l,1:NTX)*fddl+Tmp(1:NTX)*fddp
                tmomdnl(xymoms,l,1:NTX) = tmom(xymoms,l,1:NTX)*fddl+ &
                     tmomp(xymoms,  1:NTX)*fddp
                dtmr    (l,1:NTX) = dtmr    (l,1:NTX)-fddl *tm    (l,1:NTX)
                dtmomr(:,l,1:NTX) = dtmomr(:,l,1:NTX)-fddl *tmom(:,l,1:NTX)
                Tmp         (1:NTX) = Tmp         (1:NTX)*fleft
                tmomp(xymoms,1:NTX) = tmomp(xymoms,1:NTX)*fleft
#endif
                do K=1,KMAX
                  UMDNL(K,L)=.5*(ETADN*UMP(K)+DDRAFT*U_0(K,L))
                  UMP(K)=UMP(K)*FLEFT
                  DUM(K,L)=DUM(K,L)-.5*DDRAFT*U_0(K,L)
                  VMDNL(K,L)=.5*(ETADN*VMP(K)+DDRAFT*V_0(K,L))
                  VMP(K)=VMP(K)*FLEFT
                  DVM(K,L)=DVM(K,L)-.5*DDRAFT*V_0(K,L)
                end do

              end if    !   Buoyancy Test
            end if    !   More Than 1 Layer Check
            !****
            !**** COMPUTE CUMULUS UPDRAFT SPEED BASED ON GREGORY (2001, QJRMS)
            !****

            W2TEM=.16667D0*GRAV*BUOY(L)-WCU(L-1)*WCU(L-1)* &
                 (.66667D0*DET(L)+ENT(L))
            HDEP=AIRM(L)*TL(L)*RGAS/(GRAV*PL(L))
            WCU2(L)=WCU2(L-1)+2.*HDEP*W2TEM    ! use WCU2
            !     WCU(L)=WCU(L-1)+HDEP*W2TEM
            WCU(L)=0.
            if (WCU2(L).gt.0.D0) WCU(L)=sqrt(WCU2(L))
            if (WCU(L).ge.0.D0) WCU(L)=min(50.D0,WCU(L))
            if (WCU(L).lt.0.D0) WCU(L)=max(-50.D0,WCU(L))
#ifdef SCM
            if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
              !     save cumulus updraft speed and plume temperature
              WCUALL(L,IC,LMIN) = WCU(L)
              MPLUMEALL(L,IC,LMIN) = CCM(L-1)
              ENTALL(L,IC,LMIN) = 1000.D0*ENT(L)
              TPALL(L,IC,LMIN) = TP   ! Save plume temp ????? THIS IS NOT CORRECT
              !     write(iu_scm_prt,885) LMIN,ic,L,WCU(L), WCUALL(L,IC,LMIN)
              !885  format(1x,'mstcnv  lmin ic l wcu wcuall ',
              !    *      i5,i5,i5,2(f12.6))
            endif
#endif
            !**** UPDATE ALL QUANTITIES CARRIED BY THE PLUME
            !      SVLATL(L)=VLAT(L)
            SMPMAX=SMP
            SMOMPMAX(xymoms) =  SMOMP(xymoms)
            QMPMAX=QMP
            QMOMPMAX(xymoms) =  QMOMP(xymoms)
#ifdef TRACERS_ON
            !**** Tracers at top of plume
            TMPMAX(1:NTX) = TMP(1:NTX)
            TMOMPMAX(xymoms,1:NTX) = TMOMP(xymoms,1:NTX)
#endif
            MPMAX=MPLUME
            LMAX = LMAX + 1
            if(WCU2(L).lt.0.D0) exit CLOUD_TOP  !  use WCU(L)*WCU(L)

            !**** CUMULUS MICROPHYSICS
            WCUFRZ=0 ; if(LFRZ>0) WCUFRZ=WCU(LFRZ)
            call CONVECTIVE_MICROPHYSICS(PL(L),WCU(L),DWCU,LFRZ,WCUFRZ,TP, &
                 TI,FITMAX,PLAND,CN0,CN0I,CN0G,FLAMW,FLAMG,FLAMI,RHOIP,RHOG, &
                 ITMAX,TL(LMIN),TL(LMIN+1),WMAX, &
#ifdef CLD_AER_CDNC
                 TL(L),RCLD_C,MCDNCW,CONDMU,CONDPC(L), &
#endif
                 CONDP(L),CONDP1(L),CONDIP(L),CONDGP(L))

            !**** convert condp to the same units as cond
            CONDP(L)=.01d0*CONDP(L)*CCM(L-1)*TL(L)*RGAS/PL(L)
            CONDP1(L)=.01d0*CONDP1(L)*CCM(L-1)*TL(L)*RGAS/PL(L)

            !**** To test code with no precipitation from deep convection
            !**** set CONDP(L)=0. here.
            !**** To test code with no vertical advection of condensate
            !**** set CONDP1(L)=COND(L) here.

            !#ifdef SCM
            !     if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
            !         write(iu_scm_prt,339)  LMIN,IC,L,COND(L),CONDP1(L),CONDP(L),
            !    *               CONDGP(L),CONDIP(L)
            !339      format(1x,'LMIN IC L COND p1 p gp ip ',3(i5),5(f10.6))
            !     endif
            !#endif

            !**** CALCULATE VERTICALLY TRANSPORTED PART OF CONDENSATE AND REMOVE
            !**** FROM REMAINDER OF CONDENSATE (PRECIPITATING PART REMOVED LATER)
            if(CONDP1(L).gt.COND(L)) CONDP1(L)=COND(L)
            if(CONDP(L).gt.CONDP1(L)) CONDP(L)=CONDP1(L)
            CONDV(L)=COND(L)-CONDP1(L)       ! part of COND transported up
            COND(L)=COND(L)-CONDV(L)         ! CONDP1(L)

#ifdef TRACERS_WATER
            FQCONDV=CONDV(L)/((COND(L)+CONDV(L))+teeny)
            TRCONDV(:,L)=FQCONDV*TRCOND(:,L)
            TRCOND (:,L)=TRCOND(:,L)-TRCONDV(:,L)
#endif
          end do  CLOUD_TOP     !   End Loop 3  (L)

          if(LMIN.eq.LMAX) then
#ifdef TRACERS_ON
            call reset_tracer_work_arrays(lmin,lmin)
#endif
            cycle CLOUD_TYPES
          endif

        end do AREA_PARTITION    ! end area partition loop (NPPL)

        !**** UPDATE CHANGES CARRIED BY THE PLUME IN THE TOP CLOUD LAYER
        TAUMCL(LMIN:LMAX)=TAUMCL(LMIN:LMAX)+TAUMC1(LMIN:LMAX)

#ifdef SCM
        if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          if(PLE(LMIN)-PLE(LMAX+1).ge.450.) then
            !       write(iu_scm_prt,887) LMIN,LMAX,IC
            !887    format(1x,'mstcnv -- deep lmin lmax ic ',i5,i5,i5)
            do L=LMIN,LMAX
              WCUDEEP(L,IC) = WCU(L)
              MPLUMEDEEP(L,IC) = CCM(L-1)
              ENTDEEP(L,IC) = 1000.D0*ENT(L)
              if(IC.eq.1) then
                SAVWL(L)=WCU(L)
              else
                SAVWL1(L)=WCU(L)
              end if
              !       write(iu_scm_prt,888) ic,L,WCUDEEP(L,ic)
              !888    format(1x,'mstcnv--  deep  ic l wcudeep ',
              !    *                     i5,i5,f10.4)
            end do
          end if
        endif
#endif

        if(PL(LMIN).lt.850.d0) then
          LHX1=LHE
          if(TL(LMIN).lt.TF) LHX1=LHS
          U00L(LMIN)=1.d0-2.*(U00b*2.d-4/QSAT(TL(LMIN),LHX1,PL(LMIN)))
        else
          do L=1,LMIN
            LHX1=LHE
            if(TL(L).lt.TF) LHX1=LHS
            U00L(L)=1.d0-2.*(U00b*2.d-4/QSAT(TL(L),LHX1,PL(L)))
          end do
        end if
        if(TPSAV(LMAX).ge.TF) LFRZ=LMAX
        U00L(LMAX)=U00a

        !**** Add plume characteristics at LMAX
        DM(LMAX)=DM(LMAX)+MPMAX
        DSM(LMAX)=DSM(LMAX)+SMPMAX
        DSMOM(xymoms,LMAX)=DSMOM(xymoms,LMAX) + SMOMPMAX(xymoms)
        DQM(LMAX)=DQM(LMAX)+QMPMAX
        DQMOM(xymoms,LMAX)=DQMOM(xymoms,LMAX) + QMOMPMAX(xymoms)
#ifdef TRACERS_ON
       DO N=1,NTX
        select case (trname(ntix(n)))
        case('NH3')
      DTM(LMIN,N) = DTM(LMIN,N) + TMPMAX(N)
      DTMOM(xymoms,LMIN,N) = DTMOM(xymoms,LMIN,N) + TMOMPMAX(xymoms,N)

        case default
      DTM(LMAX,N) = DTM(LMAX,N) + TMPMAX(N)
      DTMOM(xymoms,LMAX,N) = DTMOM(xymoms,LMAX,N) + TMOMPMAX(xymoms,N)
        end select
       end do
#endif
        CCM(LMAX)=0.
        do K=1,KMAX
          DUM(K,LMAX)=DUM(K,LMAX)+UMP(K)
          DVM(K,LMAX)=DVM(K,LMAX)+VMP(K)
        end do
        CDHM=0.
        if(MINLVL.gt.LMIN) MINLVL=LMIN
        if(MAXLVL.lt.LMAX) MAXLVL=LMAX
        if(LMCMIN.eq.0) LMCMIN=LMIN
        if(LMCMAX.lt.MAXLVL) LMCMAX=MAXLVL

        !****
        !**** DOWNDRAFT DESCENT AND TRANSPORT LOOP (4)
        !****
        LDMIN=LDRAFT-1
        LLMIN=LDMIN
        EDRAFT=0.

        if(ETADN.gt.1d-10) then ! downdraft possible

          !**** DOWNDRAFT MASS, TEMPERATURE, HUMIDITY, MOMENTUM
          DDRAFT=DDR(LDRAFT)
          DDROLD=DDRAFT
          SMDN=SMDNL(LDRAFT)
          QMDN=QMDNL(LDRAFT)
          SMOMDN(xymoms)=SMOMDNL(xymoms,LDRAFT)
          QMOMDN(xymoms)=QMOMDNL(xymoms,LDRAFT)
          do K=1,KMAX
            UMDN(K)=UMDNL(K,LDRAFT)
            VMDN(K)=VMDNL(K,LDRAFT)
          end do
#ifdef TRACERS_ON
          TMDN(:)=TMDNL(LDRAFT,:)
          TMOMDN(xymoms,:)=TMOMDNL(xymoms,LDRAFT,:)
#endif

          !**** LOOP FROM TOP DOWN OVER POSSIBLE DOWNDRAFTS
          !****
          DOWNDRAFT: do L=LDRAFT,1,-1
            LHX=VLAT(L)               ! LHX consistency
            SLH=LHX*BYSHA
            TNX1=SMDN*PLK(L)/DDRAFT   ! save for tracers

            call get_dq_evap(smdn,qmdn,plk(l),ddraft,lhx,pl(l),cond(l) &
                 ,dqsum,fqcond1)

            !**** EVAPORATE CONVECTIVE CONDENSATE IN DOWNDRAFT AND UPDATE DOWNDRAFT
            !**** TEMPERATURE AND HUMIDITY; CURRENTLY ALL CONDENSATE IS ALLOWED TO
            !**** EVAPORATE IF POSSIBLE (FDDRT = 1)
            DQEVP=FDDRT*COND(L)       ! limit evap from condensate to fddrt of amount
            if(DQEVP.gt.DQSUM) DQEVP=DQSUM           ! limit evaporation
            if(DQEVP.gt.SMDN*PLK(L)/SLH) DQEVP=SMDN*PLK(L)/SLH
            if (L.lt.LMIN) DQEVP=0.

            FSEVP = 0
            if (PLK(L)*SMDN.gt.teeny) FSEVP = SLH*DQEVP/(PLK(L)*SMDN)
            SMDN=SMDN-SLH*DQEVP/PLK(L)
            SMOMDN(xymoms)=SMOMDN(xymoms)*(1.-FSEVP)

            FQEVP = 0
            if (COND(L).gt.0.) FQEVP = DQEVP/COND(L)
            QMDN=QMDN+DQEVP

            !**** REMOVE EVAPORATED WATER FROM CONVECTIVE CONDENSATE AMOUNT
            COND(L)=COND(L)-DQEVP
            TAUMCL(L)=TAUMCL(L)-DQEVP*FMC1
            CDHEAT(L)=CDHEAT(L)-DQEVP*SLH
            EVPSUM=EVPSUM+DQEVP*SLH

#ifdef TRACERS_WATER
            !**** RE-EVAPORATION OF TRACERS IN DOWNDRAFTS
            !**** (If 100% evaporation, allow all tracers to evaporate completely.)
            if(FQEVP.eq.1.) then  ! total evaporation
              TMDN(1:NTX)     = TMDN(1:NTX) + TRCOND(1:NTX,L)
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) &
                   trdvap_mc(l,1:ntx)=trdvap_mc(l,1:ntx)+trcond(1:ntx,l)
#endif
              TRCOND(1:NTX,L) = 0.d0
            else            ! otherwise, tracers evaporate dependent on type of tracer
              call GET_EVAP_FACTOR(NTX,TNX1,LHX,.false.,1d0,FQEVP,FQEVPT,ntix)
              dtr(1:ntx) = fqevpt(1:ntx)*trcond(1:ntx,l)
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) trdvap_mc(l,1:ntx)=trdvap_mc(l,1:ntx) &
                   +dtr(1:ntx)
#endif
              TMDN(1:NTX)     = TMDN(1:NTX)     + DTR(1:NTX)
              TRCOND(1:NTX,L) = TRCOND(1:NTX,L) - DTR(1:NTX)
            end if
#endif
            !**** ENTRAINMENT INTO DOWNDRAFTS
            if(L.lt.LDRAFT.and.L.gt.1) then
              DDRUP=DDRAFT
              DDRAFT=DDRAFT+DDROLD*ETAL(L)   ! add in entrainment
              if(DDRUP.gt.DDRAFT) DDRAFT=DDRUP
              SMIX=SMDN/(DDRUP+teeny)
              QMIX=QMDN/(DDRUP+teeny)
              SVMIX=SMIX*PLK(L-1)              ! *(1.+DELTX*QMIX)
              SVM1=SM1(L-1)*BYAM(L-1)*PLK(L-1) ! *(1.+DELTX*QM1(L-1)*BYAM(L-1))

              if ((SVMIX-SVM1).ge.DTMIN1) then
                DDRAFT=FDDET*DDRUP             ! detrain downdraft if buoyant
              end if

              !**** LIMIT SIZE OF DOWNDRAFT IF NEEDED
              if(DDRAFT.gt..95d0*(AIRM(L-1)+DMR(L-1))) &
                   DDRAFT=.95d0*(AIRM(L-1)+DMR(L-1))
              EDRAFT=DDRAFT-DDRUP

              !**** ENTRAIN INTO DOWNDRAFT, UPDATE TEMPERATURE AND HUMIDITY
              if (EDRAFT.gt.0) then  ! usual case, entrainment into downdraft
                FENTRA=EDRAFT*BYAM(L)
                SENV=SM(L)*BYAM(L)
                QENV=QM(L)*BYAM(L)
                SMDN=SMDN+EDRAFT*SENV
                QMDN=QMDN+EDRAFT*QENV
                SMOMDN(xymoms)= SMOMDN(xymoms)+ SMOM(xymoms,L)*FENTRA
                QMOMDN(xymoms)= QMOMDN(xymoms)+ QMOM(xymoms,L)*FENTRA
                DSMR(L)=DSMR(L)-EDRAFT*SENV
                DSMOMR(:,L)=DSMOMR(:,L)-SMOM(:,L)*FENTRA
                DQMR(L)=DQMR(L)-EDRAFT*QENV
                DQMOMR(:,L)=DQMOMR(:,L)-QMOM(:,L)*FENTRA
                DMR(L)=DMR(L)-EDRAFT
                !**** EFFECT OF ENTRAINMENT AND CONVECTIVE PRESSURE GRADIENT ON
                !**** HORIZONTAL MOMENTUM OF DOWNDRAFT AIR
                do K=1,KMAX        ! add in momentum entrainment
                  UMTEMP=PGRAD*DDRAFT**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
                  VMTEMP=PGRAD*DDRAFT**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
                  UMDN(K)=UMDN(K)+FENTRA*UM(K,L)-UMTEMP
                  VMDN(K)=VMDN(K)+FENTRA*VM(K,L)-VMTEMP
                  DUM(K,L)=DUM(K,L)-FENTRA*UM(K,L)+UMTEMP
                  DVM(K,L)=DVM(K,L)-FENTRA*VM(K,L)+VMTEMP
                end do
#ifdef TRACERS_ON
                Tenv(1:NTX)=tm(l,1:NTX)/airm(l)
                TMDN(1:NTX)=TMDN(1:NTX)+EDRAFT*Tenv(1:NTX)
                TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)+ TMOM(xymoms,L &
                     ,1:NTX)*FENTRA
                DTMR(L,1:NTX)=DTMR(L,1:NTX)-EDRAFT*TENV(1:NTX)
                DTMOMR(:,L,1:NTX)=DTMOMR(:,L,1:NTX)-TMOM(:,L,1:NTX)*FENTRA
#endif
              else  ! occasionally detrain into environment if ddraft too big
                FENTRA=EDRAFT/(DDRUP+teeny)  ! < 0
                DSM(L)=DSM(L)-FENTRA*SMDN
                DSMOM(xymoms,L)=DSMOM(xymoms,L)-SMOMDN(xymoms)*FENTRA
                DQM(L)=DQM(L)-FENTRA*QMDN
                DQMOM(xymoms,L)=DQMOM(xymoms,L)-QMOMDN(xymoms)*FENTRA
                SMDN=SMDN*(1+FENTRA)
                QMDN=QMDN*(1+FENTRA)
                SMOMDN(xymoms)= SMOMDN(xymoms)*(1+FENTRA)
                QMOMDN(xymoms)= QMOMDN(xymoms)*(1+FENTRA)
                DM(L)=DM(L)-EDRAFT       ! EDRAFT is negative here
                do K=1,KMAX          ! add in momentum detrainment
                  UMTEMP=PGRAD*DDRAFT**2*(U_0(K,L+1)-U_0(K,L))/(PL(L)-PL(L+1))
                  VMTEMP=PGRAD*DDRAFT**2*(V_0(K,L+1)-V_0(K,L))/(PL(L)-PL(L+1))
                  UMDN(K)=UMDN(K)*(1.+FENTRA)-UMTEMP
                  VMDN(K)=VMDN(K)*(1.+FENTRA)-VMTEMP
                  DUM(K,L)=DUM(K,L)-FENTRA*UMDN(K)+UMTEMP
                  DVM(K,L)=DVM(K,L)-FENTRA*VMDN(K)+VMTEMP
                end do
#ifdef TRACERS_ON
                DTM(L,1:NTX)=DTM(L,1:NTX)-FENTRA*TMDN(1:NTX)
                DTMOM(xymoms,L,1:NTX)=DTMOM(xymoms,L,1:NTX)- &
                     TMOMDN(xymoms,1:NTX)*FENTRA
                TMDN(1:NTX)=TMDN(1:NTX)*(1.+FENTRA)
                TMOMDN(xymoms,1:NTX)= TMOMDN(xymoms,1:NTX)*(1.+FENTRA)
#endif
              end if
            end if

            LDMIN=L
            LLMIN=LDMIN        ! save LDMIN for diagnostics
            !**** ALLOW DOWNDRAFT TO DESCEND BELOW CLOUD BASE IF IT IS NEGATIVELY BUOYANT
            if (L.gt.1) then
              SMIX=SMDN/(DDRAFT+teeny)
              QMIX=QMDN/(DDRAFT+teeny)
              SVMIX=SMIX*PLK(L-1)              ! *(1.+DELTX*QMIX)
              SVM1=SM1(L-1)*BYAM(L-1)*PLK(L-1) ! *(1.+DELTX*QM1(L-1)*BYAM(L-1))
              if (L.le.LMIN.and.SVMIX.ge.SVM1) exit
              DDM(L-1)=DDRAFT
              DDROLD=DDRAFT
              DDRAFT=DDRAFT+DDR(L-1)    ! add in downdraft one layer below
              SMDN=SMDN+SMDNL(L-1)
              QMDN=QMDN+QMDNL(L-1)
              SMOMDN(xymoms)=SMOMDN(xymoms)+SMOMDNL(xymoms,L-1)
              QMOMDN(xymoms)=QMOMDN(xymoms)+QMOMDNL(xymoms,L-1)
              do K=1,KMAX
                UMDN(K)=UMDN(K)+UMDNL(K,L-1)
                VMDN(K)=VMDN(K)+VMDNL(K,L-1)
              end do
#ifdef TRACERS_ON
              do N=1,NTX
                TMDN(N)=TMDN(N)+TMDNL(L-1,N)
                TMOMDN(xymoms,N)=TMOMDN(xymoms,N)+TMOMDNL(xymoms,L-1,N)
              end do
#endif
            end if

          end do DOWNDRAFT
          !**** end of loop over downdrafts

          !**** UPDATE PROPERTIES OF LOWEST LAYER TO WHICH DOWNDRAFT PENETRATES
          DSM(LDMIN)=DSM(LDMIN)+SMDN
          DSMOM(xymoms,LDMIN)=DSMOM(xymoms,LDMIN) + SMOMDN(xymoms)
          DQM(LDMIN)=DQM(LDMIN)+QMDN

          DQMOM(xymoms,LDMIN)=DQMOM(xymoms,LDMIN) + QMOMDN(xymoms)
          TDNL(LDMIN)=SMDN*PLK(LDMIN)/(DDRAFT+teeny)
          QDNL(LDMIN)=QMDN/(DDRAFT+teeny)
#ifdef TRACERS_ON
          DTM(LDMIN,1:NTX) = DTM(LDMIN,1:NTX) + TMDN(1:NTX)
          DTMOM(xymoms,LDMIN,1:NTX) = DTMOM(xymoms,LDMIN,1:NTX) + &
               TMOMDN(xymoms,1:NTX)
          TRDNL(1:NTX,LDMIN)=TMDN(1:NTX)/(DDRAFT+teeny)
#endif
          do K=1,KMAX
            DUM(K,LDMIN)=DUM(K,LDMIN)+UMDN(K)
            DVM(K,LDMIN)=DVM(K,LDMIN)+VMDN(K)
          end do
          DM(LDMIN)=DM(LDMIN)+DDRAFT
        end if

        !****
        !**** SUBSIDENCE AND MIXING LOOP (5)
        !****
        !**** Calculate vertical mass fluxes (Note CM for subsidence is defined
        !**** in opposite sense than normal (positive is down))
        if(LDMIN.gt.LMIN) LDMIN=LMIN    ! some loops require LMIN to LMAX
        do L=0,LDMIN-1
          CM(L) = 0.
        end do
        do L=LDMIN,LMAX
          CM(L) = CM(L-1) - DM(L) - DMR(L)
          SMT(L)=SM(L)    ! Save profiles for diagnostics
          QMT(L)=QM(L)
        end do
        do L=LMAX,LM
          CM(L) = 0.
        end do
        !**** simple upwind scheme for momentum
        do K=1,KMAX
          SUMU(K)=sum(UM(K,LDMIN:LMAX))
          SUMV(K)=sum(VM(K,LDMIN:LMAX))
        end do
        SUMDP=sum(AIRM(LDMIN:LMAX))
        ALPHA=0.
        do L=LDMIN,LMAX
          CLDM=CCM(L)
          if(L.lt.LDRAFT.and.L.ge.LLMIN.and.ETADN.gt.1d-10) &
               CLDM=CCM(L)-DDM(L)
          if(MC1) VSUBL(L)=100.*CLDM*RGAS*TL(L)/(PL(L)*GRAV*DTsrc)
          BETA=CLDM*BYAM(L+1)
          if(CLDM.lt.0.) BETA=CLDM*BYAM(L)
          BETAU=BETA
          ALPHAU=ALPHA
          if(BETA.lt.0.) BETAU=0.
          if(ALPHA.lt.0.) ALPHAU=0.
          do K=1,KMAX
            UM(K,L)= &
                 UM(K,L)+RA(K)*(-ALPHAU*UM(K,L)+BETAU*UM(K,L+1)+DUM(K,L))
            VM(K,L)= &
                 VM(K,L)+RA(K)*(-ALPHAU*VM(K,L)+BETAU*VM(K,L+1)+DVM(K,L))
          end do
          ALPHA=BETA
        end do
        do K=1,KMAX
          SUMU1(K)=sum(UM(K,LDMIN:LMAX))
          SUMV1(K)=sum(VM(K,LDMIN:LMAX))
        end do
        do K=1,KMAX                          ! momentum adjustment
          UM(K,LDMIN:LMAX)=UM(K,LDMIN:LMAX)-(SUMU1(K)-SUMU(K))* &
               AIRM(LDMIN:LMAX)/SUMDP
          VM(K,LDMIN:LMAX)=VM(K,LDMIN:LMAX)-(SUMV1(K)-SUMV(K))* &
               AIRM(LDMIN:LMAX)/SUMDP
        end do

        !****

        ! Determine the number of subsidence sub-timesteps such that
        ! courant numbers in the QUS do not exceed 1
        ksub = 1
        do l=ldmin,lmax-1
          if(    +cm(l) > airm(l+1)+dmr(l+1)) then
            ksub = max(ksub, 1+int((+cm(l)-dmr(l+1))/airm(l+1)) )
          elseif(-cm(l) > airm(l  )+dmr(l  )) then
            ksub = max(ksub, 1+int((-cm(l)-dmr(l  ))/airm(l  )) )
          endif
        enddo
        ksub = min(ksub,2) ! max 2 iterations allowed currently
        !      ksub = 2 ! non-interactive default

        byksub = 1d0/ksub
        nsub = lmax-ldmin+1
        cmneg=0.          ! initialization
        cmneg(ldmin:lmax-1) = -cm(ldmin:lmax-1)*byksub
        cmneg(lmax) = 0.  ! avoid roundoff error (esp. for qlimit)

        !**** Subsidence uses Quadratic Upstream Scheme for QM and SM
        do ITER=1,KSUB ! subsidence sub-timesteps to improve stability
          ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)*BYKSUB
          SM(LDMIN:LMAX) =  SM(LDMIN:LMAX)  +  DSMR(LDMIN:LMAX)*BYKSUB
          SMOM(:,LDMIN:LMAX)=SMOM(:,LDMIN:LMAX)+DSMOMR(:,LDMIN:LMAX)*BYKSUB
          call adv1d(sm(ldmin),smom(1,ldmin), f(ldmin),fmom(1,ldmin), &
               ml(ldmin),cmneg(ldmin), nsub,.false.,1, zdir,ierrt,lerrt)
          SM(LDMIN:LMAX) =   SM(LDMIN:LMAX) +   DSM(LDMIN:LMAX)*BYKSUB
          SMOM(:,LDMIN:LMAX)=SMOM(:,LDMIN:LMAX)+DSMOM(:,LDMIN:LMAX)*BYKSUB
          ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
          !****
          ML(LDMIN:LMAX) = AIRM(LDMIN:LMAX) +   DMR(LDMIN:LMAX)*BYKSUB
          QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +  DQMR(LDMIN:LMAX)*BYKSUB
          QMOM(:,LDMIN:LMAX)=QMOM(:,LDMIN:LMAX)+DQMOMR(:,LDMIN:LMAX)*BYKSUB
          call adv1d(qm(ldmin),qmom(1,ldmin), f(ldmin),fmom(1,ldmin), &
               ml(ldmin),cmneg(ldmin), nsub,.true.,1, zdir,ierrt,lerrt)
          QM(LDMIN:LMAX) =   QM(LDMIN:LMAX) +   DQM(LDMIN:LMAX)*BYKSUB
          QMOM(:,LDMIN:LMAX)=QMOM(:,LDMIN:LMAX)+DQMOM(:,LDMIN:LMAX)*BYKSUB
          ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
#ifdef TRACERS_ON
          !**** Subsidence of tracers by Quadratic Upstream Scheme
          do N=1,NTX
            ML(LDMIN:LMAX) =  AIRM(LDMIN:LMAX) +    DMR(LDMIN:LMAX)*BYKSUB
            TM(LDMIN:LMAX,N) =  TM(LDMIN:LMAX,N) + DTMR(LDMIN:LMAX,N)*BYKSUB
            TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N)+DTMOMR(:,LDMIN:LMAX,N) &
                 *BYKSUB
            call adv1d(tm(ldmin,n),tmom(1,ldmin,n), f(ldmin),fmom(1,ldmin), &
                 ml(ldmin),cmneg(ldmin), nsub,t_qlimit(n),1, zdir,ierrt,lerrt)
            TM(LDMIN:LMAX,N) = TM(LDMIN:LMAX,N) +   DTM(LDMIN:LMAX,N)*BYKSUB
            TMOM(:,LDMIN:LMAX,N) = TMOM(:,LDMIN:LMAX,N) +DTMOM(:,LDMIN:LMAX,N) &
                 *BYKSUB
            ierr=max(ierrt,ierr) ; lerr=max(lerrt+ldmin-1,lerr)
          end do
#endif
        end do        ! end of sub-timesteps for subsidence
        !**** Check for v. rare negative humidity error condition
        do L=LDMIN,LMAX
          if(QM(L).lt.0.d0) then
            write(6,*) ' Q neg: it,i,j,l,q,cm',itime,i_debug,j_debug,l &
                 ,qm(l),cmneg(l)
            !**** reduce subsidence post hoc.
            LM1=max(1,L-1)
            if (QM(LM1)+QM(L).lt.0) then
              write(6,*) "Q neg cannot be fixed!",L,QM(LM1:L)
            else
              QM(L-1)=QM(L-1)+QM(L)
              QM(L)=0.
#ifdef TRACERS_WATER
              !**** corresponding water tracer adjustment
              do N=1,NTX
                if (tr_wd_type(n) .eq. nWater) then
                  TM(L-1,N)=TM(L-1,N)+TM(L,N)
                  TM(L,N)=0.
                end if
              end do
#endif
            end if
          end if
        end do
#ifdef TRACERS_ON
        !**** check for independent tracer errors
        do N=1,NTX
          if (.not.t_qlimit(n)) cycle
#ifdef TRACERS_WATER
          if (tr_wd_type(n) .eq. nWater) cycle ! water tracers already done
#endif
          do L=LDMIN,LMAX
!            if (TM(L,N).lt.0.) then
!              write(6,*) trname(n),' neg: it,i,j,l,tr,cm',itime,i_debug &
!                   ,j_debug,l,tm(l,n),cmneg(l)
!              !**** reduce subsidence post hoc.
!              LM1=max(1,L-1)
!              if (TM(LM1,N)+TM(L,N).lt.0) then
!                write(6,*) trname(n)," neg cannot be fixed!",L,TM(LM1:L,N)
!              else
!                TM(L-1,N)=TM(L-1,N)+TM(L,N)
!                TM(L,N)=0.
!              end if
!            end if
            if (tm(l,n).lt.0.) then ! borrow mass from below
              vsum = tm(l,n)
              lborrow1 = l
              do while(vsum.lt.0. .and. lborrow1.gt.1)
                lborrow1 = lborrow1 - 1
                vsum = vsum + tm(lborrow1,n)
              enddo
              if(vsum.lt.0.) then
                write(6,*) trname(n)," neg cannot be fixed!",L,TM(1:L,N)
              else
                if(l-lborrow1.gt.1) then
                  write(6,*) trname(n),' nonlocal borrow: it,i,j,l,tr,cm', &
                       itime,i_debug,j_debug,l,tm(lborrow1:l,n),cmneg(l)
                else
                  write(6,*) trname(n),' neg: it,i,j,l,tr,cm', &
                       itime,i_debug,j_debug,l,tm(l,n),cmneg(l)
                endif
                ! note: borrowing from more than one layer is done by
                ! multiplication rather than subtraction
                tm(lborrow1:l-1,n)=tm(lborrow1:l-1,n)*(vsum/(vsum-tm(l,n)))
                tm(l,n)=0.
              endif
            endif
          end do
        end do
#endif
        !**** diagnostics
        do L=LDMIN,LMAX
          FCDH=0.
          if(L.eq.LMAX) FCDH=CDHSUM-CDHSUM1+CDHM
          FCDH1=0.
          if(L.eq.LLMIN) FCDH1=CDHSUM1-EVPSUM
          MCFLX(L)=MCFLX(L)+CCM(L)*FMC1
          DGDSM(L)=DGDSM(L)+(PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
          DGDQM(L)=DGDQM(L)+SLHE*(QM(L)-QMT(L))*FMC1
          DQMTOTAL(L)=DQMTOTAL(L)+(QM(L)-QMT(L))*BYAM(L)*FMC1
          if(PLE(LMAX+1).gt.700.d0) then
            DGSHLW(L)=DGSHLW(L)+ &
                 (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
            DQMSHLW(L)=DQMSHLW(L)+(QM(L)-QMT(L))*BYAM(L)*FMC1
          endif
          if(PLE(LMIN)-PLE(LMAX+1).ge.450.d0) then
            DGDEEP(L)=DGDEEP(L)+ &
                 (PLK(L)*(SM(L)-SMT(L))-FCDH-FCDH1)*FMC1
            DQMDEEP(L)=DQMDEEP(L)+(QM(L)-QMT(L))*BYAM(L)*FMC1
          endif
          DTOTW(L)=DTOTW(L)+SLHE*(QM(L)-QMT(L)+COND(L))*FMC1
          DDMFLX(L)=DDMFLX(L)+DDM(L)*FMC1
#ifdef SCM
          if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
            CUMFLX(L) = 100.*MCFLX(L)*bygrav/dtsrc
            DWNFLX(L) = 100.*DDMFLX(L)*bygrav/dtsrc
            !           write(iu_scm_prt,*) 'L CUMFLX DWNFLX ',
            !    &            L,CUMFLX(L),DWNFLX(L)
          endif
#endif
        end do
        !**** save new 'environment' profile for static stability calc.
        do L=1,LM
          SM1(L)=SM(L)
          QM1(L)=QM(L)
        end do
#ifdef TRACERS_ON
        TM1(:,1:NTX) = TM(:,1:NTX)
#endif

        !****
        !**** Partition condensate into precipitation and cloud water
        !****
        COND(LMAX)=COND(LMAX)+CONDV(LMAX)
#ifdef TRACERS_WATER
        TRCOND(:,LMAX)=TRCOND(:,LMAX)+TRCONDV(:,LMAX)
#endif
        if(PLE(LMIN)-PLE(LMAX+1).ge.450.) then
          do L=LMAX,LMIN,-1
            if(COND(L).lt.CONDP(L)) CONDP(L)=COND(L)
            FCLW=0.
            if (COND(L).gt.0) FCLW=(COND(L)-CONDP(L))/COND(L)
#ifdef SCM
            if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
              !ccc  in g/m3
              !            PRCCDEEP(L,IC,LMIN) = FMC1*1.d5*CONDP(L) *         ! save COND
              !    *             BYAM(L)*PL(L)/(RGAS*TL(L))
              !            NPRCCDEEP(L,IC,LMIN) = FMC1*1.d5*(COND(l)-CONDP(L)) *
              !    *             BYAM(L)*PL(L)/(RGAS*TL(L))
              !ccc  in kg/kg
              PRCCDEEP(L,IC,LMIN) = FMC1*CONDP(L)*BYAM(L)       ! save COND
              NPRCCDEEP(L,IC,LMIN) = FMC1*(COND(l)-CONDP(L))*BYAM(L)
              PRCCGRP(L,IC,LMIN) = FMC1*CONDGP(L)*BYAM(L)
              PRCCICE(L,IC,LMIN) = FMC1*CONDIP(L)*BYAM(L)
              DETRAINDEEP(L,IC,LMIN) = FCLW*COND(L)*BYAM(L)*FMC1
              !            write(iu_scm_prt,889) lmin,ic,FMC1,L,
              !    *              PRCCDEEP(L,IC,LMIN)*1000.,
              !    *              NPRCCDEEP(L,IC,LMIN)*1000.,TPALL(L,IC,LMIN),
              !    *              PRCCGRP(L,IC,LMIN)*1000.,
              !    *              PRCCICE(L,IC,LMIN)*1000.
              !889         format(1x,
              !    *        'mc--dp lmin ic fmc1 l prcc nprcc tp grp ice',
              !    *              2(i4),f8.3,i4,f9.5,f9.5,f8.2,f9.5,f9.5)
            endif
#endif
            !**** check in case saved condensate is different phase
            if (SVLATL(L).gt.0 .and. SVLATL(L).ne.VLAT(L)) HEAT1(L)= &
                 HEAT1(L)+(SVLATL(L)-VLAT(L))*svwmxl(l)*airm(l)*BYSHA/FMC1
            SVLATL(L)=VLAT(L)         ! moved from above
            SVWMXL(L)=SVWMXL(L)+FCLW*COND(L)*BYAM(L)*FMC1
            COND(L)=CONDP(L)

#ifdef TRACERS_WATER
            !**** Apportion cloud tracers and condensation
            !**** Note that TRSVWML is in mass units unlike SVWMX
            TRSVWML(1:NTX,L) = TRSVWML(1:NTX,L) + FCLW*TRCOND(1:NTX,L) &
                 *FMC1
#ifdef TRDIAG_WETDEPO
            if (diag_wetdep == 1) &
                 trflcw_mc(l,1:ntx)=trflcw_mc(l,1:ntx)+fclw*trcond(1:ntx,l)
#endif
            TRCOND(1:NTX,L) = (1.-FCLW)*TRCOND(1:NTX,L)
#endif
          end do
        end if

        !****
        !**** RE-EVAPORATION AND PRECIPITATION LOOP (6); INCLUDES CALCULATION OF
        !**** CONVECTIVE CLOUD FRACTION
        !****

#ifdef SCM
        if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          do L=1,LM
            MCCOND(L,IC,LMIN) = COND(L)*FMC1*BYAM(L)
          enddo
        endif
#endif
        PRCP=COND(LMAX)
        PRHEAT=CDHEAT(LMAX)

        !**** check whether environment is the same phase as cond
        TOLD=SMOLD(LMAX)*PLK(LMAX)*BYAM(LMAX)
        lhp(lmax)=lhe
        if (TOLD.le.TF) lhp(lmax)=lhs

        !**** adjust phase of precip to that of environment
        if ((TOLD.gt.TF .and. vlat(Lmax).eq.lhs) .or. (TOLD.le.TF .and. &
             vlat(lmax).eq.lhe)) then
          FSSUM = 0
          if (abs(PLK(LMAX)*SM(LMAX)).gt.teeny .and. ((lhp(lmax)-vlat(lmax &
               ))*PRCP*BYSHA).lt.0) FSSUM = -(lhp(lmax)-vlat(lmax))*PRCP &
               *BYSHA/(PLK(LMAX)*SM(LMAX))
          if (debug) print*,"cnv0",lmax,(lhp(lmax)-vlat(lmax))*PRCP*BYSHA &
               /PLK(LMAX),lhp(lmax),vlat(lmax)
          SM(LMAX)=SM(LMAX)+(lhp(lmax)-vlat(lmax))*PRCP*BYSHA/PLK(LMAX)
          SMOM(:,LMAX) =  SMOM(:,LMAX)*(1.-FSSUM)
        end if

        !**** add in heat1 from lmax
        if (heat1(lmax).ne.0) print*,"cnvA",i_debug,j_debug,lmax &
             ,heat1(lmax),vlat(lmax),lhp(lmax),prcp

        !      FSSUM = 0
        !      IF (ABS(PLK(LMAX)*SM(LMAX)).gt.teeny .and. HEAT1(lmax).gt.0) FSSUM
        !     *     = -heat1(lmax)/(PLK(LMAX)*SM(LMAX))
        !      SM(LMAX)=SM(LMAX)-heat1(lmax)/PLK(LMAX)
        !      SMOM(:,LMAX) =  SMOM(:,LMAX)*(1.-FSSUM)


#ifdef TRACERS_WATER
        !**** Tracer precipitation
        ! Note that all of the tracers that condensed do not precipitate here,
        ! since a fraction (FCLW) of TRCOND was removed above.
        TRPRCP(1:NTX) = TRCOND(1:NTX,LMAX)
#ifdef TRDIAG_WETDEPO
        if (diag_wetdep == 1) &
             trprcp_mc(lmax,1:ntx)=trprcp_mc(lmax,1:ntx) &
             +trcond(1:ntx,lmax)
#endif
#endif
        DPHASE(LMAX)=DPHASE(LMAX)+(CDHSUM-CDHSUM1+ &
             CDHM)*FMC1
        if(PLE(LMAX+1).gt.700.d0) DPHASHLW(LMAX)=DPHASHLW(LMAX)+ &
             (CDHSUM-CDHSUM1+CDHM)*FMC1
        if(PLE(LMIN)-PLE(LMAX+1).ge.450.d0) DPHADEEP(LMAX)= &
             DPHADEEP(LMAX)+(CDHSUM-CDHSUM1+CDHM)*FMC1

        !**** Loop down from top of plume
        EVAP_PRECIP: do L=LMAX-1,1,-1
#ifdef SCM
          !     if (SCM_ATURB_FLAG.eq.0) then  ! commented out by yao
          !         for SCM running with Dry Convection instead of ATURB WTURB
          !         not filled -- set TRATIO=1
          !         TRATIO = 1.d0
          !     else
          !         for SCM running with ATURB
          !         TTURB=HPBL/WTURB(L)
          !         TRATIO=TTURB/DTsrc
          !     endif
#else
          !     TTURB=HPBL/WTURB(L)
          !     TRATIO=TTURB/DTsrc
#endif

          !**** CALCULATE CONVECTIVE CLOUD FRACTION
          call MC_CLOUD_FRACTION(L,TL(L),PL(L),CCM(L),WCU(L),PLE(LMIN), &
               PLE(L+2),PLE(LMAX+1),CCM(LMIN),WCU(LMIN),LMIN,LMAX, &
               CCMUL,CCMUL2,FCLOUD)

          !**** PRECIPITATION IS ALLOWED TO EVAPORATE FULLY INTO A FRACTION OF
          !**** THE GRIDBOX HALF AS LARGE AS THE FRACTION OF GRIDBOX MASS THAT
          !**** CONVECTS
          FEVAP=.5*CCM(L)*BYAM(L+1)
          if(L.lt.LMIN) FEVAP=.5*CCM(LMIN)*BYAM(LMIN+1)
          if(FEVAP.gt..5) FEVAP=.5
          CLDMCL(L+1)=min(CLDMCL(L+1)+FCLOUD*FMC1,FMC1)
          CLDREF=CLDMCL(L+1)
          if(PLE(LMAX+1).gt.700..and.CLDREF.gt.CLDSLWIJ) &
               CLDSLWIJ=CLDREF
          if(PLE(LMIN)-PLE(LMAX+1).ge.450..and.CLDREF.gt.CLDDEPIJ) &
               CLDDEPIJ=CLDREF
          TOLD=SMOLD(L)*PLK(L)*BYAM(L)
          TOLD1=SMOLD(L+1)*PLK(L+1)*BYAM(L+1)

          !**** FORWARD STEP COMPUTES HUMIDITY CHANGE BY RECONDENSATION
          !**** Q = Q + F(TOLD,PRHEAT,QOLD+EVAP)
          PRECNVL(L+1)=PRECNVL(L+1)+PRCP*BYGRAV

          !**** CALCULATE CONVECTIVE PRECIP PHASE
          call MC_PRECIP_PHASE(L,AIRM(L),FEVAP, &
               LHP(L+1),PRCP,TOLD,TOLD1,VLAT(L),COND(L),LMIN, &
               LHP(L),MCLOUD,HEAT1(L))

          !**** set phase of precip based on local environment temperature
          LHX=LHP(L)

          DQSUM=0.
          FPRCP=0.
          SLH=LHX*BYSHA

          if (PRCP.gt.0.) then

            if (mcloud.gt.0) call get_dq_evap(smold(l),qmold(l),plk(l),airm(l) &
                 ,lhx,pl(l),prcp*AIRM(L)/MCLOUD,dqsum,fprcp)
            dqsum=dqsum*MCLOUD*BYAM(L)

            PRCP=PRCP-DQSUM
            QM(L)=QM(L)+DQSUM

          end if

          !**** UPDATE TEMPERATURE DUE TO NET REEVAPORATION IN CLOUDS
          FSSUM = 0
          if (abs(PLK(L)*SM(L)).gt.teeny .and. (SLH*DQSUM+HEAT1(L)).gt.0) &
               FSSUM = (SLH*DQSUM+HEAT1(L))/(PLK(L)*SM(L))
          if (debug) print*,"cnv4",l,SLH*DQSUM,HEAT1(L)
          SM(L)=SM(L)-(SLH*DQSUM+HEAT1(L))/PLK(L)
          SMOM(:,L) =  SMOM(:,L)*(1.-FSSUM)

          FCDH1=0.
          if(L.eq.LLMIN) FCDH1=CDHSUM1-EVPSUM
          DPHASE(L)=DPHASE(L)-(SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
          DQCOND(L)=DQCOND(L)-SLH*DQSUM*FMC1
          DQCTOTAL(L)=DQCTOTAL(L)-DQSUM*BYAM(L)*FMC1
          if(PLE(LMAX+1).gt.700.d0) then
            DPHASHLW(L)=DPHASHLW(L)- &
                 (SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
            DQCSHLW(L)=DQCSHLW(L)-DQSUM*BYAM(L)*FMC1
          endif
          if(PLE(LMIN)-PLE(LMAX+1).ge.450.d0) then
            DPHADEEP(L)=DPHADEEP(L)- &
                 (SLH*DQSUM-FCDH1+HEAT1(L))*FMC1
            DQCDEEP(L)=DQCDEEP(L)-DQSUM*BYAM(L)*FMC1
          endif

#ifdef TRACERS_WATER
          if (PRCP+DQSUM.gt.0.) then
            !**** Tracer net re-evaporation
            !**** (If 100% evaporation, allow all tracers to evaporate completely.)
            BELOW_CLOUD = L.lt.LMIN
            if(FPRCP.eq.1.) then      !total evaporation
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) trnvap_mc(l,1:ntx) = trnvap_mc(l,1:ntx) &
                   +trprcp(1:ntx)
#endif
              do N=1,NTX
                TM(L,N)   = TM(L,N)  + TRPRCP(N)
                !         if (debug .and.n.eq.1) print*,"cld2",L,TM(L,N),TRPRCP(N),2
                !     *         *FEVAP
                TRPRCP(N) = 0.d0
              end do
            else ! otherwise, tracers evaporate dependent on type of tracer
              !**** estimate effective humidity
              if (below_cloud) then
                TNX1=(SM(L)*PLK(L)-SLH*DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L)
                HEFF=min(1d0,(QM(L)+DQSUM*(1./(2.*MCLOUD)-1.))*BYAM(L) &
                     /QSAT(TNX1,LHX,PL(L)))
              else
                heff=1.
              end if
              call GET_EVAP_FACTOR(NTX,TOLD,LHX,BELOW_CLOUD,HEFF,FPRCP,FPRCPT,ntix)
              dtr(1:ntx) = fprcpt(1:ntx)*trprcp(1:ntx)
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) trnvap_mc(l,1:ntx)=trnvap_mc(l,1:ntx) &
                   +dtr(1:ntx)
#endif
              TM(L,1:NTX) = TM(L,1:NTX)     + DTR(1:NTX)
              !          if (debug .and.n.eq.1) print*,"cld3",L,TM(L,N),FPRCP
              !     *         ,FPRCPT(N),TRPRCP(N)
              TRPRCP(1:NTX) = TRPRCP(1:NTX) - DTR(1:NTX)
            end if
#ifndef NO_WASHOUT_IN_CLOUDS
            if (.not. below_cloud .and. prcp > teeny) then
              !**** Washout of tracers in cloud
              wmxtr=prcp*byam(l)
              precip_mm=prcp*100.*bygrav
              TM_dum(:) = TM(L,:)
              call GET_WASH_factor( &
                   ntx,fplume,precip_mm,fwasht,told,lhx, &
                   wmxtr,fplume,tm_dum,trprcp,thwash,pl(l),ntix,.true. &
#ifdef TRACERS_TOMAS
      ,i_debug,j_debug,L &
#endif
           )
              dtr(1:ntx) = fwasht(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) trwash_mc(l,1:ntx)=trwash_mc(l,1:ntx) &
                   +dtr(1:ntx)+thwash(1:ntx)
#endif
              do igas=1,gases_count
                n = gases_list(igas)
                if (tm(l,n) > teeny) then
                  tmfac(n)=thwash(n)/tm(l,n)
                else
                  tmfac(n)=0.
                end if
              enddo
              do n=1,ntx
                trprcp(n)=dtr(n)+trprcp(n)+thwash(n)
                tm(l,n)=tm(l,n)*(1.-fwasht(n))-thwash(n)
                tmom(xymoms,l,n)=tmom(xymoms,l,n)*(1.-fwasht(n)-tmfac(n))
              end do
            end if
#endif
            if (below_cloud .and. prcp > teeny) then
              !**** WASHOUT of TRACERS BELOW CLOUD
              WMXTR = PRCP*BYAM(L)
              precip_mm = PRCP*100.*bygrav
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
              WA_VOL= precip_mm*DXYPIJ

              call GET_SULFATE(L,TOLD,FPLUME,WA_VOL,WMXTR,SULFIN, &
                   SULFINC,SULFOUT,TR_LEFT,TM,TRPRCP,AIRM,LHX, &
                   DT_SULF_MC(1,L),CLDSAVT)

              do iaqch=1,aqchem_count
                n = aqchem_list(iaqch)
                TRPRCP(N)=TRPRCP(N)*(1.+SULFINC(N))
                TM(L,N)=TM(L,N)*(1.+SULFIN(N))
                TMOM(xymoms,L,N)=TMOM(xymoms,L,N) *(1.+SULFIN(N))
                TRCOND(N,L) = TRCOND(N,L)+SULFOUT(N)
              enddo
#endif

              !dmk Here I took out GET_COND, since we are below cloud.
              !dmk GET_WASH now has gas dissolution, extra arguments
              TM_dum(:) = TM(L,:)
              call GET_WASH_FACTOR(NTX,fplume,precip_mm,FWASHT, &
                   TOLD,LHX,WMXTR,FPLUME,TM_dum,TRPRCP,THWASH,pl(l),ntix,.true. &
#ifdef TRACERS_TOMAS
       ,i_debug,j_debug,L &
#endif
            )
              dtr(1:ntx) = fwasht(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
              if (diag_wetdep == 1) trwash_mc(l,1:ntx)=trwash_mc(l,1:ntx) &
                   +dtr(1:ntx)+thwash(1:ntx)
#endif
              do igas=1,gases_count
                n = gases_list(igas)
                if (TM(L,N).gt.teeny) then
                  TMFAC(N)=THWASH(N)/TM(L,N)
                else
                  TMFAC(N)=0.
                end if
              enddo
              do N=1,NTX
                TRPRCP(N) = DTR(N)+TRPRCP(N)+THWASH(N)
                TM(L,N)=TM(L,N)*(1.-FWASHT(N))-THWASH(N)
                !          if (debug .and.n.eq.1) print*,"cld4",L,TM(L,N),FWASHT(N)
                !     *         ,THWASH(N)
                TMOM(xymoms,L,N)=TMOM(xymoms,L,N) * &
                     (1.-FWASHT(N)-TMFAC(N))
              end do
            end if

          end if
#endif

          !**** ADD PRECIPITATION AND LATENT HEAT BELOW
          PRHEAT=CDHEAT(L)+SLH*PRCP
          if (debug) print*,"cnv5",l,prcp,cond(l)
          PRCP=PRCP+COND(L)
#ifdef TRACERS_WATER
          TRPRCP(1:NTX) = TRPRCP(1:NTX) + TRCOND(1:NTX,L)
#ifdef TRDIAG_WETDEPO
          if (diag_wetdep == 1) &
               trprcp_mc(l,1:ntx)=trprcp_mc(l,1:ntx)+trcond(1:ntx,l)
#endif
#ifdef TRACERS_SPECIAL_O18
          !**** Isotopic equilibration of liquid precip with water vapour
          if (LHX.eq.LHE .and. PRCP.gt.0) then
            do N=1,NTX
              call ISOEQUIL(NTIX(N),TOLD,.true.,QM(L),PRCP,TM(L,N),TRPRCP(N) &
                   ,0.5d0)
            end do
          end if
#endif
#endif
          !**** end of loop down from top of plume
        end do EVAP_PRECIP
        !****
        if(PRCP.gt.0.) then
          if(PLE(LMIN)-PLE(LMAX+1).lt.450.) then
            CLDMCL(1)=min(CLDMCL(1),FMC1)
          else
            RHO=PL(1)/(RGAS*TL(1))
            CLDMCL(1)=min(CLDMCL(1)+FMC1*CCM(LMIN)/(RHO*GRAV*WCU(LMIN)* &
                 DTsrc+teeny),FMC1)
          end if
        end if
        PRCPMC=PRCPMC+PRCP*FMC1
#ifdef TRACERS_WATER
        TRPRMC(1:NTX) = TRPRMC(1:NTX) + TRPRCP(1:NTX)*FMC1
#endif
        if(LMCMIN.gt.LDMIN) LMCMIN=LDMIN
        !****
        !**** END OF INNER LOOP (2) OVER CLOUD TYPES
        !****
        MC1=.false.

#ifdef TRACERS_ON
        call reset_tracer_work_arrays(ldmin,lmax)
#endif

      end do CLOUD_TYPES

      !****
      !**** END OF OUTER LOOP (1) OVER CLOUD BASE
      !****
    end do CLOUD_BASE

    if(LMCMIN.gt.0) then

      !**** set fssl array
      do L=1,LM
        FSSL(L)=1.-FMC1
      end do
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
      if (diag_wetdep == 1) then
        do n=1,ntx
          do l=1,lmcmax
            trcond_mc(l,n)=trcond_mc(l,n)*fmc1
            trdvap_mc(l,n)=trdvap_mc(l,n)*fmc1
            trflcw_mc(l,n)=trflcw_mc(l,n)*fmc1
            trprcp_mc(l,n)=trprcp_mc(l,n)*fmc1
            trnvap_mc(l,n)=trnvap_mc(l,n)*fmc1
            trwash_mc(l,n)=trwash_mc(l,n)*fmc1
          enddo
        enddo
      end if
#endif

      !**** ADJUSTMENT TO CONSERVE CP*T DURING SUBSIDENCE
      SUMAJ=0.
      SUMDP=0.
      do L=LMCMIN,LMCMAX
        SUMDP=SUMDP+AIRM(L)*FMC1
        SUMAJ=SUMAJ+DGDSM(L)
      end do
      do L=LMCMIN,LMCMAX
        DGDSM(L)=DGDSM(L)-SUMAJ*AIRM(L)*FMC1/SUMDP
        SM(L)=SM(L)-SUMAJ*AIRM(L)/(SUMDP*PLK(L))
        if (debug) print*,"cnv6",l,SUMAJ*AIRM(L)/SUMDP
      end do

      !**** LOAD MASS EXCHANGE ARRAY FOR GWDRAG
#ifdef CUBED_SPHERE
      ! maximum at any level
      AIRXL = maxval(MCFLX(LMCMIN:LMCMAX))
#else
      ! sum over levels
      AIRXL = 0.
      do L=LMCMIN,LMCMAX
        AIRXL = AIRXL+MCFLX(L)
      end do
#endif
    end if

    !****
    !**** CALCULATE CONVECTIVE CLOUD OPTICAL THICKNESS
    !****
    WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
    WMSUM=0.
#ifdef SCM
    if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
      SCM_LWP_MC = 0.d0
      SCM_IWP_MC = 0.d0
      SCM_WM_MC = 0.d0
    endif
#endif
#ifdef CLD_AER_CDNC
    WMCLWP=0.  ; WMCTWP=0. ; ACDNWM=0. ; ACDNIM=0.
    AREWM=0.   ; AREIM=0.  ; ALWWM=0.  ; ALWIM=0.
    NMCW=0     ; NMCI=0
#endif
    OPTICAL_THICKNESS: do L=1,LMCMAX
      TL(L)=(SM(L)*BYAM(L))*PLK(L)
      TEMWM=(TAUMCL(L)-SVWMXL(L)*AIRM(L))*1.d2*BYGRAV
      if(TL(L).ge.TF) WMSUM=WMSUM+TEMWM ! pick up water path
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
        SCM_WM_MC(L) = TAUMCL(L)*BYAM(L)-SVWMXL(L)
        if (TL(L).ge.TF) SCM_LWP_MC = SCM_LWP_MC + TEMWM
        if (TL(L).lt.TF) SCM_IWP_MC = SCM_IWP_MC + TEMWM
      endif
#endif



#ifdef CLD_AER_CDNC
      WMCTWP=WMCTWP+TEMWM
      if(TL(L).ge.TF) WMCLWP=WMCLWP+TEMWM
#endif
      !**** DEFAULT OPTICAL THICKNESS = 8 PER 100 MB CLOUD DEPTH, BUT 2 PER
      !**** 100 MB INSTEAD FOR DETRAINMENT LEVEL OF SHALLOW/MIDLEVEL
      !**** CONVECTION.  ALSO, AN OPTICAL THICKNESS OF 2 PER 100 MB IS
      !**** ASSUMED BELOW CLOUD BASE FOR RAIN FALLING FROM DEEP CONVECTION
      if(CLDMCL(L).gt.0.) then
        TAUMCL(L)=AIRM(L)*COETAU
        if(L.eq.LMCMAX .and. PLE(LMCMIN)-PLE(LMCMAX+1).lt.450.) &
             TAUMCL(L)=AIRM(L)*.02d0
        if(L.le.LMCMIN .and. PLE(LMCMIN)-PLE(LMCMAX+1).ge.450.) &
             TAUMCL(L)=AIRM(L)*.02d0
      end if
      SVLAT1(L) = SVLATL(L)   ! used in large-scale clouds
      if(SVLATL(L).eq.0.) then
        SVLATL(L)=LHE
        if ( (TPSAV(L).gt.0. .and. TPSAV(L).lt.TF) .or. &
             (TPSAV(L).eq.0. .and. TL(L).lt.TF) ) SVLATL(L)=LHS
      end if

      !**** FOR DEEP CONVECTIVE ANVILS USE DETRAINED CLOUD WATER AND
      !**** EFFECTIVE RADIUS BASED ON FIXED CLOUD DROPLET NUMBER
      !**** CONCENTRATION TO CALCULATE OPTICAL THICKNESS (ANALOGOUS TO
      !**** CALCULATION FOR STRATIFORM CLOUDS)
      if(SVWMXL(L).gt.0.) then
        FCLD=CLDMCL(L)+1.E-20
        TEM=1.d5*SVWMXL(L)*AIRM(L)*BYGRAV
        WTEM=1.d5*SVWMXL(L)*PL(L)/(FCLD*TL(L)*RGAS)
        if(SVLATL(L).eq.LHE.and.SVWMXL(L)/FCLD.ge.WCONST*1.d-3) &
             WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
        if(WTEM.lt.1.d-10) WTEM=1.d-10
        !**   Set CDNC for moist conv. clds (const at present)
        MNdO = 59.68d0/(RWCLDOX**3)
        MNdL = 174.d0
        MNdI = 0.06417127d0
#ifdef CLD_AER_CDNC
#ifdef ALT_CDNC_INPUTS
        ! (min+max)/2 of values from all updrafts at this level.
        ! could average over all updrafts with appropriate weighting instead
        MNdO=.5*(MNdO_max(L)+MNdO_min(L))
        MNdL=.5*(MNdL_max(L)+MNdL_min(L))
#else
        ! values from the last updraft computation at its detrainment level
        MNdO=MCDNO1
        MNdL=MCDNL1
#endif
#endif
        MCDNCW=MNdO*(1.-PEARTH)+MNdL*PEARTH
        !          if(MCDNCW.gt.0.) write(6,*)"CDNC MC cld",MNdO,MNdL,l
        MCDNCI=MNdI

        !**** COMPUTE ANVIL OPTICAL THICKNESS AND DROPLET RADIUS
        call ANVIL_OPTICAL_THICKNESS(SVLATL(L),RCLDX,MCDNCW,MCDNCI, &
             RIMAX,BYBR,FCLD,TEM,WTEM,RCLD,TAUMCL(L))

        RCLDE=RCLD/BYBR       !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
        !** Using the Liu and Daum paramet, Nature, 2002, Oct 10, Vol 419
        !** for spectral dispersion effects on droplet size distribution
        !** central value of 0.003 for alfa:Rotstayn&Daum, J.Clim,2003,16,21, Nov 2003.
        Repsi=1.d0 - 0.7d0*exp(-0.003d0*MCDNCW)
        Repsis=Repsi*Repsi
        Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**by3)
        RCLDE=RCLD*Rbeta
        TAUMCL(L)=1.5*TEM/(FCLD*RCLDE+1.E-20) ! evaluate TAUMCL again
        if(TAUMCL(L).gt.100.) TAUMCL(L)=100.
        !     write(6,*)"RCLD",RCLDE,RCLD,Rbeta,WTEM,L,MCDNCW
#endif
        CSIZEL(L)=RCLDE           !  effective droplet radius in anvil
#ifdef CLD_AER_CDNC
        if (FCLD.gt.1.d-5.and.SVLATL(L).eq.LHE) then
          ACDNWM(L)= MCDNCW
          AREWM(L) = RCLDE
          ALWWM(L)= WTEM      ! cld water density in g m-3
          NMCW  = NMCW+1
          !           if(MCDNCW.gt.0.) write(6,*)"MST CDNC b",ACDNWM(L),MCDNCW,L
        elseif(FCLD.gt.1.d-5.and.SVLATL(L).eq.LHS) then
          ACDNIM(L)= MCDNCI
          AREIM(L) = RCLDE
          ALWIM(L) = WTEM
          NMCI  = NMCI+1
        end if
#endif
      end if
      if(TAUMCL(L).lt.0..and.CLDMCL(L).le.0.) TAUMCL(L)=0.
    end do OPTICAL_THICKNESS

    if(LMCMAX.le.1) then
      do L=1,LM
        if(PL(L).lt.850.d0) U00L(L)=0.
      end do
    end if

    return

#ifdef TRACERS_ON

  contains

    subroutine reset_tracer_work_arrays(l1,l2)
      integer :: l1,l2
      integer :: l,n
      tpold(1:l2) = 0.
      do n=1,ntx
        do l=l1,l2
          dtm(l,n) = 0.
          dtmom(:,l,n) = 0.
          dtmr(l,n) = 0.
          dtmomr(:,l,n) = 0.
          tmdnl(l,n) = 0.
          tmomdnl(:,l,n) = 0.
        enddo
      enddo
#ifdef TRACERS_WATER
      do l=1,l2
        trcond(1:ntx,l) = 0.
        trcondv(1:ntx,l) = 0.
      enddo
#endif
    end subroutine reset_tracer_work_arrays
#endif

  end subroutine MSTCNV

  ! ****************************************************************************************

  subroutine LSCOND(IERR,WMERR,LERR,i_debug,j_debug)
!@sum  LSCOND column physics of large scale condensation
!@auth M.S.Yao/A. Del Genio (modularisation by Gavin Schmidt)
!@calls CTMIX,QSAT,DQSATDT,THBAR
    implicit none

!@var IERR,WMERR,LERR error reporting
    integer, intent(OUT) :: IERR,LERR
    real*8, intent(OUT) :: WMERR
    integer, intent(IN) :: i_debug,j_debug
    logical :: debug_out
    real*8 LHX

    !**** functions
    real*8 :: QSAT, DQSATDT,ERMAX
!@var QSAT saturation humidity
!@var DQSATDT dQSAT/dT

!@param CM00 upper limit for autoconversion rate
!@param AIRM0 scaling factor for computing rain evaporation
!@param HEFOLD e-folding length for computing HDEP
!@param COEFM coefficient for ratio of cloud water amount and WCONST
!@param COEFT coefficient used in computing PRATM
!@param COESIG coefficient for equ. 23 of Del Genio et al. (1996)
!@param COEEC coefficient for computing cloud evaporation
!@param ERP exponential power for computing ER
    !**** Adjust COEFT and COEFM to change proportion of super-cooled rain
    !**** to snow. Increasing COEFT reduces temperature range of super
    !**** -cooled rain, increasing COEFM enhances probability of snow.
    real*8, parameter :: AIRM0=100.d0, GbyAIRM0=GRAV/AIRM0
    real*8 CM00,TEM1,RCLDE1
    real*8, parameter :: HEFOLD=500.,COEFM=10.,COEFT=2.5
    real*8, parameter :: COESIG=1d-3,COEEC=1000.
    integer, parameter :: ERP=2
    real*8, dimension(KMAX) :: UMO1,UMO2,UMN1,UMN2 !@var dummy variables
    real*8, dimension(KMAX) :: VMO1,VMO2,VMN1,VMN2 !@var dummy variables
!@var Miscellaneous vertical arrays
    real*8, dimension(LM) :: &
         QSATL,RHF,ATH,SQ,ER,QHEAT,WMPR, &
         CLEARA,PREP,RH00,EC,WMXM
!@var QSATL saturation water vapor mixing ratio
!@var RHF environmental relative humidity
!@var ATH change in potential temperature
!@var SQ ERMAX dummy variables
!@var ER evaporation of precip
!@var QHEAT change of latent heat
!@var CLEARA fraction of clear region
!@var PREP precip conversion rate
!@var RH00 threshold relative humidity
!@var EC cloud evaporation rate
!@var WMXM cloud water mass (mb)
    real*8, dimension(LM+1) :: PREBAR,PREICE
!@var PREBAR,PREICE precip entering layer top for total, snow

#ifdef TRACERS_WATER
!@var TRPRBAR tracer precip entering layer top for total (kg)
    real*8, dimension(NTM,LM+1) :: TRPRBAR
!@var DTER change of tracer by evaporation (kg)
!@var FWTOQ fraction of CLW that goes to water vapour
!@var FPR fraction of CLW that precipitates
!@var FER fraction of precipitate that evaporates
    real*8 TWMTMP,FWTOQ,FPR,FER
!@var DTPRT tracer-specific change of tracer by precip (kg)
!@var DTERT tracer-specific change of tracer by evaporation (kg)
!@var DTWRT tracer-specific change of tracer by washout (kg)
!@var DTQWT tracer-specific change of tracer by condensation (kg)
!@var FWTOQT tracer-specific fraction of tracer in CLW that evaporates
!@var FQTOWT tracer-specific fraction of gas tracer condensing in CLW
!@var FPRT tracer-specific fraction of tracer in CLW that precipitates
!@var FERT tracer-specific fraction of tracer in precipitate evaporating
!@var FQCONDT fraction of tracer that condenses
    real*8, dimension(NTM) :: &
         FQTOWT,FQCONDT,FERT,FWTOQT,DTERT,DTPRT,DTQWT
    real*8 :: DTWRT,FPRT,PRLIQ
!@var BELOW_CLOUD logical- is the current level below cloud?
!@var CLOUD_YET logical- in L loop, has there been any cloud so far?
    logical BELOW_CLOUD,CLOUD_YET
!@var FWASHT  fraction of tracer scavenged by below-cloud precipitation
    real*8 :: FWASHT(NTM),TM_dum(NTM), DTR(NTM)
!@var WMXTR available water mixing ratio for tracer condensation ( )?
!@var b_beta_DT precipitating gridbox fraction from lowest precipitating
!@+   layer. The name was chosen to correspond to Koch et al. p. 23,802.
!@var precip_mm precipitation (mm) from the grid box above for washout
    real*8 WMXTR, b_beta_DT, precip_mm
    ! for tracers in general, added by Koch
    real*8, dimension(ntm) :: THLAW,TR_LEF,TR_LEFT
    real*8 THWASH(NTM),TMFAC(NTM),TMFAC2(NTM),CLDSAVT
    integer :: IGAS
!@var TR_LEF limits precurser dissolution following sulfate formation
!@var THLAW Henry's Law determination of amount of tracer dissolution
!@var THWASH Henry's Law for below cloud dissolution
!@var TMFAC,TMFAC2 used to adjust tracer moments
!@var CLDSAVT is present cloud fraction, saved for tracer use
!@var cldprec cloud fraction at lowest precipitating level
    real*8 :: cldprec
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
    ! for sulfur chemistry
!@var WA_VOL Cloud water volume (L). Used by GET_SULFATE.
    real*8 WA_VOL
    real*8, dimension(NTM) ::SULFOUT,SULFIN,SULFINC
    integer :: IAQCH
#endif
#endif

    real*8 AIRMR,BETA,BMAX &
         ,CBF,CBFC0,CK,CKIJ,CK1,CK2,CKM,CKR,CM,CM0,CM1,DFX,DQ,DQSDT &
         ,DQSUM,DQUP,DRHDT,DSE,DSEC,DSEDIF,DWDT,DWDT1,DWM,ECRATE,EXPST &
         ,FCLD,FMASS,FMIX,FPLUME,FPMAX,FQTOW,FRAT,FUNI,dqsum1,fqcond1 &
         ,FUNIL,FUNIO,HCHANG,HDEP,HPHASE,OLDLAT,OLDLHX,PFR,PMI,PML &
         ,HPBL,PRATIO,QCONV,QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2,QNEW &
         ,QNEWU,QOLD,QOLDU,QSATC,QSATE,RANDNO,RCLDE,RHI,RHN,RHO,RHT1 &
         ,RHW,SEDGE,SIGK,SLH,SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1 &
         ,THT2,TLT1,TNEW,TNEWU,TOLD,TOLDU,TOLDUP,VDEF,WCONST,WMN1,WMN2 &
         ,WMNEW,WMO1,WMO2,WMT1,WMT2,WMX1,WTEM,VVEL,RCLD,FCOND &
         ,PRATM,SMN12,SMO12,QF
    real*8 SNdO,SNdL,SNdI,SCDNCW,SCDNCI
#ifdef CLD_AER_CDNC
!@auth Menon  - storing var for cloud droplet number
    integer, parameter :: SNTM=31
    real*8 Repsis,Repsi,Rbeta,CDNL1,QAUT,DSU(SNTM),QCRIT &
         ,CDNL0,NEWCDN,OLDCDN,SNd
    real*8 dynvis(LM),DSGL(LM,SNTM),DSS(SNTM),r6,r6c
    real*8 NEWCLD,SAVCLD
    real*8, dimension(lm) :: vvel_sv,CLDSAV0
    real*8, dimension(sntm,lm) :: dsu_sv
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
    real*8 DPP,TEMPR,RHODK,PPRES,PRS        ! for 3 hrly diag
    real*8 D3DL(LM),CWCON(LM)               ! for 3 hrly diag
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
    integer,parameter         :: mkx=1   ! lm
    real*8,parameter          :: mw0 = 2.094395148947515E-15
    real*8,parameter          :: mi0 = 2.094395148947515E-15
    logical, parameter        :: lSCM=.false.
    logical, parameter        :: wSCM=.false.
    real(8), parameter        :: tiny = 1.0D-30
    real*8,dimension(mkx)     :: tk0,qk0,pk0,w0,v0,r0,rablk
    real*8,dimension(mkx)     :: tk0new,qk0new
    real*8,dimension(mkx)     :: ndrop,mdrop,ncrys,mcrys
    real*8,dimension(mkx)     :: nrain,mrain,mtau,rtau,ctau,dtau
    real*8                    :: qrArray(5),row,mr0,piby6
    real*8,dimension(mkx)     :: vrain
    real*8,dimension(mkx)     :: ndrop_old,mdrop_old
    real*8,dimension(mkx)     :: ndrop_new,mdrop_new
    real*8,dimension(mkx)     :: ndrop_blk,mdrop_blk
    real*8,dimension(mkx)     :: ndrop_res,mdrop_res
    real*8,dimension(mkx)     :: ncrys_old,mcrys_old
    real*8,dimension(mkx)     :: ncrys_new,mcrys_new
    real*8,dimension(mkx)     :: ncrys_blk,mcrys_blk
    real*8,dimension(mkx)     :: ncrys_res,mcrys_res
    real*8,dimension(mkx)     :: npccn,nprc,nnuccc,nnucci
    real*8,dimension(mkx)     :: mpccn,mprc,mnuccc,mnucci
    real*8,dimension(mkx)     :: nnuccd,nnucmd,nnucmt
    real*8,dimension(mkx)     :: mnuccd,mnucmd,mnucmt
    real*8,dimension(mkx)     :: nc_tnd,qc_tnd,ni_tnd,qi_tnd
    real*8,dimension(mkx)     :: nc_tot,qc_tot,ni_tot,qi_tot

    real*8                    :: DTB2M,QAUT_B2M
    real*8                    :: NEWCDNC,OLDCDNC
    real*8,dimension(LM)      :: DCLD
    logical                   :: ldummy=.false.
    character*8               :: sname='lscond: '
    integer                   :: nm,iuo=801
#ifdef TRACERS_AMP
    real*8                    :: naero (mkx,nmodes)
    !     real*8,dimension(lm,nmodes)   :: nactc
#endif
#ifdef TRACERS_TOMAS
      REAL*8,dimension(mkx)     :: nactl
!Can
!Can Droplet parameterization quantities
!Can
      INTEGER :: NCCNMx,NCC,NSECi
      PARAMETER (NCCNMx=100, NCC=10)
      REAL*8 SULFI, BOXVL, TOTi,TOT_MI, &
           TPi(NCCNMx), MLi(NCCNMx), SLVL(NCC), CCON(NCC), &
           NACTEarth, NACTOcean, NACT, NACTBL, &
           SMAXEarth, SMAXOcean, SMAX, &
           REFFEarth, REFFOcean, REFF, REFFBL, REFFGISS, &
           CLDTAUEarth, CLDTAUOcean, CLDTAU, CLDTAUBL, CLDTAULIQ, &
           CLDTAUICE, TPARC,PPARC, &
           WPARCOcean, WPARCEarth, WPARC, &
           RHOSI,QLWC,EPSILON,AUTO(6),DIFFLWMR,DIFFEPS

      INTEGER ITYP
      LOGICAL EX
!c$$$      REAL*8 CldLiqTauNS(IM,JM,LM), CldLiqTauBL(IM,JM,LM),
!c$$$     &                 CldLiqTauGS(IM,JM,LM), CldIceTauGS(IM,JM,LM)
!c$$$      COMMON /MICROPH/ CldLiqTauNS, CldLiqTauBL, CldLiqTauGS,
!c$$$     &                 CldIceTauGS
!Can
!Can

#endif
#endif

!@var BETA,BMAX,CBFC0,CKIJ,CK1,CK2,PRATM dummy variabls
!@var SMN12,SMO12 dummy variables
!@var AIRMR
!@var CBF enhancing factor for precip conversion
!@var CK ratio of cloud top jumps in moist static energy and total water
!@var CKM CTEI threshold in MacVean and Mason theory
!@var CKR CTEI threshold in Randall theory
!@var CM conversion rate for large cloud water content
!@var CM0,CM1 limiting autoconversion rate for large cloud water content
!@var DFX iteration increment
!@var DQ condensed water vapor
!@var DQSDT derivative of saturation vapor pressure w.r.t. temperature
!@var DQSUM,DQUP dummy variables
!@var DRHDT time change of relative humidity
!@var DSE moist static energy jump at cloud top
!@var DSEC critical DSE for CTEI to operate
!@var DSFDIF DSE-DSEC
!@var DWDT,DWDT1 time change rates of cloud water
!@var ECRATE cloud droplet evaporation rate
!@var EXPST exponential term in determining the fraction in CTEI
!@var FCLD cloud fraction
!@var FMIX,FRAT fraction of mixing in CTEI
!@var FMASS mass of mixing in CTEI
!@var FPLUME fraction of mixing in CTEI
!@var FPMAX max fraction of mixing in CTEI
!@var FQTOW fraction of water vapour that goes to CLW
!@var FUNI the probablity for ice cloud to form
!@var FUNIL,FUNIO FUNI over land, ocean
!@var HCHANG,HPHASE latent heats for changing phase
!@var HDEP,HPBL layer depth (m)  (note change of unit km-->m)
!@var OLDLAT,OLDLHX previous LHX
!@var PFR PROBABLITY OF GLACIATION OF SUPER-COOLED WATER
!@var PMI icy precip entering the layer top
!@var PML layer's cloud water devided by GRAV
!@var PRATIO PMI/PML
!@var QCONV convergence of latent heat
!@var QHEATC,QLT1,QLT2,QMN1,QMN2,QMO1,QMO2 dummy variables
!@var QNEW,QNEWU updated specific humidity
!@var QOLD,QOLDU previous specific humidity
!@var QSATC saturation vapor mixing ratio
!@var QSATE saturation vapor mixing ratio w.r.t. water
!@var RANDNO random number
!@var RCLD,RCLDE cloud particle's radius, effective radius
!@var RHI relative humidity w.r.t. ice
!@var RHN dummy variable
!@var RHO air density
!@var RHT1,RHW,SIGK dummy variables
!@var SEDGE potential temperature at layer edge
!@var SMN1,SMN2,SMO1,SMO2,TEM,TEMP,TEVAP,THT1,THT2,TLT1 dummy variables
!@var TNEW,TNEWU updated tempertures
!@var TOLD,TOLDU,TOLDUP previous temperatures
!@var VDEF = VVEL - VSUB
!@var WCONST,WMN1,WMN2,WMO1,WMO2,WMT1,WMT2,WMX1 dummy variables
!@var WMNEW updated cloud water mixing ratio
!@var WTEM cloud water density (g m**-3)
!@var VVEL vertical velocity (cm/s)
!@var FCOND QF dummy variables
    integer LN,ITER !@var LN,ITER loop variables
    logical BANDF  !@var BANDF true if Bergero-Findeisen proc. occurs
    logical FORM_CLOUDS !@var FORM_CLOUDS true if clouds are formed

    integer K,L,N  !@var K,L,N loop variables

    real*8 THBAR !@var THBAR potential temperature at layer edge

    !****
    !**** LARGE-SCALE CLOUDS AND PRECIPITATION
    !**** THE LIQUID WATER CONTENT IS PREDICTED
    !****
    IERR=0
    !****
    debug_out=.false.
    PRCPSS=0.
    HCNDSS=0.
    CKIJ=1.
    RCLDX=radius_multiplier
    RTEMP=funio_denominator
    CMX=autoconv_multiplier
    WMUIX=wmui_multiplier
    !**** initialise vertical arrays
    ER=0.
    EC=0.
    PREP=0.
    PREBAR=0.
    LHP=0.
#ifdef SCM
    LHPSAV=0.
    PRESAV=0.
#endif
    QHEAT=0.
    CLDSSL=0
    TAUSSL=0
    WMPR=0.
#ifdef TRACERS_WATER
    TRPRSS = 0.
    TRPRBAR = 0.
    BELOW_CLOUD=.false.
    CLOUD_YET=.false.
    CLDSAVT=0.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
    DT_SULF_SS(1:NTM,:)=0.
#endif
    TR_LEF(:)=1. ! currently, only aqchem_list elements vary with layer
    thlaw(:) = 0.  ! nonzero only for gas tracers
    thwash(:) = 0. ! nonzero only for gas tracers
    tmfac(:) = 0.  ! nonzero only for gas tracers
    tmfac2(:) = 0. ! nonzero only for gas tracers
    fwasht(:) = 0. ! nonzero only for aerosols
    fqcondt(:) = 0.
    fqtowt(:) = 0.
#endif
#ifdef CLD_AER_CDNC
    CDN3DL=0.
    CRE3DL=0.
    SMLWP=0.
    !      AERTAU=0.
    DSGL(:,1:SNTM)=0.
#endif
#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
    CTEML=0.
    CD3DL=0.
    CL3DL=0.
    CI3DL=0.
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
    WMXICE(:)=0.
    !      print *,sname,'WMX, WMXICE = ', WMX, WMXICE
#endif
#if (defined TRACERS_WATER) && (defined TRDIAG_WETDEPO)
    if (diag_wetdep == 1) then
      !**** initialize diagnostic arrays
      trwash_ls=0.D0
      trprcp_ls=0.D0
      trclwc_ls=0.D0
      trclwe_ls=0.D0
      trcond_ls=0.D0
    end if
#endif
    do L=1,LP50
      CLEARA(L)=1.-CLDSAVL(L)
      if(WMX(L).le.0.) CLEARA(L)=1.
#ifdef CLD_AER_CDNC
      CLDSAV0(L) = 1.-CLEARA(L)
#endif
    end do
    DQUP=0.
    TOLDUP=TL(LP50)
    PREICE(LP50+1)=0.
    WCONST=WMU*(1.-PEARTH)+WMUL*PEARTH
    SSHR=0.
    DCTEI=0.
    !****
    !**** MAIN L LOOP FOR LARGE-SCALE CONDENSATION, PRECIPITATION AND CLOUDS
    !****
    do L=LP50,1,-1
      TOLD=TL(L)
      QOLD=QL(L)
      OLDLHX=SVLHXL(L)
      OLDLAT=SVLATL(L)
      !**** COMPUTE VERTICAL VELOCITY IN CM/S
      TEMP=100.*RGAS*TL(L)/(PL(L)*GRAV)
      if(L.eq.1)  then
        VVEL=-SDL(L+1)*TEMP
      else if(L.eq.LM)  then
        VVEL=-SDL(L)*TEMP
      else
        VVEL=-.5*(SDL(L)+SDL(L+1))*TEMP
      end if
      VDEF=VVEL-VSUBL(L)
#ifdef CLD_AER_CDNC
      vvel_sv(l) = vvel ! save for opt. depth calc.
#endif

      FCLD=(1.-CLEARA(L))*FSSL(L)+teeny

      !**** COMPUTE THE PROBABILITY OF ICE FORMATION, FUNI, AND
      !**** THE PROBABLITY OF GLACIATION OF SUPER-COOLED WATER, PFR
      !**** DETERMINE THE PHASE MOISTURE CONDENSES TO
      !**** DETERMINE THE POSSIBILITY OF B-F PROCESS
      BANDF=.false.
      LHX=LHE
      CBF=1. + exp(-((TL(L)-258.16d0)/10.)**2)

      if (TL(L).le.238.16) then     ! below -35C: force ice
        LHX=LHS
      elseif (TL(L).ge.TF) then ! above freezing: force water
        LHX=LHE
      else                      ! in between: compute probability
        if(TL(L).gt.269.16) then ! OC/SI/LI clouds: water above -4
          FUNIO=0.
        else
          FUNIO=1.-exp(-((TL(L)-269.16d0)/RTEMP)**4)
        end if
        if(TL(L).gt.263.16) then ! land clouds water: above -10
          FUNIL=0.
        else
          FUNIL=1.-exp(-((TL(L)-263.16d0)/RTEMP)**4)
        end if
        FUNI=FUNIO*(1.-PEARTH)+FUNIL*PEARTH
        RANDNO=RNDSSL(1,L)       !  RANDNO=RANDU(XY)
        if(RANDNO.lt.FUNI) LHX=LHS

        if (OLDLHX.eq.LHS.and.TL(L).lt.TF) LHX=LHS   ! keep old phase
        if (OLDLHX.eq.LHE.and.TL(L).gt.269.16d0) LHX=LHE ! keep old phase
        !**** special case 1) if ice previously then stay as ice (if T<Tf)
        !       IF((OLDLHX.EQ.LHS.OR.OLDLAT.EQ.LHS).AND.TL(L).LT.TF) THEN
        if(OLDLAT.eq.LHS.and.TL(L).lt.TF.and.SVLAT1(L).gt.0.) then
          if(LHX.eq.LHE) BANDF=.true.
          LHX=LHS
        end if
        if (debug) print*,"ls0",l,oldlhx,oldlat,lhx,lhp(l)

        if (L.lt.LP50) then
          !**** Decide whether precip initiates B-F process
          PML=WMX(L)*AIRM(L)*BYGRAV
          PMI=PREICE(L+1)*DTsrc
          RANDNO=RNDSSL(2,L)     !  RANDNO=RANDU(XY)
          !**** Calculate probability of ice precip seeding a water cloud
          if (LHX.eq.LHE.and.PMI.gt.0) then
            PRATIO=min(PMI/(PML+1.E-20),10d0)
            CM00=3.d-5           ! reduced by a factor of 3
            if(ROICE.gt..1d0) CM00=3.d-4
            CM0=CM00
            if(VDEF.gt.0.) CM0=CM00*10.**(-0.2*VDEF)
            CBFC0=.5*CM0*CBF*DTsrc
            PFR=(1.-exp(-(PRATIO*PRATIO)))*(1.-exp(-(CBFC0*CBFC0)))
            if(PFR.gt.RANDNO) then
              BANDF=.true.
              LHX=LHS
            end if
          end if
          !**** If liquid rain falls into an ice cloud, B-F must occur
          if (LHP(L+1).eq.LHE .and. LHX.eq.LHS .and. PML.gt.0.) &
               BANDF=.true.
        end if
      end if
      if(LHX.eq.LHS .and. (OLDLHX.eq.LHE.or.OLDLAT.eq.LHE)) BANDF=.true.

      !**** COMPUTE THE LIMITING AUTOCONVERSION RATE FOR CLOUD WATER CONTENT
      CM00=1.d-4
      if(LHX.eq.LHS.and.SVWMXL(L).le.0d0) CM00=1.d-3
      if(LHX.eq.LHE) then                 ! reduced by a factor of 3
        CM00=3.d-5
        if(ROICE.gt..1d0) CM00=3.d-4
      end if
      CM0=CM00
      if(VDEF.gt.0.) CM0=CM00*10.**(-0.2*VDEF)

      !**** COMPUTE RELATIVE HUMIDITY
      QSATL(L)=QSAT(TL(L),LHX,PL(L))
      RH1(L)=QL(L)/QSATL(L)
      if(LHX.eq.LHS.and.WMX(L).le.0d0) then          ! Karcher and Lohmann formula
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=(2.583d0-TL(L)/207.83)*(QSAT(TL(L),LHS,PL(L))/QSATE)
        if(TL(L).lt.238.16) RH1(L)=QL(L)/(QSATE*RHW)
      end if
      !**** PHASE CHANGE OF CLOUD WATER CONTENT
      HCHANG=0.
      if(OLDLHX.eq.LHE.and.LHX.eq.LHS) HCHANG= WML(L)*LHM
      if(OLDLHX.eq.LHS.and.LHX.eq.LHE) HCHANG=-WML(L)*LHM
      if(OLDLAT.eq.LHE.and.LHX.eq.LHS) HCHANG=HCHANG+SVWMXL(L)*LHM
      if(OLDLAT.eq.LHS.and.LHX.eq.LHE) HCHANG=HCHANG-SVWMXL(L)*LHM
      if (debug) print*,"ls1",l,hchang

      SVLHXL(L)=LHX
      TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
      TH(L)=TL(L)/PLK(L)
#ifdef SCM
      !     preserving T for difference from before updating with ARM data
      if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
        if (NRINIT.ne.0) then
          TH(L) = SCM_SAVE_T(L)*PLK(L)+SCM_DEL_T(L)+ &
               HCHANG/(SHA*FSSL(L)+teeny)
          TH(L) = TH(L)/PLK(L)
        endif
      endif
#endif

      ATH(L)=(TH(L)-TTOLDL(L))*BYDTsrc
      !**** COMPUTE RH IN THE CLOUD-FREE AREA, RHF
      RHI=QL(L)/QSAT(TL(L),LHS,PL(L))
      ! this formulation is used for consistency with current practice
      RH00(L)=U00a
      if(PL(L).lt.850.d0) then

        RH00(L) = RH00(L)/(RH00(L) + (1.-RH00(L))*AIRM(L)/35.)

        if(VDEF.gt..2d0.and.LMCMAX.le.1) RH00(L)= &
             RH00(L)*min(sqrt(.2d0/VDEF),.5d0) ! dependece on vertical velocity
      end if
      if(U00L(L).gt.RH00(L)) RH00(L)=U00L(L)
      !**** Option to treat boundary layer differently
      if (do_blU00.eq.1) then
        if (L.le.DCL) then      ! boundary layer clouds
          !**** calculate total pbl depth
          HPBL=0.
          do LN=1,DCL
            HPBL=HPBL+AIRM(LN)*TL(LN)*RGAS/(GRAV*PL(LN))
          end do
          !**** Scale HPBL by HRMAX to provide tuning control for PBL clouds
          HDEP = min(HPBL,HRMAX*(1.-exp(-HPBL/HEFOLD)))
          !**** Special conditions for boundary layer contained wholly in layer 1
          if (DCL.le.1) then
            if (RIS.gt.1.) HDEP=10d0
            if (RIS.le.1..and.RI1.gt.1.) HDEP=50d0
            if (RIS.le.1..and.RI1.le.1..and.RI2.gt.1.) HDEP=100d0
          end if
          !**** Estimate critical rel. hum. based on parcel lifting argument
          RH00(L)=1.-GAMD*LHE*HDEP/(RVAP*TS*TS)
          if(RH00(L).lt.0.) RH00(L)=0.
        end if
      end if
      !****
      if(RH00(L).lt.0.) RH00(L)=0.
      if(RH00(L).gt.1.) RH00(L)=1.
      RHF(L)=RH00(L)+(1.-CLEARA(L))*(1.-RH00(L))
      !**** Set precip phase to be the same as the cloud, unless precip above
      !**** is ice and temperatures after ice melt would still be below TFrez
      LHP(L)=LHX
      if (LHP(L+1).eq.LHS .and. &
           TL(L).lt.TF+DTsrc*LHM*PREICE(L+1)*GRAV*BYAM(L)*BYSHA) &
           LHP(L)=LHP(L+1)
#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
!@auth Menon  saving aerosols mass for CDNC prediction
      do N=1,SNTM
        DSS(N)=1.d-10
        DSGL(L,N)=1.d-10
      end do
      do N=1,NTX
        select case (trname(ntix(n)))
        case('SO4')
          DSGL(L,1)=tm(l,n)     !n=4
          DSS(1) = DSGL(L,1)
        case('seasalt1')
          DSGL(L,2)=tm(l,n)     !n=6
          DSS(2) = DSGL(L,2)
        case('seasalt2')
          DSGL(L,3)=tm(l,n)     !n=7
          DSS(3) = DSGL(L,3)
        case('OCIA')
          DSGL(L,4)=tm(l,n)     !n=12
          DSS(4) = DSGL(L,4)
        case('OCB')
          DSGL(L,5)=tm(l,n)     !n=13
          DSS(5) = DSGL(L,5)
        case('BCIA')
          DSGL(L,6)=tm(l,n)     !n=9
          DSS(6) = DSGL(L,6)
        case('BCB')
          DSGL(L,7)=tm(l,n)     !n=10
          DSS(7) = DSGL(L,7)
        case('OCII')
          DSGL(L,8)=tm(l,n)     !n=11
          DSS(8) = DSGL(L,8)
        case('BCII')
          DSGL(L,9)=tm(l,n)     !n=8
          DSS(9) = DSGL(L,9)
#ifdef TRACERS_DUST
        case('Clay')
          DSGL(L,10)=tm(l,n)    !n=23
          DSS(10) = DSGL(L,10)
        case('Silt1')
          DSGL(L,11)=tm(l,n)    !n=23
          DSS(11) = DSGL(L,11)
        case('Silt2')
          DSGL(L,12)=tm(l,n)    !n=23
          DSS(12) = DSGL(L,12)
        case('Silt3')
          DSGL(L,13)=tm(l,n)    !n=23
          DSS(13) = DSGL(L,13)
#endif
#ifdef TRACERS_NITRATE
        case('NO3p')
          DSGL(L,14)=tm(l,n)    !n=23
          DSS(14) = DSGL(L,14)
#endif
#ifdef TRACERS_HETCHEM
          !**** Here are dust particles coated with sulfate
        case('SO4_d1')
          DSGL(L,15)=tm(l,n)    !n=20
          DSS(15) = DSGL(L,15)
        case('SO4_d2')
          DSGL(L,16)=tm(l,n)    !n=21
          DSS(16) = DSGL(L,16)
        case('SO4_d3')
          DSGL(L,17)=tm(l,n)    !n=22
          DSS(17) = DSGL(L,17)
#endif
#ifdef TRACERS_AEROSOLS_SOA
        case('isopp1a')
          DSGL(L,18)=tm(l,n)
          DSS(18) = DSGL(L,18)
        case('isopp2a')
          DSGL(L,19)=tm(l,n)
          DSS(19) = DSGL(L,19)
#ifdef TRACERS_TERP
        case('apinp1a')
          DSGL(L,20)=tm(l,n)
          DSS(20) = DSGL(L,20)
        case('apinp2a')
          DSGL(L,21)=tm(l,n)
          DSS(21) = DSGL(L,21)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
        case('OCocean')
          DSGL(L,22)=tm(l,n)
          DSS(22) = DSGL(L,22)
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_AEROSOLS_VBS
        case('vbsAm2')
          DSGL(L,23)=tm(l,n)
          DSS(23) = DSGL(L,23)
        case('vbsAm1')
          DSGL(L,24)=tm(l,n)
          DSS(24) = DSGL(L,24)
        case('vbsAz')
          DSGL(L,25)=tm(l,n)
          DSS(25) = DSGL(L,25)
        case('vbsAp1')
          DSGL(L,26)=tm(l,n)
          DSS(26) = DSGL(L,26)
        case('vbsAp2')
          DSGL(L,27)=tm(l,n)
          DSS(27) = DSGL(L,27)
        case('vbsAp3')
          DSGL(L,28)=tm(l,n)
          DSS(28) = DSGL(L,28)
        case('vbsAp4')
          DSGL(L,29)=tm(l,n)
          DSS(29) = DSGL(L,29)
        case('vbsAp5')
          DSGL(L,30)=tm(l,n)
          DSS(30) = DSGL(L,30)
        case('vbsAp6')
          DSGL(L,31)=tm(l,n)
          DSS(31) = DSGL(L,31)
#endif  /* TRACERS_AEROSOLS_VBS */
        end select
      end do      !end of n loop for tracers
#endif   /* tracerpart and cld-aer part */
      !***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
      WMUI=WMUIX*.001         ! .0001
      WMUSI=0.1
#if (defined CLD_AER_CDNC) && (defined TRACERS_AEROSOLS_Koch)
      call GET_CDNC(L,LHX,WCONST,WMUI,AIRM(L),WMX(L),DXYPIJ, &
           FCLD,CLEARA(L),CLDSAVL(L),DSS,PL(L),TL(L), &
           OLDCDL(L),VVEL,SME(L),DSU,CDNL0,CDNL1)
      DSU_SV(:,L) = DSU(:) ! save for opt. depth calc.
      !     write(6,*)"Where is",DSU(L),l
      SNd=CDNL1
      !C** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
      !     if(SNd.gt.20.) write(6,*)"SM 11 CDNC",NEWCDN,OLDCDN,L
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
      ! Microphysical time step
      dtB2M=DTsrc
      ! Set all tendencies to zero
      ldummy=execute_bulk2m_driver('all','2zero',mkx)
      ! Set thermodynamics
      tk0=TL(L)                 ! Temperature, [K]
      qk0=QL(L)                 ! Water vapor mixing ratio, [kq/kq]
      pk0=PL(L)                 ! Pressure, [hPa]
      w0=VVEL*1.d-02            ! Large-scale velocity, [m/s]
      v0=WTURB(L) !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
      r0=0.0  !RDTNDL(L)        ! tendency due to radiation, [K/s], not needed
      !      print*,"rad tendency",w0,v0,r0,L
      piby6=4.*atan(1.0)/6.0; row=1.e+03
      qrArray(1)=842.0e+00                      ! [m^(1-b)s]
      qrArray(2)=  0.8e+00                      ! [unitless]
      qrArray(3)=500.0e-06                      ! assumed mean rain diameter [m]
      qrArray(4)=piby6*row*qrArray(3)**3        ! mean rain mass [kg]
      qrArray(5)=1000.0e-06                     ! assumed rain conc [No/m^3]
      !      write(ou,*) 'k,   prebar(k)   vrain(k)    mrain(k)    nrain(k)'
      vrain=min(qrArray(1)*qrArray(3)**qrArray(2),9.2)  ! [m/s]
      !      write(6,*)"VRAIN",vrain,l
      mrain=100.d0*prebar(L)/vrain             ! [kq/m^3]
      nrain=mrain/qrArray(4)                ! [No/m^3]
      RHO=1d5*PL(L)/(RGAS*TL(L))
      mrain=1.d3*mrain/RHO                  !kg water/kg air
      nrain=1.d3*nrain/RHO                  !number/kg air
      ldummy=execute_bulk2m_driver('all' &
           ,tk0,qk0,pk0,w0,v0,r0)
      ! Set microphysics
      if(LHX.eq.LHE)  then
        mdrop=WMX(L)            ! drop content, [kg water/kg air]
        ndrop =OLDCDL(L)*1.d6   !convert from cm-3 to m-3
        if (WMX(L).eq.0.) ndrop=0.d0
        ncrys=0.d0;mcrys=0.0d0
      else
        WMXICE(L) = WMX(L)
        mcrys=WMXICE(L)         ! crys content, [kg water/kg air]
        ncrys=OLDCDI(L)*1.d6     ! convert cm-3 to m-3; set at 0.1 l-1 = 1.d-4 cm-3
        if (WMX(L).eq.0.) ncrys=0.d0
        ndrop=0.0d0;mdrop=0.0d0
      endif
      !
      ndrop_old=ndrop;mdrop_old=mdrop;ncrys_old=ncrys;mcrys_old=mcrys
      ndrop_new=0.0d0;mdrop_new=0.0d0;ncrys_new=0.0d0;mcrys_new=0.0d0
      nc_tnd=0.0d0;qc_tnd=0.0d0;ni_tnd=0.0d0;qi_tnd=0.0d0
      nc_tot=0.0d0;qc_tot=0.0d0;ni_tot=0.0d0;qi_tot=0.0d0
      !
      !** Convert from l-1 to cm-3====>  1 l^-1 = 10^-3 cm^-3 = 10^3 m^-3
      !      ldummy=execute_bulk2m_driver('all'
      !    *           ,ndrop,mdrop,ncrys,mcrys,'end')
#ifdef TRACERS_AMP
      do nm=1,nmodes
        naero(mkx,nm)=nactc(l,nm)
        !       if(l.eq.1) then
        !       if(nactc(l,nm).gt.1.)write(6,*)"Callmatrix",naero(mkx,nm)*1.e-6
        !        endif
      enddo
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,naero,nmodes,'end',qr0=mrain, &
           nr0=nrain)
#endif
#ifdef TRACERS_TOMAS
!CCC
!Can *************************************************************************
!Can      CLOUD DROPLET CALCULATION
!Can *************************************************************************
!CCC
!CCC *** Input properties for parameterization
!CCC
      TOT_MI    = 0d0
      WPARC      = 0d0
      SMAX       = 0d0
      NACT       = 0d0
      REFF       = 0d0
!      CLDTAU     = 0d0
!      CLDTAUBL   = 0d0 ! I don't account BL case- yhl
!c$$$      QautP6     = 0d0
!c$$$      QautKK     = 0d0
!c$$$      QautMC     = 0d0
!c$$$      QautBH     = 0d0
!c$$$      QautGI     = 0d0
!c$$$      QautNS     = 0d0
!c$$$C
!C Get CCN properties
!C
!      avol(l) = axyp(i_debug,j_debug)*MA(i_debug,j_debug,l)/mair*
!!!   byam(l) = [m2/kg of air]
      boxvl = DXYPIJ*airm(l)*mb2kg*rgas*TL(L)  &
          /100./PL(L)

      CALL getCCN (I_debug,J_debug,L,BOXVL,TOT_MI,TOTi,TPi,MLi, &
         NCCNMx,NSECi)        ! Get CCN properties
!C
!C Call cloud microphysics
!C
      IF (LHX.EQ.LHE) THEN         ! Liquid clouds present
!CCC
!CCC *** Nenes & Seinfeld parameterization - calcuilate droplet number
!CCC

!Two options for Updrate velocity

!1. fixed as a constant

!         WPARCOcean = 0.15d0       ! Fix Ocean and terrestrial updrafts for now
!         WPARCEarth = 0.3d0
!         WPARC      = (1.d0-PEARTH)*WPARCOcean + PEARTH*WPARCEarth

!2. computed using EGCM

        WPARC=v0(mkx) !wturb=sqrt(0.6667*EGCM(l,i,j))
!End of updrate velocity option.


         QLWC=WMX(L)/(FCLD + teeny)     !in-cloud dimensionless LWC
         QLWC=MIN(QLWC, 3.d-03)    !(upper limit for the QLWC)

         RHO=1.d5*PL(L)/(RGAS*TL(L))
         RHOSI = RHO*1.d-3
         if(rhosi.eq.0.) print*,'zero rho',rho,pl(l),tl(l)

         TPARC=tk0(mkx)
         PPARC=pk0(mkx)*100.d0  ! mbar to Pa
         IF (TOTi.GT.6.d7.and.WPARC.gt.0.) THEN  ! more than 60 particles per cc, call droplet activation
            CALL CALCNd (TPARC,PPARC,TPi,MLi,NSECi,WPARC,NACT & ! Activate droplets
                 ,SMAX ,RHOSI,QLWC,EPSILON,AUTO,DIFFLWMR,DIFFEPS,pearth)
         ELSE
!YUNHA- I don't know how to set the minimum NACT. For now, I put the same mininum as Nenes.
            NACT = 40.d6 !ndrop(mkx) ! 40.d6      ! Minimum droplet number [#/m3]
            SMAX = 0.0001    ! Minimum supersaturation
         ENDIF
       NACTL(mkx)=NACT
       CDNC_NENES(L)=NACT !FOR DIAGNOSTICS
       ENDIF

       ldummy=execute_bulk2m_driver('all' &
            ,ndrop,mdrop,ncrys,mcrys,nactl,'end',qr0=mrain, &
           nr0=nrain)

#endif
#ifdef TRACERS_AEROSOLS_Koch
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,'end',qr0=mrain,nr0=nrain)
#endif
      ! Make calls to calculate hydrometeors' growth rates due to microphysical
      ! processes
      ! content      :       [kg water/kg air/s]
      ! concentration:       [No/kg air/s]
      ! Activation of cloud droplets: prescribed AP spectrum
      !*** For the originial HM scheme with fixed distributions for amm. sulfate
      !       ldummy=execute_bulk2m_driver('hugh','drop_nucl',dtB2M,mkx)

      !*** Use this if using the Lohmann or Gultepe scheme  for mass to number
#ifdef TRACERS_AEROSOLS_Koch
      OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
      NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
      ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx, &
           OLDCDNC,NEWCDNC)
#endif
#ifdef TRACERS_AMP
      !*** Using the AMP_actv interface from MATRIX
      ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#endif
#ifdef TRACERS_TOMAS
      !*** Using the TOMAS_actv interface from TOMAS
        ldummy=execute_bulk2m_driver('toma','drop_nucl',dtB2M,mkx)
#endif
      ! Droplets' autoconversion: Beheng (concentration and content)
      !       ldummy=execute_bulk2m_driver('hugh','drop_auto',dtB2M,mkx)
      ! Droplets' autoconversion: Seifert and Beheng (concentration and content)
      !       ldummy=execute_bulk2m_driver('beheng','drop_auto',dtB2M,mkx)
      ! Freezing of cloud droplets (contact and immersion)
      ldummy=execute_bulk2m_driver('hugh','drop_frzn',dtB2M,mkx)
      ! Crystal nucleation: Ncrys=anuc(k)*wef**bnuc(k)
      ldummy=execute_bulk2m_driver('hugh','crys_nucl',dtB2M,mkx)
      ! Numerous processes of water-water,water-ice, ice-ice interaction,
      ! condensation/evaporation/deposition/sublimation, ice multiplication
      ! and sedimention are ready to be called. There are only a few examples:
      !        ldummy=execute_bulk2m_driver('hugh','drop_rain',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','drop_snow',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_auto',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_snow',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_cond',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','snow_melt',dtB2M,mkx)
      !
      ! In this chain of events the very last call, which applies saturation
      ! adjustment to keep environment about water saturation, is supposed to be
      !        ldummy=execute_bulk2m_driver('hugh','drop_cond',dtB2M,mkx)
      !
      ! To ensure calculated growth rates don't lead to negative contents:
      ! Previous call ('drop_cond') HAS to be used.
      !        ldummy=execute_bulk2m_driver('all','make_balance',mkx)
      ! Otherwise, the caller is responsible to handle the problem.
      !
      ! To get particular growth rate:
      !        rablk = execute_bulk2m_driver('get','rate','mprc')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To get tendencies of temperature, water vapor mixing ratio or
      ! tendencies of concentration/content of particular hydrometeor:
      !        rablk = execute_bulk2m_driver('get','tnd','qc_tnd')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To get parameters of hydrometeors size distributions:
      !        rablk = execute_bulk2m_driver('get','val','ec')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To update concentration and contents due to uncommented processes
      ! If "make_balance" was used:
      !        ndrop=ndrop+dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
      !        mdrop=mdrop+dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
      !        ncrys=ncrys+dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
      !        mcrys=mcrys+dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')

      !
      ! This is our case
      ! Without "make_balance":
      ! Droplet concentration
      !       if(l.eq.1) then
      !        npccn=execute_bulk2m_driver('get','npccn')
      !        nprc=execute_bulk2m_driver('get','nprc')
      !        nnuccc=execute_bulk2m_driver('get','nnuccc')
      !        nnucci=execute_bulk2m_driver('get','nnucci')
      !         if(npccn(l).gt.1.)write(6,*)"check BLK",ndrop*1.e-6,
      !    *npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
      !       endif

      ndrop=ndrop+( &
                                !       npccn              ! change n droplets activation
           +execute_bulk2m_driver('get','npccn') &
                                !       nprc               ! change n autoconversion of droplets:
           -execute_bulk2m_driver('get','nprc') &
                                !       nnuccc             ! change n due to con droplets freez
           -execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci             ! change n due to imm droplets freez
           -execute_bulk2m_driver('get','nnucci') &
                                !
           )*dtB2M
      !
      ! Droplet content
      mdrop=mdrop+( &
                                !       mpccn              ! change q droplets activation
           +execute_bulk2m_driver('get','mpccn') &
                                !       mprc               ! change q autoconversion of droplets:
           -execute_bulk2m_driver('get','mprc') &
                                !       mnuccc             ! change q due to con droplets freez
           -execute_bulk2m_driver('get','mnuccc') &
                                !       mnucci             ! change q due to imm droplets freez
           -execute_bulk2m_driver('get','mnucci') &
                                !
           )*dtB2M
      !
      ! Crystal concentration
      ncrys=ncrys+( &
                                !       nnuccc             ! change n due to contact droplets freez
           +execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci             ! change n due to immersion droplets freez
           +execute_bulk2m_driver('get','nnucci') &
                                !       nnuccd             ! change n freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','nnuccd') &
           )*dtB2M &
                                !      nnucmd              ! change n cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmd') &
                                !      nnucmt              ! change n cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmt')
      !
      ! Crystal content
      mcrys=mcrys+( &
                                !       mnuccc             ! change q due to con droplets freez
           +execute_bulk2m_driver('get','mnuccc') &
                                !       mnucci             ! change q due to imm droplets freez
           +execute_bulk2m_driver('get','mnucci') &
                                !       mnuccd             ! change q freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','mnuccd') &
           )*dtB2M &
                                !      mnucmd              ! change q cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','mnucmd') &
                                !      mnucmt              ! change q cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','mnucmt')
      !
      if_balance: if( (ndrop(mkx) .lt. 0) .or. (mdrop(mkx) .lt. 0)  .or. &
           (ncrys(mkx) .lt. 0) .or. (mcrys(mkx) .lt. 0)) then
        !
        if(lSCM) then
          write(6,*)"stop BLK: ndrop_old,mdrop_old,ncrys_old,mcrys_old" &
               ,l,ndrop_old*1.e-6,mdrop_old*1.e+3,ncrys_old*1.e-3,mcrys_old*1.e+3
          !
          write(6,*)"stop BLK: ndrop,mdrop,ncrys,mcrys" &
               ,l,ndrop*1.e-6,mdrop*1.e+3,ncrys*1.e-3,mcrys*1.e+3
        endif
        !
        ! No/m^3
        !
        npccn  =             & ! change n droplets activation
             +execute_bulk2m_driver('get','npccn')*dtB2M
        nprc   =             & ! change n autoconversion of droplets:
             -execute_bulk2m_driver('get','nprc')*dtB2M
        nnuccc =             & ! change n due to con droplets freez
             -execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci =             & ! change n due to imm droplets freez
             -execute_bulk2m_driver('get','nnucci')*dtB2M

        nc_tot = npccn + nprc + nnuccc + nnucci
        !
        ! No/cc
        !
        if(lSCM) then
          write(6,*)"stop BLK: ndrop_old,nc_tot,ndrop" &
               ,l,ndrop_old*1.e-6,nc_tot*1.e-6,ndrop*1.e-6
          write(6,*)"stop BLK: npccn,nprc,nnuccc,nnucci" &
               ,l,npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
        endif
        !
        ! kg/kg
        !
        mpccn  =              & ! change q droplets activation
             execute_bulk2m_driver('get','mpccn')*dtB2M
        mprc   =              & ! change q autoconversion of droplets:
             -execute_bulk2m_driver('get','mprc')*dtB2M
        mnuccc =              & ! change q due to con droplets freez
             -execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci =              & ! change q due to imm droplets freez
             -execute_bulk2m_driver('get','mnucci')*dtB2M

        qc_tot = mpccn + mprc + mnuccc + mnucci
        !
        ! g/kg
        !
        if(lSCM) then
          write(6,*)"stop BLK: mdrop_old,qc_tot,mdrop" &
               ,l,mdrop_old*1.e+3,qc_tot*1.e+3,mdrop*1.e+3
          write(6,*)"stop BLK: mpccn,mprc,mnuccc,mnucci" &
               ,l,mpccn*1.e+3,mprc*1.e+3,mnuccc*1.e+3,mnucci*1.e+3
        endif
        !
        ! No/m^3
        !
        nnuccc  =             & ! change n due to contact droplets freez
             +execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci  =             & ! change n due to immersion droplets freez
             +execute_bulk2m_driver('get','nnucci')*dtB2M
        nnuccd  =             & ! change n freezing aerosol (prim ice nuc)
             +execute_bulk2m_driver('get','nnuccd') ! *dtB2M
        nnucmd =              & ! change n cond freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','nnucmd') ! *dtB2M
        nnucmt =              & ! change n cont freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','nnucmt') ! *dtB2M

        ni_tot = nnuccc + nnucci + nnuccd + nnucmd + nnucmt
        !
        ! No/l
        !
        if(lSCM) then
          write(6,*)"stop BLK: ncrys_old,ni_tot,ncrys" &
               ,l,ncrys_old*1.e-3,ni_tot*1.e-3,ncrys*1.e-3
          write(6,*)"stop BLK: nnuccc,nnucci,nnuccd,nnucmd,nnucmt" &
               ,l,nnuccc*1.e-3,nnucci*1.e-3,nnuccd*1.e-3,nnucmd*1.e-3 &
               ,nnucmt*1.e-3
        endif
        !
        ! kg/kg
        !
        mnuccc  =             & ! change q due to contact droplets freez
             +execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci  =             & ! change q due to immersion droplets freez
             +execute_bulk2m_driver('get','mnucci')*dtB2M
        mnuccd  =             & ! change q freezing aerosol (prim ice nuc)
             +execute_bulk2m_driver('get','mnuccd') ! *dtB2M
        mnucmd =              & ! change q cond freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','mnucmd') ! *dtB2M
        mnucmt =              & ! change q cont freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','mnucmt') ! *dtB2M

        qi_tot = mnuccc + mnucci + mnuccd + mnucmd + mnucmt
        !
        ! g/m^3
        !
        if(lSCM) then
          write(6,*)"stop BLK: mcrys_old,qi_tot,mcrys" &
               ,l,mcrys_old*1.e+3,qi_tot*1.e+3,mcrys*1.e+3
          write(6,*)"stop BLK: mnuccc,mnucci,mnuccd,mnucmd,mnucmt" &
               ,l,mnuccc*1.e+3,mnucci*1.e+3,mnuccd*1.e+3,mnucmd*1.e+3 &
               ,mnucmt*1.e+3
        endif
        !
        ! balanced tendecies:
        !
        ldummy=execute_bulk2m_driver('all','make_balance',mkx)
        !
        nc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:00: nc_tnd",l,nc_tnd*1.e-6
        endif
        !
        npccn  =             & ! change n droplets activation
             +execute_bulk2m_driver('get','npccn')*dtB2M
        nprc   =             & ! change n autoconversion of droplets:
             -execute_bulk2m_driver('get','nprc')*dtB2M
        nnuccc =             & ! change n due to con droplets freez
             -execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci =             & ! change n due to imm droplets freez
             -execute_bulk2m_driver('get','nnucci')*dtB2M

        nc_tnd = npccn + nprc + nnuccc + nnucci
        if(lSCM) then
          write(6,*)"stop BLK:01: nc_tnd",l,nc_tnd*1.e-6
        endif
        !
        qc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:00: qc_tnd",l,qc_tnd*1.e+3
        endif
        !
        mpccn  =              & ! change q droplets activation
             execute_bulk2m_driver('get','mpccn')*dtB2M
        mprc   =              & ! change q autoconversion of droplets:
             -execute_bulk2m_driver('get','mprc')*dtB2M
        mnuccc =              & ! change q due to con droplets freez
             -execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci =              & ! change q due to imm droplets freez
             -execute_bulk2m_driver('get','mnucci')*dtB2M

        qc_tnd = mpccn + mprc + mnuccc + mnucci
        if(lSCM) then
          write(6,*)"stop BLK:01: qc_tnd",l,qc_tnd*1.e+3
        endif
        !
        ni_tnd=dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:01: ni_tnd",l,ni_tnd*1.e-3
        endif
        !
        qi_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:01: qi_tnd",l,qi_tnd*1.e+3
        endif
        !
        ndrop_new = ndrop_old + nc_tnd
        mdrop_new = mdrop_old + qc_tnd
        ncrys_new = ncrys_old + ni_tnd
        mcrys_new = mcrys_old + qi_tnd
        !
        ndrop_blk = execute_bulk2m_driver('get','val','nc')
        mdrop_blk = execute_bulk2m_driver('get','val','qc')
        ncrys_blk = execute_bulk2m_driver('get','val','ni')
        mcrys_blk = execute_bulk2m_driver('get','val','qi')
        !
        ndrop_res = ndrop_blk + nc_tnd
        mdrop_res = mdrop_blk + qc_tnd
        ncrys_res = ncrys_blk + ni_tnd
        mcrys_res = mcrys_blk + qi_tnd
        !
        if(wSCM) then
          write(6,*) &
               "stop BLK: ndrop_old,nc_tnd,ndrop_new" &
               ,l,ndrop_old*1.e-6,nc_tnd*1.e-6,ndrop_new*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_old,qc_tnd,mdrop_new" &
               ,l,mdrop_old*1.e+3,qc_tnd*1.e3,mdrop_new*1.e+3
          !
          write(6,*) &
               "stop BLK: ncrys_old,ni_tnd,ncrys_new" &
               ,l,ncrys_old*1.e-3,ni_tnd*1.e-3,ncrys_new*1.e-3
          !
          write(6,*) &
               "stop BLK: mcrys_old,qi_tnd,mcrys_new" &
               ,l,mcrys_old*1.e+3,qi_tnd*1.e3,mcrys_new*1.e+3
          !
          write(6,*) &
               "stop BLK: ndrop_blk,nc_tnd,ndrop_res" &
               ,l,ndrop_blk*1.e-6,nc_tnd*1.e-6,ndrop_res*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_blk,qc_tnd,mdrop_res" &
               ,l,mdrop_blk*1.e+3,qc_tnd*1.e3,mdrop_res*1.e+3
          !
          write(6,*) &
               "stop BLK: ndrop_old,ndrop,ndrop_new,nc_tot,nc_tnd" &
               ,l,ndrop_old*1.e-6,ndrop_new*1.e-6,ndrop*1.e-6 &
               ,nc_tot*1.e-6,nc_tnd*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_old,mdrop,mdrop_new,qc_tot,qc_tnd" &
               ,l,mdrop_old*1.e+3,mdrop_new*1.e+3,mdrop*1.e+3 &
               ,qc_tot*1.e+3,qc_tnd*1.e+3
        endif
        !
        ! output for standalone internal variables
        !
        if(lSCM) then
          write(iuo,*) 'l dtB2M'
          write(iuo,*) l
          write(iuo,*) dtB2M

          write(iuo,*) 'wmx wmxice tl ql pl svlhxl lhx '
          write(iuo,*) wmx
          write(iuo,*) wmxice
          write(iuo,*) tl
          write(iuo,*) ql
          write(iuo,*) pl
          write(iuo,*) svlhxl
          write(iuo,*) lhx

          write(iuo,*) 'ndrop_old mdrop_old ncrys_old mcrys_old'
          write(iuo,*) ndrop_old
          write(iuo,*) mdrop_old
          write(iuo,*) ncrys_old
          write(iuo,*) mcrys_old

          write(iuo,*) 'nc_tnd qc_tnd ni_tnd qi_tnd'
          write(iuo,*) nc_tnd
          write(iuo,*) qc_tnd
          write(iuo,*) ni_tnd
          write(iuo,*) qi_tnd

          write(iuo,*) 'ndrop_new mdrop_new ncrys_new mcrys_new'
          write(iuo,*) ndrop_new
          write(iuo,*) mdrop_new
          write(iuo,*) ncrys_new
          write(iuo,*) mcrys_new

          write(iuo,*) 'ndrop mdrop ncrys mcrys'
          write(iuo,*) ndrop
          write(iuo,*) mdrop
          write(iuo,*) ncrys
          write(iuo,*) mcrys

          write(iuo,*) 'ndrop_res mdrop_res ncrys_res mcrys_res'
          write(iuo,*) ndrop_res
          write(iuo,*) mdrop_res
          write(iuo,*) ncrys_res
          write(iuo,*) mcrys_res

          write(iuo,*) 'npccn nprc nnuccc nnucci nc_tot'
          write(iuo,*) npccn
          write(iuo,*) nprc
          write(iuo,*) nnuccc
          write(iuo,*) nnucci
          write(iuo,*) nc_tot

          write(iuo,*) 'mpccn mprc mnuccc mnucci qc_tot'
          write(iuo,*) mpccn
          write(iuo,*) mprc
          write(iuo,*) mnuccc
          write(iuo,*) mnucci
          write(iuo,*) qc_tot

          write(iuo,*) 'nnuccd nnucmd nnucmt ni_tot'
          write(iuo,*) nnuccd
          write(iuo,*) nnucmd
          write(iuo,*) nnucmt
          write(iuo,*) ni_tot

          write(iuo,*) 'mnuccd mnucmd mnucmt qi_tot'
          write(iuo,*) mnuccd
          write(iuo,*) mnucmd
          write(iuo,*) mnucmt
          write(iuo,*) qi_tot
          !
        endif
        !
        if( (ndrop_res(mkx) .lt. 0) .or. (mdrop_res(mkx) .lt. 0)  .or. &
             (ncrys_res(mkx) .lt. 0) .or. (mcrys_res(mkx) .lt. 0)) then
          !        call stop_model("BLK2MOM: Negative conc/cont...", 255)
          !         write(6,*)"We reached -ve con.",ndrop_res(mkx),mdrop_res(mkx),
          !     * ncrys_res(mkx), mcrys_res(mkx),l
          ndrop_res(mkx)=20.*1.d06
          mdrop_res(mkx)=1*1.d-06
          ncrys_res(mkx)=1*1.d-06
          mcrys_res(mkx)=1*1.d02
        else
          ndrop=ndrop_res;mdrop=mdrop_res;ncrys=ncrys_res;mcrys=mcrys_res
        endif
        !
      endif if_balance
      !
      ! To calculate "new" temperature and vapor mixing ratio:
      ldummy=execute_bulk2m_driver('tkqv','tk_qv',mkx)
      tk0new=tk0+dtB2M*execute_bulk2m_driver('get','tnd','tk_tnd')
      qk0new=qk0+dtB2M*execute_bulk2m_driver('get','tnd','qv_tnd')
      !
      ! At this point you have 2 phases separately.
      ! Almost all processes are switched off, but you can calculate also
      ! accreation of droplets/crystal  by rain/snow, for example, and use
      ! rain/snow as diagnostic variables. But you need one additional
      ! long-storage array to keep ice crystal hydrometeor content as a minimum
      !
      !     IF(LHX.EQ.LHE)  THEN
      !        WMX(L)=mdrop(mkx)
      !      ELSE
      !        WMX(L)=mcrys(mkx)
      !      ENDIF
      !
      ! GCM logics...........  SNd0, SNdL [No/cc; ]SNdI Units also in /cc
      !
      !      SNdI=ncrys(mkx)*1.0d-6          ! ncrys, [No/m^3]
      SNdI = 0.06417127d0
      !      if(SNdI.gt.0.) write(6,*)"ICE CRY",SNdI, SNdI/dtB2M
      if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
      SNd=ndrop(mkx)*1.d-6                 ! ndrop, [No/m^3]
      !      if(SNd.gt.20.) write(6,*)"SM 12 CDNC",SNd   ,l
      !**** Old treatment to get CDNC for cloud changes within time steps
      DCLD(L) = FCLD-CLDSAVL(L) ! cloud fraction change
      !** If previous time step is clear sky
      if(CLDSAVL(L).eq.0.) then
        SNd=SNd
      elseif (DCLD(L).le.0.d0) then
        SNd=OLDCDL(L)
      elseif(DCLD(L).gt.0.d0) then
        SNd=( (OLDCDL(L)*CLDSAVL(L)) + (SNd*DCLD(L)) )/FCLD
      endif
      !* If using an alternate definition for QAUT
      rablk=execute_bulk2m_driver('get','mprc')
      QAUT_B2M=rablk(mkx)
#endif  /* (cld-aer and BLK_2MOM) */
#ifdef CLD_AER_CDNC
      SCDNCW=SNd      ! we have already passed the grid box value
      SCDNCI=SNdI
      if (SCDNCI.le.0.0d0) SCDNCI=teeny         !set min ice crystal, do we need this, please check
      if (SCDNCW.le.20.d0) SCDNCW=20.d0         !set min CDNC, sensitivity test
      !     if(SCDNCW.gt.2000.) write(6,*)"PROBLEM",SCDNCW,L
      if (SCDNCW.ge.1400.d0) SCDNCW=1400.d0     !set max CDNC, sensitivity test
      !     if (SNd.gt.20.) write(6,*)"CDNC LSS",SCDNCW,SNd,L
#endif
      !**** COMPUTE THE AUTOCONVERSION RATE OF CLOUD WATER TO PRECIPITATION
      if(WMX(L).gt.0.) then
        RHO=1d5*PL(L)/(RGAS*TL(L))
        TEM=RHO*WMX(L)/(WCONST*FCLD+teeny)
        if(LHX.eq.LHS ) TEM=RHO*WMX(L)/(WMUI*FCLD+teeny)
        !       IF(LHX.EQ.LHE.AND.ROICE.GT..1d0) TEM=RHO*WMX(L)
        !    *    /(WMUSI*FCLD+teeny)
        TEM=TEM*TEM
        if(TEM.gt.10.) TEM=10.
        if(VDEF.gt.0..and.RHO*WMX(L).ge.10.0d0) CM0=CM00
        CM1=CM0
        if(BANDF) CM1=CM0*CBF
        if(LHX.eq.LHS) CM1=CM0
        CM=CM1*(1.-1./exp(TEM*TEM))+100.*(PREBAR(L+1)+ &
             PRECNVL(L+1)*BYDTsrc)
        !C#ifdef CLD_AER_CDNC
        !** Choice of 2 different routines to get the autoconversion rate
        !C#ifdef BLK_2MOM
        !*** using an alternate QAUT definition based on Beheng (1994)
        !C    if(FCLD.gt.teeny) then
        !C      CM=QAUT_B2M/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
        !C   *  PRECNVL(L+1)*BYDTsrc)
        !       if (QAUT_B2M.lt.0.) write(6,*)"QAUT BLK_2M",QAUT_B2M,CM,WMX(L),L
        !       if(L.eq.1) write(6,*)"4th check BLK_2M",CM,QAUT_B2M,WMX(L)
        !C      if(CM.gt.1.d-03) CM=1.d-03
        !C    else
        !C      CM=0.d0
        !C    endif
        !C#else
        !** Use Qaut definition based on Rotstayn and Liu (2005, GRL)
        !C         WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
        !C          IF(LHX.EQ.LHE)  THEN
        !C            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
        !C          ELSE
        !C            RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
        !    *         *(1.+pl(l)*xRICld)
        !C          END IF
        !C       CALL GET_QAUT(L,PL(L),TL(L),FCLD,WMX(L),SCDNCW,RCLD,RHOW,
        !C     *r6,r6c,QCRIT,QAUT)
        !** Can also use other Qaut definitions if BLK_2MOM is not defined by switching to this call
        !     CALL GET_QAUT(L,TL(L),FCLD,WMX(L),SCDNCW,RHO,QCRIT,QAUT)
        !      CALL GET_QAUT(L,FCLD,WMX(L),SCDNCW,RHO,QAUT)
        !*** If 6th moment of DSD is greater than critical radius r6c start QAUT
        !     if (r6.gt.r6c) then
        !C      if ((WMX(L)/(FCLD+teeny)).GT.QCRIT) then
        !C        CM=QAUT/(WMX(L)+1.d-20)+1.d0*100.d0*(PREBAR(L+1)+
        !C     *     PRECNVL(L+1)*BYDTsrc)
        !C      else
        !C        CM=0.d0
        !C      endif
        !** end routine for QAUT as a function of N,LWC
        !C#endif
        !C#endif
        CM=CM*CMX
        if(CM.gt.BYDTsrc) CM=BYDTsrc
        PREP(L)=WMX(L)*CM
#ifdef SCM
        if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
          PRESAV(L)=PREP(L)*DTsrc
        endif
#endif
        if(TL(L).lt.TF.and.LHX.eq.LHE) then ! check snowing pdf
          PRATM=1d5*COEFM*WMX(L)*PL(L)/(WCONST*FCLD*TL(L)*RGAS+teeny)
          PRATM=min(PRATM,1d0)*(1.-exp(max(-1d2,(TL(L)-TF)/COEFT)))
          if(PRATM.gt.RNDSSL(3,L)) LHP(L)=LHS
        end if
      else
        CM=0.
      end if
      !**** DECIDE WHETHER TO FORM CLOUDS
      !**** FORM CLOUDS ONLY IF RH GT RH00
      if (RH1(L).lt.RH00(L)) then
        FORM_CLOUDS=.false.
      else   ! COMPUTE THE CONVERGENCE OF AVAILABLE LATENT HEAT
        SQ(L)=LHX*QSATL(L)*DQSATDT(TL(L),LHX)*BYSHA
        TEM=-LHX*DPDT(L)/PL(L)
        QCONV=LHX*AQ(L)-RH(L)*SQ(L)*SHA*PLK(L)*ATH(L) &
             -TEM*QSATL(L)*RH(L)
        FORM_CLOUDS= (QCONV.gt.0. .or. WMX(L).gt.0.)
      end if
      !****
      ERMAX=LHX*PREBAR(L+1)*GRAV*BYAM(L)
      if (FORM_CLOUDS) then
        !**** COMPUTE EVAPORATION OF RAIN WATER, ER
        RHN=min(RH(L),RHF(L))
        if(WMX(L).gt.0.)  then
          ER(L)=(1.-RHN)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        else                    !  WMX(l).le.0.
          if(PREICE(L+1).gt.0..and.TL(L).lt.TF)  then
            ER(L)=(1.-RHI)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
          else
            ER(L)=(1.-RH(L))**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
          end if
        end if
        ER(L)=max(0d0,min(ER(L),ERMAX))
        !**** COMPUTATION OF CLOUD WATER EVAPORATION
        if (CLEARA(L).gt.0.) then
          WTEM=1d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
          if(LHX.eq.LHE.and.WMX(L)/FCLD.ge.WCONST*1d-3) &
               WTEM=1d2*WCONST*PL(L)/(TL(L)*RGAS)
          if(WTEM.lt.1d-10) WTEM=1d-10
          if(LHX.eq.LHE)  then
            !           RCLD=1d-6*(RWCLDOX*10.*(1.-PEARTH)+7.*PEARTH)*(WTEM*4.)**BY3
            RCLD=RCLDX*1d-6*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
          else
            !           RCLD=25.d-6*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
            RCLD=RCLDX*100.d-6*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
            RCLD=min(RCLD,RIMAX)
          end if
          CK1=1000.*LHX*LHX/(2.4d-2*RVAP*TL(L)*TL(L))
          CK2=1000.*RGAS*TL(L)/(2.4d-3*QSATL(L)*PL(L))
          TEVAP=COEEC*(CK1+CK2)*RCLD*RCLD
          WMX1=WMX(L)-PREP(L)*DTsrc
          WMPR(L)=PREP(L)*DTsrc             ! precip water
          ECRATE=(1.-RHF(L))/(TEVAP*FCLD+teeny)
          if(ECRATE.gt.BYDTsrc) ECRATE=BYDTsrc
          EC(L)=WMX1*ECRATE*LHX
        end if
        !**** COMPUTE NET LATENT HEATING DUE TO STRATIFORM CLOUD PHASE CHANGE,
        !**** QHEAT, AND NEW CLOUD WATER CONTENT, WMNEW
        DRHDT=2.*CLEARA(L)*CLEARA(L)*(1.-RH00(L))*(QCONV+ER(L))/LHX/ &
             (WMX(L)/(FCLD+teeny)+2.*CLEARA(L)*QSATL(L)*(1.-RH00(L)) &
             +teeny)
        if(ER(L).eq.0.and.WMX(L).le.0.) DRHDT=0.
        QHEAT(L)=FSSL(L)*(QCONV-LHX*DRHDT*QSATL(L))/(1.+RH(L)*SQ(L))
        if (debug) print*,"ls2",l,qheat(l)
        DWDT=QHEAT(L)/LHX-PREP(L)+CLEARA(L)*FSSL(L)*ER(L)/LHX
        WMNEW =WMX(L)+DWDT*DTsrc
        if(WMNEW.lt.0.) then
          WMNEW=0.
          QHEAT(L)=(-WMX(L)*BYDTsrc+PREP(L))*LHX-CLEARA(L)*FSSL(L)*ER(L)
          if (debug) print*,"ls3",l,qheat(l)
        end if
      else
        !**** UNFAVORABLE CONDITIONS FOR CLOUDS TO EXIT, PRECIP OUT CLOUD WATER
        QHEAT(L)=0.
        if (WMX(L).gt.0.) then

          call get_dq_evap(tl(l)*rh00(l)/plk(l),ql(l)*rh00(l),plk(l) &
               ,rh00(l),lhx,pl(l),wmx(l)/(fssl(l)*rh00(l)),dqsum,fqcond1 &
               )
          DWDT=DQSUM*RH00(L)*FSSL(L)

          !**** DWDT is amount of water going to vapour, store LH (sets QNEW below)
          QHEAT(L)=-DWDT*LHX*BYDTsrc
          if (debug) print*,"ls4",l,qheat(l)
          PREP(L)=max(0d0,(WMX(L)-DWDT)*BYDTsrc) ! precip out cloud water
          WMPR(L)=PREP(L)*DTsrc ! precip water (for opt. depth calculation)
#ifdef SCM
          if (i_debug.eq.I_TARG .and. j_debug.eq.J_TARG) then
            PRESAV(L)=PREP(L)*DTsrc
          endif
#endif
        end if
        ER(L)=(1.-RH(L))**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        if(PREICE(L+1).gt.0..and.TL(L).lt.TF) &
             ER(L)=(1.-RHI)**ERP*LHX*PREBAR(L+1)*GbyAIRM0 ! GRAV/AIRM0
        ER(L)=max(0d0,min(ER(L),ERMAX))
        QHEAT(L)=QHEAT(L)-CLEARA(L)*FSSL(L)*ER(L)
        if (debug) print*,"ls5",l,qheat(l),lhx
        WMNEW=0.
      end if

      !**** PHASE CHANGE OF PRECIPITATION, FROM ICE TO WATER
      !**** This occurs if 0 C isotherm is crossed, or ice is falling into
      !**** a super-cooled water cloud (that has not had B-F occur).
      !**** Note: on rare occasions we have ice clouds even if T>0
      !**** In such a case, no energy of phase change is needed.
      HPHASE=0.
      if (LHP(L+1).eq.LHS.and.LHP(L).eq.LHE.and.PREICE(L+1).gt.0) then
        HPHASE=HPHASE+LHM*PREICE(L+1)*GRAV*BYAM(L)
        PREICE(L+1)=0.
        if (debug) print*,"ls6",l,hphase,lhp(l),lhp(l+1),lhe
      end if
      !**** PHASE CHANGE OF PRECIP, FROM WATER TO ICE
      if (LHP(L+1).eq.LHE.and.LHP(L).eq.LHS.and.PREBAR(L+1).gt.0) then
        HPHASE=HPHASE-LHM*PREBAR(L+1)*GRAV*BYAM(L)
        if (debug) print*,"ls7",l,hphase,lhp(l),lhp(l+1),lhe
      end if

      !**** Make sure energy is conserved for transfers between P and CLW
      if (LHP(L).ne.LHX) then
        HPHASE=HPHASE+(ER(L)*CLEARA(L)*FSSL(L)/LHX-PREP(L))*LHM
        if (debug) print*,"ls8",l,hphase
      end if
      !**** COMPUTE THE PRECIP AMOUNT ENTERING THE LAYER TOP
      if (ER(L).eq.ERMAX) then ! to avoid round off problem
        PREBAR(L)=PREBAR(L+1)*(1.-CLEARA(L)*FSSL(L))+ &
             AIRM(L)*PREP(L)*BYGRAV
      else
        PREBAR(L)=max(0d0,PREBAR(L+1)+ &
             AIRM(L)*(PREP(L)-ER(L)*CLEARA(L)*FSSL(L)/LHX)*BYGRAV)
      end if
      !**** UPDATE NEW TEMPERATURE AND SPECIFIC HUMIDITY
      QNEW =QL(L)-DTsrc*QHEAT(L)/(LHX*FSSL(L)+teeny)
      if(QNEW.lt.0.) then
        QNEW=0.
        QHEAT(L)=QL(L)*LHX*BYDTsrc*FSSL(L)
        if (debug) print*,"ls9",l,qheat(l)
        DWDT1=QHEAT(L)/LHX-PREP(L)+CLEARA(L)*FSSL(L)*ER(L)/LHX
        WMNEW=WMX(L)+DWDT1*DTsrc
        !**** IF WMNEW .LT. 0., THE COMPUTATION IS UNSTABLE
        if(WMNEW.lt.0.) then
          IERR=1
          LERR=L
          WMERR=WMNEW
          WMNEW=0.
        end if
      end if
      !**** Only Calculate fractional changes of Q to W
#ifdef TRACERS_WATER
      FPR=0.
      if (WMX(L).gt.0.) FPR=PREP(L)*DTsrc/WMX(L)              ! CLW->P
      FPR=min(1d0,FPR)
      FER=0.
      if (PREBAR(L+1).gt.0.) FER=CLEARA(L)*FSSL(L)*ER(L)*AIRM(L)/ &
           (GRAV*LHX*PREBAR(L+1))                             ! P->Q
      FER=min(1d0,FER)
      FWTOQ=0.                                                ! CLW->Q
#endif
      FQTOW=0.                                                ! Q->CLW
      if (FSSL(L).gt.0) then
        if (QHEAT(L)+CLEARA(L)*FSSL(L)*ER(L).gt.0) then
          if (LHX*QL(L)+DTsrc*CLEARA(L)*ER(L).gt.0.) FQTOW=(QHEAT(L &
               )+CLEARA(L)*FSSL(L)*ER(L))*DTsrc/((LHX*QL(L)+DTsrc*CLEARA(L) &
               *ER(L))*FSSL(L))
#ifdef TRACERS_WATER
        else
          if (WMX(L)-PREP(L)*DTsrc.gt.0.) FWTOQ=-(QHEAT(L) &
               +CLEARA(L)*FSSL(L)*ER(L))*DTsrc/(LHX*(WMX(L)-PREP(L)*DTsrc))
          FWTOQ=min(1d0,FWTOQ)
#endif
        end if
      end if
      QL(L)=QNEW
      !**** adjust gradients down if Q decreases
      QMOM(:,L)= QMOM(:,L)*(1.-FQTOW)
      WMX(L)=WMNEW
      !     if(abs(DTsrc*(QHEAT(L)-HPHASE)).gt.100.*SHA*FSSL(L)) then ! warn
      !       write(0,*) 'it,i,j,l,tlold,dtl,qht,hph',itime,i_debug,j_debug,
      !    *    L,TL(L),DTsrc*(QHEAT(L)-HPHASE)/(SHA*FSSL(L)+teeny),
      !    *    QHEAT(L),HPHASE
      !       debug_out=.true.
      !     end if
      TL(L)=TL(L)+DTsrc*(QHEAT(L)-HPHASE)/(SHA*FSSL(L)+teeny)
      if (debug) print*,"lsA",l,qheat(l),hphase,lhx

      TH(L)=TL(L)/PLK(L)
      QSATC=QSAT(TL(L),LHX,PL(L))
      RH(L)=QL(L)/QSATC
#ifdef TRACERS_WATER
      !**** update tracers from cloud formation (in- and below-cloud
      !****    precipitation, evaporation, condensation, and washout)
      ! CLDSAVT is current FCLD
      if(RH(L).le.1.) then
        if (RH00(L).lt.1.) then
          CLDSAVT=1.-DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        else
          CLDSAVT=0.
        end if
      end if
      if(CLDSAVT.lt.0.) CLDSAVT=0.
      if(RH(L).gt.1.) CLDSAVT=1.
      if (CLDSAVT.gt.1.) CLDSAVT=1.
      if (WMX(L).le.0.) CLDSAVT=0.
      CLDSAVT=CLDSAVT*FSSL(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      WA_VOL=0.
      if (WMNEW.gt.teeny) then
        WA_VOL=WMNEW*AIRM(L)*1.D2*BYGRAV*DXYPIJ
      end if
      WMXTR = WMX(L)
      if (BELOW_CLOUD.and.WMX(L).lt.teeny) then
        precip_mm = PREBAR(L+1)*100.*DTsrc
        if (precip_mm.lt.0.) precip_mm=0.
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (wmxtr.lt.0.) wmxtr=0.
        WA_VOL=precip_mm*DXYPIJ
      end if

      call GET_SULFATE(L,TL(L),FCLD,WA_VOL &
           ,WMXTR,SULFIN,SULFINC,SULFOUT,TR_LEFT,TM,TRWML(1,l),AIRM,LHX &
           ,DT_SULF_SS(1,L),CLDSAVT)

      do iaqch=1,aqchem_count
        n = aqchem_list(iaqch)
        TR_LEF(n)=TR_LEFT(N)
        TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
        TM(L,N)=TM(L,N)*(1.+SULFIN(N))
        TMOM(:,L,N)  = TMOM(:,L,N)*(1. +SULFIN(N))
        if (WMX(L).lt.teeny.and.BELOW_CLOUD) then
          TRPRBAR(N,L+1)=TRPRBAR(N,L+1)+SULFOUT(N)
        else
          TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
        endif
      enddo

#endif

      ! precip. tracer evap
      FPRT = FPR
      DTPRT(1:NTX) = FPRT  *TRWML(1:NTX,L)

      if(fer.ne.0.) then
        call GET_EVAP_FACTOR(NTX,TL(L),LHX,.false.,1d0,FER,FERT,ntix)
        DTERT(1:NTX) = FERT(1:NTX)  *TRPRBAR(1:NTX,L+1)
      else
        FERT(1:NTX) = 0.
        DTERT(1:NTX) = 0.
      endif

      if(fwtoq.ne.0.) then
        call GET_EVAP_FACTOR(NTX,TL(L),LHX,.false.,1d0,FWTOQ,FWTOQT,ntix)
        FPRT=FPR
        do N=1,NTX
          DTQWT(N) = -FWTOQT(N)*TRWML(N,L)*(1.-FPRT)
        enddo
      else
        FWTOQT(1:NTX) = 0.
        DTQWT(1:NTX) = 0.
      endif

      TM_dum(:) = TM(L,:)
!TOMAS DEBUG
            DO N=1,NTM
              if(TM_dum(n).lt.0.) print*,'TM_dum<0 2',TM_dum(n),trname(n)
            ENDDO
!TOMAS DEBUG
      if(BELOW_CLOUD.and.WMX(L).lt.teeny) then
        FQTOWT(:)=0.
        THLAW(gases_list)=0.
        precip_mm = PREBAR(L+1)*100.*dtsrc
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (precip_mm.lt.0.) precip_mm=0.
        if (wmxtr.lt.0.) wmxtr=0.
        call GET_WASH_FACTOR(NTX,b_beta_DT,precip_mm,FWASHT, &
             tl(l),LHX,WMXTR,cldprec,TM_dum,TRPRBAR(:,l), &
             THWASH,pl(l),ntix,.true. &
#ifdef TRACERS_TOMAS
       ,i_debug,j_debug,L &
#endif
            )!washout
      else
        !         b_beta_DT is needed at the lowest precipitating level,
        !         so saving it here for below cloud case:
        b_beta_DT = cldsavt  !*CM*dtsrc
        !         saves cloud fraction at lowest precipitating level for washout
        cldprec=cldsavt
        !dmk added arguments above; THLAW added below (no way to factor this)
        WMXTR = WMX(L)
        call GET_COND_FACTOR( &
             NTX,WMXTR,TL(L),TL(L),LHX,FCLD,FQTOW &
             ,FQTOWT,.false.,TRWML(:,L),TM_dum,THLAW,TR_LEF,PL(L) &
             ,ntix,CLDSAVT)
#ifdef TRACERS_AEROSOLS_OCEAN
        if (trm(i_debug,j_debug,L,n_ococean) .gt. 0.d0) then
          do n=1,ntx
            select case (trname(ntix(n)))
            case ('seasalt1', 'OCocean')
              fqcondt(n)=fqcondt(n)* &
                   trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1) &
                   /(trm(i_debug,j_debug,L,n_ococean)/trpdens(n_ococean)+ &
                   trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1))
            end select
          enddo
        endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
        do N=1,NTX
          DTQWT(N) = DTQWT(N) &
               +FQTOWT(N)*TR_LEF(N)*(TM_dum(N)+DTERT(N))
        enddo
#ifdef NO_WASHOUT_IN_CLOUDS
        THWASH(:)=0.
        FWASHT(:)=0.
#endif
      endif

#ifndef NO_WASHOUT_IN_CLOUDS
      !**** washout in clouds
      ! apply certain removal processes before calculating washout in clouds
      tm_dum(1:ntx)=tm_dum(1:ntx)-dtqwt(1:ntx)-thlaw(1:ntx)
      if(.not.(BELOW_CLOUD.and.WMX(L).lt.teeny)) then
        precip_mm = PREBAR(L+1)*100.*dtsrc
        WMXTR = PREBAR(L+1)*grav*BYAM(L)*dtsrc
        if (precip_mm.lt.0.) precip_mm=0.
        if (wmxtr.lt.0.) wmxtr=0.
        call GET_WASH_FACTOR(NTX,b_beta_DT,precip_mm,FWASHT, &
             tl(l),LHX,WMXTR,cldprec,TM_dum,TRPRBAR(:,l), &
             THWASH,pl(l),ntix,.true. &
#ifdef TRACERS_TOMAS
       ,i_debug,j_debug,L &
#endif
            )!washout
      end if
#endif

#ifdef TRDIAG_WETDEPO
      if (diag_wetdep == 1) then
        FPRT=FPR
        do N=1,NTX
          trevap_ls(l,n)=dtert(n)
          trwash_ls(l,n)=fwasht(n)*tm_dum(n)+thwash(n)
          trclwc_ls(l,n)=fqtowt(n)*tr_lef(n)*(tm(l,n)+dtert(n))+thlaw(n)
          trprcp_ls(l,n)=dtprt(n)
          trclwe_ls(l,n)=fwtoqt(n)*trwml(n,l)*(1.-fprt)
        enddo
      endif
#endif

      do igas=1,gases_count
        n = gases_list(igas)
        if (TM(L,N).gt.teeny) then
          TMFAC(N)=THLAW(N)/TM(L,N)
        else
          TMFAC(N)=0.
        end if
#ifndef NO_WASHOUT_IN_CLOUDS
        if (tm_dum(n).gt.teeny) then
          TMFAC2(N)=THWASH(N)/tm_dum(n) ! use tm instead?
        else
          TMFAC2(N)=0.
        end if
#endif
      enddo

      do N=1,NTX

        dtwrt=fwasht(n)*tm_dum(n)

        ! ---------------------- apply fluxes ------------------------
        FPRT=FPR
        TRWML(N,L) = TRWML(N,L)*(1.-FPRT)  + DTQWT(N)+THLAW(N)

        TM(L,N) = max(0., TM(L,N) &
             + DTERT(N) - DTWRT - DTQWT(N) - THLAW(N) - THWASH(N) )

        TRPRBAR(N,L)=TRPRBAR(N,L+1)*(1.-FERT(N)) &
             +DTPRT(N)+DTWRT+THWASH(N)

        TMOM(:,L,N)  = TMOM(:,L,N)* &
             (1. - FQTOWT(N) - FWASHT(N) - TMFAC(N) - TMFAC2(N))

#ifdef TRACERS_SPECIAL_O18
        !**** Isotopic equilibration of the CLW and water vapour
        if (LHX.eq.LHE .and. WMX(L).gt.0) then  ! only if liquid
          call ISOEQUIL(NTIX(N),TL(L),.true.,QL(L)*FSSL(L),WMX(L), &
               TM(L,N),TRWML(N,L),1d0)
        end if
        !**** Isotopic equilibration of Precip (if liquid) and water vapour
        !**** Note that precip is either all water or all ice
        if (LHP(L).eq.LHE .and. PREBAR(L).gt.0 .and. QL(L).gt.0) then
          PRLIQ=PREBAR(L)*DTSrc*BYAM(L)*GRAV
          call ISOEQUIL(NTIX(N),TL(L),.true.,QL(L)*FSSL(L),PRLIQ, &
               TM(L,N),TRPRBAR(N,L),1d0)
        end if
#endif
      end do
      if (PREBAR(L).eq.0) TRPRBAR(1:ntx,L)=0. ! remove round off error
      if (WMX(L).eq.0)    TRWML(1:ntx,L)=0.   ! remove round off error
#endif
      !**** CONDENSE MORE MOISTURE IF RELATIVE HUMIDITY .GT. 1
      RH1(L)=QL(L)/QSATC
      if(LHX.eq.LHS) then
        if(RH(L).le.1.) then
          if (RH00(L).lt.1.) then
            CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
          else
            CLEARA(L)=1.
          end if
        end if
        if(CLEARA(L).gt.1.) CLEARA(L)=1.
        if(RH(L).gt.1.) CLEARA(L)=0.
        if(WMX(L).le.0.) CLEARA(L)=1.
        QF=(QL(L)-QSATC*(1.-CLEARA(L)))/(CLEARA(L)+teeny)
        if(QF.lt.0.) write(6,*) 'L CA QF Q QSA=',L,CLEARA(L),QF,QL(L), &
             QSATC
        QSATE=QSAT(TL(L),LHE,PL(L))
        RHW=(2.583d0-TL(L)/207.83)*(QSAT(TL(L),LHS,PL(L))/QSATE)
        if(TL(L).lt.238.16.and.WMX(L).le.0d0) RH1(L)=QF/(QSATE*RHW)
      end if
      if(RH1(L).gt.1.) then    ! RH was used in old versions
        SLH=LHX*BYSHA

        call get_dq_cond(tl(l),ql(l),1d0,1d0,lhx,pl(l),dqsum,fcond)

        if(DQSUM.gt.0.) then
          if (debug) print*,"lsB",l,slh*dqsum

          TL(L)=TL(L)+SLH*DQSUM
          QL(L)=QL(L)-DQSUM
          WMX(L)=WMX(L)+DQSUM*FSSL(L)
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
          WA_VOL=0.
          if (WMX(L).gt.teeny) then
            WA_VOL=WMX(L)*AIRM(L)*1.D2*BYGRAV*DXYPIJ
          end if
#endif
          !**** adjust gradients down if Q decreases
          QMOM(:,L)= QMOM(:,L)*(1.-FCOND)
#ifdef TRACERS_WATER
          !**** CONDENSING MORE TRACERS
          WMXTR = WMX(L)
          if(RH(L).le.1.) then
            if (RH00(L).lt.1.) then
              CLDSAVT=1.-DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
            else
              CLDSAVT=0.
            end if
          end if
          if(CLDSAVT.lt.0.) CLDSAVT=0.
          if(RH(L).gt.1.) CLDSAVT=1.
          if (CLDSAVT.gt.1.) CLDSAVT=1.
          if (WMX(L).le.0.) CLDSAVT=0.
          CLDSAVT=CLDSAVT*FSSL(L)
          !dmks  I took out some code above this that was for below cloud
          !   processes - this should be all in-cloud
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

          call GET_SULFATE(L,TL(L),FCLD,WA_VOL,WMXTR,SULFIN, &
               SULFINC,SULFOUT,TR_LEFT,TM,TRWML(1,L),AIRM,LHX, &
               DT_SULF_SS(1,L),CLDSAVT)

          do iaqch=1,aqchem_count
            n = aqchem_list(iaqch)
            TRWML(N,L)=TRWML(N,L)*(1.+SULFINC(N))
            TM(L,N)=TM(L,N)*(1.+SULFIN(N))
            TMOM(:,L,N) =TMOM(:,L,N)*(1.+SULFIN(N))
            TRWML(N,L) = TRWML(N,L)+SULFOUT(N)
            TR_LEF(N)=TR_LEFT(N)
          enddo

#endif
          ! below TR_LEFT(N) limits the amount of available tracer in gridbox
          !dmkf and below, extra arguments for GET_COND, addition of THLAW
          TM_dum(:) = TM(L,:)
!TOMAS DEBUG
            DO N=1,NTM
              if(TM_dum(n).lt.0.) print*,'TM_dum<0 3',TM_dum(n),trname(n)
            ENDDO
!TOMAS DEBUG
          call GET_COND_FACTOR(NTX,WMXTR,TL(L),TL(L),LHX,FCLD,FCOND &
               ,FQCONDT,.false.,TRWML(:,L),TM_dum,THLAW,TR_LEF,pl(l) &
               ,ntix,CLDSAVT)
#ifdef TRACERS_AEROSOLS_OCEAN
          if (trm(i_debug,j_debug,L,n_ococean) .gt. 0.d0) then
            do n=1,ntx
              select case (trname(ntix(n)))
              case ('seasalt1', 'OCocean')
                fqcondt(n)=fqcondt(n)* &
                     trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1) &
                     /(trm(i_debug,j_debug,L,n_ococean)/trpdens(n_ococean)+ &
                     trm(i_debug,j_debug,L,n_seasalt1)/trpdens(n_seasalt1))
              end select
            enddo
          endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
          dtr(1:ntx) = fqcondt(1:ntx)*tm_dum(1:ntx)
#ifdef TRDIAG_WETDEPO
          if (diag_wetdep == 1) trcond_ls(l,1:ntx) = &
               dtr(1:ntx)+thlaw(1:ntx)
#endif
          do igas=1,gases_count
            n = gases_list(igas)
            if (TM_dum(N).gt.teeny) then
              TMFAC(N)=THLAW(N)/TM_dum(N)
            else
              TMFAC(N)=0.
            end if
          enddo

          do N=1,NTX
            TRWML(N,L)  =TRWML(N,L)+ DTR(N)+THLAW(N)
            TM(L,N)     =TM(L,N)    *(1.-FQCONDT(N))   -THLAW(N)
            TMOM(:,L,N) =TMOM(:,L,N)*(1.-FQCONDT(N) - TMFAC(N))
          end do
#endif
        end if
        RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
        TH(L)=TL(L)/PLK(L)
        !     if (debug_out) write(0,*) 'after condensation: l,tlnew,',l,tl(l)
      end if

      if(RH(L).le.1.) then
        if (RH00(L).lt.1.) then
          CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        else
          CLEARA(L)=1.
        end if
      end if
      if(CLEARA(L).gt.1.) CLEARA(L)=1.
      if(RH(L).gt.1.) CLEARA(L)=0.
      if(WMX(L).le.0.) CLEARA(L)=1.
      if(CLEARA(L).lt.0.) CLEARA(L)=0.
      RHF(L)=RH00(L)+(1.-CLEARA(L))*(1.-RH00(L))
      if(RH(L).le.RHF(L).and.RH(L).lt..999999.and.WMX(L).gt.0.) then
        !**** PRECIP OUT CLOUD WATER IF RH LESS THAN THE RH OF THE ENVIRONMENT
#ifdef TRACERS_WATER
        TRPRBAR(1:NTX,L) = TRPRBAR(1:NTX,L) + TRWML(1:NTX,L)
#ifdef TRDIAG_WETDEPO
        if (diag_wetdep == 1) &
             trprcp_ls(l,1:ntx)=trprcp_ls(l,1:ntx)+trwml(1:ntx,l)
#endif
        TRWML(1:NTX,L) = 0.
#endif
        PREBAR(L)=PREBAR(L)+WMX(L)*AIRM(L)*BYGRAV*BYDTsrc
        if(LHP(L).eq.LHS .and. LHX.eq.LHE) then
          HCHANG=WMX(L)*LHM
          if (debug) print*,"lsC",l,hchang
          TL(L)=TL(L)+HCHANG/(SHA*FSSL(L)+teeny)
          !         if(debug_out) write(0,*) 'after rain out: l,tlnew',l,tl(l)
          TH(L)=TL(L)/PLK(L)
        end if
        WMX(L)=0.
      end if
      prebar1(l)=prebar(l)
      !**** set phase of condensation for next box down
      PREICE(L)=0.
      if (PREBAR(L).gt.0 .and. LHP(L).eq.LHS) PREICE(L)=PREBAR(L)
      if (PREBAR(L).le.0) LHP(L)=0.
      !**** COMPUTE THE LARGE-SCALE CLOUD COVER
      if(RH(L).le.1.) then
        if (RH00(L).lt.1.) then
          CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        else
          CLEARA(L)=1.
        end if
      end if
      if(CLEARA(L).gt.1.) CLEARA(L)=1.
      if(RH(L).gt.1.) CLEARA(L)=0.
      if(WMX(L).le.teeny) WMX(L)=0.
      if(WMX(L).le.0.) CLEARA(L)=1.
      if(CLEARA(L).lt.0.) CLEARA(L)=0.
      CLDSSL(L)=FSSL(L)*(1.-CLEARA(L))
      CLDSAVL(L)=1.-CLEARA(L)
#ifdef TRACERS_WATER
      if(CLDSSL(L).gt.0.) CLOUD_YET=.true.
      if(CLOUD_YET.and.CLDSSL(L).eq.0.) BELOW_CLOUD=.true.
#endif
      TOLDUP=TOLD
      !**** ACCUMULATE SOME DIAGNOSTICS
      HCNDSS=HCNDSS+FSSL(L)*(TL(L)-TOLD)*AIRM(L)
      SSHR(L)=SSHR(L)+FSSL(L)*(TL(L)-TOLD)*AIRM(L)
      DQLSC(L)=DQLSC(L)+FSSL(L)*(QL(L)-QOLD)
    end do  ! end of loop over L

    PRCPSS=max(0d0,PREBAR(1)*GRAV*DTsrc) ! fix small round off err
#ifdef TRACERS_WATER
    do n=1,ntx
      TRPRSS(n)=TRPRBAR(n,1)
      if(t_qlimit(n)) TRPRSS(n)=max(0d0,TRPRSS(n))
    enddo
#endif
    !****
    !**** CLOUD-TOP ENTRAINMENT INSTABILITY
    !****
    do L=LP50-1,1,-1
      SM(L)=TH(L)*AIRM(L)
      QM(L)=QL(L)*AIRM(L)
      WMXM(L)=WMX(L)*AIRM(L)
      SM(L+1)=TH(L+1)*AIRM(L+1)
      QM(L+1)=QL(L+1)*AIRM(L+1)
      WMXM(L+1)=WMX(L+1)*AIRM(L+1)
      if(WMX(L+1).gt.teeny) cycle
      TOLD=TL(L)
      TOLDU=TL(L+1)
      QOLD=QL(L)
      QOLDU=QL(L+1)
      FCLD=(1.-CLEARA(L))*FSSL(L)+teeny
      if(CLEARA(L).eq.1. .or. (CLEARA(L).lt.1..and.CLEARA(L+1).lt.1.)) &
           cycle
      !       IF(WMX(L).EQ.0. .OR. (WMX(L).GT.0..AND.WMX(L+1).GT.0.)) CYCLE
      SEDGE=THBAR(TH(L+1),TH(L))
      DSE=(TH(L+1)-SEDGE)*PLK(L+1)+(SEDGE-TH(L))*PLK(L)+ &
           SLHE*(QL(L+1)-QL(L))
      DWM=QL(L+1)-QL(L)+(WMX(L+1)-WMX(L))/FCLD
      DQSDT=DQSATDT(TL(L),LHE)*QL(L)/(RH(L)+1d-30)
      BETA=(1.+BYMRAT*TL(L)*DQSDT)/(1.+SLHE*DQSDT)
      CKM=(1.+SLHE*DQSDT)*(1.+(1.-DELTX)*TL(L)/SLHE)/ &
           (2.+(1.+BYMRAT*TL(L)/SLHE)*SLHE*DQSDT)
      CKR=TL(L)/(BETA*SLHE)
      CK=DSE/(SLHE*DWM+teeny)
      SIGK=0.
      if(CKR.gt.CKM) cycle
      if(CK.gt.CKR) SIGK=COESIG*((CK-CKR)/((CKM-CKR)+teeny))**5
      EXPST=exp(-SIGK*DTsrc)
      if(L.le.1) CKIJ=EXPST
      DSEC=DWM*TL(L)/BETA
      if(CK.lt.CKR) cycle
      FPMAX=min(1d0,1.-EXPST)       ! full CTE strength
      if(FPMAX.le.0.) cycle
      if(DSE.ge.DSEC) cycle
      !**** MIXING TO REMOVE CLOUD-TOP ENTRAINMENT INSTABILITY
      if (debug) print*,"lse",l,wmxm(l:l+1),svlhxl(l:l+1)

      AIRMR=(AIRM(L+1)+AIRM(L))*BYAM(L+1)*BYAM(L)
      SMO1=SM(L)
      QMO1=QM(L)
      WMO1=WMXM(L)
      SMO2=SM(L+1)
      QMO2=QM(L+1)
      WMO2=WMXM(L+1)
      SMO12=SMO1*PLK(L)+SMO2*PLK(L+1)
      do K=1,KMAX
        UMO1(K)=UM(K,L)
        VMO1(K)=VM(K,L)
        UMO2(K)=UM(K,L+1)
        VMO2(K)=VM(K,L+1)
      end do
      FPLUME=FPMAX*FSSL(L)
      DFX=FPLUME ! not FPMAX
      do ITER=1,9
        DFX=DFX*0.5
        FMIX=FPLUME*FCLD
        FMASS=FMIX*AIRM(L)
        FMASS=min(FMASS,(AIRM(L+1)*AIRM(L))/(AIRM(L+1)+AIRM(L)))
        FMIX=FMASS*BYAM(L)
        FRAT=FMASS*BYAM(L+1)
        SMN1=SMO1*(1.-FMIX)+FRAT*SMO2
        QMN1=QMO1*(1.-FMIX)+FRAT*QMO2
        WMN1=WMO1*(1.-FMIX)+FRAT*WMO2
        SMN2=SMO2*(1.-FRAT)+FMIX*SMO1
        QMN2=QMO2*(1.-FRAT)+FMIX*QMO1
        WMN2=WMO2*(1.-FRAT)+FMIX*WMO1
        THT1=SMN1*BYAM(L)/PLK(L)
        QLT1=QMN1*BYAM(L)
        TLT1=THT1*PLK(L)
        LHX=SVLHXL(L)
        RHT1=QLT1/(QSAT(TLT1,LHX,PL(L)))
        WMT1=WMN1*BYAM(L)
        THT2=SMN2*BYAM(L+1)/PLK(L+1)
        QLT2=QMN2*BYAM(L+1)
        WMT2=WMN2*BYAM(L+1)
        SEDGE=THBAR(THT2,THT1)
        DSE=(THT2-SEDGE)*PLK(L+1)+(SEDGE-THT1)*PLK(L)+SLHE*(QLT2-QLT1)
        DWM=QLT2-QLT1+(WMT2-WMT1)/FCLD
        DQSDT=DQSATDT(TLT1,LHE)*QLT1/(RHT1+1d-30)
        BETA=(1.+BYMRAT*TLT1*DQSDT)/(1.+SLHE*DQSDT)
        CKM=(1.+SLHE*DQSDT)*(1.+(1.-DELTX)*TLT1/SLHE)/ &
             (2.+(1.+BYMRAT*TLT1/SLHE)*SLHE*DQSDT)
        DSEC=DWM*TLT1/BETA
        DSEDIF=DSE-DSEC
        if (debug) print*,"lsf",iter,dsedif,fplume,dfx,smn1,smn2
        if(DSEDIF.gt.1d-3) FPLUME=FPLUME-DFX
        if(DSEDIF.lt.-1d-3) FPLUME=FPLUME+DFX
        if(abs(DSEDIF).le.1d-3.or.FPLUME.gt.FPMAX*FSSL(L)) exit
      end do
      !        IF (FPLUME.GT.FPMAX) print*,"lsH",i_debug,j_debug,fplume,fpmax
      !     $       ,fssl(l),dfx,iter

      !**** UPDATE TEMPERATURE, SPECIFIC HUMIDITY AND MOMENTUM DUE TO CTEI
      SMN12=SMN1*PLK(L)+SMN2*PLK(L+1)
      SMN1=SMN1-(SMN12-SMO12)*AIRM(L)/((AIRM(L)+AIRM(L+1))*PLK(L))
      SMN2=SMN2-(SMN12-SMO12)*AIRM(L+1)/((AIRM(L)+AIRM(L+1))*PLK(L+1))
      TH(L)=SMN1*BYAM(L)
      TL(L)=TH(L)*PLK(L)
      !       if(debug_out) write(0,*) 'after CTEI: l,tlnew',l,tl(l)
      QL(L)=QMN1*BYAM(L)
      LHX=SVLHXL(L)
      RH(L)=QL(L)/QSAT(TL(L),LHX,PL(L))
      WMX(L)=WMN1*BYAM(L)
      TH(L+1)=SMN2*BYAM(L+1)
      QL(L+1)=QMN2*BYAM(L+1)
      WMX(L+1)=WMN2*BYAM(L+1)
      call CTMIX (SM(L),SMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
      call CTMIX (QM(L),QMOM(1,L),FMASS*AIRMR,FMIX,FRAT)
      !****
#ifdef TRACERS_ON
      do N=1,NTX
        call CTMIX (TM(L,N),TMOM(1,L,N),FMASS*AIRMR,FMIX,FRAT)
#ifdef TRACERS_WATER
        !**** mix cloud liquid water tracers as well
        TWMTMP      = TRWML(N,L  )*(1.-FMIX)+FRAT*TRWML(N,L+1)
        TRWML(N,L+1)= TRWML(N,L+1)*(1.-FRAT)+FMIX*TRWML(N,L  )
        TRWML(N,L)  = TWMTMP
#endif
      end do
#endif
      do K=1,KMAX
        UMN1(K)=(UMO1(K)*(1.-FMIX)+FRAT*UMO2(K))
        VMN1(K)=(VMO1(K)*(1.-FMIX)+FRAT*VMO2(K))
        UMN2(K)=(UMO2(K)*(1.-FRAT)+FMIX*UMO1(K))
        VMN2(K)=(VMO2(K)*(1.-FRAT)+FMIX*VMO1(K))
        UM(K,L)=UM(K,L)+(UMN1(K)-UMO1(K))*RA(K)
        VM(K,L)=VM(K,L)+(VMN1(K)-VMO1(K))*RA(K)
        UM(K,L+1)=UM(K,L+1)+(UMN2(K)-UMO2(K))*RA(K)
        VM(K,L+1)=VM(K,L+1)+(VMN2(K)-VMO2(K))*RA(K)
      end do
      !**** RE-EVAPORATION OF CLW IN THE UPPER LAYER
      QL(L+1)=QL(L+1)+WMX(L+1)/(FSSL(L)+teeny)
      if (wmxm(l+1).gt.0. .and. svlhxl(l+1).gt.0 .and. svlhxl(l+1).ne. &
           svlhxl(l)) print*,"lsT",i_debug,j_debug,wmxm(l)*svlhxl(l) &
           +wmxm(l+1)*svlhxl(l+1)-wmx(l)*airm(l)*svlhxl(l)-wmx(l+1) &
           *airm(l+1)*svlhxl(l+1),wmxm(l:l+1),svlhxl(l:l+1)
      !**** assumes that wmx(l+1) is same phase as wmx(l)?
      !**** energy fix?
      TH(L+1)=TH(L+1)-((SVLHXL(L+1)-SVLHXL(L))*BYSHA)*WMXM(L+1)/AIRM(L &
           +1)/(PLK(L+1)*FSSL(L)+teeny)

      TH(L+1)=TH(L+1)-(LHX*BYSHA)*WMX(L+1)/(PLK(L+1)*FSSL(L)+teeny)
      if (debug) print*,"lsD",l+1,LHX*WMX(L+1)
      TL(L+1)=TH(L+1)*PLK(L+1)
      !       if(debug_out) write(0,*) 'after re-evap: l,tlnew',l+1,tl(l+1)
      RH(L+1)=QL(L+1)/QSAT(TL(L+1),LHX,PL(L+1))
      WMX(L+1)=0.
#ifdef TRACERS_WATER
      TM(L+1,1:NTX)=TM(L+1,1:NTX)+TRWML(1:NTX,L+1)
#ifdef TRDIAG_WETDEPO
      if (diag_wetdep == 1) &
           trclwe_ls(l+1,1:ntx)=trclwe_ls(l+1,1:ntx)+trwml(1:ntx,l+1)
#endif
      TRWML(1:NTX,L+1)=0.
#endif
      if(RH(L).le.1.) then
        if (RH00(L).lt.1.) then
          CLEARA(L)=DSQRT((1.-RH(L))/((1.-RH00(L))+teeny))
        else
          CLEARA(L)=1.
        end if
      end if
      if(CLEARA(L).gt.1.) CLEARA(L)=1.
      if(RH(L).gt.1.) CLEARA(L)=0.
      CLDSSL(L)=FSSL(L)*(1.-CLEARA(L))
      CLDSAVL(L)=1.-CLEARA(L)
      TNEW=TL(L)
      TNEWU=TL(L+1)
      QNEW=QL(L)
      QNEWU=QL(L+1)
      HCNDSS=HCNDSS+FSSL(L)*(TNEW-TOLD)*AIRM(L)+ &
           FSSL(L+1)*(TNEWU-TOLDU)*AIRM(L+1)
      SSHR(L)=SSHR(L)+FSSL(L)*(TNEW-TOLD)*AIRM(L)
      SSHR(L+1)=SSHR(L+1)+FSSL(L+1)*(TNEWU-TOLDU)*AIRM(L+1)
      DQLSC(L)=DQLSC(L)+FSSL(L)*(QNEW-QOLD)
      DQLSC(L+1)=DQLSC(L+1)+FSSL(L+1)*(QNEWU-QOLDU)
      DCTEI(L)=DCTEI(L)+FSSL(L)*(QNEW-QOLD)*AIRM(L)*LHX*BYSHA
      DCTEI(L+1)=DCTEI(L+1)+FSSL(L+1)*(QNEWU-QOLDU)*AIRM(L+1)*LHX*BYSHA
    end do

    !**** COMPUTE CLOUD PARTICLE SIZE AND OPTICAL THICKNESS
    WMSUM=0.
#ifdef SCM
    if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
      SCM_LWP_SS = 0.d0
      SCM_IWP_SS = 0.d0
    endif
#endif
#ifdef CLD_AER_CDNC
    ACDNWS=0.
    ACDNIS=0.
    AREWS=0.
    AREIS=0.
    ALWWS=0.
    ALWIS=0.
    NLSW = 0
    NLSI = 0
    CDNC_TOMAS=0.
#endif
    do L=1,LP50
      FCLD=CLDSSL(L)+teeny
      WTEM=1.d5*WMX(L)*PL(L)/(FCLD*TL(L)*RGAS+teeny)
      LHX=SVLHXL(L)
      if(WTEM.lt.1d-10) WTEM=1.d-10
      !***Setting constant values of CDNC over land and ocean to get RCLD=f(CDNC,LWC)
      SNdO = 59.68d0/(RWCLDOX**3)
      SNdL = 174.d0
      SNdI = 0.06417127d0
      SCDNCW=SNdO*(1.-PEARTH)+SNdL*PEARTH
      SCDNCI=SNdI
#if (defined CLD_AER_CDNC ) || (defined BLK_2MOM)
#ifdef ALT_CDNC_INPUTS
      VVEL = VVEL_sv(l)    ! retrieve value at this level
      DSU(:) = DSU_SV(:,L) ! retrieve
      NEWCLD = 1.-CLEARA(L)
      SAVCLD = CLDSAV0(L) ! from prev. timestep
#else
      ! These choices always produce DCLD <= 0.
      NEWCLD = CLDSSL(L)  ! = (1-CLEARA)*FSSL (updated)
      SAVCLD = CLDSAVL(L) ! = (1-CLEARA)      (already updated)
#endif
!@auth Menon for CDNC prediction
#ifdef TRACERS_AEROSOLS_Koch
      call GET_CDNC_UPD(L,LHX,WCONST,WMUI,WMX(L),FCLD,NEWCLD, &
           SAVCLD,VVEL,SME(L),DSU,OLDCDL(L), &
           CDNL0,CDNL1)
      OLDCDL(L) = CDNL1
      SNd=CDNL1
      !** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
      !     if (L.eq.1)write(6,*)"BLK_2M NUPD",NEWCDN,OLDCDN
#endif
#ifdef TRACERS_AMP
      OLDCDL(L)=SNd
      OLDCDI(L)=SNdi
#endif
#ifdef TRACERS_TOMAS
       OLDCDL(L)=SNd
       OLDCDI(L)=SNdi
#endif
#endif
#if (defined CLD_AER_CDNC) || (defined BLK_2MOM)
      ! Update thermo if environment was changed
      tk0=TL(L)                 ! Temperature, [K]
      qk0=QL(L)                 ! Water vapor mixing ratio, [kq/kq]
      pk0=PL(L)                 ! Pressure, [hPa]
      w0=VVEL*1.d-02           ! Large-scale velocity, [m/s]
      v0=WTURB(L) !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
      r0= 0.0  !RDTNDL(L)               ! T tendency due to radiation, [K/s]
      ldummy=execute_bulk2m_driver('all' &
           ,tk0,qk0,pk0,w0,v0,r0)
      ! Update micro if contents were changed
      !        mdrop=WMX(L)            ! drop content, [kg water/kg air]
      !        ndrop=mdrop/mw0         ! drop concent, [No/m3]
      !        mcrys=WMXICE(L)         ! crys content, [kg water/kg air]
      !        ncrys=mcrys/mi0         ! crys concent, [No/m3]
      if(LHX.eq.LHE)  then
        mdrop =WMX(L)
        ndrop= OLDCDL(L)*1.d6  !mdrop/mw0         ! drop concent, [No/m3]
        if(WMX(L).eq.0.) ndrop=0.0
      else
        mcrys =WMX(L)
        WMXICE(L) = WMX(L)
        ncrys= OLDCDI(L)*1.d6  !mcrys/mi0         ! crystal concent, [No/m3]
        if(WMX(L).eq.0.) ncrys=0.0
      endif
      !      if(L.eq.1)write(6,*)"5th check BLK_2M",
      !    *WMX(L),OLDCDL(L),OLDCDI(L)
      !
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,'end')
      ! Get new drop & crys concentration
      !     ldummy=execute_bulk2m_driver('surabi','GET_CDNC_UPD',dtB2M,mkx)
#ifdef TRACERS_AEROSOLS_Koch
      !*** Call Lohmann's or Gultepe's scheme for CDNC
      OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
      NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
      ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx, &
           OLDCDNC,NEWCDNC)
#endif
#ifdef TRACERS_AMP
      !*** Using the AMP_actv interface from MATRIX
      ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#endif
#ifdef TRACERS_TOMAS
      !*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('toma','drop_nucl',dtB2M,mkx)
#endif
      rablk=execute_bulk2m_driver('get','value','nc') + ( &
           +execute_bulk2m_driver('get','npccn') &
           )*dtB2M
      SNd=rablk(mkx)*1.0d-6            ! ndrop, [No/m^3], SNdL, [No/cc]
      !     if(SNd.gt.20.) write(6,*)"Finally out",SNd, l
      !**** Old treatment to get CDNC for cloud changes within time steps
      DCLD(L) = NEWCLD-SAVCLD ! cloud fraction change
      !** If previous time step is clear sky
      if(SAVCLD.eq.0.) then
        SNd=SNd
      elseif (DCLD(L).le.0.d0) then
        SNd=OLDCDL(L)
      elseif(DCLD(L).gt.0.d0) then
        SNd=( (OLDCDL(L)*SAVCLD) + (SNd*DCLD(L)) )/NEWCLD
      endif
      rablk=execute_bulk2m_driver('get','value','ni') + ( &
                                !       nnuccc             ! change n due to contact droplets freez
           +execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci             ! change n due to immersion droplets freez
           +execute_bulk2m_driver('get','nnucci') &
                                !       nnuccd             ! change n freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','nnuccd') &
           )*dtB2M &
                                !      nnucmd              ! change n cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmd') &
                                !      nnucmt              ! change n cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmt')

      !      SNdI=rablk(mkx)*1.0d-6             ! from ncrys [No/m^3] to SNdI in [No/cc]
      SNdI = 0.06417127d0
      if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
      OLDCDL(L) = SNd
      OLDCDI(L) = SNdI
#ifdef TRACERS_AMP
      nactc(l,1:nmodes) =  naero(mkx,1:nmodes)
      !      do nm=1,nmodes
      !        if(nactc(l,nm).gt.0.)
      !    *   write(6,*)"NMOD1",nactc(l,nm),l,nm
      !       enddo
#endif
      !      if(L.eq.1) write(6,*)"6_LO check BLK_2M",SNd,SNdI
      ! To get effective radii in micron
      rablk=execute_bulk2m_driver('get','value','ec')  ! [micron]
#endif
#ifdef CLD_AER_CDNC
      SCDNCW=SNd
      SNdI = 0.06417127d0
      SCDNCI=SNdI
      if (SCDNCW.le.20.d0) SCDNCW=20.d0   !set min CDNC sensitivity test
      !     If (SCDNCI.le.0.06d0) SCDNCI=0.06417127d0   !set min ice crystal
      if (SCDNCI.le.0.0d0) SCDNCI=teeny           !set min ice crystal
      if(SCDNCW.gt.1400.d0) SCDNCw=1400.d0
      !     if (SCDNCW.gt.20.) write(6,*) "SCND CDNC",SCDNCW,OLDCDL(l),l
#endif

      if(LHX.eq.LHE) then

        !         RCLD=(RWCLDOX*10.*(1.-PEARTH)+7.0*PEARTH)*(WTEM*4.)**BY3
        RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCW))**BY3
        QHEATC=(QHEAT(L)+FSSL(L)*CLEARA(L)*(EC(L)+ER(L)))/LHX
        if(RCLD.gt.RWMAX.and.PREP(L).gt.QHEATC) RCLD=RWMAX
        RCLDE=RCLD/BYBR
#ifdef CLD_AER_CDNC
        !** Using the Liu and Daum paramet
        !** for spectral dispersion effects on droplet size distribution
        Repsi=1.d0 - 0.7d0*exp(-0.003d0*SCDNCW)
        Repsis=Repsi*Repsi
        Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**0.333d0)
        !     write(6,*)"RCLD",Rbeta,RCLD,SCDNCW,Repsis
        RCLDE=RCLD*Rbeta
!@auth Menon    end of addition  comment out the RCLDE definition below
#endif
#ifdef BLK_2MOM
        !        if(l.eq.1) write(6,*)"7th check BLK_2M",RCLDE,RCLD
        !    *   ,SCDNCW,WTEM
        !        rablk=execute_bulk2m_driver('get','value','ec')  ! [micron]
        !        RCLDE=rablk(mkx)
        !        if(l.eq.1) write(6,*)"8th check BLK_2M",RCLDE
#endif
      else
        !         RCLD=25.0*(WTEM/4.2d-3)**BY3 * (1.+pl(l)*xRICld)
        RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*SCDNCI))**BY3
        RCLD=min(RCLD,RIMAX)
        RCLDE=RCLD/BYBR
#ifdef BLK_2MOM
        !       if(L.eq.1)  write(6,*)"9th check BLK_2M",RCLDE
        !        rablk=execute_bulk2m_driver('get','value','ei')  ! [micron]
        !        RCLDE=rablk(mkx)
        !        if(l.eq.1) write(6,*)"10th check BLK_2M",RCLDE
#endif
      end if
      RCLDE1=5.*RCLDE          ! for precip optical thickness
      CSIZEL(L)=RCLDE
      IF(FCLD.LE.teeny.AND.CSIZEL(L).GT.25.d0) CSIZEL(L)=25.d0
#ifdef CLD_AER_CDNC  /* save for diag purposesi */
      if (FCLD.gt.1.d-5.and.LHX.eq.LHE) then
        ACDNWS(L)= SCDNCW
        AREWS(L) = RCLDE
        ALWWS(L) = WTEM
        CDN3DL(L)=SCDNCW
        CRE3DL(L)=RCLDE
        NLSW  = NLSW + 1
#ifdef TRACERS_TOMAS
        CDNC_TOMAS(L)=CDNC_NENES(L)
#endif
        !      if(ACDNWS(L).gt.20.d0) write(6,*)"INWCLD",ACDNWS(L),
        !    * SCDNCW,NLSW,AREWS(L),RCLDE,LHX
      elseif(FCLD.gt.1.d-5.and.LHX.eq.LHS) then
        ACDNIS(L)= SCDNCI
        AREIS(L) = RCLDE
        ALWIS(L) = WTEM
        CDN3DL(L)=SCDNCI
        CRE3DL(L)=RCLDE
        NLSI  = NLSI + 1
        !      if(ACDNIS(L).gt.0.d0)    write(6,*)"INICLD",ACDNIS(L),
        !    * SCDNCI,NLSI,AREIS(L),RCLDE,LHX
      end if
#endif
      TEM=AIRM(L)*WMX(L)*1.d2*BYGRAV
      TAUSSL(L)=1.5d3*TEM/(FCLD*RCLDE+teeny)
      TEM1=AIRM(L)*WMPR(L)*1.d2*BYGRAV      ! precip contribution
      TAUSSL(L)=TAUSSL(L)+1.5d3*TEM1/(FCLD*RCLDE1+teeny)
      if(FCLD.le.teeny) TAUSSL(L)=0.
      if(TAUSSL(L).gt.100.) TAUSSL(L)=100.
      if(LHX.eq.LHE) WMSUM=WMSUM+TEM      ! pick up water path
#ifdef SCM
      if (i_debug.eq.I_TARG.and.j_debug.eq.J_TARG) then
        if (LHX.eq.LHE) SCM_LWP_SS = SCM_LWP_SS + TEM
        if (LHX.eq.LHS) SCM_IWP_SS = SCM_IWP_SS + TEM
      endif
#endif
#ifdef CLD_AER_CDNC
      SMLWP=WMSUM
#endif
    end do

    !**** CALCULATE OPTICAL THICKNESS
    do L=1,LP50
      CLDSV1(L)=CLDSSL(L)
      if(WMX(L).le.0.) SVLHXL(L)=0.
      if(TAUMCL(L).eq.0..or.CKIJ.ne.1.) then
        BMAX=1.-exp(-(CLDSV1(L)/.3d0))
        if(CLDSV1(L).ge..95d0) BMAX=CLDSV1(L)
        if(L.eq.1.or.L.le.DCL) then
          CLDSSL(L)=min(CLDSSL(L)+(BMAX-CLDSSL(L))*CKIJ,FSSL(L))
          TAUSSL(L)=TAUSSL(L)*CLDSV1(L)/(CLDSSL(L)+teeny)
        end if
        if(TAUSSL(L).le.0.) CLDSSL(L)=0.
        if(L.gt.DCL .and. TAUMCL(L).le.0.) then
          CLDSSL(L)=min(CLDSSL(L)**(2.*BY3),FSSL(L))
          TAUSSL(L)=TAUSSL(L)*CLDSV1(L)**BY3
        end if
      end if
      if(TAUSSL(L).lt.0.) then
        write(6,*) 'negative TAUSS CLDSS WMX I J L=',TAUSSL(L), &
             CLDSSL(L),WMX(L),i_debug,j_debug,L
        TAUSSL(L)=0.
        CLDSSL(L)=0.
        WMX(L)=0.
      end if
    end do


#if (defined CLD_AER_CDNC) || (defined CLD_SUBDD)
    !Save variables for 3 hrly diagnostics
    !     AAA=1
    !     DO L=LP50,1,-1
    !       if (CLDSSL(L).gt.0.d0.and.CLDMCL(L).gt.0.d0) then
    !        if (TAUSSL(L).gt.0.d0.and.TAUMCL(L).gt.0.d0) then
    !         if (CLDSSL(L).LE.randu(xx))GO TO 7
    !          AERTAUSS=TAUSSL(L)
    !          AERTAU(L)=AERTAUSS
    !   7     if (CLDMCL(L).LE.randu(xx))GO TO 8
    !          AERTAUMC=TAUMCL(L)
    !          IF(AERTAUMC.GT.AERTAU(L)) AERTAU(L)=AERTAUMC
    !   8      SUMTAU=SUMTAU+AERTAU(L)
    !**
    !**   DETECT AT WHAT LAYER SUMTAU (COLUMN OPTICAL THICKNESS)
    !**   IS ABOVE THOLD (THRESHOLD SET AT .1 TAU)
    !**   THIS LAYER BECOMES CLOUD TOP LAYER (ILTOP)
    !**
    !         THOLD=(.1)
    !         IF(SUMTAU.GT.THOLD)THEN
    !         IFLAG=(AAA)
    !           IF(IFLAG.EQ.1)THEN
    !             ILTOP=L
    !             AAA=0
    !           ENDIF
    !         ENDIF
    !        ENDIF
    !       ENDIF
    !     ENDDO

    do L=1,LP50
      PRS = (PL(1)-PTOP)/SIG(1)

      if (L.ge.ls1) then
        PPRES = (SIG(L)*(PSF-PTOP)+PTOP)         !in hPa
        DPP= (SIGE(L+1)-SIGE(L))*(PSF-PTOP)      !in hPa
        TEMPR=(TL(L)/PLK(L))*(SIG(L)*(PSF-PTOP)+PTOP)**KAPA
      else
        PPRES= (SIG(L)*PRS+PTOP)
        DPP= (SIGE(L+1)-SIGE(L))*PRS
        TEMPR=(TL(L)/PLK(L))*(SIG(L)*PRS+PTOP)**KAPA
      endif

      CTEML(L)=TEMPR                                        ! Cloud temperature(K)
      D3DL(L)=DPP/PPRES*TEMPR/GRAV*(gasc*1.d03)/mair        ! For Cloud thickness (m)
      if(CLDSSL(L).gt.0.d0) CD3DL(L)=-1.d0*D3DL(L)*CLDSAVL(L)/CLDSSL(L)
      RHODK=100.d0*PPRES/(gasc*1.d03)/TEMPR*mair
      !       IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHODK*CD3DL(L) ! cld water kg m-2
      !       IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHODK*CD3DL(L) ! ice water kg m-2
      if (SVLHXL(L).eq.LHE) CL3DL(L) = WMX(L)*RHODK          ! cld water kg m-3
      if (SVLHXL(L).eq.LHS) CI3DL(L) = WMX(L)*RHODK          ! ice water kg m-3
      !      write(6,*)"CT",L,WMX(L),CD3DL(l),CL3DL(L),CI3DL(L)

      !      CTEML(L)=TL(L)                            ! Cloud temperature(K)
      !      D3DL(L)=AIRM(L)*TL(L)*bygrav*rgas/PL(L)   ! For Cloud thickness (m)
      !      IF(CLDSSL(L).GT.0.d0) CD3DL(L)= 1.d0*D3DL(L)*CLDSAVL(L)/CLDSSL(L)
      !      RHO=1d2*PL(L)/(RGAS*TL(L))
      !      IF (SVLHXL(L).EQ.LHE) CL3DL(L) = WMX(L)*RHO*CD3DL(L) ! cld water kg m-2
      !      IF (SVLHXL(L).EQ.LHS) CI3DL(L) = WMX(L)*RHO*CD3DL(L) ! ice water kg m-2
      !      write(6,*)"CT",L,WMX(L),CD3DL(l),CL3DL(L),CI3DL(L)
    end do
#endif

    return
  end subroutine LSCOND

  !**************************************************************************************

  subroutine ANVIL_OPTICAL_THICKNESS(SVLATL,                         & ! input
       RCLDX,MCDNCW,MCDNCI,RIMAX,BYBR,FCLD,TEM,WTEM,                    & ! input
       RCLD,TAUMC)                                                      ! output

!@sum  MC_OPTICAL_THICKNESS calculates anvil cloud optical thickness and droplet radius
!@auth M.S. Yao

    use CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6 &
         ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
    implicit none

!@var SVLATL latent heat
    real*8, intent(IN) :: SVLATL
!@var RCLDX RCLD multiplier
    real*8, intent(IN) :: RCLDX
!@var MCDNCW,MCDNCI,RIMAX cloud droplet for liquid and ice, max ice RCLD
    real*8, intent(IN) :: MCDNCW,MCDNCI,RIMAX
!@var BYBR factor for converting cloud particle radius to effective radius
    real*8, intent(IN) :: BYBR
!@var RCLD,TAUMC cloud droplet radius and optical thickness
    real*8, intent(OUT) :: RCLD
    real*8, intent(INOUT) :: TAUMC
!@var RCLDE effective cloud droplet radius
!@var FCLD,TEM,WTEM cloud cover and dummy variables
    real*8 FCLD,TEM,WTEM,RCLDE

    !**** CALCULATE ANVIL CLOUD OPTICAL THICKNESS AND DROPLET RADIUS
    if(SVLATL.eq.LHE)  then
      RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCW))**BY3
    else
      RCLD=RCLDX*100.d0*(WTEM/(2.d0*BY3*TWOPI*MCDNCI))**BY3
      RCLD=min(RCLD,RIMAX)
    end if
    RCLDE=RCLD/BYBR       !  effective droplet radius in anvil
    TAUMC=1.5*TEM/(FCLD*RCLDE+1.E-20)
    if(TAUMC.gt.100.) TAUMC=100.

    return
  end subroutine ANVIL_OPTICAL_THICKNESS

  !*****************************************************************************************

  subroutine MC_CLOUD_FRACTION(L,TL,PL,CCM,WCU,PLEMIN,        & ! input
       PLEL2,PLEMAX,CCMMIN,WCUMIN,LMIN,LMAX,CCMUL,CCMUL2,        & ! input
       FCLOUD)                                                   ! output

!@sum  MC_CLOUD_FRACTION calculates convective cloud fraction
!@auth M.S. Yao

    use CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6 &
         ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
    use MODEL_COM, only : dtsrc
    implicit none

!@var L index for layer
    integer, intent(IN) :: L
!@var TL,PL temperature and pressure
    real*8, intent(IN) :: TL,PL
!@var CCM,WCU convective mass and cumulus updraft speed
    real*8, intent(IN) :: CCM,WCU
!@var PLEMIN,PLEL2,PLEMAX pressure at LMIN, L+2 and LMAX+1
    real*8, intent(IN) :: PLEMIN,PLEL2,PLEMAX
!@var CCMMIN,WCUMIN CCM and WCU at LMIN
    real*8, intent(IN) :: CCMMIN,WCUMIN
!@var LMIN,LMAX convective base and maximum level
    integer, intent(IN) :: LMIN,LMAX
!@var CCMUL,CCMUL2 cloud cover multipliers
    real*8, intent(IN) :: CCMUL,CCMUL2
!@var FCLOUD cloud fraction
    real*8, intent(OUT) :: FCLOUD
!@var RHO air density
    real*8 RHO

    !**** CONVECTIVE CLOUD FRACTION IS CALCULATED BASED ON RATIO OF MASS
    !**** FLUX TO UPDRAFT SPEED; AREAL FRACTION BELOW ANVIL CLOUD BASE IS
    !**** ASSUMED TO BE TWICE THE AREA OF ACTIVE UPDRAFT CORES
    RHO=PL/(RGAS*TL)
    FCLOUD=CCMUL*CCM/(RHO*GRAV*WCU*DTsrc+teeny)

    !**** DEEP ANVIL CLOUD FRACTION IS ASSUMED TO BE 5 TIMES THE NOMINAL
    !**** CONVECTIVE CLOUD FRACTION; WILL BE OVERRIDDEN IN LSCOND IF RH
    !**** IS ABOVE THRESHOLD FOR STRATIFORM CLOUD FORMATION
    if(PLEMIN-PLEL2.ge.450.) FCLOUD=5.d0*FCLOUD

    !**** "CLOUDY" AREA BELOW CLOUD BASE INCLUDED TO REPRESENT VIRGA
    if(L.lt.LMIN) then
      FCLOUD=CCMUL*CCMMIN/(RHO*GRAV*WCUMIN*DTsrc+teeny)
    end if

    !**** CLOUD FRACTION AT DETRAINMENT LEVEL OF SHALLOW/MIDLEVEL
    !**** CONVECTION IS ASSUMED TO BE 3 TIMES NOMINAL CLOUD FRACTION
    if(PLEMIN-PLEMAX.lt.450.) then
      if(L.eq.LMAX-1) then
        FCLOUD=CCMUL2*CCM/(RHO*GRAV*WCU*DTsrc+teeny)
      end if

      !**** NO VIRGA ASSUMED BELOW SHALLOW/MIDLEVEL CONVECTION
      if(L.lt.LMIN) FCLOUD=0.
    end if
    if(FCLOUD.gt.1.) FCLOUD=1.

    return
  end subroutine MC_CLOUD_FRACTION

  !***************************************************************************************

  subroutine CONVECTIVE_MICROPHYSICS(PL,WCU,DWCU,LFRZ,WCUFRZ,TP,         &
       TI,FITMAX,PLAND,CN0,CN0I,CN0G,FLAMW,FLAMG,FLAMI,RHOIP,            &
       RHOG,ITMAX,TLMIN,TLMIN1,WMAX,                                     &
#ifdef CLD_AER_CDNC
       TL,RCLD_C,MCDNCW,CONDMU,CONDPC,                                   &
#endif
       CONDP,CONDP1,CONDIP,CONDGP)                                       ! output

!@sum  CONVECTIVE_MICROPHYSICS calculates convective microphysics and precip condensate
!@auth M.S. Yao

    use CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6 &
         ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
    implicit none

    real*8, intent(IN) :: PL,WCU
!@var PL,WCU pressure and cumulus updraft speed
    integer, intent(IN) :: LFRZ,ITMAX
!@var LFRZ,ITMAX freezing level and maximum iterations
    real*8, intent(IN) :: DWCU,WCUFRZ
!@var DWCU,WCUFRZ increment for WCU and WCU at freezing level
    real*8, intent(IN) :: TP,TI
!@var TP,TI convective plume termperature and a constant 233.16
    real*8, intent(IN) :: TLMIN,TLMIN1
!@var TLMIN,TLMIN1 temperature at LMIN and LMIN+1
    real*8, intent(IN) :: FITMAX
!@var FITMAX 1.d0/ITMAX
    real*8, intent(IN) :: PLAND,WMAX
!@var PLAND,WMAX land fraction and a constant 50.
    real*8, intent(IN) :: CN0,CN0I,CN0G
!@var CN0,CN0I,CN0G intercepts of Marshall-Palmer particle size dist.
    real*8, intent(IN) :: FLAMW,FLAMG,FLAMI
!@var FLAMW,FLAMG,FLAMI Marshall-Palmer lambda for water, graupel and ice
    real*8, intent(IN) :: RHOIP,RHOG
!@var RHOG,RHOIP density of graupel and ice particles
    real*8, intent(OUT) :: CONDP,CONDP1
!@var CONDP precipitating part of convective condensate
!@var CONDP1 precipitating and detrained parts of convective condensate
    real*8 WV,VT
!@var WV,VT convective updraft speed and precip terminal velocity
    real*8 FG,FI
!@var FG,FI fraction for graupel and ice
    real*8 DCI,DCW,DCG
!@var DCI,DCW,DCG critical cloud particle sizes for onset of precip
    real*8 TIG,DDCW
!@var TIG,DDCW dummy variables
    real*8 CONDIP,CONDGP
!@var CONDIP,CONDGP precipitating part of convective condensate
    integer ITER
!@var ITER index for iteration
#ifdef CLD_AER_CDNC
    real*8 TL,RCLD_C,MCDNCW,CONDMU,CONDPC,RHO
!@var TL temperature
!@var RHO air density
!@var CONDMU,CONDPC convective condensate
!@var RCLD_C,MCDNCW
#endif
    !****
    !**** CUMULUS MICROPHYSICS IS BASED ON DEL GENIO ET AL. (2005, J. CLIM.).  THE
    !**** PARAMETERIZATION ASSUMES A MARSHALL-PALMER PARTICLE SIZE DISTRIBUTION
    !**** AND EMPIRICAL SIZE-FALLSPEED RELATIONS FOR LIQUID, GRAUPEL AND FLUFFY ICE.
    !**** COMPARISON OF THE CUMULUS UPDRAFT SPEED TO THE FALLSPEEDS OF PARTICLES
    !**** OF DIFFERENT SIZES IS USED TO PARTITION THE PARTICLE SIZE DISTRIBUTION
    !**** INTO PARTS THAT PRECIPITATE, DETRAIN AT THE CURRENT LEVEL, OR ARE
    !**** CARRIED UPWARD TO THE NEXT LEVEL (PUBLISHED VERSION ONLY HAD A 2-WAY
    !**** SEPARATION BETWEEN PRECIPITATING AND DETRAINING PARTS).
    !****

    !**** PARTICLES WITH FALLSPEEDS WITHIN +/- DWCU OF THE UPDRAFT SPEED ARE
    !**** ASSUMED TO DETRAIN AT THE CURRENT LEVEL; FIRST FIND SMALL PARTICLE END OF
    !**** THIS RANGE (LOWER CRITICAL SIZE) FOR GRAUPEL AND FLUFFY ICE

    WV=WCU-DWCU          ! apply the calculated cumulus updraft speed
    if(WV.lt.0.d0) WV=0.d0

    !**** LOWER CRITICAL SIZE FOR GRAUPEL (EQ. 4B OF DEL GENIO ET AL. (2005))
    DCG=((WV/19.3)*(PL/1000.)**.4)**2.7
    DCG=min(DCG,1.D-2)

    !**** LOWER CRITICAL SIZE FOR FLUFFY ICE (MODIFIED FROM EQ. 4C OF DEL GENIO ET AL.
    !**** (2005) TO REPRESENT UNRIMED RADIATING ASSEMBLAGES RATHER THAN DENDRITES)
    DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
    DCI=min(DCI,1.D-2)

    !**** ABOVE MELTING LEVEL CONDENSATE IS ASSUMED TO BE FROZEN, TRANSITIONING FROM
    !**** PARTLY GRAUPEL TO PURE FLUFFY ICE OVER A RANGE OF TEMPERATURES THAT INCREASES
    !**** AS MELTING LEVEL UPDRAFT SPEED STRENGTHENS
    if(LFRZ.eq.0) then
      TIG=TF           ! initialization
    else
      TIG=TF-4.*WCUFRZ
    end if
    if(TIG.lt.TI-10.d0) TIG=TI-10.d0

    !**** LOWER CRITICAL SIZE FOR LIQUID WATER (EQ. 4A OF DEL GENIO ET AL. (2005))
    if (TP.ge.TF) then ! water phase
      DDCW=6d-3*FITMAX
      if(PLAND.lt..5) DDCW=1.5d-3*FITMAX
      DCW=0.
      do ITER=1,ITMAX-1
        VT=(-.267d0+DCW*(5.15D3-DCW*(1.0225D6-7.55D7*DCW)))* &
             (1000./PL)**.4d0
        if(VT.ge.0..and.abs(VT-WV).lt..3) exit
        if(VT.gt.WMAX) exit
        DCW=DCW+DDCW
      end do

      !**** PART OF LIQUID CONDENSATE NOT ADVECTED UP (EQ. 5, DEL GENIO ET AL. (2005))
      CONDP1 = PRECIP_MP (RHOW,FLAMW,DCW,CN0)

      !**** UPPER CRITICAL SIZE FOR LIQUID WATER
      WV=WCU+DWCU          ! apply the calculated cumulus updraft speed
      DCW=0.
      do ITER=1,ITMAX-1
        VT=(-.267d0+DCW*(5.15D3-DCW*(1.0225D6-7.55D7*DCW)))* &
             (1000./PL)**.4d0
        if(VT.ge.0..and.abs(VT-WV).lt..3) exit
        if(VT.gt.WMAX) exit
        DCW=DCW+DDCW
      end do

      !**** PRECIPITATING PART OF LIQUID CONDENSATE (RAIN)
      CONDP = PRECIP_MP (RHOW,FLAMW,DCW,CN0)

#ifdef CLD_AER_CDNC
      RHO=1d2*PL/(RGAS*TL) !air density in kg/m3
      !** Here we calculate condensate converted to precip
      !** only if drop size (Reff) exceeds 14 um
      !** Conver Rvol to m and CDNC to m from um and cm, respectively
      !** Precip condensate is simply the LWC
      CONDPC=4.d0*by3*PI*RHOW*((RCLD_C*1.d-6)**3)*1.d6*MCDNCW/RHO
      if(CONDMU.gt.CONDPC)  then
        CONDP=CONDMU-CONDPC !
      else
        CONDP=0.d0
      end if
      !     if (CONDP(L).lt.0.d0)
      !    *write(6,*)"Mup",CONDP(L),CONDPC(L),CONDMU,DCW,MCDNCW,RCLD_C,L
      !      write(6,*)"Mup",CONDP(L),DCW,FLAMW,CONDMU,L
#endif

      !**** ABOVE MIXED PHASE REGION, PART OF FLUFFY ICE NOT ADVECTED UP
    else if (TP.le.TIG) then ! pure ice phase, use TIG
      CONDP1 = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)

      !**** NOW FIND LARGE PARTICLE END OF RANGE FOR DETRAINMENT AT CURRENT
      !**** LEVEL FOR EACH FROZEN PARTICLE TYPE (UPPER CRITICAL SIZE)
      WV=WCU+DWCU

      !**** UPPER CRITICAL SIZE FOR FLUFFY ICE
      DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
      DCI=min(DCI,1.D-2)

      !**** PRECIPITATING PART OF FLUFFY ICE (SNOW)
      CONDP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)

      !**** WITHIN MIXED PHASE REGION, CALCULATE FRACTION OF FROZEN CONDENSATE
      !**** THAT IS GRAUPEL VS. FLUFFY ICE
    else ! mixed phase
      FG=(TP-TIG)/((TF-TIG)+teeny)
      if(FG.gt.1d0) FG=1d0
      if(FG.lt.0d0) FG=0d0
      if(TLMIN.le.TF.or.TLMIN1.le.TF) FG=0.d0
      FI=1.-FG

      !**** PART OF FROZEN CONDENSATE NOT ADVECTED UP IN MIXED PHASE REGION
      CONDIP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)
      CONDGP = PRECIP_MP (RHOG, FLAMG,DCG,CN0G)
      CONDP1=FG*CONDGP+FI*CONDIP
      WV=WCU+DWCU

      !**** UPPER CRITICAL SIZE FOR FLUFFY ICE
      DCI=((WV/11.72)*(PL/1000.)**.4)**2.439
      DCI=min(DCI,1.D-2)

      !**** UPPER CRITICAL SIZE FOR GRAUPEL
      DCG=((WV/19.3)*(PL/1000.)**.4)**2.7
      DCG=min(DCG,1.D-2)

      !**** PRECIPITATING PART OF FROZEN CONDENSATE (HAIL + SNOW) IN MIXED PHASE REGION
      CONDIP = PRECIP_MP (RHOIP,FLAMI,DCI,CN0I)
      CONDGP = PRECIP_MP (RHOG, FLAMG,DCG,CN0G)
      CONDP=FG*CONDGP+FI*CONDIP
    end if

    return
  end subroutine CONVECTIVE_MICROPHYSICS

  !**************************************************************************************

  subroutine MC_PRECIP_PHASE(L,AIRM,FEVAP,LHP1,                     & ! input
       PRCP,TOLD,TOLD1,VLAT,COND,LMIN,                                 & ! input
       LHP,MCLOUD,HEAT1)                                               ! output

!@sum  MC_PRECIP_PHASE calculates the convective precip phase
!@auth M.S. Yao

    use CONSTANT, only :rgas,grav,lhe,lhs,lhm,sha,bysha,pi,by6 &
         ,by3,tf,bytf,rvap,bygrav,deltx,bymrat,teeny,gamd,rhow,twopi
    implicit none

!@var L,LMIN index for layer and the convective base layer
    integer, intent(IN) :: L,LMIN
!@var AIRM air mass
    real*8, intent(IN) :: AIRM
!@var FEVAP,LHP1 fraction for evaporation, precip phase at L+1
    real*8, intent(IN) :: FEVAP,LHP1
!@var TOLD,TOLD1 old temperature at L and L+1
    real*8, intent(IN) :: TOLD,TOLD1
!@var VLAT,COND save LHX temporarily, condensate
    real*8, intent(IN) :: VLAT,COND
!@var LHP precip phase
    real*8, intent(OUT) :: LHP
!@var MCLOUD air mass for re-evaporation of precip
    real*8, intent(OUT) :: MCLOUD
!@var HEAT1 heating for phase change
    real*8, intent(OUT) :: HEAT1
!@var PRCP precipitation
    real*8, intent(INOUT) :: PRCP

    !**** DETERMINE THE AIR MASS FOR RAIN EVAPORATION TO TAKE PLACE
    MCLOUD=0.
    if(L.le.LMIN) MCLOUD=2.*FEVAP*AIRM
    if(MCLOUD.gt.AIRM) MCLOUD=AIRM
    lhp=lhp1

    !**** phase of precip is defined by original environment (i.e. TOLD)
    !**** phase of condensate is defined by VLAT (based on plume T).

    !**** decide whether to melt frozen precip
    if (lhp.eq.lhs .and. TOLD.gt.TF.and.TOLD1.le.TF) then
      if (debug) print*,"cnv1",l,lhp,told,told1,LHM*PRCP*BYSHA
      HEAT1=HEAT1+LHM*PRCP*BYSHA
      lhp=lhe
    end if
    !**** and deal with possible inversions and re-freezing of rain
    if (lhp.eq.lhe .and. TOLD.le.TF.and.TOLD1.gt.TF) then
      if (debug) print*,"cnv2",l,lhp,told,told1,-LHM*PRCP*BYSHA
      HEAT1=HEAT1-LHM*PRCP*BYSHA
      lhp=lhs
    end if
    !**** check for possible inconsistency with cond at this level
    if (LHP.ne.VLAT.and. COND.gt.0) then ! convert to cond to precip phase
      HEAT1=HEAT1+(VLAT-LHP)*COND*BYSHA
      if (debug) print*,"cnv3",l,lhp,vlat,cond,(VLAT-LHP &
           )*COND*BYSHA
    end if

    return
  end subroutine MC_PRECIP_PHASE

  !***************************************************************************************

  function  PRECIP_MP (RHO, FLAM, DC, CN)
    !
    !     Calculates the mass density of the precipitating condensate part of the
    !     Marshall-Palmer distribution for specified phase (liquid,graupel,ice)
    !
    !     INPUTS:
    !        RHO  -  Density pf Condenssed Phase
    !        FLAM -  Marshall-Palmer Slope Parameter
    !        DC   -  Critical Size for Precipitation
    !        CN   -  Marshall-Palmer Intercept
    !
    use CONSTANT, only : PI, BY6
    implicit none
    !
    real*8   PRECIP_MP, RHO, FLAM, DC, CN
    !
    PRECIP_MP = RHO*(PI*BY6)*CN*exp(-FLAM*DC)* &
         (DC*DC*DC/FLAM+3.*DC*DC/(FLAM*FLAM)+ &
         6.*DC/(FLAM*FLAM*FLAM)+6./FLAM**4.)
    !
    return
  end function PRECIP_MP

  !***************************************************************************************

  subroutine MASS_FLUX ( FPLUME, FMP2, DQSUM, LMIN, &
       LHX, QMO1, QMO2, SLH, SMO1, SMO2, WMDN, WMUP, WMEDG )
    !
    use CONSTANT, only : deltx
    !         deltx = coeff. of humidity in virtual temperature defn. (0.6078)
    !
    implicit none
    !
    ! ****************************************************************************
    ! * MASS_FLUX  - PERFORMS ITERATION TO FIND FPLUME MASS WHICH RESTORES CLOUD *
    ! *              BASE TO NEUTRAL STATE                                       *
    ! ****************************************************************************
    !
    ! **************
    ! **  INPUTS  **
    ! **************
    !        REAL*8
    !             AIRM                        the layer's pressure depth (mb)
    !             BYAM                        1./AIRM
    !             PL                          layer pressure (mb)
    !             PLK                         PL**KAPA
    !             SM,QM                       Vertical profiles of T/Q
    integer &
         LMIN                     !  base layer of a convective event
    real*8 &
         LHX,                     & !  latent heat of evaporation (J/Kg)
         QMO1,QMO2, SMO1,SMO2,    & !  LMIN & LMIN+1 elements of QM & SM
         SLH,                     & !  LHX/SHA
         WMDN,WMUP,WMEDG          !  LMIN & LMIN+1 elements of WML & middle
    !
    ! ***************
    ! **  OUTPUTS  **
    ! ***************
    real*8 &
         FPLUME,                  & !  fraction of convective plume
         FMP2,                    & !  FPLUME*AIRM(LMIN)
         DQSUM                    !  SUM OF DQs DURING ITERATION
    !
    ! **************
    ! **  LOCAL   **
    ! **************
    integer ITER        !   number for iteration
    real*8 &
         DMSE1,         & !   difference in moist static energy
         DFP,           & !   an iterative increment
         FEVAP,         & !   fraction of air available for precip evaporation
         MCLOUD,        & !   air mass available for re-evaporation of precip
         SMP, QMP,      & !   plume's SM, QM
         TP,            & !   plume's temperature (K)
         QSATC,         & !   saturation vapor mixing ratio
                                !             QNX,TNX        !   Dummy Variables in /CLDPRV/
         FRAT1,FRAT2,QMN1,QMN2,SMN1,SMN2,   & !   Dummy Variables
         GAMA,DQ,QDN,QUP,SDN,SUP,SEDGE,     & !   Dummy Variables
         SVDN,SVUP,SVEDG,QEDGE              !   Dummy Variables
    !
    ! *************
    !   Functions
    !     QSATMP  plume's saturation vapor mixing ratio (function)
    !     QSAT    saturation specific humidity
    !     DQSATDT dQSAT/dT
    !     THBAR   virtual temperature at layer edge
    real*8 QSATMP, QSAT, DQSATDT, THBAR
    !
    !
    FPLUME=.25
    DFP = .25
    do ITER=1,8
      DFP=DFP*0.5
      FMP2=FPLUME*AIRM(LMIN)
      FRAT1=FMP2*BYAM(LMIN+1)
      FRAT2=FMP2*BYAM(LMIN+2)
      SMN1=SMO1*(1.-FPLUME)+FRAT1*SMO2
      QMN1=QMO1*(1.-FPLUME)+FRAT1*QMO2
      SMN2=SMO2*(1.-FRAT1)+FRAT2*SM(LMIN+2)
      QMN2=QMO2*(1.-FRAT1)+FRAT2*QM(LMIN+2)
      SMP=SMO1*FPLUME
      QMP=QMO1*FPLUME
      TP=SMO1*PLK(LMIN+1)*BYAM(LMIN)
      QSATMP=FMP2*QSAT(TP,LHX,PL(LMIN+1))
      GAMA=SLH*QSATMP*DQSATDT(TP,LHX)/FMP2
      DQSUM=(QMP-QSATMP)/(1.+GAMA)
      if(DQSUM.gt.0.)  then
        FEVAP=.5*FPLUME
        MCLOUD=FEVAP*AIRM(LMIN+1)
        TNX=SMO2*PLK(LMIN+1)*BYAM(LMIN+1)
        QNX=QMO2*BYAM(LMIN+1)
        QSATC=QSAT(TNX,LHX,PL(LMIN+1))
        DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
        if(DQ.gt.DQSUM) DQ=DQSUM
        SMN2=SMN2-SLH*DQ/PLK(LMIN+1)
        QMN2=QMN2+DQ
        DQSUM=DQSUM-DQ
        if(DQSUM.gt.0.)  then
          MCLOUD=FEVAP*AIRM(LMIN)
          TNX=SMO1*PLK(LMIN)*BYAM(LMIN)
          QNX=QMO1*BYAM(LMIN)
          QSATC=QSAT(TNX,LHX,PL(LMIN))
          DQ=MCLOUD*(QSATC-QNX)/(1.+SLH*QSATC*DQSATDT(TNX,LHX))
          if(DQ.gt.DQSUM) DQ=DQSUM
          SMN1=SMN1-SLH*DQ/PLK(LMIN)
          QMN1=QMN1+DQ
        end if
      end if
      SDN=SMN1*BYAM(LMIN)
      SUP=SMN2*BYAM(LMIN+1)
      SEDGE=THBAR(SUP,SDN)
      QDN=QMN1*BYAM(LMIN)
      QUP=QMN2*BYAM(LMIN+1)
      SVDN=SDN*(1.+DELTX*QDN-WMDN)
      SVUP=SUP*(1.+DELTX*QUP-WMUP)
      QEDGE=.5*(QUP+QDN)
      SVEDG=SEDGE*(1.+DELTX*QEDGE-WMEDG)
      DMSE1=(SVUP-SVEDG)*PLK(LMIN+1)+(SVEDG-SVDN)*PLK(LMIN)+ &
           SLHE*(QSAT(SUP*PLK(LMIN+1),LHX,PL(LMIN+1))-QDN)
      if (abs(DMSE1).le.1.d-3) exit ! finish iteration
      if(DMSE1.gt. 1.d-3) FPLUME=FPLUME-DFP
      if(DMSE1.lt.-1.d-3) FPLUME=FPLUME+DFP
    end do
    !
    return
    !
  end subroutine MASS_FLUX

  !****************************************************************************************

end module CLOUDS

!----------

subroutine ISCCP_CLOUD_TYPES(sunlit,pfull &
     ,phalf,qv,cc,conv,dtau_s,dtau_c,skt,at,dem_s,dem_c,itrop &
     ,fq_isccp,meanptop,meantaucld,boxtau,boxptop,nbox,jerr)
!@sum  ISCCP_CLOUD_TYPES calculate isccp cloud diagnostics in a column
!@auth Gavin Schmidt
!@ver  2.0 (from isccp version 3.5)
!@     GISS modifications:
!@       1) used modelE random number generator
!@       2) real-->real*8, e-2 -> d-2 etc.
!@       3) ncol, nlev taken from modules, add jerr flag
!@       4) emsfc_lw,top_height,overlap,npoints fixed
!@       5) sub in for constants
!@       6) use pre-calculated tropopause
!@       7) tautab/invtau from module
!@       8) removed boxtau,boxptop from output
!@       9) added back nbox for backwards compatibility
!@@@@@ added boxtau,boxptop back to output      Audrey Wolf
  ! *****************************COPYRIGHT*******************************
  ! (c) COPYRIGHT Steve Klein and Mark Webb 2004, All Rights Reserved.
  ! Steve Klein klein21@mail.llnl.gov
  ! Mark Webb mark.webb@metoffice.gov.uk markwebb@mail.com
  ! ISCCP SIMULATOR icarus-scops version 3.5
  !
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version.
  !
  ! This library is distributed in the hope that it will be useful,
  ! but WITHOUT ANY WARRANTY; without even the implied warranty of
  ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  ! Lesser General Public License for more details.
  !
  ! You should have received a copy of the GNU Lesser General Public
  ! License along with this library; if not, write to the Free Software
  ! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  ! See also http://www.gnu.org/copyleft/lesser.txt
  !
  ! The Met Office hereby disclaims all copyright interest in
  ! the ISCCP Simulator ( a library which converts model clouds into
  ! direct equivalents to the satellite observations of the ISCCP)
  ! written by Steve Klein and Mark Webb
  !
  ! Catherine Senior, August 2004
  ! Manager Climate Sensitivy Group
  ! Met Office   Hadley Centre for Climate Prediction and Research
  ! FitzRoy Road  Exeter  EX1 3PB  United Kingdom
  ! *****************************COPYRIGHT*******************************
  use CONSTANT, only : wtmair=>mair,Navo=>avog,bygrav,bymrat
  use RANDOM, only : randu
  use RESOLUTION, only : nlev=>lm
  use MODEL_COM, only : qcheck
  use CLOUDS, only : tautab,invtau
  use CLOUDS_COM, only : ncol
  implicit none
!@var  emsfc_lw    longwave emissivity of surface at 10.5 microns
  real*8, parameter :: emsfc_lw=0.99d0
!@var pc1bylam Planck constant c1 by wavelength (10.5 microns)
  real*8, parameter :: pc1bylam = 1.439d0/10.5d-4
!@var boxarea fractional area of each sub-grid scale box
  real*8, parameter :: boxarea = 1d0/ncol
!@var t0 ave temp  (K)
  real*8, parameter :: t0 = 296.
!@var t0bypstd ave temp by sea level press
  real*8, parameter :: t0bypstd = t0/1.013250d6
  real*8, parameter :: bywc = 1./2.56d0 , byic= 1./2.13d0
  real*8, parameter :: isccp_taumin = 0.3d0  ! used to be 0.1

  !  number of model points in the horizontal (FIXED FOR ModelE CODE)
  integer, parameter :: npoints=1

  !     -----
  !     Input
  !     -----

  !      INTEGER nlev          !  number of model levels in column
  !      INTEGER ncol          !  number of subcolumns

  integer sunlit(npoints) !  1 for day points, 0 for night time

  !      INTEGER seed(npoints)
  !  seed values for marsaglia  random number generator
  !  It is recommended that the seed is set
  !  to a different value for each model
  !  gridbox it is called on, as it is
  !  possible that the choice of the same
  !  seed value every time may introduce some
  !  statistical bias in the results, particularly
  !  for low values of NCOL.

  real*8 pfull(npoints,nlev)
  !  pressure of full model levels (Pascals)
  !  pfull(npoints,1) is top level of model
  !  pfull(npoints,nlev) is bot of model

  real*8 phalf(npoints,nlev+1)
  !  pressure of half model levels (Pascals)
  !  phalf(npoints,1) is top of model
  !  phalf(npoints,nlev+1) is the surface pressure

  real*8 qv(npoints,nlev)
  !  water vapor specific humidity (kg vapor/ kg air)
  !         on full model levels

  real*8 cc(npoints,nlev)
  !  input cloud cover in each model level (fraction)
  !  NOTE:  This is the HORIZONTAL area of each
  !         grid box covered by clouds

  real*8 conv(npoints,nlev)
  !  input convective cloud cover in each model
  !   level (fraction)
  !  NOTE:  This is the HORIZONTAL area of each
  !         grid box covered by convective clouds

  real*8 dtau_s(npoints,nlev)
  !  mean 0.67 micron optical depth of stratiform
  !  clouds in each model level
  !  NOTE:  this the cloud optical depth of only the
  !  cloudy part of the grid box, it is not weighted
  !  with the 0 cloud optical depth of the clear
  !         part of the grid box

  real*8 dtau_c(npoints,nlev)
  !  mean 0.67 micron optical depth of convective
  !  clouds in each
  !  model level.  Same note applies as in dtau_s.

  !      INTEGER overlap           !  overlap type
  !  1=max
  !  2=rand
  !  3=max/rand

  !      INTEGER top_height        !  1 = adjust top height using both a computed
  !  infrared brightness temperature and the visible
  !  optical depth to adjust cloud top pressure. Note
  !  that this calculation is most appropriate to compare
  !  to ISCCP data during sunlit hours.
  !  2 = do not adjust top height, that is cloud top
  !  pressure is the actual cloud top pressure
  !  in the model
  !  3 = adjust top height using only the computed
  !  infrared brightness temperature. Note that this
  !  calculation is most appropriate to compare to ISCCP
  !  IR only algortihm (i.e. you can compare to nighttime
  !  ISCCP data with this option)
  integer, parameter :: top_height=1, overlap=3


  !      REAL*8 tautab(0:255)      !  ISCCP table for converting count value to
  !  optical thickness

  !      INTEGER invtau(-20:45000) !  ISCCP table for converting optical thickness
  !  to count value
  !
  !     The following input variables are used only if top_height = 1 or top_height = 3
  !
  real*8 skt(npoints)       !  skin Temperature (K)
  !      REAL*8 emsfc_lw           !  10.5 micron emissivity of surface (fraction)

  real*8 at(npoints,nlev)   !  temperature in each model level (K)
  real*8 dem_s(npoints,nlev) !  10.5 micron longwave emissivity of stratiform
  !  clouds in each
  !  model level.  Same note applies as in dtau_s.
  real*8 dem_c(npoints,nlev) !  10.5 micron longwave emissivity of convective
  !  clouds in each
  !  model level.  Same note applies as in dtau_s.
  integer itrop(npoints)    ! model level containing tropopause (WMO defn)

  !     ------
  !     Output
  !     ------

  real*8 fq_isccp(npoints,7,7) !  the fraction of the model grid box covered by
  !  each of the 49 ISCCP D level cloud types

  real*8 totalcldarea(npoints) !  the fraction of model grid box columns
  !  with cloud somewhere in them.  This should
  !  equal the sum over all entries of fq_isccp
  integer nbox(npoints) !  the number of model grid box columns
  !  with cloud somewhere in them.

  ! The following three means are averages over the cloudy areas only.  If no
  ! clouds are in grid box all three quantities should equal zero.

  real*8 meanptop(npoints)  !  mean cloud top pressure (mb) - linear averaging
  !  in cloud top pressure.

  real*8 meantaucld(npoints) !  mean optical thickness
  !  linear averaging in albedo performed.

  real*8 boxtau(npoints,ncol) !  optical thickness in each column

  real*8 boxptop(npoints,ncol) !  cloud top pressure (mb) in each column

  integer JERR  ! error flag

  !
  !     ------
  !     Working variables added when program updated to mimic Mark Webb's PV-Wave code
  !     ------

  real*8 frac_out(npoints,ncol,nlev) ! boxes gridbox divided up into
  ! Equivalent of BOX in original version, but
  ! indexed by column then row, rather than
  ! by row then column

  real*8 tca(npoints,0:nlev) ! total cloud cover in each model level (fraction)
  ! with extra layer of zeroes on top
  ! in this version this just contains the values input
  ! from cc but with an extra level
  real*8 cca(npoints,nlev)  ! convective cloud cover in each model level (fraction)
  ! from conv

  real*8 threshold(npoints,ncol) ! pointer to position in gridbox
  real*8 maxocc(npoints,ncol) ! Flag for max overlapped conv cld
  real*8 maxosc(npoints,ncol) ! Flag for max overlapped strat cld

  real*8 boxpos(npoints,ncol) ! ordered pointer to position in gridbox

  real*8 threshold_min(npoints,ncol) ! minimum value to define range in with new threshold
  ! is chosen

  real*8 dem(npoints,ncol),bb(npoints) !  working variables for 10.5 micron longwave
  !  emissivity in part of
  !  gridbox under consideration

  real*8 ran(npoints)       ! vector of random numbers
  real*8 ptrop(npoints)
  real*8 attrop(npoints)
  !      REAL*8 attropmin (npoints)
  real*8 atmax(npoints)
  real*8 atmin(npoints)
  real*8 btcmin(npoints)
  real*8 transmax(npoints)

  integer i,j,ilev,ibox
  integer ipres(npoints)
  integer itau(npoints),ilev2
  integer acc(nlev,ncol)
  integer match(npoints,nlev-1)
  integer nmatch(npoints)
  integer levmatch(npoints,ncol)

  !variables needed for water vapor continuum absorption
  real*8 fluxtop_clrsky(npoints),trans_layers_above_clrsky(npoints)
  real*8 taumin(npoints)
  real*8 dem_wv(npoints,nlev)
  real*8 press(npoints), dpress(npoints), atmden(npoints)
  real*8 rvh20(npoints), wk(npoints), rhoave(npoints)
  real*8 rh20s(npoints), rfrgn(npoints)
  real*8 tmpexp(npoints),tauwv(npoints)

  character*1 cchar(6),cchar_realtops(6)
  integer icycle
  real*8 tau(npoints,ncol)
  logical box_cloudy(npoints,ncol)
  real*8 tb(npoints,ncol)
  real*8 ptop(npoints,ncol)
  real*8 emcld(npoints,ncol)
  real*8 fluxtop(npoints,ncol)
  real*8 trans_layers_above(npoints,ncol)
  real*8 fluxtopinit(npoints),tauir(npoints)
  real*8 meanalbedocld(npoints)
  real*8 albedocld(npoints,ncol)
  !      integer debug       ! set to non-zero value to print out inputs
  !                    ! with step debug
  !      integer debugcol    ! set to non-zero value to print out column
  !                    ! decomposition with step debugcol
  integer rangevec(npoints),rangeerror

  real*8 tauchk,xx

  !      character*10 ftn09

  data cchar / ' ','-','1','+','I','+'/
  data cchar_realtops / ' ',' ','1','1','I','I'/

  JERR=0

  !      INTEGER irand,i2_16,huge32,overflow_32  ! variables for RNG
  !      PARAMETER(huge32=2147483647)
  !      i2_16=65536

  tauchk = -1.*log(0.9999999)

  !      ncolprint=0
  !
  !      if ( debug.ne.0 ) then
  !          j=1
  !          write(6,'(a10)') 'j='
  !          write(6,'(8I10)') j
  !          write(6,'(a10)') 'debug='
  !          write(6,'(8I10)') debug
  !          write(6,'(a10)') 'debugcol='
  !          write(6,'(8I10)') debugcol
  !          write(6,'(a10)') 'npoints='
  !          write(6,'(8I10)') npoints
  !          write(6,'(a10)') 'nlev='
  !          write(6,'(8I10)') nlev
  !          write(6,'(a10)') 'ncol='
  !          write(6,'(8I10)') ncol
  !          write(6,'(a10)') 'top_height='
  !          write(6,'(8I10)') top_height
  !          write(6,'(a10)') 'overlap='
  !          write(6,'(8I10)') overlap
  !          write(6,'(a10)') 'emsfc_lw='
  !          write(6,'(8f10.2)') emsfc_lw
  !          write(6,'(a10)') 'tautab='
  !          write(6,'(8f10.2)') tautab
  !          write(6,'(a10)') 'invtau(1:100)='
  !          write(6,'(8i10)') (invtau(i),i=1,100)
  !        do j=1,npoints,debug
  !          write(6,'(a10)') 'j='
  !          write(6,'(8I10)') j
  !          write(6,'(a10)') 'sunlit='
  !          write(6,'(8I10)') sunlit(j)
  !          write(6,'(a10)') 'seed='
  !          write(6,'(8I10)') seed(j)
  !          write(6,'(a10)') 'pfull='
  !          write(6,'(8f10.2)') (pfull(j,i),i=1,nlev)
  !          write(6,'(a10)') 'phalf='
  !          write(6,'(8f10.2)') (phalf(j,i),i=1,nlev+1)
  !          write(6,'(a10)') 'qv='
  !          write(6,'(8f10.3)') (qv(j,i),i=1,nlev)
  !          write(6,'(a10)') 'cc='
  !          write(6,'(8f10.3)') (cc(j,i),i=1,nlev)
  !          write(6,'(a10)') 'conv='
  !          write(6,'(8f10.2)') (conv(j,i),i=1,nlev)
  !          write(6,'(a10)') 'dtau_s='
  !          write(6,'(8g12.5)') (dtau_s(j,i),i=1,nlev)
  !          write(6,'(a10)') 'dtau_c='
  !          write(6,'(8f10.2)') (dtau_c(j,i),i=1,nlev)
  !          write(6,'(a10)') 'skt='
  !          write(6,'(8f10.2)') skt(j)
  !          write(6,'(a10)') 'at='
  !          write(6,'(8f10.2)') (at(j,i),i=1,nlev)
  !          write(6,'(a10)') 'dem_s='
  !          write(6,'(8f10.3)') (dem_s(j,i),i=1,nlev)
  !          write(6,'(a10)') 'dem_c='
  !          write(6,'(8f10.3)') (dem_c(j,i),i=1,nlev)
  !        enddo
  !      endif

  !     ---------------------------------------------------!

  !     assign 2d tca array using 1d input array cc

  do j=1,npoints
    tca(j,0)=0
  enddo

  do ilev=1,nlev
    do j=1,npoints
      tca(j,ilev)=cc(j,ilev)
    enddo
  enddo

  !     assign 2d cca array using 1d input array conv

  do ilev=1,nlev
    do j=1,npoints
      cca(j,ilev)=conv(j,ilev)
    enddo
  enddo

  !      if (ncolprint.ne.0) then
  !      do j=1,npoints,1000
  !        write(6,'(a10)') 'j='
  !        write(6,'(8I10)') j
  !        write (6,'(a)') 'seed:'
  !        write (6,'(I3.2)') seed(j)
  !
  !        write (6,'(a)') 'tca_pp_rev:'
  !        write (6,'(8f5.2)')
  !     &   ((tca(j,ilev)),
  !     &      ilev=1,nlev)
  !
  !        write (6,'(a)') 'cca_pp_rev:'
  !        write (6,'(8f5.2)')
  !     &   ((cca(j,ilev),ibox=1,ncolprint),ilev=1,nlev)
  !      enddo
  !      endif

  if (top_height .eq. 1 .or. top_height .eq. 3) then

    do j=1,npoints ! use pre-computed tropopause level
      ptrop(j) = pfull(j,itrop(j))
      attrop(j) = at(j,itrop(j))
      atmin(j) = 400.
      atmax(j) = 0.
    enddo

    do ilev=1,nlev
      do j=1,npoints
        if (at(j,ilev) .gt. atmax(j)) atmax(j)=at(j,ilev)
        if (at(j,ilev) .lt. atmin(j)) atmin(j)=at(j,ilev)
      enddo
    end do

  end if

    !     -----------------------------------------------------!

    !     ---------------------------------------------------!

  do ilev=1,nlev
    do j=1,npoints

      rangevec(j)=0

      if (cc(j,ilev) .lt. 0. .or. cc(j,ilev) .gt. 1.) then
        !           error = cloud fraction less than zero
        !           error = cloud fraction greater than 1
        rangevec(j)=rangevec(j)+1
      endif

      if (conv(j,ilev) .lt. 0. .or. conv(j,ilev) .gt. 1.) then
        !           ' error = convective cloud fraction less than zero'
        !           ' error = convective cloud fraction greater than 1'
        rangevec(j)=rangevec(j)+2
      endif

      if (dtau_s(j,ilev) .lt. 0.) then
        !           ' error = stratiform cloud opt. depth less than zero'
        rangevec(j)=rangevec(j)+4
      endif

      if (dtau_c(j,ilev) .lt. 0.) then
        !           ' error = convective cloud opt. depth less than zero'
        rangevec(j)=rangevec(j)+8
      endif

      if (dem_s(j,ilev) .lt. 0. .or. dem_s(j,ilev) .gt. 1.) then
        !             ' error = stratiform cloud emissivity less than zero'
        !             ' error = stratiform cloud emissivity greater than 1'
        rangevec(j)=rangevec(j)+16
      endif

      if (dem_c(j,ilev) .lt. 0. .or. dem_c(j,ilev) .gt. 1.) then
        !             ' error = convective cloud emissivity less than zero'
        !             ' error = convective cloud emissivity greater than 1'
        rangevec(j)=rangevec(j)+32
      endif
    enddo

    rangeerror=0
    do j=1,npoints
      rangeerror=rangeerror+rangevec(j)
    enddo

    if (rangeerror.ne.0) then
      write (6,*) 'Input variable out of range'
      write (6,*) 'rangevec:'
      write (6,*) rangevec
      call sys_flush(6)
      JERR=1 ; return  !STOP
    endif
  enddo

  do ibox=1,ncol
    do j=1,npoints
      boxpos(j,ibox)=(ibox-.5)*boxarea
    enddo
  enddo

  !     ---------------------------------------------------!
  !     Initialise working variables
  !     ---------------------------------------------------!

  !     Initialised frac_out to zero

  do ilev=1,nlev
    do ibox=1,ncol
      do j=1,npoints
        frac_out(j,ibox,ilev)=0.0
      enddo
    enddo
  enddo

  !      if (ncolprint.ne.0) then
  !        write (6,'(a)') 'frac_out_pp_rev:'
  !          do j=1,npoints,1000
  !          write(6,'(a10)') 'j='
  !          write(6,'(8I10)') j
  !          write (6,'(8f5.2)')
  !     &     ((frac_out(j,ibox,ilev),ibox=1,ncolprint),ilev=1,nlev)
  !
  !          enddo
  !        write (6,'(a)') 'ncol:'
  !        write (6,'(I3)') ncol
  !      endif
  !      if (ncolprint.ne.0) then
  !        write (6,'(a)') 'last_frac_pp:'
  !          do j=1,npoints,1000
  !          write(6,'(a10)') 'j='
  !          write(6,'(8I10)') j
  !          write (6,'(8f5.2)') (tca(j,0))
  !          enddo
  !      endif

  !     ---------------------------------------------------!
  !     ALLOCATE CLOUD INTO BOXES, FOR NCOLUMNS, NLEVELS
  !     frac_out is the array that contains the information
  !     where 0 is no cloud, 1 is a stratiform cloud and 2 is a
  !     convective cloud

  !loop over vertical levels
  ilev_loop: do ilev = 1,nlev

    !     Initialise threshold

    if (ilev.eq.1) then
      ! If max overlap
      if (overlap.eq.1) then
        ! select pixels spread evenly
        ! across the gridbox
        do ibox=1,ncol
          do j=1,npoints
            threshold(j,ibox)=boxpos(j,ibox)
          enddo
        enddo
      else
        do ibox=1,ncol
          !                include 'congvec.f'
          do j=1,npoints
            ran(j)=randu(xx)
          end do
          ! select random pixels from the non-convective
          ! part the gridbox ( some will be converted into
          ! convective pixels below )
          do j=1,npoints
            threshold(j,ibox)= &
                 cca(j,ilev)+(1-cca(j,ilev))*ran(j)
          enddo
        enddo
      end if
      !            IF (ncolprint.ne.0) then
      !              write (6,'(a)') 'threshold_nsf2:'
      !                do j=1,npoints,1000
      !                write(6,'(a10)') 'j='
      !                write(6,'(8I10)') j
      !                write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
      !                enddo
      !            END IF
    end if

    !        IF (ncolprint.ne.0) then
    !            write (6,'(a)') 'ilev:'
    !            write (6,'(I2)') ilev
    !        END IF

    do ibox=1,ncol

      ! All versions
      do j=1,npoints
        if (boxpos(j,ibox).le.cca(j,ilev)) then
          maxocc(j,ibox) = 1
        else
          maxocc(j,ibox) = 0
        end if
      enddo

      ! Max overlap
      if (overlap.eq.1) then
        do j=1,npoints
          threshold_min(j,ibox)=cca(j,ilev)
          maxosc(j,ibox)=1
        enddo
      endif

      ! Random overlap
      if (overlap.eq.2) then
        do j=1,npoints
          threshold_min(j,ibox)=cca(j,ilev)
          maxosc(j,ibox)=0
        enddo
      endif

      ! Max/Random overlap
      if (overlap.eq.3) then
        do j=1,npoints
          threshold_min(j,ibox)=max(cca(j,ilev), &
               min(tca(j,ilev-1),tca(j,ilev)))
          if (threshold(j,ibox) &
               .lt.min(tca(j,ilev-1),tca(j,ilev)) &
               .and.(threshold(j,ibox).gt.cca(j,ilev))) then
            maxosc(j,ibox)= 1
          else
            maxosc(j,ibox)= 0
          end if
        enddo
      endif

      ! Reset threshold

      !          include 'congvec.f'
      do j=1,npoints
        ran(j)=randu(xx)
      end do

      do j=1,npoints
        threshold(j,ibox)= &
                                !if max overlapped conv cloud
             maxocc(j,ibox) * ( &
             boxpos(j,ibox) &
             ) + &
                                !else
             (1-maxocc(j,ibox)) * ( &
                                !if max overlapped strat cloud
             (maxosc(j,ibox)) * ( &
                                !threshold=boxpos
             threshold(j,ibox) &
             ) + &
                                !else
             (1-maxosc(j,ibox)) * ( &
                                !threshold_min=random[thrmin,1]
             threshold_min(j,ibox)+ &
             (1-threshold_min(j,ibox))*ran(j) &
             ) &
             )
      enddo

    end do ! ibox

    !          Fill frac_out with 1's where tca is greater than the threshold

    do ibox=1,ncol
      do j=1,npoints
        if (tca(j,ilev).gt.threshold(j,ibox)) then
          frac_out(j,ibox,ilev)=1
        else
          frac_out(j,ibox,ilev)=0
        end if
      enddo
    end do

    !         Code to partition boxes into startiform and convective parts
    !         goes here

    do ibox=1,ncol
      do j=1,npoints
        if (threshold(j,ibox).le.cca(j,ilev)) then
          ! = 2 IF threshold le cca(j)
          frac_out(j,ibox,ilev) = 2
        else
          ! = the same IF NOT threshold le cca(j)
          frac_out(j,ibox,ilev) = frac_out(j,ibox,ilev)
        end if
      enddo
    end do

    !         Set last_frac to tca at this level, so as to be tca
    !         from last level next time round

    !          if (ncolprint.ne.0) then
    !
    !            do j=1,npoints ,1000
    !            write(6,'(a10)') 'j='
    !            write(6,'(8I10)') j
    !            write (6,'(a)') 'last_frac:'
    !            write (6,'(8f5.2)') (tca(j,ilev-1))
    !
    !            write (6,'(a)') 'cca:'
    !            write (6,'(8f5.2)') (cca(j,ilev),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'max_overlap_cc:'
    !            write (6,'(8f5.2)') (maxocc(j,ibox),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'max_overlap_sc:'
    !            write (6,'(8f5.2)') (maxosc(j,ibox),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'threshold_min_nsf2:'
    !            write (6,'(8f5.2)') (threshold_min(j,ibox),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'threshold_nsf2:'
    !            write (6,'(8f5.2)') (threshold(j,ibox),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'frac_out_pp_rev:'
    !            write (6,'(8f5.2)')
    !     &       ((frac_out(j,ibox,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
    !          enddo
    !          endif

  end do ilev_loop

      !
      !     ---------------------------------------------------!


      !
      !     ---------------------------------------------------!
      !     COMPUTE CLOUD OPTICAL DEPTH FOR EACH COLUMN and
      !     put into vector tau

      !initialize tau and albedocld to zero
  do ibox=1,ncol
    do j=1,npoints
      tau(j,ibox)=0.
      albedocld(j,ibox)=0.
      boxtau(j,ibox)=0.
      boxptop(j,ibox)=0.
      box_cloudy(j,ibox)=.false.
    enddo
  end do

        !compute total cloud optical depth for each column
  do ilev=1,nlev
    !increment tau for each of the boxes
    do ibox=1,ncol
      do j=1,npoints
        if (frac_out(j,ibox,ilev).eq.1) then
          tau(j,ibox)=tau(j,ibox) &
               + dtau_s(j,ilev)
        endif
        if (frac_out(j,ibox,ilev).eq.2) then
          tau(j,ibox)=tau(j,ibox) &
               + dtau_c(j,ilev)
        end if
      enddo
    enddo ! ibox
  enddo ! ilev
  !          if (ncolprint.ne.0) then
  !
  !              do j=1,npoints ,1000
  !                write(6,'(a10)') 'j='
  !                write(6,'(8I10)') j
  !                write(6,'(i2,1X,8(f7.2,1X))')
  !     &          ilev,
  !     &          (tau(j,ibox),ibox=1,ncolprint)
  !              enddo
  !          endif
  !
  !     ---------------------------------------------------!


  !
  !     ---------------------------------------------------!
  !     COMPUTE INFRARED BRIGHTNESS TEMPERATURES
  !     AND CLOUD TOP TEMPERATURE SATELLITE SHOULD SEE
  !
  !     again this is only done if top_height = 1 or 3
  !
  !     fluxtop is the 10.5 micron radiance at the top of the
  !              atmosphere
  !     trans_layers_above is the total transmissivity in the layers
  !             above the current layer
  !     fluxtop_clrsky(j) and trans_layers_above_clrsky(j) are the clear
  !             sky versions of these quantities.

  if (top_height .eq. 1 .or. top_height .eq. 3) then


    !----------------------------------------------------------------------
    !
    !             DO CLEAR SKY RADIANCE CALCULATION FIRST
    !
    !compute water vapor continuum emissivity
    !this treatment follows Schwarkzopf and Ramasamy
    !JGR 1999,vol 104, pages 9467-9499.
    !the emissivity is calculated at a wavenumber of 955 cm-1,
    !or 10.47 microns
    !        wtmair = 28.9644
    !        wtmh20 = 18.01534
    !        Navo = 6.023E+23
    !        grav = 9.806650E+02
    !        pstd = 1.013250E+06
    !        t0 = 296.
    !        if (ncolprint .ne. 0)
    !     &         write(6,*)  'ilev   pw (kg/m2)   tauwv(j)      dem_wv'
    do ilev=1,nlev
      do j=1,npoints
        !press and dpress are dyne/cm2 = Pascals *10
        press(j) = pfull(j,ilev)*10.
        dpress(j) = (phalf(j,ilev+1)-phalf(j,ilev))*10
        !atmden = g/cm2 = kg/m2 / 10
        ! minor GISS changes correct for unit difference in bygrav, bymrat,t0bypstd
        atmden(j) = dpress(j)*bygrav*0.01d0
        rvh20(j) = qv(j,ilev)*bymrat    !wtmair/wtmh20
        wk(j) = rvh20(j)*Navo*atmden(j)/wtmair
        !               rhoave(j) = (press(j)/pstd)*(t0/at(j,ilev))
        rhoave(j) = (press(j)/at(j,ilev))*t0bypstd
        rh20s(j) = rvh20(j)*rhoave(j)
        rfrgn(j) = rhoave(j)-rh20s(j)
        tmpexp(j) = exp(-0.02d0*(at(j,ilev)-t0))
        tauwv(j) = wk(j)*1.d-20*( &
             (0.0224697d0*rh20s(j)*tmpexp(j)) + &
             (3.41817d-7*rfrgn(j)) )*0.98d0
        dem_wv(j,ilev) = 1. - exp( -1. * tauwv(j))
      enddo
      !               if (ncolprint .ne. 0) then
      !               do j=1,npoints ,1000
      !               write(6,'(a10)') 'j='
      !               write(6,'(8I10)') j
      !               write(6,'(i2,1X,3(f8.3,3X))') ilev,
      !     &           qv(j,ilev)*(phalf(j,ilev+1)-phalf(j,ilev))*bygrav,
      !     &           tauwv(j),dem_wv(j,ilev)
      !               enddo
      !             endif
    end do

    !initialize variables
    do j=1,npoints
      fluxtop_clrsky(j) = 0.
      trans_layers_above_clrsky(j)=1.
    enddo

    do ilev=1,nlev
      do j=1,npoints

        ! Black body emission at temperature of the layer
        ! 1307.27 = Planck constant c1 by wavelength 10.5
        bb(j)=1 / ( exp(pc1bylam/at(j,ilev)) - 1. )

        ! increase TOA flux by flux emitted from layer
        ! times total transmittance in layers above

        fluxtop_clrsky(j) = fluxtop_clrsky(j) &
             + dem_wv(j,ilev)*bb(j)*trans_layers_above_clrsky(j)

        ! update trans_layers_above with transmissivity
        ! from this layer for next time around loop

        trans_layers_above_clrsky(j)= &
             trans_layers_above_clrsky(j)*(1.-dem_wv(j,ilev))


      enddo
      !            if (ncolprint.ne.0) then
      !             do j=1,npoints ,1000
      !              write(6,'(a10)') 'j='
      !              write(6,'(8I10)') j
      !              write (6,'(a)') 'ilev:'
      !              write (6,'(I2)') ilev
      !
      !              write (6,'(a)')
      !     &        'emiss_layer,100.*bb(j),100.*f,total_trans:'
      !              write (6,'(4(f7.2,1X))') dem_wv(j,ilev),100.*bb(j),
      !     &             100.*fluxtop_clrsky(j),trans_layers_above_clrsky(j)
      !             enddo
      !            endif

    enddo   !loop over level

    do j=1,npoints
      !add in surface emission
      bb(j)=1/( exp(pc1bylam/skt(j)) - 1. )
      !bb(j)=5.67d-8*skt(j)**4

      fluxtop_clrsky(j) = fluxtop_clrsky(j) + emsfc_lw * bb(j) &
           * trans_layers_above_clrsky(j)
    enddo

    !        if (ncolprint.ne.0) then
    !        do j=1,npoints ,1000
    !          write(6,'(a10)') 'j='
    !          write(6,'(8I10)') j
    !          write (6,'(a)') 'id:'
    !          write (6,'(a)') 'surface'
    !
    !          write (6,'(a)') 'emsfc,100.*bb(j),100.*f,total_trans:'
    !          write (6,'(4(f7.2,1X))') emsfc_lw,100.*bb(j),
    !     &      100.*fluxtop_clrsky(j),
    !     &       trans_layers_above_clrsky(j)
    !        enddo
    !      endif


    !
    !           END OF CLEAR SKY CALCULATION
    !
    !----------------------------------------------------------------



    !        if (ncolprint.ne.0) then
    !
    !        do j=1,npoints ,1000
    !            write(6,'(a10)') 'j='
    !            write(6,'(8I10)') j
    !            write (6,'(a)') 'ts:'
    !            write (6,'(8f7.2)') (skt(j),ibox=1,ncolprint)
    !
    !            write (6,'(a)') 'ta_rev:'
    !            write (6,'(8f7.2)')
    !     &       ((at(j,ilev2),ibox=1,ncolprint),ilev2=1,nlev)
    !
    !        enddo
    !        endif
    !loop over columns
    do ibox=1,ncol
      do j=1,npoints
        fluxtop(j,ibox)=0.
        trans_layers_above(j,ibox)=1.
      enddo
    enddo

    do ilev=1,nlev
      do j=1,npoints
        ! Black body emission at temperature of the layer

        bb(j)=1 / ( exp(pc1bylam/at(j,ilev)) - 1. )
        !bb(j)= 5.67d-8*at(j,ilev)**4
      enddo

      do ibox=1,ncol
        do j=1,npoints

          ! emissivity for point in this layer
          if (frac_out(j,ibox,ilev).eq.1) then
            dem(j,ibox)= 1. - &
                 ( (1. - dem_wv(j,ilev)) * (1. -  dem_s(j,ilev)) )
          else if (frac_out(j,ibox,ilev).eq.2) then
            dem(j,ibox)= 1. - &
                 ( (1. - dem_wv(j,ilev)) * (1. -  dem_c(j,ilev)) )
          else
            dem(j,ibox)=  dem_wv(j,ilev)
          end if


          ! increase TOA flux by flux emitted from layer
          ! times total transmittance in layers above

          fluxtop(j,ibox) = fluxtop(j,ibox) &
               + dem(j,ibox) * bb(j) &
               * trans_layers_above(j,ibox)

          ! update trans_layers_above with transmissivity
          ! from this layer for next time around loop

          trans_layers_above(j,ibox)= &
               trans_layers_above(j,ibox)*(1.-dem(j,ibox))

        enddo ! j
      enddo ! ibox

      !            if (ncolprint.ne.0) then
      !              do j=1,npoints,1000
      !              write (6,'(a)') 'ilev:'
      !              write (6,'(I2)') ilev
      !
      !              write(6,'(a10)') 'j='
      !              write(6,'(8I10)') j
      !              write (6,'(a)') 'emiss_layer:'
      !              write (6,'(8f7.2)') (dem(j,ibox),ibox=1,ncolprint)
      !
      !              write (6,'(a)') '100.*bb(j):'
      !              write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
      !
      !              write (6,'(a)') '100.*f:'
      !              write (6,'(8f7.2)')
      !     &         (100.*fluxtop(j,ibox),ibox=1,ncolprint)
      !
      !              write (6,'(a)') 'total_trans:'
      !              write (6,'(8f7.2)')
      !     &          (trans_layers_above(j,ibox),ibox=1,ncolprint)
      !            enddo
      !          endif

    enddo ! ilev


    do j=1,npoints
      !add in surface emission
      bb(j)=1/( exp(pc1bylam/skt(j)) - 1. )
      !bb(j)=5.67d-8*skt(j)**4
    end do

    do ibox=1,ncol
      do j=1,npoints

        !add in surface emission

        fluxtop(j,ibox) = fluxtop(j,ibox) &
             + emsfc_lw * bb(j) &
             * trans_layers_above(j,ibox)

      end do
    end do

    !        if (ncolprint.ne.0) then
    !
    !          do j=1,npoints ,1000
    !          write(6,'(a10)') 'j='
    !          write(6,'(8I10)') j
    !          write (6,'(a)') 'id:'
    !          write (6,'(a)') 'surface'
    !
    !          write (6,'(a)') 'emiss_layer:'
    !          write (6,'(8f7.2)') (dem(1,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') '100.*bb(j):'
    !          write (6,'(8f7.2)') (100.*bb(j),ibox=1,ncolprint)
    !
    !          write (6,'(a)') '100.*f:'
    !          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
    !          end do
    !      endif

    !now that you have the top of atmosphere radiance account
    !for ISCCP procedures to determine cloud top temperature

    !account for partially transmitting cloud recompute flux
    !ISCCP would see assuming a single layer cloud
    !note choice here of 2.13, as it is primarily ice
    !clouds which have partial emissivity and need the
    !adjustment performed in this section
    !
    !If it turns out that the cloud brightness temperature
    !is greater than 260K, then the liquid cloud conversion
    !factor of 2.56 is used.
    !
    !Note that this is discussed on pages 85-87 of
    !the ISCCP D level documentation (Rossow et al. 1996)

    do j=1,npoints
      !compute minimum brightness temperature and optical depth
      btcmin(j) = 1. /  ( exp(pc1bylam/(attrop(j)-5.)) - 1. )
    enddo
    do ibox=1,ncol
      do j=1,npoints
        transmax(j) = (fluxtop(j,ibox)-btcmin(j)) &
             /(fluxtop_clrsky(j)-btcmin(j))
        !note that the initial setting of tauir(j) is needed so that
        !tauir(j) has a realistic value should the next if block be
        !bypassed
        tauir(j) = tau(j,ibox) * byic
        taumin(j) = -1. * log(max(min(transmax(j),0.9999999d0),0.001d0))

      enddo

      if (top_height .eq. 1) then
        do j=1,npoints
          if (transmax(j) .gt. 0.001 .and. &
               transmax(j) .le. 0.9999999) then
            fluxtopinit(j) = fluxtop(j,ibox)
            tauir(j) = tau(j,ibox) *byic
          endif
        enddo
        do icycle=1,2
          do j=1,npoints
            if (tau(j,ibox) .gt. (tauchk            )) then
              if (transmax(j) .gt. 0.001 .and. &
                   transmax(j) .le. 0.9999999) then
                emcld(j,ibox) = 1. - exp(-1. * tauir(j)  )
                fluxtop(j,ibox) = fluxtopinit(j) - &
                     ((1.-emcld(j,ibox))*fluxtop_clrsky(j))
                fluxtop(j,ibox)=max(1.d-06, &
                     (fluxtop(j,ibox)/emcld(j,ibox)))
                tb(j,ibox)= pc1bylam &
                     / (log(1. + (1./fluxtop(j,ibox))))
                if (tb(j,ibox) .gt. 260.) then
                  tauir(j) = tau(j,ibox) * bywc
                end if
              end if
            end if
          enddo
        enddo

      endif

      do j=1,npoints
        if (tau(j,ibox) .gt. (tauchk            )) then
          !cloudy box
          tb(j,ibox)= pc1bylam/ (log(1. + (1./fluxtop(j,ibox))))
          if (top_height.eq.1.and.tauir(j).lt.taumin(j)) then
            tb(j,ibox) = attrop(j) - 5.
            tau(j,ibox) = 2.13d0*taumin(j)
          end if
        else
          !clear sky brightness temperature
          tb(j,ibox) = pc1bylam/(log(1.+(1./fluxtop_clrsky(j))))
        end if
      enddo ! j
    enddo ! ibox

    !        if (ncolprint.ne.0) then
    !
    !          do j=1,npoints,1000
    !          write(6,'(a10)') 'j='
    !          write(6,'(8I10)') j
    !
    !          write (6,'(a)') 'attrop:'
    !          write (6,'(8f7.2)') (attrop(j))
    !
    !          write (6,'(a)') 'btcmin:'
    !          write (6,'(8f7.2)') (btcmin(j))
    !
    !          write (6,'(a)') 'fluxtop_clrsky*100:'
    !          write (6,'(8f7.2)')
    !     &      (100.*fluxtop_clrsky(j))
    !
    !          write (6,'(a)') '100.*f_adj:'
    !          write (6,'(8f7.2)') (100.*fluxtop(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'transmax:'
    !          write (6,'(8f7.2)') (transmax(ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'tau:'
    !          write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'emcld:'
    !          write (6,'(8f7.2)') (emcld(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'total_trans:'
    !          write (6,'(8f7.2)')
    !     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'total_emiss:'
    !          write (6,'(8f7.2)')
    !     &        (1.0-trans_layers_above(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'total_trans:'
    !          write (6,'(8f7.2)')
    !     &        (trans_layers_above(j,ibox),ibox=1,ncolprint)
    !
    !          write (6,'(a)') 'ppout:'
    !          write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
    !          enddo ! j
    !      endif

  end if

  !     ---------------------------------------------------!

  !
  !     ---------------------------------------------------!
  !     DETERMINE CLOUD TOP PRESSURE
  !
  !     again the 2 methods differ according to whether
  !     or not you use the physical cloud top pressure (top_height = 2)
  !     or the radiatively determined cloud top pressure (top_height = 1 or 3)
  !

  !compute cloud top pressure
  do ibox=1,ncol
    !segregate according to optical thickness
    if (top_height .eq. 1 .or. top_height .eq. 3) then
      !find level whose temperature
      !most closely matches brightness temperature
      do j=1,npoints
        nmatch(j)=0
      enddo
      do ilev=1,nlev-1
        !cdir nodep
        do j=1,npoints
          if ((at(j,ilev)   .ge. tb(j,ibox) .and. &
               at(j,ilev+1) .lt. tb(j,ibox)) .or. &
               (at(j,ilev) .le. tb(j,ibox) .and. &
               at(j,ilev+1) .gt. tb(j,ibox))) then

            nmatch(j)=nmatch(j)+1
            if(abs(at(j,ilev)-tb(j,ibox)) .lt. &
                 abs(at(j,ilev+1)-tb(j,ibox))) then
              match(j,nmatch(j))=ilev
            else
              match(j,nmatch(j))=ilev+1
            end if
          end if
        enddo
      end do

      do j=1,npoints
        if (nmatch(j) .ge. 1) then
          ptop(j,ibox)=pfull(j,match(j,nmatch(j)))
          levmatch(j,ibox)=match(j,nmatch(j))
        else
          if (tb(j,ibox) .lt. atmin(j)) then
            ptop(j,ibox)=ptrop(j)
            levmatch(j,ibox)=itrop(j)
          end if
          if (tb(j,ibox) .gt. atmax(j)) then
            ptop(j,ibox)=pfull(j,nlev)
            levmatch(j,ibox)=nlev
          end if
        end if
      enddo ! j

    else ! if (top_height .eq. 1 .or. top_height .eq. 3)

      do j=1,npoints
        ptop(j,ibox)=0.
      enddo
      do ilev=1,nlev
        do j=1,npoints
          if ((ptop(j,ibox) .eq. 0. ) &
               .and.(frac_out(j,ibox,ilev) .ne. 0)) then
            ptop(j,ibox)=pfull(j,ilev)
            levmatch(j,ibox)=ilev
          end if
        end do
      end do
    end if

    do j=1,npoints
      if (tau(j,ibox) .le. (tauchk            )) then
        ptop(j,ibox)=0.
        levmatch(j,ibox)=0
      endif
    enddo

  end do

  !
  !
  !     ---------------------------------------------------!


  !
  !     ---------------------------------------------------!
  !     DETERMINE ISCCP CLOUD TYPE FREQUENCIES
  !
  !     Now that ptop and tau have been determined,
  !     determine amount of each of the 49 ISCCP cloud
  !     types
  !
  !     Also compute grid box mean cloud top pressure and
  !     optical thickness.  The mean cloud top pressure and
  !     optical thickness are averages over the cloudy
  !     area only. The mean cloud top pressure is a linear
  !     average of the cloud top pressures.  The mean cloud
  !     optical thickness is computed by converting optical
  !     thickness to an albedo, averaging in albedo units,
  !     then converting the average albedo back to a mean
  !     optical thickness.
  !

  !compute isccp frequencies

  !reset frequencies
  do ilev=1,7
    do ilev2=1,7
      do j=1,npoints !
        fq_isccp(j,ilev,ilev2)=0.
      enddo
    end do
  end do

  !reset variables need for averaging cloud properties
  do j=1,npoints
    totalcldarea(j) = 0.
    meanalbedocld(j) = 0.
    meanptop(j) = 0.
    meantaucld(j) = 0.
  enddo ! j

  !      boxarea = 1./real(ncol,kind=8)

  do ibox=1,ncol
    do j=1,npoints

      !          if (tau(j,ibox) .gt. (tauchk            )
      !     &      .and. ptop(j,ibox) .gt. 0.) then
      !              box_cloudy(j,ibox)=.true.
      !          endif

      box_cloudy(j,ibox)= (tau(j,ibox) .gt. tauchk .and. ptop(j,ibox &
           ) .gt. 0.)

      if (box_cloudy(j,ibox)) then

        ! totalcldarea always diagnosed day or night
        totalcldarea(j) = totalcldarea(j) + boxarea

        if (sunlit(j).eq.1) then

          ! tau diagnostics only with sunlight

          boxtau(j,ibox) = tau(j,ibox)

          !convert optical thickness to albedo
          albedocld(j,ibox) &
               =real(invtau(min(nint(100.*tau(j,ibox)),45000)))

          !contribute to averaging
          meanalbedocld(j) = meanalbedocld(j) &
               +albedocld(j,ibox)*boxarea

        endif

      endif

      if (sunlit(j).eq.1 .or. top_height .eq. 3) then

        !convert ptop to millibars
        ptop(j,ibox)=ptop(j,ibox) / 100.

        !save for output cloud top pressure and optical thickness
        boxptop(j,ibox) = ptop(j,ibox)

        if (box_cloudy(j,ibox)) then

          meanptop(j) = meanptop(j) + ptop(j,ibox)*boxarea

          !reset itau(j), ipres(j)
          itau(j) = 0
          ipres(j) = 0

          !determine optical depth category
          if (tau(j,ibox) .lt. isccp_taumin) then
            itau(j)=1
          else if (tau(j,ibox) .ge. isccp_taumin &
               .and. tau(j,ibox) .lt. 1.3) then
            itau(j)=2
          else if (tau(j,ibox) .ge. 1.3 &
               .and. tau(j,ibox) .lt. 3.6) then
            itau(j)=3
          else if (tau(j,ibox) .ge. 3.6 &
               .and. tau(j,ibox) .lt. 9.4) then
            itau(j)=4
          else if (tau(j,ibox) .ge. 9.4 &
               .and. tau(j,ibox) .lt. 23.) then
            itau(j)=5
          else if (tau(j,ibox) .ge. 23. &
               .and. tau(j,ibox) .lt. 60.) then
            itau(j)=6
          else if (tau(j,ibox) .ge. 60.) then
            itau(j)=7
          end if

          !determine cloud top pressure category
          if (    ptop(j,ibox) .gt. 0. &
               .and.ptop(j,ibox) .lt. 180.) then
            ipres(j)=1
          else if(ptop(j,ibox) .ge. 180. &
               .and.ptop(j,ibox) .lt. 310.) then
            ipres(j)=2
          else if(ptop(j,ibox) .ge. 310. &
               .and.ptop(j,ibox) .lt. 440.) then
            ipres(j)=3
          else if(ptop(j,ibox) .ge. 440. &
               .and.ptop(j,ibox) .lt. 560.) then
            ipres(j)=4
          else if(ptop(j,ibox) .ge. 560. &
               .and.ptop(j,ibox) .lt. 680.) then
            ipres(j)=5
          else if(ptop(j,ibox) .ge. 680. &
               .and.ptop(j,ibox) .lt. 800.) then
            ipres(j)=6
          else if(ptop(j,ibox) .ge. 800.) then
            ipres(j)=7
          end if

          !update frequencies
          if(ipres(j) .gt. 0.and.itau(j) .gt. 0) &
               fq_isccp(j,itau(j),ipres(j))= &
               fq_isccp(j,itau(j),ipres(j))+ boxarea
        end if

      end if

    enddo ! j
  end do ! ibox

  !compute mean cloud properties
  do j=1,npoints
    nbox(j)=totalcldarea(j)*ncol    ! for backwards compatibility
    if (totalcldarea(j) .gt. 0.) then
      meanptop(j) = meanptop(j) / totalcldarea(j)
      if (sunlit(j).eq.1) then
        meanalbedocld(j) = meanalbedocld(j) / totalcldarea(j)
        meantaucld(j) = tautab(min(255,max(1,nint(meanalbedocld(j)))))
      end if
    end if
  enddo ! j

  !
  !     ---------------------------------------------------!

  !     ---------------------------------------------------!
  !     OPTIONAL PRINTOUT OF DATA TO CHECK PROGRAM
  !
  !      if (debugcol.ne.0) then
  !
  !         do j=1,npoints,debugcol
  !
  !            !produce character output
  !            do ilev=1,nlev
  !              do ibox=1,ncol
  !                   acc(ilev,ibox)=0
  !              enddo
  !            enddo
  !
  !            do ilev=1,nlev
  !              do ibox=1,ncol
  !                   acc(ilev,ibox)=frac_out(j,ibox,ilev)*2
  !                   if (levmatch(j,ibox) .eq. ilev)
  !     &                 acc(ilev,ibox)=acc(ilev,ibox)+1
  !              enddo
  !            enddo
  !
  !             !print test
  !
  !          write(ftn09,11) j
  !11        format('ftn09.',i4.4)
  !          open(9, FILE=ftn09, FORM='FORMATTED')
  !
  !             write(9,'(a1)') ' '
  !             write(9,'(10i5)')
  !     &                  (ilev,ilev=5,nlev,5)
  !             write(9,'(a1)') ' '
  !
  !             do ibox=1,ncol
  !               write(9,'(40(a1),1x,40(a1))')
  !     &           (cchar_realtops(acc(ilev,ibox)+1),ilev=1,nlev)
  !     &           ,(cchar(acc(ilev,ibox)+1),ilev=1,nlev)
  !             end do
  !             close(9)
  !
  !             if (ncolprint.ne.0) then
  !               write(6,'(a1)') ' '
  !                    write(6,'(a2,1X,5(a7,1X),a50)')
  !     &                  'ilev',
  !     &                  'pfull','at',
  !     &                  'cc*100','dem_s','dtau_s',
  !     &                  'cchar'

  !               do 4012 ilev=1,nlev
  !                    write(6,'(60i2)') (box(i,ilev),i=1,ncolprint)
  !                   write(6,'(i2,1X,5(f7.2,1X),50(a1))')
  !     &                  ilev,
  !     &                  pfull(j,ilev)/100.,at(j,ilev),
  !     &                  cc(j,ilev)*100.0,dem_s(j,ilev),dtau_s(j,ilev)
  !     &                  ,(cchar(acc(ilev,ibox)+1),ibox=1,ncolprint)
  !4012           continue
  !               write (6,'(a)') 'skt(j):'
  !               write (6,'(8f7.2)') skt(j)
  !
  !               write (6,'(8I7)') (ibox,ibox=1,ncolprint)
  !
  !               write (6,'(a)') 'tau:'
  !               write (6,'(8f7.2)') (tau(j,ibox),ibox=1,ncolprint)
  !
  !               write (6,'(a)') 'tb:'
  !               write (6,'(8f7.2)') (tb(j,ibox),ibox=1,ncolprint)
  !
  !               write (6,'(a)') 'ptop:'
  !               write (6,'(8f7.2)') (ptop(j,ibox),ibox=1,ncolprint)
  !             endif
  !
  !        enddo
  !
  !      end if

  return
end

subroutine get_dq_cond(sm,qm,plk,mass,lhx,pl,   & ! input
     dqsum,fcond)         ! output
!@sum get_dq calculation of condensation from vapour
  use CONSTANT, only : bysha
  implicit none
!@var SM,QM heat and water content
  real*8, intent(IN) :: sm,qm
!@var MASS air mass
!@var PLK P**K for mid point
!@var LHX relevant latent heat
!@var PL mid point pressure
  real*8, intent(IN) :: plk,mass,lhx,pl
!@var DQSUM amount of water change (+ve is condensation)
!@var FCOND fractional amount of vapour that condenses
  real*8, intent(OUT) :: dqsum,fcond
  real*8 TP,QST,QSAT,DQSATDT,DQ,QMT,SLH
  integer N

  DQSUM=0.
  FCOND=0.
  if (QM.gt.0) then
    SLH=LHX*BYSHA
    QMT=QM
    TP=SM*PLK/MASS
    do N=1,3
      QST=QSAT(TP,LHX,PL)
      DQ=(QMT-MASS*QST)/(1.+SLH*QST*DQSATDT(TP,LHX))
      TP=TP+SLH*DQ/MASS
      QMT=QMT-DQ
      DQSUM=DQSUM+DQ
    end do
    !**** condensing (DQSUM > 0)
    DQSUM=max(0d0,min(DQSUM,QM))
    FCOND=DQSUM/QM
  end if

  return
end subroutine get_dq_cond

subroutine get_dq_evap(sm,qm,plk,mass,lhx,pl,cond,   & ! input
     dqsum,fevp) ! output
!@sum get_dq calculation of evaporation from condensate
  use CONSTANT, only : bysha
  implicit none
!@var SM,QM heat and water content
  real*8, intent(IN) :: sm,qm
!@var MASS air mass
!@var PLK P**K for mid point
!@var LHX relevant latent heat
!@var PL mid point pressure
!@var COND amount of condensate
  real*8, intent(IN) :: plk,mass,lhx,pl,cond
!@var DQSUM amount of water change (+ve is evaporation)
!@var FEVP fractional amount of condensate that evaporates
  real*8, intent(OUT) :: dqsum,fevp
  real*8 TP,QST,QSAT,DQSATDT,DQ,QMT,SLH
  integer N

  DQSUM=0.
  FEVP=0.
  if (COND.gt.0) then
    SLH=LHX*BYSHA
    QMT=QM
    TP=SM*PLK/MASS
    do N=1,3
      QST=QSAT(TP,LHX,PL)
      DQ=(QMT-MASS*QST)/(1.+SLH*QST*DQSATDT(TP,LHX))
      TP=TP+SLH*DQ/MASS
      QMT=QMT-DQ
      DQSUM=DQSUM-DQ
    end do
    !**** evaporating (DQ < 0, DQSUM > 0)
    DQSUM=max(0d0,min(DQSUM,COND))
    FEVP=DQSUM/COND
  end if

  return
end subroutine get_dq_evap
