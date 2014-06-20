#include "rundeck_opts.h"
#ifdef SKIP_TRACERS_RAD
#undef TRACERS_ON
#endif
!@sum RAD_DRV contains drivers for the radiation related routines
!@ver  2009/05/11
!@cont init_RAD, RADIA
C**** semi-random cloud overlap (computed opt.d+diagn)
C**** to be used with R99E or later radiation  routines.  carbon/2
C****

      SUBROUTINE CALC_ZENITH_ANGLE
!@sum calculate zenith angle for current time step
!@auth Gavin Schmidt (from RADIA)
      USE CONSTANT, only : twopi
      USE MODEL_COM, only : itime,nday,dtsrc, calendar
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE RAD_COM, only : cosz1
      USE RAD_COSZ0, only : coszt
      USE TimeInterval_mod
      IMPLICIT NONE
      INTEGER JTIME
      REAL*8 ROT1,ROT2
      type (TimeInterval) :: sPerDay

      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
      sPerDay = calendar%getSecondsPerDay()
      ROT2=ROT1+TWOPI*DTsrc/sPerDay%convertToReal()
      CALL COSZT (ROT1,ROT2,COSZ1)

      END SUBROUTINE CALC_ZENITH_ANGLE

      SUBROUTINE init_RAD(istart)
!@sum  init_RAD initialises radiation code
!@auth Original Development Team
!@calls RADPAR:RCOMP1, ORBPAR
      USE FILEMANAGER
      USE Dictionary_mod
      USE CONSTANT, only : grav,bysha,twopi,planet_name
      USE RESOLUTION, only : jm,lm,psf
      USE ATM_COM, only : t,pk,kradia,lm_req
      USE MODEL_COM, only : dtsrc,iyear1,modelEclock,master_yr
      USE MODEL_COM, only: orbit
      USE ATM_COM, only : pednl00
      USE DOMAIN_DECOMP_ATM, only : grid, write_parallel, am_i_root
     &     ,readt_parallel, getDomainBounds
#ifndef CUBED_SPHERE
      USE GEOM, only : lat_dg
#endif
      USE RADPAR, only : !rcomp1,writer,writet       ! routines
     &      PTLISO ,KTREND ,LMR=>NL, PLB, LS1_loc
     &     ,planck_tmin,planck_tmax
     &     ,transmission_corrections
     *     ,KCLDEM,KSIALB,KSOLAR, SHL, snoage_fac_max, KZSNOW
     *     ,KYEARS,KJDAYS,MADLUV, KYEARG,KJDAYG,MADGHG
     *     ,KYEARO,KJDAYO,MADO3M, KYEARA,KJDAYA,MADAER
     *     ,KYEARD,KJDAYD,MADDST, KYEARV,KJDAYV,MADVOL
     *     ,KYEARE,KJDAYE,MADEPS, KYEARR,KJDAYR
!g95     *     ,FSXAER,FTXAER    ! scaling (on/off) for default aerosols
     *     ,ITR,NTRACE        ! turning on options for extra aerosols
     *     ,FS8OPX,FT8OPX, TRRDRY,KRHTRA,TRADEN,REFDRY
     *     ,rcomp1, writer, writet
     *     ,FSTASC,FTTASC
#ifdef ALTER_RADF_BY_LAT
     *     ,FS8OPX_orig,FT8OPX_orig
#endif
      USE RAD_COM, only : rqt, s0x, co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx
     *     ,CH4X_RADoverCHEM
     *     ,s0_yr,s0_day,ghg_yr,ghg_day,volc_yr,volc_day,aero_yr,O3_yr
     *     ,H2ObyCH4,dH2O,h2ostratx,O3x,RHfix,CLDx,ref_mult,COSZ1
     *     ,CC_cdncx,OD_cdncx,cdncl,pcdnc,vcdnc
     *     ,cloud_rad_forc,aer_rad_forc
     *     ,PLB0,shl0  ! saved to avoid OMP-copyin of input arrays
     *     ,albsn_yr,dALBsnX,nradfrc
     *     ,rad_interact_aer,clim_interact_chem,rad_forc_lev,ntrix,wttr
     *     ,nrad_clay,variable_orb_par,orb_par_year_bp,orb_par,nrad
     *     ,radiationSetOrbit
#ifdef TRACERS_ON
     &     ,nTracerRadiaActive,tracerRadiaActiveFlag
#endif
#ifdef TRACERS_SPECIAL_Shindell
     *     ,maxNtraceFastj
#endif
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
#ifdef CHL_from_SeaWIFs
     *     ,iu_CHL
#endif
#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
     *     ,wfac
#endif
#ifdef OBIO_RAD_coupling
      USE RAD_COM, only : DIRVIS,FSRDIF,DIRNIR,DIFNIR
#endif
      use RAD_COSZ0, only : cosz_init
      USE CLOUDS_COM, only : llow
      USE DIAG_COM, only : iwrite,jwrite,itwrite,save3dAOD
#ifdef TRACERS_ON
      USE TRACER_COM, only: ntm
      USE TRACER_COM, only: n_BCIA, n_BCB, n_NO3p
      USE TRACER_COM, only: n_Clay, n_Silt1, n_Silt2, n_Silt3, n_Silt4
      USE TRACER_COM, only: n_SO4, n_Seasalt1, n_Seasalt2
      USE TRACER_COM, only: n_OCB, n_OCIA, n_Isopp1a, n_SO4, ntm_dust
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only: n_ASO4, n_ANACL, n_AECOB, n_AECIL,
     &     n_AOCOB, n_AOCIL, n_ADUST
#endif
#endif
#ifdef TRACERS_AMP
      USE AERO_CONFIG, only: nmodes
      USE TRACER_COM, only:
     *     n_N_AKK_1 ,n_N_ACC_1 ,n_N_DD1_1 ,n_N_DS1_1 ,n_N_DD2_1,
     *     n_N_DS2_1, n_N_SSA_1, n_N_SSC_1, n_N_OCC_1, n_N_BC1_1,
     *     n_N_BC2_1 ,n_N_BC3_1,
     *     n_N_DBC_1, n_N_BOC_1, n_N_BCS_1, n_N_MXX_1
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only: icomp
#endif
      use AerParam_mod, only : aermix
      use AerParam_mod, only: depoBC,depoBC_1990

      use AbstractOrbit_mod, only: AbstractOrbit
      IMPLICIT NONE

      integer, intent(in) :: istart

      INTEGER L,LR,n1,n,nn,iu2 ! LONR,LATR
      REAL*8 PLBx(LM+1),pyear
!@var NRFUN indices of unit numbers for radiation routines
      INTEGER NRFUN(14),IU,DONOTREAD
!@var RUNSTR names of files for radiation routines
      CHARACTER*5 :: RUNSTR(14) = (/"RADN1","RADN2","RADN3",
     *     "RADN4","RADN5","RADN6","RADN7","RADN8",
     *     "RADN9","RADNA","RADNB","RADNC","RADND",
     *     "RADNE"/)
!@var QBIN true if files for radiation input files are binary
      LOGICAL :: QBIN(14)=(/.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.
     *     ,.TRUE.,.TRUE.,.FALSE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)

#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
      integer, parameter :: nlt=33
      real*8 :: aw(nlt), bw(nlt), saw, sbw
      real*8 :: b0, b1, b2, b3, a0, a1, a2, a3, expterm, tlog, fac, rlam
      integer :: nl,ic , iu_bio, lambda, lam(nlt)
      character title*50
      data a0,a1,a2,a3 /0.9976d0, 0.2194d0,  5.554d-2,  6.7d-3 /
      data b0,b1,b2,b3 /5.026d0, -0.01138d0, 9.552d-6, -2.698d-9/
#endif

      character(len=300) :: out_line
      character*6 :: skip

      INTEGER :: I,J
      INTEGER :: I_0,I_1,J_0,J_1

C**** sync radiation parameters from input
      call sync_param( "NRAD", NRAD ) !!
      if (is_set_param("variable_orb_par")) then
        call get_param( "variable_orb_par", variable_orb_par )
      else
        if (master_yr == 0) then
          variable_orb_par=1
        else
          variable_orb_par=0
        endif
      endif
      if (is_set_param("orb_par_year_bp")) then
        call get_param( "orb_par_year_bp", orb_par_year_bp )
      else
        if (master_yr == 0) then
          orb_par_year_bp=0
        else
          orb_par_year_bp=1950-master_yr
        endif
      endif
      call sync_param( "orb_par", orb_par, 3 )
      call sync_param( "S0X", S0X )
      call sync_param( "CO2X", CO2X )
      call sync_param( "N2OX", N2OX )
      call sync_param( "CH4X", CH4X )
      call sync_param( "CH4X_RADoverCHEM", CH4X_RADoverCHEM )
      call sync_param( "CFC11X", CFC11X )
      call sync_param( "CFC12X", CFC12X )
      call sync_param( "XGHGX", XGHGX )
      call sync_param( "H2OstratX", H2OstratX )
      call sync_param( "O3X", O3X )
      call sync_param( "CLDX", CLDX )
      call sync_param( "H2ObyCH4", H2ObyCH4 )
      call get_param( "S0_yr", S0_yr, default=master_yr )
      if (is_set_param("S0_day")) then
        call get_param( "S0_day", S0_day )
      else
        if (s0_yr==0) s0_day=0 ! else use default value
      endif
      call get_param( "ghg_yr", ghg_yr, default=master_yr )
      if (is_set_param("ghg_day")) then
        call get_param( "ghg_day", ghg_day )
      else
        if (ghg_yr==0) ghg_day=0 ! else use default value
      endif
      call get_param( "volc_yr", volc_yr, default=master_yr )
      if (is_set_param("volc_day")) then
        call get_param( "volc_day", volc_day )
      else
        if (volc_yr==0) volc_day=0 ! else use default value
      endif
      call get_param( "aero_yr", aero_yr, default=master_yr )
      call sync_param( "madaer", madaer )
      call sync_param( "dALBsnX", dALBsnX )
      call get_param( "albsn_yr", albsn_yr, default=master_yr )
      call sync_param( "aermix", aermix , 13 )
      call sync_param( "REFdry", REFdry , 8 )
      call sync_param( "FS8OPX", FS8OPX , 8 )
      call sync_param( "FT8OPX", FT8OPX , 8 )
      call sync_param( "RHfix", RHfix )
      call sync_param( "CC_cdncx", CC_cdncx )
      call sync_param( "OD_cdncx", OD_cdncx )
      call get_param( "O3_yr", O3_yr, default=master_yr )
      if(planet_name.ne.'Earth') then
        PTLISO = .015d0*psf ! reasonable default
      endif
      call sync_param( "PTLISO", PTLISO )
      call sync_param( "KSOLAR", KSOLAR )
      call sync_param( "KSIALB", KSIALB )
      call sync_param( "KZSNOW", KZSNOW )
      call sync_param( "snoage_fac_max", snoage_fac_max )
      call sync_param( "nradfrc", nradfrc )
      if(snoage_fac_max.lt.0. .or. snoage_fac_max.gt.1.) then
        write(out_line,*) 'set 0<snoage_fac_max<1, not',snoage_fac_max
        call write_parallel(trim(out_line),unit=6)
        call stop_model('init_RAD: snoage_fac_max out of range',255)
      end if
      call sync_param( "rad_interact_aer", rad_interact_aer )
      call sync_param( "clim_interact_chem", clim_interact_chem )
      call sync_param( "rad_forc_lev", rad_forc_lev )
      call sync_param( "cloud_rad_forc", cloud_rad_forc )
      call sync_param( "aer_rad_forc", aer_rad_forc )
      call sync_param( "ref_mult", ref_mult )
      call sync_param( "save3dAOD", save3dAOD)
      REFdry = REFdry*ref_mult

      if(is_set_param('planck_tmin')) then
        call get_param('planck_tmin',planck_tmin)
      endif
      if(is_set_param('planck_tmax')) then
        call get_param('planck_tmax',planck_tmax)
      endif

      call getDomainBounds(grid,
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

      if(istart==2) then ! replace with cold vs warm start logic
C**** SET RADIATION EQUILIBRIUM TEMPERATURES FROM LAYER LM TEMPERATURE
        DO J=J_0,J_1
        DO I=I_0,I_1
          RQT(:,I,J)=T(I,J,LM)*PK(LM,I,J)
        ENDDO
        ENDDO
      endif


      call radiationSetOrbit(orbit)
      call cosz_init

C****
C**** SET THE CONTROL PARAMETERS FOR THE RADIATION (need mean pressures)
C****
      LMR=LM+LM_REQ
      PLB(1:LMR+1)=PEDNL00(1:LMR+1)
      DO L=1,LM
        PLBx(L)=PLB(L)           ! needed for CH4 prod. H2O
      END DO
      PLBx(LM+1)=0.
      DO LR=LM+1,LMR
        PLB0(LR-LM) = PLB(LR+1)
      END DO
      cdncl=0 ; call reterp(vcdnc,pcdnc,7, cdncl,plb,llow+2)

      KTREND=1   !  GHgas trends are determined by input file
!note KTREND=0 is a possible but virtually obsolete option
C****
C             Model Add-on Data of Extended Climatology Enable Parameter
C     MADO3M  = -1   Reads                      Ozone data the GCM way
C     MADAER  =  1   Reads   Tropospheric Aerosol climatology 1850-2050
C     MADAER  =  3   uses Koch,Bauer 2008 aerosol climatology 1890-2000
C     MADDST  =  1   Reads   Dust-windblown mineral climatology   RFILE6
C     MADVOL  =  1   Reads   Volcanic 1950-00 aerosol climatology RFILE7
C     MADEPS  =  1   Reads   Epsilon cloud heterogeniety data     RFILE8
C     MADLUV  =  1   Reads   Lean's SolarUV 1882-1998 variability RFILE9
C**** Radiative forcings are either constant = obs.value at given yr/day
C****    or time dependent (year=0); if day=0 an annual cycle is used
C****                                         even if the year is fixed
      KYEARS=s0_yr   ; KJDAYS=s0_day ;  MADLUV=1   ! solar 'constant'
      KYEARG=ghg_yr  ; KJDAYG=ghg_day              ! well-mixed GHGases
#ifndef ALTER_RADF_BY_LAT
      if(ghg_yr.gt.0)  MADGHG=0                    ! skip GHG-updating
#endif
      KYEARO=O3_yr   ; KJDAYO=0 ;       MADO3M=-1  ! ozone (ann.cycle)
      if(KYEARO.gt.0) KYEARO=-KYEARO              ! use ONLY KYEARO-data
      KYEARA=Aero_yr ; KJDAYA=0 ! MADAER=1 or 3, trop.aeros (ann.cycle)
      if(KYEARA.gt.0) KYEARA=-KYEARA              ! use ONLY KYEARA-data
      if(file_exists('TAero_SSA')) MADAER=3   ! one of the TAero_XXX set
      KYEARV=Volc_yr ; KJDAYV=Volc_day
      if(file_exists('RADN7')) MADVOL=1   ! Volc. Aerosols

!***  KYEARV=0 : use current year
!***  KYEARV<0 : use long term mean stratospheric aerosols (use -1)
!     Hack: KYEARV= -2000 and -2010 were used for 2 specific runs that
!           ended in 2100 and repeated some 20th century volcanos
!***  KYEARV=-2000: use volcanos from 100 yrs ago after 2000
!***  KYEARV=-2010: repeat 2nd half, then first half of 20th century
      if(KYEARV.le.-2000) KYEARV=0   ! use current year (before 2000)
C**** NO time history (yet), except for ann.cycle, for forcings below;
C****  if KJDAY?=day0 (1->365), data from that day are used all year
      KYEARD=0       ; KJDAYD=0 ;
      KYEARE=0       ; KJDAYE=0 ;
      KYEARR=0       ; KJDAYR=0           ! surf.reflectance (ann.cycle)
      KCLDEM=1                  ! 0:old 1:new LW cloud scattering scheme

      if(file_exists('DUSTaer'))   MADDST=1   ! Desert dust
      if(file_exists('RADN8'))     MADEPS=1   ! cloud Epsln - KCLDEP
      transmission_corrections = file_exists('RADN4')

C**** Aerosols:
C**** Currently there are five different default aerosol controls
C****   1:total 2:background+tracer 3:Climatology 4:dust 5:volcanic
C**** By adjusting FSXAER,FTXAER you can remove the default
C**** aerosols and replace them with your version if required
C**** (through TRACER in RADIA).
C**** FSXAER is for the shortwave,    FTXAER is for the longwave effects
caer  FSXAER = (/ 1.,1.,1.,1.,1. /) ; FTXAER = (/ 1.,1.,1.,1.,1. /)

C**** climatology aerosols are grouped into 6 types from 13 sources:
C****  Pre-Industrial+Natural 1850 Level  Industrial Process  BioMBurn
C****  ---------------------------------  ------------------  --------
C****   1    2    3    4    5    6    7    8    9   10   11   12   13
C****  SNP  SBP  SSP  ANP  ONP  OBP  BBP  SUI  ANI  OCI  BCI  OCB  BCB
C**** using the following default scaling/tuning factors  AERMIX(1-13)
C****  1.0, 1.0, .26, 1.0, 2.5, 2.5, 1.9, 1.0, 1.0, 2.5, 1.9, 2.5, 1.9
C**** The 8 groups are (adding dust and volcanic aerosols as 7. and 8.)
C**** 1. Sulfates (industr and natural), 2. Sea Salt, 3. Nitrates
C**** 4. Organic Carbons, 5. industr Black Carbons(BC), 6. Biomass BC
C**** 7. Dust aerosols, 8. Volcanic aerosols
C**** use FS8OPX and FT8OPX to enhance the optical effect; defaults:
caer  FS8OPX = (/1., 1., 1., 1., 2., 2.,    1.   ,   1./)     solar
caer  FT8OPX = (/1., 1., 1., 1., 1., 1.,    1.3d0,   1./)     thermal
!!!!! Note: FS|T8OPX(7-8) makes FS|TXAER(4-5) redundant.
C**** Particle sizes of the first 4 groups have RelHum dependence

C**** To add up to 8 further aerosols:
C****  1) set NTRACE to the number of extra aerosol fields
C****  2) ITR defines which set of Mie parameters get used, choose
C****     from the following:
C****     1 SO4,  2 seasalt, 3 nitrate, 4 OCX organic carbons
C****     5 BCI,  6 BCB,     7 dust,    8 H2SO4 volc
C****  2b) set up the indexing array NTRIX to map the RADIATION tracers
C****      to the main model tracers
C****  2c) set up the weighting array WTTR to weight main model tracers
C****
C****  3) Use FSTOPX/FTTOPX(1:NTRACE) to scale them in RADIA
C****  4) Set TRRDRY to dry radius
C****  5) Set KRHTRA=1 if aerosol has RH dependence, 0 if not
C**** Note: whereas FSXAER/FTXAER are global (shared), FSTOPX/FTTOPX
C****       have to be reset for each grid box to allow for the way it
C****       is used in RADIA (TRACERS_AEROSOLS_Koch)
caer   NTRACE = 0
caer   ITR = (/ 0,0,0,0, 0,0,0,0 /)
caer   TRRDRY=(/ .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0, .1d0/)
caer   KRHTRA=(/1,1,1,1,1,1,1,1/)

#if (defined TRACERS_AMP) || (defined TRACERS_AMP_M1)
      if (rad_interact_aer > 0) then
C                  SO4    SEA    NO3    OCX    BCI    BCB    DST   VOL
        FS8OPX = (/0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0 , 1d0/)
        FT8OPX = (/0d0,   0d0,   0d0,   0d0,   0d0,   0d0,   0d0,  1d0/)
      end if
      NTRACE=nmodes
      NTRIX(1:NMODES)=
     *     (/ n_N_AKK_1 ,n_N_ACC_1 ,n_N_DD1_1 ,n_N_DS1_1 ,n_N_DD2_1,
     *        n_N_DS2_1, n_N_SSA_1, n_N_SSC_1, n_N_OCC_1, n_N_BC1_1,
     *        n_N_BC2_1 ,n_N_BC3_1,
     *        n_N_DBC_1, n_N_BOC_1, n_N_BCS_1, n_N_MXX_1/)
#endif
#ifdef TRACERS_TOMAS
!TOMAS does not include NO3 AND VOL, which use its default radiation.
      if (rad_interact_aer > 0) then
C                  SO4    SEA    NO3    OCX    BCI    BCB    DST   VOL
        FS8OPX = (/0d0,   0d0,   1d0,   0d0,   0d0,   0d0,   0d0,  1d0/)
        FT8OPX = (/0d0,   0d0,   1d0,   0d0,   0d0,   0d0,   0d0,  1d0/)
      end if
      NTRACE=icomp-2
! ANUM(1) for internal-mixing case
! Others(ncomp-1) for external-mixing case.
      NTRIX(1:ntrace)=
     *     (/n_ASO4(1), n_ANACL(1), n_AECOB(1), n_AECIL(1),
     &     n_AOCOB(1), n_AOCIL(1), n_ADUST(1)/)
#endif
#ifdef TRACERS_AEROSOLS_Koch
      if (rad_interact_aer > 0) then  ! if BC's sol.effect are doubled:
#ifdef SULF_ONLY_AEROSOLS
        NTRACE=1
        FS8OPX(1:NTRACE) = (/0d0/)
        FT8OPX(1:NTRACE) = (/0d0/)
#elif (defined TRACERS_AEROSOLS_VBS)
        NTRACE=7
        FS8OPX(1:NTRACE) = (/0d0, 0d0, 1d0, 0d0, 0d0,  1d0 , 1d0/)
        FT8OPX(1:NTRACE) = (/0d0, 0d0, 1d0, 0d0, 0d0, 1.3d0, 1d0/)
#else
c       FS8OPX = (/0d0, 0d0, 1d0, 0d0, 2d0, 2d0,  1d0 , 1d0/)
        FS8OPX = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0,  1d0 , 1d0/)
        FT8OPX = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1.3d0, 1d0/)
#endif
      end if
#ifndef TRACERS_NITRATE
#ifdef SULF_ONLY_AEROSOLS
      NTRACE=1
      TRRDRY(1:NTRACE)=(/.15d0/)
      ITR(1:NTRACE) = (/1/)
      KRHTRA(1:NTRACE)=(/1/)
      NTRIX(1:NTRACE)=(/n_SO4/)
#else /* SULF aerosol allowed */
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_AEROSOLS_VBS
      NTRACE=7
      TRRDRY(1:NTRACE)=
     &(/.15d0, .44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=
     &(/1.d0,   1.d0, 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=(
     &/n_sO4,n_seasalt1,n_seasalt2,n_vbsAm2,n_isopp1a,n_BCIA,n_BCB/)
#else
      NTRACE=8
      TRRDRY(1:NTRACE)=
     &(/.15d0, .44d0, 1.7d0, .2d0, .2d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=
     &(/1.d0,   1.d0, 1.d0 , 1.d0, 1.d0, 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=(
     &/n_sO4,n_seasalt1,n_seasalt2,n_OCIA,n_OCB,n_isopp1a,n_BCIA,n_BCB/)
#endif /* TRACERS_AEROSOLS_VBS */
#else /* OFF: TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_VBS
      NTRACE=6
      TRRDRY(1:NTRACE)=(/.15d0, .44d0, 1.7d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=(/1.d0,   1.d0, 1.d0 , 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     *(/ n_sO4, n_seasalt1, n_seasalt2, n_vbsAm2, n_BCIA, n_BCB/)
#else
      NTRACE=7
      TRRDRY(1:NTRACE)=(/.15d0, .44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0/)
c augment BC by 50%
      FSTASC(1:NTRACE)=(/1.d0,   1.d0, 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4,4, 5,6/)
      KRHTRA(1:NTRACE)=(/1,1,1,1,1, 0,0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     *(/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_BCIA, n_BCB/)
#endif /* TRACERS_AEROSOLS_VBS */
#endif /* TRACERS_AEROSOLS_SOA */
#endif /* NOT defined: SULF_ONLY_AEROSOLS */
C**** define weighting (only used for clays so far)
      WTTR(1:NTRACE) = 1d0

#else /* TRACERS_NITRATE ON */

#ifdef SULF_ONLY_AEROSOLS
      call stop_model('SULF_ONLY_AEROSOLS and TRACERS_NITRATE on',255)
#else /* SULF aerosol allowed */
#ifdef TRACERS_AEROSOLS_SOA
#ifdef TRACERS_AEROSOLS_VBS
      NTRACE=8
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_vbsAm2, n_isopp1a,
     *   n_BCIA, n_BCB, n_NO3p/)
#else
      NTRACE=9
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .2d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
      FSTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     * (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2, 4,4,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1, 1,1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_isopp1a,
     *   n_BCIA, n_BCB, n_NO3p/)
#endif /* TRACERS_AEROSOLS_VBS */
#else /* OFF: TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_VBS
      NTRACE=7
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
      FSTASC(1:NTRACE)=
     *        (/1.d0,   1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     *        (/1.d0,   1.d0 , 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_vbsAm2, n_BCIA, n_BCB,
     *   n_NO3p/)
#else
      NTRACE=8
      if (rad_interact_aer > 0) then ! turn off default nitrate
        FS8OPX(3) = 0. ; FT8OPX(3) = 0.
      end if
      TRRDRY(1:NTRACE)=
     * (/.15d0,.44d0, 1.7d0, .2d0, .2d0, .08d0, .08d0,0.15d0/)
c augment BC by 50% (solar)
      FSTASC(1:NTRACE)=
     *        (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.5d0, 1.5d0,1.0d0/)
      FTTASC(1:NTRACE)=
     *        (/1.d0,   1.d0 , 1.d0 , 1.d0, 1.d0, 1.0d0, 1.0d0,1.0d0/)
cc tracer 1 is sulfate, tracers 2 and 3 are seasalt
      ITR(1:NTRACE) = (/ 1,2,2,4,4, 5,6,3/)
      KRHTRA(1:NTRACE)=(/1,1,1,1,1, 0,0,1/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics
      NTRIX(1:NTRACE)=
     * (/ n_sO4, n_seasalt1, n_seasalt2, n_OCIA, n_OCB, n_BCIA, n_BCB,
     *   n_NO3p/)
#endif /* TRACERS_AEROSOLS_VBS */
#endif /* TRACERS_AEROSOLS_SOA */
C**** define weighting (only used for clays so far)
      WTTR(1:NTRACE) = 1d0
#endif /* OFF: SULF_ONLY_AEROSOLS */
#endif /* TRACERS_NITRATE ON */
#endif /* TRACERS_AEROSOLS_Koch ON */

#ifdef TRACERS_DUST
C**** add dust optionally to radiatively active aerosol tracers
C**** should also work if other aerosols are not used
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      nrad_clay=n1
      NTRACE=NTRACE+ntm_dust+3  ! add dust tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = 7
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for dust

      SELECT CASE (ntm_dust)
      CASE (4)
C**** effective radii for dust
        TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,1.386D0,
     &       2.773D0,5.545D0/)
C**** Particle density of dust
        TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,
     &       2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
        NTRIX(n1:NTRACE)=(/n_clay,n_clay,n_clay,n_clay,n_silt1,n_silt2,
     &       n_silt3/)
C**** define weighting for different clays
        WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,1D0,1D0,1D0/)
      CASE (5)
        TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,1.386D0,
     &       2.773D0,5.545D0,11.090D0/)
        TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,
     &       2.65D0,2.65D0/)
        NTRIX(n1:NTRACE)=(/n_clay,n_clay,n_clay,n_clay,n_silt1,n_silt2,
     &       n_silt3,n_silt4/)
        WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,1D0,1D0,1D0,
     &       1D0/)
      END SELECT
#else
#ifdef TRACERS_MINERALS
C**** add minerals optionally to radiatively active aerosol tracers
C**** so far all minerals have the properties of far traveled Saharan
C**** dust - to be changed soon
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      NTRACE = NTRACE + ntm_minerals ! add mineral tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = (/7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
     &     7,7,7,7,7,7,7,7,7,7,7,7,7,7/)
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for minerals
C**** effective radii for minerals
      TRRDRY(n1:NTRACE)=(/0.132D0,0.23D0,0.416D0,0.766D0,0.132D0,0.23D0,
     &     0.416D0,0.766D0,0.132D0,0.23D0,0.416D0,0.766D0,0.132D0,
     &     0.23D0,0.416D0,0.766D0,0.132D0,0.23D0,0.416D0,0.766D0,
     &     1.386D0,1.386D0,1.386D0,1.386D0,1.386D0,2.773D0,2.773D0,
     &     2.773D0,2.773D0,2.773D0,5.545D0,5.545D0,5.545D0,5.545D0,
     &     5.545D0/)
C**** Particle density of dust
      TRADEN(n1:NTRACE)=(/2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,
     &     2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,2.5D0,
     &     2.5D0,2.5D0,2.5D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,
     &     2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,2.65D0,
     &     2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
      NTRIX(n1:NTRACE)=(/n_clayilli,n_clayilli,n_clayilli,n_clayilli,
     &     n_claykaol,n_claykaol,n_claykaol,n_claykaol,n_claysmec,
     &     n_claysmec,n_claysmec,n_claysmec,n_claycalc,n_claycalc,
     &     n_claycalc,n_claycalc,n_clayquar,n_clayquar,n_clayquar,
     &     n_clayquar,n_sil1quar,n_sil1feld,n_sil1calc,n_sil1hema,
     &     n_sil1gyps,n_sil2quar,n_sil2feld,n_sil2calc,n_sil2hema,
     &     n_sil2gyps,n_sil3quar,n_sil3feld,n_sil3calc,n_sil3hema,
     &     n_sil3gyps/)
C**** define weighting for different clays
      WTTR(n1:NTRACE)=(/0.009D0,0.081D0,0.234D0,0.676D0,0.009D0,0.081D0,
     &     0.234D0,0.676D0,0.009D0,0.081D0,0.234D0,0.676D0,0.009D0,
     &     0.081D0,0.234D0,0.676D0,0.009D0,0.081D0,0.234D0,0.676D0,1D0,
     &     1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0,1D0/)
#endif
#ifdef TRACERS_QUARZHEM
C**** add quartz/hematite to radiatively active aerosol tracers
C**** so far all minerals have the properties of far traveled Saharan
C**** dust - to be changed soon
      if (rad_interact_aer > 0) then ! turn off default dust
        FS8OPX(7) = 0. ; FT8OPX(7) = 0.
      end if
      n1=NTRACE+1
      NTRACE = NTRACE + ntm_quarzhem ! add quartz/hematite aggregate tracers
c tracer 7 is dust
      ITR(n1:NTRACE) = (/7,7,7/)
      KRHTRA(n1:NTRACE)= 0.  ! no deliq for quartz/hematite aggregates
C**** effective radii for quartz/hematite
      TRRDRY(n1:NTRACE)=(/1.386D0,2.773D0,5.545D0/)
C**** Particle density of quartz/hematite
      TRADEN(n1:NTRACE)=(/2.65D0,2.65D0,2.65D0/)
C**** Define indices to map model tracer arrays to radiation arrays
C**** for the diagnostics. Adjust if number of dust tracers changes.
      NTRIX(n1:NTRACE)=(/n_sil1quhe,n_sil2quhe,n_sil3quhe/)
C**** define weighting
      WTTR(n1:NTRACE)=(/1.D0,1.D0,1.D0/)
#endif
#endif

#ifdef TRACERS_SPECIAL_Shindell
      if(NTRACE > maxNtraceFastj)
     &call stop_model("NTRACE > maxNtraceFastj in init_Rad",13)
#endif

      if (ktrend.ne.0) then
C****   Read in time history of well-mixed greenhouse gases
        call openunit('GHG',iu,.false.,.true.)
        call ghghst(iu)
        call closeunit(iu)
        if(file_exists('dH2O').and.H2ObyCH4.ne.0..and.Kradia.le.0) then
C****     Read in dH2O: H2O prod.rate in kg/m^2 per day and ppm_CH4
          call openunit('dH2O',iu,.false.,.true.)
#if defined(CUBED_SPHERE)
          call read_qma(iu,plbx)
#else
          call getqma(iu,lat_dg,plbx,dh2o,lm,jm)
#endif
          call closeunit(iu)
        else
          H2ObyCH4 = 0.
        end if
      end if
      if(dalbsnX.ne.0.) then
        call updBCd(1990) ; depoBC_1990 = depoBC
      endif
C**** set up unit numbers for 14 more radiation input files
      donotread = -9999
      nrfun(:) = 0 ! green light
      nrfun(12:13) = donotread ! not used in GCM
      nrfun(10:11) = donotread ! obsolete O3 data
      nrfun(6)     = donotread ! dust read externally now
      if(.not.transmission_corrections) nrfun(4) = donotread
      if(madvol == 0) nrfun(7) = donotread
      if(madeps == 0) nrfun(8) = donotread
      if(ksolar < 0)  nrfun(9) = donotread
      DO IU=1,14
        if(nrfun(iu) == donotread) cycle
        call openunit(RUNSTR(IU),NRFUN(IU),QBIN(IU),.true.)
      END DO

#ifdef CHL_from_SeaWIFs
C**** open chlorophyll data
      call openunit("CHL_DATA",iu_CHL,.true.,.true.)
#endif

      LS1_loc=1  ! default
C***********************************************************************
C     Main Radiative Initializations
C     ------------------------------------------------------------------
      CALL RCOMP1 (NRFUN)
      if (am_i_root()) CALL WRITER(6,0)  ! print rad. control parameters
C***********************************************************************
      DO IU=1,14
        if(nrfun(iu) == donotread) cycle
        call closeunit(NRFUN(IU))
      END DO
C**** Save initial (currently permanent and global) Q in rad.layers
      do LR=1,LM_REQ
        shl0(LR) = shl(LM+LR)
      end do
      write(out_line,*) 'spec.hum in rad.equ.layers:',shl0
      call write_parallel(trim(out_line),unit=6)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FS8OPX_orig(:)=FS8OPX(:); FT8OPX_orig(:)=FT8OPX(:) ! aerosols

C**** Read in the factors used for alterations:
      call openunit('ALT_GHG_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,13D8.3)') skip,(FULGAS_lat(nn,n),nn=1,13)
      enddo
      call closeunit(iu2)
      call openunit('ALT_AER_LAT',iu2,.false.,.true.)
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FS8OPX_lat(nn,n),nn=1,8)
      enddo
      read(iu2,*) ! skip first line
      do n=1,46
        read(iu2,'(a6,8D8.3)') skip,(FT8OPX_lat(nn,n),nn=1,8)
      enddo
      call closeunit(iu2)
#endif

#if (defined OBIO_RAD_coupling) || (defined CHL_from_SeaWIFs)
      call openunit('cfle1',iu_bio,.false.,.true.)
      do ic = 1,6
        read(iu_bio,'(a50)')title
      enddo
      do nl = 1,nlt
        read(iu_bio,20) lambda,saw,sbw
        lam(nl) = lambda
        aw(nl) = saw
        bw(nl) = sbw
        if (lam(nl) .lt. 900) then
          expterm = exp(-(aw(nl)+0.5*bw(nl)))
          tlog = dlog(1.0D-36+expterm)
          fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
          wfac(nl) = max(0d0,min(fac,1d0))
        else
          rlam = float(lam(nl))
          fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
          wfac(nl) = max(fac,0d0)
        endif
      enddo
      print*,'RAD_DRV, wfac initializ= ', wfac
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)
#endif

#ifdef TRACERS_ON
c**** set tracerRadiaActiveFlag for radiatively active tracer
      do n=1,ntrace
        if (ntrix(n) > 0) tracerRadiaActiveFlag(ntrix(n))=.true.
      end do
      nTracerRadiaActive=count(tracerRadiaActiveFlag)
#endif

c transplanted from main().  needs reviving
c      USE RAD_COM, only : dimrad_sv
c      CHARACTER aDATE*14
c      if (Kradia.ne.0 .and. Kradia<10) then
c        write(aDATE(1:7),'(a3,I4.4)') aMON(1:3),Jyear
c        if (Kradia.gt.0) aDATE(4:7)='    '
c        call openunit(trim('RAD'//aDATE(1:7)),iu_RAD,.true.,.false.)
c        if (Kradia.lt.0) call io_POS(iu_RAD,Itime-1,2*dimrad_sv,Nrad)
c      end if

      RETURN
      END SUBROUTINE init_RAD

      subroutine SETATM             ! dummy routine in gcm
      end subroutine SETATM

      subroutine GETVEG(LONR,LATR)  ! dummy routine in gcm
      integer LONR,LATR
      end subroutine GETVEG

      SUBROUTINE daily_RAD(end_of_day)
!@sum  daily_RAD sets radiation parameters that change every day
!@auth G. Schmidt
!@calls RADPAR:RCOMPT
      USE CONSTANT, only : by12
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : am_I_root,GRID,REWIND_PARALLEL
     *     ,READT_PARALLEL, getDomainBounds
      USE RESOLUTION, only : im,jm
      use model_com, only: modelEclock, calendar
      USE GEOM, only : imaxj
      USE RADPAR, only : FULGAS,JYEARR=>JYEAR,JDAYR=>JDAY
     *     ,xref,KYEARV
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_orig
#endif
      USE RADPAR, only : rcompt,writet
      USE RAD_COM, only : co2x,n2ox,ch4x,cfc11x,cfc12x,xGHGx,h2ostratx
     *     ,o3x,o3_yr,ghg_yr,co2ppm,Volc_yr,albsn_yr,dalbsnX
#ifdef CHL_from_SeaWIFs
     *     ,iu_CHL,achl,echl1,echl0,bchl,cchl
      USE FLUXES, only : focean,atmocn
#endif
      use DIAG_COM, only : iwrite,jwrite,itwrite
      use CalendarMonth_mod
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
!@var TEMP_LOCAL stores ACHL+ECHL1 to avoid the use of common block
      REAL*8 :: TEMP_LOCAL(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
!@var IMON0 current month for CHL climatology reading
      INTEGER, SAVE :: IMON0 = 0
      INTEGER :: LSTMON,I,J
      REAL*8 TIME

      INTEGER :: J_0,J_1, I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      integer :: year, month, dayOfYear, date
      type (CalendarMonth) :: cMonth

      call modelEclock%get(year=year, month=month,
     &     dayOfYear=dayOfYear, date=date)

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1,
     &         HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Update time dependent radiative parameters each day
!     Get black carbon deposition data for the appropriate year
!     (does nothing except at a restart or the beginning of a new year)
      if(dalbsnX.ne.0.) then
        if (albsn_yr.eq.0) then
          call updBCd (year)
        else
          call updBCd (albsn_yr)
        end if
      endif
!     Hack: 2 specific volc. eruption scenarios for 2000-2100 period
      if(volc_yr.eq.-2010) then              ! repeat some old volcanos
         KYEARV=YEAR
         if(YEAR.GT.2010) KYEARV=YEAR-100  ! go back 100 years
      end if
      if(volc_yr.eq.-2000) then
         KYEARV=YEAR
         if(YEAR.GT.2000) KYEARV=YEAR-50   ! go back 50 years til 2050
         if(YEAR.GT.2050) KYEARV=YEAR-150  ! then go back 150 years
      end if

      JDAYR=dayOfYear
      JYEARR=YEAR
      CALL RCOMPT
!     FULGAS(2:) is set only in the first call to RCOMPT unless ghg_yr=0
!     Optional scaling of the observed value only in case it was (re)set
      if(ghg_yr.eq.0 .or. .not. end_of_day) then
         FULGAS(2)=FULGAS(2)*CO2X
         FULGAS(6)=FULGAS(6)*N2OX
         FULGAS(7)=FULGAS(7)*CH4X
         FULGAS(8)=FULGAS(8)*CFC11X
         FULGAS(9)=FULGAS(9)*CFC12X
         FULGAS(11)=FULGAS(11)*XGHGX
      end if
      IF(.not. end_of_day .and. H2OstratX.GE.0.)
     *     FULGAS(1)=FULGAS(1)*H2OstratX
      IF(.not. end_of_day .or. O3_yr==0.) FULGAS(3)=FULGAS(3)*O3X

C**** write trend table for forcing 'itwrite' for years iwrite->jwrite
C**** itwrite: 1-2=GHG 3=So 4-5=O3 6-9=aerosols: Trop,DesDust,Volc,Total
      if (am_i_root() .and.
     *  jwrite.gt.1500) call writet (6,itwrite,iwrite,jwrite,1,0)

#ifdef ALTER_RADF_BY_LAT
C**** Save initial rad forcing alterations:
      FULGAS_orig(:)=FULGAS(:) ! GHGs
#endif

C**** Define CO2 (ppm) for rest of model
      co2ppm = FULGAS(2)*XREF(1)

#ifdef CHL_from_SeaWIFs
C**** Read in Seawifs files here
      IF (month.NE.IMON0) THEN
      IF (IMON0==0) THEN
C**** READ IN LAST MONTH'S END-OF-MONTH DATA
        LSTMON=month-1
        if (lstmon.eq.0) lstmon = 12
        CALL READT_PARALLEL
     *       (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,LSTMON)
        ECHL0 = TEMP_LOCAL(:,:,2)
      ELSE
C**** COPY END-OF-OLD-MONTH DATA TO START-OF-NEW-MONTH DATA
        ECHL0=ECHL1
      END IF
C**** READ IN CURRENT MONTHS DATA: MEAN AND END-OF-MONTH
      IMON0=month
      if (month.eq.1) CALL REWIND_PARALLEL( iu_CHL )
      CALL READT_PARALLEL
     *     (grid,iu_CHL,NAMEUNIT(iu_CHL),TEMP_LOCAL,1)
      ACHL  = TEMP_LOCAL(:,:,1)
      ECHL1 = TEMP_LOCAL(:,:,2)

C**** FIND INTERPOLATION COEFFICIENTS (LINEAR/QUADRATIC FIT)
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          BCHL(I,J)=ECHL1(I,J)-ECHL0(I,J)
          CCHL(I,J)=3.*(ECHL1(I,J)+ECHL0(I,J))-6.*ACHL(I,J)
        END DO
      END DO
      END IF
C**** Calculate CHL for current day
      cMonth = calendar%getCalendarMonth(month, year)
      TIME=(DATE-.5)/cMonth%daysinMonth -.5 ! -.5<TIME<.5
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          IF (FOCEAN(I,J).gt.0) THEN
C**** CHL always uses quadratic fit
            atmocn%CHL(I,J)=ACHL(I,J)+BCHL(I,J)*TIME
     *           +CCHL(I,J)*(TIME**2-BY12)
            if (atmocn%CHL(I,J).lt. 0) atmocn%CHL(I,J)=0. ! just in case
          END IF
        END DO
      END DO
C**** REPLICATE VALUES AT POLE
      IF(HAVE_NORTH_POLE) then
       if (FOCEAN(1,JM).gt.0) atmocn%CHL(2:IM,JM)=atmocn%CHL(1,JM)
      ENDIF
      IF(HAVE_SOUTH_POLE) then
       if (FOCEAN(1, 1).gt.0) atmocn%CHL(2:IM, 1)=atmocn%CHL(1, 1)
      ENDIF
#endif

      RETURN
      END SUBROUTINE daily_RAD

      SUBROUTINE DAILY_orbit(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@calls constant:orbit
      use model_com, only: modelEclock
      USE RAD_COM, only : RSDIST,COSD,SIND,COSZ_day,SUNSET,
     *     variable_orb_par,orb_par_year_bp, useOrbit => orbit
      USE DOMAIN_DECOMP_ATM, only : am_I_root
      use RAD_COSZ0, only : daily_cosz
      use BaseTime_mod
      use TimeInterval_mod
      use Rational_mod
      IMPLICIT NONE
      REAL*8 :: SUNLON,SUNLAT,LAM,EDPY,VEDAY,PYEAR
      LOGICAL, INTENT(IN) :: end_of_day
      integer :: year, dayOfYear
      type (BaseTime) :: t
      real*8 :: declinationAngle
      type (TimeInterval) :: halfDay

      call modelEclock%get(year=year, dayOfYear=dayOfYear)

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
C**** This is for noon (GMT) for new day.

C**** The orbital calculation will need to vary depending on the kind
C**** of calendar adopted (i.e. a generic 365 day year, or a transient
C**** calendar including leap years etc.).  For transient calendars the
C**** dayOfYear passed to orbit needs to be adjusted to represent the number
C**** of days from Jan 1 2000AD.
c      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
c      dayOfYear => dayOfYear + 365 * (YEAR-2000) + appropriate number of leaps
C**** Default calculation (no leap, VE=Mar 21 hr 0)
c      EDPY=365d0 ; VEDAY=79d0           ! Generic year
C**** PMIP calculation (no leap, VE=Mar 21 hr 12)
      EDPY=365d0 ; VEDAY=79.5d0           ! Generic year
C**** Update orbital parameters at start of year
      if (variable_orb_par == 1.and.dayOfYear == 1) then
        pyear = YEAR - orb_par_year_bp ! bp=before present model year
        call useOrbit%setYear(pYear)
        if (am_I_root()) then
          write(6,*) 'Set orbital parameters for year ',pyear,' (CE)'
          if (orb_par_year_bp.ne.0) write(6,*) 'offset by',
     *      orb_par_year_bp,' years from model year'
          write(6,*) "   Eccentricity: ", useOrbit%getEccentricity()
          write(6,*) "   Obliquity (degs): ",useOrbit%getObliquity()
          write(6,*) "   Precession (degs from ve): ",
     *         useOrbit%getLongitudeAtPeriapsis()
        end if
      end if

      ! Use time for the _middle_ of the day to compute 
      ! zenith angle:

      halfDay = TimeInterval(useOrbit%getMeanDay() / 2)
      t = newBaseTime(modelEClock%getCurrentTime() + halfDay)

      sinD = useOrbit%getSinDeclinationAngle(t)
      cosD = sqrt(1-sinD**2)
      rsdist = useOrbit%getDistance(t)**2

      call daily_cosz(sind,cosd,cosz_day,sunset)

      RETURN
      END SUBROUTINE DAILY_orbit

      SUBROUTINE DAILY_ch4ox(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@vers 2013/03/27
!@auth Original Development Team
!@calls constant:orbit
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : q
      use model_com, only: modelEclock
      USE MODEL_COM, only : itime
      USE GEOM, only : axyp,imaxj,lat2d
      USE ATM_COM, only : byMA
      USE RADPAR, only : ghgam,ghgyr2,ghgyr1
      USE RAD_COM, only : dh2o,H2ObyCH4,ghg_yr
#ifdef TRACERS_WATER
      use OldTracer_mod, only: tr_wd_type, nWATER,tr_H2ObyCH4, itime_tr0
      USE TRACER_COM, only: trm,NTM
#endif
      USE DIAG_COM, only : aj=>aj_loc,j_h2och4,ftype,ntype
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, am_I_root
      IMPLICIT NONE
      REAL*8 :: xCH4,xdH2O
      INTEGER i,j,l,iy,it
      LOGICAL, INTENT(IN) :: end_of_day
#ifdef TRACERS_WATER
      INTEGER n
#endif
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      integer :: year, month

      call modelEclock%get(year=year, month=month)

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      IF (.not.end_of_day) RETURN

C**** Tasks to be done at end of day only
      if (H2ObyCH4.gt.0) then
C****   Add obs. H2O generated by CH4(*H2ObyCH4) using a 2 year lag
        iy = year - 2 - ghgyr1 + 1
        if (ghg_yr.gt.0) iy = ghg_yr - 2 - ghgyr1 + 1
        if (iy.lt.1) iy=1
        if (iy.gt.ghgyr2-ghgyr1+1) iy=ghgyr2-ghgyr1+1
        xCH4=ghgam(3,iy)*H2ObyCH4
c        If (AM_I_ROOT())
c     &    write(6,*) 'add in stratosphere: H2O gen. by CH4(ppm)=',xCH4

        do l=1,lm
        do j=J_0,J_1
        do i=I_0,imaxj(j)
#ifdef CUBED_SPHERE
          call lat_interp_qma(lat2d(i,j),l,month,xdH2O)
#else
          xdH2O = dH2O(j,l,month)
#endif
          q(i,j,l)=q(i,j,l)+xCH4*xdH2O*byMA(l,i,j)
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
          do n=1,ntm
            if (itime_tr0(n).le.itime) then
              select case (tr_wd_type(n))
              case (nWater)    ! water: add CH4-sourced water to tracers
                trm(i,j,l,n) = trm(i,j,l,n) +
     +                tr_H2ObyCH4(n)*xCH4*xdH2O*axyp(i,j)
              end select
            end if
          end do
#endif
          do it=1,ntype
            call inc_aj(i,j,it,j_h2och4,xCH4*xdH2O*ftype(it,i,j))
          end do
        end do
        end do
        If (HAVE_NORTH_POLE) q(2:im,jm,l)=q(1,jm,l)
        If (HAVE_SOUTH_POLE) q(2:im, 1,l)=q(1, 1,l)
#ifdef TRACERS_WATER
        do n=1,ntm
          If (HAVE_SOUTH_POLE) trm(2:im, 1,l,n)=trm(1, 1,l,n)
          If (HAVE_NORTH_POLE) trm(2:im,jm,l,n)=trm(1,jm,l,n)
        end do
#endif
        end do
      end if

      RETURN
      END SUBROUTINE DAILY_ch4ox

      SUBROUTINE RADIA
!@sum  RADIA adds the radiation heating to the temperatures
!@vers 2013/03/27
!@auth Original Development Team
!@calls tropwmo,coszs,coszt, RADPAR:rcompx ! writer,writet
      USE CONSTANT, only : lhe,lhs,twopi,tf,stbo,rhow,mair,grav
     *     ,bysha,pi,radian,areag
      USE RESOLUTION, only : pmtop
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : kradia,lm_req,p,t,q,iu_rad,req_fac_d
      USE MODEL_COM
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE ATM_COM, only : byaml00
      USE GEOM, only : imaxj, axyp, byaxyp
     &     ,lat2d,lon2d
      USE RADPAR
     &  , only :  ! routines
     &           lx  ! for threadprivate copyin common block
     &          ,tauwc0,tauic0 ! set in radpar block data
     &          ,writer,rcompx,updghg
C     INPUT DATA         ! not (i,j) dependent
     X          ,S00WM2,RATLS0,S0,JYEARR=>JYEAR,JDAYR=>JDAY,FULGAS
     &          ,use_tracer_chem,FS8OPX,FT8OPX,use_o3_ref,KYEARG,KJDAYG
     &          ,planck_tmin,planck_tmax
#ifdef ALTER_RADF_BY_LAT
     &          ,FS8OPX_orig,FT8OPX_orig,FULGAS_orig
#endif
C     INPUT DATA  (i,j) dependent
     &             ,JLAT46=>JLAT,ILON72=>ILON,JGCM,IGCM
     &             ,L1,LMR=>NL, PLB ,TLB,TLM ,SHL,RHL
     &             ,ltopcl,TAUWC ,TAUIC ,SIZEWC ,SIZEIC, kdeliq
     &             ,POCEAN,PEARTH,POICE,PLICE,PLAKE,COSZ,PVT
     &             ,TGO,TGE,TGOI,TGLI,TSL,WMAG,WEARTH
     &             ,AGESN,SNOWE,SNOWOI,SNOWLI,dALBsn, ZSNWOI,ZOICE
     &             ,zmp,fmp,flags,LS1_loc,snow_frac,zlake
     *             ,TRACER,NTRACE,FSTOPX,FTTOPX,chem_IN
     *             ,FTAUC,LOC_CHL,FSTASC,FTTASC
#ifdef HEALY_LM_DIAGS
     *             ,VTAULAT
#endif

C     OUTPUT DATA
     &          ,TRDFLB ,TRNFLB ,TRUFLB, TRFCRL ,chem_out
     &          ,SRDFLB ,SRNFLB ,SRUFLB, SRFHRL
     &          ,PLAVIS ,PLANIR ,ALBVIS ,ALBNIR ,FSRNFG
     &          ,SRRVIS ,SRRNIR ,SRAVIS ,SRANIR ,SRXVIS ,SRDVIS
     &          ,BTEMPW ,TTAUSV ,SRAEXT ,SRASCT ,SRAGCB
     &          ,SRDEXT ,SRDSCT ,SRDGCB ,SRVEXT ,SRVSCT ,SRVGCB
     &          ,aesqex,aesqsc,aesqcb
     &          ,SRXNIR,SRDNIR
      USE RAD_COM, only : modrd,nrad
      USE RAD_COM, only : rqt,srhr,trhr,fsf,cosz1,s0x,rsdist,nradfrc
     *     ,CH4X_RADoverCHEM
     *     ,plb0,shl0,tchg,alb,fsrdir,srvissurf,srdn,cfrac,rcld
     *     ,chem_tracer_save,rad_interact_aer,kliq,RHfix,CLDx
     *     ,ghg_yr,CO2X,N2OX,CH4X,CFC11X,CFC12X,XGHGX,rad_forc_lev,ntrix
     *     ,wttr,cloud_rad_forc,CC_cdncx,OD_cdncx,cdncl,nrad_clay
     *     ,dALBsnX,rad_to_chem,trsurf,dirvis
     *     ,FSRDIF,DIRNIR,DIFNIR,aer_rad_forc,clim_interact_chem
     *     ,TAUSUMW,TAUSUMI
#ifdef mjo_subdd
     *     ,SWHR,LWHR,SWHR_cnt,LWHR_cnt,OLR_acc,OLR_cnt
     *     ,swu_avg,swu_cnt
#endif
#ifdef ALTER_RADF_BY_LAT
     *     ,FULGAS_lat,FS8OPX_lat,FT8OPX_lat
#endif
#ifdef TRACERS_DUST
     &     ,srnflb_save,trnflb_save
#endif
#ifdef TRACERS_ON
     &     ,ttausv_sum,ttausv_sum_cs,ttausv_count,ttausv_save
     &     ,ttausv_cs_save,aerAbs6SaveInst
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ,ttausv_ntrace
#endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &     ,stratO3_tracer_save
#endif
      USE RANDOM
      USE CLOUDS_COM, only : tauss,taumc,svlhx,rhsav,svlat,cldsav,
     *     cldmc,cldss,csizmc,csizss,llow,lmid,lhi,fss
#ifdef SCM
      USE SCMCOM, only : SCM_SURF_ALBEDO_FLAG,iu_scm_prt
      USE SCMDIAG, only : SRDFLBTOP,SRNFLBTOP,SRUFLBTOP,TRUFLBTOP,
     &                    SRDFLBBOT,SRNFLBBOT,SRUFLBBOT,TRUFLBBOT,
     &                    TRDFLBBOT,TRDFLBTOP,SRFHRLCOL,TRFCRLCOL,
     &                    CSSRNTOP,CSTRUTOP,CSSRNBOT,CSTRNBOT,
     &                    CSSRDBOT,TRNFLBBOT,dTradlw,dTradsw
#endif
      USE DIAG_COM, only : ia_rad,jreg,aij=>aij_loc,aijl=>aijl_loc
     &     ,ntype,ftype,itocean,itlake,itearth,itlandi,itoice,itlkice
     *     ,adiurn=>adiurn_loc,ndiuvar,ia_rad_frc,
#ifndef NO_HDIURN
     *     hdiurn=>hdiurn_loc,
#endif
     *     iwrite,jwrite,itwrite,ndiupt,j_pcldss,j_pcldmc,ij_pmccld,
     *     j_clddep,j_pcld,ij_cldcv,ij_pcldl,ij_pcldm,ij_pcldh,
     *     ij_cldtppr,j_srincp0,j_srnfp0,j_srnfp1,j_srincg,
     *     j_srnfg,j_brtemp,j_trincg,j_hsurf,j_hatm,ij_trnfp0,
     *     j_plavis, j_planir, j_albvis, j_albnir,
     *     j_srrvis, j_srrnir, j_sravis, j_sranir,
     *     ij_srnfp0,ij_srincp0,ij_srnfg,ij_srincg,ij_btmpw,ij_srref
     *     ,ij_srvis,ij_rnfp1,j_clrtoa,j_clrtrp,j_tottrp,ijl_rc
     *     ,ijdd,idd_cl7,idd_ccv,idd_isw,idd_palb,idd_galb, idd_aot
     *     ,idd_aot2
     *     ,idd_absa,jl_srhr,jl_trcr,jl_totcld,jl_sscld,jl_mccld
     *     ,ij_frmp,jl_wcld,jl_icld,jl_wcod,jl_icod,jl_wcsiz,jl_icsiz
     *     ,jl_wcldwt,jl_icldwt
     *     ,ij_clr_srincg,ij_CLDTPT,ij_cldt1t,ij_cldt1p,ij_cldcv1
     *     ,ij_wtrcld,ij_icecld,ij_optdw,ij_optdi,ij_swcrf,ij_lwcrf
     *     ,AFLX_ST, hr_in_day,hr_in_month,ij_srntp,ij_trntp
     *     ,ij_clr_srntp,ij_clr_trntp,ij_clr_srnfg,ij_clr_trdng
     *     ,ij_clr_sruptoa,ij_clr_truptoa,ijl_cf
     *     ,ij_swdcls,ij_swncls,ij_lwdcls,ij_swnclt,ij_lwnclt, NREG
     *     ,adiurn_dust,j_trnfp0,j_trnfp1,ij_srvdir, ij_srvissurf
     *     ,ij_chl, ij_swaerrf, ij_lwaerrf,ij_swaersrf,ij_lwaersrf
     *     ,ij_swaerrfnt,ij_lwaerrfnt,ij_swaersrfnt,ij_lwaersrfnt
     *     ,ij_swcrf2,ij_lwcrf2, ij_siswd, ij_siswu, save3dAOD
#ifdef ACCMIP_LIKE_DIAGS
     *     ,ij_fcghg ! array
#endif
#ifdef HEALY_LM_DIAGS
     *     ,j_vtau,j_ghg
#endif

      USE ATM_COM, only : pk,pedn,pmid,pdsig,ltropo,MA,byMA
      USE SEAICE, only : rhos,ace1i,rhoi
      USE SEAICE_COM, only : si_atm
      USE GHY_COM, only : snowe_com=>snowe,snoage,fr_snow_rad_ij,fearth
#ifdef USE_ENT
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
     &                    ,N_COVERTYPES !YKIM-temp hack
      use ent_drv, only : map_ent2giss  !YKIM-temp hack
#else
      USE VEG_COM, only : vdata
#endif
      USE LAKES_COM, only : flake,mwl
      USE FLUXES, only : asflx4,atmocn,atmice,atmgla,atmlnd,atmsrf
     &     ,flice,fland,focean
      USE DOMAIN_DECOMP_ATM, ONLY: grid, write_parallel
      USE DOMAIN_DECOMP_ATM, ONLY: GLOBALSUM, getDomainBounds
      USE RAD_COSZ0, only : COSZT,COSZS

#ifdef TRACERS_ON
      use OldTracer_mod, only: trname, trpdens
      USE TRACER_COM, only: NTM
     *     ,n_Ox,trm,n_OCB,n_BCII,n_BCIA
     *     ,n_OCIA,N_OCII,n_so4_d2,n_so4_d3,n_SO4,n_stratOx
     *     ,n_N_AKK_1
#ifdef TRACERS_NITRATE
      use OldTracer_mod, only: tr_mm
      use TRACER_COM, only: n_NH4,n_NO3p
#endif
#ifdef TRACERS_AEROSOLS_SOA
      use TRACER_COM, only: n_isopp1a,n_isopp2a
#ifdef TRACERS_TERP
      use TRACER_COM, only: n_apinp1a,n_apinp2a
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
      use TRACER_COM, only: n_ococean
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_AEROSOLS_Koch
c          use TRACER_COM, only: SNFST0,TNFST0
#endif  /* TRACERS_AEROSOLS_Koch */
#ifdef TRACERS_AEROSOLS_VBS
      use TRACERS_VBS, only: vbs_tr
#endif
      USE TRDIAG_COM, only: taijs=>taijs_loc,taijls=>taijls_loc,ijts_fc
     *     ,ijts_tau,ijts_tausub,ijts_fcsub,ijlt_3dtau,ijlt_3daaod
     *     ,ijts_sqex
     *     ,ijts_sqexsub,ijts_sqsc,ijts_sqscsub,ijts_sqcb,ijts_sqcbsub
     *     ,diag_rad
#ifdef AUXILIARY_OX_RADF
     *     ,ijts_auxfc
#endif /* AUXILIARY_OX_RADF */
#ifdef BC_ALB
     *     ,ijts_alb
#endif  /* BC_ALB */
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: Lmax_rad_O3,Lmax_rad_CH4
#endif /* TRACERS_SPECIAL_Shindell */
#endif /* TRACERS_ON */
#ifdef TRACERS_AMP
      USE AERO_CONFIG, only: nmodes
      USE AMP_AEROSOL, only: AMP_DIAG_FC
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only: icomp,TOMAS_DIAG_FC
      USE TRACER_COM, only : n_ANUM
#endif
      use AerParam_mod, only: dCDNC_est
      use AerParam_mod, only: depoBC,depoBC_1990
      USE TimerPackage_mod, only: startTimer => start, stopTimer => stop
      USE Dictionary_mod, only : get_param, is_set_param
      IMPLICIT NONE
C
C     INPUT DATA   partly (i,j) dependent, partly global
      REAL*8 U0GAS,taulim, xdalbs,sumda,tauda,fsnow
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     sumda_psum,tauda_psum

      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     COSZ2,COSZA,TRINCG,BTMPW,WSOIL,fmp_com
      REAL*8, DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                    grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFS,TNFS
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFSCRF,TNFSCRF,SNFSCRF2,TNFSCRF2
      REAL*8, DIMENSION(18,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     SNFSAERRF,TNFSAERRF
#ifdef ACCMIP_LIKE_DIAGS
!@var snfs_ghg,tnfs_ghg like SNFS/TNFS but with reference GHG for
!@+   radiative forcing calculations. TOA only.
!@+   index 1=CH4, 2=N2O, 3=CFC11, 4=CFC12.
      real*8,dimension(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                   snfs_ghg,tnfs_ghg
      real*8,dimension(4) :: sv_fulgas_ref,sv_fulgas_now
      integer :: nf,GFrefY,GFrefD,GFnowY,GFnowD
!@var nfghg fulgas( ) index of radf diag ghgs:
      integer, dimension(4) :: nfghg=(/7,6,8,9/)
#endif
#ifdef HEALY_LM_DIAGS
C  GHG Effective forcing relative to 1850
       real*8 :: ghg_totforc, CO2I=285.2,N2OI=.2754,CH4I=.791  ! 1850  GHG's
       real*8 ::              CO2R=337.9,N2OR=.3012,CH4R=1.547 ! RAD's 1979 Reference values
       real*8 ::              FCO2,FN2O,FCH4                   ! Current Model GHG
       real*8 :: Fe !! Function
#endif
#ifdef TRACERS_ON
!@var StauL sum of ttausv over L
      real*8 :: StauL
!@var SNFST,TNFST like SNFS/TNFS but with/without specific tracers for
!@+   radiative forcing calculations
      REAL*8,DIMENSION(2,NTRACE,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                          grid%J_STRT_HALO:grid%J_STOP_HALO)::
     *     SNFST,TNFST
!@var SNFST_o3ref,TNFST_o3ref like snfst,tnfst for special case ozone for
!@+   which ntrace fields are not defined. Indicies are :
!@+   1=LTROPO,reference, 2=TOA,reference; not saving surface forcing.
!@+   3=LTROPO,auxiliary, 4=TOA,auxiliary; not saving surface forcing.
      REAL*8,DIMENSION(4,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     SNFST_o3ref,TNFST_o3ref
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &    ,snfst_stratOx,tnfst_stratOx
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#ifdef BC_ALB
      REAL*8,DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                 grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     ALBNBC,NFSNBC
#endif /* BC_ALB */
#endif /* TRACERS_ON */
      REAL*8, DIMENSION(LM_REQ,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                         grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     TRHRS,SRHRS
      REAL*8, DIMENSION(0:LM+LM_REQ,
     *     grid%I_STRT_HALO:grid%I_STOP_HALO,
     *     grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *     TRHRA,SRHRA ! for adj.frc
      REAL*8, DIMENSION(LM) :: TOTCLD,dcc_cdncl,dod_cdncl
      INTEGER I,J,L,K,KR,LR,JR,IH,IHM,INCH,JK,IT,iy,iend,N,onoff_aer
     *     ,onoff_chem,LFRC,JTIME,n1,moddrf
      REAL*8 ROT1,ROT2,PLAND,CSS,CMC,DEPTH,QSS,TAUSSL,RANDSS
     *     ,TAUMCL,ELHX,CLDCV,X,OPNSKY,CSZ2,tauup,taudn,ptype4(4)
     *     ,taucl,wtlin,MSTRAT,STRATQ,STRJ,MSTJ,optdw,optdi,rsign_aer
     *     ,rsign_chem,tauex5,tauex6,tausct,taugcb,dcdnc
     *     ,QR(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &            grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,CLDinfo(LM,3,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8 tmpS(8),tmpT(8)
      REAL*8 RANDXX ! temporary
      REAL*8 QSAT
#ifdef BC_ALB
      REAL*8 dALBsn1
#endif
      LOGICAL NO_CLOUD_ABOVE, set_clayilli,set_claykaol,set_claysmec,
     &     set_claycalc,set_clayquar
C
      REAL*8  RDSS(LM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                grid%J_STRT_HALO:grid%J_STOP_HALO)
     *     ,RDMC(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &           grid%J_STRT_HALO:grid%J_STOP_HALO)

      REAL*8 :: TMP(NDIUVAR)
      INTEGER, PARAMETER :: NLOC_DIU_VAR = 8
      INTEGER :: idx(NLOC_DIU_VAR)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      INTEGER, PARAMETER :: NLOC_DIU_VARb = 5
#else
      INTEGER, PARAMETER :: NLOC_DIU_VARb = 3
#endif
      INTEGER :: idxb(NLOC_DIU_VARb)

      integer :: aj_alb_inds(8)
      real*8, dimension(lm_req) :: bydpreq

c     INTEGER ICKERR,JCKERR,KCKERR
      INTEGER :: J_0, J_1, I_0, I_1
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      character(len=300) :: out_line

      integer :: nij_before_j0,nij_after_j1,nij_after_i1
      integer :: initial_GHG_setup

#ifdef USE_ENT
      real*8 :: PVT0(N_COVERTYPES)
#endif
#ifdef TRACERS_NITRATE
      real*8 :: nh4_on_no3
#endif
#ifdef TRACERS_TOMAS
      real*8 :: qcb_col(6,ICOMP-2)
#endif

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,POND_MELT
      LOGICAL, DIMENSION(:,:), POINTER :: FLAG_DSWS
      integer :: year, dayOfYear, hour, date

      call modelEclock%get(year=year, dayOfYear=dayOfYear,
     *     hour=hour, date=date)

      RSI => SI_ATM%RSI
      MSI => SI_ATM%MSI
      SNOWI => SI_ATM%SNOWI
      POND_MELT => SI_ATM%POND_MELT
      FLAG_DSWS => SI_ATM%FLAG_DSWS

C
C****
      call startTimer('RADIA()')

      idx = (/ (IDD_CL7+i-1,i=1,7), IDD_CCV /)
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      idxb = (/ IDD_PALB, IDD_GALB, IDD_ABSA, idd_aot, idd_aot2 /)
#else
      idxb = (/ IDD_PALB, IDD_GALB, IDD_ABSA /)
#endif
      call getDomainBounds(grid, HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      J_0S = grid%J_STRT_SKP
      J_1S = grid%J_STOP_SKP


C****
C**** FLAND     LAND COVERAGE (1)
C**** FLICE     LAND ICE COVERAGE (1)
C****
C**** GTEMPR RADIATIVE TEMPERATURE ARRAY OVER ALL SURFACE TYPES (K)
C****   RSI  RATIO OF OCEAN ICE COVERAGE TO WATER COVERAGE (1)
C****
C**** VDATA  1-11 RATIOS FOR THE 11 VEGETATION TYPES (1)
C****

C**** limit optical cloud depth from below: taulim
      taulim=min(tauwc0,tauic0) ! currently both .001
      tauwc0 = taulim ; tauic0 = taulim
C**** Calculate mean cosine of zenith angle for the current physics step
      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
c      ROT2=ROT1+TWOPI*DTsrc/SECONDS_PER_DAY
c      CALL COSZT (ROT1,ROT2,COSZ1)
      CALL CALC_ZENITH_ANGLE   ! moved to main loop

      if (kradia.gt.0) then    ! read in all rad. input data (frc.runs)
        iend = 1
        it = itime-1           ! make sure, at least 1 record is read
        do while (mod(itime-it,NDAY*INT_DAYS_PER_YEAR).ne.0)
          read(iu_rad,end=10,err=10) it
C****   input data:          WARNINGS
C****        1 - any changes here also go in later (look for 'iu_rad')
C****        2 - keep "dimrad_sv" up-to-date:         dimrad_sv=IM*JM*{
     *     ,T,RQT,atmsrf%TsAvg                         ! LM+LM_REQ+1+
     *     ,QR,P,CLDinfo,rsi,msi                       ! LM+1+3*LM+1+1+
!     *     ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm)    ! (4+)
     *     ,wsoil,atmsrf%wsavg,snowi,atmgla%snow,snowe_com ! 1+1+1+1+1+
     *     ,snoage,fmp_com,flag_dsws,ltropo            ! 3+1+.5+.5+
     *     ,fr_snow_rad_ij,mwl,flake                   ! 2+1+1
C****   output data: really needed only if kradia=2
     *     ,srhra,trhra                                ! 2(LM+LM_REQ+1)}
C****   total: dimrad_sv= IM*JM*(7*LM + 3*LM_REQ + 24 (+4)) => RAD_COM.f
     *     ,iy
          if (qcheck) then
            write(out_line,*) 'reading RADfile at Itime',Itime,it,iy
            call write_parallel(trim(out_line),unit=6)
          endif
        end do
        iend = 0
   10   if (it.ne.iy.or.iend.eq.1) then
          write(out_line,*) 'RAD input file bad or too short:',itime,
     *                       it,iy,iend
          call write_parallel(trim(out_line),unit=6)
          call stop_model('RADIA: input file bad or too short',255)
        end if
C****   Find arrays derived from P : PEdn and PK (forcing experiments)
        call calc_ampk(lm)
      end if

      IF (MODRD.NE.0) GO TO 900
      IDACC(ia_rad)=IDACC(ia_rad)+1
      moddrf=1    ! skip rad.forcing diags if nradfrc.le.0
      if (nradfrc>0) moddrf=mod(itime-itimei,nrad*nradfrc)
C****
      if (moddrf==0) IDACC(ia_rad_frc)=IDACC(ia_rad_frc)+1
C**** Interface with radiation routines, done only every NRAD time steps
C****
C**** Calculate mean cosine of zenith angle for the full radiation step
      ROT2=ROT1+TWOPI*NRAD*DTsrc/SECONDS_PER_DAY
      CALL COSZS (ROT1,ROT2,COSZ2,COSZA)
      JDAYR=dayOfYear
      JYEARR=YEAR

      S0=S0X*S00WM2*RATLS0/RSDIST

c**** find scaling factors for surface albedo reduction
      if(dalbsnX.ne.0.) then
      IF (HAVE_SOUTH_POLE) THEN
         sumda_psum(:,1)=axyp(1,1)
         tauda_psum(:,1)=axyp(1,1)*depobc_1990(1,1)
      End If
      do j=J_0S,J_1S
      do i=I_0,I_1
c ilon72, jlat46 are indices w.r.t 72x46 grid
c      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5)   ! slightly more general
c      ILON72=INT(.5+(I-.5)*72./IM+.5)
        ilon72 = 1 + int( 72d0*lon2d(i,j)/twopi )
        jlat46 = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
        fsnow = flice(i,j) + rsi(i,j)*(1-fland(i,j))
        if(SNOWE_COM(I,J).gt.0.) fsnow = fsnow+fearth(i,j)
        sumda_psum(i,j) = axyp(i,j)*fsnow
        tauda_psum(i,j) = axyp(i,j)*fsnow*depobc_1990(i,j)
      end do
      end do
      IF (HAVE_NORTH_POLE) THEN
         sumda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)
         tauda_psum(:,JM)=axyp(1,jm)*rsi(1,jm)*depobc_1990(1,jm)
      END IF
      CALL GLOBALSUM(grid, sumda_psum,sumda,all=.true.)
      CALL GLOBALSUM(grid, tauda_psum,tauda,all=.true.)

      xdalbs=-dalbsnX*sumda/tauda
      IF(QCHECK) write(6,*) 'coeff. for snow alb reduction',xdalbs
      endif ! dalbsnX not zero

      if(kradia.le.0) then
      IF (QCHECK) THEN
C****   Calculate mean strat water conc
        STRATQ=0.
        MSTRAT=0.
        DO J=J_0,J_1
          STRJ=0.
          MSTJ=0.
          DO I=I_0,IMAXJ(J)
            DO L=LTROPO(I,J)+1,LM
              STRJ = STRJ + Q(I,J,L)*MA(L,I,J)*AXYP(I,J)
              MSTJ = MSTJ + MA(L,I,J)*AXYP(I,J)
            END DO
          END DO
          IF (J.eq.1 .or. J.eq.JM) THEN
            STRJ=STRJ*IM
            MSTJ=MSTJ*IM
          END IF
          STRATQ=STRATQ+STRJ
          MSTRAT=MSTRAT+MSTJ
        END DO
        PRINT*,"Strat water vapour (ppmv), mass (mb)",1d6*STRATQ*mair
     *       /(18.*MSTRAT),PMTOP+1d-2*GRAV*MSTRAT/AREAG
      END IF

C**** Get the random numbers outside openMP parallel regions
C**** but keep MC calculation separate from SS clouds
C**** To get parallel consistency also with mpi, force each process
C**** to generate random numbers for all latitudes (using BURN_RANDOM)

C**** MC clouds are considered as a block for each I,J grid point

      CALL BURN_RANDOM(nij_before_j0(J_0))

      DO J=J_0,J_1                    ! complete overlap
      CALL BURN_RANDOM((I_0-1))
      DO I=I_0,IMAXJ(J)
        RDMC(I,J) = RANDU(X)
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1))
      END DO

      CALL BURN_RANDOM((nij_after_j1(J_1)))

C**** SS clouds are considered as a block for each continuous cloud
      CALL BURN_RANDOM(nij_before_j0(j_0)*LM)

      DO J=J_0,J_1                    ! semi-random overlap
      CALL BURN_RANDOM((I_0-1)*LM)
      DO I=I_0,IMAXJ(J)
        NO_CLOUD_ABOVE = .TRUE.
        DO L=LM,1,-1
          IF(TAUSS(L,I,J).le.taulim) CLDSS(L,I,J)=0.
          IF(TAUMC(L,I,J).le.taulim) CLDMC(L,I,J)=0.
          RANDXX = RANDU(X)
          IF(CLDSS(L,I,J).GT.0.) THEN
            IF (NO_CLOUD_ABOVE) THEN
              RANDSS = RANDXX
              NO_CLOUD_ABOVE = .FALSE.
            END IF
          ELSE
            RANDSS = 1.
            NO_CLOUD_ABOVE = .TRUE.
          END IF
          RDSS(L,I,J) = RANDSS
        END DO
      END DO
      CALL BURN_RANDOM(nij_after_i1(I_1)*LM)
      END DO

      CALL BURN_RANDOM(nij_after_j1(j_1)*LM)

      end if                    ! kradia le 0

#ifdef TRACERS_ON
      ttausv_count=ttausv_count+1.d0
#endif
#ifdef ACCMIP_LIKE_DIAGS
! because of additional updghg calls, these factors won't apply:
      if(CO2X.ne.1.)  call stop_model('CO2x.ne.1 accmip diags',255)
      if(N2OX.ne.1.)  call stop_model('N2Ox.ne.1 accmip diags',255)
      if(CH4X.ne.1.)  call stop_model('CH4x.ne.1 accmip diags',255)
      if(CFC11X.ne.1.)call stop_model('CFC11x.ne.1 accmip diags',255)
      if(CFC12X.ne.1.)call stop_model('CFC12x.ne.1 accmip diags',255)
      if(XGHGX.ne.1.) call stop_model('XGHGx.ne.1 accmip diags',255)
      GFrefY=1850; GFrefD=182     ! ghg forcing refrnce year, day
      GFnowY=JyearR; GFnowD=JdayR ! ghg current desired year, day
      if(KJDAYG > 0) GFnowD=KJDAYG ! unless presribed in deck
      if(KYEARG > 0) GFnowY=KYEARG !           "
      call updghg(GFrefY,GFrefD)
      sv_fulgas_ref(1:4)=fulgas(nfghg(1:4))
      call updghg(GFnowY,GFnowD)
      sv_fulgas_now(1:4)=fulgas(nfghg(1:4))
#endif
#ifdef HEALY_LM_DIAGS
      FCO2=FULGAS(2)*CO2R
      FN2O=FULGAS(6)*N2OR
      FCH4=FULGAS(7)*CH4R
c
c      write(6,*) 'RJH: GHG: CONC=',
c     * FCO2,FN2O,FCH4
      ghg_totforc=5.35d0*log(FCO2/CO2I)
     *  +0.036d0*(sqrt(FCH4)-sqrt(CH4I))
     *  -(Fe(FCH4,N2OI)-Fe(CH4I,N2OI))
     *  + 0.12d0*(sqrt(FN2O)-sqrt(N2OI))
     *  -(Fe(CH4I,FN2O)-Fe(CH4I,N2OI))
c      write(6,*) 'RJH: GHG: FORC=',ghg_totforc
#endif

      aj_alb_inds = (/ J_PLAVIS, J_PLANIR, J_ALBVIS, J_ALBNIR,
     &                 J_SRRVIS, J_SRRNIR, J_SRAVIS, J_SRANIR /)

C****
C**** MAIN J LOOP
C****
      DO J=J_0,J_1

c     ICKERR=0
c     JCKERR=0
c     KCKERR=0

C****
C**** MAIN I LOOP
C****
      DO I=I_0,IMAXJ(J)
C**** Radiation input files use a 72x46 grid independent of IM and JM
C**** (ilon72,jlat46) is the 4x5 box containing the center of box (i,j)
c      JLAT46=INT(1.+(J-1.)*45./(JM-1.)+.5)  !  lat_index w.r.to 72x46 grid
c      JLAT46=INT(1.+(J-1.)*0.25*DLAT_DG+.5) ! slightly more general
c      ILON72=INT(.5+(I-.5)*72./IM+.5)  ! lon_index w.r.to 72x46 grid
      igcm = i
      jgcm = j
      ilon72 = 1 + int( 72d0*lon2d(i,j)/twopi )
      jlat46 = 1 + int( 45d0*(lat2d(i,j)+92d0*radian)/pi )
#ifdef ALTER_RADF_BY_LAT
      FULGAS(:)=FULGAS_orig(:)*FULGAS_lat(:,JLAT46)
      FS8OPX(:)=FS8OPX_orig(:)*FS8OPX_lat(:,JLAT46)
      FT8OPX(:)=FT8OPX_orig(:)*FT8OPX_lat(:,JLAT46)
#endif
      L1 = 1                         ! lowest layer above ground
      LMR=LM+LM_REQ                  ! radiation allows var. # of layers
      JR=JREG(I,J)
C**** DETERMINE FRACTIONS FOR SURFACE TYPES AND COLUMN PRESSURE
      PLAND=FLAND(I,J)
      POICE=RSI(I,J)*(1.-PLAND)
      POCEAN=(1.-PLAND)-POICE
      PLAKE=FLAKE(I,J)
      PLICE=FLICE(I,J)
      PEARTH=FEARTH(I,J)
      ptype4(1) = pocean ! open ocean and open lake
      ptype4(2) = poice  ! ocean/lake ice
      ptype4(3) = plice  ! glacial ice
      ptype4(4) = pearth ! non glacial ice covered soil

C**** CHECK SURFACE TEMPERATURES
      DO IT=1,4
        IF(ptype4(IT) > 0.) then
          IF(int(asflx4(it)%GTEMPR(I,J)) .LT. planck_tmin .OR.
     &       int(asflx4(it)%GTEMPR(I,J)) .GE. planck_tmax ) then
            WRITE(6,*) 'In Radia: Time,I,J,IT,TG1',ITime,I,J,IT
     *         ,asflx4(it)%GTEMPR(I,J)
CCC         STOP 'In Radia: Grnd Temp out of range'
c           ICKERR=ICKERR+1
          END IF
        END IF
      END DO

#if (defined CHL_from_SeaWIFs) || (defined TRACERS_OceanBiology)
C**** Set Chlorophyll concentration
      if (POCEAN.gt.0) then
          LOC_CHL = atmocn%chl(I,J)
          AIJ(I,J,IJ_CHL)=AIJ(I,J,IJ_CHL)+atmocn%CHL(I,J)*FOCEAN(I,J)
!         write(*,'(a,3i5,e12.4)')'RAD_DRV:',
!    .    itime,i,j,chl(i,j)
      endif
#endif

      LS1_loc=LTROPO(I,J)+1  ! define stratosphere for radiation
C**** kradia>1: adjusted forcing, i.e. T adjusts in L=LS1_loc->LM+3
      if(kradia>1) LS1_loc=LS1_loc+2-kradia     ! favorite:kradia=3
      if(kradia>3) LS1_loc=1                    ! favorite:kradia=3
      kdeliq=0   ! initialize mainly for L>LM
      if (kradia.gt.0) then     ! rad forcing model
        do l=1,lm
          tlm(l) = T(i,j,l)*pk(l,i,j)
          shl(l) = QR(l,i,j) ; if(shl(l)<0) shl(l)=0
          tauwc(l) = cldx*CLDinfo(l,1,i,j)
          tauic(l) = cldx*CLDinfo(l,2,i,j)
          SIZEWC(L)= CLDinfo(l,3,i,j)
          SIZEIC(L)= SIZEWC(L)
        end do
      else                      ! full model
C****
C**** DETERMINE CLOUDS (AND THEIR OPTICAL DEPTHS) SEEN BY RADIATION
C****
      CSS=0. ; CMC=0. ; CLDCV=0. ; DEPTH=0. ; OPTDW=0. ; OPTDI=0.
      if(cc_cdncx.ne.0. .or. od_cdncx.ne.0.) then
        call dCDNC_EST(i,j,pland, dCDNC)
      else
        dCDNC = 0.
      endif
      dCC_CDNCL = CC_cdncx*dCDNC*CDNCL
      dOD_CDNCL = OD_cdncx*dCDNC*CDNCL
      DO L=1,LM
        if(q(i,j,l)<0) then
           WRITE(6,*)'In Radia: Time,I,J,L,Q<0',ITime,I,J,L,Q,'->0'
           Q(I,J,L)=0.
        end if
        QSS=Q(I,J,L)/(RHSAV(L,I,J)+1.D-20)
        shl(L)=QSS
        IF(FSS(L,I,J)*CLDSAV(L,I,J).LT.1.)
     *       shl(L)=(Q(I,J,L)-QSS*FSS(L,I,J)*CLDSAV(L,I,J))/
     /              (1.-FSS(L,I,J)*CLDSAV(L,I,J))
        TLm(L)=T(I,J,L)*PK(L,I,J)
        TAUSSL=0.
        TAUMCL=0.
        TAUWC(L)=0.
        TAUIC(L)=0.
        SIZEWC(L)=0.
        SIZEIC(L)=0.
        TOTCLD(L)=0.
C**** Determine large scale and moist convective cloud cover for radia
        IF (CLDSS(L,I,J)*(1.+dcc_cdncl(l)).GT.RDSS(L,I,J)) THEN
          TAUSSL=TAUSS(L,I,J)*(1.+dod_cdncl(l))
          shl(L)=QSS
          CSS=1.
          call inc_ajl(i,j,l,jl_sscld,css)
        END IF
        IF (CLDMC(L,I,J).GT.RDMC(I,J)) THEN
          CMC=1.
          call inc_ajl(i,j,l,jl_mccld,cmc)
          DEPTH=DEPTH+PDSIG(L,I,J)
          IF(TAUMC(L,I,J).GT.TAUSSL) THEN
            TAUMCL=TAUMC(L,I,J)
            ELHX=LHE
            IF(TLm(L).LE.TF) ELHX=LHS
            shl(L)=QSAT(TLm(L),ELHX,PMID(L,I,J))
          END IF
        END IF
        IF(TAUSSL+TAUMCL.GT.0.) THEN
             CLDCV=1.
          TOTCLD(L)=1.
          call inc_ajl(i,j,l,jl_totcld,1d0)
C**** save 3D cloud fraction as seen by radiation
          if(cldx>0) AIJL(I,J,L,IJL_CF)=AIJL(I,J,L,IJL_CF)+1.
          IF(TAUMCL.GT.TAUSSL) THEN
            SIZEWC(L)=CSIZMC(L,I,J)
            SIZEIC(L)=CSIZMC(L,I,J)
            IF(SVLAT(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUMCL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
            ELSE
              TAUIC(L)=cldx*TAUMCL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
            END IF
          ELSE
            SIZEWC(L)=CSIZSS(L,I,J)
            SIZEIC(L)=CSIZSS(L,I,J)
            IF(SVLHX(L,I,J).EQ.LHE) THEN
              TAUWC(L)=cldx*TAUSSL
              OPTDW=OPTDW+TAUWC(L)
              call inc_ajl(i,j,l,jl_wcld,1d0)
              call inc_ajl(i,j,l,jl_wcldwt,pdsig(l,i,j))
            ELSE
              TAUIC(L)=cldx*TAUSSL
              OPTDI=OPTDI+TAUIC(L)
              call inc_ajl(i,j,l,jl_icld,1d0)
              call inc_ajl(i,j,l,jl_icldwt,pdsig(l,i,j))
            END IF
          END IF
          call inc_ajl(i,j,l,jl_wcod,tauwc(l))
          call inc_ajl(i,j,l,jl_icod,tauic(l))
          call inc_ajl(i,j,l,jl_wcsiz,sizewc(l)*tauwc(l))
          call inc_ajl(i,j,l,jl_icsiz,sizeic(l)*tauic(l))
        END IF
C**** save some radiation/cloud fields for wider use
        RCLD(L,I,J)=TAUWC(L)+TAUIC(L)
      END DO
      CFRAC(I,J) = CLDCV    ! cloud fraction consistent with radiation
C**** effective cloud cover diagnostics
         OPNSKY=1.-CLDCV
         DO IT=1,NTYPE
           call inc_aj(i,j,it,J_PCLDSS,CSS  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLDMC,CMC  *FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_CLDDEP,DEPTH*FTYPE(IT,I,J))
           call inc_aj(i,j,it,J_PCLD  ,CLDCV*FTYPE(IT,I,J))
         END DO
         call inc_areg(i,j,jr,J_PCLDSS,CSS  )
         call inc_areg(i,j,jr,J_PCLDMC,CMC  )
         call inc_areg(i,j,jr,J_CLDDEP,DEPTH)
         call inc_areg(i,j,jr,J_PCLD  ,CLDCV)
         AIJ(I,J,IJ_PMCCLD)=AIJ(I,J,IJ_PMCCLD)+CMC
         AIJ(I,J,IJ_CLDCV) =AIJ(I,J,IJ_CLDCV) +CLDCV
         DO L=1,LLOW
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDL)=AIJ(I,J,IJ_PCLDL)+1.
           exit
         end do
         DO L=LLOW+1,LMID
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDM)=AIJ(I,J,IJ_PCLDM)+1.
           exit
         end do
         DO L=LMID+1,LHI
           IF (TOTCLD(L).NE.1.) cycle
           AIJ(I,J,IJ_PCLDH)=AIJ(I,J,IJ_PCLDH)+1.
           exit
         end do

         TAUSUMW(I,J) = OPTDW
         TAUSUMI(I,J) = OPTDI
         if(optdw.gt.0.) then
            AIJ(I,J,IJ_optdw)=AIJ(I,J,IJ_optdw)+optdw
            AIJ(I,J,IJ_wtrcld)=AIJ(I,J,IJ_wtrcld)+1.
         end if
         if(optdi.gt.0.) then
            AIJ(I,J,IJ_optdi)=AIJ(I,J,IJ_optdi)+optdi
            AIJ(I,J,IJ_icecld)=AIJ(I,J,IJ_icecld)+1.
         end if

         DO KR=1,NDIUPT
           IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
C**** Warning: this replication may give inaccurate results for hours
C****          1->(NRAD-1)*DTsrc (ADIURN) or skip them (HDIURN)
             TMP(IDD_CL7:IDD_CL7+6)=TOTCLD(1:7)
             TMP(IDD_CCV)=CLDCV
             DO INCH=1,NRAD
               IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
               IH=IHM
               IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
               ADIURN(IDX(:),KR,IH)=ADIURN(IDX(:),KR,IH)+TMP(IDX(:))
#ifndef NO_HDIURN
               IHM = IHM+(DATE-1)*HR_IN_DAY
               IF(IHM.LE.HR_IN_MONTH) THEN
                 HDIURN(IDX(:),KR,IHM)=HDIURN(IDX(:),KR,IHM)+TMP(IDX(:))
               ENDIF
#endif
             ENDDO
           END IF
         END DO
      end if ! kradia le 0 (full model)
C****
C**** SET UP VERTICAL ARRAYS OMITTING THE I AND J INDICES
C****
C**** EVEN PRESSURES
#ifdef TRACERS_TOMAS
      TTAUSV(:,:)=0.
      aesqex(:,:,:)=0.0
      aesqsc(:,:,:)=0.0
      aesqcb(:,:,:)=0.0
#endif
      PLB(LM+1)=PEDN(LM+1,I,J)
      DO L=1,LM
        PLB(L)=PEDN(L,I,J)
C**** TEMPERATURES
C---- TLm(L)=T(I,J,L)*PK(L,I,J)     ! already defined
        IF(int(TLm(L)) .LT. planck_tmin .OR.
     &     int(TLm(L)) .GE. planck_tmax ) THEN
          WRITE(6,*) 'In Radia: Time,I,J,L,TL',ITime,I,J,L,TLm(L)
          WRITE(6,*) 'GTEMPR:',
     &         asflx4(1)%GTEMPR(I,J),asflx4(2)%GTEMPR(I,J),
     &         asflx4(3)%GTEMPR(I,J),asflx4(4)%GTEMPR(I,J)
CCC       STOP 'In Radia: Temperature out of range'
c         ICKERR=ICKERR+1
        END IF
C**** MOISTURE VARIABLES
C---- shl(L)=Q(I,J,L)        ! already defined and reset to 0 if <0
c       if(shl(l).lt.0.) then
c         WRITE(0,*)'In Radia: Time,I,J,L,QL<0',ITime,I,J,L,shl(L),'->0'
c         KCKERR=KCKERR+1
c         shl(l)=0.
c       end if
        RHL(L) = shl(L)/QSAT(TLm(L),LHE,PMID(L,I,J))
        if(RHfix.ge.0.) RHL(L)=RHfix
C**** Extra aerosol data
C**** For up to NTRACE aerosols, define the aerosol amount to
C**** be used (kg/m^2)
C**** Only define TRACER is individual tracer is actually defined.
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
C**** loop over tracers that are passed to radiation.
C**** Some special cases for black carbon, organic carbon, SOAs where
C**** more than one tracer is lumped together for radiation purposes
      do n=1,NTRACE
        if (NTRIX(n).gt.0) then
          select case (trname(NTRIX(n)))
          case ("OCIA", 'vbsAm2')
            TRACER(L,n)=(
#ifdef TRACERS_AEROSOLS_VBS
     *           sum(trm(i,j,l,vbs_tr%iaer))
#else
     *           trm(i,j,l,n_OCII)+trm(i,j,l,n_OCIA)
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_OCEAN
     *          +trm(i,j,l,n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
     *                  )*BYAXYP(I,J)
          case ("OCB")
#ifdef TRACERS_AEROSOLS_VBS
           TRACER(L,n)=0.d0
#else
           TRACER(L,n)=trm(i,j,l,n_OCB)*BYAXYP(I,J)
           ! The only reason this is a special case is, following
           ! how OCIA is done, it has no wttr factor like default does.
#endif  /* TRACERS_AEROSOLS_VBS */
#ifdef TRACERS_AEROSOLS_SOA
          case ("isopp1a")
           TRACER(L,n)=( trm(i,j,l,n_isopp1a)+trm(i,j,l,n_isopp2a)
#ifdef TRACERS_TERP
     &                  +trm(i,j,l,n_apinp1a)+trm(i,j,l,n_apinp2a)
#endif /* TRACERS_TERP */
     &                 )*BYAXYP(I,J)
#endif /* TRACERS_AEROSOLS_SOA */
          case ("BCIA")
           TRACER(L,n)=(trm(i,j,l,n_BCII)+trm(i,j,l,n_BCIA))*BYAXYP(I,J)
          case default
#ifdef TRACERS_NITRATE
! assume full neutralization of NO3p, if NH4 suffice
           select case (trname(NTRIX(n)))
           case ("NO3p")
            if (trm(i,j,l,NTRIX(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(nh4_on_no3+trm(i,j,l,NTRIX(n)))/
     *                trm(i,j,l,NTRIX(n))
            endif
           case ("SO4")
            if (trm(i,j,l,NTRIX(n)) > 0.d0) then
              nh4_on_no3=min(trm(i,j,l,n_NO3p)*(tr_mm(n_NO3p)+
     *            tr_mm(n_NH4))/tr_mm(n_NO3p)-trm(i,j,l,n_NO3p),
     *                       trm(i,j,l,n_NH4))
              wttr(n)=(trm(i,j,l,n_NH4)-nh4_on_no3+trm(i,j,l,NTRIX(n)))/
     *                trm(i,j,l,NTRIX(n))
            endif
           end select
#endif
           TRACER(L,n)=wttr(n)*trm(i,j,l,NTRIX(n))*BYAXYP(I,J)
          end select
        end if
      end do
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/QUARZHEM/OM_SP */

#ifdef TRACERS_AMP
      CALL SETAMP_LEV(i,j,l)
#endif
#ifdef TRACERS_TOMAS
      CALL SETTOMAS_LEV(i,j,l)
#endif
      END DO
C**** Radiative Equilibrium Layer data
      DO K=1,LM_REQ
        IF(int(RQT(K,I,J)) .LT. planck_tmin .OR.
     &     int(RQT(K,I,J)) .GE. planck_tmax ) THEN
        WRITE(6,*) 'In RADIA: Time,I,J,L,TL',ITime,I,J,LM+K,RQT(K,I,J)
CCC     STOP 'In Radia: RQT out of range'
c       JCKERR=JCKERR+1
        END IF
        TLm(LM+K)=RQT(K,I,J)
        PLB(LM+k+1) = PLB0(k)
        shl(LM+k)    = shl0(k)
        RHL(LM+k) = shl(LM+k)/QSAT(TLm(LM+k),LHE,
     *                                .5d0*(PLB(LM+k)+PLB(LM+k+1)) )
        tauwc(LM+k) = 0.
        tauic(LM+k) = 0.
        sizewc(LM+k)= 0.
        sizeic(LM+k)= 0.
#ifdef TRACERS_ON
C**** set radiative equilibirum extra tracer amount to zero
        IF (NTRACE.gt.0) TRACER(LM+k,1:NTRACE)=0.
#endif
      END DO
      if (kradia.gt.1) then
        do l=1,lm+lm_req
          tlm(l) = tlm(l) + Tchg(l,i,j)
          AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+Tchg(L,I,J)
        end do
      end if
C**** Zenith angle and GROUND/SURFACE parameters
      COSZ=COSZA(I,J)
      TGO =atmocn%GTEMPR(I,J)
      TGOI=atmice%GTEMPR(I,J)
      TGLI=atmgla%GTEMPR(I,J)
      TGE =atmlnd%GTEMPR(I,J)
      TSL=atmsrf%TSAVG(I,J)
      SNOWOI=SNOWI(I,J)
      SNOWLI=atmgla%SNOW(I,J)
      SNOWE=SNOWE_COM(I,J)                    ! snow depth (kg/m**2)
      snow_frac(:) = fr_snow_rad_ij(:,i,j)    ! snow cover (1)
      AGESN(1)=SNOAGE(3,I,J)    ! land         ! ? why are these numbers
      AGESN(2)=SNOAGE(1,I,J)    ! ocean ice        so confusing ?
      AGESN(3)=SNOAGE(2,I,J)    ! land ice
c      print*,"snowage",i,j,SNOAGE(1,I,J)
C**** set up parameters for new sea ice and snow albedo
      zsnwoi=snowoi/rhos
      if(dalbsnX.ne.0.) then
        dALBsn = xdalbs*depobc(i,j)
      else
        dALBsn = 0.
      endif
c to use on-line tracer albedo impact, set dALBsnX=0. in rundeck
#if (defined BC_ALB) &&\
    ((defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS))
      call GET_BC_DALBEDO(i,j,dALBsn1)
      dALBsn=dALBsn1
#endif
      if (poice.gt.0.) then
        zoice=(ace1i+msi(i,j))/rhoi
        flags=flag_dsws(i,j)
        if (kradia .le. 0) then
          fmp=min(1.6d0*sqrt(pond_melt(i,j)/rhow),1d0)
             AIJ(I,J,IJ_FRMP) = AIJ(I,J,IJ_FRMP) + fmp*POICE
        else
          fmp = fmp_com(i,j)
        end if
        zmp=min(0.8d0*fmp,0.9d0*zoice)
      else
        zoice=0. ; flags=.FALSE. ; fmp=0. ; zmp=0.
      endif
C**** set up new lake depth parameter to incr. albedo for shallow lakes
      zlake=0.
      if (plake.gt.0) then
        zlake = MWL(I,J)/(RHOW*PLAKE*AXYP(I,J))
      end if
C****
      if (kradia .le. 0) then
        !WEARTH=(WEARTH_COM(I,J)+AIEARTH(I,J))/(WFCS(I,J)+1.D-20)
        WEARTH=atmlnd%bare_soil_wetness(i,j)
        if (wearth.gt.1.) wearth=1.
      else                            ! rad.frc. model
        wearth = wsoil(i,j)
      end if
#ifdef USE_ENT
      if ( fearth(i,j) > 0.d0 ) then
        call ent_get_exports( entcells(i,j),
     &       vegetation_fractions=PVT0 )
        call map_ent2giss(PVT0,PVT) !temp hack: ent pfts->giss veg
      else
        PVT(:) = 0.d0  ! actually PVT is not supposed to be used in this case
      endif
#else
      DO K=1,12
        PVT(K)=VDATA(I,J,K)
      END DO
#endif
      WMAG=atmsrf%WSAVG(I,J)
C****
C**** Radiative interaction and forcing diagnostics:
C**** If no radiatively active tracers are defined, nothing changes.
C**** Currently this works for aerosols and ozone but should be extended
C**** to cope with all trace gases.
C****
      FSTOPX(:)=1. ; FTTOPX(:)=1. ; FTAUC=1. ! deflt (aeros/clouds on)
      use_tracer_chem(:) = 0 ! by default use climatological ozone/ch4
C**** Set level for inst. rad. forc. calcs for aerosols/trace gases
C**** This is set from the rundeck.
      LFRC=LM+LM_REQ+1          ! TOA
      if (rad_forc_lev.gt.0) LFRC=LTROPO(I,J) ! TROPOPAUSE
#ifdef ACCMIP_LIKE_DIAGS
      if(rad_forc_lev.gt.0)call stop_model
     &('ACCMIP_LIKE_DIAGS desires TOA RADF diags',255)
#endif
C**** The calculation of the forcing is slightly different.
C**** depending on whether full radiative interaction is turned on
C**** or not.
      onoff_aer=0; onoff_chem=0
      if (rad_interact_aer > 0) onoff_aer=1
      if (clim_interact_chem > 0) onoff_chem=1
      use_o3_ref=0

C YUNHA LEE - took the shindell outside of the Koch/dust directives.
#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone and Methane:
      CHEM_IN(1,1:LM)=chem_tracer_save(1,1:LM,I,J)
      CHEM_IN(2,1:LM)=chem_tracer_save(2,1:LM,I,J)*CH4X_RADoverCHEM
      if (clim_interact_chem > 0) then
        use_tracer_chem(1)=Lmax_rad_O3  ! O3
        use_tracer_chem(2)=Lmax_rad_CH4 ! CH4
      endif
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
      if(clim_interact_chem<=0)
     &call stop_model("stratOx RADF on, clim_interact_chem<=0",255)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
C**** Aerosols incl. Dust:        set up for radiative forcing diagnostics
      if (NTRACE>0 .and. moddrf==0) then
        FSTOPX(:)=onoff_aer ; FTTOPX(:)=onoff_aer
        set_clayilli=.FALSE.
        set_claykaol=.FALSE.
        set_claysmec=.FALSE.
        set_claycalc=.FALSE.
        set_clayquar=.FALSE.
        do n=1,NTRACE
          IF (ntrix(n) > 0) THEN
            IF (trname(NTRIX(n)).eq."seasalt2") CYCLE ! not for seasalt2
            IF (trname(ntrix(n)) == 'ClayIlli' .AND. set_clayilli) cycle
            IF (trname(ntrix(n)) == 'ClayKaol' .AND. set_claykaol) cycle
            IF (trname(ntrix(n)) == 'ClaySmec' .AND. set_claysmec) cycle
            IF (trname(ntrix(n)) == 'ClayCalc' .AND. set_claycalc) cycle
            IF (trname(ntrix(n)) == 'ClayQuar' .AND. set_clayquar) cycle
            FSTOPX(n)=1-onoff_aer ; FTTOPX(n)=1-onoff_aer ! turn on/off tracer
C**** Warning: small bit of hardcoding assumes that seasalt2 immediately
C****          succeeds seasalt1 in NTRACE array
            IF (trname(NTRIX(n)).eq."seasalt1") THEN          !add seasalt2
              FSTOPX(n+1)=1-onoff_aer;FTTOPX(n+1)=1-onoff_aer !to seasalt1
            END IF
C**** Do radiation calculations for all clay classes at once
C**** Assumes that 4 clay tracers are adjacent in NTRACE array
            SELECT CASE (trname(ntrix(n)))
            CASE ('ClayIlli')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_clayilli=.true.
            CASE ('ClayKaol')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claykaol=.true.
            CASE ('ClaySmec')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claysmec=.true.
            CASE ('ClayCalc')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_claycalc=.true.
            CASE ('ClayQuar')
              fstopx(n+1:n+3)=1-onoff_aer; fttopx(n+1:n+3)=1-onoff_aer
              set_clayquar=.true.
            END SELECT
            kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
            CALL RCOMPX  ! tr.aero.Koch/dust/miner./quarz/om_sp
            SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
            TNFST(1,n,I,J)=TRNFLB(1)
            SNFST(2,n,I,J)=SRNFLB(LFRC)
            TNFST(2,n,I,J)=TRNFLB(LFRC)
            FSTOPX(n)=onoff_aer ; FTTOPX(n)=onoff_aer   ! back to default
            IF (trname(NTRIX(n)).eq."seasalt1") THEN    ! also for seasalt2
              FSTOPX(n+1)=onoff_aer ; FTTOPX(n+1)=onoff_aer
            END IF
            SELECT CASE (trname(ntrix(n)))           ! also for clays
            CASE ('ClayIlli','ClayKaol','ClaySmec','ClayCalc',
     &           'ClayQuar')
              fstopx(n+1:n+3)=onoff_aer ; fttopx(n+1:n+3)=onoff_aer
            END SELECT
          END IF
        end do
      end if
#endif /* TRACERS_AEROSOLS_Koch/DUST/MINERALS/QUARZHEM/OM_SP */

      if (moddrf==0) then
#ifdef TRACERS_SPECIAL_Shindell
C**** Ozone:
! ozone rad forcing diags now use a constant reference year
! for this first call. And no tracer values...
        use_o3_ref=1 ; use_tracer_chem(1)=0
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
        SNFST_o3ref(1,I,J)=SRNFLB(LTROPO(I,J)) ! tropopause
        TNFST_o3ref(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_o3ref(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(2,I,J)=TRNFLB(LM+LM_REQ+1)
#ifdef AUXILIARY_OX_RADF
! if needed, also save the auxiliary ozone field (i.e. climatology
! if tracer is used in final call, tracers if climatology is used.)
        use_o3_ref=0 ; use_tracer_chem(1)=(1-onoff_chem)*Lmax_rad_O3
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX        ! tr_Shindell Ox tracer
        SNFST_o3ref(3,I,J)=SRNFLB(LTROPO(I,J)) ! tropopause
        TNFST_o3ref(3,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_o3ref(4,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_o3ref(4,I,J)=TRNFLB(LM+LM_REQ+1)
#endif
! ... after which it can use either climatological or tracer O3:
        use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
! Optional intermediate call with stratOx tracer:
!NEED chem_IN(1,1:LM)=stratO3_tracer_save(1:LM,I,J)
!NEED kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
!NEED CALL RCOMPX        ! stratOx diag tracer
        SNFST_stratOx(1,I,J)=SRNFLB(LTROPO(I,J)) ! tropopause
        TNFST_stratOx(1,I,J)=TRNFLB(LTROPO(I,J))
        SNFST_stratOx(2,I,J)=SRNFLB(LM+LM_REQ+1) ! T.O.A.
        TNFST_stratOx(2,I,J)=TRNFLB(LM+LM_REQ+1)
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS */
        chem_IN(1,1:LM)=chem_tracer_save(1,1:LM,I,J)  ! Ozone
        chem_IN(2,1:LM)=chem_tracer_save(2,1:LM,I,J)*CH4X_RADoverCHEM  ! Methane
#ifdef ACCMIP_LIKE_DIAGS
! TOA GHG rad forcing: nf=1,4 are CH4, N2O, CFC11, and CFC12:
! Initial calls are reference year/day:
      do nf=1,4
        if(nf==1)then ! CH4 reference call must not use tracer
          use_tracer_chem(2)=0
        else ! N2O and CFC call's CH4 should match final call
          use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
        endif
        fulgas(nfghg(nf))=sv_fulgas_ref(nf)
        kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
        CALL RCOMPX
        SNFS_ghg(nf,I,J)=SRNFLB(LM+LM_REQ+1)
        TNFS_ghg(nf,I,J)=TRNFLB(LM+LM_REQ+1)
        fulgas(nfghg(nf))=sv_fulgas_now(nf)
      enddo
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_SPECIAL_Shindell */
      end if ! moddrf=0
#ifdef TRACERS_SPECIAL_Shindell
! final (main) RCOMPX call can use tracer methane (or not):
      use_tracer_chem(2)=onoff_chem*Lmax_rad_CH4
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
          use_tracer_chem(2)=0  ! special case; model outputs climatology
        end if
      end if
#endif /* TRACERS_SPECIAL_Shindell */


      if (moddrf==0) then
#ifdef BC_ALB
        dalbsn=0.d0
        CALL RCOMPX
        NFSNBC(I,J)=SRNFLB(LM+LM_REQ+1)
c       NFSNBC(I,J)=SRNFLB(LFRC)
        ALBNBC(I,J)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
c set for BC-albedo effect
        dALBsn=dALBsn1
#endif
#ifdef TRACERS_AMP
        IF (AMP_DIAG_FC == 2) THEN
          Do n = 1,nmodes
            FSTOPX(n) = 1-onoff_aer !turns off online tracer
            FTTOPX(n) = 1-onoff_aer !
            if (n.eq.1) FSTOPX(:) = 1-onoff_aer
            if (n.eq.1) FTTOPX(:) = 1-onoff_aer
            CALL RCOMPX
            SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
            TNFST(1,n,I,J)=TRNFLB(1)
            SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
            TNFST(2,n,I,J)=TRNFLB(LFRC)
            FSTOPX(:) = onoff_aer !turns on online tracer
            FTTOPX(:) = onoff_aer !
          ENDDO
        ELSE
           n = 1
          FSTOPX(:) = 1-onoff_aer !turns off online tracer
          FTTOPX(:) = 1-onoff_aer !
          CALL RCOMPX
          SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
          TNFST(1,n,I,J)=TRNFLB(1)
          SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
          TNFST(2,n,I,J)=TRNFLB(LFRC)
          FSTOPX(:) = onoff_aer !turns on online tracer
          FTTOPX(:) = onoff_aer !
        ENDIF
#endif
#ifdef TRACERS_TOMAS
        IF (TOMAS_DIAG_FC == 2) THEN
            FSTOPX(:) = onoff_aer !turns on online tracer
            FTTOPX(:) = onoff_aer !
          Do n = 1,icomp-2
            FSTOPX(n) = 1-onoff_aer !turns off online tracer
            FTTOPX(n) = 1-onoff_aer !
            CALL RCOMPX
            SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
            TNFST(1,n,I,J)=TRNFLB(1)
            SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
            TNFST(2,n,I,J)=TRNFLB(LFRC)
            FSTOPX(n) = onoff_aer !turns on online tracer
            FTTOPX(n) = onoff_aer !
          ENDDO
        ELSE
           n = 1
          FSTOPX(:) = 1-onoff_aer !turns off online tracer
          FTTOPX(:) = 1-onoff_aer !
          CALL RCOMPX
          SNFST(1,n,I,J)=SRNFLB(1) ! surface forcing
          TNFST(1,n,I,J)=TRNFLB(1)
          SNFST(2,n,I,J)=SRNFLB(LFRC) ! Tropopause forcing
          TNFST(2,n,I,J)=TRNFLB(LFRC)
          FSTOPX(:) = onoff_aer !turns on online tracer
          FTTOPX(:) = onoff_aer !
        ENDIF
#endif
C**** Optional calculation of CRF using a clear sky calc.
        if (cloud_rad_forc.gt.0) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX          ! cloud_rad_forc>0 : clr sky
          SNFSCRF(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA
#ifdef SCM
          CSSRNTOP = SRNFLB(LM+LM_REQ+1)*COSZ2(I,J)
          CSTRUTOP = TRUFLB(LM+LM_REQ+1)
          CSSRNBOT = SRNFLB(1)*COSZ2(I,J)
          CSTRNBOT = TRNFLB(1)
          CSSRDBOT = SRDFLB(1)*COSZ2(I,J)
#endif
C         BEGIN AMIP
          AIJ(I,J,IJ_SWDCLS)=AIJ(I,J,IJ_SWDCLS)+SRDFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_SWNCLS)=AIJ(I,J,IJ_SWNCLS)+SRNFLB(1)*COSZ2(I,J)
          AIJ(I,J,IJ_LWDCLS)=AIJ(I,J,IJ_LWDCLS)+TRDFLB(1)
          AIJ(I,J,IJ_SWNCLT)=AIJ(I,J,IJ_SWNCLT)+SRNFLB(LM+LM_REQ+1)
     *     *COSZ2(I,J)
          AIJ(I,J,IJ_LWNCLT)=AIJ(I,J,IJ_LWNCLT)+TRNFLB(LM+LM_REQ+1)
C       END AMIP
        end if
        FTAUC=1.     ! default: turn on cloud tau


C**** 2nd Optional calculation of CRF using a clear sky calc. without aerosols and Ox
        if (cloud_rad_forc.gt.0) then
          FTAUC=0.   ! turn off cloud tau (tauic +tauwc)
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
c Including turn off of aerosols and Ox during crf calc.+++++++++++++++++++
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=1 ; use_tracer_chem(1)=0  !turns off ozone
#endif
       FSTOPX(:) = 1-onoff_aer !turns off aerosol tracer
       FTTOPX(:) = 1-onoff_aer !
        CALL RCOMPX          ! cloud_rad_forc>0 : clr sky
       FSTOPX(:) = onoff_aer !turns on aerosol tracer
       FTTOPX(:) = onoff_aer !
#ifdef TRACERS_SPECIAL_Shindell
       use_o3_ref=0 ; use_tracer_chem(1)=onoff_chem*Lmax_rad_O3 ! turns on ozone tracers
#endif
          SNFSCRF2(I,J)=SRNFLB(LM+LM_REQ+1)   ! always TOA
          TNFSCRF2(I,J)=TRNFLB(LM+LM_REQ+1)   ! always TOA

c          AIJ(I,J,IJ_SWDCLS2)=AIJ(I,J,IJ_SWDCLS2)+SRDFLB(1)*COSZ2(I,J)
c          AIJ(I,J,IJ_SWNCLS2)=AIJ(I,J,IJ_SWNCLS2)+SRNFLB(1)*COSZ2(I,J)
c          AIJ(I,J,IJ_LWDCLS2)=AIJ(I,J,IJ_LWDCLS2)+TRDFLB(1)
c          AIJ(I,J,IJ_SWNCLT2)=AIJ(I,J,IJ_SWNCLT2)+SRNFLB(LM+LM_REQ+1)
c     *     *COSZ2(I,J)
c          AIJ(I,J,IJ_LWNCLT2)=AIJ(I,J,IJ_LWNCLT2)+TRNFLB(LM+LM_REQ+1)
        end if
        FTAUC=1.     ! default: turn on cloud tau


C**** Optional calculation of the impact of default aerosols
        if (aer_rad_forc.gt.0) then
C**** first, separate aerosols
          DO N=1,8
          tmpS(N)=FS8OPX(N)   ; tmpT(N)=FT8OPX(N)
          FS8OPX(N)=0.     ; FT8OPX(N)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX           ! aer_rad_forc>0 : no aerosol N
          SNFSAERRF(N,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(N,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(N+8,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(N+8,I,J)=TRNFLB(1) ! SURF
          FS8OPX(N)=tmpS(N)   ; FT8OPX(N)=tmpT(N)
          END DO
C**** second, net aerosols
          tmpS(:)=FS8OPX(:)   ; tmpT(:)=FT8OPX(:)
          FS8OPX(:)=0.     ; FT8OPX(:)=0.
          kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)
          CALL RCOMPX             ! aer_rad_forc>0 : no aerosols
          SNFSAERRF(17,I,J)=SRNFLB(LM+LM_REQ+1) ! TOA
          TNFSAERRF(17,I,J)=TRNFLB(LM+LM_REQ+1) ! TOA
          SNFSAERRF(18,I,J)=SRNFLB(1) ! SURF
          TNFSAERRF(18,I,J)=TRNFLB(1) ! SURF
          FS8OPX(:)=tmpS(:)   ; FT8OPX(:)=tmpT(:)
        end if
      end if  ! moddrf=0
C**** End of initial computations for optional forcing diagnostics

C**** Localize fields that are modified by RCOMPX
      kdeliq(1:lm,1:4)=kliq(1:lm,1:4,i,j)

C*****************************************************
C     Main RADIATIVE computations, SOLAR and THERM(A)L
      CALL RCOMPX
C*****************************************************

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)

#ifdef TRACERS_AMP
      NTRACE = nmodes
#endif
#ifdef TRACERS_TOMAS
      NTRACE = icomp-2
#endif
C**** Save optical depth diags
      do n=1,NTRACE
        IF (ntrix(n) > 0) THEN
          SELECT CASE (trname(ntrix(n)))
          CASE ('Clay')
            n1=n-nrad_clay+1
            IF (diag_rad /= 1) THEN
              IF (ijts_tausub(1,ntrix(n),n1) > 0)
     &             taijs(i,j,ijts_tausub(1,ntrix(n),n1))
     &             =taijs(i,j,ijts_tausub(1,ntrix(n),n1))
     &             +SUM(ttausv(1:Lm,n))
              IF (ijts_tausub(2,ntrix(n),n1) > 0)
     &             taijs(i,j,ijts_tausub(2,ntrix(n),n1))
     &             =taijs(i,j,ijts_tausub(2,ntrix(n),n1))
     &             +SUM(ttausv(1:Lm,n))*OPNSKY
            END IF
            if (ijlt_3Daaod(NTRIX(n)).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Daaod(NTRIX(n)))
     *           =taijls(i,j,1:lm,ijlt_3Daaod(NTRIX(n)))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3Dtau(NTRIX(n)).gt.0)
     *           taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))
     *           =taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))+TTAUSV(1:lm,n)
            IF (diag_rad == 1) THEN
              DO kr=1,6
                IF (ijts_sqexsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqexsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqexsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqexsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqexsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqexsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqex(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqscsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqscsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqscsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqsc(1:Lm,kr,n))
                IF (ijts_sqscsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqscsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqscsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqsc(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqcbsub(1,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqcbsub(1,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqcbsub(1,kr,ntrix(n),n1))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
                IF (ijts_sqcbsub(2,kr,ntrix(n),n1) > 0)
     &               taijs(i,j,ijts_sqcbsub(2,kr,ntrix(n),n1))
     &               =taijs(i,j,ijts_sqcbsub(2,kr,ntrix(n),n1))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
              END DO
            END IF
          CASE DEFAULT

            IF (diag_rad /= 1) THEN
              if (ijts_tau(1,NTRIX(n)).gt.0)
     &             taijs(i,j,ijts_tau(1,NTRIX(n)))
     &             =taijs(i,j,ijts_tau(1,NTRIX(n)))+SUM(TTAUSV(1:lm,n))
              if (ijts_tau(2,NTRIX(n)).gt.0)
     &             taijs(i,j,ijts_tau(2,NTRIX(n)))
     &             =taijs(i,j,ijts_tau(2,NTRIX(n)))
     &             +SUM(TTAUSV(1:lm,n))*OPNSKY
            END IF
            if (ijlt_3Daaod(NTRIX(n)).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Daaod(NTRIX(n)))
     &           =taijls(i,j,1:lm,ijlt_3Daaod(NTRIX(n)))+
     *            (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
            if (ijlt_3Dtau(NTRIX(n)).gt.0)
     &           taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))
     &           =taijls(i,j,1:lm,ijlt_3Dtau(NTRIX(n)))+TTAUSV(1:lm,n)
            IF (diag_rad == 1) THEN
              DO kr=1,6
c                 print*,'SUSA  diag',SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqex(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqex(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqex(1,kr,ntrix(n)))
     &               +SUM(aesqex(1:Lm,kr,n))
                IF (ijts_sqex(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqex(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqex(2,kr,ntrix(n)))
     &               +SUM(aesqex(1:Lm,kr,n))*OPNSKY
                IF (ijts_sqsc(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqsc(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqsc(1,kr,ntrix(n)))
     &               +SUM(aesqsc(1:Lm,kr,n))
                IF (ijts_sqsc(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqsc(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqsc(2,kr,ntrix(n)))
     &               +SUM(aesqsc(1:Lm,kr,n))*OPNSKY
#ifndef TRACERS_TOMAS
                IF (ijts_sqcb(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
                IF (ijts_sqcb(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               +SUM(aesqcb(1:Lm,kr,n))
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
#else
                qcb_col(kr,n)=0.d0
                do l=1,lm
                  qcb_col(kr,n)=qcb_col(kr,n)+aesqcb(l,kr,n)*
     *                 aesqsc(l,kr,n)
                enddo

                IF (ijts_sqcb(1,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(1,kr,ntrix(n)))
     &               +qcb_col(kr,n)
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)
                IF (ijts_sqcb(2,kr,ntrix(n)) > 0)
     &               taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               =taijs(i,j,ijts_sqcb(2,kr,ntrix(n)))
     &               +qcb_col(kr,n)
     &               /(SUM(aesqsc(1:Lm,kr,n))+1.D-10)*OPNSKY
#endif
              END DO
            END IF
          END SELECT
        END IF
      end do

#endif

#ifdef TRACERS_ON
! this to accumulate daily SUM of optical thickness for
! each active tracer. Will become average in DIAG.f.
! Also saving the aerosol absorption (band 6) 3D, summed over species.
      aerAbs6SaveInst(i,j,:)=0.d0
      do n=1,NTRACE
        if(ntrix(n) > 0) then
          StauL=sum(ttausv(1:LM,n))
          ttausv_sum(i,j,ntrix(n))=ttausv_sum(i,j,ntrix(n))+StauL
          ttausv_sum_cs(i,j,ntrix(n))=ttausv_sum_cs(i,j,ntrix(n))
     &         +StauL*OPNSKY
          aerAbs6SaveInst(i,j,1:lm)=aerAbs6SaveInst(i,j,1:lm) +
     &    (aesqex(1:lm,6,n)-aesqsc(1:lm,6,n))
        endif
      enddo
      IF (adiurn_dust == 1 .or. save3dAOD == 1) THEN
        ttausv_save(i,j,:,:)=0.D0
        DO n=1,NTRACE
          IF (ntrix(n) > 0) THEN
            do k=1,LM
              ttausv_save(i,j,ntrix(n),k)=ttausv_save(i,j,ntrix(n),k)
     &             +ttausv(k,n)
              ttausv_cs_save(i,j,ntrix(n),k)
     &             =ttausv_cs_save(i,j,ntrix(n),k)+ttausv(k,n)*OPNSKY
            end do
          END IF
        END DO
      END IF
#ifdef TRACERS_SPECIAL_Shindell
      do n=1,NTRACE
        if(ntrix(n) > 0) then
          do k=1,LM
            ttausv_ntrace(i,j,n,k)=ttausv(k,n)
          end do
        end if
      end do
#endif
#endif /* TRACERS_ON */

      IF (I.EQ.IWRITE .and. J.EQ.JWRITE) CALL WRITER(6,ITWRITE)
      CSZ2=COSZ2(I,J)
      do L=1,LM
        rad_to_chem(:,L,i,j)=chem_out(L,:)
        rad_to_chem(4,L,i,j)=chem_out(L,4)/CH4X_RADoverCHEM
        do k=1,4
          kliq(L,k,i,j)=kdeliq(L,k) ! save updated flags
        end do
      end do
      if (kradia.gt.0) then  ! rad. forc. model; acc diagn
        do L=1,LM+LM_REQ+1
          AFLX_ST(L,I,J,1)=AFLX_ST(L,I,J,1)+SRUFLB(L)*CSZ2
          AFLX_ST(L,I,J,2)=AFLX_ST(L,I,J,2)+SRDFLB(L)*CSZ2
          AFLX_ST(L,I,J,3)=AFLX_ST(L,I,J,3)+TRUFLB(L)
          AFLX_ST(L,I,J,4)=AFLX_ST(L,I,J,4)+TRDFLB(L)
        end do
        if(kradia.eq.1) then
          tauex6=0. ; tauex5=0. ; tausct=0. ; taugcb=0.
          do L=1,LM
            AFLX_ST(L,I,J,5)=AFLX_ST(L,I,J,5)+1.d2*RHL(L)
            tauex6=tauex6+SRAEXT(L,6)+SRDEXT(L,6)+SRVEXT(L,6)
            tauex5=tauex5+SRAEXT(L,5)+SRDEXT(L,5)+SRVEXT(L,5)
            tausct=tausct+SRASCT(L,6)+SRDSCT(L,6)+SRVSCT(L,6)
            taugcb=taugcb+SRASCT(L,6)*SRAGCB(L,6)+
     +        SRDSCT(L,6)*SRDGCB(L,6)+SRVSCT(L,6)*SRVGCB(L,6)
          end do
          AFLX_ST(LM+1,I,J,5)=AFLX_ST(LM+1,I,J,5)+tauex5
          AFLX_ST(LM+2,I,J,5)=AFLX_ST(LM+2,I,J,5)+tauex6
          AFLX_ST(LM+3,I,J,5)=AFLX_ST(LM+3,I,J,5)+tausct
          AFLX_ST(LM+4,I,J,5)=AFLX_ST(LM+4,I,J,5)+taugcb
          cycle
        end if
        do l=LS1_loc,lm
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     +      (-trfcrl(l)-trhra(l,i,j)) )*nrad*DTsrc*bysha*byMA(l,i,j)
        end do
        do l=lm+1,lm+lm_req
          tchg(l,i,j) = tchg(l,i,j) + ( srfhrl(l)*csz2-srhra(l,i,j) +
     +      (-trfcrl(l)-trhra(l,i,j)) )*nrad*DTsrc*bysha*byaml00(l)
        end do
        cycle
      else if (kradia.lt.0) then ! save i/o data for frc.runs
        fmp_com(i,j) = fmp                  ! input data
        wsoil(i,j) = wearth
        do L=1,LM
          QR(L,I,J) = shl(L)
          CLDinfo(L,1,I,J) = tauwc(L)
          CLDinfo(L,2,I,J) = tauic(L)
          CLDinfo(L,3,I,J) = sizeic(L)  ! sizeic=sizewc currently
        end do
        SRHRA(0,I,J)=SRNFLB(1)*CSZ2      ! output data (for adj frc)
        TRHRA(0,I,J)=-TRNFLB(1)
        do L=1,LM+LM_REQ
          SRHRA(L,I,J)=SRFHRL(L)*CSZ2
          TRHRA(L,I,J)=-TRFCRL(L)
        end do
      end if
C****
C**** Save relevant output in model arrays
C****
C**** (some generalisation and coherence needed in the rad surf type calc)
      FSF(1,I,J)=FSRNFG(1)   !  ocean
      FSF(2,I,J)=FSRNFG(3)   !  ocean ice
      FSF(3,I,J)=FSRNFG(4)   !  land ice
      FSF(4,I,J)=FSRNFG(2)   !  soil
      SRHR(0,I,J)=SRNFLB(1)
      TRHR(0,I,J)=STBO*(
     &      POCEAN*atmocn%GTEMPR(I,J)**4
     &     + POICE*atmice%GTEMPR(I,J)**4
     &     + PLICE*atmgla%GTEMPR(I,J)**4
     &     +PEARTH*atmlnd%GTEMPR(I,J)**4)
     &     -TRNFLB(1)
      TRSURF(1,I,J) = STBO*atmocn%GTEMPR(I,J)**4  !  ocean
      TRSURF(2,I,J) = STBO*atmice%GTEMPR(I,J)**4  !  ocean ice
      TRSURF(3,I,J) = STBO*atmgla%GTEMPR(I,J)**4  !  land ice
      TRSURF(4,I,J) = STBO*atmlnd%GTEMPR(I,J)**4  !  soil
      DO L=1,LM
        SRHR(L,I,J)=SRFHRL(L)
        TRHR(L,I,J)=-TRFCRL(L)
      END DO
      DO LR=1,LM_REQ
        SRHRS(LR,I,J)= SRFHRL(LM+LR)
        TRHRS(LR,I,J)=-TRFCRL(LM+LR)
      END DO
C**** Save fluxes at four levels surface, P0, P1, LTROPO
      SNFS(1,I,J)=SRNFLB(1)     ! Surface
      TNFS(1,I,J)=TRNFLB(1)
      SNFS(2,I,J)=SRNFLB(LM+1)  ! P1
      TNFS(2,I,J)=TRNFLB(LM+1)
      SNFS(3,I,J)=SRNFLB(LM+LM_REQ+1) ! P0 = TOA
      TNFS(3,I,J)=TRNFLB(LM+LM_REQ+1)
      SNFS(4,I,J)=SRNFLB(LTROPO(I,J)) ! LTROPO
      TNFS(4,I,J)=TRNFLB(LTROPO(I,J))

#ifdef SCM
      do L=1,LM
        SRFHRLCOL(L) = SRFHRL(L) * COSZ1(I,J)
        TRFCRLCOL(L) = TRFCRL(L)
      enddo
      SRNFLBBOT = SRNFLB(1) * COSZ1(I,J)           ! Surface
      SRNFLBTOP = SRNFLB(LM+LM_REQ+1) * COSZ1(I,J) ! P0 = TOA
      SRDFLBBOT = SRDFLB(1) * COSZ1(I,J)
      SRDFLBTOP = SRDFLB(LM+LM_REQ+1) * COSZ1(I,J)
      SRUFLBBOT = SRUFLB(1) * COSZ1(I,J)
      SRUFLBTOP = SRUFLB(LM+LM_REQ+1) * COSZ1(I,J)
      TRUFLBTOP = TRUFLB(LM+LM_REQ+1)
      TRDFLBTOP = TRDFLB(LM+LM_REQ+1)
      TRUFLBBOT = TRUFLB(1)
      TRDFLBBOT = TRDFLB(1)
#endif

C****
      TRINCG(I,J)=TRDFLB(1)
      BTMPW(I,J)=BTEMPW-TF
      ALB(I,J,1)=SRNFLB(1)/(SRDFLB(1)+1.D-20)
      ALB(I,J,2)=PLAVIS
      ALB(I,J,3)=PLANIR
      ALB(I,J,4)=ALBVIS
      ALB(I,J,5)=ALBNIR
      ALB(I,J,6)=SRRVIS
      ALB(I,J,7)=SRRNIR
      ALB(I,J,8)=SRAVIS
      ALB(I,J,9)=SRANIR

#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN
        srnflb_save(i,j,1:lm)=srnflb(1:lm)
        trnflb_save(i,j,1:lm)=trnflb(1:lm)
      END IF
#endif
#ifdef mjo_subdd
      swu_avg(I,J)=swu_avg(I,J)+SRUFLB(1)*CSZ2
#endif

      SRDN(I,J) = SRDFLB(1)     ! save total solar flux at surface
C**** SALB(I,J)=ALB(I,J,1)      ! save surface albedo (pointer)
      FSRDIR(I,J)=SRXVIS        ! direct visible solar at surface **coefficient
      SRVISSURF(I,J)=SRDVIS     ! total visible solar at surface
      DIRVIS(I,J)=SRXVIS*SRDVIS  ! direct visible solar at surface
      FSRDIF(I,J)=SRDVIS*(1-SRXVIS) ! diffuse visible solar at surface

      DIRNIR(I,J)=SRXNIR*SRDNIR     ! direct beam nir solar at surface
      DIFNIR(I,J)=SRDNIR*(1-SRXNIR) ! diffuse     nir solar at surface

cdiag write(*,'(a,2i5,6e12.4)')'RAD_DRV: ',
cdiag.    I,J,FSRDIR(I,J),SRVISSURF(I,J),FSRDIF(I,J),
cdiag.        DIRNIR(I,J),SRDNIR,DIFNIR(I,J)
C**** Save clear sky/tropopause diagnostics here
      AIJ(I,J,IJ_CLR_SRINCG)=AIJ(I,J,IJ_CLR_SRINCG)+OPNSKY*
     *     SRDFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_SRNFG)=AIJ(I,J,IJ_CLR_SRNFG)+OPNSKY*
     *     SRNFLB(1)*CSZ2
      AIJ(I,J,IJ_CLR_TRDNG)=AIJ(I,J,IJ_CLR_TRDNG)+OPNSKY*TRHR(0,I,J)
      AIJ(I,J,IJ_CLR_SRUPTOA)=AIJ(I,J,IJ_CLR_SRUPTOA)+OPNSKY*
     *     SRUFLB(LM+LM_REQ+1)*CSZ2
      AIJ(I,J,IJ_CLR_TRUPTOA)=AIJ(I,J,IJ_CLR_TRUPTOA)+OPNSKY*
     *     TRUFLB(LM+LM_REQ+1)
      AIJ(I,J,IJ_CLR_SRNTP)=AIJ(I,J,IJ_CLR_SRNTP)+OPNSKY*
     *     SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_CLR_TRNTP)=AIJ(I,J,IJ_CLR_TRNTP)+OPNSKY*
     *     TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SRNTP)=AIJ(I,J,IJ_SRNTP)+SRNFLB(LTROPO(I,J))*CSZ2
      AIJ(I,J,IJ_TRNTP)=AIJ(I,J,IJ_TRNTP)+TRNFLB(LTROPO(I,J))
      AIJ(I,J,IJ_SISWD)=AIJ(I,J,IJ_SISWD)+POICE*SRDFLB(1)*CSZ2
      AIJ(I,J,IJ_SISWU)=AIJ(I,J,IJ_SISWU)+
     *     POICE*(SRDFLB(1)-FSRNFG(3))*CSZ2

      DO IT=1,NTYPE
         call inc_aj(i,j,it,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *        *CSZ2-TRNFLB(LM+LM_REQ+1))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
         call inc_aj(i,j,it,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *        *CSZ2-TRNFLB(LTROPO(I,J)))*FTYPE(IT,I,J))
      END DO
      call inc_areg(i,j,jr,J_CLRTOA,OPNSKY*(SRNFLB(LM+LM_REQ+1)
     *     *CSZ2-TRNFLB(LM+LM_REQ+1)))
      call inc_areg(i,j,jr,J_CLRTRP,OPNSKY*(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
      call inc_areg(i,j,jr,J_TOTTRP,(SRNFLB(LTROPO(I,J))
     *     *CSZ2-TRNFLB(LTROPO(I,J))))
C**** Save cloud top diagnostics here
      if (CLDCV.le.0.) go to 590
      AIJ(I,J,IJ_CLDTPPR)=AIJ(I,J,IJ_CLDTPPR)+plb(ltopcl+1)
      AIJ(I,J,IJ_CLDTPT)=AIJ(I,J,IJ_CLDTPT)+(tlb(ltopcl+1) - tf)
C**** Save cloud tau=1 related diagnostics here (opt.depth=1 level)
      tauup=0.
      DO L=LM,1,-1
         taucl=tauwc(l)+tauic(l)
         taudn=tauup+taucl
         if (taudn.gt.1.) then
            aij(i,j,ij_cldcv1)=aij(i,j,ij_cldcv1)+1.
            wtlin=(1.-tauup)/taucl
            aij(i,j,ij_cldt1t)=aij(i,j,ij_cldt1t)+( tlb(l+1)-tf +
     +           (tlb(l)-tlb(l+1))*wtlin )
            aij(i,j,ij_cldt1p)=aij(i,j,ij_cldt1p)+( plb(l+1)+
     +           (plb(l)-plb(l+1))*wtlin )
            go to 590
         end if
      end do
 590  continue

      END DO
C****
C**** END OF MAIN LOOP FOR I INDEX
C****

      END DO
C****
C**** END OF MAIN LOOP FOR J INDEX
C****

#ifdef mjo_subdd
       swu_cnt=swu_cnt+1.
#endif

      if(kradia.gt.0) then
         call stopTimer('RADIA()')
         return
      end if
C**** Stop if temperatures were out of range
C**** Now only warning messages are printed for T,Q errors
c     IF(ICKERR.GT.0)
c     &     call stop_model('In Radia: Temperature out of range',11)
c     IF(JCKERR.GT.0)  call stop_model('In Radia: RQT out of range',11)
c     IF(KCKERR.GT.0)  call stop_model('In Radia: Q<0',255)
C**** save all input data to disk if kradia<0
      if (kradia.lt.0) write(iu_rad) itime
     &     ,T,RQT,atmsrf%TsAvg   ! LM+LM_REQ+1+
     &     ,QR,P,CLDinfo,rsi,msi ! LM+1+3*LM+1+1+
!     &     ,(((GTEMPR(k,i,j),k=1,4),i=1,im),j=1,jm) ! (4+)
     &     ,wsoil,atmsrf%wsavg,snowi,atmgla%snow,snowe_com ! 1+1+1+1+1+
     &     ,snoage,fmp_com,flag_dsws,ltropo ! 3+1+.5+.5+
     &     ,fr_snow_rad_ij,mwl,flake ! 2+1+1
C**** output data: really needed only if kradia=2
     &     ,srhra,trhra         ! 2(LM+LM_REQ+1)
     &     ,itime
C****
C**** ACCUMULATE THE RADIATION DIAGNOSTICS
C****
      bydpreq(:) = 1d0/(req_fac_d(:)*pmtop)
      DO 780 J=J_0,J_1
         DO 770 I=I_0,IMAXJ(J)
            do l=1,lm
              call inc_ajl(i,j,l,jl_srhr,SRHR(L,I,J)*COSZ2(I,J))
              call inc_ajl(i,j,l,jl_trcr,TRHR(L,I,J))
            enddo
            CSZ2=COSZ2(I,J)
            JR=JREG(I,J)
            DO LR=1,LM_REQ
              call inc_asjl(i,j,lr,3,bydpreq(lr)*SRHRS(LR,I,J)*CSZ2)
              call inc_asjl(i,j,lr,4,bydpreq(lr)*TRHRS(LR,I,J))
            END DO
            DO KR=1,NDIUPT
            IF (I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
#if (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
              TMP(idd_aot) =SUM(aesqex(1:Lm,6,1:NTRACE))!*OPNSKY
              TMP(idd_aot2) =SUM(aesqsc(1:Lm,6,1:NTRACE))!*OPNSKY
#endif
              TMP(IDD_PALB)=(1.-SNFS(3,I,J)/S0)
              TMP(IDD_GALB)=(1.-ALB(I,J,1))
              TMP(IDD_ABSA)=(SNFS(3,I,J)-SRHR(0,I,J))*CSZ2
              DO INCH=1,NRAD
                IHM=1+(JTIME+INCH-1)*HR_IN_DAY/NDAY
                IH=IHM
                IF(IH.GT.HR_IN_DAY) IH = IH - HR_IN_DAY
                ADIURN(IDXB(:),KR,IH)=ADIURN(IDXB(:),KR,IH)+TMP(IDXB(:))
#ifndef NO_HDIURN
                IHM = IHM+(DATE-1)*HR_IN_DAY
                IF(IHM.LE.HR_IN_MONTH) THEN
                  HDIURN(IDXB(:),KR,IHM)=HDIURN(IDXB(:),KR,IHM)+
     &                 TMP(IDXB(:))
                ENDIF
#endif
              ENDDO
            END IF
            END DO

      DO IT=1,NTYPE
         call inc_aj(I,J,IT,J_SRINCP0,(S0*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRINCG ,(SRHR(0,I,J)*CSZ2/(ALB(I,J,1)+1.D
     *        -20))*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_BRTEMP ,BTMPW(I,J) *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRINCG ,TRINCG(I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HSURF  ,-(TNFS(3,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP0 ,-TNFS(3,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_TRNFP1 ,-TNFS(2,I,J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_SRNFP1 ,SNFS(2,I,J)*CSZ2*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,J_HATM   ,-(TNFS(2,I,J)-TNFS(1,I,J))
     *        *FTYPE(IT,I,J))
#ifdef HEALY_LM_DIAGS
         call inc_aj(I,J,IT,j_vtau ,10.*-20*VTAULAT(J)*FTYPE(IT,I,J))
         call inc_aj(I,J,IT,j_ghg ,10.*ghg_totforc*FTYPE(IT,I,J))
#endif

      END DO
C**** Note: confusing because the types for radiation are a subset
         call inc_aj(I,J,ITOCEAN,J_SRNFG,(FSF(1,I,J)*CSZ2)*FOCEAN(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITLAKE ,J_SRNFG,(FSF(1,I,J)*CSZ2)* FLAKE(I,J)
     *        *(1.-RSI(I,J)))
         call inc_aj(I,J,ITEARTH,J_SRNFG,(FSF(4,I,J)*CSZ2)*FEARTH(I,J))
         call inc_aj(I,J,ITLANDI,J_SRNFG,(FSF(3,I,J)*CSZ2)* FLICE(I,J))
         call inc_aj(I,J,ITOICE ,J_SRNFG,(FSF(2,I,J)*CSZ2)*FOCEAN(I,J)
     *        *RSI(I,J))
         call inc_aj(I,J,ITLKICE,J_SRNFG,(FSF(2,I,J)*CSZ2)* FLAKE(I,J)
     *        *RSI(I,J))
C****
         call inc_areg(I,J,JR,J_SRINCP0,(S0*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP0 ,(SNFS(3,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRNFP1 ,(SNFS(2,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_SRINCG ,(SRHR(0,I,J)*CSZ2/
     *        (ALB(I,J,1)+1.D-20)))
         call inc_areg(I,J,JR,J_HATM,-(TNFS(2,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_SRNFG,(SRHR(0,I,J)*CSZ2))
         call inc_areg(I,J,JR,J_HSURF,-(TNFS(3,I,J)-TNFS(1,I,J)))
         call inc_areg(I,J,JR,J_BRTEMP, BTMPW(I,J))
         call inc_areg(I,J,JR,J_TRINCG,TRINCG(I,J))
         call inc_areg(I,J,JR,J_TRNFP0,-TNFS(3,I,J))
         call inc_areg(I,J,JR,J_TRNFP1,-TNFS(2,I,J))
         DO K=2,9
           JK=AJ_ALB_INDS(K-1) ! accumulate 8 radiation diags.
           DO IT=1,NTYPE
             call inc_aj(I,J,IT,JK,(S0*CSZ2)*ALB(I,J,K)*FTYPE(IT,I,J))
           END DO
           call inc_areg(I,J,JR,JK,(S0*CSZ2)*ALB(I,J,K))
         END DO
         AIJ(I,J,IJ_SRINCG) =AIJ(I,J,IJ_SRINCG) +(SRHR(0,I,J)*CSZ2/
     /        (ALB(I,J,1)+1.D-20))
         AIJ(I,J,IJ_SRNFG) =AIJ(I,J,IJ_SRNFG) +(SRHR(0,I,J)*CSZ2)
         AIJ(I,J,IJ_BTMPW) =AIJ(I,J,IJ_BTMPW) +BTMPW(I,J)
         AIJ(I,J,IJ_SRREF) =AIJ(I,J,IJ_SRREF) +S0*CSZ2*ALB(I,J,2)
         AIJ(I,J,IJ_SRVIS) =AIJ(I,J,IJ_SRVIS) +S0*CSZ2*ALB(I,J,4)
         AIJ(I,J,IJ_TRNFP0)=AIJ(I,J,IJ_TRNFP0)-TNFS(3,I,J)
         AIJ(I,J,IJ_SRNFP0)=AIJ(I,J,IJ_SRNFP0)+(SNFS(3,I,J)*CSZ2)
         AIJ(I,J,IJ_RNFP1) =AIJ(I,J,IJ_RNFP1) +(SNFS(2,I,J)*CSZ2
     *                                         -TNFS(2,I,J))
         AIJ(I,J,ij_srvdir)=AIJ(I,J,ij_srvdir)
     &        + FSRDIR(I,J)*SRVISSURF(I,J)
         AIJ(I,J,IJ_SRVISSURF)=AIJ(I,J,IJ_SRVISSURF)+SRVISSURF(I,J)
#ifdef mjo_subdd
         OLR_acc(I,J) = OLR_acc(I,J) -TNFS(3,I,J)
#endif
C**** CRF diags if required
         if (moddrf .ne. 0) go to 770
         if (cloud_rad_forc > 0) then
c    CRF diagnostics
           AIJ(I,J,IJ_SWCRF)=AIJ(I,J,IJ_SWCRF)+
     +          (SNFS(3,I,J)-SNFSCRF(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF)=AIJ(I,J,IJ_LWCRF)-
     -          (TNFS(3,I,J)-TNFSCRF(I,J))
c    CRF diagnostics without aerosols and Ox
           AIJ(I,J,IJ_SWCRF2)=AIJ(I,J,IJ_SWCRF2)+
     +          (SNFS(3,I,J)-SNFSCRF2(I,J))*CSZ2
           AIJ(I,J,IJ_LWCRF2)=AIJ(I,J,IJ_LWCRF2)-
     -          (TNFS(3,I,J)-TNFSCRF2(I,J))
         end if

C**** AERRF diags if required
         if (aer_rad_forc > 0) then
           do N=1,8
             AIJ(I,J,IJ_SWAERRF+N-1)=AIJ(I,J,IJ_SWAERRF+N-1)+
     *            (SNFS(3,I,J)-SNFSAERRF(N,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERRF+N-1)=AIJ(I,J,IJ_LWAERRF+N-1)-
     *            (TNFS(3,I,J)-TNFSAERRF(N,I,J))
             AIJ(I,J,IJ_SWAERSRF+N-1)=AIJ(I,J,IJ_SWAERSRF+N-1)+
     *            (SNFS(1,I,J)-SNFSAERRF(N+8,I,J))*CSZ2
             AIJ(I,J,IJ_LWAERSRF+N-1)=AIJ(I,J,IJ_LWAERSRF+N-1)-
     *            (TNFS(1,I,J)-TNFSAERRF(N+8,I,J))
           end do
           AIJ(I,J,IJ_SWAERRFNT)=AIJ(I,J,IJ_SWAERRFNT)+
     *          (SNFS(3,I,J)-SNFSAERRF(17,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERRFNT)=AIJ(I,J,IJ_LWAERRFNT)-
     *          (TNFS(3,I,J)-TNFSAERRF(17,I,J))
           AIJ(I,J,IJ_SWAERSRFNT)=AIJ(I,J,IJ_SWAERSRFNT)+
     *          (SNFS(1,I,J)-SNFSAERRF(18,I,J))*CSZ2
           AIJ(I,J,IJ_LWAERSRFNT)=AIJ(I,J,IJ_LWAERSRFNT)-
     *          (TNFS(1,I,J)-TNFSAERRF(18,I,J))
         end if

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_DUST) ||\
    (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
C**** Generic diagnostics for radiative forcing calculations
C**** Depending on whether tracers radiative interaction is turned on,
C**** diagnostic sign changes (for aerosols)
         rsign_aer=1. ; rsign_chem=-1.
         if (rad_interact_aer > 0) rsign_aer=-1.
C**** define SNFS/TNFS level (TOA/TROPO) for calculating forcing
         LFRC=3                 ! TOA
         if (rad_forc_lev.gt.0) LFRC=4 ! TROPOPAUSE
#ifdef  TRACERS_AMP
         IF (AMP_DIAG_FC == 2) THEN
            NTRACE = nmodes
         ELSE
            NTRACE = 1
            NTRIX(1) = 1
         ENDIF
#endif /* TRACERS_AMP */
#ifdef  TRACERS_TOMAS
         IF (TOMAS_DIAG_FC == 2) THEN
            NTRACE = icomp-2
         ELSE
            NTRACE = 1
            NTRIX(1) = 1
         ENDIF
#endif /* TRACERS_TOMAS */
         if (ntrace > 0) then
#ifdef BC_ALB
      if (ijts_alb(1).gt.0)
     * TAIJS(I,J,ijts_alb(1))=TAIJS(I,J,ijts_alb(1))
     *   + 100.d0*(ALBNBC(I,J)-ALB(I,J,1))
      if (ijts_alb(2).gt.0)
     & taijs(i,j,ijts_alb(2))
     &     =taijs(i,j,ijts_alb(2))
     &         +(SNFS(3,I,J)-NFSNBC(I,J))*CSZ2
#endif /* BC_ALB */
#ifdef TRACERS_AEROSOLS_Koch
c          snfst0(:,:,i,j)=0.D0
c          tnfst0(:,:,i,j)=0.D0
#endif /* TRACERS_AEROSOLS_Koch */
c     ..........
c     accumulation of forcings for tracers for which ntrace fields are
c     defined
c     ..........
           set_clayilli=.FALSE.
           set_claykaol=.FALSE.
           set_claysmec=.FALSE.
           set_claycalc=.FALSE.
           set_clayquar=.FALSE.
           do n=1,ntrace
             IF (ntrix(n) > 0) THEN
               SELECT CASE (trname(ntrix(n)))
               CASE ('Clay')
                 n1=n-nrad_clay+1
c shortwave forcing (TOA or TROPO) of Clay sub size classes
                 if (ijts_fcsub(1,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(1,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(1,ntrix(n),n1))
     &                +rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2
c longwave forcing  (TOA or TROPO) of Clay size sub classes
                 if (ijts_fcsub(2,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(2,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(2,ntrix(n),n1))
     &                -rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
c shortwave forcing (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(5,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(5,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(5,ntrix(n),n1))
     &                +rsign_aer*(snfst(2,n,i,j)-snfs(lfrc,i,j))*csz2
     &                *(1.D0-cfrac(i,j))
c longwave forcing  (TOA or TROPO) clear sky of Clay sub size classes
                 if (ijts_fcsub(6,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(6,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(6,ntrix(n),n1))
     &                -rsign_aer*(tnfst(2,n,i,j)-tnfs(lfrc,i,j))
     &                *(1.D0-cfrac(i,j))
c shortwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(3,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(3,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(3,ntrix(n),n1))
     &                +rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))*csz2
c longwave forcing at surface (if required) of Clay sub size classes
                 if (ijts_fcsub(4,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(4,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(4,ntrix(n),n1))
     &                -rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))
c shortwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(7,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(7,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(7,ntrix(n),n1))
     &                +rsign_aer*(snfst(1,n,i,j)-snfs(1,i,j))*csz2
     &                *(1.D0-cfrac(i,j))
c longwave forcing at surface clear sky (if required) of Clay sub size classes
                 if (ijts_fcsub(8,ntrix(n),n1) > 0)
     &                taijs(i,j,ijts_fcsub(8,ntrix(n),n1))
     &                =taijs(i,j,ijts_fcsub(8,ntrix(n),n1))
     &                -rsign_aer*(tnfst(1,n,i,j)-tnfs(1,i,j))
     &                *(1.D0-cfrac(i,j))
               CASE DEFAULT
                 SELECT CASE (trname(ntrix(n)))
                 CASE ('seasalt2')
                   CYCLE
                 CASE ('ClayIlli')
                   IF (set_clayilli) CYCLE
                   set_clayilli=.TRUE.
                 CASE ('ClayKaol')
                   IF (set_claykaol) CYCLE
                   set_claykaol=.TRUE.
                 CASE ('ClaySmec')
                   IF (set_claysmec) CYCLE
                   set_claysmec=.TRUE.
                 CASE ('ClayCalc')
                   IF (set_claycalc) CYCLE
                   set_claycalc=.TRUE.
                 CASE ('ClayQuar')
                   IF (set_clayquar) CYCLE
                   set_clayquar=.TRUE.
                 END SELECT
c shortwave forcing (TOA or TROPO)
                 if (ijts_fc(1,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(1,ntrix(n)))
     &                =taijs(i,j,ijts_fc(1,ntrix(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
c longwave forcing  (TOA or TROPO)
                 if (ijts_fc(2,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(2,ntrix(n)))
     &                =taijs(i,j,ijts_fc(2,ntrix(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
c shortwave forcing (TOA or TROPO) clear sky
                 if (ijts_fc(5,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(5,ntrix(n)))
     &                =taijs(i,j,ijts_fc(5,ntrix(n)))
     &                +rsign_aer*(SNFST(2,N,I,J)-SNFS(LFRC,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing  (TOA or TROPO) clear sky
                 if (ijts_fc(6,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(6,ntrix(n)))
     &                =taijs(i,j,ijts_fc(6,ntrix(n)))
     &                -rsign_aer*(TNFST(2,N,I,J)-TNFS(LFRC,I,J))
     &                *(1.d0-CFRAC(I,J))
c shortwave forcing at surface (if required)
                 if (ijts_fc(3,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(3,ntrix(n)))
     &                =taijs(i,j,ijts_fc(3,ntrix(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
c longwave forcing at surface (if required)
                 if (ijts_fc(4,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(4,ntrix(n)))
     &                =taijs(i,j,ijts_fc(4,ntrix(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
c shortwave forcing at surface clear sky (if required)
                 if (ijts_fc(7,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(7,ntrix(n)))
     &                =taijs(i,j,ijts_fc(7,ntrix(n)))
     &                +rsign_aer*(SNFST(1,N,I,J)-SNFS(1,I,J))*CSZ2
     &                *(1.d0-CFRAC(I,J))
c longwave forcing at surface clear sky (if required)
                 if (ijts_fc(8,ntrix(n)).gt.0)
     &                taijs(i,j,ijts_fc(8,ntrix(n)))
     &                =taijs(i,j,ijts_fc(8,ntrix(n)))
     &                -rsign_aer*(TNFST(1,N,I,J)-TNFS(1,I,J))
     &                *(1.d0-CFRAC(I,J))
               END SELECT
#ifdef  TRACERS_AMP
         IF (AMP_DIAG_FC == 2) THEN
         ELSE
         NTRIX(1)=  n_N_AKK_1
         ENDIF
#endif /* TRACERS_AMP */
#ifdef  TRACERS_TOMAS
!I don't know why this is used..
         IF (TOMAS_DIAG_FC == 2) THEN
         ELSE
         NTRIX(1)=  n_ANUM(1)
         ENDIF
#endif /* TRACERS_TOMAS */
#ifdef TRACERS_AEROSOLS_Koch
c              SNFST0(1,ntrix(n),I,J)=SNFST0(1,ntrix(n),I,J)
c    &              +rsign_aer*(SNFST(2,n,I,J)-SNFS(LFRC,I,J))*CSZ2
c              SNFST0(2,ntrix(n),I,J)=SNFST0(2,ntrix(n),I,J)
c    &              +rsign_aer*(SNFST(1,n,I,J)-SNFS(1,I,J))*CSZ2
c              TNFST0(1,ntrix(n),I,J)=TNFST0(1,ntrix(n),I,J)
c    &              -rsign_aer*(TNFST(2,n,I,J)-TNFS(LFRC,I,J))
c              TNFST0(2,ntrix(n),I,J)=TNFST0(2,ntrix(n),I,J)
c    &              -rsign_aer*(TNFST(1,n,I,J)-TNFS(1,I,J))
#endif /* TRACERS_AEROSOLS_Koch */
             END IF   ! ntrix(n)>0
           end do     ! n=1,ntrace
         end if       ! ntrace>0

c ..........
c accumulation of forcings for special case ozone (ntrace fields
c not defined) Warning: indicies used differently, since we don't
c need CS or Surface, but are doing both TOA and Ltropo:
c ..........
         if(n_Ox > 0) then ! ------ main Ox tracer -------
c shortwave forcing at tropopause
           if (ijts_fc(1,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(1,n_Ox))=taijs(i,j,ijts_fc(1,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(1,I,J)-SNFS(4,I,J))*CSZ2
c longwave forcing at tropopause
           if (ijts_fc(2,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(2,n_Ox))=taijs(i,j,ijts_fc(2,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(1,I,J)-TNFS(4,I,J))
c shortwave forcing at TOA
           if (ijts_fc(3,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(3,n_Ox))=taijs(i,j,ijts_fc(3,n_Ox))
     &     +rsign_chem*(SNFST_o3ref(2,I,J)-SNFS(3,I,J))*CSZ2
c longwave forcing at TOA
           if (ijts_fc(4,n_Ox).gt.0)
     &     taijs(i,j,ijts_fc(4,n_Ox))=taijs(i,j,ijts_fc(4,n_Ox))
     &     -rsign_chem*(TNFST_o3ref(2,I,J)-TNFS(3,I,J))
         endif
#ifdef AUXILIARY_OX_RADF
c shortwave forcing at tropopause
         if (ijts_auxfc(1)>0)
     &   taijs(i,j,ijts_auxfc(1))=taijs(i,j,ijts_auxfc(1))
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_o3ref(3,I,J))*CSZ2
c longwave forcing at tropopause
         if (ijts_auxfc(2)>0)
     &   taijs(i,j,ijts_auxfc(2))=taijs(i,j,ijts_auxfc(2))
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_o3ref(3,I,J))
c shortwave forcing at TOA
         if (ijts_auxfc(3)>0)
     &   taijs(i,j,ijts_auxfc(3))=taijs(i,j,ijts_auxfc(3))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_o3ref(4,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_auxfc(4)>0)
     &   taijs(i,j,ijts_auxfc(4))=taijs(i,j,ijts_auxfc(4))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_o3ref(4,I,J))
#endif /* AUXILIARY_OX_RADF */
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
                      ! ------ diag stratOx tracer -------
! note for now for this diag, there is a failsafe that stops model
! if clim_interact_chem .le. 0 when the below would be wrong:
c shortwave forcing at tropopause
         if (ijts_fc(1,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(1,n_stratOx))=taijs(i,j,ijts_fc(1,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(1,I,J)-SNFST_stratOx(1,I,J))*CSZ2
c longwave forcing at tropopause
         if (ijts_fc(2,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(2,n_stratOx))=taijs(i,j,ijts_fc(2,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(1,I,J)-TNFST_stratOx(1,I,J))
c shortwave forcing at TOA
         if (ijts_fc(3,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(3,n_stratOx))=taijs(i,j,ijts_fc(3,n_stratOx))
     &   +rsign_chem*(SNFST_o3ref(2,I,J)-SNFST_stratOx(2,I,J))*CSZ2
c longwave forcing at TOA
         if (ijts_fc(4,n_stratOx).gt.0)
     &   taijs(i,j,ijts_fc(4,n_stratOx))=taijs(i,j,ijts_fc(4,n_stratOx))
     &   -rsign_chem*(TNFST_o3ref(2,I,J)-TNFST_stratOx(2,I,J))
#endif /* SHINDELL_STRAT_EXTRA && ACCMIP_LIKE_DIAGS*/
#endif /* any of various tracer groups defined */

#ifdef ACCMIP_LIKE_DIAGS
         do nf=1,4 ! CH4, N2O, CFC11, and CFC12:
c shortwave GHG forcing at TOA
           if(ij_fcghg(1,nf).gt.0)aij(i,j,ij_fcghg(1,nf))=
     &     aij(i,j,ij_fcghg(1,nf))+(SNFS(3,I,J)-SNFS_ghg(nf,I,J))
     &     *CSZ2
c longwave GHG forcing at TOA
           if(ij_fcghg(2,nf).gt.0)aij(i,j,ij_fcghg(2,nf))=
     &     aij(i,j,ij_fcghg(2,nf))+(TNFS_ghg(nf,I,J)-TNFS(3,I,J))
         enddo
#endif /* ACCMIP_LIKE_DIAGS */

  770    CONTINUE
  780    CONTINUE

#ifdef mjo_subdd
         OLR_cnt = OLR_cnt + 1.
#endif

      DO J=J_0,J_1
      DO I=I_0,I_1
      DO L=1,LM
        AIJL(i,j,l,IJL_RC)=AIJL(i,j,l,IJL_RC)+
     &       (SRHR(L,I,J)*COSZ2(I,J)+TRHR(L,I,J))
#ifdef mjo_subdd
       SWHR(I,J,L)=SWHR(I,J,L)+
     *             SRHR(L,I,J)*COSZ2(I,J)*bysha*byMA(L,I,J)
       LWHR(I,J,L)=LWHR(I,J,L)+TRHR(L,I,J)*bysha*byMA(L,I,J)
#endif
      END DO
      END DO
      END DO
#ifdef mjo_subdd
      SWHR_cnt=SWHR_cnt+1
      LWHR_cnt=LWHR_cnt+1
#endif

C****
C**** Update radiative equilibrium temperatures
C****
      DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO LR=1,LM_REQ
            RQT(LR,I,J)=RQT(LR,I,J)+(SRHRS(LR,I,J)*COSZ2(I,J)
     &           +TRHRS(LR,I,J))*NRAD*DTsrc*bysha*byaml00(lr+lm)
          END DO
        END DO
      END DO
C****
C**** Update other temperatures every physics time step
C****
  900 DO J=J_0,J_1
        DO I=I_0,IMAXJ(J)
          DO L=1,LM
            T(I,J,L)=T(I,J,L)+(SRHR(L,I,J)*COSZ1(I,J)+TRHR(L,I,J))*
     *           DTsrc*bysha*byMA(l,i,j)/PK(L,I,J)
#ifdef SCM
            dTradlw(L)=TRHR(L,I,J)*DTsrc*bysha*byMA(l,i,j)/PK(L,I,J)
            dTradsw(L)=SRHR(L,I,J)*COSZ1(I,J)
     &                *DTsrc*bysha*byMA(L,I,J)/PK(L,I,J)
#endif
          END DO
          AIJ(I,J,IJ_SRINCP0)=AIJ(I,J,IJ_SRINCP0)+(S0*COSZ1(I,J))
        END DO
      END DO

C**** daily diagnostics
      IH=1+modelEclock%getHour()
      IHM = IH+(modelEclock%getDate()-1)*24
      DO KR=1,NDIUPT
        I = IJDD(1,KR)
        J = IJDD(2,KR)
        IF ((J >= J_0) .AND. (J <= J_1) .AND.
     &      (I >= I_0) .AND. (I <= I_1)) THEN
          ADIURN(IDD_ISW,KR,IH)=ADIURN(IDD_ISW,KR,IH)+S0*COSZ1(I,J)
#ifndef NO_HDIURN
          HDIURN(IDD_ISW,KR,IHM)=HDIURN(IDD_ISW,KR,IHM)+S0*COSZ1(I,J)
#endif
        ENDIF
      ENDDO

      call stopTimer('RADIA()')
      RETURN
      END SUBROUTINE RADIA

      SUBROUTINE RESET_SURF_FLUXES(I,J,ITYPE_OLD,ITYPE_NEW,FTYPE_ORIG,
     *     FTYPE_NOW)
!@sum set incident solar and upward thermal fluxes appropriately
!@+   as fractions change to conserve energy, prevent restart problems
!@auth Gavin Schmidt
      use rad_com, only : fsf,trsurf
      implicit none
!@var itype_old, itype_new indices for the old type turning to new type
      integer, intent(in) :: i,j,itype_old,itype_new
!@var ftype_orig, ftype_now original and current fracs of the 'new' type
      real*8, intent(in) :: ftype_orig, ftype_now
      real*8 :: delf ! change in fraction from old to new

      if (( (ITYPE_OLD==1 .and. ITYPE_NEW==2)
     *  .or.(ITYPE_OLD==2 .and. ITYPE_NEW==1))
     *  .and. (FTYPE_NOW .le. 0. .or. FTYPE_NOW .gt. 1.)) then
        write (6,*) ' RESET_SURF_FLUXES:
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW = ',
     *    I, J, ITYPE_OLD, ITYPE_NEW, FTYPE_ORIG, FTYPE_NOW
        call stop_model('RESET_SURF_FLUXES: INCORRECT RESET',255)
      end if

      delf = FTYPE_NOW-FTYPE_ORIG
C**** Constrain fsf_1*ftype_1+fsf_2*ftype_2 to be constant
      FSF(ITYPE_NEW,I,J)=(FSF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     FSF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

C**** Same for upward thermal
      TRSURF(ITYPE_NEW,I,J)=(TRSURF(ITYPE_NEW,I,J)*FTYPE_ORIG+
     *     TRSURF(ITYPE_OLD,I,J)*DELF)/FTYPE_NOW

      RETURN
      END SUBROUTINE RESET_SURF_FLUXES

      SUBROUTINE GHGHST(iu)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy

      use domain_decomp_atm, only : write_parallel
      USE RADPAR, only : nghg,ghgyr1,ghgyr2,ghgam
      USE RAD_COM, only : ghg_yr
      IMPLICIT NONE
      INTEGER :: iu,n,k,nhead=4,iyr
      CHARACTER*80 title
      character(len=300) :: out_line

      write(out_line,*)  ! print header lines and first data line
      call write_parallel(trim(out_line),unit=6)
      do n=1,nhead+1
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
        nhead=5
      end if

!**** find range of table: ghgyr1 - ghgyr2
      read(title,*) ghgyr1
      do ; read(iu,'(a)',end=20) title ; end do
   20 read(title,*) ghgyr2
      rewind iu  !   position to data lines
      do n=1,nhead ; read(iu,'(a)') ; end do

      allocate (ghgam(nghg,ghgyr2-ghgyr1+1))
      do n=1,ghgyr2-ghgyr1+1
        read(iu,*) iyr,(ghgam(k,n),k=1,nghg)
        do k=1,nghg ! replace -999. by reasonable numbers
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
        if(ghg_yr>0 .and. abs(ghg_yr-iyr).le.1) then
          write(out_line,'(i5,6f10.4)') iyr,(ghgam(k,n),k=1,nghg)
          call write_parallel(trim(out_line),unit=6)
        endif
      end do
      write(out_line,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      call write_parallel(trim(out_line),unit=6)
      return
      end SUBROUTINE GHGHST

#if defined(CUBED_SPHERE)
      subroutine read_qma (iu,plb)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
      use domain_decomp_atm, only : write_parallel
      use rad_com, only : dH2O,jma=>jm_dh2o,lat_dh2o
      use resolution, only : lm
      use constant, only : radian
      use TimeConstants_mod, only: DAYS_PER_YEAR
      implicit none
      integer, parameter:: lma=24
      integer m,iu,j,l,ll,ldn(lm),lup(lm)
      real*8 :: plb(lm+1)
      real*4 pb(0:lma+1),h2o(jma,0:lma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,dh,fracl
      character(len=300) :: out_line

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (lat_dh2o(j),j=1,jma)
      lat_dh2o(:) = lat_dh2o(:)*radian

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) vertically
        do j=1,jma
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+h2o(j,ll)*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/DAYS_PER_YEAR !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine read_qma

      subroutine lat_interp_qma (rlat,lev,mon,dh2o_interp)
!@sum  interpolate CH4->H2O production rates in latitude
!@auth R. Ruedy
      use rad_com, only : jma=>jm_dh2o,xlat=>lat_dh2o,dh2o
      implicit none
      real*8 :: rlat ! input latitude (radians)
      integer :: lev,mon ! input level, month
      real*8 :: dh2o_interp  ! output
      real*8 w1,w2
      integer :: j1,j2

C**** Interpolate (extrapolate) horizontally
      j2 = 2+(jma-1)*(rlat-xlat(1))/(xlat(jma)-xlat(1)) ! first guess
      j2 = min(max(2,j2),jma)
      j1 = j2-1
      if(rlat.gt.xlat(j2)) then ! j guess was too low
         do while (j2.lt.jma .and. rlat.gt.xlat(j2))
          j2 = j2+1
         end do
         j1 = j2-1
      elseif(rlat.lt.xlat(j1)) then ! j guess was too high
         do while (j1.gt.1 .and. rlat.lt.xlat(j1))
          j1 = j1-1
         end do
         j2 = j1+1
      endif
C**** coeff. for latitudinal linear inter/extrapolation
      w1 = (xlat(j2)-rlat)/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
      if(w1.gt.1.) w1=.5+.5*w1
      if(w1.lt.0.) w1=.5*w1
      w2 = 1.-w1
      dh2o_interp = w1*dh2o(j1,lev,mon)+w2*dh2o(j2,lev,mon)
      return
      end subroutine lat_interp_qma

#endif /* CUBED_SPHERE */

      subroutine getqma (iu,dglat,plb,dh2o,lm,jm)
!@sum  reads H2O production rates induced by CH4 (Tim Hall)
!@auth R. Ruedy
      use domain_decomp_atm, only : grid,getDomainBounds,write_parallel
      use TimeConstants_mod, only: DAYS_PER_YEAR
      implicit none
      integer, parameter:: jma=18,lma=24
      integer m,iu,jm,lm,j,j1,j2,l,ll,ldn(lm),lup(lm)
      real*8 PLB(lm+1),dH2O(grid%j_strt_halo:grid%j_stop_halo,lm,12)
     &     ,dglat(jm)
      real*4 pb(0:lma+1),h2o(jma,0:lma),xlat(jma),z(lma),dz(0:lma)
      character*100 title
      real*4 pdn,pup,w1,w2,dh,fracl
      integer :: j_0,j_1
      character(len=300) :: out_line
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

C**** read headers/latitudes
      read(iu,'(a)') title
      write(out_line,'(''0'',a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
      write(out_line,'(1x,a100)') title
      call write_parallel(trim(out_line),unit=6)
      read(iu,'(a)') title
c      write(6,'(1x,a100)') title
      read(title(10:100),*) (xlat(j),j=1,jma)

C**** read heights z(km) and data (kg/km^3/year)
      do m=1,12
        read(iu,'(a)') title
        write(out_line,'(1x,a100)') title
        call write_parallel(trim(out_line),unit=6)
        do l=lma,1,-1
          read(iu,'(a)') title
c          write(6,'(1x,a100)') title
          read(title,*) z(l),(H2O(j,l),j=1,jma)
        end do
        do j=1,jma
          h2o(j,0) = 0.
        end do

C**** Find edge heights and pressures
        dz(0) = 0.
        dz(1) = z(2)-z(1)
        do l=2,lma-1
           dz(l)=.5*(z(l+1)-z(l-1))
        end do
        dz(lma) = z(lma)-z(lma-1)

        pb(0) = plb(1)
        do l=1,lma
           Pb(l)=1000.*10.**(-(z(l)-.5*dz(l))/16.)
        end do
C**** extend both systems vertically to p=0
        pb(lma+1)=0.
        plb(lm+1)=0.

C**** Interpolate vertical resolution to model layers
        ldn(:) = 0
        do l=1,lm
          do while (pb(ldn(l)+1).ge.plb(l) .and. ldn(l).lt.lma)
            ldn(l)=ldn(l)+1
          end do
          lup(l)=ldn(l)
          do while (pb(lup(l)+1).gt.plb(l+1) .and. lup(l).lt.lma)
            lup(l)=lup(l)+1
          end do
        end do

C**** Interpolate (extrapolate) horizontally and vertically
        j2=2
        do j=j_0,j_1
C**** coeff. for latitudinal linear inter/extrapolation
          do while (j2.lt.jma .and. dglat(j).gt.xlat(j2))
            j2 = j2+1
          end do
          j1 = j2-1
          w1 = (xlat(j2)-dglat(j))/(xlat(j2)-xlat(j1))
C**** for extrapolations, only use half the slope
          if(w1.gt.1.) w1=.5+.5*w1
          if(w1.lt.0.) w1=.5*w1
          w2 = 1.-w1
          do l=1,lm
            dh = 0.
            pdn = plb(l)
            if (lup(l).gt.0) then
              do ll=ldn(l),lup(l)
                pup = max(REAL(pb(ll+1),KIND=8),plb(l+1))
                fracl= (pdn-pup)/(pb(ll)-pb(ll+1))
                dh = dh+(w1*h2o(j1,ll)+w2*h2o(j2,ll))*fracl*dz(ll)
                pdn = pup
              end do
            end if
            dh2o(j,l,m) = 1.d-6*dh/1.74d0/DAYS_PER_YEAR !->(kg/m^2/ppm_CH4/day)
          end do
        end do
      end do
      return
      end subroutine getqma

