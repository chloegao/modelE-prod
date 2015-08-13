E4TdusM20.R GISS Model E  2004 modelE    jan perlwitz 10/27/2009
 Template for setting up simulations using the E4M20 model with dust tracers

E4TdusM20.R: E4M20 + dust tracers

E4M20: modelE as frozen (or not yet) in July 2009 without gravity wave drag
modelE4 4x5 hor. grid with 20 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 2000
ocean data: prescribed, 1993-2002 climatology  (see OSST/SICE)
uses turbulence scheme, simple strat.drag (not grav.wave drag)
time steps: dynamics 7.5 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define USE_ENT                  ! include dynamic vegetation model
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! default wet deposition or water tracers 
#define TRACERS_DUST             ! include dust tracers
#define TRACERS_DUST_Silt4       ! include 4th silt size class
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define NO_HDIURN                ! switches off hdiurn diagnostics
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL20 STRAT_DUM          ! vertical resolution is 20 layers -> 0.1mb
DIAG_RES_M FFT72 
MODEL_COM GEOM_B IORSF              ! model variables and geometry
MODELE                              ! Main and model overhead
ALLOC_DRV                           ! allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
QUS3D                               ! advection of Q and tracer gases
TRACER_COM TRACERS_DRV              ! configurable tracer code 
TRACERS                             ! generic tracer code 
TRDRYDEP                            ! dry deposition of tracers 
TRDUST_COM TRDUST TRDUST_DRV        ! dust tracer specific code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV ! + component giss_LSM: land surface and soils
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + component Ent: new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
ATURB_E1                            ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
POUT_netcdf                         ! post-processing output

Components:
tracers Ent shared MPI_Support solvers giss_LSM

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB
OPTS_giss_LSM = USE_ENT=YES

Data input files:
    ! resolution dependent files
    ! start up from restart file of earlier run
! AIC=1DECxxxx.rsfEyyyy ! initial conditions (atm./ground), no GIC, ISTART=8
    ! or start up from observed conditions
AIC=AIC.RES_M20A.D771201.nc          ! initial conditions (atm.)      needs GIC, ISTART=2
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc   ! initial conditions (ground)
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1993-2002avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1993-2002avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1993-2002avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1993-2002avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1993-2002avg.Hadl1.1.nc
! For q-flux ocean, replace lines above by the next 2 lines & set KOCEAN=1, ISTART=8 
!! AIC=1JAN1961.rsfE4M20.MXL65m   ! = end of preliminary run with KOCEAN=0,Kvflxo=1
!! OHT=OTSPEC.E4M20.MXL65m.1956-1960 ! ocean horizontal heat transport
OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (needed for post processing)
TOPO=Z72X46N.cor4_nocasp.nc       ! topography
SOIL=S4X50093.ext.nc              ! soil bdy.conds
! VEG=V72X46.1.cor2   ! or:       ! vegetation fractions  (sum=1), need crops_yr=-1
VEG=V72X46.1.cor2_no_crops.ext.nc  ! veg. fractions
CROPS=CROPS2007_72X46N.cor4_nocasp.nc ! crops history
soil_textures=soil_textures_top30cm
SOILCARB_global=soilcarb_top30cm_4x5.nc
CDN=CD4X500S.ext.nc               ! surf.drag coefficient
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
TOP_INDEX=top_index_72x46_a.ij.ext.nc  ! only used if #define DO_TOPMODEL_RUNOFF
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution

    ! resolution independent files
RADN1=sgpgxg.table8               ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02_hdr
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr      ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
MSU_wts=MSU.RSS.weights.data

#include "rad_72x46_input_files"

!------- Needed for dry deposition ---------
VEGTYPE=chem_files/vegtype.global
OLSON=chem_files/drydep.table
DRYCOEFF=chem_files/drydep.coef
LAI01=chem_files/lai01.global
LAI02=chem_files/lai02.global
LAI03=chem_files/lai03.global
LAI04=chem_files/lai04.global
LAI05=chem_files/lai05.global
LAI06=chem_files/lai06.global
LAI07=chem_files/lai07.global
LAI08=chem_files/lai08.global
LAI09=chem_files/lai09.global
LAI10=chem_files/lai10.global
LAI11=chem_files/lai11.global
LAI12=chem_files/lai12.global

! sulfate+black carbon files:
!  MADAER=1 (default) needs:
!TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
!TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
!TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
!TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
! updated aerosols need MADAER=3
TAero_SUL=SUL_Koch2008_kg_m2_72x46x20_1890-2000h.nc
TAero_SSA=SSA_Koch2008_kg_m2_72x46x20h.nc
TAero_NIT=NIT_Bauer2008_kg_m2_72x46x20_1890-2000h.nc
TAero_OCA=OCA_Koch2008_kg_m2_72x46x20_1890-2000h.nc
TAero_BCA=BCA_Koch2008_kg_m2_72x46x20_1890-2000h.nc
TAero_BCB=BCB_Koch2008_kg_m2_72x46x20_1890-2000h.nc

! files for dust tracers
VTRSH=vtr-mod-o0.mean-pb       ! threshold wind speeds for old emission scheme
FRCLAY=claygcm-f               ! clay fraction for old emission scheme
FRSILT=siltgcm-f               ! silt fraction for old emission scheme
DRYHR=text5hr-f                ! threshold of dry hours for old emission scheme
ERS=ERS1_1993_MONTHLY.72x46.threshold-13 ! ERS data
DSRC=Ginoux2001_source_VegMask_72x46     ! preferred dust sources
                                         ! optimized use: prefDustSources=0
                                         ! (only choice so far)
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities
dust_bin1=DUST_bin1_2000_new.nc  ! for AEROCOM
dust_bin2=DUST_bin2_2000_new.nc  ! for AEROCOM
dust_bin3=DUST_bin3_2000_new.nc  ! for AEROCOM
dust_bin4=DUST_bin4_2000_new.nc  ! for AEROCOM

Label and Namelist:
E4TdusM20 (ModelE1 4x5, 20 lyrs, 2000 atm, 1993-2002 clim ocn, with passive dust
DTFIX=300

&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0 ! ocn is prescribed
!! KOCEAN=1 ! ocn is computed
Kvflxo=0 ! set =1 after spinup to prepare for q-flux run (edit "I")

variable_lk=1 ! 1: let lakes grow or shrink in horizontal extent
wsn_max=2.   ! restrict snow depth to 2 m-h2o (if 0. snow depth is NOT restricted)

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP

PTLISO=15.  ! press(mb) above which rad. assumes isothermal layers

xCDpbl=1.
cond_scheme=2    ! more elaborate conduction scheme (GHY, Nancy Kiang)

! tuning param.: this setting works for 1850; use U00wtrX=1.28 for 1979
 
U00a=.60    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=2.3    ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.62      ! U00ice+.01 =>dBal=1.5,dPl.alb=-.9%   goals:Bal=0,plan.alb=30%
U00wtrX=1.29    ! U00wtrX+.01=>dBal=0.7,dPl.alb=-.25%  Bal=glb.ann NetHt at z0
! HRMAX=500.    ! not needed unless do_blU00=1, HRMAX up => nethtz0 down (alb up)

H2OstratX=1.
H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=-1 no land ice fixup, 1 Lacis' scheme)
KSOLAR=2
cloud_rad_forc=1 ! 1: calculate cloud radiative forcing

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=2000 ! if -1, crops in VEG-file is used! =???? , also change OSST,SICE
s0_yr=2000                                      ! =???? , also change OSST,SICE
s0_day=182
ghg_yr=2000                                     ! =???? , also change OSST,SICE
ghg_day=182
volc_yr=-1  ! 1850-1999 mean strat.aeros        ! =???? , also change OSST,SICE
volc_day=182
aero_yr=2000                                    ! =???? , also change OSST,SICE
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (was .0036)
albsn_yr=2000                                   ! =???? , also change OSST,SICE
dalbsnX=.024
o3_yr=2000                                      ! =???? , also change OSST,SICE
CO2X=1.
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

variable_orb_par=0
orb_par_year_bp=-50  !  BP i.e. 1950-orb_par_year AD = 2000 AD

!--------- general aerosol parameters-------------
aer_rad_forc=0     ! 1: calculate aerosol radiative forcing
rad_forc_lev=1     ! 0: for TOA, 1: for tropopause for rad forcing diags.
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: additional radiation diagnostics
diag_wetdep=1      ! 1: additional wet deposition diagnostics
to_conc=0,1,1,1,1,1 ! 1: taijln diags as concentration, 0: as mixing ratio (default), needs to be exactly set according to the number and order of tracers, if used

!--------- sulfate and carbon aerosol parameters -----
madaer=3           ! 1: default sulfate and carbon aerosol, 3: updated aerosols

!--------- dust aerosol parameters----------------
imDust=0              ! 0: PDF emission scheme, 1: AEROCOM, 2: legacy emission scheme
adiurn_dust=0         ! 1: daily dust diagnostics at certain grid points
prefDustSources=0     ! 0: Ginoux 2001 w/ vegetation mask
                      ! 1-4: No files and optimization available yet
                      ! >4: Free choice of emis. parameters
!fracClayPDFscheme=1. ! Frac. clay emis, only effective for prefDustSources > 4
!fracSiltPDFscheme=1. ! Frac. silt emis, only effective for prefDustSources > 4
                      ! set internally for 0-4

!-------------------------------------------------

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

DTsrc=1800.     ! cannot be changed after a run has been started
! parameters that may have to be changed in emergencies:
DT=450.
NIsurf=1        ! increase as layer 1 gets thinner

! parameters that affect at most diagn. output:
Ndisk=1920
SUBDD=' '       ! no sub-daily frequency diags
NSUBDD=0        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
KCOPY=2         ! saving acc + rsf   =3 to also save "oda"-files
isccp_diags=1   ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13        ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48         ! to get daily energy history use nda4=24*3600/DTsrc
nssw=2          ! until diurnal diagn. are fixed, nssw should be even

&&END_PARAMETERS

 &INPUTZ
   YEARI=1999,MONTHI=12,DATEI=1,HOURI=0, IYEAR1=1875 ! or earlier
   YEARE=2005,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1999,MONTHE=12,DATEE=1,HOURE=1,
/
