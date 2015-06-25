E001tr.R GISS Model E with sample tracers               jal 10/01

WARNING: This rundeck is being used as a template to exercise normally
unused bits of code. It should not therefore be used for production
runs. Please look at the other templates for standard
configurations. If example tracer code is required, the
tracer-relevant bits used here (pre-processing options, extra TR object
modules, and input files), should be copied to a more functional rundeck.

E001tr: ModelE1 (3.0) (based on B402A with sample tracers)
 Air mass, SF6, RN222, CO2, 14CO2, CFC-11, CH4, N2O, linearizedO3, SF6_c, Water
23 lyrs, top at 85km - 1979 atmosphere/ocean
gravity wave drag;     uses dry convection (rather than turbulence)
Sdrag: weak linear strat. drag in top layer, near poles down to 20 mb
       lost ang.mom is added in below 150 mb
sealevel pressure filter applied every hour
6-band oice albedo; Hogstrom(1984) pbl drag
Note: Many of these choices may be changed using the PARAMETERs below.
Further Note: This has been extended to 23 layers and some plug and play
options are different to provide a wider test suite for testing. (i.e.
4th order momentum, dummy ice dynamics and DRYCNV instead of ATRUB).

Preprocessor Options
#define TRACERS_ON          ! include tracers code
#define TRACERS_SPECIAL_Lerner ! also activate TRACER_SPECIAL_Lerner in Obj.modules !!
#define TRACERS_WATER      ! include water tracers code
#define NEW_IO
End Preprocessor Options

Run Options
STACKSIZE=131072

Object modules: (in order of decreasing priority)
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL23                      ! vertical resolution is 23 layers -> 0.1mb
MODEL_COM GEOM_B IO_DRV             ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN4TH          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q and tracer gases
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT TRDIAG        ! tracer diagnostic printout
! use next line if #define TRACERS_SPECIAL_Lerner
TRACER_SPECIAL_Lerner               ! routines called when TRACERS_SPECIAL_Lerner is activated
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                      ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices ATURB or DRYCNV
! ATURB                               ! turbulence in whole atmosphere
DRYCNV                              ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DUM  ! or ICEDYN_DRV ICEDYN  ! dynamic ice modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_RES_M                          ! diagnostics (resolution dependent)
      FFT72                         ! utilities
POUT                                ! post-processing output

Components:
tracers MPI_Support shared dd2d

Data input files:
    ! the first 4 files are specific to prescribed ocean runs
AIC=AIC.RES_M23.D771201.nc           ! initial conditions (atm.)
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc ! initial conditions (ground)
! prescr. climatological ocean (1 yr of data)
OSST=OST4X5.B.1975-84avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1975-84avg.Hadl1.1.nc
! prescr. climatological sea ice
SICE=SICE4X5.B.1975-84avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1975-84avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1975-84avg.Hadl1.1.nc
    ! if the prescr. ocean varies from year to year use instead:
! OSST=OST4X5.B.1950.M02.Hadl1.1  ! ocean data   Feb 1950 - 1999
! SICE=SICE4X5.B.1950.M02.Hadl1.1 ! ocean data   Feb 1950 - 1999
    ! the next 3 files are specific to q-flux ocean runs
! AIC=E001/1JAN1956.rsfE001.O250D      ! AIC/OHT made by aux/mkOTSPEC.E001.M250D
! OHT=E001/OTSPEC.E001.M250D.1951-1955 ! horizontal ocean heat transport
OCNML=Z1O.B4X5.cor.nc                ! mixed layer depth (use for post processing)
    ! files needed for all models
CDN=CD4X500S.ext.nc
VEG=V72X46.1.cor2_no_crops.ext.nc    ! vegetation fractions
SOIL=S4X50093.ext.nc
TOPO=Z72X46N.cor4_nocasp.nc   ! soil/topography bdy.conds
REG=REG4X5                        ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
ZVAR=ZVAR4X5.nc         ! topographic variation for gwdrag
#include "rad_input_files"
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
#include "rad_72x46_input_files"
TOP_INDEX=top_index_72x46_a.ij.ext.nc
MSU_wts=MSU.RSS.weights.data
GLMELT=GLMELT_4X5.OCN.nc   ! glacial melt distribution
CO2_IC=CO2ijl_IC_Jan1_scale334_M23  !wofsy+B140TQaM9
CO2_FOS_FUEL=CO2_sources/gcm_data/CO2FOS_MRL_4X5
CO2_FERT=CO2_sources/gcm_data/CO2fert01_4X5
CO2_REGROWTH=CO2_sources/gcm_data/CO2_Nforest_4X5
CO2_LAND_USE=CO2_sources/gcm_data/CO2DEF_HOU_4X5
CO2_VEG=CO2_sources/gcm_data/CO2VEG_MON_4X5          ! Monthly source
CO2_OCEAN=CO2_sources/gcm_data/CO2_4X5_Ocean_flux02  ! Monthly source
14CO2_IC_DATA=workshop.14co2                         ! for 14CO2 Oct. 1963
LINOZ_TABLE=O3_linoz_coeff                !linoz coefficients for stratosphere
LO3_Trop_loss=linoz/LOx_IJ_M23_trop_L11   !Troposphere ozone chemical loss (Harvard
LO3_Trop_prod=linoz/POx_IJ_M23_trop_L11   !Troposphere ozone chemical production (Harvard)
LINOZ_Dep_vel=linoz/O3dv_IJ.bin           !Deposition velocity for O3 (Harvard)
N2O_TABLE=N2Oloss.table                     ! Stratosphere tracer forcing
CFC11_TABLE=F11loss.table                   ! Stratosphere tracer forcing
CH4_TABLE=CH4chem.table                     ! Stratosphere tracer forcing
CH4_TROP_FRQ=CLIM.RUN.OHCH4.FRQ      !tropo loss frequency table (9 layers, n-grid)
N2O_IC=N2O_Shindell_Jan9293_M23             !initial conditions
CH4_IC=Wofsy_data_CH4       !wofsy jl initial conditions
CH4_ANIMALS=methane/gcm_data/CH4ANIMLS_4X5      ! Annual
CH4_COALMINE=methane/gcm_data/CH4COAL_4X5       ! Annual
CH4_GASLEAK=methane/gcm_data/CH4GASLEAK_4X5     ! Annual
CH4_GASVENT=methane/gcm_data/CH4GASVENT_4X5     ! Annual
CH4_CITYDUMP=methane/gcm_data/CH4MSW_4X5        ! Annual
CH4_SOIL_ABS=methane/gcm_data/CH4SOILABS_4X5    ! Annual
CH4_TERMITES=methane/gcm_data/CH4TRMITE_4X5     ! Annual
CH4_COALBURN=methane/gcm_data/COAL_BURN_BY_POP84_4X5   ! Annual
CH4_BURN=methane/gcm_data/CH4BURN_4X5           ! Monthly
CH4_RICE=methane/gcm_data/CH4RICEC_4X5          ! Monthly
CH4_WETL=methane/gcm_data/CH4WETL+TUNDRA_4X5    ! Monthly

Label and Namelist:
E001tr (ModelE1 (3.0) based on B402A, uses dry adiab. adjustment; tracers)
R=00BG/B
DTFIX=90
&&PARAMETERS
! parameters set for prescribed ocean runs:
KOCEAN=0        ! ocn is prescribed
Kvflxo=0        ! save VFLXO (daily) if ocn prescribed

! parameters usually not changed when switching to q-flux ocean:

X_SDRAG=.00025,.000025  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.1     ! lin. sdrag above .1mb (top 2 layers) except near poles
PP_SDRAG=1.    ! lin. sdrag above 1.mb near poles (top 4 layers)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=1000.     ! maximum wind velocity in sdrag; default=200 when GW drag not used
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000035 !the default is 15d-6
PCONPEN=500.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=1.     ! Shear drag coefficient
CMTN=0.25

xCDpbl=1.
cond_scheme=1   ! 2 = more elaborate conduction scheme (GHY, Nancy Kiang)
 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.61      ! U00ice up  => nethtz0 down (alb down) goals: nethtz0=0 (ann.
U00wtrX=1.259   ! U00wtrX up => nethtz0 up   (alb down)           global mean)
! HRMAX=550.    ! HRMAX up   => nethtz0 down (alb up  )        plan.alb 30%

RWCLDOX=1.5  !  wtr cld particle size *3/2 over ocean
RICLDX=.3333 !  ice cld particle size * 1(at 0mb)->1/3(at 1000mb)

H2ObyCH4=1.     ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the atmospheric composition
! if set to 0, the current (day/) year is used: transient run
s0_yr=1979
s0_day=182
ghg_yr=1979
ghg_day=182
volc_yr=1979
volc_day=182
aero_yr=1979
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1979
dalbsnX=.015
o3_yr=-1979

! parameters that control the Shapiro filter
DT_XUfilter=180. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=180. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

! parameters that may have to be changed in emergencies:
LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DT=180.,             ! from default: DTsrc=3600.,
NIsurf=4,            ! number of surface time steps

! parameters that affect at most diagn. output:
Ndisk=24        ! use =240 on halem
SUBDD=' '       ! save SLP at sub-daily frequency
NSUBDD=0        ! saving sub-daily diags 12hrly
KCOPY=2         ! saving acc + rsf
isccp_diags=0   ! use =0 to save cpu time
nda5d=1         ! use =7 to save cpu time
nda5s=1         ! use =7 to save cpu time

itime_tr0=8016,8016,8016,8016,10920,15696,8016,8016,8016,8016
to_volume_MixRat=1,1,1,1,1,1,1,1,1,1   ! for tracer printout
nstrtc=12                    ! Number of layers for Prather stratosphere chemistry (LM=23)

&&END_PARAMETERS

 &INPUTZ
   YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! IYEAR1=YEARI (default)
   YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=13*0,
   ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/
