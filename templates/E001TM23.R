E001TM23.R GISS Model E                                 gas 06/00

E001TM23: Differences from original TRACERS_SPECIAL_Lerner include:
  Tropospheric ozone scheme from Loretta Mickley with deposition in layer 1
  New tracer (SF6_c) with constant source
  parameter NSTRTC replaces ltropo
  Better initial distribution of CFC-11
  Air tracer removed
Based on E001M23
 SF6, Rn222, CO2, N2O, CFC-11, 14CO2, CH4, linearizedO3, SF6_c

Preprocessor Options
#define TRACERS_ON                  ! include tracers code
#define TRACERS_SPECIAL_Lerner
End Preprocessor Options

Object modules: (in order of decreasing priority)
Atm72x46                   ! horizontal resolution is 72x46 -> 4x5deg
AtmL23                      ! vertical resolution is 23 layers -> 0.1mb
DIAG_RES_M FFT72   
MODEL_COM GEOM_B IORSF              ! model variables and geometry
TRIDIAG                             ! tridiagonal matrix solver
MODELE                              ! Main and model overhead
                                    ! parameter database
              ALLOC_DRV             ! domain decomposition, allocate global distributed arrays
ATMDYN_COM ATMDYN MOMEN2ND          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
STRATDYN STRAT_DIAG                 ! strospheric dynamics (incl. gw drag)
QUS_COM QUSDEF QUS_DRV              ! advection of tracers
TQUS_DRV                            ! advection of Q
TRACER_COM TRACERS_DRV              ! configurable tracer code
TRACERS                             ! generic tracer code
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
! use next line if #define TRACERS_SPECIAL_Lerner
TRACER_SPECIAL_Lerner               ! routines called when TRACERS_SPECIAL_Lerner is activated
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES                              ! surface calculation and fluxes
GHY_COM GHY_DRV GHY GHY_H           ! land surface and soils
VEG_DRV VEG_COM VEGETATION          ! vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
! pick exactly one of the next 2 choices: ATURB or DRYCNV
ATURB                               ! turbulence in whole atmosphere
! DRYCNV                            ! drycnv
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV                 ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
OCEAN OCNML                         ! ocean modules
SNOW_DRV SNOW                       ! snow model
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS ALBEDO READ_AERO          ! radiation and albedo
DIAG_COM DIAG DEFACC DIAG_PRT       ! diagnostics
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
                                    ! utilities
POUT                                ! post-processing output

Components:
tracers MPI_Support shared

Data input files:
AIC=1OCT1951.rsfE464jalTM23_noTrace
GIC=GIC.E046D3M20A.1DEC1955.ext_1.nc
OSST=OST4X5.B.1876-85avg.Hadl1.1.nc
OSST_eom=OST4X5.B.1876-85avg.Hadl1.1.nc
SICE=SICE4X5.B.1876-85avg.Hadl1.1.nc
SICE_eom=SICE4X5.B.1876-85avg.Hadl1.1.nc
ZSIFAC=SICE4X5.B.1876-85avg.Hadl1.1.nc
OCNML=Z1O.B4X5.cor.nc   ! needed only for postprocessing
!                                             (end of section 1 of data input files)
TOPO=Z72X46N.cor4_nocasp.nc
SOIL=S4X50093.ext.nc              ! bdy.cond
VEG=V72X46.1.cor2_no_crops.ext.nc
CROPS=CROPS2007_72X46N.cor4_nocasp.nc
CDN=CD4X500S.ext.nc
REG=REG4X5           ! special regions-diag
RVR=RD_modelE_M.nc                ! river direction file
NAMERVR=RD_modelE_M.names.txt     ! named river outlets
TOP_INDEX=top_index_72x46_a.ij.ext.nc
ZVAR=ZVAR4X5.nc         ! topographic variation for gwdrag
!                                             (end of section 2 of data input files)
RADN1=sgpgxg.table8    ! rad.tables
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800              ! rad.tables and history files
RADN4=LWCorrTables33k              ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_Ma_2008
RADN3=miescatpar.abcdv2
! RADNA,RADNB are no longer used
TAero_PRE=dec2003_PRE_Koch_kg_m2_ChinSEA_Liao_1850 ! pre-industr trop. aerosols
TAero_SUI=sep2003_SUI_Koch_kg_m2_72x46x9_1875-1990 ! industrial sulfates
TAero_OCI=sep2003_OCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial organic carbons
TAero_BCI=sep2003_BCI_Koch_kg_m2_72x46x9_1875-1990 ! industrial black carbons
RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-1999.Apr02
RADN8=cloud.epsilon4.72x46
RADN9=solar.lean02.ann.uvflux_hdr     ! need KSOLAR=2
RADNE=topcld.trscat8
ISCCP=ISCCP.tautables
#include "rad_72x46_input_files"
GHG=GHG.Mar2004.txt
dH2O=dH2O_by_CH4_monthly
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
CFCic_Lerner=CFCic_Lerner_72x46             ! 1st layer initial distribution
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
E001TM23 (Ejal1T+New ozone; add SF6_c; NSTRTC replaces ltropo)
R=00BG/B

&&PARAMETERS
X_SDRAG=.0005,.00005  ! used for lin. sdrag above P_SDRAG mb
C_SDRAG=0.     ! no constant sdrag
P_SDRAG=.01     ! lin. sdrag above p_sdrag mb (top layer for M23) except near poles
PP_SDRAG=4.6   ! lin. sdrag above  pp_sdrag mb near poles (top 5 layers for M23)
ANG_SDRAG=1    ! if =1: sdrag conserves ang mom.
WMAX=1000.     ! maximum wind velocity in sdrag; default=200 when GW drag not used

PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000030 !the default is 15d-6
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=1.      ! Shear drag coefficient
CMTN=0.25      ! default is 0.5
CDEF=1.5       ! deformation drag coefficient

KOCEAN=0
Kvflxo=0        ! saving VFLXO (daily)

variable_lk=1
wsn_max=2.

 
U00a=.55    ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; then tune this to get rad.balance
! U00a,U00b replace the U00 parameters below - U00ice/U00wtrX are kept only for the _E1 version
U00ice=.62   ! tune this first to get reas.alb/cldcvr (range: .4-.6), then
u00wtrx = 1.29  ! 1880 conditions

cond_scheme=2  ! more elaborate conduction scheme (GHY, Nancy Kiang)
H2ObyCH4=1.  ! activates strat.H2O generated by CH4
KSIALB=0        ! 6-band albedo (Hansen) (=1 A.Lacis orig. 6-band alb)
KSOLAR=2

! parameters that control the Shapiro filter
DT_XUfilter=450. ! Shapiro filter on U in E-W direction; usually same as DT (below)
DT_XVfilter=450. ! Shapiro filter on V in E-W direction; usually same as DT (below)
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

LMCM=16              ! max level of moist convection
XCDNST=300.,10000.   ! strat. gw drag parameters
DTsrc = 1800.        ! half-hour physics time step (default: DTsrc=3600.)
DT=450.,             ! dynamic time step
NIsurf=1,            ! number of surface time steps

! The next 4 parameters correspond to the half-hour physics time step
NDISK = 480,
ndaa = 13,
nda5k = 13,
nda4  = 48,

NSUBDD=0        ! saving sub-daily diags
KCOPY=2         ! saving acc + rsf
isccp_diags=0

to_volume_MixRat=1,1,1,1,1,1,1,1,1   ! for tracer printout
itime_tr0=30624,30624,30624,30624,39360,31344,30624,30624,30624
nstrtc=12                    ! Number of layers for Prather stratosphere chemistry (LM=23)
! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
crops_yr=1880 ! if -1, crops in VEG-file is used
s0_yr=1880
ghg_yr=1880
ghg_day=182
s0_day=182
volc_yr=1880
volc_day=182
aero_yr=1880
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.0036    ! include 2nd indirect effect
albsn_yr=1880
dalbsnX=.024
o3_yr=-1880
&&END_PARAMETERS

 &INPUTZ
   YEARI=1951,MONTHI=10,DATEI=1,HOURI=0,
   IYEAR1=1950,
   YEARE=1951,MONTHE=12, DATEE=1,HOURE=0,     KDIAG=0,2,2,9*0,9,
   YEARE=1951,MONTHE=11, DATEE=1,HOURE=0,     KDIAG=0,2,2,9*0,9,
   ISTART=5,IRANDI=0, YEARE=1951,MONTHE=10,DATEE=1,HOURE=1,

/


