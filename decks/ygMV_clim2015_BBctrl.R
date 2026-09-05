ygMV_clim2015_BBctrl.R GISS ModelE Lat-Lon Atmosphere Model, climatological ocn/atm MATRIX tracers

ygMV_clim2015_BBctrl: E6TmatrixF40clim but year 2000 settings
               (e.g. 9-year averages centered around nominal date)

Lat-lon: 2x2.5 degree horizontal resolution
F40: 40 vertical layers with standard hybrid coordinate, top at .1 mb
Atmospheric composition for year 2015
Ocean climatology prescribed from years 2003-2017, CMIP6
Uses turbulence scheme (no dry conv), grav.wave drag
Time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
Filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)

Preprocessor Options
#define STDHYB                   ! standard hybrid vertical coordinate
#define ATM_LAYERING L40         ! 40 layers, top at .1 mb
#define NEW_IO                   ! new I/O (netcdf) on
#define IRRIGATION_ON
#define MODIS_LAI
#define NEW_BCdalbsn
!---> generic tracers code start
#define TRAC_ADV_CPU             ! timing index for tracer advection on
#define TRACERS_ON               ! include tracers code
#define TRACERS_WATER            ! wet deposition and water tracer
#define TRACERS_DRYDEP           ! default dry deposition
#define TRDIAG_WETDEPO           ! additional wet deposition diags for tracers
#define TRACERS_AIR
!  OFF #define CALCULATE_LIGHTNING ! Calculate lightning flash rates when NOx is not needed
!  OFF #define AUTOTUNE_LIGHTNING  ! Automatically generate lightning tuning parameters (present-day only)
!<--- generic tracers code end
!---> chemistry start
#define TRACERS_SPECIAL_Shindell    ! includes drew's chemical tracers
!  OFF #define AUXILIARY_OX_RADF ! radf diags for climatology or tracer Ozone
#define TRACERS_TERP                ! include terpenes in gas-phase chemistry
#define BIOGENIC_EMISSIONS       ! turns on interactive isoprene emissions
!  OFF #define WATER_MISC_GRND_CH4_SRC ! adds lake, ocean, misc. ground sources for CH4
!  OFF #define INTERACTIVE_WETLANDS_CH4 ! turns on interactive CH4 wetland source
!<--- chemistry end
!---> MATRIX start
#define TRACERS_AMP
#define TRACERS_AMP_M9
#define TRACERS_AEROSOLS_VBS     ! Volatility-Basis Set
!<--- MATRIX end
#define BC_ALB                    !optional tracer BC affects snow albedo
#define CLD_AER_CDNC              !aerosol-cloud interactions
#define BLK_2MOM                  !aerosol-cloud interactions
!  OFF #define NUDGE_ON                 ! nudge the meteorology
#define CACHED_SUBDD
End Preprocessor Options

Object modules:
     ! resolution-specific source codes
Atm144x90                           ! horizontal resolution is 144x90 -> 2x2.5deg
AtmLayering                         ! vertical resolution
FFT144                              ! Fast Fourier Transform

IO_DRV                              ! new i/o

     ! GISS dynamics with gravity wave drag
ATMDYN MOMEN2ND                     ! atmospheric dynamics
QUS_DRV QUS3D                       ! advection of Q/tracers
STRATDYN STRAT_DIAG                 ! stratospheric dynamics (incl. gw drag)

    ! lat-lon grid specific source codes
AtmRes
GEOM_B                              ! model geometry
DIAG_ZONAL GCDIAGb                  ! grid-dependent code for lat-circle diags
DIAG_PRT POUT                       ! diagn/post-processing output
MODEL_COM                           ! calendar, timing variables
MODELE_DRV                          ! ModelE cap
MODELE                              ! initialization and main loop
ATM_COM                             ! main atmospheric variables
ATM_DRV                             ! driver for atmosphere-grid components
ATMDYN_COM                          ! atmospheric dynamics
ATM_UTILS                           ! utilities for some atmospheric quantities
QUS_COM QUSDEF                      ! T/Q moments, 1D QUS
CLOUDS2 CLOUDS2_DRV CLOUDS_COM      ! clouds modules
SURFACE SURFACE_LANDICE FLUXES      ! surface calculation and fluxes
GHY_COM GHY_DRV    ! + giss_LSM     ! land surface and soils + snow model
VEG_DRV                             ! vegetation
! VEG_COM VEGETATION                ! old vegetation
ENT_DRV  ENT_COM   ! + Ent          ! new vegetation
PBL_COM PBL_DRV PBL                 ! atmospheric pbl
IRRIGMOD                            ! irrigation module
ATURB                               ! turbulence in whole atmosphere
LAKES_COM LAKES                     ! lake modules
SEAICE SEAICE_DRV                   ! seaice modules
LANDICE LANDICE_COM LANDICE_DRV     ! land ice modules
ICEDYN_DRV ICEDYN                   ! ice dynamics modules
RAD_COM RAD_DRV RADIATION           ! radiation modules
RAD_UTILS GHGMOD ALBEDO READ_AERO ocalbedo ! radiation and albedo
DIAG_COM DIAG DEFACC                ! diagnostics
OCN_DRV                             ! driver for ocean-grid components
OCEAN OCNML                         ! ocean modules

! Codes common to atmospheric tracer sets
TRACER_COM                          ! configurable tracer code
TRACERS_DRV |-O0|                   ! O0 speeds compilation and no difference for this file
initTracerGriddedData		    ! grid independent info
initTracerMetadata		    ! misc initialization that used to be mixed in with metadata
TRACERS                             ! generic tracer code
TRDRYDEP                            ! dry deposition of tracers
TRDIAG_COM TRACER_PRT               ! tracer diagnostic printout
lightning                           ! lightning module
MiscTracersMetadata
sharedTracersMetadata
ATMCOL_COMtr                        ! basic atmospheric state seen in tracer_3Dsource
! ---TRACER SPECIFIC CODES----------
TRACERS_SPECIAL_Shindell            ! routines specific to drew's 15-tracers
TRCHEM_Shindell_COM                 ! Drew Shindell's tracers common
TRCHEM_calc                         ! chemical reaction calculations
TRCHEM_init                         ! chemistry initialization, I/O
TRCHEM_family                       ! tracer family chemistry
TRCHEM_fastj2                       ! used for trop+strat chem version
TRCHEM_master                       ! trop chem "driver"/strat prescrioption
BIOGENIC_EMISSIONS                  ! old N.Unger interactive isoprene
ShindellTracersMetadata

! ----------------------------------
TRDUST_COM TRDUST TRDUST_DRV               ! dust tracer specific code
TRACERS_AEROSOLS_SEASALT                   ! seasalt
TRACERS_AEROSOLS_Koch_e4                   ! BC/OC/sulfate
TRACERS_AEROSOLS_VBS                ! Volatility-Basis Set
! Aerosol MicroPhysics
TRAMP_drv        | $(EXTENDED_SOURCE) |
TRAMP_mod        | $(EXTENDED_SOURCE) |
TRAMP_actv       | $(EXTENDED_SOURCE) |
TRAMP_diam       | $(EXTENDED_SOURCE) |
TRAMP_nomicrophysics | $(EXTENDED_SOURCE) |
TRAMP_subs       | $(EXTENDED_SOURCE) |
TRAMP_coag       | $(EXTENDED_SOURCE) |
TRAMP_depv       | $(EXTENDED_SOURCE) |
TRAMP_param_GISS | $(EXTENDED_SOURCE) |
TRAMP_config
TRAMP_dicrete    | $(EXTENDED_SOURCE) |
TRAMP_init       | $(EXTENDED_SOURCE) |
TRAMP_quad       | $(EXTENDED_SOURCE) |
TRAMP_matrix     | $(EXTENDED_SOURCE) |
TRAMP_setup      | $(EXTENDED_SOURCE) |
TRAMP_npf        | $(EXTENDED_SOURCE) |
TRAMP_rad        | $(EXTENDED_SOURCE) |
! Enable next two lines when using EQSAM Thermodynamics
TRACER_drv_eqsam | $(EXTENDED_SOURCE) |    ! Nitrate aerosol
TRAMP_eqsam_v03d                           ! EQSAM module for inorganic aerosol thermodynamic equilibrium
! Enable next four lines when using ISORROPIA Thermodynamics
!TRACER_drv_isorropia | $(EXTENDED_SOURCE) |
!TRAMP_isocom2
!TRAMP_isofwd2
!TRAMP_isorev2
! ----------------------------------
AmpTracersMetadata
VbsTracersMetadata
TRDIAG                              ! new i/o
SUBDD
! NUDGE
CLD_AEROSOLS_Menon_MBLK_MAT_E29q BLK_DRV ! aerosol-cloud interactions
CLD_AER_CDNC                        ! aerosol-cloud interactions wrapper

Components:
shared MPI_Support solvers giss_LSM 
dd2d
tracers
Ent

Component Options:
OPTS_Ent = ONLINE=YES PS_MODEL=FBB PFT_MODEL=ENT 
OPTS_dd2d = NC_IO=PNETCDF

Data input files:
    ! start from the restart file of an earlier run ...                 ISTART=8
! AIC=1....rsfE... ! initial conditions, no GIC needed, use
!! AIC=1JAN1961.rsfE4F40.MXL65m   ! end of run with KOCEAN=0

    ! start from observed conditions AIC(,OIC), model ground data GIC   ISTART=2
! AIC=AIC.RES_F40.D771201.nc      ! observed initial conditions for F40 1977/12/01
! AIC=AIC_144x90_DEC01_L96.nc     ! observed initial conditions for F96 1977/12/01
AIC=NCARIC.144x90.D7712010_ext.nc ! AIC for automatic relayering to model vertical grid
GIC=GIC.144X90.DEC01.1.ext_1.nc   ! initial ground conditions

OSST=OST_144x90.2003-2017avg.HadISST.nc     ! climatological ocean temperature
SICE=SICE_144x90.2003-2017avg.HadISST.nc    ! climatological sea ice cover
ZSIFAC=ZSIfac_144x90.2003-2017avg.HadISST.nc  ! climatological sea ice thickness
TOPO=Z2HX2fromZ1QX1N.BS1.nc               ! ocean fraction and surface topography
!! Q-flux ocean: use the next line instead, set KOCEAN=1
!! OHT=OTSPEC.E4F40.MXL65m.1956-1960         ! ocean horizontal heat transports
!! OCNML=Z1O.B144x90.nc                      ! mixed layer depth for Q-flux model
RVR=RD_Fd.nc             ! river direction file
NAMERVR=RD_Fd.names.txt  ! named river outlets

CDN=CD144X90.ext.nc
VEG=V144x90_EntMM16_lc_max_trimmed_scaled_nocrops.ext.nc
LAIMAX=V144x90_EntMM16_lai_max_trimmed_scaled_ext.nc
HITEent=V144x90_EntMM16_height_trimmed_scaled_ext.nc
LAI=V144x90_EntMM16_lai_trimmed_scaled_ext.nc
CROPS=CROPS_and_pastures_Pongratz_to_Hurtt_144X90N_nocasp.nc
IRRIG=Irrig144x90_1848to2100_FixedFuture_v3.nc
SOIL=S144X900098M.ext.nc
TOP_INDEX=top_index_144x90_a.ij.ext.nc
ZVAR=ZVAR2X25A.nc             ! topographic variation for gwdrag
! probably need these (should convert to 144x90)
soil_textures=soil_textures_top30cm_2x2.5
SOILCARB_global=soilcarb_top30cm_2x2.5.nc
GLMELT=GLMELT_144X90_gas.OCN.nc
RADN1=sgpgxg.table8                           ! rad.tables and history files
RADN2=LWTables33k_lowH2O_CO2_O3_planck_1-800  ! rad.tables and history files
RADN4=LWCorrTables33k                         ! rad.tables and history files
RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies H2O continuum table
! other available H2O continuum tables:
!    RADN5=H2Ocont_Ma_2000
!    RADN5=H2Ocont_Ma_2004
!    RADN5=H2Ocont_Roberts
!    RADN5=H2Ocont_MT_CKD  ! Mlawer/Tobin_Clough/Kneizys/Davies
RADN3=miescatpar.abcdv2

RH_QG_Mie=oct2003.relhum.nr.Q633G633.table
RADN7=STRATAER.VOL.1850-2014_CMIP6_hdr  ! needs MADVOL=2
RADN8=cloud.epsilon4.72x46
!RADN9=solar.lean2015.ann1610-2014.nc ! need KSOLAR=2
RADN9=solar.CMIP6official.ann1850-2299_with_E3_fastJ.nc ! need KSOLAR=2
RADNE=topcld.trscat8

ISCCP=ISCCP.tautables
GHG=GHG.CMIP6.1-2014.txt  !  GreenHouse Gases for CMIP6 runs up to 2014
CO2profile=CO2profile.Jul16-2017.txt ! scaling of CO2 in stratosphere
dH2O=dH2O_by_CH4_monthly

! Begin NINT E2.1 input files

BCdalbsn=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCdalbsn
DUSTaer=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/DUST
TAero_SUL=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SUL
TAero_SSA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/SSA
TAero_NIT=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/NIT
TAero_OCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/OCA
TAero_BCA=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCA
TAero_BCB=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/BCB

O3file=cmip6_nint_inputs_E14TomaOCNf10_4av_decadal/O3
Ox_ref=o3_2010_shindell_144x90x49_April1850.nc

! End NINT E2.1 input files
!-----------------------------------------------
!  resolution-independent chemistry input files:
!-----------------------------------------------
JPLRX=chem_files/JPL2011_CMIP6_E2.1_fastterp_newRX2
JPLPH=chem_files/ds4_photlist_T25
RATJ=chem_files/ratj.giss_25
SPECFJ=chem_files/jv_spec_AV_X68d.dat
ATMFJ=chem_files/jv_atms.dat
LNOxCDF=lightning/light_dist.ott2010.dat
! these to scale various initial conditions:
trICratOth=giss2/IC/trICratOth.nc
trICratN=giss2/IC/trICratN.nc
trICratN2O=giss2/IC/trICratN2O.nc
trICratCOt=giss2/IC/trICratCOt.nc
trICratCOs=giss2/IC/trICratCOs.nc
trICratCFC=giss2/IC/trICratCFC.nc
trICch4=giss2/IC/trICch4.nc
!-----------------------------------------------
!-----------------------------------------------
!  3D chemistry input files:
!-----------------------------------------------
N2O_IC=giss2/IC/N2O_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
CFC_IC=giss2/IC/CFC_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
CH4_IC=giss2/IC/CH4_IC_M23_4x5_6.17_conc_2x2.5_conc.nc
Ox_IC=giss2/IC/Ox_init_cond_M23_4x5_conc_2x2.5_conc.nc
CO_IC=giss2/IC/CO_init_cond_M23_conc_2x2.5_conc.nc
ALB_IC=giss2/ALBIJ1_IC_144x90_from_ANN1854_of_M40cadiOK0.nc
!! OFFLINE_AERO=giss2/uncoupledAero_E14TomaOCNf10_4av_decadal_F40 ! only for coupled_chem.eq.0
! files for dust tracers
ERS=ERS1_1993_MONTHLY.144x90.threshold.nc ! ERS data
DSRC=Ginoux_source_v2009_VegMask_144x90.nc   ! preferred dust sources
! alternative preferred dust source files:
!  Ginoux2001_source_VegMask_144x90
!  Ginoux_source_v2009_NoVegMask_144x90
!  GriniZender_DustSources_144x90
!  Tegen_DustSources_144x90
LKTAB=log_dust_emission_60ms-1 ! look up table for emission calculations
LKTAB1=table_wspdf             ! look up table for wind speed probabilities
! mineral fractions of dust aerosols at emission (also needed for default
! one-type dust if size size distribution of mineral fraction version is
! used for dust emission; in this case set imDust=4)
MINFR=mineralfractionsatemission_AMFmethod_NASAGISS_201412_v2_144x90.nc
! for AeroCom year 2000 simulations (set imDust=1) or if AeroCom size
! distribution is used for dust emission (set imDust=3) or if AeroCom source
! distribution by location is used as additional mask for model calculated
! emission (set imDust=5)
!dust_bin1=DUST_bin1_2000_2x2.5.nc
!dust_bin2=DUST_bin2_2000_2x2.5.nc
!dust_bin3=DUST_bin3_2000_2x2.5.nc
!dust_bin4=DUST_bin4_2000_2x2.5.nc
!---------- start chemistry emissions files --------------
CO_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
CO_09=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
NOx_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
NOx_09=emisnc/F/NAT/NOx_Soil_GEIA_2x2.5_HALF_h.nc ! half because we have ag source
NOx_10=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
! Note that the Isoprene emis file is ignored when BIOGENIC_EMISSIONS
! directive is on. But I am commenting it anyway.
! if BIOGENIC_EMISSIONS or PS_BVOC are defined, and only one
! Terpenes file is available, the Isoprene file is needed too
! Isoprene_01=ORCHIDEE_Isoprene_1990_2x2.5_h
Terpenes_01=emisnc/F/OTHER/ORCHIDEE_Terpenes_1990_2x2.5_h.nc
Terpenes_02=emisnc/F/OTHER/ORCHIDEE_ORVOC_1990_2x2.5_h.nc
! === If you add any of your own Alkenes or Paraffin emissions ===
! === please remember they are to be emitted with molecular wt ===
! === of 1.0. In other words in Kmole m-2 s-1 units            ===
Alkenes_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
Alkenes_09=emis/Alkenes_vegetation_GEIA_2x2.5_sname.nc
Alkenes_10=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
Paraffin_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
Paraffin_09=emis/Paraffin_vegetation_GEIA_2x2.5_sname.nc
Paraffin_10=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
codirect_01=emisnc/F/OTHER/HTAP_codirect_emissions_2x2.5_h.nc
!------------ end chemistry emissions files --------------
!-------Aerosol inputs--------------
!----oxidants needed if not coupled to chem ----
!---these are place-holders, right now must run coupled ----
!OFFLINE_CHEM=giss2/uncoupledChem_E14TomaOCNf10_4av_decadal_F40
!-----------------------------------------------
!       AMP Radiation Input
AMP_MIE_TABLES=AMP_MIE_Q_G_A_S.nc
AMP_CORESHELL_TABLES=AMP_CORESHELL_TABLES.nc
!---------- start aerosol emissions files --------------
PSREF=ANN1960.E70F40pi.prsurf.nc  ! time avg. surf. pres. on model grid
SO2_VOLCANO=SO2_volc_2000_1x1_pres.nc
DMS_SEA=DMS_Kettle_Andreae_2x2.5.nc
NH3_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
NH3_09=emisnc/F/OTHER/NH3hCON_OCEANflux_Jan10_2x2.5_h.nc
NH3_10=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_BC1_BC_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_BC1_BC_09=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM2_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM1_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCM0_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP1_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP2_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP3_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP4_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP5_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_AIRC=emis/M_OCC/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_01=emis/M_OCC/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_02=emis/M_OCC/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_03=emis/M_OCC/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_04=emis/M_OCC/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_05=emis/M_OCC/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_06=emis/M_OCC/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_07=emis/M_OCC/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_08=emis/M_OCC/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
M_OCC_OCP6_09=emis/M_OCC/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
SO2_AIRC=emis/CMIP6_AIR_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_01=emis/CMIP6_IND_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_02=emis/CMIP6_TRA_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_03=emis/CMIP6_WST_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_04=emis/CMIP6_RCO_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_05=emis/CMIP6_SHP_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_06=emis/CMIP6_SLV_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_07=emis/CMIP6_ENE_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_08=emis/CMIP6_AGR_2017-05-18_2.0x2.5_CLIMATOLOGY
SO2_09=emis/CMIP6_BBURN_v1.2_2.0x2.5_CLIMATOLOGY
!------------ end aerosol emissions files --------------

MSU_wts=MSU_SSU_RSS_weights.txt      ! MSU-diag
REG=REG2X2.5                      ! special regions-diag

Label and Namelist:  (next 2 lines)
ygMV_clim2015_BBctrl (climatological prescribed ocean atmospheric tracer model with MATRIX and Shindell chemistry)

&&PARAMETERS
! parameters set for choice of ocean model:
KOCEAN=0        ! ocean is prescribed
!! KOCEAN=1        ! ocean is computed
Kvflxo=0        ! usually set to 1 only during a prescr.ocn run by editing "I"
!  Kvflxo=1     ! saves VFLXO files to prepare for q-flux runs (mkOTSPEC)

variable_lk=1   ! variable lakes

! drag params if grav.wave drag is not used and top is at .01mb
X_SDRAG=.002,.0002  ! used above P(P)_sdrag mb (and in top layer)
C_SDRAG=.0002       ! constant SDRAG above PTOP=150mb
P_sdrag=1.          ! linear SDRAG only above 1mb (except near poles)
PP_sdrag=1.         ! linear SDRAG above PP_sdrag mb near poles
P_CSDRAG=1.         ! increase CSDRAG above P_CSDRAG to approach lin. drag
Wc_JDRAG=30.        ! crit.wind speed for J-drag (Judith/Jim)
ANG_sdrag=1     ! if 1: SDRAG conserves ang.momentum by adding loss below PTOP
! vsdragl is a tuning coefficient for SDRAG starting at LS1
! layer:   24    25    26    27   28    29    30    31   32   33     34   35   36  37  38   39 40
vsdragl=0.000,0.000,0.000,0.000,0.00,0.000,0.000,0.000,0.00,0.00,  0.00,0.00,0.00,0.3,0.6,0.83,1.

! Gravity wave parameters
PBREAK = 200.  ! The level for GW breaking above.
DEFTHRESH=0.000055  ! threshold (1/s) for triggering deformation waves
PCONPEN=400.   ! penetrating convection defn for GWDRAG
CMC = 0.0000002 ! parameter for GW Moist Convective drag
CSHEAR=10.     ! Shear drag coefficient
CMTN=0.1       ! default is 0.5
CDEF=1.6       ! tuning factor for deformation -> momentum flux
XCDNST=400.,10000.   ! strat. gw drag parameters
QGWMTN=1 ! mountain waves ON
QGWDEF=1 ! deformation waves ON
QGWSHR=0 ! shear drag OFF
QGWCNV=0 ! convective drag OFF


! cond_scheme=2   ! newer conductance scheme (N. Kiang) ! not used with Ent

! The following two lines are only used when aerosol/radiation interactions are off
FS8OPX=1.,1.,1.,1.,1.5,1.5,1.,1.
FT8OPX=1.,1.,1.,1.,1.,1.,1.3,1.

! Increasing U00a decreases the high cloud cover; increasing U00b decreases net rad at TOA
U00a=0.625  ! above 850mb w/o MC region;  tune this first to get 30-35% high clouds
U00b=1.00   ! below 850mb and MC regions; tune this last  to get rad.balance
WMUI_multiplier = 2.
use_vmp=1
radius_multiplier=1.1

PTLISO=0.        ! pressure(mb) above which radiation assumes isothermal layers
H2ObyCH4=0.      ! if =1. activates stratospheric H2O generated by CH4 without interactive chemistry
KSOLAR=2         ! 2: use long annual mean file ; 1: use short monthly file

! parameters that control the atmospheric/boundary conditions
! if set to 0, the current (day/) year is used: transient run
master_yr=2015
!crops_yr=2015  ! if -1, crops in VEG-file is used
!s0_yr=2015
!s0_day=182
!ghg_yr=2015
!ghg_day=182
!irrig_yr=1850
volc_yr=-1
!volc_day=182
!aero_yr=2015
od_cdncx=0.        ! don't include 1st indirect effect
cc_cdncx=0.        ! don't include 2nd indirect effect (used 0.0036)
!albsn_yr=2015
dalbsnX=1.
!o3_yr=-2015
!aer_int_yr=2015    !select desired aerosol emissions year or 0 to use JYEAR
! atmCO2=368.6          !uatm for year 2000 - enable for CO2 tracer runs

!variable_orb_par=0
!orb_par_year_bp=100  !  BP i.e. 1950-orb_par_year_bp AD = 1850 AD
MADVOL=2
!!!!!!!!!!!!!!!!!!!!!!!
! Please note that making o3_yr non-zero tells the model
! to override the transient chemistry tracer emissions'
! use of model year and use o3_yr instead!
!!!!!!!!!!!!!!!!!!!!!!!
!--------- general aerosol parameters-------------
rad_forc_lev=0     ! 0: for TOA, 1: for tropopause for rad forcing diags.
rad_interact_aer=1 ! 1: couples aerosols to radiation, 0: use climatology
prather_limits=1   ! 1: to avoid some negative tracers in sub-gridscale
diag_rad=1         ! 1: save ext/sct/asf for 6 bands; only ext for band6 otherwise
diag_aod_3d=3      ! 0: off; 1: save 3d all-sky tau/aaod; 2: save clear-sky; 3: 1+2 combo
                   ! negative values work the same but save total aod only, not speciated
diag_fc=1          ! 2=one radiation call per tracer (slow) || 1=one radiation call || 0=no radiation calls
diag_wetdep=1      ! 1: additional wet deposition diagnostics
VBSemifactFF=0.03d0,0.06d0,0.09d0,0.14d0,0.18d0,0.30d0,0.40d0,0.50d0,0.80d0 
VBSemifactBB=0.03d0,0.06d0,0.09d0,0.14d0,0.18d0,0.30d0,0.40d0,0.50d0,0.80d0
!--- number of biomass burning sources (per tracer) Remember to list those sources last!
NH3_nBBsources=1
SO2_nBBsources=1
M_BC1_BC_nBBsources=1
M_OCC_OCM2_nBBsources=1
M_OCC_OCM1_nBBsources=1
M_OCC_OCM0_nBBsources=1
M_OCC_OCP1_nBBsources=1
M_OCC_OCP2_nBBsources=1
M_OCC_OCP3_nBBsources=1
M_OCC_OCP4_nBBsources=1
M_OCC_OCP5_nBBsources=1
M_OCC_OCP6_nBBsources=1
!------------------  AMP parameters
AMP_RAD_KEY = 1    ! 1=Volume Mixing || 2=Core-Shell || 3=Maxwell Garnett
!--------- dust aerosol parameters----------------
imDust=0                     ! 0: PDF emission scheme, 1: prescr. (AeroCom)
                             ! 3: PDF emission with AeroCom size distr.
                             ! 4: PDF emission with OMA-mineral size distr.
                             ! 5: as 4, but with AeroCom source mask
adiurn_dust=0                ! 1: daily dust diagnostics for selected grid boxes
!to_conc_soildust=1   ! three-dimensional dust output as concentration [kg/m^3]

!! MATRIX with model winds
! MATRIX w/ VMP clouds:
! for imDust=0:
scaleDustEmission=1.0
fracClayPDFscheme=0.0142993763417 ! clay emission parameter from calibration
fracSiltPDFscheme=0.0421175030543 ! silt emission parameter from calibration
   !fracClayPDFscheme = 0.0136705318754[1] * 1.046[2]
   !fracSiltPDFscheme = 0.0402652992871[1] * 1.046[2]
   ![1] emission parameters used for calibration run E20170907TmatrixF40climIM0_000
   ![2] total emitted mass factor derived from calibration run
! for imDust=3
!scaleDustEmission=0.144776901493 ! scales total dust emission
  !scaleDustEmission = 0.15517352786[1] * 0.933[2]
  ![1]: emission parameter used for calibration run E20170907TmatrixF40climIM3_000
  ![2]: total emitted mass factor derived from calibration run
! for imDust=4
!scaleDustEmission=0.145680277394 ! scales total dust emission         
  !scaleDustEmission = 0.148049062392[1] * 0.984[2]
  ![1]: emission parameter used for calibration run E20170907TmatrixF40climIM4_000
  ![2]: total emitted mass factor derived from calibration run
!for imDust=3-5 (also set as default in the model)
  !fracClayPDFscheme=1.0
  !fracSiltPDFscheme=1.0
!-------------------------------------------------
!! MATRIX nudged toward NCEP winds:
! for imDust=0:
!fracClayPDFscheme=0.00678986343007 ! clay emission parameter from calibration
!fracSiltPDFscheme=0.0205586498967  ! silt emission parameter from calibration
   !fracClayPDFscheme = 0.00683772752273[1] * 0.993[2]
   !fracSiltPDFscheme = 0.0207035749211 [1] * 0.993[2]
   ![1] emission parameters used for calibration run E2p1_20190525TmatNcepF40IM0_003
   ![2] total emitted mass factor derived from calibration
! for imDust=3
!scaleDustEmission=0.067761297736 ! scales total dust emission
  !scaleDustEmission = 0.0681018067698[1] * 0.995[2]
  ![1]: emission parameter used for calibration run E2p1_20190525TmatNcepF40IM3_003
  ![2]: total emitted mass factor derived from calibration run
! for imDust=4
!scaleDustEmission=0.0671151632357 ! scales total dust emission
  !scaleDustEmission = 0.0709462613485[1] * 0.946[2]
  ![1]: emission parameter used for calibration run E2p1_20190525TmatNcepF40IM4_003
!for imDust=3-5 (also set as default in the model)
  !fracClayPDFscheme=1.0
  !fracSiltPDFscheme=1.0
!------------------------------------------------------
!! all versions:
! for imDust=1 (prescribed AeroCom emissions):
!scaleDustEmission=1.0
!-------------------------------------------------
!-----------------------------------------------
!  Start tracer code parameters:
!-----------------------------------------------
! Lightning parameterization
lightning_param=1          ! 1 = Cloud Top Height; 2 = Upward Convective Mass Flux; 3 = Convective Precipitation
! Lightning tuning parameters
tune_lt_land=2.5d0
tune_lt_sea= 5.8d0
FLASH_PERTURB=1.0d0        ! 1.1 = 10% increase in lightning flash rate globally

! -----------------------------------
! chemistry debug output
!print_reaction_lists=0
!print_reaction_rates=0
!print_chemical_changes=0
!coords_for_print_chemical=17,35,1 ! ilon,jlat,klev
! -----------------------------------
! Pressure above which Ox, NOx, BrOx, and ClOx will be
! overwritten with climatology based on NINT ozone input.
PltOx=0.0
! NH and SH polar stratospheric cloud formation temperature offsets:
Tpsc_offset_N=0.d0
Tpsc_offset_S=0.d0
! -----------------------------------
! 40-layer model CMIP6 chemistry tuning parameters: 
windowN2Ocorr=1.0
windowO2corr=1.15
reg1Power_SpherO2andN2Ocorr=2.5
reg1TopPres_SpherO2andN2Ocorr=50.
reg2Power_SpherO2andN2Ocorr=2.5
reg2TopPres_SpherO2andN2Ocorr=10.
reg3Power_SpherO2andN2Ocorr=0.5
reg3TopPres_SpherO2andN2Ocorr=1.
reg4Power_SpherO2andN2Ocorr=0.5
! -----------------------------------
COUPLED_CHEM=1     ! to couple chemistry and aerosols
use_sol_Ox_cycle=0 ! (=1) apply ozone changes in radiation, based on solar cycle
clim_interact_chem=1 ! 1=use calculated Ox/CH4 in radiation, 0=use climatology
                   ! If = 0, consider turning on AUXILIARY_OX_RADF CPP directive.
                   ! Note: 0 also turns off chemistry(H2O)-->Q(humidity) feedback
                   ! if you want humidity feedback on but radiation feedback off
                   ! you could do: clim_interact_chem=1, Lmax_rad_{O3,CH4}=0...
! Lmax_rad_O3=0    ! Ox levels used in rad code default is LM
use_rad_n2o=1      ! use the radiation code's N2O
use_rad_cfc=1      ! use rad code cfc11+cfc12, adjusted
rad_FL=1           ! use rad code insolation getting fastj2 photon flux
which_trop=0       ! choose tropopause for chemistry purposes:
                   ! 0=LTROPO(I,J), 1=LS1-1
!--- number of biomass burning sources (per tracer) Remember to list those sources last!
Alkenes_nBBsources=1
CO_nBBsources=1
NOx_nBBsources=1
Paraffin_nBBsources=1
! -----------------------------------
! Lightning NOx yield per flash
FLASH_YIELD_MIDLAT=290.0d0 ! NOx yield per flash (moles N/flash), applied poleward of 23deg N/S
FLASH_YIELD_TROPIC=290.0d0 ! NOx yield per flash (moles N/flash), applied equatorward of 23deg N/S
! -----------------------------------
! Tune Isoprene emissions (biogenic_emissions version, not MEGAN):
base_isopreneX=0.6

! --- CH4 tracer specific settings: ---
!    important:
use_rad_ch4=1      ! use rad code CH4, shut off sfc sources
! OFF FOR NOW: CH4_nBBsources=1
!    rarer settings:
! Lmax_rad_CH4=0   ! CH4 levels used in rad code default is LM
fix_CH4_chemistry=0    ! for setting fixed methane value for chemistry:
scale_ch4_IC_file=1.d0 ! multiplicative factor on CH4 IC file (fix_CH4_chemistry=-1)
! --- end of CH4 tracers specific settings ---

DTsrc=1800.      ! cannot be changed after a run has been started
DT=225.
! parameters that control the Shapiro filter
DT_XUfilter=225. ! Shapiro filter on U in E-W direction; usually same as DT
DT_XVfilter=225. ! Shapiro filter on V in E-W direction; usually same as DT
DT_YVfilter=0.   ! Shapiro filter on V in N-S direction
DT_YUfilter=0.   ! Shapiro filter on U in N-S direction

NIsurf=2         ! surface interaction computed NIsurf times per source time step
NRAD=5           ! radiation computed NRAD times per source time step
! parameters that affect at most diagn. output:  standard if DTsrc=1800. (sec)
aer_rad_forc=0   ! if set =1, radiation is called numerous times - slow !!
cloud_rad_forc=1 ! calls radiation twice; use =0 to save cpu time, 2= calculates crf_toa2
SUBDD=' '        ! no sub-daily frequency diags
NSUBDD=0         ! saving sub-daily diags every NSUBDD-th physics time step (1/2 hr)
KCOPY=1          ! 0: no output; 1: save .acc; 2: unused; 3: include ocean data
KRSF=12          ! 0: no output; X: save rsf at the beginning of every X month
isccp_diags=1    ! use =0 to save cpu time, but you lose some key diagnostics
nda5d=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
nda5s=13         ! use =1 to get more accurate energy cons. diag (increases CPU time)
ndaa=13
nda5k=13
nda4=48          ! to get daily energy history use nda4=24*3600/DTsrc

Nssw=2           ! until diurnal diags are fixed, Nssw has to be even
Ndisk=960        ! write fort.1.nc or fort.2.nc every NDISK source time step
&&END_PARAMETERS

&INPUTZ
 YEARI=2010,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=2021,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=2010,MONTHE=12,DATEE=1,HOURE=1,
/
### Information below describes your run. Do not delete! ###
Wed Mar 25 08:44:05 EDT 2020
ifort version 19.1.0.166
Git reference: bfffdd83f27acced35476567f16519058b53af78
