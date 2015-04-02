#include "rundeck_opts.h"
      MODULE tracers_dust
!@sum  tracers_dust dust/mineral tracer parameter and variable declarations
!@auth Reha Cakmur, Jan Perlwitz, Ina Tegen
!@ver 2.0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      USE constant,ONLY : By6
      USE resolution,ONLY : Im,Jm,Lm
      use tracer_com, only: ntm_dust
      use OldTracer_mod, only: trName, MAX_LEN_NAME

      IMPLICIT NONE

!@param By8 0.25d0/2d0
      REAL*8,PARAMETER :: By8=0.25D0/2D0
!@param By4 1D0/4D0
      REAL*8,PARAMETER :: By4=1D0/4D0

!@param n_soilDust index of first soil dust aerosol tracer
      integer :: n_soilDust
!@param dust_names names of soil dust aerosol tracers
      character(len=MAX_LEN_NAME),allocatable,dimension(:) ::
     &   dust_names(:)

c**** rundeck parameter to switch between different emission schemes
c****
!@dbparam imDust: 0: scheme using PDF of wind speed (default)
!@+               1: prescribed AEROCOM emissions
!@+               2: legacy emission scheme using third power of wind speeds
!@+                  (only works with 72x46 horizontal resolution)
      INTEGER :: imDust=0

c**** legacy emission code (Tegen, I. and R. Miller, JGR (1998))
c**** declarations for emission scheme using third power of wind speed
c****
!@param CWiCub uplift factor [kg*s**2/m**5] for all size classes
      REAL*8,PARAMETER :: CWiCub=52.D-9
!@param FClWiCub fraction [1] of uplifted clay for scheme using cubes of
!@+     wind speed
!@param FSiWiCub fractions [1] of uplifted silt for scheme using cubes of
!@+     wind speed
      REAL*8 :: FClWiCub=By6,FSiWiCub=By8
!@var hbaij accumulated precipitation - evaporation balance  [kg/m^2]
!@var ricntd no. of hours with negative precipitation - evaporation balance [1]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: hbaij,ricntd
!@var dryhr  number of hours with evaporation-precipitation greater Zero
!@+          to allow dust emission
!@var frclay fraction of clay
!@var frsilt fraction of silt
!@var vtrsh  threshold wind speed above which dust emis. is allowed [m/s]
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: dryhr,frclay,frsilt,vtrsh

c**** for legacy wet deposition
c**** declaration for simple wet deposition scheme
!@var prelay distributed array with some prec info needed for simple wet
!@+          deposition scheme
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: prelay

c**** current default emission scheme (Cakmur, R. et al. (2004))
c**** declarations for emission scheme using probability density function of
c**** wind speed

!@dbparam prefDustSources  rundeck parameter to choose preferred source file
!+                         and according parameter settings for the fractions
!+                         of uplifted clay and silt
!+                         0: Ginoux 2001 with overlaid annual vegetation mask
!+                         1: Ginoux 2009 with overlaid annual vegetation mask
!+                         2: Ginoux 2009 w/o overlaid vegetation mask
!+                         3: Grini/Zender preferred sources
!+                         4: Tegen preferred sources
!+                        >4: Free choice of file
      integer :: prefDustSources=1

!@param numDustSourceOpt number of options for preferred dust sources
      integer,parameter :: numDustSourceOpt = 5
!@param numResolutionOpt number of resolution options for preferred dust sources
      integer,parameter :: numResolutionOpt = 5
!@param dustSourceFile default files with preferred dust sources for output
!@+                    info
      character(len=71),parameter
     &    ,dimension(0:numDustSourceOpt-1,numResolutionOpt) ::
     &    dustSourceFile = reshape((/
     &    'Ginoux2001_source_VegMask_72x46    (optimized for old model)'
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'Ginoux2001_source_VegMask_144x90                            '
     &   ,'Ginoux_source_v2009_VegMask_144x90                          '
     &   ,'Ginoux_source_v2009_NoVegMask_144x90                        '
     &   ,'GriniZender_source_144x90                                   '
     &   ,'Tegen_source_144x90                                         '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'No file available yet                                       '
     &   ,'Ginoux2001_source_VegMask_144x90_C90                        '
     &   ,'Ginoux_source_v2009_VegMask_0.25x0.25_to_C90                '
     &   ,'Ginoux_source_v2009_NoVegMask_144x90_C90                    '
     &   ,'GriniZender_DustSources_C90_from_144x90                     '
     &   ,'Tegen_DustSources_C90_from_0.5x0.5                          '
     &   /),(/numDustSourceOpt,numResolutionOpt/))

!@param CWiPdf uplift factor [kg*s**2/m**5] for all size classes of soil dust
      REAL*8,PARAMETER :: CWiPdf=12.068996D-9
!@dparam FracClayPDFscheme fraction [1] of uplifted clay
!@dparam FracSiltPDFscheme fractions [1] of uplifted silt
      real(kind=8) :: fracClayPDFscheme = 1.D0, fracSiltPDFscheme = 1.D0
!@var ers_data field of ERS data
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: ers_data
!@var dustSourceFunction distribution of preferred dust sources
      real(kind=8),allocatable,dimension(:,:) :: dustSourceFunction

      ! note that Lim,Ljm happen to match kim,kjm in wspdf_mod
      INTEGER,PARAMETER :: Lim=234,Ljm=234,Lkm=22
!@var x1,x2,x3 indices of lock up table for emission
      REAL*8 :: x1(Lim),x2(Ljm),x3(Lkm)
!@var table array of lock up table for emission local to each grid box
      real*8 table(Lim,Ljm,Lkm)

      !REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: table
      !ALLOCATE(table(Lim,Ljm,Lkm),STAT=ier)

!@var wsubtke_com distributed array of subscale turbulent term
!@var wsubwd_com distributed array of subscale dry convective term
!@var wsubwm_com distributed array of subscale moist convective term
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: wsubtke_com,wsubwd_com,
     &     wsubwm_com

c****
c**** declarations for prescribed daily dust emissions
c****
!@param nAerocomDust Number of AEROCOM dust classes
      integer,parameter :: nAerocomDust = 4
!@var d_dust prescribed daily dust emissions [kg] (e.g. AEROCOM)
      REAL*8,ALLOCATABLE,DIMENSION(:,:,:,:) :: d_dust

c**** additional declarations for dust tracers with mineralogical composition
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!@param Mtrac number of different fields with tracer fractions in grid box
!@+           5 clay; 5 silt
      INTEGER,PARAMETER :: Mtrac=10
!@param DensityQuartz  particle density of Quartz
!@+            (measured; http://www.mindat.org/min-3337.html)
      real(kind=8), parameter :: DensityQuartz = 2.655d3 ! 
!@param DensityHematite  particle density of average of Hematite and Goethite
!@+            (Hematite measured: 5.26d3; http://www.mindat.org/min-1856.html)
!@+            (Goethite measured; 4.28d3; http://www.mindat.org/min-1719.html)
      real(kind=8), parameter :: DensityHematite = 4.77d3
!@dbparam brittleFactor  one-dimension array of size Mtrac with factors to
!@+                      weight fractionation of minerals during emission
      real(kind=8), dimension(Mtrac) :: brittleFactor = 1.d0
#ifdef TRACERS_QUARZHEM
      REAL*8 :: FreeFe=0.5D0 ! not used anywhere in code, currently
!@dbparam frHemaInQuarAggr fraction of hematite in quartz/hematite aggregate
      real(kind=8) :: frHemaInQuarAggr = 0.1d0 ! arbitrary assumption
!@dbparam purebyTotalHematite fraction of pure hematite minerals to all
!@+                           hematite in minerals (pure + aggregate)
      real(kind=8) :: pureByTotalHematite = 0.1d0 ! arbitrary assumption
!@dbparam calcMineralAggrProb 1: calculate mineral aggregation probabilities
!@+                           from mineral fractions in soil; 0: use
!@+                           prescribed pureByTotalHematite
      integer :: calcMineralAggrProb = 0
#endif
!@var minfr distribution of mtrac mineral fractions in soils
       REAL*8,ALLOCATABLE,DIMENSION(:,:,:) :: minfr
!@var mineralFractions distribution of mineral fractions in soils for each
!@+   mineralogical soil dust tracer
      real(kind=8), allocatable, dimension(:,:,:) :: mineralFractions
#endif

c**** Parameters for dust/mineral tracer specific diagnostics
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@param nDustEmij index of dust emission in ijts_isrc
      INTEGER,PARAMETER :: nDustEmij=1
!@param nDustEm2ij index of dust emission according to cubic scheme
!@+                in ijts_isrc
      INTEGER,PARAMETER :: nDustEm2ij=2
!@param nDustTurbij index of dust dry turbulent deposition in ijts_isrc
      INTEGER,PARAMETER :: nDustTurbij=3  ! not used?
!@param nDustEv1ij index of number of dust events below threshold wind
!@+                in ijts_spec
!@param nDustEv2ij index of number of dust events above threshold wind
!@+                in ijts_spec
!@param nDustWthij index of threshold velocity in ijts_spec
      INTEGER,PARAMETER :: nDustEv1ij=1,nDustEv2ij=2,nDustWthij=3
#endif
!@dbparam to_conc_soildust: For printout of 3D soil dust aerosol concentration
!@+   in kg/m3
!@+   to_conc_soildust = 0: printout is as defined by to_volume_MixRat
!@+   to_conc_soildust = 1: printout is in kg/m3
      integer :: to_conc_soildust = 0

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
!@param nDustEmjl index of dust emission in jls_source
      INTEGER,PARAMETER :: nDustEmjl=1
!@param nDustEm2jl index of dust emission according to cubic scheme
!@+                in jls_source
      INTEGER,PARAMETER :: nDustEm2jl=2
!@param nDustTurbjl index of dust dry turbulent deposition in jls_source
      INTEGER,PARAMETER :: nDustTurbjl=3
!@param nDustEv1jl index of number of dust events below threshold wind
!@+                in jls_spec
!@param nDustEv2jl index of number of dust events above threshold wind
!@+                in jls_spec
!@param nDustWthjl index of threshold velocity in ijts_spec
      INTEGER,PARAMETER :: nDustEv1jl=1,nDustEv2jl=2,nDustWthjl=3
#endif

c**** Variables for specific subdaily soil dust aerosol diagnostics
!@var dustDiagSubdd structured type for specific subdaily dust
!+                  aerosol diagnostics
!@var  %dustEmission   soil dust aerosol emission [kg/m^2/s]
!@var  %dustEmission2  soil dust emission calculated using cubed formula
!@+                    [kg/m^2/s] (only for diagnostic purposes)
!@var  %dustDepoTurb   turbulent soil dust aerosol deposition  [kg/m^2/s]
!@var  %dustDepoGrav   gravitational settling of soil dust aerosols [kg/m^2/s]
!@var  %dustMassInPrec wet deposition of soil dust aerosols [k]
!@var  %dustSurfMixR   surface mixing mixing ratio of soil dust aerosols [kg/kg]
!@var  %dustSurfConc   surface concentration of soil dust aerosols [kg/m^3]
!@var  %dustMass       mass of soil dust aerosol [kg]
!@var  %dustConc       dust concentration [kg/m] (later divided by grid box area)
      type dustDiagSubdd
      real(kind=8),allocatable,dimension(:,:,:) :: dustEmission
     &     ,dustEmission2,dustDepoTurb,dustDepoGrav,dustMassInPrec
     &     ,dustSurfMixR,dustSurfConc
      real(kind=8),allocatable,dimension(:,:,:,:) :: dustMass,dustConc
      end type dustDiagSubdd
!@var dustDiagSubdd_acc structured variable to accumulate data for
!+                      subdaily dust aerosol diagnostics
      type(dustDiagSubdd) :: dustDiagSubdd_acc

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP || TRACERS_TOMAS*/

      END MODULE tracers_dust

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

      SUBROUTINE alloc_dust(grid)
!@sum  alloc_dust allocates dust/mineral tracer arrays
!@auth Jan Perlwitz

      USE domain_decomp_atm, ONLY : dist_grid
      USE resolution,ONLY : Lm
      USE tracer_com,ONLY : Ntm_dust
      use tracers_dust
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      use AbstractCalendar_mod
      use MODEL_COM, only: calendar
      IMPLICIT NONE

      TYPE(DIST_GRID),INTENT(IN) :: grid

      INTEGER :: i_0h,i_1h,j_1h,j_0h
      INTEGER :: ier
      LOGICAL,SAVE :: qfirst=.TRUE.
      integer :: maxDaysInYear

      IF (.NOT. qfirst) RETURN
      qfirst=.FALSE.

      i_0h=grid%i_strt_halo
      i_1h=grid%i_stop_halo
      j_0h=grid%j_strt_halo
      j_1h=grid%j_stop_halo

      allocate(dust_names(ntm_dust))

      maxDaysInYear = calendar%getMaxDaysInYear()

      ALLOCATE(hbaij(i_0h:i_1h,j_0h:j_1h),ricntd(i_0h:i_1h,j_0h:j_1h),
     &     dryhr(i_0h:i_1h,j_0h:j_1h),frclay(i_0h:i_1h,j_0h:j_1h),
     &     frsilt(i_0h:i_1h,j_0h:j_1h),vtrsh(i_0h:i_1h,j_0h:j_1h),
     &     dustSourceFunction(i_0h:i_1h,j_0h:j_1h),
     &     ers_data(i_0h:i_1h,j_0h:j_1h,INT_MONTHS_PER_YEAR),
     &     wsubtke_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwd_com(i_0h:i_1h,j_0h:j_1h),
     &     wsubwm_com(i_0h:i_1h,j_0h:j_1h),
     &     prelay(i_0h:i_1h,j_0h:j_1h,LM),
     &     d_dust(i_0h:i_1h,j_0h:j_1h,nAerocomDust,maxDaysInYear),
#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
     &     minfr(i_0h:i_1h,j_0h:j_1h,Mtrac),
     &     mineralFractions( i_0h:i_1h, j_0h:j_1h, ntm_dust ),
#endif
     &     STAT=ier)

      d_dust(i_0h:i_1h,j_0h:j_1h,:,:)=0.D0

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      minfr( i_0h:i_1h, j_0h:j_1h, : ) = 0.d0
      mineralFractions( i_0h:i_1h, j_0h:j_1h, : ) = 0.d0
#endif

      allocate(dustDiagSubdd_acc%dustEmission(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustEmission2(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustDepoTurb(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustDepoGrav(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustMassInPrec(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustSurfMixR(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustSurfConc(i_0h:i_1h,j_0h:j_1h
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustMass(i_0h:i_1h,j_0h:j_1h,Lm
     &     ,Ntm_dust))
      allocate(dustDiagSubdd_acc%dustConc(i_0h:i_1h,j_0h:j_1h,Lm
     &     ,Ntm_dust))

      dustDiagSubdd_acc%dustEmission = 0.D0
      dustDiagSubdd_acc%dustEmission2 = 0.D0
      dustDiagSubdd_acc%dustDepoTurb = 0.D0
      dustDiagSubdd_acc%dustDepoGrav = 0.D0
      dustDiagSubdd_acc%dustMassInPrec = 0.D0
      dustDiagSubdd_acc%dustSurfMixR = 0.D0
      dustDiagSubdd_acc%dustSurfConc = 0.D0
      dustDiagSubdd_acc%dustMass = 0.D0
      dustDiagSubdd_acc%dustConc = 0.D0

      RETURN
      END SUBROUTINE alloc_dust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP || TRACERS_TOMAS*/
