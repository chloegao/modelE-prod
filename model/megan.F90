#include "rundeck_opts.h"

! Please see ../doc/megan_suggested_todo.txt for further notes and suggestions
! for improvement that have been extracted from this program's comments.

!TODO: save at *least* some running average stuff to the rsf files 
!      (see flammability code for guide)
!
!TODO: When writing the flammability code, exporting LAI from Ent, I noted:
!      "I guess that is the LAI from the *last surface timestep only*?"
!      Igor said that was OK, as LAI is only computed once per day. However,
!      now we are using ent_get_exports all over the place, so
!      we need to revisit this worry about time steps.

!TODO: find out from someone whether my "save" or lack thereof will allow
!      objects like SAT to be persistent between calls. (e.g. not just with respect
!      to restarts but to program scope.)

!TODO: Tune the CCE parameter (see noted below)
!
! MORE TODOS are below (Throughout)


module megan_objects_mod
use constant, only: undef

implicit none
private
public :: runningAverage
public :: biogenicSpecies
integer, parameter, public :: nMeganPFT=16 ! number of MEGAN plant functional types

type runningAverage
  integer :: stepsPerDay=0 ! expected accumulation steps each day (e.g. 48 for DTsrc=1800.)
  integer :: daysPerPeriod=0 ! expected accumulation days in averaging period
  logical :: laggedValue=.false. ! if true, instead of the runningAverage array holding a
                      ! running average, it holds the value from daysPerPeriod days ago.
  ! (:,:) arrays here are to be I,J dimensions.
  ! (:,:,:) are either I,J,stepsPerDay(stepSave) or IM,JM,daysPerPeriod (dayAvg)
  logical, allocatable, dimension(:,:):: first ! Whether in first averaging period
  real*8, allocatable, dimension(:,:):: step ! Saves number of local accum. in current day
  real*8, allocatable, dimension(:,:,:):: stepSave ! saves step() through day
  real*8, allocatable, dimension(:,:):: day ! Saves number of day accum. in first period
  real*8, allocatable, dimension(:,:,:):: dayAvg ! daily avg saved each day in period
  real*8, allocatable, dimension(:,:):: runningAverage ! The current running average
  real*8, allocatable, dimension(:,:):: periodRunningSum ! The current running sum
  integer, allocatable, dimension(:,:):: marker ! Current position in runningAverage
end type runningAverage

type biogenicSpecies
  character*8 :: itsname='_unknown' ! name of species, which should match trname
  real*8 :: cceo=undef ! Coefficient for temperature activity factor in gamma_tld routine
  real*8 :: ct1=undef ! A temperature needed for the gamma_tld routine
  real*8 :: tdf_prm=undef ! a temperature-dependent parameter needed for gamma_tli routine
  real*8 :: ldf=undef ! light dependant fraction, used for relative weighting between
                      ! gamma_PPFD*gamma_tld (ldf) and gamma_tli (1-ldf)
  integer :: aindx=-1 ! an index to position in arrays Anew, Agro, Amat, Aold for aging
                      ! gamma routine. Related to relative emission acitivity?
  real*8, allocatable, dimension(:,:):: source ! holds kg m-2 s-1 emissions source for exporting
  real*8, dimension(nMeganPFT) :: ef ! emissions factors by MEGAN PFT
end type biogenicSpecies

end module megan_objects_mod


module megan
!@sum Contains MEGAN routines, including gamma calculations from MEGAN2.1
!@+ and modelE routines to use them for biogenic emissions calculations.
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
! Tan Sakulyanontvittaya, Alex Guenther
! See: http://www.lar.wsu.edu/megan/
! These routines mostly come from the MEGAN2.1 code: EMPROC/gamma_etc.f

! Note, throughout, "G 2012" means:
! Guenther, A. B., X. Jiang, C. L. Heald, et al. 2012The Model of Emissions of
!    Gases and Aerosols from Nature Version 2.1 (MEGAN2.1): An Extended and
!    Updated Framework for Modeling Biogenic Emissions. Geoscientific Model
!    Development 5(6): 1471–1492.
!
! and any references to "hammoz" is referring to the ECHAM-HAMMOZ 
! implementation of MEGAN referenced in the online MEGAN2.1 documentation. 
! Specifically the single file: mo_hammoz_emi_biogenic.f90 obtained from 
! Alexandra Henrot (see her paper:
! http://www.geosci-model-dev-discuss.net/gmd-2016-248/gmd-2016-248.pdf)
!
! Note that the MEGAN gamma routines used to loop over ncols and nrows within.
! We now call from a driver that loops over i,j, passes in variables needed
! from the GCM and returns the gamma for the current i,j conditions. Along the 
! eay, changing to double-precision, using GCM constants when obvious, etc.
!
! Most notes that look like this next section are copied from MEGAN directly
! for reference, though ACTUAL USAGE MAY DIVERGE IN SPECIFICS:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Scientific algorithm
!
!             Emission = [EF][GAMMA][RHO]
!           where [EF]    = emission factor (ug/m2h)
!                 [GAMMA] = emission activity factor (non-dimension)
!                 [RHO]   = production and loss within plant canopies
!                           (non-dimensino)
!                 Assumption: [RHO] = 1 (11/27/06) (See PDT_LOT_CP.EXT)
!
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM]
!           where [GAMMA_CE]  = canopy correction factor
!                 [GAMMA_age] = leaf age correction factor
!                 [GAMMA_SM]  = soil moisture correction factor
!                 Assumption: [GAMMA_SM]  = 1 (11/27/06)
!
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T]
!           where [GAMMA_LAI] = leaf area index factor
!                 [GAMMA_P]   = PPFD emission activity factor
!                 [GAMMA_T]   = temperature response factor
!
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!        Derivation:
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P]
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] }
!             Emission = [EF][GAMMA_ect]{ (1-LDF) + [LDF][GAMMA_P] }

!     For ISOPRENE
!                 Assumption: LDF = 1 for isoprene            (11/27/06)
!
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age][GAMMA_SM]
!
!     For NON-ISOPRENE
!        Final Equation
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age][GAMMA_SM]*
!                        { (1-LDF) + [LDF][GAMMA_P] }
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use megan_objects_mod, only: runningAverage, biogenicSpecies

implicit none

!TODO: are these "save"s needed? Elsewhere? (I think the flammability code has a global save in the module...)
type(runningAverage), save :: SAT, T, LAI, PPFD
type(biogenicSpecies), save :: isoprene

! next two lines are from MEGAN2.1 canopy.f, but other models
! may use 4.6 or 4.55 not separated by Shade and Sun (See
! G 2012 paper.) The idea here is that the "Sun" one is more
! appropriate for direct sunlight hitting a leaf and "Shade"
! one is better for diffuse light hitting shaded canopy leaves.
! Units are [micro mole photons] per [Joule]:
real*8, parameter :: ConvertShadePPFD = 4.6d0
real*8, parameter :: ConvertSunPPFD = 4.0d0


end module megan


subroutine biogenicEmissions_drv(i,j)
!@sum calculate biogenic emissions of atmospheric constituents using
!@+ MEGAN model 2.1 and fill in source array
!@auth Greg Faluvegi (intial modelE implementation)

use megan
use resolution, only: IM
use model_com, only: modelEclock,itime
use fluxes, only: atmsrf
use ghy_com, only: fearth
use ent_com, only: entcells,n_covertypes
use ent_mod, only: ent_get_exports
use rad_com, only: cosz1, CO2ppm, CO2X, FSRDIR, SRVISSURF
use constant, only: radian, undef, tf
use TimeConstants_mod, only: HOURS_PER_DAY, SECONDS_PER_HOUR
use megan_objects_mod, only: runningAverage, biogenicSpecies, nMeganPFT
use OldTracer_mod, only: nBBsources,trname,do_megan,itime_tr0
use tracer_com, only: ntm, ntsurfsrc, sfc_src
use Tracer_mod, only: ntsurfsrcmax

implicit none

! The *_megan variables here are quantities extracted from the GCM for
! input to megan gamma routines.
! The gamma_* variables are unitless values returned from the calls
! to megan gamma routines.
real*8 :: CO2_megan,LAI_megan,cosSZA_megan,PPFD_megan,PPFD_daily_megan
real*8 :: LAI_current_megan, LAI_previous_megan, T_daily_megan, T_megan
real*8 :: SAT_daily_megan, SAT_megan
integer :: JDAY_megan
real*8 :: gamma_CO2, gamma_LAI, gamma_PPFD, gamma_AGE, gamma_SM 
real*8 :: gamma_TLD, gamma_TLI, bulk_EF
real*8 :: par_total, par_direct, par_diffuse ! assumed components of ppfd below
real*8, parameter :: radianToDegree=1.d0/radian
!@param convertUnits to convert from emission factor in microgram m-2 hr-1
!@+ from megan to kg m-2 s-1 for GCM
real*8, parameter :: convertUnits=1.d-9/SECONDS_PER_HOUR
!@param CCE G 2012 canopy environment coefficient to be set
!@+ such that total gamma of 1 results during standard conditions (e.g. 
!@+ 0.3 or 0.57 in CLM4 and WRF-AQ models, respectivesly)
!@+ Is this the same as RHO in MEGAN 2.1?
!TODO: Tune this CCE parameter:
real*8, parameter :: CCE=1.d0
real*8, dimension(n_covertypes) :: pvt0,hvt0 ! ent types and heights
real*8, dimension(nMeganPFT) :: pvt ! locat fraction of MEGAN PFTs
integer, intent(IN) :: i,j
integer :: n, localTimeIndex, hour, dayOfYear, nTracer
integer :: ipft
integer, parameter :: nMeganSpecies=1
type(biogenicSpecies), dimension(nMeganSpecies) :: species
character*80 :: message

! Fill in Megan species we're using, until I learn how to iterate the objects better:
species(1)=Isoprene

call modelEclock%get(dayOfYear=dayOfYear, hour=hour)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Preparations. Gather needed information for MEGAN:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Get near-surface CO2 value:
! ---------------------------

! Getting CO2 from rad code exporting (one option GHY_DRV uses; see
! land_CO2_bc_flag for others. We could, e.g. hook up with CO2 tracer).
CO2_megan=CO2ppm*CO2X ! remains in-routine in case eventually is I,J-dependent


! Get the local Cosine of the Solar Zenith angle:
! ---------------------------------------------

! First trying to use the COSZ1 from the radiation code.
! COSZ1 is the cosine of the solar zenith angle
! averaged over the physics time step. Note I believe they set this to 0 
! when sza > 90 deg, so there will be no "twilight" SZA when sun is below
! the horizon. However, I think this is OK because within the gamma routine,
! gamma is set to zero when COS(SZA) is negative (or zero).
! Note: beta is confusingly called both the SZA and the elevation angle (SZA
! complement) in the megan routines. I believe sinbeta (sine of the elevation angle)
! equals cos(sza). So, to avoid changing the gamma routine, keeping this called
! sin(beta) in there, but passing it cos(sza) here:
cosSZA_megan = cosz1(i,j)

! Get the local julian day:
! -------------------------

! Note that if there is another way to get Ptoa from the rad code or elsewhere
! to send to the gamma_p routine, then we wouldn't need this JDAY_megan.
! Here, localTimeIndex is intended to be an integer ranging from 1 to
! INT_HOURS_PER_DAY. When it is beyond that range, increment (or decrement)
! the local julian day:
localTimeIndex=(hour+1)+NINT((i-(IM+1)/2.)*HOURS_PER_DAY/float(IM))
if(localTimeIndex > HOURS_PER_DAY) then
  JDAY_megan=dayOfYear+1
else if(localTimeIndex < 1) then
  JDAY_megan=dayOfYear-1
else
  JDAY_megan=dayOfYear
end if

! Get the Surface Air Temperature, both instantaneous & running average:
! ----------------------------------------------------------------------

! atmsrf%tsavg should already be in K and is the current average over
! surface types in this (i,j):
SAT_megan=atmsrf%tsavg(i,j)

! For the running average, attempting to pass in the running average
! object and the current value/location for updating.
! If, for this location, the running average has not finished it's
! first averaging period, it will remain undefined and use the instantaneous
! value of the SAT for megan instead:
call running_average( SAT, SAT_megan, i, j)
if(SAT%runningAverage(i,j)==undef)then
  ! i.e. in first averaging period use instantaneous value instead:
  SAT_daily_megan=SAT_megan
else
  SAT_daily_megan=SAT%runningAverage(i,j)
end if

! Get the near-top-of-canopy temperature, running average:
! ----------------------------------------------------------------------
call get_canopy_temperaure(T_megan, i, j) ! from GHY_DRV

! Note that the 'instantaneous' value returned there is only used to get the
! running average (not used on it's own in a call to gamma routines, *except*
! during the first averaging period. See below.)

! This canopy temperature is defined as -1.d30 (not 'undef') when there's no
! vegetation. Hence we have to deal with missing values here. That is very
! tricky when dealing with running averages, as it might be -1.d30 at one
! time but not later... So for safety, I replace canopy temperature with SAT
! when the former is missing (since SAT should be available everywhere):
if(T_megan == -1.d30 .or. T_megan == undef) then
  T_megan=SAT_megan
else
  T_megan=T_megan+tf ! degC --> K
end if

call running_average( T, T_megan, i, j)
if(T%runningAverage(i,j)==undef)then
  ! e.g. in first averaging period use instantaneous value instead:
  T_daily_megan=T_megan
else
  T_daily_megan=T%runningAverage(i,j)
end if

! Get LAI, instantaneous (current) & "previous" value.
! ----------------------------------------------------
! Below I used the "bulk" LAI from Ent. Unclear to me if anything can be
! gained by instead doing a fraction-weighted average of gamma_LAI, by feeding
! the gamma routine PFT-specific LAIs, instead of feeding it bulk LAI and getting
! a gamma from that. E.g. getting a "bulk gamma_LAI" by accumulation, instead
! of starting with the bulk LAI. Things would have to be moved inside loops if
! you try that. 

! Get LAI leaf area index from current and previous time step (This is weird;
! see notes below on questionable methodology and TSTLEN set to 30 days...)
!
! For now, for LAI, the running_average subroutine returns the LAI saved from
! 30 days ago, and this will act as the "previous" value. (It knows to return
! this, instead of the running average because LAI%laggedValue=.true.):

if(fearth(i,j)>0.d0) then
  call ent_get_exports( entcells(i,j),leaf_area_index=LAI_megan)
else
  LAI_megan=0.d0
end if
LAI_current_megan=LAI_megan

call running_average( LAI, LAI_megan, i, j)
if(LAI%runningAverage(i,j)==undef)then
  ! In first averaging period use instantaneous value instead (=no aging yet):
  LAI_previous_megan=LAI_current_megan
else
  ! Then use the value from start of period (not an actual average in this case):
  LAI_previous_megan=LAI%runningAverage(i,j)
end if


! Obtain photosynthetic photon flux density, instantaneous & running average:
! ---------------------------------------------------------------------------
! This is an 'educated guess' obtained from GHY_DRV subroutine earth, but what's
! the 0.82? :
par_total=SRVISSURF(i,j)*cosSZA_megan*0.82d0 ! visible (400-700nm) rad from GHY_DRV
par_direct=par_total*FSRDIR(i,j)
par_diffuse=par_total-par_direct

! Next line takes direct and diffuse PAR (which Ent says are in W m-2
! for the 400-700 nm range) and converts to PPFD in micro-mol(photons) m-2 s-1:
PPFD_megan=par_direct*ConvertSunPPFD+par_diffuse*ConvertShadePPFD

call running_average( PPFD, PPFD_megan, i, j)
if(PPFD%runningAverage(i,j)==undef)then
  ! e.g. in first averaging period use instantaneous value instead:
  PPFD_daily_megan=PPFD_megan
else
  PPFD_daily_megan=PPFD%runningAverage(i,j)
end if

! Map the Ent PFT fractions to the MEGAN PFTs fractions in the box:
! -----------------------------------------------------------------
if(fearth(i,j)>0.d0)then
  call ent_get_exports(entcells(i,j),vegetation_fractions=PVT0)
  call ent_get_exports(entcells(i,j),vegetation_heights=HVT0)
  call map_ent_pfts_to_megan_pfts(pvt0,hvt0,pvt,i,j)
else
  pvt(:)=0.d0
end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the gammas from MEGAN for current conditions:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Gamma for Leaf Area Index (independent of species properties):
call get_gamma_lai( LAI_megan, gamma_LAI )

! Gamma for photosynthetic photon flux density activity
! (independent of species properties):
call get_gamma_p( JDAY_megan, cosSZA_megan, PPFD_megan, &
                & PPFD_daily_megan, gamma_PPFD)

! Gamma for CO2 Inhibition (independent of species properties):
call get_gamma_CO2(CO2_megan, gamma_CO2)

! Gamma for Soil Moisture (independent of species properties):
! Right now just returns a gamma of 1.0:
call get_gamma_s(gamma_SM)

! begin loop over species objects. I.e. below gammas are species-dependant:
tracers_loop: do nTracer=1,ntm

  ! skip if tracer not turned on yet or not intended for megan use:
  if(itime < itime_tr0(nTracer) .or. do_megan(nTracer) <= 0) cycle

  ! try to match tracer with megan-defined species, otherwise skip:
  species_loop: do n=1,size(species)
    if(trim(trname(nTracer))==trim(species(n)%itsname)) then

      ! G 2012 says that CO2 gamma and soil moisture gamma should be non-unity
      ! only for Isoprene. So overwrite here for non-Isoprene species:
      if (trim(species(n)%itsname) .ne. 'Isoprene')then
        gamma_CO2=1.d0
        gamma_SM=1.d0
      end if

      ! Gamma for Aging (species-dependant):
      !   TODO: Potentially-important: the hammoz model noted that the algorithm here is
      !   questionalbe because it assumes the timestep of a month (see TSTLEN parameter set to
      !   30 days) is the same as the timestep used for LAI. This is the reason I set the "previous"
      !   LAI in the call below to be the 30-day lagged value of the Ent LAI. Someone needs to look
      !   into it to see, e.g., if we can use the previous model timestep's LAI instead and then change
      !   TSTLEN (=t) in gamma_a routine to DTsrc (or whatever).

      call get_gamma_a( LAI_previous_megan, LAI_current_megan, T_daily_megan, &
                      & species(n), gamma_AGE)
      ! note that it looks like hammoz passes a daily and monthly LAI (instead of latest
      ! instantaneous one and month-old one)...

      ! Light-dependant temperature gamma: (only one used for Isoprene):
      call get_gamma_tld(SAT_megan, SAT_daily_megan, species(n), gamma_tld)
      ! Light-independant temperature gamma:
      call get_gamma_tli(SAT_megan, species(n), gamma_tli)

      ! Calculate the emissions flux, to be exported and applied elsewhere:

      ! I believe that with the following check, we don't have to treat Isoprene as a
      ! special case of the emissions formula a few lines down -- because the light-
      ! independent portion will drop out! But it does worry me why the hammoz model at least
      ! treated Isoprene separately...
      if (trim(species(n)%itsname) == 'Isoprene') then
        if(species(n)%ldf .ne. 1.d0) call stop_model( &
        & 'Isoprene MEGAN LDF .ne. 1.',255)
      end if

      ! Calculate the bulk emission factor for this species in a loop over fractions
      ! of *MEGAN* (not Ent) plant functional types:
      bulk_EF=0.d0
      do ipft=1,nMeganPFT
        bulk_EF=bulk_EF+species(n)%EF(ipft)*pvt(ipft)
      end do

      ! Calculate the Emissions:

      ! I am aiming for kg m-2 s-1 units for "source". Since EF is in microGram m-2 hr-1
      ! and the gammas are unitless, conversion to kg m-2 s-1 is 1.d-9/DTsrc (see
      ! convertUnits param):
      sfc_src(i,j,nTracer,do_megan(nTracer))= &
      & convertUnits*CCE*bulk_EF*gamma_LAI*gamma_AGE*gamma_SM*gamma_CO2&
      & * ( (1.d0-species(n)%ldf) * gamma_tli + &
      & species(n)%ldf * gamma_PPFD*gamma_tld)

      cycle tracers_loop ! done with this tracer

    end if ! matching megan species to tracer name

  end do species_loop

  write(message,*)'MEGAN species '//trim(trname(nTracer))// 'not found.'
  call stop_model(trim(message),255)

end do tracers_loop


end subroutine biogenicEmissions_drv


subroutine alloc_megan(grid)
!@sum init_megan initialize some properties for MEGAN
!@+ at startup and allocate array dimensions
!@auth Greg Faluvegi

use megan
use megan_objects_mod, only: runningAverage
! Would like to use this, but it's not set yet at this point in code:
! use model_com, only: nday
use domain_decomp_atm, only: dist_grid, getDomainBounds
use constant, only: undef
use model_com, only: dtsrc, calendar
use Rational_mod, only: nint
use dictionary_mod, only : get_param, is_set_param

implicit none

integer :: nday_local
integer :: ier, J_1H, J_0H, I_1H, I_0H
real*8 :: dtsrc_local
type (dist_grid), intent(in) :: grid
call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
call getDomainBounds( grid , I_STRT_HALO=I_0H, I_STOP_HALO=I_1H )

! since MODELE.f didn't set nday yet when this code is called, must
! calculate it here for now: Note that the nint here is used from
! Rational_mod. Not simply the fortran intrinsic:
dtsrc_local = dtsrc
if(is_set_param("DTsrc"))call get_param("DTsrc",dtsrc_local)
nday_local=2*nint(calendar%getSecondsPerDay()/(dtsrc_local*2))


! For the type(runningAverage) :: SAT, T, LAI, PPFD, etc. set the
! daysPerPeriod and stepsPerDay before allocate statements, so they
! can be used as a dimension:
SAT%daysPerPeriod=1  ! we want daily running average
T%daysPerPeriod=1    ! we want daily running average
PPFD%daysPerPeriod=1 ! we want daily running average
LAI%daysPerPeriod=30 ! we want 30-day lagged value?
SAT%stepsPerDay=nday_local ! e.g. DTsrc timesteps in a day
PPFD%stepsPerDay=nday_local ! e.g. DTsrc timesteps in a day
LAI%stepsPerDay=nday_local ! e.g. DTsrc timesteps in a day
T%stepsPerDay=nday_local ! e.g. DTsrc timesteps in a day

! Allocations:

! Allocate running average stuff:

allocate( SAT%first(I_0H:I_1H,J_0H:J_1H) )
allocate( SAT%step(I_0H:I_1H,J_0H:J_1H) )
allocate( SAT%stepSave(I_0H:I_1H,J_0H:J_1H,SAT%stepsPerDay ) )
allocate( SAT%day(I_0H:I_1H,J_0H:J_1H) )
allocate( SAT%dayAvg(I_0H:I_1H,J_0H:J_1H,SAT%daysPerPeriod) )
allocate( SAT%runningAverage(I_0H:I_1H,J_0H:J_1H) )
allocate( SAT%periodRunningSum(I_0H:I_1H,J_0H:J_1H) )
allocate( SAT%marker(I_0H:I_1H,J_0H:J_1H) )

allocate( PPFD%first(I_0H:I_1H,J_0H:J_1H) )
allocate( PPFD%step(I_0H:I_1H,J_0H:J_1H) )
allocate( PPFD%stepSave(I_0H:I_1H,J_0H:J_1H,PPFD%stepsPerDay) )
allocate( PPFD%day(I_0H:I_1H,J_0H:J_1H) )
allocate( PPFD%dayAvg(I_0H:I_1H,J_0H:J_1H,PPFD%daysPerPeriod) )
allocate( PPFD%runningAverage(I_0H:I_1H,J_0H:J_1H) )
allocate( PPFD%periodRunningSum(I_0H:I_1H,J_0H:J_1H) )
allocate( PPFD%marker(I_0H:I_1H,J_0H:J_1H) )

allocate( LAI%first(I_0H:I_1H,J_0H:J_1H) )
allocate( LAI%step(I_0H:I_1H,J_0H:J_1H) )
allocate( LAI%stepSave(I_0H:I_1H,J_0H:J_1H,LAI%stepsPerDay) )
allocate( LAI%day(I_0H:I_1H,J_0H:J_1H) )
allocate( LAI%dayAvg(I_0H:I_1H,J_0H:J_1H,LAI%daysPerPeriod) )
allocate( LAI%runningAverage(I_0H:I_1H,J_0H:J_1H) )
allocate( LAI%periodRunningSum(I_0H:I_1H,J_0H:J_1H) )
allocate( LAI%marker(I_0H:I_1H,J_0H:J_1H) )

allocate( T%first(I_0H:I_1H,J_0H:J_1H) )
allocate( T%step(I_0H:I_1H,J_0H:J_1H) )
allocate( T%stepSave(I_0H:I_1H,J_0H:J_1H,T%stepsPerDay) )
allocate( T%day(I_0H:I_1H,J_0H:J_1H) )
allocate( T%dayAvg(I_0H:I_1H,J_0H:J_1H,T%daysPerPeriod) )
allocate( T%runningAverage(I_0H:I_1H,J_0H:J_1H) )
allocate( T%periodRunningSum(I_0H:I_1H,J_0H:J_1H) )
allocate( T%marker(I_0H:I_1H,J_0H:J_1H) )

! Initializing running average stuff (values may be overwritten by reading from restart files):

SAT%laggedValue=.false. ! do actual running average
SAT%first = .true.
SAT%step = 0.d0
SAT%stepSave = undef
SAT%day = 0.d0
SAT%dayAvg = undef
SAT%runningAverage = undef
SAT%periodRunningSum = undef ! starts at 0 end of first averaging period
SAT%marker = 0

PPFD%laggedValue=.false. ! do actual running average
PPFD%first = .true.
PPFD%step = 0.d0
PPFD%stepSave = undef
PPFD%day = 0.d0
PPFD%dayAvg = undef
PPFD%runningAverage = undef
PPFD%periodRunningSum = undef ! starts at 0 end of first averaging period
PPFD%marker = 0

LAI%laggedValue=.true. ! want lagged value, not running average
LAI%first = .true.
LAI%step = 0.d0
LAI%stepSave = undef
LAI%day = 0.d0
LAI%dayAvg = undef
LAI%runningAverage = undef
LAI%periodRunningSum = undef ! starts at 0 end of first averaging period
LAI%marker = 0

T%laggedValue=.false. ! do actual running average
T%first = .true.
T%step = 0.d0
T%stepSave = undef
T%day = 0.d0
T%dayAvg = undef
T%runningAverage = undef
T%periodRunningSum = undef ! starts at 0 end of first averaging period
T%marker = 0

! Initializing biogenic species stuff (nothing to allocate currently):

isoprene%itsname='Isoprene' ! .le. 8 characters and matching trname please.
isoprene%cceo=2.0d0 ! Coefficient for temperature activity factor in gamma_tld routine
isoprene%ct1=95.0d0 ! A temperature needed for the gamma_tld routine
isoprene%tdf_prm=0.13d0 ! a temperature-dependent parameter needed for gamma_tli routine (beta in G 2012)
isoprene%ldf=1.0d0 ! light dependant fraction, used for relative weighting gammas
isoprene%aindx=5 ! an index to position in arrays Anew, Agro, Amat, Aold for aging gamma
isoprene%ef=(/ 600.d0,     1.d0,  3000.d0, 7000.d0, 10000.d0, & ! emission factors by MEGAN PFT from ...
  &           7000.d0, 10000.d0, 11000.d0, 2000.d0,  4000.d0, & ! ... MGN2MECH/INCLDIR/EFS_PFT.EXT.womap
  &           4000.d0,  1600.d0,   800.d0,  200.d0,    50.d0, &
  &              1.d0  /)
! IMPORTANT: in entering more megan ef values, always confirm index 2 vs. 3 as there
!            was some bug in the code that swapped the two vs. the paper their based on?
!            (see Tables 2 and 3 in G 2012 vs. MGN2MECH/INCLDIR/EFS_PFT.EXT.womap)

!TODO: if a species we are using in the model doesn't happen to line up with
! one of the 20 megan categories, then we'll have to run a mechanism translation to get it
! and this is a major task not yet programmed. E.g. see stuff in the MGN2MECH/ megan dir.

#ifdef PS_BVOC
call stop_model('DO_MEGAN + PS_BVOC conflict',255)
#endif
#ifdef BIOGENIC_EMISSIONS
call stop_model('DO_MEGAN + BIOGENIC_EMISSIONS conflict',255)
#endif

end subroutine alloc_megan


subroutine running_average(this, val, i, j)
!@sum running_average keeps a running average of model quantities
!@+ needed for MEGAN input. In practice, this does hourly and daily
!@+ running averages and uses those to maintain the period-long running
!@+ average, to avoid saving a huge array. This multi-day functionality
!@+ was put in in case we need it, but right now in MEGAN I think only
!@+ daily averages are requested, so the period will be a single day.
! Note that there is a second option here, controlled by the "laggedValue"
! property of the runningAverage type object, that will return instead
! of a running average, the value (val) at i,j from X-days ago...
!@auth Greg Faluvegi

use megan
use constant, only: undef
use megan_objects_mod, only: runningAverage

implicit none

!@var this current pointed-to running average object
!@var val passed in value to be used to update the running average
!@var i longitude index to update the running average
!@var j latitude index to update the running average
type(runningAverage), intent(inout) :: this
real*8, intent(in) :: val
integer, intent(in) :: i,j
!@var byStepsPerDay = reciprocal of expected calls per day
!@var byDaysPerPeriod = reciprocal of expected days per period
!@var currentDayAverage temporary working variable
real*8 :: byStepsPerDay, byDaysPerPeriod, currentDayAverage
integer :: n

!TODO: Greg: sanity check by printing things out to confirm that the whole
! globe (all i,j's) running average are "in sync". Meaning their 'step(i,j)'s
! are incrementing each call all over and whether their 'first(i,j)'s turn
! to false at the same time. If so, I guess such arrays could be made 
! scalar (though then you'd have to increment them and set them high up in
! the model - outseide the I,J loop.)

! Make sure we didn't run outside of expected calls per day:
if(nint(this%step(i,j)) < 0 .or. &
 & nint(this%step(i,j)) > this%stepsPerDay) then
  write(6,*) "i,j,step,max=",i,j,nint(this%step(i,j)),this%stepsPerDay
  call stop_model('step problem in running_average MEGAN',255)
end if

byStepsPerDay=1.d0/dble(this%stepsPerDay)
byDaysPerPeriod=1.d0/dble(this%daysPerPeriod)

! Increment the "step", and save current value in an array that
! saves all steps for later:
this%step(i,j) = this%step(i,j) + 1.d0
this%stepSave(i,j,nint(this%step(i,j)))=val

! If it is not at the end of the days' worth of accumulations,
! we are done! Otherwise, must continue:
if(nint(this%step(i,j)) == this%stepsPerDay) then

  ! Reset the steps:
  this%step(i,j) = 0.d0

  ! If we're still in the first averaging period, update the
  ! daily running average and the period running sum, but
  ! return an undefined running average and let the calling routine
  ! deal with that:
  if(this%first(i,j)) then

    ! increment the day:
    this%day(i,j) = this%day(i,j) + 1.d0

    ! calculate today's average, storing it in daily average array:
    currentDayAverage=0.d0
    do n=1,this%stepsPerDay
      currentDayAverage=currentDayAverage+this%stepSave(i,j,n) &
      & * byStepsPerDay
    end do
    this%dayAvg(i,j,nint(this%day(i,j)))=currentDayAverage

    ! we are not ready to define this yet:
    this%runningAverage(i,j)=undef

    ! if today is the end of the first period do special things:
    if(nint(this%day(i,j)) == this%daysPerPeriod) then

      this%periodRunningSum(i,j) = 0.d0
      do n=1,this%daysPerPeriod
        this%periodRunningSum(i,j)=this%periodRunningSum(i,j) &
        & + this%dayAvg(i,j,n)
      end do

      ! first time defining running average:
      if(.NOT.this%laggedValue)then
        ! normal running average case:
        this%runningAverage(i,j) = this%periodRunningSum(i,j) &
      & *byDaysPerPeriod
      else
        ! lagged value case (return oldest daily value):
        this%runningAverage(i,j)=this%dayAvg(i,j,1)
      end if

      ! no longer first averaging period, so record that fact in
      ! 'first' logical, initialize the marker to 0 (will increment
      ! to 1 the next time step; marker keeps track of position in the
      ! current period's running sum) and reset the day counter to 0
      ! (not sure that is needed):
      this%first(i,j)=.false.
      this%day(i,j)=0.d0
      this%marker(i,j)=0

    end if ! end of special stuff at end of first period

  else ! no longer in first averaging period; update the running average:

    ! move the marker (pointer to position in array that is current day):
    this%marker(i,j)=this%marker(i,j)+1

    ! if past the end of period, cycle marker around:
    if(this%marker(i,j) == this%daysPerPeriod+1) this%marker(i,j)=1

    ! calculate today's average:
    currentDayAverage=0.d0
    do n=1,this%stepsPerDay
      currentDayAverage=currentDayAverage+this%stepSave(i,j,n) &
      &*byStepsPerDay 
    end do

    ! Only for special case where we want to return a lagged 
    ! (daily average) value instead of a running average, define
    ! that now:
    if(this%laggedValue)then
      this%runningAverage(i,j)=this%dayAvg(i,j,this%marker(i,j))
    end if
    ! now continuing on...

    ! remove previous daily average that was stored in today's position
    ! from the period running sum:
    this%periodRunningSum(i,j) = this%periodRunningSum(i,j) &
    & - this%dayAvg(i,j,this%marker(i,j))

    ! Replace daily average stored in today's position with today's value.
    ! I.e. saving it so it can be removed next time we come arround to the
    ! above line of code:
    this%dayAvg(i,j,this%marker(i,j))=currentDayAverage

    ! add today's average to the period running sum:
    this%periodRunningSum(i,j) = this%periodRunningSum(i,j) &
    & + currentDayAverage

    ! Finally, update the running average (but ignore for the
    ! special case of desiring a lagged value instead):
    if(.NOT.this%laggedValue)then
      this%runningAverage(i,j) = this%periodRunningSum(i,j) &
      & *byDaysPerPeriod
    end if

  end if ! whether or not first averaging period

end if ! whether of not at the end of the day

end subroutine running_average


! Moving on to gamma routines:

! MEGAN notes:
!-----------------------------------------------------------------------
!.....1) Calculate GAM_L (GAMMA_LAI)
!-----------------------------------------------------------------------
!                            0.49[LAI]
!             GAMMA_LAI = ----------------    (non-dimension)
!                         (1+0.2LAI^2)^0.5
!
!     SUBROUTINE GAMMA_LAI returns the GAMMA_LAI values
!-----------------------------------------------------------------------

subroutine get_gamma_lai(lai,gam_l)
!@sum Calculate gamma for leaf area index from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
implicit none
real*8,intent(IN) :: lai
real*8,intent(OUT) :: gam_l

gam_l = (0.49d0*lai) / SQRT(1.d0 + 0.2d0*(lai**2))

return
end subroutine get_gamma_lai


! MEGAN notes:
!-----------------------------------------------------------------------
!.....2) Calculate GAM_P (GAMMA_P)
!-----------------------------------------------------------------------
!             GAMMA_P = 0.0         a<=0, a>=180, sin(a) <= 0.0
!
!             GAMMA_P = sin(a)[ 2.46*(1+0.0005(Pdaily-400))*PHI - 0.9*PHI^2 ]
!                                   0<a<180, sin(a) > 0.0
!           where PHI    = above canopy PPFD transmission (non-dimension)
!                 Pdaily = daily average above canopy PPFD (umol/m2s)
!                 a      = solar angle (degree)
!
!                 Note: AAA = 2.46*BBB*PHI - 0.9*PHI^2
!                       BBB = (1+0.0005(Pdaily-400))
!                       GAMMA_P = sin(a)*AAA
!
!                       Pac
!             PHI = -----------
!                   sin(a)*Ptoa
!           where Pac  = above canopy PPFD (umol/m2s)
!                 Ptoa = PPFD at the top of atmosphere (umol/m2s)
!
!             Pac =  SRAD * 4.766 mmmol/m2-s * 0.5
!
!             Ptoa = 3000 + 99*cos[2*3.14-( DOY-10)/365 )]
!           where DOY = day of year
!
!     SUBROUTINE GAMMA_P returns the GAMMA_P values
!-----------------------------------------------------------------------

subroutine get_gamma_p( local_jday, sinbeta, ppfd, d_ppfd, gam_p )
!@sum Calculate gamma for photosynthetic photon flux density activity 
!@+ from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
! This routine used to calculate the zenith angle itself; we want 
! to use the GCM's instead.

use constant, only : twopi, pi, radian
use TimeConstants_mod, only: INT_DAYS_PER_YEAR

implicit none 

!@var local_jday integer julian day at current lat, lon. DAY in MEGAN2.1
integer, intent(IN) :: local_jday
!@var ppfd instantaneous photosynthetic photon flux density
!@var d_ppfd daily photosynthetic photon flux density
real*8, intent(IN) :: ppfd, d_ppfd
real*8, intent(OUT) :: gam_p
!@var sinbeta sine of solar elevation angle (=COS(SZA))
!@var Ptoa photosynthetic photon flux density at top of atmosphere (umol m-2 s-1)
!@var Pac photosynthetic photon flux density above the canopy (umol m-2 s-1)
!@var phi non-dimensional above canopy photosynthetic photon flux density transmission
real*8 :: aaa,bbb,sinbeta,Ptoa=0.d0,Pac,phi=0.d0

Pac = PPFD


if(sinbeta <= 0.d0) then
  gam_p = 0.d0
else if(sinbeta > 0.d0) then
  ! TODO: See if any more of these "magic numbers" are actually physical
  !       constants that could be replaced (instead of just parameters of
  !       parameterization. COS( ) looks to be below is COS(hourAngle).
  Ptoa = 3000.d0 + 99.d0 * COS(twopi*(local_jday-10)/INT_DAYS_PER_YEAR)
  phi = Pac/(sinbeta*Ptoa)
  bbb = 1.d0 + 0.0005d0*( d_ppfd-400.d0 )
  aaa = ( 2.46d0 * bbb * phi ) - ( 0.9d0 * phi**2 )
  gam_p = sinbeta * aaa
else
  call stop_model('invalid sza in MEGAN gamm_p routine',255)
end if

! MEGAN says: "Screening the unforced errors": If solar elevation
! angle is less than 1, gamma_p can not be greater than 0.1.
! (but I think they mean zenith, not elevation angle):
if(asin(sinbeta)/radian < 1.d0 .and. gam_p > 0.1d0) gam_p = 0.d0

return
end subroutine get_gamma_p


! MEGAN Notes:
!-----------------------------------------------------------------------
!.....3) Calculate GAM_T (GAMMA_T) for isoprene
!-----------------------------------------------------------------------
!                          Eopt*CT2*exp(CT1*x)
!             GAMMA_T =  ------------------------
!                        [CT2-CT1*(1-exp(CT2*x))]
!           where x      = [ (1/Topt)-(1/Thr) ] / 0.00831
!                 Eopt   = 1.75*exp(0.08(Tdaily-297)
!                 CT1    = 80
!                 CT2    = 200
!                 Thr    = hourly average air temperature (K)
!                 Tdaily = daily average air temperature (K)
!                 Topt   = 313 + 0.6(Tdaily-297)
!
!                 Note: AAA = Eopt*CT2*exp(CT1*x)
!                       BBB = [CT2-CT1*(1-exp(CT2*x))]
!                       GAMMA_T = AAA/BBB
!
!     SUBROUTINE GAMMA_TLD returns the GAMMA_T value for isoprene
!-----------------------------------------------------------------------

subroutine get_gamma_tld(temp,d_temp,this,gam_t)
!@sum Calculate gamma temperature response factor for Isopene 
!@+ (tld=light dependent?) from MEGAN2.1
! (This is more complex in G 2012 incorporating also a 240-hour T)
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
! In MEGAN it uses INCLUDE 'EACO.EXT' and uses its INDEX1 function to look
! those parameters up. Here, we pass in 'this' biogenicSpecies object
! containing the needed information.

use megan_objects_mod, only: biogenicSpecies
use constant, only: bygasc
implicit none
!@var this current pointed-to species object
type(biogenicSpecies), intent(inout) :: this
!@var temp instantaneous surface air temperature (K)
!@var d_temp daily average surface air temperature (K)
real*8, intent(IN) :: temp,d_temp
real*8, intent(OUT) :: gam_t
real*8, parameter :: ct2=200.d0 ! some reference temperature (K)?
real*8, parameter :: tempS=297.d0 ! standard leaf temperature
real*8 :: Eopt, Topt, X, aaa, bbb

Eopt = this%cceo * exp(0.08d0*(d_temp-tempS))
Topt = 313.d0 + ( 0.6d0*(d_temp-tempS) )
X = ( (1.d0/Topt)-(1.d0/temp) ) *bygasc*1.d3 ! bygasc*1.d3 was 1/0.00831d0 (i.e. mol K / kJ)
aaa = Eopt*ct2*exp(this%ct1*X)
bbb = ( ct2-this%ct1*( 1.d0-exp(ct2*X) ) )
gam_t = aaa/bbb

return
end subroutine get_gamma_tld


! MEGAN Notes:
!-----------------------------------------------------------------------
!.....4) Calculate GAM_T (GAMMA_T) for non-isoprene
!-----------------------------------------------------------------------
!
!             GAMMA_T =  exp[TDP_FCT*(T-Ts)]
!           where TDP_FCT = temperature dependent parameter ('beta')
!                 Ts     = standard temperature (normally 303K, 30C)
!
!     SUBROUTINE GAMMA_TLI returns the GAMMA_T value for non-isoprene
!-----------------------------------------------------------------------

subroutine get_gamma_tli(temp,this,gam_t)
!@sum Calculate gamma temperature response factor for non-Isopene 
!@+ species. (tli=light independent?) from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
use megan_objects_mod, only: biogenicSpecies
implicit none
!@var this current pointed-to species object
type(biogenicSpecies), intent(inout) :: this
!@var temp instantaneous surface air temperature (K)
!@param Ts standard temperature (K)
! Note this%tdf_prm is a temperature dependent parameter for current
! species. In MEGAN it's gotten using INCLUDE 'EACO.EXT' and the INDEX1 
! function. Here, we pass in the biogenicSpecies object containing that
! information. Note, I think this is called beta in G 2012.
real*8, intent(IN) :: temp
real*8, intent(OUT) :: gam_t
real*8, parameter :: Ts = 303.0d0

gam_t = exp( this%tdf_prm*(temp-Ts) )

return
end subroutine get_gamma_tli
!-----------------------------------------------------------------------


! MEGAN Notes:
!-----------------------------------------------------------------------
!.....5) Calculate GAM_A (GAMMA_age)
!-----------------------------------------------------------------------
!
!             GAMMA_age = Fnew*Anew + Fgro*Agro + Fmat*Amat + Fold*Aold
!           where Fnew = new foliage fraction
!                 Fgro = growing foliage fraction
!                 Fmat = mature foliage fraction
!                 Fold = old foliage fraction
!                 Anew = relative emission activity for new foliage
!                 Agro = relative emission activity for growing foliage
!                 Amat = relative emission activity for mature foliage
!                 Aold = relative emission activity for old foliage
!
!             For foliage fraction
!             Case 1) LAIc = LAIp
!             Fnew = 0.0  , Fgro = 0.1  , Fmat = 0.8  , Fold = 0.1
!
!             Case 2) LAIp > LAIc
!             Fnew = 0.0  , Fgro = 0.0
!             Fmat = 1-Fold
!             Fold = (LAIp-LAIc)/LAIp
!
!             Case 3) LAIp < LAIc
!             Fnew = 1-(LAIp/LAIc)                       t <= ti
!                  = (ti/t) * ( 1-(LAIp/LAIc) )          t >  ti
!
!             Fmat = LAIp/LAIc                           t <= tm
!                  = (LAIp/LAIc) +
!                      ( (t-tm)/t ) * ( 1-(LAIp/LAIc) )  t >  tm
!
!             Fgro = 1 - Fnew - Fmat
!             Fold = 0.0
!
!           where
!             ti = 5 + (0.7*(300-Tt))                   Tt <= 303
!                = 2.9                                  Tt >  303
!             tm = 2.3*ti
!
!             t  = length of the time step (days)
!             ti = number of days between budbreak and the induction of
!                  emission
!             tm = number of days between budbreak and the initiation of
!                  peak emissions rates
!             Tt = average temperature (K) near top of the canopy during
!                  current time period (daily ave temp for this case)
!
!
!             For relative emission activity
!             Case 1) Constant
!             Anew = 1.0  , Agro = 1.0  , Amat = 1.0  , Aold = 1.0
!
!             Case 2) Monoterpenes
!             Anew = 2.0  , Agro = 1.8  , Amat = 0.95 , Aold = 1.0
!
!             Case 3) Sesquiterpenes
!             Anew = 0.4  , Agro = 0.6  , Amat = 1.075, Aold = 1.0
!
!             Case 4) Methanol
!             Anew = 3.0  , Agro = 2.6  , Amat = 0.85 , Aold = 1.0
!
!             Case 5) Isoprene
!             Anew = 0.05 , Agro = 0.6  , Amat = 1.125, Aold = 1.0
!
!     SUBROUTINE GAMMA_A returns GAMMA_A
!-----------------------------------------------------------------------

subroutine get_gamma_a(LAIp,LAIc,Tt,this,gam_a)
!@sum Calculate gamma foliage aging factor from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
! MEGAN uses INCLUDE 'EACO.EXT' and function INDEX1 to look 
! the relative emissions activity parameter up. See REA_INDEX( )
use megan_objects_mod, only: biogenicSpecies
implicit none
!@var this current pointed-to species object
type(biogenicSpecies), intent(inout) :: this
real*8, parameter :: TSTLEN=30.d0 ! from MEGAN 30 day time step
!@var t MEGAN time step. Set to TSTLEN to make origin clear.
!@param Ts standard temperature (K)
! next line could be a coincidence. It was hardcoded 303 in the aging routine,
! but called "standard temperature" in the gamma_tli routine:
real*8, parameter :: Ts = 303.0d0
real*8 :: t ! for some reason was an integer in MEGAN
!@var LAIp leaf area index at previous time step
!@var LAIc leaf area index at current time step
!@var k the index to hold the current species aindx, the relative emission 
!@+ acitivity index for current species, used to choose Anew, Agro, Amat, Aold
!@var Tt MEGAN daily average temperature (K). Did not bother with intermediate d_temp var
real*8, intent(IN) :: Tt,LAIp,LAIc
real*8, intent(OUT) :: gam_a
!@var ti number of days between budbreak and induction of emission
!@var tm number of days between budbreak and initiation of peak emissions rates
!@var Fnew new foliage fraction
!@var Fgro growing foliage fraction
!@var Fmat mature foliage fraction
!@var Fold old foliage fraction
real*8 :: Fnew, Fgro, Fmat, Fold, ti, tm

! MEGAN notes and values copied from EACO.EXT:
! C=======================================================================
! C  REL_EM_ACT.EXT
! C  This include file contains "produciton and loss within canopy"
! C  factors.
! C
! C
! C  MEGAN v2.02
! C  INPUT version 200
! C
! C  History:
! C  Who          When       What
! C  ---------------------------------------------------------------------
! C  Tan          12/02/06 - Creates this file
! C  Tan          08/14/07 - Move from MEGAN v2.0 to MEGAN v2.02 with no update.
! C=======================================================================
!@var Anew relative emission activity for new foliage
!@var Agro relative emission activity for growing foliage
!@var Amat relative emission activity for mature foliage
!@var Aold relative emission activity for old foliage
integer, parameter :: N_CAT  = 5
real*8, dimension(N_CAT) :: Anew, Agro, Amat, Aold
integer :: k
! Note that, in the comments from MEGAN above, it's the Aold column that's
! all 1.00's, but in the below (which I got from EACO.EXT) it's the Amat that
! are all 1.00's. I think hammoz was updated to G 2012 paper, and has 6 instead
! of 5 categories, and has the Amat all 1's as well, as also stated in that paper:
data Anew(1),Agro(1),Amat(1),Aold(1) /1.00d0, 1.00d0, 1.00d0, 1.00d0/
data Anew(2),Agro(2),Amat(2),Aold(2) /2.00d0, 1.80d0, 1.00d0, 1.05d0/
data Anew(3),Agro(3),Amat(3),Aold(3) /0.40d0, 0.60d0, 1.00d0, 0.95d0/
data Anew(4),Agro(4),Amat(4),Aold(4) /3.50d0, 3.00d0, 1.00d0, 1.20d0/
data Anew(5),Agro(5),Amat(5),Aold(5) /0.05d0, 0.60d0, 1.00d0, 0.90d0/
! end of section taken from EACO.EXT

t = TSTLEN

!... Calculate foliage fraction

if(LAIp < LAIc) then ! i.e. growing:

!        Calculate ti and tm
  if(Tt <= Ts) then
    ti = 5.0d0 + 0.7d0*(300.d0-Tt)
  else
    ti = 2.9d0
  end if
  tm = 2.3d0*ti

!       Calculate Fnew and Fmat, then Fgro and Fold

!       Fnew
  if(LAIc==0.d0)call stop_model("LAIc div by 0, MEGAN gamma_a",255)
  ! note, currently t is set to non-zero TSTLEN so can not cause problem 
  ! in divides by t below:
  if(ti >= t) then
    Fnew = 1.d0 - (LAIp/LAIc)
  else
    Fnew = (ti/t) * ( 1.d0-(LAIp/LAIc) )
  end if

!       Fmat
  if(tm >= t) then
    Fmat = LAIp/LAIc
  else
    Fmat = (LAIp/LAIc) + ( (t-tm)/t ) * ( 1.d0-(LAIp/LAIc) )
  end if

  Fgro = 1.d0 - Fnew - Fmat
  Fold = 0.d0

else if(LAIp == LAIc) then ! i.e. mature

  Fnew = 0.d0
  Fgro = 0.1d0
  Fmat = 0.8d0
  Fold = 0.1d0

else ! LAIp > LAIc          i.e. senescing

  if(LAIp==0.d0)call stop_model("LAIp div by 0, MEGAN gamma_a",255)
  Fnew = 0.d0
  Fgro = 0.d0
  Fold = ( LAIp-LAIc ) / LAIp
  Fmat = 1.d0-Fold

end if

! Ready to calculate the gamma for aging:
k=this%aindx
gam_a=Fnew*Anew(k)+Fgro*Agro(k)+Fmat*Amat(k)+Fold*Aold(k)

return
end subroutine get_gamma_a


! MEGAN Notes:
!-----------------------------------------------------------------------
!.....6) Calculate GAM_SMT (GAMMA_SM)
!-----------------------------------------------------------------------
!
!             GAMMA_SM =     1.0   (non-dimension)
!
!
!     SUBROUTINE GAMMA_S returns the GAMMA_SM values
!-----------------------------------------------------------------------

subroutine get_gamma_s(gam_s)
!@sum Calculate gamma soil moisture response factor
!@+ from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
implicit none
real*8, intent(OUT) :: gam_s

gam_s = 1.d0

! While MEGAN 2.1 sets this to unity, as above, we could parameterize like G 2012.
! We would need to pass in the volumetric soil moisture (m3 m-3), theta, and
! define theta_w, the wilting point soil moisture (below which plants can't extract
! moisture). Theta_1 and delta_theta_1 would follow like the code commented below.
! Note that the wilting point should be set at what's appropriate for *our* ground
! hydrology model/Ent:
!     subroutine gamma_s(theta,gam_s)
!     real*8, intent(IN) :: theta ! volumetric soil moisture
!     ! Next line are placeholders from hamoz, G 2012 has 0.04 for delta_theta_1
!     real*8, parameter :: theta_w=0.35d0, delta_theta_1=0.06d0
!     real*8, parameter :: theta_1=theta_w+delta_theta_1
!     real*8 :: theta_1
!     if(theta >= theta_1) then
!       gam_s=1.d0
!     else if (theta <= theta_w) then
!       gam_s=0.d0
!     else
!       gam_s=(theta-theta_w)/delta_theta_1
!     end if

return
end subroutine get_gamma_s


!  MEGAN Notes:
!-----------------------------------------------------------------------
!.....7) Calculate GAM_CO2(GAMMA_CO2)
!-----------------------------------------------------------------------
!
!             GAMMA_CO2 =     1.0   (non-dimension)
!             When CO2 =400ppm
!
!     SUBROUTINE GAM_CO2 returns the GAMMA_CO2 values
!    Xuemei Wang-2009-06-22 
!-----------------------------------------------------------------------

subroutine get_gamma_CO2(CO2,gam_CO2)
!@sum Calculate gamma CO2 factor from MEGAN2.1
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
implicit none
real*8, intent(IN) :: CO2 ! I think ppmv
real*8, intent(OUT):: gam_CO2
real*8 :: Ci
real*8, parameter :: ISmax = 1.344d0, h=1.4614d0, Cstar =585.d0

Ci = 0.7d0*CO2
if(CO2==400.d0)then
  gam_CO2 = 1.0d0
else
  gam_CO2 = ISmax - ((ISmax*Ci**h)/(Cstar**h+Ci**h))
end if

return
end subroutine get_gamma_CO2


!-----------------------------------------------------------------------
!.....8) Calculate GAMMA_LAIbidir(gam_LAIbidir,LAI)
!-----------------------------------------------------------------------
!From Alex Guenther 2010-01-26
!If lai < 2 Then
!gammaLAIbidir= 0.5 * lai
!ElseIf lai <= 6 Then
!gammaLAIbidir= 1 - 0.0625 * (lai - 2)
!Else
!gammaLAIbidir= 0.75
!End If
!-----------------------------------------------------------------------
subroutine get_gamma_LAIbidir(lai,gam_l)
!@sum Calculate gamma for leaf area index from MEGAN2.1; but for bidirectional VOCs
!@auth MEGAN team, initial modelE implementation by Greg Faluvegi
! It's not clear why this is never called in MEGAN2.1 (nor hammoz), but
! including it here, in case we decide to call it for VOCs classified as bidirectional:
implicit none
real*8,intent(IN) :: lai
real*8,intent(OUT) :: gam_l

if(lai < 2.d0) then
  gam_l = 0.5d0*LAI
else if(lai <= 6.d0 .and. lai >= 2.d0) then
  gam_l = 1.d0 - 0.0625d0 * (lai - 2.d0)
else
  gam_l = 0.75d0
end if

return
end subroutine get_gamma_LAIbidir


! END OF GAMMA ROUTINES


subroutine map_ent_pfts_to_megan_pfts(v_ent,h_ent,v_megan,i,j)
!@sum Map Ent plant functional type coverage fractions onto MEGAN types.
!@+ Similar "hack" by Y. Kim was in radiation code and flammability
!@+ codes to map to GISS vegetation.
!@auth initial modelE implementation by Greg Faluvegi
use ent_com, only: n_covertypes
use megan_objects_mod, only: nMeganPFT
use geom, only: lat2D_dg
implicit none
real*8, dimension(n_covertypes), intent(in) :: h_ent,v_ent
real*8, dimension(nMeganPFT), intent(out) :: v_megan
integer, intent(in) :: i,j
logical :: tropical, temperate, boreal, arctic

! TODO: I am simply putting in the most basic, 'semi-reasonable' mapping
! here. The whole routine should be looked at by an expert in the vegetation.
! E.g. things should be decided by temperature climatologies [or possibly
! MEGAN ecoregion maps], not by latitude as I have done here. Note that h_ent
! (vegetation heights) are passed in, in case useful in refining this, but not
! used in current mapping:
! The person who refines this might want to look at the CLM4 model
! (http://www.cesm.ucar.edu/) pft definitions, as that's what the MEGAN
! ones (basically) are.
! Crops are particularly sketchy in my quick mapping. For example hammoz has
! crop1 and then crop2, though the value for EF of 1 for isoprene seems like 16
! should be crop1?

! Only mapped for a single configuration so far:
if(n_covertypes.ne.18)call stop_model( &
& 'n_covertypes.ne.18 in map_ent_pfts_to_megan_pfts',255)
if(nMeganPFT.ne.16)call stop_model( &
& 'nMeganPFT.ne.16 in map_ent_pfts_to_megan_pfts',255)

! Determine latitude zone as a very rough(!) determination of 
! temperate, tropical, boreal & arctic vegetation:
! Please make only one of the four true!
if(abs(lat2D_dg(i,j))>=60.d0)then
  arctic=.true.  ; boreal=.false. ; temperate=.false. ; tropical=.false.
else if(abs(lat2D_dg(i,j))<60.d0.and.abs(lat2D_dg(i,j))>=50.d0)then
  arctic=.false. ; boreal=.true.  ; temperate=.false. ; tropical=.false.
else if(abs(lat2D_dg(i,j))<50.d0.and.abs(lat2D_dg(i,j))>=23.d0)then
  arctic=.false. ; boreal=.false. ; temperate=.true.  ; tropical=.false.
else if(abs(lat2D_dg(i,j))<23.d0)then
  arctic=.false. ; boreal=.false. ; temperate=.false. ; tropical=.true.
else
  call stop_model( &
  & 'bad lat band definition in map_ent_pfts_to_megan_pfts.',255)
end if

!  I *believe* the 16 Megan types are:
!  =====================================
!  1 needleleaf evergreen temperate tree
!  2 needleleaf deciduous boreal tree   - swapped order compared to G 2012
!  3 needleleaf evergreen boreal tree   - swapped order compared to G 2012
!  4 broadleaf evergreen tropical tree
!  5 broadleaf evergreen temperate tree
!  6 broadleaf deciduous tropical tree
!  7 broadleaf deciduous temperate tree
!  8 broadleaf deciduous boreal tree
!  9 broadleaf evergreen temperate shrub
! 10 broadleaf deciduous temperate shrub
! 11 broadleaf deciduous boreal shrub
! 12 arctic C3 grass
! 13 cool C3 grass
! 14 warm C4 grass (though I think megan *code* says C3 warm... *paper* says this)
! 15 crop1 (?) (*paper* has no 16 and 15=Crop1)
! 16 crop2 (corn??)

! I *believe* the 18 Ent types are:
! ===================================
!  1 evergreen broadleaf early successional
!  2 evergreen broadleaf late successional
!  3 evergreen needleleaf early successional
!  4 evergreen needleleaf late successional
!  5 cold deciduous broadleaf early successional
!  6 cold deciduous broadleaf late successional
!  7 drought deciduous broadleaf
!  8 deciduous needleleaf
!  9 cold adapted shrub
! 10 arid adapted shrub
! 11 C3 grass perennial
! 12 C4 grass
! 13 C3 grass - annual
! 14 arctic C3 grass
! 15 C4 crops
! 16 crops broadleaf woody
! 17 like sand
! 18 like dirt

v_megan(:)=0.d0

! Depending on latitude band, placing Ent:
!  1 evergreen broadleaf early successional
!  2 evergreen broadleaf late successional
! in Megan categories:
!  4 broadleaf evergreen tropical tree
!  5 broadleaf evergreen temperate tree

if(tropical)then
  v_megan(4)=v_megan(4)+v_ent(1)+v_ent(2)
else ! extra-tropical
  v_megan(5)=v_megan(5)+v_ent(1)+v_ent(2)
end if

! Depending on latitude band, placing Ent:
!  3 evergreen needleleaf early successional
!  4 evergreen needleleaf late successional
! in Megan categories:
!  1 needleleaf evergreen temperate tree
!  3 needleleaf evergreen boreal tree
if(boreal .or. arctic)then
  v_megan(3)=v_megan(3)+v_ent(3)+v_ent(4)
else ! temperate or tropical
  v_megan(1)=v_megan(1)+v_ent(3)+v_ent(4)
end if

! Depending on latitude band, placing Ent:
!  5 cold deciduous broadleaf early successional
!  6 cold deciduous broadleaf late successional
!  7 drought deciduous broadleaf
! Into Megan categories:
!  6 broadleaf deciduous tropical tree
!  7 broadleaf deciduous temperate tree
!  8 broadleaf deciduous boreal tree
if(boreal .or. arctic)then
  v_megan(8)=v_megan(8)+v_ent(5)+v_ent(6)+v_ent(7)
else if(temperate)then
  v_megan(7)=v_megan(7)+v_ent(5)+v_ent(6)+v_ent(7)
else if(tropical)then
  v_megan(6)=v_megan(6)+v_ent(5)+v_ent(6)+v_ent(7)
end if

! Place Ent: {8 deciduous needleleaf} into Megan {2 needleleaf deciduous boreal tree}
v_megan(2)=v_megan(2)+v_ent(8)

! Depending on latitude band, placing Ent:
!  9 cold adapted shrub
! Into Megan categories:
! 10 broadleaf deciduous temperate shrub
! 11 broadleaf deciduous boreal shrub
if(boreal .or. arctic)then
  v_megan(11)=v_megan(11)+v_ent(9)
else
  v_megan(10)=v_megan(10)+v_ent(9)
end if

! Place Ent: {10 arid adapted shrub} into Megan:
!  9 broadleaf evergreen temperate shrub
! though this is hardly obvious:
v_megan(9)=v_megan(9)+v_ent(10)

! Place Ent: 
! 11 C3 grass perennial
! 13 C3 grass - annual
! Into Megan category:
! 13 cool C3 grass
! OK?
v_megan(13)=v_megan(13)+v_ent(11)+v_ent(13)

! Place Ent: {12 C4 grass} into Megan:
! 14 warm C4 grass (though I think megan notes says C3 warm...?)
v_megan(14)=v_megan(14)+v_ent(12)

! Place Ent: {14 arctic C3 grass} into Megan:
! 12 arctic C3 grass
v_megan(12)=v_megan(12)+v_ent(14)

! Place Ent: {15 C4 crops} into Megan:
! 16 crop2 (corn??). Simply because I think corn is C4. 
v_megan(16)=v_megan(16)+v_ent(15)

! Place Ent: {16 crops broadleaf woody} into Megan:
! 15 crop1 (?). But for no good reason. 
! Vegetation height criterion could be used to determine these woody crops.
v_megan(15)=v_megan(15)+v_ent(16)

! End categories 17 and 18 are bare sand/dirt, so no
! accumulation into a MEGAN type. So we'ere done.
end subroutine map_ent_pfts_to_megan_pfts
