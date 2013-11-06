!     SCM_COM.f
!@sum set up variables for SCM forcing data 
!@auth  Audrey Wolf
!
C--------------------------------------------------------------------------------
      Module SCMCOM
      USE RESOLUTION , ONLY : LM
      IMPLICIT NONE
      SAVE

C SCM DATA as provided from ARM variational analysis 
!@var SG_P Pressure at GCM sigma levels (mb)
      REAL*8 SG_P(LM)
      REAL*8 SGE_P(LM+1)
!@var SG_HGT Height of layer at GCM sigma levels (m)
      REAL*8 SG_HGT(LM)
!@var SG_T ARM Temperature at GCM sigma levels (K)
      REAL*8 SG_T(LM)
!@var SG_Q ARM Specific Humidity at GCM sigma levels (g/kg)
      REAL*8 SG_Q(LM)
!@var SG_U ARM U wind at GCM Sigma levels (m/s)
      REAL*8 SG_U(LM)
!@var SG_V ARM V wind at GCM sigma levels (m/s)
      REAL*8 SG_V(LM)
!@var SG_OMEGA ARM Omega at GCM sigma levels (mb/hr)
      REAL*8 SG_OMEGA(LM)
!@var SG_WINDIV ARM Wind Divergence at GCM sigma levels (1/s)
      REAL*8 SG_WINDIV(LM)
!@var SG_CONV   convergence as calculated in FCONV subr from windiv for clouds   
      REAL*8 SG_CONV(LM)   
!@var SG_HOR_TMP_ADV ARM Horizontal temperature advection at GCM sigma levels (K/s)
      REAL*8 SG_HOR_TMP_ADV(LM)
!@var SG_VER_S_ADV ARM Vertical S advection at GCM sigma levels (K/s)
      REAL*8 SG_VER_S_ADV(LM)
!@var SG_HOR_Q_ADV ARM Horizontal Q advection at GCM sigma levels( kg/kg/s)
      REAL*8 SG_HOR_Q_ADV(LM)
!@var SG_VER_Q_ADV ARM Vertical Q advection at GCM sigma levels (kg/kg/s) 
      REAL*8 SG_VER_Q_ADV(LM)
!@var SG_ARSCL ARSCL Cloud amounts at GCM sigma levels (%)
      REAL*8 SG_ARSCL(LM)
!@var ASTIME ARM TIME STAMP - surface data 
      REAL*8 ASTIME
!@var ALTIME ARM TIME STAMP - layer data
      REAL*8 ALTIME
!@var APREC ARM precipitation (mm/hour)
      REAL*8 APREC
!@var ALH Latent Heat Flux (W/m**2)
      REAL*8 ALH
!@var ASH Sensible Heat Flux (W/m**2)
      REAL*8 ASH
!@var AMEANPS ARM Mean Surface Pressure (mb) 
      REAL*8 AMEANPS
!@var ATSAIR ARM surface air temperature (C)
      REAL*8 ATSAIR
!@var ATSKIN surface skin temperature (C)
      REAL*8 ATSKIN
!@var ARHSAIR ARM Surface Air Relative Humidity (%)
      REAL*8 ARHSAIR
!@var ASWINDSPD 10m  wind speed (m/s)
      REAL*8 ASWINDSPD
!@var AQS 2m water vapor mixing ration(kg/kg)
      REAL*8 AQS
!@var AUS 10m u component (m/s)
      REAL*8 AUS
!@var AVS 10m v component (m/s)
      REAL*8 AVS 
!@var AUSRF surface u wind (m/s)
      REAL*8 AUSRF
!@var AVSRF surface v wind (m/s) 
      REAL*8 AVSRF
!@var ALWP  ARM MWR Cloud Liquid Water Path (cm)
      REAL*8 ALWP
!@var ADWDT ARM d(Column H20)/dt  (mm/hr)
      REAL*8 ADWDT
!@var ADWADV ARM Column_H20_Advection_ (mm/hr)
      REAL*8 ADWADV
!@var ATLWUP ARM TOA LW UP (W/m**2)
      REAL*8 ATLWUP
!@var ATSWDN ARM TOA SW DOWN (W/m**2)
      REAL*8 ATSWDN
!@var ATSWIN ARM TOA SW INS (W/m**2)
      REAL*8 ATSWIN
!@var ASRFALBEDO Surface Albedo (initless)
      REAL*8 ASRFALBEDO
!@var ARMDATE Julian day fraction from ARM
      REAL*8 ARMDATE




!@var ARMFAC   factor to take into consideration the difference in size 
!              between the area of the ARM site and the GCM grid box area
!              for the Wind Divergence    
      REAL*8 ARMFAC
!@var ARM_ELEV terrain hgt at ARM site in m
      REAL*8 ARM_ELEV
  
!@var SCM_SURFACE_FLAG 0-use GCM calculated surface fluxes,
!                      1-use SCM prescribed surface fluxes  
!                      2-use SCM surf temps and GCM calculated surface fluxes
      INTEGER SCM_SURFACE_FLAG
!@var SCM_ATURB_FLAG   0-run with dry convection routine
!                      1-run with aturb routine
      INTEGER SCM_ATURB_FLAG
!@var SCM_SURF_ALBEDO_FLAG    0-run with GCM calculated surface albedo
!                             1-run with SCM ARM prescribed surface albedo
      INTEGER SCM_SURF_ALBEDO_FLAG
!@var SCM_RELAX_FORCING_FLAG  0 - run with ARM forcings as given
!                             1 - run with a relaxing over time of the ARM forcings
      INTEGER SCM_RELAX_FORCING_FLAG


!@var NARM #of GCM time steps per ARM time step
      INTEGER NARM
!@var NRINIT #of GCM time steps between reinitializing T,Q 
      INTEGER NRINIT
!@var TAUARM starting TAU of ARM DATA  (should this be starting date and time)
      INTEGER TAUARM
!@var IKT index to arm data interpolated to time steps
      INTEGER IKT
      INTEGER iu_scm_prt,iu_scm_diag,iu_scm_seed
      INTEGER jrandscm    

!@var IFLRESET,NRAMP,IRESET  used for doing updating with a ramp, then saving only
!     time steps after ramp, then backtracking and ramping again before the next saved
!     time steps
!     IFLRESET = flag
!     NRAMP = length of ramp in time steps
!     NRESET = length of ramp + buffer
!     IRESET = index to output buffers
      INTEGER IFLRESET,NRAMP,NRESET,IRESET

      parameter (NRESET=72)

!**** Target Coordinates for SCM
      REAL*8 :: LON_TARG,LAT_TARG
c      INTEGER*4 :: I_TARG,J_TARG   !TWP I=125,J=39  set targets in parameter list
      INTEGER*4 :: NSTEPSCM=0      !Time step counter for SCM
   
      
      INTEGER MCT
      INTEGER NTOTSCM

      parameter (MCT=3000)

      REAL*8 HTA_HR(LM,MCT)    
      REAL*8 VSA_HR(LM,MCT)
      REAL*8 HQA_HR(LM,MCT)
      REAL*8 VQA_HR(LM,MCT)
    
      REAL*4 STMSTEP(MCT)
      REAL*4 STMSTEPL(MCT)
      REAL*4 AMPS(MCT)
      REAL*4 ATSK(MCT)
      REAL*4 ATSA(MCT)
      REAL*4 ARHHR(MCT)
      REAL*4 THR(LM,MCT)
      REAL*4 QHR(LM,MCT)
      REAL*4 APRCHR(MCT)
      REAL*4 ALHHR(MCT)
      REAL*4 ASHHR(MCT)
      REAL*4 AWSHR(MCT)
      REAL*4 AQSHR(MCT)
      REAL*4 AUSHR(MCT)
      REAL*4 AVSHR(MCT)
      REAL*4 UHR(LM,MCT)
      REAL*4 VHR(LM,MCT)
      REAL*4 OMGHR(LM,MCT)
      REAL*4 WDHR(LM,MCT)
      REAL*4 ACLDHR(LM,MCT)
      REAL*4 ALWPHR(MCT)
      REAL*4 ADWDTHR(MCT)
      REAL*4 ADWADVHR(MCT)
      REAL*4 ATLWUPHR(MCT) 
      REAL*4 ATSWDNHR(MCT)
      REAL*4 ATSWINHR(MCT)
      REAL*4 ASRFALBHR(MCT)

      REAL*8 SCM_SAVE_T(LM),SCM_SAVE_Q(LM),SCM_DEL_T(LM),
     &       SCM_DEL_Q(LM)

c
c     add buffers for running ramp/reset to store cloud variables
c
      real*8 CBTTOLD(LM,0:MCT),CBQTOLD(LM,0:MCT),CBWM(LM,0:MCT),
     *       CBPTOLD(0:MCT),CBSVLHX(LM,0:MCT),
     *       CBRHSAV(LM,0:MCT),CBCLDSAV(LM,0:MCT),
     *       CBCLDSAV1(LM,0:MCT)

      end module SCMCOM
c
c    
      subroutine ALLOC_SCM_COM()
   
      USE SCMCOM, only : SCM_SURFACE_FLAG,SCM_ATURB_FLAG,
     &            SCM_SURF_ALBEDO_FLAG,SCM_RELAX_FORCING_FLAG,
     &            ARMFAC,ARM_ELEV,IFLRESET,IRESET,NRAMP,
     &            iu_scm_prt


!@var SCM_SURFACE FLAG   0-run with GCM calculated surface fluxes
!                        1-run with ARM prescribed surface fluxes
!                        2-RUN WITH ARM srf tmps and GCM calc srf fluxes
      SCM_SURFACE_FLAG = 1     

!@var SCM_ATURB_FLAG     0-run with DRYCNV dry convection routine 
!                        1-run with ATURB turbulence routine    
      SCM_ATURB_FLAG = 1
 
!@var SCM_SURF_ALBEDO_FLAG   0-run with GCM calculated surface albedo
!                            1-run with SCM ARM prescribed surface albedo
      SCM_SURF_ALBEDO_FLAG = 0

!@var SCM_RELAX_FORCING_FLAG   0-run with ARM forcings as given
!                              1-run with relaxation of forcings
      SCM_RELAX_FORCING_FLAG = 0

!@var ARMFAC
!     variable to scale windiv for difference in area of ARM site and
!     GCM grid box    see subroutine FCONV in ATMDYN_SCM.f
c     for now set ARMFAC TO 1.0   to be determined when setting up run.
      ARMFAC=1.0

!@var ARM_ELEV   set terrain hgt of ARM site in m
!                SGP = 320.m
      ARM_ELEV = 320.


!@var  IFLRESET          0-run without ramp/reset/updating
!                        1-run with ramp/reset/updating
      IFLRESET = 0
      NRAMP = 24
      IRESET = 0
      if (IFLRESET.eq.1) write(0,*) 'run with ramps iflreset nramp ',
     &                   iflreset,nramp

      return

      end subroutine ALLOC_SCM_COM

      MODULE GEOM
!@sum  GEOM contains geometric variables and arrays
!@auth M. Kelley

      IMPLICIT NONE
      SAVE

      real*8, dimension(1,1) ::
!@var  lat2d latitude of mid point of primary grid box (radians)
!@var  lat2d_dg latitude of mid point of primary grid box (degrees)
     &     lat2d,lat2d_dg
!@var  lon2d longitude of mid point of primary grid box (radians)
!@var  lon2d_dg longitude of mid point of primary grid box (degrees)
     &     ,lon2d,lon2d_dg
     &     ,lon_dg,lat_dg
     &     ,sinlat2d, coslat2d

!@var  axyp,byaxyp area of grid box (+inverse) (m^2)
!@+    Note these should not play any role in a single-column
!@+    model. For clarity, references to them will be replaced
!@+    by usage of axyp_nominal in parts of the model which
!@+    apply a scale-aware parameterization approach.
      real*8, dimension(1,1) :: axyp, byaxyp, axyp_nominal

      integer, dimension(1) :: imaxj

      CONTAINS

      SUBROUTINE GEOM_ATM
      USE CONSTANT, only : PI,TWOPI,radian
      use Dictionary_mod, only : get_param,sync_param
      use scmcom, only : lon_targ, lat_targ
      implicit none

      ! mandatory rundeck parameters: lon and lat of target point
      call get_param('lon_targ',lon_targ)
      call get_param('lat_targ',lat_targ)

      if(abs(lon_targ).gt.180d0 .or. abs(lat_targ).gt.90d0)
     &     call stop_model(
     &       'geom_atm: invalid lon_targ,lat_targ in rundeck',255)

      lon2d_dg(1,1) = lon_targ
      lat2d_dg(1,1) = lat_targ

      axyp_nominal(1,1) = 1.
      call sync_param('nominal_area',axyp_nominal(1,1))

      axyp(1,1) = axyp_nominal(1,1)

      byaxyp(1,1) = 1d0/axyp(1,1)

      lon2d(1,1) = lon2d_dg(1,1)*radian
      lat2d(1,1) = lat2d_dg(1,1)*radian

      sinlat2d(1,1) = sin(lat2d(1,1))
      coslat2d(1,1) = cos(lat2d(1,1))
      lon2d(1,1) = lon2d(1,1) + pi ! IDL has a value of zero
      if(lon2d(1,1) .lt. 0.) lon2d(1,1)= lon2d(1,1) + twopi

      imaxj = 1

      lon_dg = lon2d_dg
      lat_dg = lat2d_dg

      return
      END SUBROUTINE GEOM_ATM

      subroutine lonlat_to_ij(ll,ij)
      implicit none
      real*8, intent(in) :: ll(2)
      integer, intent(out) :: ij(2)
      ij = (/ 1, 1 /)
      return
      end subroutine lonlat_to_ij

      END MODULE GEOM
