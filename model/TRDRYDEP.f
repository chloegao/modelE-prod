#include "rundeck_opts.h"

#ifdef CUBED_SPHERE
#define BIN_OLSON
#endif

!#undef BIN_OLSON

      MODULE tracers_DRYDEP

!@sum  tracers_DRYDEP tracer dry deposition from Harvard CTM.
!@+    Current version only calculates the "bulk surface resistance
!@+    to deposition" component of the deposition velocity.
!@auth D.J. Jacob and Y.H. Wang, modularized by G.M. Gardner, 
!@+    adapted for GISS GCM by D. Koch, modelEified by G. Faluvegi
!@+    Harvard version 3.1: 12/17/97)  
C*********************************************************************
C  Literature cited in drydep.f routines: 
C     Baldocchi, D.D., B.B. Hicks, and P. Camara, A canopy stomatal
C       resistance model for gaseous deposition to vegetated surfaces,
C       Atmos. Environ. 21, 91-101, 1987.
C     Guenther, A., and 15 others, A global model of natural volatile
C       organic compound emissions, J. Geophys. Res., 100, 8873-8892,
C       1995.  
C     Brutsaert, W., Evaporation into the Atmosphere, Reidel, 1982.
C     Dwight, H.B., Tables of integrals and other mathematical data,
C       MacMillan, 1957.
C     Hicks, B.B., and P.S. Liss, Transfer of SO2 and other reactive
C       gases across the air-sea interface, Tellus, 28, 348-354, 1976.
C     Jacob, D.J., and S.C. Wofsy, Budgets of reactive nitrogen,
C       hydrocarbons, and ozone over the Amazon forest during the wet
C       season, J.  Geophys. Res., 95, 16737-16754, 1990.
C     Jacob, D.J., and 9 others, Deposition of ozone to tundra,
C       J. Geophys. Res., 97, 16473-16479, 1992.
C     Levine, I.N., Physical Chemistry, 3rd ed., McGraw-Hill, New York,
C        1988.  
C     Munger, J.W., and 8 others, Atmospheric deposition of reactive
C       nitrogen oxides and ozone in a temperate deciduous forest and a
C       sub-arctic woodland, J. Geophys. Res., in press, 1996.
C     Walcek, C.J., R.A. Brost, J.S. Chang, and M.L. Wesely, SO2,
C        sulfate, and HNO3 deposition velocities computed using
C        regional landuse and meteorological data, Atmos. Environ.,
C        20, 949-964, 1986.   
C     Wang, Y.H., paper in preparation, 1996.   
C     Wesely, M.L, Improved parameterizations for surface resistance to
C       gaseous dry deposition in regional-scale numerical models,
C       Environmental Protection Agency Report EPA/600/3-88/025,   
C       Research Triangle Park (NC), 1988.   
C     Wesely, M.L., same title, Atmos. Environ., 23, 1293-1304, 1989.
C*********************************************************************

      IMPLICIT NONE
      SAVE

!@param NPOLY number of dry dep polynomial coefficients
!@param NENT number of ENT vegetation types
!@param NTYPE number of surface types
!@param NDRYDEP number of dry deposition surface types 
!@var IJUSE fraction of gridbox area occupied by land type
!@+          element LDT
!@var IDEP maps from Ent to dry deposition PFTs
!@var IRI read in ___ resistance
!@var IRLU Cuticular resistances per unit area of leaf
!@var IRAC read in ___ resistance
!@var IRGSS read in ___ resistance
!@var IRGSO read in ___ resistance
!@var IRCLS read in ___ resistance
!@var IRCLO read in ___ resistance
!@var IVSMAX maximum deposition velocity for aerosol
!@var XYLAI Leaf Area Index of land type element LDT
!@var DRYCOEFF polynomial fittings coeffcients  
      INTEGER, PARAMETER :: NENT    = 18 
      INTEGER, PARAMETER :: NPOLY   = 20, NTYPE = 21
      INTEGER, PARAMETER :: NDRYDEP = 11
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)  :: XYLAI,IJUSE 
      REAL*8, DIMENSION(NPOLY)               :: DRYCOEFF
      INTEGER, DIMENSION(NENT)               :: IDEP
      REAL*8, DIMENSION(NDRYDEP)             :: IRI,IRLU,IRAC,IRGSS,
     &                                      IRGSO,IRCLS,IRCLO,IVSMAX
      
      END MODULE tracers_DRYDEP



      subroutine alloc_trdrydep(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth G.Faluvegi
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use tracers_DRYDEP, only: NENT,XYLAI,IJUSE

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: I_0,I_1,J_0,J_1
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call getDomainBounds(grid)
      I_0=grid%I_STRT
      I_1=grid%I_STOP
      J_0=grid%J_STRT
      J_1=grid%J_STOP
 
      allocate(   XYLAI(I_0:I_1,J_0:J_1,NENT) )
      allocate(   IJUSE(I_0:I_1,J_0:J_1,NENT) )

      end subroutine alloc_trdrydep



#ifdef TRACERS_DRYDEP
      SUBROUTINE get_dep_vel(I,J,ITYPE,OBK,ZHH,USTARR,TEMPK,DEP_VEL,
     & stomatal_dep_vel,trnmm
#ifdef TRACERS_TOMAS
     & ,gs_vel,psurf
#endif
     &     )
!@sum  get_dep_vel computes the Bulk surface reistance to
!@+    tracer dry deposition using a resistance-in-series model
!@+    from a portion of the Harvard CTM dry deposition routine.
!@+    The deposition velocity is the reciprocal of the
!@+    Bulk surface reistance...
!@auth D.J. Jacob and Y.H. Wang, modularized by G.M. Gardner, 
!@+    adapted for GISS GCM by D. Koch modelEified by G. Faluvegi
!@+    Harvard version 3.1: 12/17/97)  
C      uses functions: BIOFIT,DIFFG
c
C**** GLOBAL parameters and variables:  
C
      USE GEOM,       only : imaxj
      USE CONSTANT,   only : tf,pi,grav,gasC     
      USE RAD_COM,     only: COSZ1,cfrac,srdn
      use OldTracer_mod, only: tr_wd_TYPE, nPart, trname
      use OldTracer_mod, only: dodrydep, F0_glob=>F0, HSTAR_glob=>HSTAR
      USE TRACER_COM, only : NTM
#ifdef TRACERS_SPECIAL_Shindell
     & , n_NOx
      USE TRCHEM_Shindell_COM, only : pNOx
#endif
      USE tracers_DRYDEP, only: NPOLY,XYLAI,
     & DRYCOEFF,IJUSE,NTYPE,IDEP,IRI,IRLU,IRAC,IRGSS,IRGSO,
     & IRCLS,IRCLO,IVSMAX,NENT
#ifdef TRACERS_TOMAS 
      USE TRACER_COM, only : NBS, NBINS, n_ASO4
#endif
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments
C
!@param XMWH2O molecular weight of water in KG/mole
      REAL*8, PARAMETER ::  XMWH2O = 18.d-3
!@var N the current tracer number
!@var LDT dummy loop variable for Ent land type
!@var k,IW dummy loop variables
!@var RLU cuticular resistance for the bulk canopy
!@var RAC transfer resistance
!@var RGSS ______ resistance
!@var RGSO ______ resistance
!@var RCLS low canopy resistance
!@var RCLO ______ resistance
!@var RSURFACE Bulk surface resistance for species K landtype LDT
!@var RI internal resistance (minimum stomatal resistance for
!@+   water vapor, per unit area of leaf)
!@var RSTOMATAL Stomatal portion of RSURFACE for ozone, landtype LDT
!@var RAD0 downward solar radiation flux at surface (w/m2)
!@var RT correction term for bulk surface resistance for gases?
!@var RIX ______ resistance
!@var GFACT,GFACI,RGSX,RCLX,CZH ?
!@var RIXX corrected RIX
!@var RLUXX corrected RLU
!@var RDC aerodynamic resistance to elements in lower part of
!@+   the canopy or structure
!@var DTMP1,DTMP2,DTMP3,DTMP4 temp arrays for recipricol resistances
!@var VDS temp deposition velocity
!@var DUMMY1,DUMMY2,DUMMY3,DUMMY4 dummy temp variables
!@var TEMPK Surface air temperature (Kelvin)
!@var TEMPC Surface air temperature (oC)  
!@var byTEMPC 1/TEMPC
!@var I,J GCM grid box horizontal position
!@var ITYPE GCM surface type 1=ocean; 2=ocean ice; 3=land ice; 4=land
!@var TOTA total landuse area, excluding water and ice dep. types.
!@var problem_point logical that is true if ITYPE=4, but there are no
!@+   non-water, non-ice portions, according to Harvard dataset.
!@var OBK Monin-Obukhov length (m) (now is the PBL lmonin variable)
!@var ZHH boundary layer depth (m) ("mixing depth") (PBL dbl variable)
!@var USTARR Friction velocity (m s-1) (PBL ustar variable)
!@var tr_mm_temp temporary variable to hold trnmm(k), etc.
!@var SUNCOS Cosine of solar zenith angle
!@var VD deposition velocity temp array   (s m-1)
!@var SVD stomatal deposition velocity temp array   (s m-1)
      REAL*8,  DIMENSION(ntm,NTYPE) :: RSURFACE
      REAL*8,  DIMENSION(NTYPE) :: RSTOMATAL
      REAL*8,  DIMENSION(ntm)   :: TOTA,VD,HSTAR,F0,trnmm
      REAL*8,  DIMENSION(NTYPE) :: RI,RLU,RAC,RGSS,RGSO,RCLS,RCLO
      REAL*8 :: RT,RAD0,RIX,GFACT,GFACI,RDC,RIXX,RLUXX,RGSX,RCLX,VDS,
     &       DTMP1,DTMP2,DTMP3,DTMP4,CZH,DUMMY1,DUMMY2,DUMMY3,DUMMY4,
     &       TEMPC,byTEMPC,BIOFIT,DIFFG,tr_mm_temp,SUNCOS,SVD
      REAL*8, INTENT(IN) :: OBK,ZHH,USTARR,TEMPK
!@var dep_vel the deposition velocity = 1/bulk sfc. res. (m/s)
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: dep_vel
      REAL*8, INTENT(OUT) :: stomatal_dep_vel
      INTEGER :: k,n,LDT,II,IW
      INTEGER, INTENT(IN) :: I,J,ITYPE  
      LOGICAL :: problem_point
#ifdef TRACERS_TOMAS 
      real*8  XNU
      integer binnum,Ni    !@var binnum size bin # that corresponds to current tracer
      real*8 rb(NBINS)  !@var rb quasilaminar sublayer resistance (s m-1)
      real*8 vs(NBINS)  !@var vs gravitational settling velocity (m s-1)
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: gs_vel
      REAL*8,INTENT(IN) :: psurf ! surface pressure
!      REAL*8, INTENT(IN), DIMENSION(NTM) :: TM
      real*8 Dp(nbins),density(nbins) !particle diameter (m)
      real*8 Dk          !@var Dk particle diffusivity (m2/s)
      real*8 mu          !@var mu air viscosity (kg/m s)
      real*8 Sc,St       !@var Sc/St  particle Schmidt and Stokes numbers
      real*8 kB           !@var kB Boltzmann constant (J/K)
      parameter (kB=1.38d-23)
      real mfp  !@var mfp mean free path of air molecule (m)
      real slipc(nbins)
#endif

C Use cosine of the solar zenith angle from the radiation code,
C ...which seems to have a minumum of 0

      SUNCOS = COSZ1(I,J)

C* Initialize VD and RSURFACE and reciprocal: 
      stomatal_dep_vel = 0.d0
      SVD = 0.d0
      RSTOMATAL(1:NTYPE) = 0.d0
      DO K = 1,ntm
        if(dodrydep(K))then
          RSURFACE(K,1:NTYPE) = 0.d0
          VD(K)               = 0.d0
          dep_vel(K)          = 0.d0
#ifdef TRACERS_TOMAS
          gs_vel(K)           = 0.d0
#endif
       end if
      END DO    

C** TEMPK and TEMPC are surface air temperatures in K and in C  
      TEMPC = TEMPK-tf ! TEMPK was BLDATA(I,J,2)
      byTEMPC = 1.D0/TEMPC    
      RAD0 = srdn(I,J)*suncos
              
#ifdef TRACERS_TOMAS

      XNU = 0.0000151*(TEMPK/tf)**1.77
      mu=2.5277e-7*tempk**0.75302
                                
!calculate particle physical and depositional properties      
 
      call dep_getdp(i,j,1,Dp,density)
       
      do k=1,nbins         
         Dk=kB*tempk/(3.0*pi*mu*Dp(k)) !S&P Table 12.1
         Sc=XNU/Dk
         mfp=2.0*mu/(psurf*100.*sqrt(8.0*0.0289/(pi*gasC*TEMPK))) !mfp of air
         Slipc(k)=1.d0+2.d0*mfp/Dp(k)*(1.257+0.4*exp(-.55d0*Dp(k)/mfp)) !S&P eqn 19.20
         vs(k)=density(k)*(Dp(k)**2)*grav/18.d0/mu*Slipc(k) !S&P eqn 19.19
         St=vs(k)*USTARR**2/grav/XNU
         rb(k)=1.d0/(USTARR     !S&P eqn 19.18
     &        *(Sc**(-2.d0/3.d0)+10.d0**(-3.d0/St)))
      enddo
#endif

C* Compute bulk surface resistance for gases.
C*   
C* Adjust external surface resistances for temperature;
C* from Wesely [1989], expression given in text on p. 1296.
C* Note: the sign of 4.0 was fixed after consulting w/ Harvard.

      RT = 1.d3*EXP(-TEMPC-4.0)

C  Check for species-dependant corrections to drydep parameters:      
      DO K = 1,ntm
        if(dodrydep(K)) then
          HSTAR(K)=HSTAR_glob(K)
          F0(K)=F0_glob(K)
#ifdef TRACERS_SPECIAL_Shindell
          ! For NOx, sum deposition for NO and NO2
          ! HSTAR: NO2 = 0.01, NO = 0.002, F0: NO2 = 0.1, NO = 0.0
          if(trname(K) == 'NOx')then
             HSTAR(K)=pNOx(1,i,j)*0.01d0+(1.-pNOx(1,i,j))*2.d-3
             F0(K)=pNOx(1,i,j)*1.d-1
          endif 
#endif
        end if
      END DO

      SELECT CASE(ITYPE)
      CASE(4)     ! LAND*************************************

C    Get surface resistances - loop over land types LDT   
C**********************************************************************
C* The land types within each grid square are defined using the Ent 
C* land-types. Each of the Ent land types is assigned a
C* corresponding "deposition land type" with characteristic values of
C* surface resistance components. There are 18 Ent land-types but
C* only 11 deposition land-types (i.e., some of the Ent land types
C* share the same deposition characteristics). Surface resistance 
C* components for the "deposition land types" are from Wesely [1989]
C* except for tropical forests [Jacob and Wofsy, 1990] 
C* and for tundra [Jacob et al., 1992]. All surface resistance
C* components are normalized to a leaf area index of unity.
C*  
C* Previously, 74 Olson land-types were used instead of Ent's 18. Olson 
C* land-types, deposition land types, and surface resistance components 
C* were previously read from file DRYTBL='drydep.table'; check that file
C* for further details.
C**********************************************************************
C* first check for any points that are ITYPE=4 but all ice/water:
        problem_point=.true.
        DO LDT=1, NENT
          IF (IJUSE(I,J,LDT).gt.0.d0) THEN
               II=IDEP(LDT)
               if(II /= 1 .and. II /= 11) problem_point=.false.
          ENDIF
        END DO

        DO LDT = 1,NENT
          IF (IJUSE(I,J,LDT).gt.0.d0) THEN
            II = IDEP(LDT)
            ! exclude ice and water (except for problem points) :
            IF((II /= 1 .and. II /= 11) .or. problem_point) THEN
   
C** Here, we should probably put some provisions that if the GCM land
C** grid box is mostly snow-covered, set II=1 <######################
C**       
C* Read the internal resistance RI (minimum stomatal resistance for
C* water vapor, per unit area of leaf) from the IRI array; a '9999'
C* value means no deposition to stomata so we impose a very large   
C* value for RI.  

            RI(LDT) = REAL(IRI(II))   
            IF (RI(LDT) >=  9999.) RI(LDT)= 1.d12

C** Cuticular resistances IRLU are per unit area of leaf; divide them 
C** by the LAI to get a cuticular resistance for the bulk canopy.  If IRLU is '9999' it
C** means there are no cuticular surfaces on which to deposit so we  
C** impose a very large value for RLU.  

            IF (IRLU(II) >=  9999 .OR. XYLAI(I,J,LDT) <= 0.) THEN
              RLU(LDT)  = 1.d6
            ELSE
              RLU(LDT)= REAL(IRLU(II))/XYLAI(I,J,LDT) + RT
            ENDIF

C** The following are the remaining resistances for the Wesely   
C** resistance-in-series model for a surface canopy   
C** (see Atmos. Environ. paper, Fig.1).

            RAC(LDT)  = MAX(DBLE(IRAC(II)), 1.d0) 
            IF (RAC(LDT)   >=  9999.) RAC(LDT)  = 1.d12 
            RGSS(LDT) = MAX(DBLE(IRGSS(II)) + RT ,1.d0)
            IF (RGSS(LDT)  >=  9999.) RGSS(LDT) = 1.d12
            RGSO(LDT) = MAX(DBLE(IRGSO(II)) + RT ,1.d0)
            IF (RGSO(LDT)  >=  9999.) RGSO(LDT) = 1.d12
            RCLS(LDT) = DBLE(IRCLS(II)) + RT
            IF (RCLS(LDT)  >=  9999.) RCLS(LDT) = 1.d12
            RCLO(LDT) = DBLE(IRCLO(II)) + RT
            IF (RCLO(LDT)  >=  9999.) RCLO(LDT) = 1.d12
 
C**********************************************************************
C*  Adjust stomatal resistances for insolation and temperature:
C*  Temperature adjustment is from Wesely [1989], equation (3).
C*
C*  Light adjustment by the function BIOFIT described by Wang [1996].
C*  It combines
C*  - Local dependence of stomal resistance on the intensity I of light
C*    impinging the leaf; this is expressed as a mutliplicative
C*    factor I/(I+b) to the stomatal resistance where b = 50 W m-2
C*    (equation (7) of Baldocchi et al. [1987])
C*  - radiative transfer of direct and diffuse radiation in the
C*    canopy using equations (12)-(16) from Guenther et al. [1995]
C*  - separate accounting of sunlit and shaded leaves using
C*    equation (12) of Guenther et al. [1995]
C*  - partitioning of the radiation @ the top of the canopy into direct
C*    and diffuse components using a parameterization to results from
C*    an atmospheric radiative transfer model [Wang, 1996]
C*  The dependent variables of the function BIOFIT are the leaf area
C*  index (XYLAI), the cosine of zenith angl (SUNCOS) and the fractional
C*  cloud cover (CFRAC).  The factor GFACI integrates the light
C*  dependence over the canopy depth; sp even though RI is input per 
C*  unit area of leaf it need not be scaled by LAI to yield a bulk
C*  canopy value because that's already done in the GFACI formulation.
C**********************************************************************

            RIX = RI(LDT)  
            IF (RIX  <  9999.) THEN
              IF (TEMPC > 0. .AND. TEMPC < 40.) THEN
                GFACT = 400.d0*byTEMPC/(40.d0-TEMPC)  
              ELSE
                GFACT = 100.d0
              END IF
              IF (RAD0 > 0. .AND. XYLAI(I,J,LDT) > 0.) THEN  
                GFACI=1./BIOFIT
     *          (DRYCOEFF,XYLAI(I,J,LDT),SUNCOS,CFRAC(I,J))
              ELSE
                GFACI = 100.d0
              ENDIF   
              RIX = RIX*GFACT*GFACI   
            END IF

C*    Compute aerodynamic resistance to lower elements in lower part
C*    of the canopy or structure, assuming level terrain -
C*    equation (5) of Wesely [1989].

            RDC = 100.d0*(1.+1000.d0/(RAD0 + 10.d0))

C*    Loop over species; species-dependent corrections to resistances
C*    are from equations (6)-(9) of Wesely [1989].

            DO K = 1,ntm
             if(dodrydep(K))then
              tr_mm_temp = trnmm(k)*1.d-3
#ifdef TRACERS_SPECIAL_Shindell
              ! For NOx, use NO2 mol. wt. to get collision diameter:
              if(trname(K) == 'NOx') tr_mm_temp = 4.4d-2
#endif
C**           For non-aerosols:
              IF(tr_wd_TYPE(K) /= nPART) THEN
                RIXX=RIX*DIFFG(TEMPK,XMWH2O)/DIFFG(TEMPK,tr_mm_temp)
     &          + 1.d0/(HSTAR(K)/3.d3+100.d0*F0(K))
                IF(RLU(LDT) < 9999.) THEN
                  RLUXX = RLU(LDT)/(1.d-5*HSTAR(K)+F0(K))
                ELSE
                  RLUXX = 1.d12
                END IF

C* To prevent virtually zero resistance to species with huge HSTAR,
C* such as HNO3, a minimum value of RLUXX needs to be set. The
C* rationality of the existence of such a minimum is demonstrated by
C* the observed relationship between Vd(NOy-NOx) and Ustar in Munger
C* et al.[1996]; Vd(HNO3) never exceeds 2 cm s-1 in observations. The
C* corresponding minimum resistance is 50 s m-1.  This correction
C* was introduced by J.Y. Liang on 7/9/95.  

                IF(RLUXX  <  50.) RLUXX = 50.d0
                RGSX = 1.d0/(1.d-5*HSTAR(K)/RGSS(LDT)+F0(K)/RGSO(LDT))
                RCLX = 1.d0/(1.d-5*HSTAR(K)/RCLS(LDT)+F0(K)/RCLO(LDT))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Note from Greg & Gavin: it might be necessary to 
C limit some of these other resistances too:
                IF(RGSX <  50.) RGSX= 50.d0
C               IF(RCLX <  50.) RCLX= 50.d0
C               IF(RIXX <  50.) RIXX= 50.d0
C               IF(RDC <  50.) RDC= 50.d0
C               IF(RAC(LDT) <  50.) RAC(LDT)= 50.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C** Get bulk surface resistance of the canopy, from the network
C** of resistances in parallel and in series (Fig. 1 of Wesely [1989])
                DTMP1=1.d0/RIXX
                DTMP2=1.d0/RLUXX   
                DTMP3=1.d0/(RAC(LDT)+RGSX)
                DTMP4=1.d0/(RDC+RCLX)
                RSURFACE(K,LDT) = 1.d0/(DTMP1+DTMP2+DTMP3+DTMP4)
                if(trname(K)=='Ox')
     &          RSTOMATAL(LDT)=RIXX
              END IF                                 ! gases (above)
            
              IF (tr_wd_TYPE(K) == nPART) THEN ! aerosols (below)
C**             Get surface deposition velocity for aerosols if needed
C**             Equations (15)-(17) of Walcek et al. [1986]            
                VDS = 0.002d0*USTARR
                IF(OBK < 0.)VDS = VDS*(1.d0+(-300.d0/OBK)**0.6667)
                IF(OBK == 0.) call stop_model('OBK=0 in TRDRYDEP',255)
                CZH  = ZHH/OBK
                IF(CZH < -30.)VDS=0.0009d0*USTARR*(-CZH)**0.6667

C* Set VDS to be less than VDSMAX (entry in input file divided by 1.D4)
C* VDSMAX is taken from Table 2 of Walcek et al. [1986].
C* Invert to get corresponding R

                RSURFACE(K,LDT)=1.d0/MIN(VDS,1.d-4*REAL(IVSMAX(II)))
              END IF ! aerosols

C* Set max, min for bulk surface (and stomatal portion) resistances:

              RSURFACE(K,LDT)=MAX(1.d0, MIN(RSURFACE(K,LDT), 9999.d0))
#ifdef TRACERS_TOMAS
              if(k.ge.n_ASO4(1))THEN 
! for size-resolved aerosol model
! use this formula for size-resolved aerosols
! Seinfeld & Pandis, eqn 19.7
                binnum=mod(K-n_ASO4(1)+1,NBINS)
                if (binnum.eq.0) binnum=NBINS
                VDS=1.d0/rb(binnum)
                gs_vel(k)=vs(binnum) !grav. settling velocity for TOMAS model
                RSURFACE(K,LDT)=MAX(1.D0, rb(binnum)) 
              endif
#endif
             end if! dodrydep
            END DO ! K loop   
            ! 1.d12 max in next line is from the ocean section below
            RSTOMATAL(LDT)=MAX(1.d0, min(1.d12, RSTOMATAL(LDT)))
            END IF ! end ice and water exclusion          
          END IF   ! IJUSE ne 0.d0
        END DO     ! LDT loop

C* Loop through the different landuse types present in the grid square.
C* IJUSE is the fraction of the grid square occupied by surface LDT.  
C* Add the contribution of surface type LDT to the bulk surface resistance:

        TOTA = 0.d0 ! total non-water, non-ice area
        DO LDT=1, NENT
          IF (IJUSE(I,J,LDT).gt.0.d0) THEN
            II = IDEP(LDT)
            ! exclude ice and water (except for problem points) :
            IF((II /= 1 .and. II /= 11) .or. problem_point) THEN
             DO K = 1,ntm
              if(dodrydep(K))then
                VD(K) = VD(K) +
     &          IJUSE(I,J,LDT)/RSURFACE(K,LDT)
                TOTA(K)=TOTA(K)+IJUSE(I,J,LDT)
                if(trname(K)=='Ox')SVD = SVD +
     &          IJUSE(I,J,LDT)/RSTOMATAL(LDT)
              end if ! dodrydep
             END DO  ! K
            END IF   ! end ice and water exclusion
          END IF     ! IJUSE gt 0.d0
        END DO       ! LDT  


C* Calculate the deposition velocity, to be returned:

        DO K = 1,ntm
          if(dodrydep(K))then
            dep_vel(K) = VD(K)/TOTA(K)
            if(trname(K)=='Ox')stomatal_dep_vel=SVD/TOTA(K) 
          end if
        END DO       
         
      CASE(1:3) ! OCEAN, OCEAN ICE, LANDICE *******************
        
C       Please see comments in land case above:     
        II = 1                  ! ice
        IF(ITYPE == 1) II  = 11 ! water
        LDT=1 ! just so we do not have to use all new variables
        RI(LDT)   = 1.d12 ! No stomatal deposition
        RLU(LDT)  = 1.d6  ! No cuticular surface deposition 
        RAC(LDT)  = 1.d0
        RGSS(LDT) = MAX(REAL(IRGSS(II)) + RT ,1.d0)
        IF (RGSS(LDT)  >=  9999.) RGSS(LDT) = 1.d12
        RGSO(LDT) = MAX(REAL(IRGSO(II)) + RT ,1.d0)
        IF (RGSO(LDT)  >=  9999.) RGSO(LDT) = 1.d12
        RCLS(LDT) = 1.d12
        RCLO(LDT) = REAL(IRCLO(II)) + RT
        IF (RCLO(LDT)  >=  9999.) RCLO(LDT) = 1.d12
        RIX = RI(LDT)  
        IF (RIX  <  9999.) THEN
          IF (TEMPC > 0. .AND. TEMPC < 40.) THEN
            GFACT = 400.d0*byTEMPC/(40.d0-TEMPC)
          ELSE
            GFACT = 100.d0
          END IF
          IF (RAD0 > 0. .AND. XYLAI(I,J,LDT) > 0.) THEN
            GFACI=1.d0/BIOFIT
     *      (DRYCOEFF,XYLAI(I,J,LDT),SUNCOS,CFRAC(I,J))
          ELSE
            GFACI = 100.d0
          ENDIF
          RIX = RIX*GFACT*GFACI
        END IF
        RDC = 100.d0*(1.+1000.d0/(RAD0 + 10.d0))
        DO K = 1,ntm
         if(dodrydep(K)) then
          tr_mm_temp = trnmm(k)*1.d-3
#ifdef TRACERS_SPECIAL_Shindell
c         For NOx, use NO2 molecular weight to get collision diameter:
          if(trname(K) == 'NOx') tr_mm_temp = 4.4d-2
#endif
          IF(tr_wd_TYPE(K) /= nPART) THEN  ! NON-AEROSOLS
            RIXX = RIX*DIFFG(TEMPK,XMWH2O)/DIFFG(TEMPK,tr_mm_temp)
     &      + 1.d0/(HSTAR(K)/3.d3+100.d0*F0(K))
            RLUXX = 1.d12
            RGSX = 1.d0 / (1.d-5*HSTAR(K)/RGSS(LDT)+F0(K)/RGSO(LDT))
            RCLX = 1.d0 / (1.d-5*HSTAR(K)/RCLS(LDT)+F0(K)/RCLO(LDT))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Note from Greg & Gavin: it might be necessary to 
C limit some of these other resistances too:
            IF(RGSX <  50.) RGSX= 50.d0
C           IF(RCLX <  50.) RCLX= 50.d0
C           IF(RIXX <  50.) RIXX= 50.d0
C           IF(RDC <  50.) RDC= 50.d0
C           IF(RAC(LDT) <  50.) RAC(LDT)= 50.d0
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            DTMP1=1.d0/RIXX
            DTMP2=1.d0/RLUXX   
            DTMP3=1.d0/(RAC(LDT)+RGSX)
            DTMP4=1.d0/(RDC+RCLX)
            RSURFACE(K,LDT) = 1.d0/(DTMP1+DTMP2+DTMP3+DTMP4)
          ELSE IF (tr_wd_TYPE(K) == nPART) THEN ! AEROSOLS    
            VDS = 0.002d0*USTARR
            IF(OBK < 0.)VDS = VDS*(1.d0+(-300.d0/OBK)**0.6667) 
            IF(OBK == 0.) call stop_model('OBK=0 in TRDRYDEP',255)
            CZH  = ZHH/OBK
            IF(CZH < -30.)VDS=0.0009d0*USTARR*(-CZH)**0.6667
            RSURFACE(K,LDT)=1.d0/MIN(VDS,1.d-4*REAL(IVSMAX(II)))
          END IF ! tracer type
          dep_vel(K)=1.d0/MAX(1.d0, MIN(RSURFACE(K,LDT), 9999.d0))
#ifdef TRACERS_TOMAS
            if(k.ge.n_ASO4(1))THEN 
!     for size-resolved aerosol model
!     use this formula for size-resolved aerosols
!     Seinfeld & Pandis, eqn 19.7
              binnum=mod(K-n_ASO4(1)+1,NBINS)
              if (binnum.eq.0) binnum=NBINS
              VDS=1.d0/rb(binnum)
              gs_vel(k)=vs(binnum) !grav. settling velocity for TOMAS model
              RSURFACE(K,LDT)=rb(binnum)
              dep_vel(K)=1.d0/MAX(1.d0, RSURFACE(K,LDT))
            endif
#endif
         end if  ! do drydep
        END DO   ! tracer loop                  

      CASE DEFAULT
        call stop_model('ITYPE error in TRDRYDEP',255)
      END SELECT

      RETURN
      END SUBROUTINE get_dep_vel 



      REAL*8 FUNCTION BIOFIT(COEFF1,XLAI1,SUNCOS1,CFRAC1)
!@sum BIOFIT Calculates the 'light correction' for the stomatal 
!@+   resistance?
!@auth Y.H. Wang
!@ver ? 

C**** GLOBAL parameters and variables:  

      USE tracers_DRYDEP, only : NPOLY
     
      IMPLICIT NONE

C**** Local parameters and variables and arguments
!@param KK number of terms
!@param ND scaling factor for each variable
!@param X0 maximum for each variable
!@var COEFF1 drydep "coefficients" of fit?
!@var TERM,REALTERM ?
!@var XLAI1 local copy of leaf area index
!@var SUNCOS1 local copy of cosine of solar zenith angle
!@var CFRAC1 local copy of the cloud fraction
!@var K,K1,K2,K3,I loop index
!@var XLOW minimum for each variable
      INTEGER, PARAMETER :: KK=4
      INTEGER, PARAMETER :: ND(KK)=(/0,55,20,11/)
      REAL*8, PARAMETER  :: X0(KK)=(/0.,11.,1.,1./)
      REAL*8, DIMENSION(NPOLY), INTENT(IN) :: COEFF1
      REAL*8, DIMENSION(KK)                :: TERM
      REAL*8, DIMENSION(NPOLY)             :: REALTERM
      REAL*8, INTENT(IN) :: XLAI1,SUNCOS1,CFRAC1
      REAL*8                               :: XLOW
      INTEGER                              :: I, K, K1, K2, K3

      TERM(1)=1.d0 
      TERM(2)=XLAI1 
      TERM(3)=SUNCOS1 
      TERM(4)=CFRAC1
C --- this section relaces SUNPARAM routine --- 
      DO I=2,KK ! NN
        TERM(I)=MIN(TERM(I),X0(I))
        IF (I /= KK) THEN ! hardcode of array position !
          XLOW=X0(I)/REAL(ND(I))
        ELSE
          XLOW= 0.d0
        END IF
        TERM(I)=MAX(TERM(I),XLOW)
        TERM(I)=TERM(I)/X0(I)
      END DO
C ---------------------------------------------
      K=0 
      DO K3=1,KK 
        DO K2=K3,KK 
          DO K1=K2,KK 
            K=K+1 
            REALTERM(K)=TERM(K1)*TERM(K2)*TERM(K3) 
          END DO 
        END DO 
      END DO 
      BIOFIT=0.d0 
      DO K=1,NPOLY 
        BIOFIT=BIOFIT+COEFF1(K)*REALTERM(K) 
      END DO 
      IF (BIOFIT < 0.1) BIOFIT=0.1d0 
      RETURN 
      END FUNCTION BIOFIT
     


      REAL*8 FUNCTION DIFFG(TK,XM)
!@sum DIFFG to calculate tracer molecular diffusivity.
!@auth ? HARVARD CTM
      USE CONSTANT, only : bygasc, gasc, pi, mair, avog
      IMPLICIT NONE
C=====================================================================
C  This function calculates the molecular diffusivity (m2 s-1) in air
C  for a gas X of molecular weight XM (kg) at temperature TK (K) and
C  pressure PRESS (Pa).
C  We specify the molecular weight of air (XMAIR) and the hard-sphere
C  molecular radii of air (RADAIR) & of the diffusing gas (RADX). The
C  molecular radius of air is given in a Table on p. 479 of Levine
C  [1988].  The Table also gives radii for some other molecules. Rather
C  than requesting the user to supply molecular radius we specify here
C  a generic value of 2.E-10 m for all molecules, which is good enough
C  in terms of calculating the diffusivity as long as molecule is not
C  too big.
C======================================================================
!@param XMAIR air molecular weight (KG/mole)
!@param RADX hard-sphere molecular radius of the diffusing gas
!@param RADAIR hard-sphere molecular radius of air
!@param PRESS pressure (kg/s2/m) used to calculate molec. diffusivities
!@var TK passed local temperature in Kelvin
!@var XM molecular weight of diffusing gas (KG/mole)
!@var Z,DIAM ?
!@var FRPATH mean free path of the gas
!@var SPEEN average speed of the gas molecules
!@var AIRDEN local air density
      REAL*8, PARAMETER :: XMAIR = mair * 1.D-3,
     &                     RADX  = 1.5D-10,
     &                     RADAIR= 1.2D-10,
     &                     PRESS=  1.0d5
      REAL*8, INTENT(IN):: TK,XM
      REAL*8 :: Z,DIAM,FRPATH,SPEED,AIRDEN 

C* Calculate air density AIRDEN:
      AIRDEN = PRESS*avog*bygasc/TK 
C* Calculate the mean free path for gas X in air: eq. 8.5 of Seinfeld
C*  [1986]; DIAM is the collision diameter for gas X with air :
      Z = XM/XMAIR
      DIAM = RADX+RADAIR
      FRPATH = 1.d0/(PI*SQRT(1.d0+Z)*AIRDEN*(DIAM**2.))
C* Calculate average speed of gas X; eq. 15.47 of Levine [1988]
      SPEED = SQRT(8.d0*gasc*TK/(PI*XM))
C* Calculate diffusion coefficient of gas X in air; eq. 8.9 of 
C* Seinfeld [1986]  :
      DIFFG = (3.d0*PI*0.03125d0)*(1.d0+Z)*FRPATH*SPEED
      RETURN
      END
         
  
      subroutine rdlai
!@sum Updates the Leaf Area Index (LAI) and land fraction 
!@+   arrays daily, formerly from the formatted 
!@+   input file 'vegtype.global'.

C**** GLOBAL parameters and variables:
      use domain_decomp_atm, only : grid, getDomainBounds
      use ent_com,           only : entcells
      use ent_mod,           only : ent_get_exports
      use ghy_com,           only : fearth
      USE tracers_DRYDEP,    only : IJUSE, XYLAI, NENT
      IMPLICIT NONE

**** Local parameters and variables and arguments
!@var I,J local lat lon index
!@var K, counter dummy loop variable
!@var frac_var, land fractions for each Ent PFT
!@var lai_var, leaf area indices for each Ent PFT

      integer :: J_0, J_1, J_1H, J_0H, I_0, I_1
      real*8, dimension(NENT) :: frac_var, lai_var
      INTEGER :: I, J, K, counter

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &                      J_STRT=J_0 , J_STOP=J_1 )
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

      ! Update LAI and land fraction arrays from Ent
      do j=J_0,J_1
      do i=I_0,I_1
          if(fearth(i,j).gt.0.d0) then
            frac_var(:) = 0.d0
            lai_var(:) = 0.d0
            call ent_get_exports(entcells(i,j),
     &           leaf_area_index_pft=lai_var)
            call ent_get_exports(entcells(i,j),
     &           vegetation_fractions=frac_var)
            do k=1, NENT
                if (frac_var(k).gt.0.d0) then
                  XYLAI(i,j,k) = lai_var(k)/frac_var(k)
                  IJUSE(i,j,k) = frac_var(k)
                end if
            end do ! k
          end if ! fearth>0
      end do ! i
      end do ! j
 
      return 
      end subroutine rdlai 


      subroutine rddrycf 
!@sum Hardcoded resistances for each deposition surface type,
!@+   maximum deposition velocity for aerosol,
!@+   polynomial dry deposition coefficients,
!@+   and the ent->drydep PFT mapping.

C**** GLOBAL parameters and variables:
      USE tracers_DRYDEP, only: IDEP,NPOLY,
     & IRI,IRLU,IRAC,IRGSS,IRGSO,IRCLS,IRCLO,IVSMAX,DRYCOEFF
      IMPLICIT NONE

! Hard coded resistances for each deposition surface type,
! formerly from 'drydep.table'
      IRI = (/9999.d0, 200.d0, 400.d0, 200.d0, 200.d0, 200.d0,
     & 200.d0, 9999.d0, 200.d0, 9999.d0, 9999.d0/)

      IRLU = (/9999.d0, 9000.d0, 9000.d0, 9000.d0, 9000.d0,
     & 1000.d0, 4000.d0, 9999.d0, 9000.d0, 9999.d0, 9999.d0/)

      IRAC = (/0.d0, 2000.d0, 2000.d0, 200.d0, 100.d0, 
     & 2000.d0, 0.d0, 0.d0, 300.d0, 100.d0, 0.d0/)

      IRGSS = (/100.d0, 500.d0, 500.d0, 150.d0, 350.d0, 
     & 200.d0, 340.d0, 1000.d0, 0.d0, 400.d0, 0.d0/)

      IRGSO = (/3500.d0, 200.d0, 200.d0, 150.d0, 200.d0, 200.d0,
     & 340.d0, 400.d0, 1000.d0, 300.d0, 2000.d0/)

      IRCLS = (/9999.d0, 2000.d0, 2000.d0, 2000.d0, 2000.d0,
     & 9999.d0, 9999.d0, 9999.d0, 2500.d0, 9999.d0, 9999.d0/)

      IRCLO = (/1000.d0, 1000.d0, 1000.d0, 1000.d0, 1000.d0,
     & 9999.d0, 9999.d0, 9999.d0, 1000.d0, 9999.d0, 9999.d0/)

! Maximum deposition velocity (Vsmax) for aerosol
      IVSMAX = (/100.d0, 100.d0, 100.d0, 100.d0, 100.d0, 100.d0,
     & 100.d0, 10.d0, 100.d0, 100.d0, 10.d0/)

! Polynomial dry deposition coefficients, formerly from 'drydep.coef'
      DRYCOEFF = (/-3.58d-1, 3.02d0, 3.85d0, -9.78d-2, -3.66d0,
     &  1.2d1, 2.52d-1, -7.8d0, 2.26d-1, 2.74d-1, 1.14d0, -2.19d0, 
     & 2.61d-1, -4.62d0, 6.85d-1, -2.54d-1, 4.37d0, -2.66d-1, -1.59d-1,
     & -2.06d-1/)

! Mapping from ENT to Dry Deposition PFTs
      IDEP = (/6, 6, 3, 3, 2, 2, 2, 2, 7, 5, 5, 5, 5, 5, 4, 4,
     & 8, 8/)

      return
      end subroutine rddrycf

#endif
