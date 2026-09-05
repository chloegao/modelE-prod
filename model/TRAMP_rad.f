#include "rundeck_opts.h"

c -----------------------------------------------------------------

      module AMP_Optics_mod
!@sum Table sizes, index conventions and lookup helpers shared by the AMP
!@+   radiation routines SETAMP, SETAMP_LEV, SETUP_RAD and GET_LW.
!@+
!@+   Everything here is mechanism independent. Per-mode quantities are keyed on
!@+   the mode name rather than on a mode number, and per-species quantities on
!@+   named species indices, so nothing in this file has to be edited when a
!@+   mechanism is added to AERO_CONFIG or its mode list is reordered.
!@auth Susanne Bauer

#ifndef USE_OFFLINE_AEROSOLS
      USE AERO_PARAM, only: RHO_NH4NO3, RHO_H2O
#endif

      IMPLICIT NONE
      SAVE

c-------------------------------------------------------------------------------------------------------
c     Radiation bands.
c-------------------------------------------------------------------------------------------------------
!@param NSW number of shortwave bands
!@param NLW number of longwave (thermal) bands
!@param IBAND_550 shortwave band used for the aerosol optical thickness at 550 nm
      INTEGER, PARAMETER :: NSW = 6
      INTEGER, PARAMETER :: NLW = 33
      INTEGER, PARAMETER :: IBAND_550 = 6

c-------------------------------------------------------------------------------------------------------
c     Axes of the internal-mixture Mie lookup tables read from AMP_MIE_TABLES by
c     SETUP_RAD: real refractive index, imaginary refractive index and effective
c     radius. AMP_EXT, AMP_SCA and AMP_ASY are (NRE,NIM,NSIZ,NSW); AMP_Q55 is
c     (NRE,NIM,NSIZ). SETUP_RAD checks these against the module arrays.
c-------------------------------------------------------------------------------------------------------
!@param NRE  number of tabulated real      refractive indices
!@param NIM  number of tabulated imaginary refractive indices
!@param NSIZ number of tabulated effective radii
      INTEGER, PARAMETER :: NRE  = 15
      INTEGER, PARAMETER :: NIM  = 17
      INTEGER, PARAMETER :: NSIZ = 23

c-------------------------------------------------------------------------------------------------------
c     Axis of the core-shell Mie lookup tables read from AMP_CORESHELL_TABLES:
c     one axis per shell volume fraction (OC, SO4 and H2O), all sharing the same
c     tabulated fractions CS_MIX. AMP_EXT_CS etc. are (NSIZ,NCS,NCS,NCS,NSW) and
c     AMP_Q55_CS is (NSIZ,NCS,NCS,NCS).
c-------------------------------------------------------------------------------------------------------
!@param NCS number of tabulated shell volume fractions
      INTEGER, PARAMETER :: NCS = 26

c-------------------------------------------------------------------------------------------------------
c     The species making up each mode, in the order used by VMass, VolFrac,
c     dry_Vf_LEV and the refractive index table Ri in SETAMP_LEV. Water is last so
c     that the dry species form a leading section, VMass(n,1:NDRY); this is what
c     lets dry_Vf_LEV be normalised over the dry species alone.
c-------------------------------------------------------------------------------------------------------
!@param SPC_SU,SPC_BC,SPC_OC,SPC_DU,SPC_SS,SPC_NO3,SPC_H2O species indices
!@param NDRY number of dry species (all but water)
!@param NSPC total number of species
      INTEGER, PARAMETER :: SPC_SU  = 1   ! sulfate
      INTEGER, PARAMETER :: SPC_BC  = 2   ! black carbon
      INTEGER, PARAMETER :: SPC_OC  = 3   ! organic carbon, summed over all organic tracers
      INTEGER, PARAMETER :: SPC_DU  = 4   ! dust
      INTEGER, PARAMETER :: SPC_SS  = 5   ! sea salt
      INTEGER, PARAMETER :: SPC_NO3 = 6   ! nitrate, together with ammonium
      INTEGER, PARAMETER :: SPC_H2O = 7   ! aerosol water
      INTEGER, PARAMETER :: NDRY = SPC_NO3
      INTEGER, PARAMETER :: NSPC = SPC_H2O

c-------------------------------------------------------------------------------------------------------
c     Densities used to convert the nitrate and water tracer masses to volumes.
c     Online these come from AERO_PARAM so there is one definition; the offline
c     path does not build AERO_PARAM, so it repeats the same two values.
c-------------------------------------------------------------------------------------------------------
!@param DENS_NO3 density of dry NH4NO3 [kg/m^3]
!@param DENS_AQ  density of aerosol water [kg/m^3]
#ifdef USE_OFFLINE_AEROSOLS
      REAL*8, PARAMETER :: DENS_NO3 = 1720.0d0   ! = 1.d3 * RHO_NH4NO3 of AERO_PARAM
      REAL*8, PARAMETER :: DENS_AQ  = 1000.0d0   ! = 1.d3 * RHO_H2O    of AERO_PARAM
#else
      REAL*8, PARAMETER :: DENS_NO3 = 1.0d3 * RHO_NH4NO3
      REAL*8, PARAMETER :: DENS_AQ  = 1.0d3 * RHO_H2O
#endif

c-------------------------------------------------------------------------------------------------------
c     Aerosol compositions understood by GET_LW, and where each one starts in the
c     22-size blocks of the RADPAR thermal tables TRUQEX/TRUQSC. Those blocks run
c     sulfate, sea salt, nitrate, water, organic, so the organic block starts at 88
c     rather than 66: the water block is deliberately skipped. Black carbon and
c     dust are not in TRUQEX at all and are handled from TRSQEX and TRDQEX, hence
c     the zero offsets. This matches the identical offsets in SETSPH of RADIATION.f.
c-------------------------------------------------------------------------------------------------------
!@param LWC_* aerosol composition codes for GET_LW
!@param NLWCLASS number of aerosol compositions
!@param LW_TRUQEX_OFFSET start of each composition's block in TRUQEX/TRUQSC
      INTEGER, PARAMETER :: LWC_SO4  = 1
      INTEGER, PARAMETER :: LWC_SEAS = 2
      INTEGER, PARAMETER :: LWC_NO3  = 3
      INTEGER, PARAMETER :: LWC_OCAR = 4
      INTEGER, PARAMETER :: LWC_BCAR = 5
      INTEGER, PARAMETER :: LWC_DUST = 6
      INTEGER, PARAMETER :: NLWCLASS = 6
      INTEGER, PARAMETER :: NREFU = 22   ! sizes per composition block of TRUQEX
      INTEGER, PARAMETER :: NREFS = 25   ! sizes in TRSQEX (black carbon)
      INTEGER, PARAMETER :: NREFD = 25   ! sizes in TRDQEX (dust)
      INTEGER, PARAMETER, DIMENSION(NLWCLASS) :: LW_TRUQEX_OFFSET =
     +  (/ 0, 22, 44, 88, 0, 0 /)

c-------------------------------------------------------------------------------------------------------
c     SPC_TO_LWC maps the species order above onto the GET_LW composition codes.
c     The two orderings differ, so this mapping must not be bypassed:
c         species index   :  1 SU   2 BC   3 OC   4 DU   5 SS   6 NO3
c         GET_LW code     :  1 SO4  5 BC   4 OC   6 DUST 2 SEAS 3 NO3
c-------------------------------------------------------------------------------------------------------
      INTEGER, PARAMETER, DIMENSION(NDRY) :: SPC_TO_LWC =
     +  (/ LWC_SO4, LWC_BCAR, LWC_OCAR, LWC_DUST, LWC_SEAS, LWC_NO3 /)

c-------------------------------------------------------------------------------------------------------
c     Tabulated axis values. These are written as default-precision literals, as
c     they were when they were DATA statements, so that the tabulated points are
c     bit for bit what they have always been.
c-------------------------------------------------------------------------------------------------------
      REAL*8, PARAMETER, DIMENSION(NSIZ) :: SIZEBINS = (/
     +  0.002, 0.005, 0.01, 0.05, 0.08, 0.1, 0.13, 0.17, 0.2, 0.25, 0.3, 0.4,
     +  0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.5, 2., 3., 5., 10./)
      REAL*8, PARAMETER, DIMENSION(NRE) :: MIE_RE = (/
     +  1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.9/)
      REAL*8, PARAMETER, DIMENSION(NIM) :: MIE_IM = (/
     +  0.0, 0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005,
     +  0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0/)
      REAL*8, PARAMETER, DIMENSION(NCS) :: CS_MIX = (/
     +  0., 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48,
     +  0.52, 0.56, 0.6, 0.64, 0.68, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0/)

c-------------------------------------------------------------------------------------------------------
c     Every mode that can appear in any mechanism, with its effective radius. This
c     is the full set of mode names of MNAME in AERO_CONFIG; the order here does
c     not matter, because modes are looked up by name in MODE_SLOT. The values are
c     the ones the per-mechanism tables used to hold, which were identical for a
c     given mode name in every mechanism that carried it.
c
c     CORE_CLASS/SHELL_CLASS used to sit beside these, assigning each mode a single
c     LW composition. The longwave section of SETAMP now mixes over the composition
c     actually present in each mode, so they are no longer used; the assignment they
c     made is recorded here for reference:
c        AKK,ACC -> SO4;  DD1,DS1,DD2,DS2,DBC,MXX -> DUST;  SSA,SSC,SSS -> SEAS;
c        OCC,OCS,BOC -> OCAR;  BC1,BC2,BC3,BCS -> BCAR;  shell class 0 throughout.
c-------------------------------------------------------------------------------------------------------
!@param NMODES_ALL number of modes that can appear in any mechanism
!@param ALL_MODE_NAME names of all possible modes
!@param ALL_REFF effective radius of each possible mode [um]
      INTEGER, PARAMETER :: NMODES_ALL = 18
      CHARACTER(LEN=3), PARAMETER, DIMENSION(NMODES_ALL) :: ALL_MODE_NAME = (/
     +  'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','SSS',
     +  'OCC','BC1','BC2','BC3','OCS','DBC','BOC','BCS','MXX'/)
      REAL*8, PARAMETER, DIMENSION(NMODES_ALL) :: ALL_REFF = (/
     +  0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00, 2.00D+00 , 0.12D+00 ,
     +  2.D+00   , 1.380D+00, 0.075D+00, 0.050D+00, 0.100D+00, 0.100D+00, 0.075D+00,
     +  0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/)

      contains

c -----------------------------------------------------------------
      INTEGER FUNCTION MODE_SLOT(MODE)
!@sum Position of a mode in ALL_MODE_NAME, or 0 if the name is not listed there.
!@+   Returning 0 rather than a default lets the caller stop the model, so a mode
!@+   added to a mechanism cannot silently pick up another mode's optical properties.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: MODE
      INTEGER :: M
      MODE_SLOT = 0
      DO M = 1, NMODES_ALL
        IF ( MODE .EQ. ALL_MODE_NAME(M) ) MODE_SLOT = M
      ENDDO
      RETURN
      END FUNCTION MODE_SLOT

c -----------------------------------------------------------------
      LOGICAL FUNCTION IS_CORE_SHELL(MODE)
!@sum True for the modes treated as a black carbon core inside a shell of OC, SO4
!@+   and H2O, and so read from the core-shell tables when AMP_RAD_KEY=2 and mixed
!@+   with the Maxwell Garnett rule when AMP_RAD_KEY=3. Every other mode is treated
!@+   as an internal mixture. This is the one place the distinction is made, so the
!@+   two treatments cannot disagree and no mode can fall between them.
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: MODE
      SELECT CASE (MODE)
      CASE ('BC1','BC2','BC3','BOC','BCS')
        IS_CORE_SHELL = .TRUE.
      CASE DEFAULT
        IS_CORE_SHELL = .FALSE.
      END SELECT
      RETURN
      END FUNCTION IS_CORE_SHELL

c -----------------------------------------------------------------
      INTEGER FUNCTION BRACKET_INDEX(X,TABLE)
!@sum Index of the first entry of an ascending TABLE that is not below X, or the
!@+   last index if X lies above the whole table.
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: X, TABLE(:)
      INTEGER :: M
      BRACKET_INDEX = SIZE(TABLE)
      DO M = 1, SIZE(TABLE)
        IF ( X .LE. TABLE(M) ) THEN
          BRACKET_INDEX = M
          RETURN
        ENDIF
      ENDDO
      RETURN
      END FUNCTION BRACKET_INDEX

c -----------------------------------------------------------------
      INTEGER FUNCTION MIE_RE_INDEX(RINDEX)
!@sum Real refractive index axis of the internal-mixture Mie tables.
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: RINDEX
      MIE_RE_INDEX = BRACKET_INDEX( DBLE(REAL(RINDEX)), MIE_RE )
      RETURN
      END FUNCTION MIE_RE_INDEX

c -----------------------------------------------------------------
      INTEGER FUNCTION MIE_IM_INDEX(RINDEX)
!@sum Imaginary refractive index axis of the internal-mixture Mie tables.
      IMPLICIT NONE
      COMPLEX*8, INTENT(IN) :: RINDEX
      MIE_IM_INDEX = BRACKET_INDEX( DBLE(AIMAG(RINDEX)), MIE_IM )
      RETURN
      END FUNCTION MIE_IM_INDEX

c -----------------------------------------------------------------
      INTEGER FUNCTION CS_INDEX(FRAC)
!@sum Shell volume fraction axis of the core-shell Mie tables.
!@+   NOTE: this reproduces exactly the index arithmetic this code has always used,
!@+   which steps back twice from the bracketing entry, so a fraction is read from
!@+   the tabulated point roughly two bins (0.08) below it. That looks like an
!@+   off-by-one, but it is left alone here because correcting it would change
!@+   results, and the intended indexing of AMP_CORESHELL_TABLES has to be
!@+   confirmed against the table itself before it can be changed.
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FRAC
      INTEGER :: M, RAW
      RAW = NCS + 1
      DO M = 1, NCS
        IF ( FRAC .LE. CS_MIX(M) ) THEN
          RAW = M
          EXIT
        ENDIF
      ENDDO
      CS_INDEX = MAX( 1, MIN( NCS, RAW-1 ) - 1 )
      RETURN
      END FUNCTION CS_INDEX

c -----------------------------------------------------------------
      SUBROUTINE SIZE_WEIGHTS(REFF,MD1,MD2,WD1,WD2)
!@sum Linear interpolation of the effective radius axis of the Mie tables.
!@+   Returns the two bracketing size indices MD1 and MD2 and their weights, so
!@+   that any tabulated quantity Q is interpolated as WD1*Q(MD1) + WD2*Q(MD2).
!@+   Below the first tabulated radius the first entry is used unweighted; above
!@+   the last, the last two entries are extrapolated, as before.
!@+   The two weights are formed separately, as b/(a+b) and a/(b+a) rather than
!@+   one as the complement of the other, so that the interpolation is bit for bit
!@+   what it was when this arithmetic was written out at each of its six sites.
      IMPLICIT NONE
      REAL*8,  INTENT(IN)  :: REFF
      INTEGER, INTENT(OUT) :: MD1, MD2
      REAL*8,  INTENT(OUT) :: WD1, WD2
      REAL*8 :: A, B
      MD2 = BRACKET_INDEX( REFF, SIZEBINS )
      IF ( MD2 .GT. 1 ) THEN
        MD1 = MD2 - 1
        B   = SIZEBINS(MD2) - REFF
        A   = REFF - SIZEBINS(MD1)
        WD1 = B / ( A + B )
        WD2 = A / ( B + A )
      ELSE
        MD1 = 1
        WD1 = 0.d0
        WD2 = 1.d0
      ENDIF
      RETURN
      END SUBROUTINE SIZE_WEIGHTS

      end module AMP_Optics_mod

c -----------------------------------------------------------------

      module AMP_Utilities_mod

      contains

      character*2 function aerosolkind(tracerName)
      character(len=*), intent(in) :: tracerName
      aerosolkind = tracerName(7:8)
      end function aerosolkind

      end module AMP_Utilities_mod

c -----------------------------------------------------------------

      SUBROUTINE SETAMP(EXT,SCT,GCB,TAB)
!@sum Calculation of extinction, asymmetry and scattering for AMP Aerosols
!@sum Calculation of absorption in the longwave
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer

      USE AMP_Optics_mod

#ifdef USE_OFFLINE_AEROSOLS
      USE OFFLINE_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV,
     +                       MIX_OC, MIX_SU, MIX_AQ, AMP_RAD_KEY,
     +                       NMODES
c     MODE_NAME is not carried by OFFLINE_AEROSOL, so it is spelled out below.
#else
      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV,
     +                       MIX_OC, MIX_SU, MIX_AQ, AMP_RAD_KEY
      USE AERO_CONFIG, only: NMODES
      USE AERO_SETUP,  only: MODE_NAME
#endif

      USE RESOLUTION,  only: lm
      USE MODEL_COM,   only: itime,itimeI
c     AMP_TAB_SPEC (RADPAR) held one LW spectrum per mode and is no longer used
c     here; it is superseded by AMP_TAB_CLASS below, which resolves composition too.
      USE RADPAR,      only: aesqex,aesqsc,aesqcb,FSTOPX,FTTOPX

      IMPLICIT NONE
      INTEGER, save:: Ifirstrad = 1
      ! Arguments: Optical Parameters dimension(lm,wavelength)
      REAL(8), INTENT(OUT) :: EXT(LM,NSW)     ! Extinction, SW
      REAL(8), INTENT(OUT) :: SCT(LM,NSW)     ! Single Scattering Albedo, SW
      REAL(8), INTENT(OUT) :: GCB(LM,NSW)     ! Asymmetry Factor, SW
      REAL(8), INTENT(OUT) :: TAB(LM,NLW)     ! Thermal absorption Cross section, LW
      REAL(8), DIMENSION(LM,NMODES) :: TTAUSV

      ! Local

      INTEGER l,n,w,s,m,MA,MB,MC,MD1,MD2,NA
      LOGICAL CORE_SHELL
      REAL*8  WD1,WD2,HELP,AMP_TAB(NLW),Vf(NDRY),Reff_mode(NMODES)
      REAL*8  AMPEXT,AMPSCA,AMPASY
c-----------------------------------------------------------------------------------------
c     Longwave absorption is mixed over the composition of each mode, rather than taken
c     from a single composition per mode (see the LW section below for why).
c     AMP_TAB_CLASS holds the LW absorption spectrum of each pure composition at each
c     mode's effective radius, precalculated once; the runtime mixing is a
c     volume-fraction weighted sum over these.
c-----------------------------------------------------------------------------------------
      REAL*8, SAVE :: AMP_TAB_CLASS(NLW,NLWCLASS,NMODES)
#ifdef USE_OFFLINE_AEROSOLS
      CHARACTER(LEN=3), PARAMETER, DIMENSION(NMODES) :: MODE_NAME = (/
     +  'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC',
     +  'OCC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/)
#endif

      EXT(:,:)    = 0.d0
      SCT(:,:)    = 0.d0
      GCB(:,:)    = 0.d0
      TAB(:,:)    = 0.d0
      TTAUSV(:,:) = 0.d0
C Longwave Pre calculate TAB: ---------------------------------------------------------------------------------------------
c     One spectrum per (composition, mode) instead of one per mode. The radius is still
c     Reff_mode(n), so only the composition is resolved here, not the size; the size
c     dependence enters through TTAUSV when the spectra are applied below.
c     GET_LW is called with a shell class of 0: the shell branch inside GET_LW is
c     commented out, and mixing is now handled by the volume weighting below.

      if ( Ifirstrad==1 ) then
        Ifirstrad = 0
c       Effective radius of each mode of this mechanism, looked up by mode name.
        DO n = 1,NMODES
          m = MODE_SLOT( MODE_NAME(n) )
          IF ( m .EQ. 0 ) THEN
            WRITE(*,*) 'SETAMP: no effective radius tabulated for mode ', MODE_NAME(n)
            WRITE(*,*) 'Add it to ALL_MODE_NAME and ALL_REFF in AMP_Optics_mod.'
            CALL stop_model('SETAMP: unknown aerosol mode name',255)
          ENDIF
          Reff_mode(n) = ALL_REFF(m)
        ENDDO
c       Vf is only read by the commented-out shell branch of GET_LW.
        Vf(:) = 0.d0
        DO n  = 1,NMODES
          DO NA = 1,NLWCLASS
            CALL GET_LW(NA,0,Reff_mode(n),AMP_TAB,Vf)
            AMP_TAB_CLASS(:,NA,n) = AMP_TAB(:)
          ENDDO   ! compositions
        ENDDO   ! modes
      endif

      if (itime.ne.itimeI) then
          IF (AMP_RAD_KEY == 1 .or. AMP_RAD_KEY ==3) THEN

c Shortwave: ---------------------------------------------------------------------------------------------

      DO l = 1,lm
      DO n = 1,NMODES

         CALL SIZE_WEIGHTS(Reff_LEV(l,n),MD1,MD2,WD1,WD2)
c---- INTERNAL MIXTURE ---------------------------------------------
         w  = IBAND_550    ! aot at 550
         MA = MIE_RE_INDEX(RindexAMP(l,n,w))
         MB = MIE_IM_INDEX(RindexAMP(l,n,w))
         TTAUSV(l,n) = NUMB_LEV(l,n) * (WD1 * AMP_Q55(MA,MB,MD1) + WD2 * AMP_Q55(MA,MB,MD2))
C----------------------------------------------------------------------
      DO w = 1,NSW  !wavelength
c---- INTERNAL MIXTURE ---------------------------------------------
          MA = MIE_RE_INDEX(RindexAMP(l,n,w))
          MB = MIE_IM_INDEX(RindexAMP(l,n,w))

          AMPEXT = (WD1 * AMP_EXT(MA,MB,MD1,w) + WD2 * AMP_EXT(MA,MB,MD2,w))
          AMPSCA = (WD1 * AMP_SCA(MA,MB,MD1,w) + WD2 * AMP_SCA(MA,MB,MD2,w))
          AMPASY = (WD1 * AMP_ASY(MA,MB,MD1,w) + WD2 * AMP_ASY(MA,MB,MD2,w))
c--------------------------------------------------------------------
          EXT(l,w) = EXT(l,w) + ( AMPEXT * TTAUSV(l,n) * FSTOPX(n))
          HELP     = ((GCB(l,w) * SCT(l,w) ) + (AMPASY * AMPSCA * TTAUSV(l,n)* FSTOPX(n)))
          SCT(l,w) = SCT(l,w) +  (AMPSCA * TTAUSV(l,n) * FSTOPX(n))
          GCB(l,w) = HELP / (SCT(l,w)+ 1.D-10)

          aesqex(l,w,n)= AMPEXT * TTAUSV(l,n)
          aesqsc(l,w,n)= AMPSCA * TTAUSV(l,n)
          aesqcb(l,w,n)= AMPASY * aesqsc(l,w,n)

      ENDDO   ! wave
      ENDDO   ! modes
      ENDDO   ! level

         ENDIF       ! AMP_RAD_KEY=1or3

c --------------------------------------------------------------------------------------------------------
c --------------------------------------------------------------------------------------------------------

         IF (AMP_RAD_KEY == 2) THEN
c Shortwave: ---------------------------------------------------------------------------------------------
c     The black carbon modes are read from the core-shell tables and every other
c     mode from the internal-mixture tables. IS_CORE_SHELL decides, so a mode that
c     is neither black carbon nor listed anywhere here still gets a treatment; when
c     the two groups were spelled out as separate case lists, a mode missing from
c     both (SSS, in mechanisms 4 and 8) was left with TTAUSV=0 and so dropped out
c     of the shortwave altogether.

      DO l = 1,lm
      DO n = 1,NMODES

         CALL SIZE_WEIGHTS(Reff_LEV(l,n),MD1,MD2,WD1,WD2)
         CORE_SHELL = IS_CORE_SHELL(MODE_NAME(n))
         w = IBAND_550    ! aot at 550

       IF (CORE_SHELL) THEN
C------ CORE SHELL -------------------------------------------------
c     The shell composition indices do not depend on wavelength, so they are set
c     here once and reused in the wavelength loop below.
         MA = CS_INDEX(MIX_OC(l,n))
         MB = CS_INDEX(MIX_SU(l,n))
         MC = CS_INDEX(MIX_AQ(l,n))
         TTAUSV(l,n) = NUMB_LEV(l,n) * (WD1 * AMP_Q55_CS(MD1,MA,MB,MC) + WD2 * AMP_Q55_CS(MD2,MA,MB,MC))
       ELSE
c---- INTERNAL MIXTURE ---------------------------------------------
         MA = MIE_RE_INDEX(RindexAMP(l,n,w))
         MB = MIE_IM_INDEX(RindexAMP(l,n,w))
         TTAUSV(l,n) = NUMB_LEV(l,n) * (WD1 * AMP_Q55(MA,MB,MD1) + WD2 * AMP_Q55(MA,MB,MD2))
       ENDIF
C----------------------------------------------------------------------
      DO w = 1,NSW  !wavelength

       IF (CORE_SHELL) THEN
C------ CORE SHELL -------------------------------------------------
          AMPEXT = (WD1 * AMP_EXT_CS(MD1,MA,MB,MC,w) + WD2 * AMP_EXT_CS(MD2,MA,MB,MC,w))
          AMPSCA = (WD1 * AMP_SCA_CS(MD1,MA,MB,MC,w) + WD2 * AMP_SCA_CS(MD2,MA,MB,MC,w))
          AMPASY = (WD1 * AMP_ASY_CS(MD1,MA,MB,MC,w) + WD2 * AMP_ASY_CS(MD2,MA,MB,MC,w))
       ELSE
c---- INTERNAL MIXTURE ---------------------------------------------
          MA = MIE_RE_INDEX(RindexAMP(l,n,w))
          MB = MIE_IM_INDEX(RindexAMP(l,n,w))
          AMPEXT = (WD1 * AMP_EXT(MA,MB,MD1,w) + WD2 * AMP_EXT(MA,MB,MD2,w))
          AMPSCA = (WD1 * AMP_SCA(MA,MB,MD1,w) + WD2 * AMP_SCA(MA,MB,MD2,w))
          AMPASY = (WD1 * AMP_ASY(MA,MB,MD1,w) + WD2 * AMP_ASY(MA,MB,MD2,w))
       ENDIF
c--------------------------------------------------------------------
          EXT(l,w) = EXT(l,w) + ( AMPEXT * TTAUSV(l,n) * FSTOPX(n))
          HELP     = ((GCB(l,w) * SCT(l,w) ) + (AMPASY * AMPSCA * TTAUSV(l,n)* FSTOPX(n)))
          SCT(l,w) = SCT(l,w) +  (AMPSCA * TTAUSV(l,n) * FSTOPX(n))
          GCB(l,w) = HELP / (SCT(l,w)+ 1.D-10)

          aesqex(l,w,n)= AMPEXT * TTAUSV(l,n)
          aesqsc(l,w,n)= AMPSCA * TTAUSV(l,n)
          aesqcb(l,w,n)= AMPASY * aesqsc(l,w,n)

      ENDDO   ! wave
      ENDDO   ! modes
      ENDDO   ! level

        ENDIF     ! AMP_RAD_KEY = 2


C Longwave: ---------------------------------------------------------------------------------------------
c     The absorption spectrum of each mode is the volume-fraction weighted average of the
c     spectra of the pure compositions it contains, using the dry volume fractions that
c     SETAMP_LEV computes for this column. This mirrors the volume mixing already applied
c     to the shortwave refractive index in SETAMP_LEV, so both parts of the spectrum see
c     the same composition.
c
c     Previously each mode used a single fixed composition (CORE_CLASS), which was
c     workable while only the carbonaceous modes carried organics. It is not workable for
c     the MATRIX-VBS mechanisms (M9, and M10 for its non-volatile organics), where every
c     mode except AKK carries organic mass: the dust and sea salt modes would otherwise
c     absorb as pure dust and pure sea salt no matter how much organic had condensed onto
c     them. See MSPCS in AERO_CONFIG for which modes carry which species.
c
c     dry_Vf_LEV(l,n,1:NDRY) is normalised over the dry species, so the weights sum to one
c     for any mode holding mass, and to zero for an empty mode (which then contributes
c     nothing, consistent with its TTAUSV also being zero).

      DO l = 1,lm
      DO n = 1,NMODES
         AMP_TAB(:) = 0.d0
         DO s = 1,NDRY
            AMP_TAB(:) = AMP_TAB(:)
     +                 + dry_Vf_LEV(l,n,s) * AMP_TAB_CLASS(:,SPC_TO_LWC(s),n)
         ENDDO   ! dry species
         TAB(l,:) = TAB(l,:) + (AMP_TAB(:) *  TTAUSV(l,n) * FTTOPX(n))
      ENDDO   ! modes
      ENDDO   ! level
      endif

      RETURN
      END SUBROUTINE SETAMP
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETAMP_LEV(i,j,l)
!@sum Calulates effective Radius and Refractive Index for Mixed Aerosols
!@sum Puts AMP Aerosols in 1 dimension CALLED in RADIA
!@auth Susanne Bauer

      USE AMP_Optics_mod

#ifdef USE_OFFLINE_AEROSOLS
      USE OFFLINE_AEROSOL, only :  DIAM,Reff_LEV, NUMB_LEV, RindexAMP,
     +  dry_Vf_LEV,MIX_OC,MIX_SU,MIX_AQ,AMP_RAD_KEY,
     +  NMODES
      USE CONSTANT,   only: rgas, pi
      use ATMCOL_COM, only: tl   ! layer temperature (K)
      use ATMCOL_COM, only: pl   ! layer pressure (mb)
      USE ATM_COM, only: MA,byMA  ! Air mass of each box (kg m-2)
      USE RESOLUTION, only: lm
      USE aeractv_streams_mod, only : actvqtys, actvqtys2

#else
      USE AMP_AEROSOL, only: DIAM,Reff_LEV, NUMB_LEV, RindexAMP,
     +  dry_Vf_LEV,MIX_OC,MIX_SU,MIX_AQ,AMP_RAD_KEY
      USE AmpTracersMetadata_mod,  only: AMP_NUMB_MAP,
     +  AMP_MODES_MAP
      USE TRACER_COM,  only: TRM, ntmAMPi,ntmAMPe
      USE AERO_CONFIG, only: NMODES
      USE AERO_SETUP,  only: SIG0, CONV_DPAM_TO_DGN   !(nmodes * npoints) lognormal parameters for each mode
      USE AERO_SETUP,  only: MODE_NAME

      USE AERO_ACTV, only: DENS_SULF, DENS_DUST,DENS_SEAS, DENS_BCAR, DENS_OCAR
      use OldTracer_mod, only: trname
      use AMP_utilities_mod, only: aerosolkind
#endif

      IMPLICIT NONE

      ! Arguments:
      INTEGER, INTENT(IN) :: i,j,l

      ! Local
      INTEGER n,w,s,nAMP,k
      REAL*8,     DIMENSION(NMODES,NSPC) :: VolFrac, VMass
      REAL*8                             :: H2O, NO3
      REAL(8), PARAMETER :: TINYNUMER = 1.0D-30
      COMPLEX*8, DIMENSION(NSW,NSPC)     :: Ri
c     Variables for Maxwell Garnett:
      REAL*8                           :: V_bc, V_host
      COMPLEX*8                        :: M_mg, M_bc, M_host
c Andies data incl Solar weighting - integral over 6 radiation band
c     One column of six shortwave values per species, in the SPC_* order.
      DATA Ri/(1.46099,    0.0764233)  ,(1.48313,  0.000516502),    !Su
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

c     &        (1.80056,     0.605467)  ,(1.68622,     0.583112),    !Bc
c     &        (1.63586,     0.551897)  ,(1.59646,     0.515333),
c     &        (1.57466,     0.484662)  ,(1.56485,     0.487992),
c only 550nm values
c     &        (1.85,     0.71)  ,(1.85,     0.71),    !Bc
c     &        (1.85,     0.71)  ,(1.85,     0.71),
c     &        (1.85,     0.71)  ,(1.85,     0.71),
cBond + Berstroem, all wavelength
     &        (2.15,     1.05)  ,(1.98,     0.79),    !Bc
     &        (1.92,     0.75)  ,(1.86,     0.73),
     &        (1.85,     0.71)  ,(1.85,     0.71),

     &        (1.46099,    0.0761930)  ,(1.48313,   0.00470000),    !Oc
     &        (1.49719,   0.00470000)  ,(1.50805,   0.00480693),
     &        (1.52000,   0.00540000)  ,(1.52775,    0.0144927),

     &        (1.47978,    0.0211233)  ,(1.50719,   0.00584169),    !Du
     &        (1.51608,   0.00378434)  ,(1.52998,   0.00178703),
     &        (1.54000,  0.000800000)  ,(1.56448,   0.00221463),

     &        (1.46390,   0.00571719)  ,(1.45000,      0.00000),    !Ss
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),

     &        (1.46099,  0.0764233)    ,(1.48313,  0.000516502),    !No3
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

     &        (1.26304,    0.0774872)  ,(1.31148,  0.000347758),    !H2O
     &        (1.32283,  0.000115835)  ,(1.32774,  3.67435e-06),
     &        (1.33059,  1.58222e-07)  ,(1.33447,  3.91074e-08)/

#ifdef USE_OFFLINE_AEROSOLS
c     Layout of the offline aerosol streams: actvqtys holds the NOFF_SPCS dry mass
c     species of each mode followed by its number and its diameter of average mass,
c     and actvqtys2 the three mode-independent extra tracers.
      INTEGER, PARAMETER :: NOFF_SPCS = SPC_SS   ! SU, BC, OC, DU, SS
      INTEGER, PARAMETER :: IQ_NUMB   = NOFF_SPCS + 1
      INTEGER, PARAMETER :: IQ_DPAM   = NOFF_SPCS + 2
      INTEGER, PARAMETER :: IQ2_NO3   = 1
      INTEGER, PARAMETER :: IQ2_NH4   = 2
      INTEGER, PARAMETER :: IQ2_H2O   = 3
      REAL(8) :: NI(NMODES)             ! number concentration for each tracer [#/m^3]
      REAL(8) :: M(NMODES,NOFF_SPCS)    ! mass   concentration for each tracer [ug/m^3]
      REAL(8) :: DG_WET(NMODES)         ! geometric mean diameter for each dry tracer [um]
      REAL(8) :: DGN(NMODES)            ! geometric mean diameter for each dry tracer [um]
      REAL(8) :: TK               ! absolute temperature [K]
      REAL(8) :: PRES             ! ambient pressure [Pa]
      REAL(8) :: AIRD(lm)          ! air density [kg/m^3]
      REAL(8) :: VOLTMP_WET
      real*8, parameter, dimension(NMODES) ::
     &    sig0=(/ 1.6d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0,
     &                 2.0d0, 2.0d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0,
     &                 1.8d0, 1.8d0, 1.8d0, 2.0d0/)
      real*8, parameter, dimension(NMODES) ::
     &    CONV_DPAM_TO_DGN=(/ 0.71795016727196403,   0.59556797724590516     ,
     &     0.59556797724590516     ,  0.59556797724590516   ,    0.59556797724590516  ,
     &     0.59556797724590516      , 0.59556797724590516   ,    0.48642160999311468  ,
     &     0.59556797724590516   ,    0.59556797724590516  ,     0.59556797724590516  ,
     &     0.59556797724590516  ,     0.59556797724590516    ,   0.59556797724590516  ,
     &     0.59556797724590516 ,      0.48642160999311468/)

      real*8, parameter, dimension(NOFF_SPCS) ::
     &    dens=(/1.77D+03,1.70D+03,1.00D+03,2.60D+03,2.165D+03/)   ! [kg/m^3]
      CHARACTER(LEN=3), PARAMETER, DIMENSION(NMODES) :: MODE_NAME = (/
     &  'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC',
     &  'OCC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/)
#endif

c     VMass must be zeroed here: not every (mode,species) pair is assigned
c     below, and the organics are accumulated over several tracers (VBS).
      VMass(:,:) = 0.d0

#ifdef USE_OFFLINE_AEROSOLS

       ! + Effective Radius [um] per Mode = geometric mass mean radius
! meteorology
        call load_atmcol(i,j)   ! remove this in E3
           TK = tl(l)                            ! [K]
           PRES= pl(l)*100.d0                    ! pmid in [hPa]
           AIRD(l) = PRES/(rgas * TK)            ! air density  [kg/m3]
!
      do n=1,NMODES ! loop over modes
           VOLTMP_WET = 1.0D-30
           do k = 1, NOFF_SPCS
             M(n,k) = actvqtys(l,n,k,i,j) * MA(l,i,j)               ! mass [kg/kg]  -> [kg/m2/layer]
             VOLTMP_WET = VOLTMP_WET + actvqtys(l,n,k,i,j)/dens(k)  ! dry volume [m3]
           enddo

             VOLTMP_WET = VOLTMP_WET
     &         + (M(n,SPC_SU)/sum(M(:,SPC_SU))) * actvqtys2(l,IQ2_NO3,i,j)/DENS_NO3  ! NO3
             VOLTMP_WET = VOLTMP_WET
     &         + (M(n,SPC_SU)/sum(M(:,SPC_SU))) * actvqtys2(l,IQ2_NH4,i,j)/DENS_NO3  ! NH4
             VOLTMP_WET = VOLTMP_WET
     &         + (M(n,SPC_SU)/sum(M(:,SPC_SU))) * actvqtys2(l,IQ2_H2O,i,j)/DENS_AQ   ! H2O

         DG_WET(n) = 1.D6 *( (6.d0/pi) * (VOLTMP_WET / actvqtys(l,n,IQ_NUMB,i,j)) ) ! dry gemoetric mean mass diameter [um]
     &                 **0.333333333333333  ! [um]
         DG_WET(n) = MIN( MAX( DG_WET(n),  0.01), 10.D0 )

         DGN(n) = DG_WET(n) * (1.0D+00 / EXP( 1.5D+00 * ( log(sig0(n)) )**2 ))

c         Reff_LEV(l,n) = DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5
         Reff_LEV(l,n) = actvqtys(l,n,IQ_DPAM,i,j)*CONV_DPAM_TO_DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5e6

          VMass(n,1:NOFF_SPCS)  = M(n,1:NOFF_SPCS)/DENS(1:NOFF_SPCS)

       enddo



       ! + Volume Fraction
       DO n=1,NMODES  ![#/m2]         pi/4     [m2]
         NI(n)  =  actvqtys(l,n,IQ_NUMB,i,j) * MA(l,i,j)     ! number [#/kg] *  [kg/m2] = [#/m2]
c         NUMB_LEV(l,n) = NI(n)* 0.7853 * (1.e-6*DG_WET(n))**2   ! [#/layer]
         NUMB_LEV(l,n) = NI(n)* 0.7853 *actvqtys(l,n,IQ_DPAM,i,j)**2

        ! NO3
        VMass(n,SPC_NO3) = VMass(n,SPC_SU) / Sum(VMass(:,SPC_SU))
     &    * ((actvqtys2(l,IQ2_NO3,i,j) + actvqtys2(l,IQ2_NH4,i,j))* MA(l,i,j) )/ DENS_NO3
        ! H2O
        VMass(n,SPC_H2O) = VMass(n,SPC_SU) / Sum(VMass(:,SPC_SU))
     &    * (actvqtys2(l,IQ2_H2O,i,j) * MA(l,i,j) )/DENS_AQ
       ENDDO

#else   /* Original MATRIX online code below */

       ! + Effective Radius [um] per Mode = geometric mass mean radius
       DO n=1,NMODES
         Reff_LEV(l,n) = DIAM(i,j,l,n)*CONV_DPAM_TO_DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5e6
       ENDDO

       ! + Mass and Number Concentration
       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
           if(trname(n) .eq.'M_NO3') NO3 =trm(i,j,l,n)
           if(trname(n) .eq.'M_H2O') H2O =trm(i,j,l,n)
           if(AMP_NUMB_MAP(nAMP).eq. 0) then  ! Volume fraction
             if (trname(n).ne.'M_NO3'.and.trname(n).ne.'M_H2O'.and.
     &           trname(n).ne.'M_NH4') then
               select case (aerosolKind(trname(n)))
               case ('SU')
                  VMass(AMP_MODES_MAP(nAMP),SPC_SU) =trm(i,j,l,n)/DENS_SULF
               case ('BC')
                  VMass(AMP_MODES_MAP(nAMP),SPC_BC) =trm(i,j,l,n)/DENS_BCAR
               case ('OC')
                  VMass(AMP_MODES_MAP(nAMP),SPC_OC) =VMass(AMP_MODES_MAP(nAMP),SPC_OC)+ trm(i,j,l,n)/DENS_OCAR
               case ('DU')
                  VMass(AMP_MODES_MAP(nAMP),SPC_DU) =trm(i,j,l,n)/DENS_DUST
               case ('SS')
                  VMass(AMP_MODES_MAP(nAMP),SPC_SS) =trm(i,j,l,n)/DENS_SEAS
               end select
             endif
           else                           ! Number
!          [ - ]                        [trm units: #/m2/layer]
           NUMB_LEV(l,AMP_NUMB_MAP(nAMP)) =trm(i,j,l,n)
          endif
       ENDDO

       ! + Volume Fraction
       DO n=1,NMODES  ![#/m2]         pi/4     [m2]
        NUMB_LEV(l,n) = NUMB_LEV(l,n)* 0.7853 * DIAM(i,j,l,n)**2
        ! NO3
        VMass(n,SPC_NO3) = VMass(n,SPC_SU) / Sum(VMass(:,SPC_SU)) * NO3/ DENS_NO3
        ! H2O
        VMass(n,SPC_H2O) = VMass(n,SPC_SU) /(Sum(VMass(:,SPC_SU)) + TINYNUMER)  * H2O /DENS_AQ
       ENDDO
#endif   /* After that code should work for all cases */

      DO s=1,NSPC  ! loop over species
        DO n=1,NMODES           ! loop over modes
          Volfrac(n,s) = VMass(n,s) / (Sum(VMass(n,:)) + TINYNUMER)
          dry_Vf_LEV(l,n,s) = VMass(n,s) / (Sum(VMass(n,1:NDRY)) + TINYNUMER)
        ENDDO
      ENDDO

      ! Core Shell Composition: core is BC, shell material is OC, SO4 and H2O.
      ! IS_CORE_SHELL picks out the same modes that SETAMP reads from the
      ! core-shell tables, so the two cannot drift apart.
      DO n=1,NMODES             ! loop over modes
        IF ( IS_CORE_SHELL(MODE_NAME(n)) ) THEN
          MIX_OC(l,n) = VMass(n,SPC_OC)
     &      / (VMass(n,SPC_SU) + VMass(n,SPC_BC) + VMass(n,SPC_OC) + VMass(n,SPC_H2O) + TINYNUMER)
          MIX_SU(l,n) = VMass(n,SPC_SU)
     &      / (VMass(n,SPC_SU) + VMass(n,SPC_BC) + VMass(n,SPC_OC) + VMass(n,SPC_H2O) + TINYNUMER)
          MIX_AQ(l,n) = VMass(n,SPC_H2O)
     &      / (VMass(n,SPC_SU) + VMass(n,SPC_BC) + VMass(n,SPC_OC) + VMass(n,SPC_H2O) + TINYNUMER)
        ENDIF
      ENDDO

      ! + Refractive Index of Aerosol mix per mode and wavelength

      RindexAMP(l,:,:) = 0.d0
      DO s=1,NSPC               ! loop over species
        DO w=1,NSW              ! loop over wavelength
          DO n=1,NMODES         ! loop over modes
            RindexAMP(l,n,w) = RindexAMP(l,n,w) + ( Volfrac(n,s) * Ri(w,s))
          ENDDO
        ENDDO
      ENDDO

      if (AMP_RAD_KEY == 3) then      ! - - - Maxwell Garnett Mixing Rule
        DO w=1,NSW              ! loop over wavelength
          DO n=1,NMODES         ! loop over modes
       IF ( IS_CORE_SHELL(MODE_NAME(n)) ) THEN
             M_bc   = Ri(w,SPC_BC)
             V_bc   = Volfrac(n,SPC_BC)
             M_host = (0.d0, 0.d0)
             V_host = 0.d0
             DO s=1,NSPC               ! loop over species other than BC
               IF ( s .NE. SPC_BC ) THEN
                 M_host = M_host + ( Volfrac(n,s) * Ri(w,s))
                 V_host = V_host + Volfrac(n,s)
               ENDIF
             ENDDO
             M_mg = M_host**2  * (M_bc**2 + 2.d0 * M_host**2 + 2.d0 * V_bc * (M_bc    - M_host   ) )
     +                         / (M_bc**2 + 2.d0 * M_host**2 -        V_host*(M_bc**2 - M_host**2) )

             RindexAMP(l,n,w) = SQRT( M_mg)
       ENDIF
          ENDDO
        ENDDO
      endif

      RETURN
      END SUBROUTINE SETAMP_LEV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETUP_RAD
!@sum Initialization for Radiation incl. Aerosol Microphysics
!@auth Susanne Bauer

      USE AMP_Optics_mod

#ifdef USE_OFFLINE_AEROSOLS
      USE OFFLINE_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       AMP_Q55, dry_Vf_LEV, RindexAMP
#else
      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA, AMP_Q55,
     +           AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +           dry_Vf_LEV, RindexAMP
#endif
      IMPLICIT NONE
      include 'netcdf.inc'
      integer start(4),count(4),count3(3),status
      integer start2(5),count2(5),count32(4)
      integer ncid, id1, id2, id3, id4,ncid2

      real*4, DIMENSION(NRE,NIM,NSIZ,NSW) ::  ASY,SCA,EXT
      real*4, DIMENSION(NRE,NIM,NSIZ) ::    QEX
      real*4, DIMENSION(NSIZ,NCS,NCS,NCS,NSW) ::  CS_ASY,CS_SCA,CS_EXT
      real*4, DIMENSION(NSIZ,NCS,NCS,NCS) ::    CS_QEX
c -----------------------------------------------------------------
c   The table sizes in AMP_Optics_mod are what this routine reads with, and what
c   SETAMP indexes with, but the arrays themselves are declared in AMP_AEROSOL
c   (or OFFLINE_AEROSOL). Check the two agree rather than trusting them to: a
c   silent mismatch here would corrupt every optical property in the model.
c -----------------------------------------------------------------
      IF ( SIZE(AMP_EXT,1)    .NE. NRE  .OR. SIZE(AMP_EXT,2)    .NE. NIM  .OR.
     +     SIZE(AMP_EXT,3)    .NE. NSIZ .OR. SIZE(AMP_EXT,4)    .NE. NSW  .OR.
     +     SIZE(AMP_Q55,1)    .NE. NRE  .OR. SIZE(AMP_Q55,2)    .NE. NIM  .OR.
     +     SIZE(AMP_Q55,3)    .NE. NSIZ .OR.
     +     SIZE(AMP_EXT_CS,1) .NE. NSIZ .OR. SIZE(AMP_EXT_CS,2) .NE. NCS  .OR.
     +     SIZE(AMP_EXT_CS,3) .NE. NCS  .OR. SIZE(AMP_EXT_CS,4) .NE. NCS  .OR.
     +     SIZE(AMP_EXT_CS,5) .NE. NSW  .OR.
     +     SIZE(RindexAMP,3)  .NE. NSW  .OR. SIZE(dry_Vf_LEV,3) .NE. NSPC ) THEN
        WRITE(*,*) 'SETUP_RAD: Mie table or species dimensions in AMP_AEROSOL do not'
        WRITE(*,*) 'match NRE,NIM,NSIZ,NCS,NSW,NSPC in AMP_Optics_mod.'
        CALL stop_model('SETUP_RAD: inconsistent aerosol optics dimensions',255)
      ENDIF
c -----------------------------------------------------------------
c   Opening of the files to be read: MIE TABLES
c -----------------------------------------------------------------
          status=NF_OPEN('AMP_MIE_TABLES',NCNOWRIT,ncid)
          status=NF_INQ_VARID(ncid,'ASYM',id1)
          status=NF_INQ_VARID(ncid,'QEXT',id2)
          status=NF_INQ_VARID(ncid,'QSCT',id3)
          status=NF_INQ_VARID(ncid,'Q55E',id4)
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start(:)=1

          count(1)=NRE
          count(2)=NIM
          count(3)=NSIZ
          count(4)=NSW
          count3(1)=NRE
          count3(2)=NIM
          count3(3)=NSIZ


          status=NF_GET_VARA_REAL(ncid,id1,start,count,ASY)
          status=NF_GET_VARA_REAL(ncid,id2,start,count,EXT)
          status=NF_GET_VARA_REAL(ncid,id3,start,count,SCA)
          status=NF_GET_VARA_REAL(ncid,id4,start,count3,QEX)

c         NF_CLOSE takes only the netCDF id; it used to be handed the file name
c         and NCNOWRIT as well, so it closed a garbage id and left the file open.
          status=NF_CLOSE(ncid)

          AMP_ASY = ASY
          AMP_EXT = EXT
          AMP_SCA = SCA
          AMP_Q55 = QEX
c -----------------------------------------------------------------
c   Opening of the files to be read: Core Shell Mie Tables
c   Core is BC, Shell Material is OC, SO4 and H2O
c -----------------------------------------------------------------
          status=NF_OPEN('AMP_CORESHELL_TABLES',NCNOWRIT,ncid2)
          status=NF_INQ_VARID(ncid2,'CS_ASYM',id1)
          status=NF_INQ_VARID(ncid2,'CS_QEXT',id2)
          status=NF_INQ_VARID(ncid2,'CS_QSCT',id3)
          status=NF_INQ_VARID(ncid2,'CS_Q55E',id4)
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start2(:)=1

          count2(1)=NSIZ
          count2(2)=NCS
          count2(3)=NCS
          count2(4)=NCS
          count2(5)=NSW
          count32(1)=NSIZ
          count32(2)=NCS
          count32(3)=NCS
          count32(4)=NCS


          status=NF_GET_VARA_REAL(ncid2,id1,start2,count2,CS_ASY)
          status=NF_GET_VARA_REAL(ncid2,id2,start2,count2,CS_EXT)
          status=NF_GET_VARA_REAL(ncid2,id3,start2,count2,CS_SCA)
          status=NF_GET_VARA_REAL(ncid2,id4,start2,count32,CS_QEX)

          status=NF_CLOSE(ncid2)

          AMP_ASY_CS = CS_ASY
          AMP_EXT_CS = CS_EXT
          AMP_SCA_CS = CS_SCA
          AMP_Q55_CS = CS_QEX
      RETURN
      END SUBROUTINE SETUP_RAD
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE GET_LW(NA,NS,AREFF,TQAB,Vf)
!@sum Calculation of LW absorption for AMP aerosols
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer

      USE AMP_Optics_mod, only: NLW, NDRY, NREFU, NREFS, NREFD,
     +     LWC_BCAR, LWC_DUST, LW_TRUQEX_OFFSET
      USE RADPAR, only: TRUQEX, TRSQEX, TRDQEX, TRUQSC, TRSQSC, TRDQSC
     *                  , REFU22, REFS25, REFD25
      IMPLICIT NONE
      INTEGER, intent(IN) :: NA,NS
      REAL*8,  intent(in) :: areff,Vf(NDRY)
      REAL*8,  intent(out):: TQAB(NLW)
      REAL*8   TQEX(NLW),TQSC(NLW),TQEX_S(NLW),TQSC_S(NLW)
      REAL*8   QXAERN(NREFS),QSAERN(NREFS)
      REAL*8   wts,wta
      INTEGER  n0,k,n,nn

c CORE
       IF(NA==0) THEN
         TQAB(:)= 0d0
       ENDIF
c     Compositions held in the 22-size blocks of TRUQEX: SO4, SEAS, NO3 and OCAR.
      IF(NA > 0 .and. NA < LWC_BCAR) THEN
        N0=LW_TRUQEX_OFFSET(NA)
        DO K=1,NLW
        DO N=1,NREFU
        NN=N0+N
        QXAERN(N)=TRUQEX(K,NN)
        QSAERN(N)=TRUQSC(K,NN)
        ENDDO
        CALL SPLINE(REFU22,QXAERN,NREFU,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,NREFU,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF

      IF(NA==LWC_BCAR) THEN           !   NA : Aerosol compositions BC
        DO K=1,NLW
        QXAERN(:)=TRSQEX(K,:)    ! 1:NREFS
        QSAERN(:)=TRSQSC(K,:)    ! 1:NREFS
        CALL SPLINE(REFS25,QXAERN,NREFS,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,NREFS,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF

      IF(NA==LWC_DUST) THEN           !   NA : Aerosol composition DST
       DO K=1,NLW
        QXAERN(:)=TRDQEX(K,:)    ! 1:NREFD
        QSAERN(:)=TRDQSC(K,:)    ! 1:NREFD
        CALL SPLINE(REFD25,QXAERN,NREFD,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,NREFD,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
      ENDDO
      ENDIF

c SHELL
c     Vf holds the shell volume fraction of each species, in the SPC_* order of
c     AMP_Optics_mod. Restoring this branch would need SPC_SU and SPC_SS from there.
c         IF(NS > 0 .and. NS < LWC_BCAR) THEN    !    NS : Aerosol compositions SO4,SEA,NO3,OC
c         N0=LW_TRUQEX_OFFSET(NS)


c         DO K=1,NLW
c         DO N=1,NREFU
c         NN=N0+N
c         IF (NS==LWC_SO4)  WTS=Vf(SPC_SU)          ! <- shell fraction of aerosol composition
c         IF (NS==LWC_SEAS) WTS=Vf(SPC_SS)          ! <- shell fraction of aerosol composition
c         WTA=1.D0-WTS
c         QXAERN(N)=TRUQEX(K,NN)
c         QSAERN(N)=TRUQSC(K,NN)
c         ENDDO
c         CALL SPLINE(REFU22,QXAERN,NREFU,AREFF,TQEX_S(K),1.D0,1.D0,1)
c         CALL SPLINE(REFU22,QSAERN,NREFU,AREFF,TQSC_S(K),1.D0,1.D0,1)
c         TQAB(K)=(TQEX(K)*WTA + TQEX_S(K)*WTS)-(TQSC(K)*WTA + TQSC_S(K)*WTS)
c         ENDDO

c      ENDIF

      RETURN
      END SUBROUTINE GET_LW
c -----------------------------------------------------------------

