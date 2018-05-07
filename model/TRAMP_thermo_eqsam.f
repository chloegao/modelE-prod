      SUBROUTINE AERO_THERMO(ASO4,ANO3,ANH4,DUST,SALT,AH2O,ApH,SSH2O,
     &                       GNH3,GHNO3,TEMP,RH,RHD,RHC)
!@sum
!@+     This routine sets up for and calls the thermodynamic module for aerosol
!@+     gas-particle partitioning.
!@+
!@+      A version of EQSAM (eqsam_v03d) is the current thermodynamic model. 
!@auth Susanne Bauer/Doug Wright

!----------------------------------------------------------------------------------------------------------------------
!     This routine sets up for and calls the thermodynamic module for aerosol
!     gas-particle partitioning.
!
!     A version of EQSAM (eqsam_v03d) is the current thermodynamic model. 
!
!     EQSAM is called with control variable IOPT=1. 
!
!     Although EQSAM takes as input the total S(VI) (H2SO4+SO4=), since the
!     aerosol model does not necessarily transfer all H2SO4 to the aerosol
!     phase (depending on configuration), we pass only the particulate SO4
!     as the total sulfate to EQSAM.
!
!     Also, this version of EQSAM takes as input the mineral cation 
!     concentrations K+, Ca++, Mg++, Na+. Given the 'well-mixed' treatment
!     of inorganic aerosol constituents, these cations are included.
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_PARAM, ONLY: WRITE_LOG, AUNIT1
      IMPLICIT NONE

      ! Arguments.
      REAL(8), INTENT(IN)    :: ASO4      ! aerosol sulfate       [ug/m^3]
      REAL(8), INTENT(INOUT) :: ANO3      ! aerosol nitrate       [ug/m^3]
      REAL(8), INTENT(INOUT) :: ANH4      ! aerosol ammonium      [ug/m^3]
      REAL(8), INTENT(IN)    :: DUST      ! dust                  [ug/m^3]
      REAL(8), INTENT(IN)    :: SALT      ! sea salt              [ug/m^3]
      REAL(8), INTENT(OUT)   :: AH2O      ! aerosol water         [ug/m^3]
      REAL(8), INTENT(OUT)   :: ApH       ! aerosol pH
      REAL(8), INTENT(OUT)   :: SSH2O     ! sea salt assoc. H2O   [ug/m^3]
      REAL(8), INTENT(INOUT) :: GNH3      ! gas-phase ammonia     [ugNH4/m^3]
      REAL(8), INTENT(INOUT) :: GHNO3     ! gas-phase nitric acid [ugNO3/m^3]
      REAL(8), INTENT(IN)    :: TEMP      ! temperature           [K]
      REAL(8), INTENT(IN)    :: RH        ! relative humidity     [0-1]
      REAL(8), INTENT(OUT)   :: RHD       ! RH of deliquescence   [0-1]
      REAL(8), INTENT(OUT)   :: RHC       ! RH of crystallization [0-1]

      ! Call parameters for the EQSAM thermodynamic model. 
      INTEGER, PARAMETER :: NCA  = 11    ! fixed number of input variables
      INTEGER, PARAMETER :: NCO  = 37    ! fixed number of output variables
      INTEGER, PARAMETER :: IOPT =  1    ! =1 selects the metastable (wet) state and history
!     INTEGER, PARAMETER :: IOPT =  2    ! =2 selects the solid      (dry) state and history
      INTEGER, PARAMETER :: LOOP =  1    ! only a single time step done
      INTEGER, PARAMETER :: IMAX =  1    ! only a single time step done

      REAL :: YI(IMAX,NCA)            ! [umol/m^3] for chemical species - input
      REAL :: YO(IMAX,NCO)            ! [umol/m^3] for chemical species - output

      ! Parameters.

      REAL(8), PARAMETER :: MW_ANH4   = 18.03850  ! [g/mol]
      REAL(8), PARAMETER :: MW_GNH3   = MW_ANH4   ! [g/mol] NH3  is passed as equivalent conc. of NH4+
      REAL(8), PARAMETER :: MW_ANO3   = 62.00494  ! [g/mol]
      REAL(8), PARAMETER :: MW_GHNO3  = MW_ANO3   ! [g/mol] HNO3 is passed as equivalent conc. of NO3-
      REAL(8), PARAMETER :: MW_ASO4   = 96.0636   ! [g/mol]
      REAL(8), PARAMETER :: MW_K      = 39.0983   ! [g/mol]
      REAL(8), PARAMETER :: MW_CA     = 40.078    ! [g/mol]
      REAL(8), PARAMETER :: MW_MG     = 24.3050   ! [g/mol]
      REAL(8), PARAMETER :: MW_NA     = 22.989768 ! [g/mol]
      REAL(8), PARAMETER :: MW_NACL   = 58.442468 ! [g/mol]

      REAL(8), PARAMETER :: MASS_FRAC_K  = 0.0028  ! From Ghan et al. (2001).
      REAL(8), PARAMETER :: MASS_FRAC_CA = 0.024   !   JGR, Vol. 106, p. 5295-5316.
      REAL(8), PARAMETER :: MASS_FRAC_MG = 0.0038  !   on p. 5296
      REAL(8), PARAMETER :: MASS_FRAC_NA = 0.014   !   "water sol. mass frac. in soil dust"

!----------------------------------------------------------------------------------------------------------------------
!     These values were used for development and testing.
!----------------------------------------------------------------------------------------------------------------------
!     REAL(8), PARAMETER :: MASS_FRAC_K  = 0.0000  
!     REAL(8), PARAMETER :: MASS_FRAC_CA = 0.0000  
!     REAL(8), PARAMETER :: MASS_FRAC_MG = 0.0000  
!     REAL(8), PARAMETER :: MASS_FRAC_NA = 0.0000  
!----------------------------------------------------------------------------------------------------------------------

      REAL(8), PARAMETER :: FRAC_DUST  = 0.1                              ! [1] fraction of dust conc. passed to EQSAM         
      REAL(8), PARAMETER :: CONV_KION  = FRAC_DUST * MASS_FRAC_K  / MW_K  ! [mol/g]
      REAL(8), PARAMETER :: CONV_CAION = FRAC_DUST * MASS_FRAC_CA / MW_CA ! [mol/g]
      REAL(8), PARAMETER :: CONV_MGION = FRAC_DUST * MASS_FRAC_MG / MW_MG ! [mol/g]
      REAL(8), PARAMETER :: CONV_NAION = FRAC_DUST * MASS_FRAC_NA / MW_NA ! [mol/g]

      REAL(8), PARAMETER :: RMW_GNH3  = 1.0 / MW_GNH3          ! [mol/g]
      REAL(8), PARAMETER :: RMW_ANH4  = 1.0 / MW_ANH4          ! [mol/g]
      REAL(8), PARAMETER :: RMW_GHNO3 = 1.0 / MW_GHNO3         ! [mol/g]
      REAL(8), PARAMETER :: RMW_ANO3  = 1.0 / MW_ANO3          ! [mol/g]
      REAL(8), PARAMETER :: RMW_ASO4  = 1.0 / MW_ASO4          ! [mol/g]
      REAL(8), PARAMETER :: RMW_NA    = 1.0 / MW_NA            ! [mol/g]
      REAL(8), PARAMETER :: RMW_NACL  = 1.0 / MW_NACL          ! [mol/g]

      REAL(8), PARAMETER :: DH2O   = 1.00D+00    ! density of water [g/cm^3]
      REAL(8), PARAMETER :: DNACL  = 2.165D+00   ! density of NaCl  [g/cm^3]
      REAL(8), PARAMETER :: CSS    = 1.08D+00    ! for sea salt ...
      REAL(8), PARAMETER :: BSS    = 1.2D+00     ! for sea salt ...              
      REAL(8), PARAMETER :: SSH2OA = (CSS*CSS*CSS*BSS-1.0D+00)*DH2O/DNACL
      REAL(8), PARAMETER :: SSH2OB = (CSS*CSS*CSS            )*DH2O/DNACL
      REAL(8), PARAMETER :: RHMAX  = 0.995D+00   ! [0-1]
      REAL(8), PARAMETER :: RHMIN  = 0.010D+00   ! [0-1]   
      REAL(8), PARAMETER :: SMALL_SO4 = 1.0D-05  ! [umol SO4/m^3] EQSAM has crashed at low RH and low sulfate conc.

      REAL(8) :: H   ! local RH, with RHMIN < H < RHMAX

      !----------------------------------------------------------------------------------------------------------------
      ! Call for the bulk non-sea salt inorganic aerosol.
      !----------------------------------------------------------------------------------------------------------------
      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(/A,2F12.3/)') 'EQSAM: TEMP[K], RH[0-1]= ', TEMP, RH
        WRITE(AUNIT1,'(A4,7A14  )') '   ','ASO4','ANO3','ANH4','AH2O', 'GNH3','GHNO3','DUST'
        WRITE(AUNIT1,'(A4,7E14.5)') 'TOP',ASO4,ANO3,ANH4,AH2O,GNH3,GHNO3,DUST
      ENDIF

      H = MAX( MIN( RH, RHMAX ), RHMIN )

      YI(1,1)  = TEMP                             ! [K]
      YI(1,2)  = H                                ! [0-1]
      YI(1,3)  = GNH3*RMW_GNH3   + ANH4*RMW_ANH4  ! from [ug/m^3] to [umol/m^3]
      YI(1,4)  =                   ASO4*RMW_ASO4  ! from [ug/m^3] to [umol/m^3]
      YI(1,5)  = GHNO3*RMW_GHNO3 + ANO3*RMW_ANO3  ! from [ug/m^3] to [umol/m^3]
      YI(1,6)  = DUST*CONV_NAION                  ! from [ug dust/m^3] to [umol Na+/m^3]
      YI(1,7)  = 0.0                              ! (HCl + Cl-)
      YI(1,8)  = DUST*CONV_KION                   ! from [ug dust/m^3] to [umol K+ /m^3]
      YI(1,9)  = DUST*CONV_CAION                  ! from [ug dust/m^3] to [umol Ca+/m^3]
      YI(1,10) = DUST*CONV_MGION                  ! from [ug dust/m^3] to [umol Mg+/m^3]
      YI(1, :) = MAX( YI(1,:), 0.0d-10 )          ! Lower limit was 1.0E-10 before 102406.
      YI(1,4)  = YI(1,4) + SMALL_SO4              ! EQSAM has crashed at low RH and low sulfate conc.

      CALL EQSAM_V03D(YI,YO,NCA,NCO,IOPT,LOOP,IMAX,AUNIT1)

      GHNO3 = MAX(YO(1, 9) * MW_GHNO3, 0.d0 )     ! from [umol/m^3] to [ug/m^3]
      GNH3  = MAX(YO(1,10) * MW_GNH3 , 0.d0 )     ! from [umol/m^3] to [ug/m^3]
      AH2O  = MAX(YO(1,12)           , 0.d0 )     ! already in [ugH2O/m^3]
      ApH   = -log10(YO(1,37)+tiny(1.e0))
      ANH4  = MAX(YO(1,19) * MW_ANH4 , 0.d0 )     ! from [umol/m^3] to [ug/m^3]
      ANO3  = MAX(YO(1,20) * MW_ANO3 , 0.d0 )     ! from [umol/m^3] to [ug/m^3]
! eqsam does not modify ASO4, so the lines below are not needed
!      ASO4  = ( YO(1,21) - SMALL_SO4 ) * MW_ASO4  ! from [umol/m^3] to [ug/m^3]
!      ASO4  = MAX( ASO4,  0.d0 )

      RHD   = 0.80D+00                            ! RHD = 0.80 for ammonium sulfate (Ghan et al., 2001).
      RHC   = 0.35D+00                            ! RHC = 0.35 for ammonium sulfate (Ghan et al., 2001).

      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(A4,7E14.5)') 'END',ASO4,ANO3,ANH4,AH2O,GNH3,GHNO3,DUST
        WRITE(AUNIT1,'(A4,7F14.5)') 'RHD',RHD
      ENDIF 

      !-------------------------------------------------------------------------
      ! Get the sea salt-associated water (only).
      !
      ! A simple parameterization provided by E. Lewis is used.
      !-------------------------------------------------------------------------
      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(A4,3A12  )') '   ','SALT','SSH2O'           
        WRITE(AUNIT1,'(A4,3F12.5)') 'TOP' ,SALT
      ENDIF

      IF ( H .GT. 0.45D+00 ) THEN     ! ... then we are above the crystallization RH of NaCl
        SSH2O = SALT * ( SSH2OA + SSH2OB / ( 1.0D+00 - H ) )
      ELSE
        SSH2O = 0.0D+00
      ENDIF

      IF ( WRITE_LOG ) THEN
        WRITE(AUNIT1,'(A4,3F12.5)') 'END',SALT,SSH2O
      ENDIF


      END SUBROUTINE AERO_THERMO
