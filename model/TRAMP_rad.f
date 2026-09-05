#include "rundeck_opts.h"

c -----------------------------------------------------------------

      module AMP_Utilities_mod

      contains

      character*2 function aerosolkind(tracerName) 
      character(len=*), intent(in) :: tracerName
      aerosolkind = tracerName(7:8)
      end function aerosolkind

      end module AMP_Utilities_mod

      SUBROUTINE SETAMP(EXT,SCT,GCB,TAB)
!@sum Calculation of extinction, asymmetry and scattering for AMP Aerosols
!@sum Calculation of absorption in the longwave
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer 
      USE domain_decomp_atm,ONLY: am_i_root


#ifdef USE_OFFLINE_AEROSOLS
      USE OFFLINE_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV,
     +                       MIX_OC, MIX_SU, MIX_AQ, AMP_RAD_KEY,
     +                       NMODES
c       -> read this in                       MODE_NAME
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
c     AMP_TAB_SPEC (RADPAR) held one LW spectrum per mode and is no longer used here;
c     it is superseded by AMP_TAB_CLASS below, which resolves composition as well.
      USE RADPAR,      only: aesqex,aesqsc,aesqcb,FSTOPX,FTTOPX

      IMPLICIT NONE
      INTEGER, save:: Ifirstrad = 1
      ! Arguments: Optical Parameters dimension(lm,wavelength)
      REAL(8), INTENT(OUT) :: EXT(LM,6)       ! Extinction, SW
      REAL(8), INTENT(OUT) :: SCT(LM,6)       ! Single Scattering Albedo, SW
      REAL(8), INTENT(OUT) :: GCB(LM,6)       ! Asymmetry Factor, SW
      REAL(8), INTENT(OUT) :: TAB(LM,33)      ! Thermal absorption Cross section, LW
      REAL(8), DIMENSION(LM,NMODES) :: TTAUSV

      ! Local
      
      INTEGER l,n,w,s,MA,MB,MC,MD,NA,NS
      REAL*8 sizebins(23), Mie_IM(17), Mie_RE(15), HELP, AMP_TAB(33)
      REAL*8 CORE_CLASS(nmodes), SHELL_CLASS(nmodes),Reff_mode(nmodes),Vf(6),CS_Mix(26)
      REAL*8 a,b,AMPEXT,AMPSCA,AMPASY
c-----------------------------------------------------------------------------------------
c     Longwave absorption is mixed over the composition of each mode, rather than taken
c     from a single composition per mode (see the LW section below for why). AMP_TAB_CLASS
c     holds the LW absorption spectrum of each pure composition at each mode's effective
c     radius, precalculated once; the runtime mixing is a volume-fraction weighted sum
c     over these. NDRY is the number of dry species tracked in dry_Vf_LEV.
c
c     SPC_TO_NA maps the species order of VMass/dry_Vf_LEV (set in SETAMP_LEV) onto the
c     aerosol composition codes NA used by GET_LW. The two orderings differ, so this
c     mapping must not be bypassed:
c         dry_Vf_LEV slot :  1 SU   2 BC   3 OC   4 DU   5 SS   6 NO3
c         GET_LW NA code  :  1 SO4  2 SEA  3 NO3  4 OC   5 BC   6 DUST
c-----------------------------------------------------------------------------------------
      INTEGER, PARAMETER :: NDRY=6, NCLASS=6
      INTEGER SPC_TO_NA(NDRY)
      REAL*8, SAVE :: AMP_TAB_CLASS(33,NCLASS,nmodes)
      DATA SPC_TO_NA /1, 5, 4, 6, 2, 3/
      DATA sizebins/0.002, 0.005,0.01,0.05,0.08,0.1,0.13,0.17,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.5,2.,3.,5.,10./
      DATA CS_Mix/0.,0.04,0.08,0.12,0.16,0.2,0.24,0.28
     +          ,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6
     +          ,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,1.0/
      DATA Mie_RE/1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.9/
      DATA Mie_IM/0.0,0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0/
#ifdef USE_OFFLINE_AEROSOLS
           CHARACTER(LEN=3), SAVE :: MODE_NAME(16)
      DATA MODE_NAME(1:16)/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/
#endif
c     One table per mechanism, holding the modes listed in IMODES in AERO_CONFIG.
c     The per-mode values are common to all mechanisms; each table below just
c     selects the modes that mechanism actually carries. The only exception is
c     REFF_mode of SSS (mechanisms 4 and 8), which appears in no other mechanism
c     and is set here to DG_SSS from TRAMP_param_GISS.
#if (defined TRACERS_AMP_M1) || (defined USE_OFFLINE_AEROSOLS)
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  BC3  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   5,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
c      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   1,   2/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,  
     +                 0.100D+00, 0.100D+00, 0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/    
#elif defined TRACERS_AMP_M2
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  OCS  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   4,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,
     +                 0.100D+00, 0.075D+00, 0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/
#elif defined TRACERS_AMP_M3
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  BOC  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   4,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,
     +                 0.100D+00, 0.100D+00, 0.100D+00/
#elif defined TRACERS_AMP_M4
c                        ACC  DD1  DS1  DD2  DS2  SSS  OCC  BC1  BC2  MXX
c                        1    2    3    4    5    6    7    8    9    10
      DATA CORE_CLASS   /1,   6,   6,   6,   6,   2,   4,   5,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00, 2.00D+00 ,
     +                 1.380D+00, 0.075D+00, 0.050D+00, 0.100D+00, 0.100D+00/
#elif defined TRACERS_AMP_M5
c                        AKK  ACC  DD1  DS1  SSA  SSC  OCC  BC1  BC2  BC3  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14
      DATA CORE_CLASS   /1,   1,   6,   6,   2,   2,   4,   5,   5,   5,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 0.12D+00 ,
     +                 2.D+00   , 0.075D+00, 0.050D+00, 0.100D+00, 0.100D+00,
     +                 0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/
#elif defined TRACERS_AMP_M6
c                        AKK  ACC  DD1  DS1  SSA  SSC  OCC  BC1  BC2  OCS  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14
      DATA CORE_CLASS   /1,   1,   6,   6,   2,   2,   4,   5,   5,   4,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 0.12D+00 ,
     +                 2.D+00   , 0.075D+00, 0.050D+00, 0.100D+00, 0.075D+00,
     +                 0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/
#elif defined TRACERS_AMP_M7
c                        AKK  ACC  DD1  DS1  SSA  SSC  OCC  BC1  BC2  BOC  MXX
c                        1    2    3    4    5    6    7    8    9    10   11
      DATA CORE_CLASS   /1,   1,   6,   6,   2,   2,   4,   5,   5,   4,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 0.12D+00 ,
     +                 2.D+00   , 0.075D+00, 0.050D+00, 0.100D+00, 0.100D+00,
     +                 0.100D+00/
#elif defined TRACERS_AMP_M8
c                        ACC  DD1  DS1  SSS  OCC  BC1  BC2  MXX
c                        1    2    3    4    5    6    7    8
      DATA CORE_CLASS   /1,   6,   6,   2,   4,   5,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.075D+00, 1.160D+00, 2.000D+00, 1.380D+00, 0.075D+00,
     +                 0.050D+00, 0.100D+00, 0.100D+00/
#elif defined TRACERS_AMP_M9
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  OCS  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   4,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
c      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   1,   1,   2/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,  
     +                 0.100D+00, 0.075D+00, 0.100D+00, 0.070D+00, 0.100D+00/
#elif defined TRACERS_AMP_M10
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  OCS  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   4,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,
     +                 0.100D+00, 0.075D+00, 0.100D+00, 0.070D+00, 0.100D+00/
#endif
  
c                   NA1= SO4  NA2=SS  NA3=NO3 NA4=OC NA5=BC NA6=DU

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
c     commented out, and mixing is now handled by the volume weighting below, which
c     supersedes the CORE_CLASS/SHELL_CLASS pair. Both arrays are kept in the mechanism
c     tables above for reference but are no longer read.

      if ( Ifirstrad==1 ) then
      Ifirstrad = 0
      DO n = 1,nmodes
      DO NA = 1,NCLASS
      Vf(:)=0.d0
      CALL GET_LW(NA,0,Reff_mode(n),AMP_TAB,Vf)
      AMP_TAB_CLASS(:,NA,n)=AMP_TAB(:)
      enddo
      enddo
      endif

      if (itime.ne.itimeI) then 
          IF (AMP_RAD_KEY == 1 .or. AMP_RAD_KEY ==3) THEN
 
c Shortwave: ---------------------------------------------------------------------------------------------    

      DO l = 1,lm
      DO n = 1,nmodes

         w = 6    ! aot at 550
         do MD = 1,23
            if (Reff_LEV(l,n) .le. sizebins(md)) goto 100            
         enddo
 100      continue  
          if (MD.gt.1) then
            MD = min(23,MD)
          b = sizebins(md) - Reff_LEV(l,n)
          a = Reff_LEV(l,n) - sizebins(md-1)
          endif
c---- INTERNAL MIXTURE ---------------------------------------------        
            do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 200
            enddo
 200        continue
            do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 201
            enddo
 201        continue
            MA = min(15,MA)              
            MB = min(17,MB)
            if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55(MA,MB,MD) + b/(a+b) * AMP_Q55(MA,MB,MD-1))
            else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55(MA,MB,MD)
            endif
C----------------------------------------------------------------------
      DO w = 1,6  !wavelength
c---- INTERNAL MIXTURE ---------------------------------------------        
         do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 401
         enddo
 401     continue
         do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 402
         enddo
 402     continue
          MA = min(15,MA)
          MB = min(17,MB)

          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT(MA,MB,MD,w)  + b/(a+b) * AMP_EXT(MA,MB,MD-1,w))
          AMPSCA = (a/(b+a)* AMP_SCA(MA,MB,MD,w)  + b/(a+b) * AMP_SCA(MA,MB,MD-1,w))
          AMPASY = (a/(b+a)* AMP_ASY(MA,MB,MD,w)  + b/(a+b) * AMP_ASY(MA,MB,MD-1,w))
          else
          AMPEXT = AMP_EXT(MA,MB,MD,w) 
          AMPSCA = AMP_SCA(MA,MB,MD,w) 
          AMPASY = AMP_ASY(MA,MB,MD,w) 
          endif
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

      DO l = 1,lm
      DO n = 1,nmodes

         w = 6    ! aot at 550
         do MD = 1,23
            if (Reff_LEV(l,n) .le. sizebins(md)) goto 500            
         enddo
 500      continue  
          if (MD.gt.1) then
            MD = min(23,MD)
          b = sizebins(md) - Reff_LEV(l,n)
          a = Reff_LEV(l,n) - sizebins(md-1)
          endif

       select case (MODE_NAME(n))
c---- INTERNAL MIXTURE ---------------------------------------------        
       case ('AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','OCS','DBC','MXX')
   
            do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 600
            enddo
 600        continue
            do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 601
            enddo
 601        continue
            MA = min(15,MA)              
            MB = min(17,MB)
            if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55(MA,MB,MD) + b/(a+b) * AMP_Q55(MA,MB,MD-1))
            else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55(MA,MB,MD)
            endif

C------ CORE SHELL -------------------------------------------------
       case ('BC1','BC2','BC3','BOC','BCS')
   
         do MA = 1,26
            if (MIX_OC(l,n) .le. CS_MIX(MA)) goto 701
         enddo
 701     continue
         do MB = 1,26
            if (MIX_SU(l,n) .le. CS_MIX(MB)) goto 702
         enddo
 702     continue
         do MC = 1,26
            if (MIX_AQ(l,n) .le. CS_MIX(MC)) goto 703
         enddo
 703     continue
         
          MA = min(26,MA-1)
          MB = min(26,MB-1)
          MC = min(26,MC-1)
          MA = max(1,MA-1)
          MB = max(1,MB-1)
          MC = max(1,MC-1)

         if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55_CS(MD,MA,MB,MC) + b/(a+b) * AMP_Q55_CS(MD-1,MA,MB,MC))
         else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55_CS(MD,MA,MB,MC)
         endif
       end select
C----------------------------------------------------------------------
      DO w = 1,6  !wavelength
c---- INTERNAL MIXTURE ---------------------------------------------        

       select case (MODE_NAME(n))

       case ('AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','OCS','DBC','MXX')
   
         do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 801
         enddo
 801     continue
         do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 802
         enddo
 802     continue
          MA = min(15,MA)
          MB = min(17,MB)

          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT(MA,MB,MD,w)  + b/(a+b) * AMP_EXT(MA,MB,MD-1,w))
          AMPSCA = (a/(b+a)* AMP_SCA(MA,MB,MD,w)  + b/(a+b) * AMP_SCA(MA,MB,MD-1,w))
          AMPASY = (a/(b+a)* AMP_ASY(MA,MB,MD,w)  + b/(a+b) * AMP_ASY(MA,MB,MD-1,w))
          else
          AMPEXT = AMP_EXT(MA,MB,MD,w) 
          AMPSCA = AMP_SCA(MA,MB,MD,w) 
          AMPASY = AMP_ASY(MA,MB,MD,w) 
          endif
C------ CORE SHELL -------------------------------------------------
       case ('BC1','BC2','BC3','BOC','BCS')
   
         do MA = 1,26
            if (MIX_OC(l,n) .le. CS_MIX(MA)) goto 805
         enddo
 805     continue
         do MB = 1,26
            if (MIX_SU(l,n) .le. CS_MIX(MB)) goto 806
         enddo
 806     continue
         do MC = 1,26
            if (MIX_AQ(l,n) .le. CS_MIX(MC)) goto 807
         enddo
 807     continue

          MA = min(26,MA-1)
          MB = min(26,MB-1)
          MC = min(26,MC-1)
          MA = max(1,MA-1)
          MB = max(1,MB-1)
          MC = max(1,MC-1)
                    
          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_EXT_CS(MD-1,MA,MB,MC,w))
          AMPSCA = (a/(b+a)* AMP_SCA_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_SCA_CS(MD-1,MA,MB,MC,w))
          AMPASY = (a/(b+a)* AMP_ASY_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_ASY_CS(MD-1,MA,MB,MC,w))
          else
          AMPEXT = AMP_EXT_CS(MD,MA,MB,MC,w)
          AMPSCA = AMP_SCA_CS(MD,MA,MB,MC,w)
          AMPASY = AMP_ASY_CS(MD,MA,MB,MC,w)
          endif
        end select
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
      DO n = 1,nmodes
         AMP_TAB(:) = 0.d0
         DO s = 1,NDRY
            AMP_TAB(:) = AMP_TAB(:)
     +                 + dry_Vf_LEV(l,n,s) * AMP_TAB_CLASS(:,SPC_TO_NA(s),n)
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
      REAL*8,     DIMENSION(nmodes,7) :: VolFrac, VMass
      REAL*8                          :: H2O, NO3 
      REAL(8), PARAMETER :: TINYNUMER = 1.0D-30 
      COMPLEX*8, DIMENSION(6,7)      :: Ri
c     Variables for Maxwell Garnett:
      REAL*8                           :: V_bc, V_host
      COMPLEX*8                        :: M_mg, M_bc, M_host
c Andies data incl Solar weighting - integral over 6 radiation band
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
      REAL(8) :: NI(16)           ! number concentration for each tracer [#/m^3]
      REAL(8) :: M(16,5)          ! mass   concentration for each tracer [ug/m^3]
      REAL(8) :: DG_WET(16)       ! geometric mean diameter for each dry tracer [um]
      REAL(8) :: DGN(16)          ! geometric mean diameter for each dry tracer [um]
      REAL(8) :: TK               ! absolute temperature [K]
      REAL(8) :: PRES             ! ambient pressure [Pa]
      REAL(8) :: AIRD(lm)          ! air density [kg/m^3]
      REAL(8) :: VOLTMP_WET 
      real*8, parameter, dimension(16) ::
     &    sig0=(/ 1.6d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0, 
     &                 2.0d0, 2.0d0, 1.8d0, 1.8d0, 1.8d0, 1.8d0, 
     &                 1.8d0, 1.8d0, 1.8d0, 2.0d0/)
      real*8, parameter, dimension(16) ::
     &    CONV_DPAM_TO_DGN=(/ 0.71795016727196403,   0.59556797724590516     ,  
     &     0.59556797724590516     ,  0.59556797724590516   ,    0.59556797724590516  ,     
     &     0.59556797724590516      , 0.59556797724590516   ,    0.48642160999311468  ,    
     &     0.59556797724590516   ,    0.59556797724590516  ,     0.59556797724590516  ,   
     &     0.59556797724590516  ,     0.59556797724590516    ,   0.59556797724590516  ,     
     &     0.59556797724590516 ,      0.48642160999311468/)

      real*8, parameter, dimension(5) ::
     &    dens=(/1.77D+03,1.70D+03,1.00D+03,2.60D+03,2.165D+03/)   ! [kg/m^3] 
           CHARACTER(LEN=3), SAVE :: MODE_NAME(16)
      DATA MODE_NAME(1:16)/'AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','BC1','BC2','BC3','DBC','BOC','BCS','MXX'/
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
      do n=1,nmodes ! loop over modes
           VOLTMP_WET = 1.0D-30
           do k = 1, 5
             M(n,k) = actvqtys(l,n,k,i,j) * MA(l,i,j)               ! mass [kg/kg]  -> [kg/m2/layer]
             VOLTMP_WET = VOLTMP_WET + actvqtys(l,n,k,i,j)/dens(k)  ! dry volume [m3]
           enddo

             VOLTMP_WET = VOLTMP_WET + (M(n,1)/sum(M(:,1))) * actvqtys2(l,1,i,j)/1720.d0  ! NO3
             VOLTMP_WET = VOLTMP_WET + (M(n,1)/sum(M(:,1))) * actvqtys2(l,2,i,j)/1720.d0  ! NH4
             VOLTMP_WET = VOLTMP_WET + (M(n,1)/sum(M(:,1))) * actvqtys2(l,3,i,j)/1000.d0  ! H2O

         DG_WET(n) = 1.D6 *( (6.d0/pi) * (VOLTMP_WET / actvqtys(l,n,6,i,j)) ) ! dry gemoetric mean mass diameter [um]
     &                 **0.333333333333333  ! [um]
         DG_WET(n) = MIN( MAX( DG_WET(n),  0.01), 10.D0 )

         DGN(n) = DG_WET(n) * (1.0D+00 / EXP( 1.5D+00 * ( log(sig0(n)) )**2 )) 
     
c         Reff_LEV(l,n) = DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5      
         Reff_LEV(l,n) = actvqtys(l,n,7,i,j)*CONV_DPAM_TO_DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5e6
       
          VMass(n,1:5)  = M(n,1:5)/DENS(1:5)           

       enddo



       ! + Volume Fraction
       DO n=1,nmodes  ![#/m2]         pi/4     [m2]
         NI(n)  =  actvqtys(l,n,6,i,j) * MA(l,i,j)           ! number [#/kg] *  [kg/m2] = [#/m2]
c         NUMB_LEV(l,n) = NI(n)* 0.7853 * (1.e-6*DG_WET(n))**2   ! [#/layer]
         NUMB_LEV(l,n) = NI(n)* 0.7853 *actvqtys(l,n,7,i,j)**2
         
        ! NO3   
        VMass(n,6) = VMass(n,1) / Sum(VMass(:,1)) * ((actvqtys2(l,1,i,j) + actvqtys2(l,2,i,j))* MA(l,i,j) )/ 1720.
        ! H2O
        VMass(n,7) = VMass(n,1) / Sum(VMass(:,1))  *(actvqtys2(l,3,i,j) * MA(l,i,j) )/1000.
       ENDDO

#else   /* Original MATRIX online code below */

       ! + Effective Radius [um] per Mode = geometric mass mean radius
       DO n=1,nmodes
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
                  VMass(AMP_MODES_MAP(nAMP),1) =trm(i,j,l,n)/DENS_SULF
               case ('BC')
                  VMass(AMP_MODES_MAP(nAMP),2) =trm(i,j,l,n)/DENS_BCAR
               case ('OC')
                  VMass(AMP_MODES_MAP(nAMP),3) =VMass(AMP_MODES_MAP(nAMP),3)+ trm(i,j,l,n)/DENS_OCAR
               case ('DU')
                  VMass(AMP_MODES_MAP(nAMP),4) =trm(i,j,l,n)/DENS_DUST
               case ('SS')
                  VMass(AMP_MODES_MAP(nAMP),5) =trm(i,j,l,n)/DENS_SEAS
               end select
             endif
           else                           ! Number
!          [ - ]                        [trm units: #/m2/layer]      
           NUMB_LEV(l,AMP_NUMB_MAP(nAMP)) =trm(i,j,l,n)
          endif
       ENDDO

       ! + Volume Fraction
       DO n=1,nmodes  ![#/m2]         pi/4     [m2]
        NUMB_LEV(l,n) = NUMB_LEV(l,n)* 0.7853 * DIAM(i,j,l,n)**2
        ! NO3   
        VMass(n,6) = VMass(n,1) / Sum(VMass(:,1)) * NO3/ 1720.
        ! H2O
        VMass(n,7) = VMass(n,1) /(Sum(VMass(:,1)) + TINYNUMER)  * H2O /1000.
       ENDDO
#endif   /* After that code should work for all cases */

      DO s=1,7  ! loop over species
        DO n=1,nmodes           ! loop over modes
          Volfrac(n,s) = VMass(n,s) / (Sum(VMass(n,:)) + TINYNUMER)
          dry_Vf_LEV(l,n,s) = VMass(n,s) / (Sum(VMass(n,1:6)) + TINYNUMER)
        ENDDO
      ENDDO

      ! Core Shell Composition: core is BC, shell material is OC, SO4 and H2O.
      ! The modes selected here must match the core-shell modes in SETAMP.
      DO n=1,nmodes             ! loop over modes
        select case (MODE_NAME(n))
        case ('BC1','BC2','BC3','BOC','BCS')
          MIX_OC(l,n) = VMass(n,3) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
          MIX_SU(l,n) = VMass(n,1) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
          MIX_AQ(l,n) = VMass(n,7) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
        end select
      ENDDO
 
      ! + Refractive Index of Aerosol mix per mode and wavelength
      
      RindexAMP(l,:,:) = 0.d0
      DO s=1,7                  ! loop over species 
        DO w=1,6                ! loop over wavelength
          DO n=1,nmodes         ! loop over modes
            RindexAMP(l,n,w) = RindexAMP(l,n,w) + ( Volfrac(n,s) * Ri(w,s))
          ENDDO
        ENDDO
      ENDDO

      if (AMP_RAD_KEY == 3) then      ! - - - Maxwell Garnett Mixing Rule
        DO w=1,6                ! loop over wavelength
          DO n=1,nmodes         ! loop over modes
       select case (MODE_NAME(n))
       case ('BC1','BC2','BC3','BOC','BCS')
             M_bc   = Ri(w,2)
             M_host = ( Volfrac(n,1) * Ri(w,1))
             DO s=3,7                  ! loop over species other that BC
             M_host = M_host + ( Volfrac(n,s) * Ri(w,s))
             ENDDO
             V_bc   = Volfrac(n,2)
             V_host = Volfrac(n,1)+Volfrac(n,3)+Volfrac(n,4)+Volfrac(n,5)+Volfrac(n,6)+Volfrac(n,7)
             M_mg = M_host**2  * (M_bc**2 + 2.d0 * M_host**2 + 2.d0 * V_bc * (M_bc    - M_host   ) ) 
     +                         / (M_bc**2 + 2.d0 * M_host**2 -        V_host*(M_bc**2 - M_host**2) )

             RindexAMP(l,n,w) = SQRT( M_mg)
       end select      
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

#ifdef USE_OFFLINE_AEROSOLS
      USE OFFLINE_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       AMP_Q55
#else
      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA, AMP_Q55,
     +           AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS  
#endif	  
	  IMPLICIT NONE
      include 'netcdf.inc'
      integer start(4),count(4),count3(3),status
      integer start2(5),count2(5),count32(4)
      integer ncid, id1, id2, id3, id4,ncid2

      real*4, DIMENSION(15,17,23,6) ::  ASY,SCA,EXT
      real*4, DIMENSION(15,17,23) ::    QEX
      real*4, DIMENSION(23,26,26,26,6) ::  CS_ASY,CS_SCA,CS_EXT
      real*4, DIMENSION(23,26,26,26) ::    CS_QEX
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
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=1

          count(1)=15
          count(2)=17
          count(3)=23
          count(4)=6
          count3(1)=15
          count3(2)=17
          count3(3)=23
 

          status=NF_GET_VARA_REAL(ncid,id1,start,count,ASY)
          status=NF_GET_VARA_REAL(ncid,id2,start,count,EXT)
          status=NF_GET_VARA_REAL(ncid,id3,start,count,SCA)
          status=NF_GET_VARA_REAL(ncid,id4,start,count3,QEX)

          status=NF_CLOSE('AMP_MIE_TABLES',NCNOWRIT,ncid)

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
          start2(1)=1
          start2(2)=1
          start2(3)=1
          start2(4)=1
          start2(5)=1

          count2(1)=23
          count2(2)=26
          count2(3)=26 
          count2(4)=26
          count2(5)=6
          count32(1)=23
          count32(2)=26
          count32(3)=26
          count32(4)=26
 

          status=NF_GET_VARA_REAL(ncid2,id1,start2,count2,CS_ASY)
          status=NF_GET_VARA_REAL(ncid2,id2,start2,count2,CS_EXT)
          status=NF_GET_VARA_REAL(ncid2,id3,start2,count2,CS_SCA)
          status=NF_GET_VARA_REAL(ncid2,id4,start2,count32,CS_QEX)

          status=NF_CLOSE('AMP_CORESHELL_TABLES',NCNOWRIT,ncid2)

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

      USE RADPAR, only: TRUQEX, TRSQEX, TRDQEX, TRUQSC, TRSQSC, TRDQSC
     *                  , REFU22, REFS25, REFD25
      INTEGER, intent(IN) :: NA,NS
      REAL*8,  intent(in) :: areff,Vf(6)
      REAL*8   TQEX(33),TQSC(33),TQAB(33),TQEX_S(33),TQSC_S(33)
      REAL*8   QXAERN(25),QSAERN(25)
      REAL*8   wts,wta
      INTEGER  n0,k,n,nn

c CORE  
       IF(NA==0) THEN
         TQAB(:)= 0d0
       ENDIF
                                      !                               1   2   3   4
      IF(NA > 0 .and. NA < 5) THEN    !    NA : Aerosol compositions SO4,SEA,NO3,OC
        N0=0
        IF(NA==2) N0=22
        IF(NA==3) N0=44
        IF(NA==4) N0=88
        DO K=1,33
        DO N=1,22
        NN=N0+N
        QXAERN(N)=TRUQEX(K,NN)
        QSAERN(N)=TRUQSC(K,NN)
        ENDDO
        CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF
                              
      IF(NA==5) THEN                  !   NA : Aerosol compositions BC
        DO K=1,33
        QXAERN(:)=TRSQEX(K,:)    ! 1:25
        QSAERN(:)=TRSQSC(K,:)    ! 1:25
        CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF

                                      !                              6
      IF(NA==6) THEN                  !   NA : Aerosol composition DST
       DO K=1,33
        QXAERN(:)=TRDQEX(K,:)    ! 1:25
        QSAERN(:)=TRDQSC(K,:)    ! 1:25
        CALL SPLINE(REFD25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
      ENDDO
      ENDIF

c SHELL
c         IF(NS > 0 .and. NS < 5) THEN    !    NS : Aerosol compositions SO4,SEA,NO3,OC
c         N0=0
c         IF(NS==2) N0=22
c         IF(NS==3) N0=44
c         IF(NS==4) N0=88


c         DO K=1,33
c         DO N=1,22
c         NN=N0+N
c         IF (NS==1) WTS=Vf(1)                      ! <- shell fraction of aerosol composition
c         IF (NS==2) WTS=Vf(5)                      ! <- shell fraction of aerosol composition
c         WTA=1.D0-WTS
c         QXAERN(N)=TRUQEX(K,NN)
c         QSAERN(N)=TRUQSC(K,NN)
c         ENDDO
c         CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX_S(K),1.D0,1.D0,1)
c         CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC_S(K),1.D0,1.D0,1)
c         TQAB(K)=(TQEX(K)*WTA + TQEX_S(K)*WTS)-(TQSC(K)*WTA + TQSC_S(K)*WTS)
c         ENDDO

c      ENDIF

      RETURN
      END SUBROUTINE GET_LW
c -----------------------------------------------------------------

