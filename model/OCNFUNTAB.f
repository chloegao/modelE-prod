C****   
C**** FUNTABLE.OCN   Table functions of ocean parameters   10/09/91
C****
      MODULE OCFUNC
!@sum  OCFUNC contains the ocean function lookup tables
!@auth Gary Russell/Gavin Schmidt      
      IMPLICIT NONE
      SAVE
C**** All look up tables take specific pot. enthalpy (J/kg), 
C**** salinity (ppt) and pressure (MPa = 10^6 Pa) as arguments
!@var VGSP specific volume (m^3/kg)
!@var TGSP potential temperature (C)
!@var HGSP specific enthalpy of seawater (J/kg)
!@var AGSP thermal expansion coefficient (kg/m**3/C)
!@var BGSP saline expansion coefficient (kg/m**3/PSU)
!@var CGS specific heat capacity of sea water (J/kg*C)
      REAL*8, DIMENSION(-2:40,0:40,0:39) :: VGSP,TGSP,HGSP,AGSP,BGSP
      REAL*8, DIMENSION(-2:40,0:40) :: CGS
      END MODULE OCFUNC

      REAL*8 FUNCTION VOLGSP (G,S,P)
C****
C**** VOLGSP returns a linearly interpolated specific volume from
C**** an input table that depends on potential specific enthalpy,
C**** salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: VOLGSP (m**3/kg) = specific volume of sea water
C****
      USE OCFUNC, only : v=>vgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      VOLGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*V(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*V(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION VOLGS (G,S)
C****
C**** VOLGS returns a linearly interpolated specific volume from
C**** an input table that depends on potential specific enthalpy
C**** and salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: VOLGS (m**3/kg) = specific volume of sea water
C****
      USE OCFUNC, only : v=>vgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS
      INTEGER IG,JS
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      VOLGS = (JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,0)
     *                 + (GG-IG  )*V(IG+1,JS  ,0))
     *      + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,0)
     *                 + (GG-IG  )*V(IG+1,JS+1,0))
      RETURN
      END

      REAL*8 FUNCTION TEMGS (G,S)
C****
C**** TEMGS returns a linearly interpolated temperature from an
C**** input table that depends on potential specific enthalpy and
C**** salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure. 
C**** 
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: TEMGS (C) = temperature of sea water
C****
      USE OCFUNC, only : t=>tgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS
      INTEGER IG,JS
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      TEMGS = (JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,0)
     *                 + (GG-IG  )*T(IG+1,JS  ,0))
     *      + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,0)
     *                 + (GG-IG  )*T(IG+1,JS+1,0))
      RETURN
      END

      REAL*8 FUNCTION SHCGS (G,S)
C****
C**** SHCGS returns a linearly interpolated specific heat capacity
C**** from an input table that depends on potential specific enthalpy
C**** and salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: SHCGS (J/kg*C) = specific heat capacity of sea water
C****
      USE OCFUNC, only : c=>cgs
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      SHCGS = (JS-SS+1)*((IG-GG+1)*C(IG  ,JS  )
     *                 + (GG-IG  )*C(IG+1,JS  ))
     *      + (SS-JS  )*((IG-GG+1)*C(IG  ,JS+1)
     *                 + (GG-IG  )*C(IG+1,JS+1))
      RETURN
      END

      REAL*8 FUNCTION GFREZS (S)
C****
C**** GFREZS returns a linearly interpolated freezing point of
C**** potential specific enthalpy that depends on salinity.  Pressure
C**** is assumed to be the normal surface ocean reference pressure.
C****
C**** Input: S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: GFREZS (J/kg) = freezing point of potential specific
C****                         enthalpy
C****
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: S
      REAL*8 F(0:40)
      DATA F /                                                0.000,
     * -232.482d0,  -461.136d0,  -688.288d0,  -914.774d0, -1141.078d0,
     *-1367.530d0, -1594.374d0, -1821.798d0, -2049.957d0, -2278.978d0,
     *-2508.971d0, -2740.030d0, -2972.237d0, -3205.667d0, -3440.386d0,
     *-3676.454d0, -3913.926d0, -4152.851d0, -4393.287d0, -4635.259d0,
     *-4878.815d0, -5123.994d0, -5370.830d0, -5619.358d0, -5869.609d0,
     *-6121.615d0, -6375.402d0, -6631.001d0, -6888.436d0, -7147.733d0,
     *-7408.917d0, -7672.010d0, -7937.037d0, -8204.019d0, -8472.976d0,
     *-8743.931d0, -9016.917d0, -9291.927d0, -9568.992d0, -9848.131d0/
      REAL*8 SS
      INTEGER JS
C****
      SS = S*1000.
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      GFREZS = (JS-SS+1)*F(JS) + (SS-JS)*F(JS+1)
      RETURN
      END

      REAL*8 FUNCTION TFREZS (SIN)
C****
C**** TFREZS calculates the freezing temperature of sea water as a
C**** function of salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: SIN (1) = salinity (kg NaCl/kg sea water), from .004 to .04
C****
C**** Output: TFREZS (C) = freezing temperature of sea water
C****
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: SIN
      REAL*8 :: A01 = -.0575d0, A02 = -2.154996D-4, A03 =1.710523D-3
      REAL*8 S,S32
C****
      S   = SIN*1.D3
      S32 = S*DSQRT(S)
      TFREZS = (A01 + A02*S)*S + A03*S32
      RETURN
      END

      REAL*8 FUNCTION HETGSP (G,S,P)
C****
C**** HETGSP returns a linearly interpolated specific enthalpy from
C**** an input table that depends on potential specific enthalpy,
C**** salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: HETGSP (J/kg) = specific enthalpy of sea water
C****
      USE OCFUNC, only : h=>hgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      HETGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*H(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*H(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*H(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*H(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*H(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*H(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*H(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*H(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION ALPHAGSP (G,S,P)
C****
C**** ALPHAGSP returns a linearly interpolated thermal expansion
C**** coefficient from an input table that depends on potential
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: ALPHAGSP (kg/m**3/C) = thermal expansion coefficient 
C****
      USE OCFUNC, only : a=>agsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      ALPHAGSP=(KP-PP+1)*((JS-SS+1)*((IG-GG+1)*A(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*A(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*A(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*A(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*A(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*A(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*A(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*A(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION BETAGSP (G,S,P)
C****
C**** BETAGSP returns a linearly interpolated saline expansion
C**** coefficient from an input table that depends on potential 
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: BETAGSP (kg/m**3/PSU) = saline expansion coefficient 
C****
      USE OCFUNC, only : b=>bgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      BETAGSP =(KP-PP+1)*((JS-SS+1)*((IG-GG+1)*B(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*B(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*B(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*B(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*B(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*B(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*B(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*B(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION TEMGSP (G,S,P)
C****
C**** TEMGSP returns a linearly interpolated in situ temperature
C**** from an input table that depends on potential 
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: TEMGSP (C) = in situ temperature
C****
      USE OCFUNC, only : t=>tgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      TEMGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*T(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*T(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*T(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*T(IG+1,JS+1,KP+1)))
      RETURN
      END
