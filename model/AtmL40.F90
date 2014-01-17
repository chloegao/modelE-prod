module VerticalRes
!@sum Vertical Resolution file, 40 layers, top at .1 mb, no GWDRAG
!@auth Original Development Team
  Implicit None
!@var LM    = number of dynamical layers
!@var LS1   = lowest layer of strtosphere
  Integer*4,Parameter :: LM=40, LS1=24

!@var MDRYA = dry atmospheric mass (kg/m^2) = 100*PSF/GRAV
!@var MTOP  = mass above dynamical top (kg/m^2) = 100*PMTOP/GRAV
!@var MFIXs = summation of MFIX (kg/m^2) = 100*(PTOP-PMTOP)/GRAV
!@var MVAR  = spatially and temporally varying column mass (kg/m^2)
!@var MSURF = MTOP + MFIXs + MVAR (kg/m^2)
!@var AM(L) = MFIX(L) + MVAR*MFRAC(L) (kg/m^2)
!@var MFIX(L)  = fixed mass in each layer (kg/m^2) = 100*[PLBOT(L)-PLBOT(L+1)]/GRAV
!@var MFRAC(L) = fraction of variable mass in each layer = DSIG(L)
  Real*8,Parameter :: MDRYA = 98400/9.80665d0, MTOP = 10/9.80665d0, MFIXs = 14990/9.80665d0, &
    MFIX(LM) = (/ 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0,  &
    0d0,0d0,0d0, &
    2200/9.80665d0, 2000/9.80665d0, 1800/9.80665d0, 1700/9.80665d0, 1600/9.80665d0, &
    1400/9.80665d0, 1200/9.80665d0, 1100/9.80665d0, 1000/9.80665d0,  438/9.80665d0, &
    246/9.80665d0,  138/9.80665d0,   78/9.80665d0, 43.8d0/9.80665d0, 24.6d0/9.80665d0, &
    13.8d0/9.80665d0,7.8d0/9.80665d0 /), &
    MFRAC(LM) = (/ 20/834d0, 22/834d0, 25/834d0, 27/834d0, 30/834d0, &
    35/834d0, 40/834d0, 45/834d0, 48/834d0, 50/834d0, &
    51/834d0, 52/834d0, 50/834d0, 48/834d0, 45/834d0, &
    42/834d0, 38/834d0, 34/834d0, 31/834d0, 28/834d0, &
    26/834d0, 24/834d0, 23/834d0, &   
    0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0, 0d0,0d0,0d0,0d0,0d0 /)

!**** Vertival resolution
!****                         ---MSURF=10034.0---    ---MSURF=5781.8----
!**** Layer   MFIX   MFRAC    MVAR     AM     SUM    MVAR     AM     SUM
!**** =====   ====   =====    ====     ==     ===    ====     ==     ===
!****   Top                                   1.0                    1.0
!****    40     .8    0        0.0     .8     1.8     0.0     .8     1.8
!****    39    1.4    0        0.0    1.4     3.2     0.0    1.4     3.2
!****    38    2.5    0        0.0    2.5     5.7     0.0    2.5     5.7
!****    37    4.5    0        0.0    4.5    10.2     0.0    4.5    10.2
!****    36    8.0    0        0.0    8.0    18.2     0.0    8.0    18.2
!****    35   14.1    0        0.0   14.1    32.2     0.0   14.1    32.2
!****
!****    25  203.9    0        0.0  203.9  1305.2     0.0  203.9  1305.2
!****    24  224.3    0        0.0  224.3  1529.6     0.0  224.3  1529.6
!****    23    0.0  23/834   234.5  234.5  1764.1   117.3  117.3  1646.8
!****    22    0.0  24/834   244.7  244.7  2008.8   122.4  122.4  1769.2
!****
!****     3    0.0  25/834   254.9  254.9  9605.7   127.5  127.5  5567.7
!****     2    0.0  22/834   224.3  224.3  9830.1   112.2  112.2  5679.8
!****     1    0.0  20/834   203.9  203.9 10034.0   102.0  102.0  5781.8
!****         ----  ------   -----  -----           -----  -----
!****       1528.6 834/834  8504.4 10033.0         4252.2 5780.8  

!@var PSF,PMTOP global mean surface, model top pressure  (mb)
!@var PTOP pressure at interface level sigma/const press coord syst (mb)
!@var PSFMPT,PSTRAT pressure due to troposhere,stratosphere
!@var PLbot pressure levels at bottom of layers (mb)
  Real*8,Parameter :: PSF=984.d0, PTOP = 150.d0, PMTOP = .1d0, PSFMPT = PSF-PTOP, &
    PSTRAT = PTOP-PMTOP, &
    PLBOT(1:LM+1) = (/ PSF,   964d0, 942d0, 917d0, 890d0, 860d0, 825d0, &  !  L=1,..   
    785d0, 740d0, 692d0, 642d0, 591d0, 539d0, 489d0, &  !  L=...    
    441d0, 396d0, 354d0, 316d0, 282d0, 251d0, 223d0, &
    197d0, 173d0,                                    &
    PTOP,                                            &  !  L=LS1    
    128d0, 108d0,  90d0,  73d0,  57d0,  43d0,  31d0, &  !  L=...    
    20d0,  10d0,5.62d0,3.16d0,1.78d0,  1.d0,        &
    .562d0,.316d0,.178d0, PMTOP /)                       !  L=..,LM+1

#ifndef SCM
!**** KEP depends on whether stratos. EP flux diagnostics are calculated
!**** If dummy EPFLUX is used set KEP=0, otherwise KEP=21
!@param KEP number of lat/height E-P flux diagnostics
  Integer*4,Parameter :: KEP = 0

!**** Based on model top, determine how much of stratosphere is resolved
!**** ISTRAT = 2:          PMTOP <   1 mb
!**** ISTRAT = 1:  1 mb <= PMTOP <  10 mb
!**** ISTRAT = 0: 10 mb <= PMTOP
  Integer*4,Parameter :: ISTRAT = 2
#endif
End Module VerticalRes

Subroutine DUMMY_STRAT
!**** Dummy routines in place of STRATDYN
!**** The vertical resolution also determines whether stratospheric wave drag will be applied or not.
!**** Hence also included here are some dummy routines for non-stratospheric models.
!@sum DUMMY dummy routines for non-stratospheric models
  Entry INIT_GWDRAG
  Entry GWDRAG
  Entry VDIFF
  Entry io_strat
  Entry ALLOC_STRAT_COM
!**** Dummy routines in place of STRAT_DIAG (EP flux calculations)
!**** Note that KEP=0 is set to zero above for the dummy versions.
  Entry EPFLUX
  Entry EPFLXI
  Entry EPFLXP
  Return
End Subroutine DUMMY_STRAT

