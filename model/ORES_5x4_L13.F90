!****   
!@sum  RES_5x4_L13 Resolution info for 2x2.5 horizontal
!!       and 13 vertical ocean layer
!@auth Gary Russel/Larissa Nazarenko
!****   
!**** RES_5x4_L13.f    2008/08/08
!****
      Module OCEANRES
      Implicit None
!****
!**** Vertical resolution of ocean
!****
      Integer*4,Parameter :: &
           IMO = 72, &
           JMO = 46, &
           LMO = 13, &!  maximum number of ocean layers in a column
       LMO_MIN =  2, &!  minimum number of ocean layers in a column
         LSRPD =  3, &!  deepest layer for penetrating solar radiation
         MAXGL =  5, &!  maximum no. layers for depositing glac. melt
        NOCEAN =  1  &!  NDYNO must be multiple of 2*NOCEAN
      , NORDER=4     !  order of Alternating Binomial Filter (must be even)
!*** 
!***
      Real*8,Parameter :: &
         AKHMIN = 1.5d8   & ! minimum horizontal diffusion
       , AKHFAC = 1d0     & ! factor to multiply built-in scaling for k_diff
       , OABFUX=0.        & ! coef. for divergence filter in X dir. 
       , OABFVX=0.        & ! coef. for vorticity filter in X dir.
       , by4tonv=0.       & ! coef. for divergence filter in Y dir.
       , by4tonu=0.         ! coef. for vorticity filter in Y dir.
      Real*8,Parameter :: &
          dZO(LMO)  =     & !  approximate thicknesses (m) of ocean layers
                      (/ 12d0, 18d0, 27d0, 40.5d0, 60.75d0, 91.125d0, &
                        136.6875d0, 205.03125d0, 307.546875d0,        & 
                        461.3203125d0, 691.98046875d0,                &
                        1037.970703125d0, 1556.9560546875d0 /) 
!****
!**** Vertical resolution of ocean
!****
!****  L  dZO    ZE   ZOC      L  dZO    ZE   ZOC
!****  =  ===    ==   ===      =  ===    ==   ===
!****  1   12    12     6      7  137   386   318
!****  2   18    30    21      8  205   591   488
!****  3   27    57    44      9  308   899   745
!****  4   41    98    77     10  461  1360  1129
!****  5   61   158   128     11  692  2052  1706
!****  6   91   249   204     12 1038  3090  2571
!****                         13 1557  4647  3868
!****
      EndModule OCEANRES
