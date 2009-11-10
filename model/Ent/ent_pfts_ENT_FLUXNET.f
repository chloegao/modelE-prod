      module ent_pfts
!@sum Ent default supported 13 vegetation types

      !use ent_pftconst
      use ent_const
      use ent_types
      implicit none

      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************

      !* 1 - evergreen broadleaf early successional
      !* 2 - evergreen broadleaf late successional
      !* 3 - evergreen needleleaf early successional
      !* 4 - evergreen needleleaf late successional
      !* 5 - cold deciduous broadleaf early successional
      !* 6 - cold deciduous broadleaf late successional
      !* 7 - drought deciduous broadleaf
      !* 8 - deciduous needleleaf
      !* 9 - cold adapted shrub
      !* 10 - arid adapted shrub
      !* 11- C3 grass perennial
      !* 12 - C4 grass
      !* 13 - C3 grass - annual
      !* 14- arctic C3 grass
      !* 15- C4 crops
      !* 16 - crops broadleaf woody

      integer, parameter :: EVGRBROADEARLY = 1
      integer, parameter :: EVGRBROADLATE = 2
      integer, parameter :: EVGRNEEDLEEARLY = 3
      integer, parameter :: EVGRNEEDLELATE = 4
      integer, parameter :: COLDDECIDBROADEARLY = 5
      integer, parameter :: COLDDECIDBROADLATE = 6
      integer, parameter :: DROUGHTDECIDBROAD =7
      integer, parameter :: DECIDNEEDLE = 8
      integer, parameter :: COLDSHRUB = 9
      integer, parameter :: ARIDSHRUB = 10
      integer, parameter :: GRASSC3PER = 11
      integer, parameter :: GRASSC4 = 12
      integer, parameter :: GRASSC3 = 13
      integer, parameter :: GRASSC3ARCTIC = 14
      integer, parameter :: CROPSC4 = 15
      integer, parameter :: CROPSWOODY = 16
!##### TEMPORARY HACK - YK #####
!to avoid the conflict in ent_prescribed_drv.f90, using CROPS!
      integer, parameter :: CROPS = 15
!##### END OF TEMPORARY HACK #####
      
      ! other parameters needed for Ent to compile
      integer, parameter :: COVEROFFSET = 0

      !pst - photosynthetic pathway, 1-C3, 2-C4
      !woody - Woody, FALSE=NO, TRUE=YES.
      !leaftype - 1=broadleaf, 2=needleleaf, 3=monocot (not crops)
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !lrage - leaf and root litter age (years) (CASA)
      !woodage - stem litter age (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content
      !* Parameter for phenology
      !phenotype - phenological types
      !          - evergreen (1), 
      !            cold deciduous (2), 
      !            drought deciduous (3), 
      !            cold/drought deciduous (4),
      !            annual  (5)
      !annual - Annual, FALSE = NO, TRUE=YES.
      !* Parameters for plant allomteries 
      !* (Albani et al. Global Change Biology 2006 & 
      !* estimated from KM67 (Santarem, Amazon) tree survey.)
      !b1Cf - para 1 for allometric relation between DBH & foliage C 
      !b2Cf - para 2 for allometric relation between DBH & foliage C
      !b1Cd - para 1 for allometric relation between DBH & structural(dead) C
      !b2Cd - para 2 for allometric relation between DBH & structural(dead) C
      !b1Ht - para 1 for allometric relation between DBH & height
      !b2Ht - para 2 for allometric relation between DBH & height

      !SOURCES:
      !soil  moisture points:  7 & 10 from Kiang (2002) dissertation.
      !                       Other values are guesses.
      !nf:  Kull&Kruijt ps cap/leaf N param. Guesses from Friend&Kiang (2005) 
      !sla:  from LSM parameterizations.
      !lrage, woodage,lit_C2N,lignin: from CASA parameterizations.
      
      !***************************************************
      !Temp values for Ent pfts (See ent_const.f for types)
      type(pftype),parameter :: pfpar(N_PFT) =         !PFT parameters
      !pst, woody,leaftype, hwilt, sstar, swilt,nf,sla,r,
      !lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype, 
      !b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht
     &     (/      
     ! !* 1 - evergreen broadleaf early successional
     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.2d0, 8.8d0, 0.5d0, 
     &     1.8d0, 41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     0.0347d0, 1.560d0, 0.0816d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 2 - evergreen broadleaf late successional
!YK - KM67, Tapajo National Forest, Santarem, Brazil
!YK - lrage=3.0 & sla=9.7
!YK     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.1d0, 8.8d0, 0.5d0, 
!YK     &     1.8d0,41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.1d0, 9.7d0, 0.5d0, 
     &     3.0d0,41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     0.0395d0, 1.560d0, 0.1017d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 3 - evergreen needleleaf early successional
     &     pftype(1,.true.,2,-153.d0, .50d0, .25d0, 0.9d0, 9.5d0, 1.2d0, 
     &     4.0d0, 42.0d0, 80.d0, 0.25d0, 0.184d0, 1,
     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0),
     ! !* 4 - evergreen needleleaf late successional
     &     pftype(1,.true.,2,-153.d0, .50d0, .25d0, 0.85d0,9.5d0, 1.2d0, 
     &     4.0d0, 42.0d0, 80.d0, 0.25d0, 0.184d0, 1,
     &     0.0450d0, 1.683d0, 0.1617d0, 2.1536d0, 22.79d0, -0.0445d0),
     ! !* 5 - cold deciduous broadleaf early successional
     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.5d0, 34.5d0,0.6d0, !34.5
!     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.5d0, 12.4d0,0.6d0, 
     &     1.2d0, 58.0d0, 57.d0, 0.2d0, 0.093d0, 2,
     &     0.0240d0, 1.860d0, 0.1480d0, 2.411d0, 25.18d0, -0.0496d0),
     ! !* 6 - cold deciduous broadleaf late successional
     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.4d0, 11.5d0,0.6d0, 
     &     0.75d0, 58.0d0, 57.d0, 0.2d0, 0.093d0, 2,
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 23.39d0, -0.0540d0),
     ! !* 7 - drought deciduous broadleaf
!     &     pftype(1,.true.,1,-500.d0, .45d0, .22d0, 1.4d0, 8.3d0, 0.5d0,
!YKIM - Oaks at Tonzi
!YKIM - sstar = 0.34 & swilt 0.28
     &     pftype(1,.true.,1,-500.d0, .34d0, .28d0, 1.4d0, 8.3d0, 0.5d0, 
     &     1.2d0,25.0d0, 60.d0, 0.2d0, 0.153d0, 3,
     &     0.0296d0, 1.560d0, 0.0621d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 8 - deciduous needleleaf !## SLA from Reich (1997) leaf longev. 1 yr
     &     pftype(1,.true.,2,-100.d0,.55d0, .25d0, 0.9d0, 10.0d0, 0.9d0, 
     &     1.8d0, 27.0d0, 50.d0, 0.2d0, 0.2d0, 2,
     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0),
     ! !* 9 - cold adapted shrub
     &     pftype(1,.true.,1,-153.d0,.50d0, .30d0, 1.4d0, 2.25d0, 0.6d0, 
     &     2.8d0, 5.5d0, 50.d0, 0.15d0, 1.40d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 10 - arid adapted shrub
     &     pftype(1,.true.,1,-2030.d0,.40d0,.22d0, 1.3d0, 3.25d0, 0.6d0, 
     &     1.0d0, 5.5d0, 65.d0, 0.2d0, 0.32d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 11 - C3 grass perennial
     &     pftype(1,.false.,3,-2030.d0,.30d0,.10d0,1.5d0, 21.6d0, 1.2d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 12 - C4 grass
     &     pftype(2,.false.,3,-2030.d0,.30d0,.10d0, 1.3d0,9.4d0, 0.6d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 13 - C3 grass - annual 
     &     pftype(1,.false.,3,-2030.d0,.30d0,.1d0, 1.5d0, 21.6d0, 1.2d0, !10->15  
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 5,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 14 - arctic C3 grass
     &     pftype(1,.false.,3,-153.d0,.60d0, .27d0, 1.4d0, 9.0d0, 0.6d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 15 - C4 crops herbaceous
     &     pftype(2,.false., 1,-153.d0,.45d0,.27d0,1.3d0, 15.d0, 0.6d0, 
     &     1.1d0, UNDEF, 52.5d0, 0.16d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 16 - crops - broadleaf woody !## COPIED FROM BROAD COLDDECID LATE ##
     &     pftype(1,.true.,1,-153.d0,.50d0, .29d0, 1.4d0, 8.3d0, 0.9d0, 
     &     1.2d0, 58.0d0, 50.d0, 0.2d0, 0.093d0, 4,
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 23.39d0, -0.0540d0)
     &     /)

      end module ent_pfts
