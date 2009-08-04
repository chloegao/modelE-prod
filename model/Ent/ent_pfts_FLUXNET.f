      module ent_pfts
!@sum GISS vegetation types for Ent model.

      use ent_const
      use ent_types
      implicit none


      !***************************************************
      !*            GISS VEGETATION TYPES                *
      !***************************************************
      !Types: These are GISS GCM types until we get full data sets for Ent.
      !1-tundra, 2-grassland, 3-shrubland, 4-savanna, 5-deciduous forest,
      !6-evergreen needleleaf forest, 7-tropical rainforest, 8-crops
      !**** Correspondence to CASA/LSM vegetation types for sla from 
      !**** Dickinson et al. (J. Clim., Nov. 1998).  See casa_bgfluxes.f and
      !**** 
      !#-GISS   = Dickinson = CASA/LSM
      !1-tundra = tundra and semidesert
      !2-grassland = avg short grass & tall grass = LSM cool C3, warm C4 grasses
      !3-shrubland = evergreen shrubs, deciduous shrubs
      !4-savanna = shrubs & broadleaf evergreen & grasses (avg 25,35,40)
      !5-deciduous forest = deciduous broadleaf & decid needleleaf
      !6-evergreen needleleaf = needleleaf evergreen
      !7-tropical rainforest = deciduous broadleaf = LSM tropical seasonal tree 
      !8-crops = crops

      integer, parameter :: COVEROFFSET = 1 !SAND in first position in GISS array
      integer, parameter :: TUNDRA = 1
      integer, parameter :: GRASSC3 = 2
      integer, parameter :: SHRUB = 3
      integer, parameter :: SAVANNA = 4
      integer, parameter :: DECIDFOREST = 5
      integer, parameter :: EVERGRNEEDLE = 6
      integer, parameter :: TROPRAINF = 7
      integer, parameter :: CROPS = 8
!##### TEMPORARY HACK - YK #####
!to avoid the conflict in phenology.f
      integer, parameter :: DROUGHTDECIDBROAD = 20
!##### TEMPORARY HACK - YK #####
      !*-----------------------------------------
      !* Veg types correspondence between models:
*       LSM:   1  2    3  4  5  6  7  8  9 10 11 12 13 14^M
*      CASA:   4  5    1  2  6  7  9 11 10 10 12 12  6  8^M
*      GISS:   6  5,6  7  5  4  2  3 x  1  1  8  8   4  x
*     NOTE:  Need to distinguish: - needleleaf decid vs. evergreen
*                                 - bare soil
*                                 - savanna mix
*            Check for:  sla, lrage, woodage, albedo
      !*-----------------------------------------
      !For dependencies to work, these constants must go in ent_const.f
!      integer,parameter :: N_PFT = 8
!      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
!      integer,parameter :: N_OTHER = 1
!      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER
       !pst - photosynthetic pathway, 1-C3, 2-C4
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !r - CASA respiration parameter (gC/gN)(CASA)
      !lrage - leaf and root litter turnover time (years) (CASA)
      !woodage - stem litter turnover time (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content
      !*for GISS pft, 
      !*phenotype, b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht are dummies!
      !phenotype - phenological type 
      !* Parameters for plant allomteries 
      !* (Albani et al. Global Change Biology 2006 & 
      !* estimated from KM67 (Santarem, Amazon) tree survey.)
      !b1Cf - para 1 for allometric relation between DBH & foliage C 
      !b2Cf - para 2 for allometric relation between DBH & foliage C
      !b1Cd - para 1 for allometric relation between DBH & structural(dead) C
      !b2Cd - para 2 for allometric relation between DBH & structural(dead) C
      !b1Ht - para 1 for allometric relation between DBH & height
      !b2Ht - para 2 for allometric relation between DBH & height

      type(pftype),parameter :: pfpar(N_PFT) =          !PFT parameters
     &!        pst,woody,leaftype,hwilt,sstar,swilt,nf,sla,
     &!        r,lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype ! 
     &!        b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht
     &     (/
     &     pftype(1,.true., 1,-153.d0,  .50d0, .30d0,  1.4d0, !tundra
     &     2.25d0, 0.6d0, 2.8d0, 5.5d0, 50.0d0,0.15d0,1.4d0,2,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(2,.false., 3,-2030.d0,  .30d0, .10d0,  1.5d0, !grassC3 !hwilt=Vaira grassland final senescence soilmp layer2
     &     21.6d0, 1.2d0, 1.5d0, UNDEF, 50.0d0, 0.1d0,0.d0,4, ! Ponca nf and SLA
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(2,.true., 1,-153.d0,  .40d0, .22d0,  1.3d0, !shrub
     &     3.25d0, 0.6d0, 1.0d0, 5.5d0, 57.5d0, 0.15d0,0.32d0,3,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(2,.true., 1,-2030.d0,  .65d0, .22d0,  1.3d0, !savanna
     &     5.1d0, 1.d0, 1.8d0, 25.d0, 50.0d0, 0.15d0,0.153d0,3,
     &     0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0),
     &     pftype(1,.true., 1,-500.d0,  .50d0, .29d0,  1.5d0, !decidforest
     &     34.d0, 0.6d0, 0.4d0, 58.0d0, 57.0d0, 0.2d0,0.093,2, 
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 23.39d0, -0.0540d0), !late successional
     &     pftype(1,.true., 2,-153.d0,  .50d0, .25d0,  0.9d0, !evergrneedle
     &     9.5d0, 1.2d0, 4.d0,42.0d0, 80.0d0,0.25d0,0.184d0,1, !SLA for Pinus sylvestris, Pensa and Sellin (2002).SLA-CLM 10. lrage WAS 5.0!!!! Pinus sylvestris 2-4 years!!!-NYK
!     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0), !early succ
     &     0.0450d0, 1.683d0, 0.1617d0, 2.1536d0, 22.79d0, -0.0445d0), !late succ
     &     pftype(1,.true., 1, -153.d0,  .60d0, .29d0,  1.1d0, !troprainf
     &     9.9d0, 0.5d0, 1.8d0, 41.0d0, 40.0d0, 0.2d0,0.075d0,1,
     &     0.0347d0, 1.560d0, 0.0816d0, 2.306d0, 34.62d0, -0.0232d0),
     &     pftype(2,.false., 1,-153.d0,  .45d0, .27d0,  1.3d0, !crops
     &     15.d0, 0.6d0, 1.1d0, UNDEF, 52.5d0, 0.16d0,0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0)
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  0.76d0)&
!     &     0.0d0,0.0d0,0.0d0,0.0d0)
     &     /)

      !**NOTE:  Above, sstar and swilt are guesses for all except grassland
      ! and savanna.


      !***************************************************
      !*         END - GISS VEGETATION TYPES             *
      !***************************************************
      end module ent_pfts
