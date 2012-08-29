      module ent_prescr_veg
!@sum ent_prescr_veg - This module should contain routines that are
!@sum primarily specific to the Matthews prescribed vegetation scheme.
      
      use ent_const
!      use ent_pfts
      !use GCM_module, only:  GCMi, GCMj !Fix to names from GCM

      implicit none
      private
      save

      public 
     &     init_params,!prescr_calcconst renamed init_params 
     &     prescr_veg_albedo,prescr_calc_rootprof,
     &     prescr_calc_lai
!      public GISS_shc  !NYK - moved to phenology.f
      public prescr_plant_cpools
      public prescr_calc_hdata
     &     ,prescr_calc_woodydiameter,prescr_get_pop,prescr_get_crownrad
     &     ,prescr_calc_initnm, prescr_calc_rootprof_all
     &     ,prescr_calc_soilcolor
      public ED_woodydiameter,popdensity,Ent_dbh

#ifdef ENT_STANDALONE_DIAG
      public print_ent_pfts
#endif
!*********************************************************************
!* Ent PFTs
!* 1.  evergreen broadleaf early successional (BROADEVERGRTREES1)
!* 2.  evergreen broadleaf late successional  (BROADEVERGRTREES2)
!* 3.  evergreen needleleaf early successional (NEEDLEEVERGRTREES1)
!* 4.  evergreen needleleaf late successional  (NEEDLEEVERGRTREES2)
!* 5.  cold deciduous broadleaf early successional (BROADCOLDDECIDTREES1)
!* 6.  cold deciduous broadleaf late successional  (BROADCOLDDECIDTREES2)
!* 7.  drought deciduous broadleaf	(BROADDRYDECIDTEE)
!* 8.  decidous needleleaf	        (NEEDLEDECIDTREE)
!* 9.  cold adapted shrub               (SHRUBCOLD)
!* 10.  arid adapted shrub              (SHRUBARID)
!* 11.  C3 grass - perennial            (GRASSC3PERENN)
!* 12.  C4 grass - perennial            (GRASSC4PERENN)
!* 13.  C3 grass - annual               (GRASSC3ANN)
!* 14.  arctic C3 grass                 (GRASSC3ARCTIC)
!* 15.  crops - C4 herbaceous           (CROPC4HERB)
!* 16.  crops - broadleaf woody         (CROPTREE)


!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt

!## alamax, alamin,and laday have been moved to ent_pfts_ENT.f and 
!##   ent_pfts_ENT_FLUXNET.f - NK

      real*8,parameter :: EDPERY=365. !GISS CONST.f

      contains

!***************************************************************************

!      subroutine prescr_calcconst() 
      subroutine init_params()
      !*    SUBROUTINE TO CALCULATE CONSTANT ARRAYS                     
      use ent_const
      use ent_pfts
      !--Local------
      integer :: n
      real*8 :: lnscl

      !* annK - Turnover time of litter and soil carbon *!
      do n = 1, N_PFT
        if(pfpar(n)%lrage.gt.0.0d0)then
          annK(n,LEAF)  = (1.0d0/(pfpar(n)%lrage*secpy))
          annK(n,FROOT) = (1.0d0/(pfpar(n)%lrage*secpy))
        else
          annK(n,LEAF)  = 1.0d-40  !CASA originally 1.0e-40
          annK(n,FROOT) = 1.0d-40
        end if

        if(pfpar(n)%woodage.gt.0.0d0)then
          annK(n,WOOD)  = 1.0d0/(pfpar(n)%woodage*secpy)
        else
          annK(n,WOOD)  = 1.0d-40
        end if
       !* iyf: 1/(turnover times) for dead pools.  Want annK in sec-1.
        annK(n,SURFMET)    = 14.8d0  /secpy
        annK(n,SURFMIC)    = 6.0d0   /secpy
        annK(n,SURFSTR)    = 3.9d0   /secpy
        annK(n,SOILMET)    = 18.5d0  /secpy
        annK(n,SOILMIC)    = 7.3d0   /secpy
        annK(n,SOILSTR)    = 4.9d0   /secpy        ! 4.8 in casa v3.0
        annK(n,CWD)        = 0.2424d0/secpy
        annK(n,SLOW)       = 0.2d0   /secpy
        annK(n,PASSIVE)    = 0.1d0 * 0.02d0  /secpy
      enddo

#ifdef DEBUG
      do n = 1,N_PFT
        write(98,*) 'pft',n
        write(98,*) 'annK(n,LEAF)',annK(n,LEAF)
        write(98,*) 'annK(n,FROOT)',annK(n,FROOT)
        write(98,*) 'annK(n,WOOD)' ,annK(n,WOOD) 
        write(98,*) 'annK(n,SURFMET)',annK(n,SURFMET)
        write(98,*) 'annK(n,SURFMIC)',annK(n,SURFMIC)
        write(98,*) 'annK(n,SURFSTR)',annK(n,SURFSTR)
        write(98,*) 'annK(n,SOILMET)',annK(n,SOILMET)
        write(98,*) 'annK(n,SOILMIC)',annK(n,SOILMIC)
        write(98,*) 'annK(n,SOILSTR)',annK(n,SOILSTR)
        write(98,*) 'annK(n,CWD)',annK(n,CWD)    
        write(98,*) 'annK(n,SLOW)',annK(n,SLOW)   
        write(98,*) 'annK(n,PASSIVE)',annK(n,PASSIVE)
      enddo
#endif

      !* solubfrac - Soluble fraction of litter*!
      !structurallignin, lignineffect - frac of structural C from lignin, effect of lignin on decomp -PK 6/29/06 
      do n = 1,N_PFT
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22 !lignin:nitrogen scalar        
        solubfract(n) = 0.85 - (0.018 * lnscl)
        structurallignin(n) = (pfpar(n)%lignin * 0.65 * 2.22) 
     &                      / (1. - solubfract(n))
        lignineffect(n) = exp(-3.0 * structuralLignin(n))
      end do

      end subroutine init_params

!**************************************************************************
      
      real*8 function prescr_calc_lai(pnum,jday,hemi ) RESULT(lai)
!@sum Returns GISS GCM leaf area index for given vegetation type, julian day
!@+   and hemisphere
      use ent_const
      use ent_pfts, only: alamax, alamin, laday
      !real*8, intent(out) :: lai !@var lai leaf area index - returned
      integer, intent(in) :: pnum !@var pnum cover type
      integer, intent(in) :: jday !@var jday julian day
      integer, intent(in) :: hemi !@var hemi =1 in N. hemisphere, =-1 S.hemi
      !-----Local variables------
      real*8 dphi

      dphi = 0
      if ( hemi < 0 ) dphi = 2d0*pi*.5d0

      !* Return lai *!
      lai =  .5d0 * (alamax(pnum) + alamin(pnum))
     $     + .5d0 * (alamax(pnum) - alamin(pnum))
     $     * cos( 2d0*pi*(laday(pnum)-jday)/dble(EDPERY) + dphi )

      end function prescr_calc_lai


!**************************************************************************
      
!     real*8 function GISS_shc(pft) Result(shc)
!!@sum Returns GISS GCM specific heat capacity for cohort.
!      use ent_const
!      use ent_pfts, only: COVEROFFSET, alamax, alamin
!      integer,intent(in) :: pft
!      !Seems like this ought to use actual LAI, too? - NK
!      !-----Local----
!      real*8 :: lai
!      integer :: anum
!
!      anum = pft+COVEROFFSET
!      lai = .5d0*(alamax(anum) + alamin(anum))
!      shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
!
!      end function GISS_shc

!*** MOVED TO phenology.f to keep clean dependencies
!      real*8 function GISS_shc(meanLAI) Result(shc)
!!@sum Returns GISS GCM specific heat capacity for cohort.
!!     meanLAI = (sum over pfts) 
!!           {.5d0*(alamax(anum) + alamin(anum)) * vfraction }/sum(vfraction)
!!     I.e. GISS GCM computes shc_entcell = shc(mean entcell LAI*vfaction)
!!     instead of shc_entcell = mean(shc(patch LAI)*vfraction)
!      use ent_const
!      real*8,intent(in) :: meanlai    
!
!      !Seems like this ought to use actual LAI, too? - NK
!      !-----Local----
!
!      shc = (.010d0+.002d0*meanLAI+.001d0*meanLAI**2)*shw*rhow
!
!      end function GISS_shc

!**************************************************************************
!## This function does not appear to be used - NK
!      real*8 function prescr_calc_shoot(pft,hdata,dbhdata)
!     &     Result(Bshoot)
!!@sum Returns GISS GCM veg shoot kg-C per plant for given vegetation type.
!!@+   From Moorcroft, et al. (2001), who takes allometry data from
!!@+   Saldarriaga et al. (1998).
!     use ent_pfts, only : COVEROFFSET
!      integer,intent(in) :: pft !@var pft vegetation type
!      real*8,intent(in) :: hdata(N_COVERTYPES), dbhdata(N_COVERTYPES)
!      !-----Local-------
!      real*8 :: wooddens
!      integer :: ncov !covertypes index
!
!      ncov = pft + COVEROFFSET
!      wooddens = wooddensity_gcm3(pft)
!      Bshoot = 0.069 * (hdata(ncov))**0.572
!     &     * (dbhdata(ncov))**1.94 * (wooddens**0.931)
!
!      end function prescr_calc_shoot

!**************************************************************************

      subroutine prescr_veg_albedo(hemi, ncov, jday, albedo)
!@sum returns albedo for vegetation of type pft 
!@+   as it is computed in GISS modelE
      !use ent_pfts, only:  COVEROFFSET, albvnd
      use ent_pfts, only:  ALBVND
      integer, intent(in) :: hemi !@hemi hemisphere (-1 south, +1 north)
      !integer, intent(in) :: pft !@var pftlike iv, plant functional type
      integer, intent(in) :: ncov !@var cover type, soil or pft, if pft then ncov=pft+COVEROFFSET
      integer, intent(in) :: jday !@jday julian day
      real*8, intent(out) :: albedo(N_BANDS) !@albedo returned albedo
      !----------Local----------
      integer, parameter :: NV=N_COVERTYPES
      !@var SEASON julian day for start of season (used for veg albedo calc)
C                      1       2       3       4
C                    WINTER  SPRING  SUMMER  AUTUMN
      real*8, parameter, dimension(4)::
     *     SEASON=(/ 15.00,  105.0,  196.0,  288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
!@+   albvnd has been moved to ent_pfts_ENT.f - NK

ccc or pass k-vegetation type, L-band and 1 or 2 for Hemisphere
      integer k,kh1,kh2,l
      real*8 seasn1,seasn2,wt2,wt1
c
c                      define seasonal albedo dependence
c                      ---------------------------------
c
      seasn1=-77.0d0
      do k=1,4
        seasn2=SEASON(k)
        if(jday.le.seasn2) go to 120
        seasn1=seasn2
      end do
      k=1
      seasn2=380.0d0
  120 continue
      wt2=(jday-seasn1)/(seasn2-seasn1)
      wt1=1.d0-wt2
      if ( hemi == -1 ) then    ! southern hemisphere
        kh1=1+mod(k,4)
        kh2=1+mod(k+1,4)
      else                      ! northern hemisphere
        kh1=1+mod(k+2,4)
        kh2=k
      endif

      do l=1,6

!        albedo(l)=wt1*ALBVND(pft+COVEROFFSET,kh1,l)
!     &        +wt2*ALBVND(pft+COVEROFFSET,kh2,l)
        albedo(l)=wt1*ALBVND(ncov,kh1,l)
     &        +wt2*ALBVND(ncov,kh2,l)
      enddo

      end subroutine prescr_veg_albedo

!**************************************************************************

      subroutine prescr_calc_rootprof(rootprof, ncov)
      !Return array rootprof of fractions of roots in soil layer
      !Cohort/patch level.
      use ent_pfts, only: COVEROFFSET, aroot, broot
      real*8 :: rootprof(:)
      integer :: ncov !plant functional type + COVEROFFSET
      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      integer :: n,l
      real*8 :: z, frup,frdn

c**** calculate root fraction afr averaged over vegetation types
      !Initialize zero
      do l=1,N_DEPTH
        rootprof(l) = 0.0
      end do
      do n=1,N_DEPTH
        if (dz_soil(n) <= 0.0) exit !Get last layer w/roots in it.
      end do
      n=n-1
      z=0.
      frup=0.
      do l=1,n
        z=z+dz_soil(l)
        frdn=aroot(ncov)*z**broot(ncov) !cumulative root distrib.
        !frdn=min(frdn,one)
        frdn=min(frdn,1d0)
        if(l.eq.n)frdn=1.
        rootprof(l) = frdn-frup
        frup=frdn
      end do
      !Return rootprof(:)
      end subroutine prescr_calc_rootprof

!**************************************************************************

      subroutine prescr_calc_rootprof_all(rootprofdata)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH) 
      !---Local--------
      integer :: ncov !plant functional type + COVEROFFSET     

      do ncov=1,N_COVERTYPES
        call prescr_calc_rootprof(rootprofdata(ncov,:), ncov)
        !Return array rootprof of fractions of roots in soil layer
        !by vegetation type.
      end do
      end subroutine prescr_calc_rootprof_all

!**************************************************************************

!      subroutine prescr_get_hdata(hdata) !Renamed
      subroutine prescr_calc_hdata(hdata)
      !* Return array parameter of GISS vegetation heights.
       !*## (Can get rid of this subroutine and just assign array vhght)
      use ent_pfts, only : vhght
      real*8 :: hdata(N_COVERTYPES) 
      !------

      !* Copy prescr code for calculating seasonal canopy height here.
      ! For prescr Model E replication, don't need to fill in an
      ! i,j array of vegetation height, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file entdata
      ! containing vegetation heights.

      !* Return hdata heights for all vegetation types
      hdata = vhght
      end subroutine prescr_calc_hdata

!**************************************************************************

      subroutine prescr_calc_initnm(nmdata)
!@sum  Mean canopy nitrogen (nmv; g/m2[leaf])
      !* ##(Can get rid of this subroutine and just assign array vhght)
      use ent_pfts, only : nmv
      real*8 :: nmdata(N_COVERTYPES)
      !-------

      !* Return intial nm for all vegetation and cover types
      nmdata = nmv
      end subroutine prescr_calc_initnm

!*************************************************************************

      subroutine prescr_get_pop(dbhdata,laimaxdata,popdata)
      !* Return array of GISS-derived vegetation population density (#/m2)
      !* Derived from Moorcroft, et al. (2001)
      use ent_pfts, only : COVEROFFSET, alamax
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(in) :: laimaxdata(N_COVERTYPES)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      !---Local-----------
      integer :: ncov,pft

      popdata(:) = 0.0 !Zero initialize, and zero bare soil.
      do pft=1,N_PFT
        ncov = pft + COVEROFFSET
        !popdata(ncov) = popdensity(pft,dbhdata(ncov),alamax(ncov))
        popdata(ncov) = popdensity(pft,dbhdata(ncov),laimaxdata(ncov))
      enddo
      end subroutine prescr_get_pop

!*************************************************************************
      real*8 function popdensity(pft,dbh,LAImax) Result(popdens)
      !* No. per m^2.  From ED
      use ent_pfts, only: pfpar, COVEROFFSET, alamax
      use allometryfn, only : wooddensity_gcm3
       integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: LAImax
      !---Local-----------
      real*8 :: Blmax !Max foliage mass per tree (kg-C/individ)
      real*8 :: wooddens !Wood density (g/cm^3)

      if (.not.pfpar(pft)%woody) then
        popdens = 10.d0       !Grass ##HACK See Stampfli et al 2008 (~25 seedlings/m2 for cover %1-10, but big range)
      else
        wooddens = wooddensity_gcm3(pft)
        Blmax = 0.0419d0 * (dbh**1.56d0) * (wooddens**0.55d0)
        !popdens = (alamax(pft+COVEROFFSET)/pfpar(pft)%sla)/Blmax 
        popdens = (LAImax/pfpar(pft)%sla)/Blmax 
        !leaf area/ground area/(leaf area/kg-C) / (kg-C/individ)
      endif
      end function popdensity

!*************************************************************************

      subroutine prescr_calc_woodydiameter(hdata, wddata)
      !* Return array of woody plant diameters at breast height (dbh, cm)
      use ent_pfts, only : COVEROFFSET
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: wddata(N_COVERTYPES)
      !----Local---------
      integer :: ncov,pft

      wddata(:) = 0.0           !Zero initialize.
      do pft = 1,N_PFT
         ncov = pft + COVEROFFSET
         wddata(ncov) = Ent_dbh(pft,hdata(ncov))
      enddo
      end subroutine prescr_calc_woodydiameter

!*************************************************************************
      real*4 function DBH_from_h (pft,h_m)
      use ent_pfts, only : pfpar
      use allometryfn, only : height2dbh
      !Calculates DBH (cm) from h (m) for woody plants
      !NK -07/02/2012
      !Original relation:  
      !  h_m = a0h + pfpar(pft)%b1Ht*(1-exp(pfpar(pft)%b2Ht*DBH_cm))
      !Inverse relation:   
      !  DBH_cm = log(1-(h_m - a0h)/pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      !This places a cap on DBH if h_m goes above pfpar(pft)%b1Ht, 
      ! which is really wrong, since it is h that is limited while
      ! DBH gets fatter -
      !  but we don't know what the upper limit of h is.
      !  Need to change function unlimited form h = a*dbh^b
      integer :: pft            !plant functional type number
      real*8 :: h_m             !height (m)
      !------

      if (pfpar(pft)%woody) then
         !if ((h_m.gt.0.0).and.(h_m.lt.1.01*a0h(pft))) then
         !   DBH_from_h = 1.01*a0h(pft)
         !elseif (h_m.gt.(a0h(pft)+0.99*pfpar(pft)%b1Ht)) then
         !   write(*,*) "pft, h, pfpar(pft)%b1Ht",pft,h_m,pfpar(pft)%b1Ht
         !   DBH_from_h = log(.01)/pfpar(pft)%b2Ht
         !else
         !   DBH_from_h = log(1-(h_m - a0h(pft))/pfpar(pft)%b1Ht(pft))/
         !     &           pfpar(pft)%b2Ht
         !endif

         DBH_from_h = height2dbh(pft,h_m)
      else
         DBH_from_h = 0.0
      endif

      end function DBH_from_h

!*************************************************************************
      real*8 function Ent_dbh(pft,h) Result(dbh)
      !* Return dbh (cm).  Checks by pft and modifies ED allometry.
      use ent_pfts, only : pfpar, TUNDRA
      integer, intent(in) :: pft
      real*8, intent(in) :: h !(m)

      if (pfpar(pft)%woody) then !Woody
         if (pft.eq.TUNDRA) then !May need separate shrub allometry
            dbh = ED_woodydiameter(pft,h)! * 20.d0  NEW PARAMS OKAY
         else                   !Most trees
            dbh = ED_woodydiameter(pft,h)
            !dbh = height2dbh(pft,h) !phenology.f function
         endif
      else
         dbh = 0.d0 !Not woody.
      endif
      end function Ent_dbh
!*************************************************************************

      real*8 function ED_woodydiameter(pft,h) Result(dbh)
      !* Return woody plant diameter (cm).
      !* From Moorcroft, et al. (2001)
      use allometryfn, only : height2dbh
      integer,intent(in) :: pft !plant functional type
      real*8,intent(in) ::  h !height (m)
      !real*8,intent(out) :: dbh !(cm)
      
      !ED version
      !dbh = ((1/2.34d0)*h)**(1/0.64d0)

      !phenology.f height2dbh version. ###!
      dbh = height2dbh(pft,h)

      end function ED_woodydiameter

!*************************************************************************

      subroutine prescr_get_crownrad(popdata,craddata)
!@sum prescr_get_crownrad - assumes closed-canopy packing of crowns
!@+      in rows and columns (not staggered).
      use ent_pfts, only : COVEROFFSET
      use allometryfn, only : Crown_rad_max_from_density
      real*8,intent(in) :: popdata(N_COVERTYPES)
      real*8,intent(out) :: craddata(N_COVERTYPES)
      !---Local----
      integer :: n, pft

      craddata(:) = 0.0 !Zero initialize.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        !craddata(n) = crown_radius_closed(popdata(n))
        craddata(n) = Crown_rad_max_from_density(popdata(n))
      end do

      end subroutine prescr_get_crownrad

!*************************************************************************

      subroutine prescr_plant_cpools(pft, lai, h, dbh, popdens
     &     , cpool )
      !* Calculate plant carbon pools for single plant (g-C/plant)
      !* After Moorcroft, et al. (2001). No assignment of LABILE pool here.
      !* Coarse root fraction is estimated from Zerihun (2007) for Pinus radiata
      !*  at h=20 m for evergr, 50% wood is carbon, and their relation:
      !*  CR(kg-C/tree) = 0.5*CR(kg/tree) = 0.5*exp(-4.4835)*(dbh**2.5064).
      !*  The ratio of CR(kg-C/tree)/HW(kg-C/tree) at h=20 m is ~0.153.
      !*  This is about the same as the ratio of their CR biomass/AG biomass.
      use ent_pfts, only: COVEROFFSET, pfpar !, alamax
      use allometryfn, only : wooddensity_gcm3, Cfol_fn
      integer,intent(in) :: pft !plant functional type
      real*8, intent(in) :: lai, h,dbh,popdens  !lai, h(m), dbh(cm),popd(/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      !----Local------
      !real*8 :: max_cpoolFOL
      real*8 :: LAmax !m2/plant

      !* Initialize
      cpool(:) = 0.d0
      
      !max_cpoolFOL = alamax(pft+COVEROFFSET)/pfpar(pft)%sla/popdens*1d3
      !max_cpoolFOL = laimax/pfpar(pft)%sla/popdens*1.d3
      !LAmax = laimax/popdens
      LAmax = .001d0*Cfol_fn(pft,dbh,h)*pfpar(pft)%sla
      cpool(FOL) = lai/pfpar(pft)%sla/popdens *1d3 !Bl
      cpool(FR) = cpool(FOL)   !Br
      !cpool(LABILE) = prognostic as diagnostic.
      if (pfpar(pft)%woody) then !Woody
	 !NOTE: Coefficient 0.128 corrects error in Moorcroft et al. (2001)
         !      See detailed notes for iqsw in phenology.f
!        cpool(SW) = 0.00128d0 * pfpar(pft)%sla
!     &      *  max_cpoolFOL * h   !Bsw = 0.00128[kgC/m3]*1d3*(sla*Bfol*1d-3)*h
        cpool(SW) = 0.128d0 *
     &       (LAmax) * h  !Bsw=0.128*LAtot(1d-3*1d3)*h 
        cpool(HW) = 0.069d0*(h**0.572d0)*(dbh**1.94d0) * 
     &       (wooddensity_gcm3(pft)**0.931d0) *1d3
        cpool(CR) =  pfpar(pft)%croot_ratio*cpool(HW) !Estimated from Zerihun (2007)
      else
        cpool(SW) = 0.d0
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif
      !write(997,*) pft, lai, h, dbh, popdens, cpool

      end subroutine prescr_plant_cpools

!*************************************************************************
!      subroutine prescr_init_Clab(pft,n,laimax,Clabile)
!@sum prescr_init_Clab - Initializes labile carbon pool
!@sum Deciduous woody: Clab = 4 x max Cfol of plant.
!@sum Evergreen woody: Clab = 0.5 x max Cfol of plant.
!@sum 4x requirement is from Bill Parton (personal communication).
!@sum 7/1/2012 - Replaced with different subroutine - NK

!! This subroutine has been revised and moved to phenology.f.      
!      subroutine prescr_init_Clab_old(pft,n,cpool)
!!@sum prescr_init_Clab - Initializes labile carbon pool to 4x mass of 
!!@sum (alamax - alamin) for woody and perennial plants. 
!!@sum 4x requirement is from Bill Parton (personal communication).
!!@sum For herbaceous annuals, assume seed provides 0.5 of alamax mass (guess).
!!@sum 7/1/2012 - Replaced with different subroutine - NK
!      use ent_pfts, only : COVEROFFSET, pfpar !, alamax, alamin
!      implicit none
!      integer, intent(in) :: pft
!      real*8, intent(in) :: n !Density (#/m^2)
!      real*8, intent(in) :: laimax, laimin
!      real*8, intent(inout) :: cpool(N_BPOOLS) !g-C/pool/plant
!      
!!      cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/pfpar(pft)%sla/n*1d3 !g-C/individ.
!
!      if (pfpar(pft)%phenotype.ne.ANNUAL) then
!        !Enough to grow peak foliage and fine roots.
!        cpool(LABILE) = (alamax(pft+COVEROFFSET)-alamin(pft+COVEROFFSET)
!     &       )*4.d0/pfpar(pft)%sla/n*1d3 !g-C/individ.
!      else
!       if( n > 0.d0 ) then
!        cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/
!     &       pfpar(pft)%sla/n*1d3 !g-C/individ.
!       else
!        cpool(LABILE) = 0.d0
!       endif
!      endif
!
!      end subroutine prescr_init_Clab_old
!*************************************************************************

cddd      real*8 function wooddensity_gcm3(pft) Result(wooddens)
cddd      use ent_pfts, only : pfpar
cddd      integer,intent(in) :: pft
cddd      !* Wood density (g cm-3). Moorcroft et al. (2001).
cddd
cddd      wooddens = max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0))
cddd
cddd      end function wooddensity_gcm3


      real*8 function wooddensity_gcm3(pft) Result(wooddens)
      !* Returns wood density in total mass per volume (total mass = dry mass = C + N + everything else)
      use ent_pfts, only : pfpar
      integer,intent(in) :: pft
      !* Wood density (g cm-3). Moorcroft et al. (2001).

      if (pfpar(pft)%leaftype.eq.NEEDLELEAF) then
        wooddens = 0.5d0
      else
        wooddens = min(
     &       max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0)), 1.05d0)       
      endif

      end function wooddensity_gcm3

!*************************************************************************

      subroutine prescr_calc_soilcolor(soil_color)
      !* Return arrays of GISS soil color and texture.
      !## Can get rid of this subroutine and replace with array assignment.
      use ent_pfts, only : soil_color_prescribed
      integer, intent(out) :: soil_color(N_COVERTYPES)
      !------

      soil_color(:) = soil_color_prescribed(:)

      end subroutine prescr_calc_soilcolor
!*************************************************************************
#ifdef ENT_STANDALONE_DIAG
      !Assume PS_MODEL=FBB if running Ent_standalone, so can print out FBBpfts.f.
      subroutine print_ent_pfts()
      use ent_const
      use ent_types
      use ent_pfts
      use FarquharBBpspar
      integer pft

      write(*,*) "ent_pfts pfpar:"
      write(*,*) "pst,woody,leaftype,hwilt,sstar,swilt,nf,sla,
     &r,lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype 
     &b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht"
       do pft = 1,N_PFT
         write(*,*) pft, pfpar(pft)
      enddo

      write(*,*) "FarquharBBpspar pftpar: "
      write(*,*) "pst, PARabsorb,Vcmax,m,b,Nleaf"
      do pft = 1,N_PFT
         write(*,*) pft, pftpar(pft)
      enddo
      end subroutine print_ent_pfts
#endif
!*************************************************************************
      end module ent_prescr_veg

