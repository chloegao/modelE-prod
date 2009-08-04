      module ent_prescr_veg

      use ent_const
      use ent_pfts
      !use GCM_module, only:  GCMi, GCMj !Fix to names from GCM

      implicit none
      private
      save

      public 
     &     prescr_veg_albedo,prescr_calc_rootprof,
     &     prescr_calcconst, prescr_calc_lai
     &     ,alamax !For temporary phenology
      public prescr_calc_shc,prescr_plant_cpools, prescr_init_Clab
      public prescr_get_hdata,prescr_get_initnm,prescr_get_rootprof,
     &     prescr_get_woodydiameter,prescr_get_pop,prescr_get_crownrad
     &     ,prescr_get_soilcolor,ED_woodydiameter,popdensity
#ifdef ENT_STANDALONE_DIAG
      public print_ent_pfts, alamin
#endif 
!*********************************************************************
!--- sand tundr  grass  shrub  trees  decid evrgr  rainf crops bdirt algae  c4grass
      real*8, parameter :: alamax(N_COVERTYPES) =
      !* Matthews LAI *!
!     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0
!     &     ,0.d0, 0.d0, 2.d0 /)
#ifndef FLUXNETINIT
      !* Revised Matthews LAI *!
     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,8.0d0,7.0d0,3.0d0
     &     ,0.d0, 0.d0, 2.d0 /)
#else
      !* FLUXNET Sites *!
     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,6.0d0,7.0d0,4.5d0
     &     ,0.d0, 0.d0, 2.d0 /)
#endif

      real*8, parameter :: alamin(N_COVERTYPES) =
      !* Matthews LAI *!
!     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0
!     &     ,0.d0, 0.d0, 1.d0 /)
#ifndef FLUXNETINIT
      !* Revised Matthews LAI *!
     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 6.0d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)
#else
      !* FLUXNET site *!
     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 4.5d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)
#endif

      integer, parameter :: laday(N_COVERTYPES) =
     $     (/ 0, 196,  196,  196,  196,  196,  196,  196,  196
     &     ,0, 0, 196 /)

      real*8,parameter :: EDPERY=365. !GISS CONST.f

      contains

!***************************************************************************

      subroutine prescr_calcconst() 
      !*    SUBROUTINE TO CALCULATE CONSTANT ARRAYS                     
      use ent_const
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
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22d0 !lignin:nitrogen scalar        
        solubfract(n) = 0.85d0 - (0.018d0 * lnscl)
        structurallignin(n) = (pfpar(n)%lignin * 0.65d0 * 2.22d0) 
     &                      / (1.d0 - solubfract(n))
        lignineffect(n) = exp(-3.0 * structuralLignin(n)) !original CASA
      end do

      end subroutine prescr_calcconst

!**************************************************************************
      
      real*8 function prescr_calc_lai(pnum,jday,hemi ) RESULT(lai)
!@sum Returns GISS GCM leaf area index for given vegetation type, julian day
!@+   and hemisphere
      use ent_const
      !real*8, intent(out) :: lai !@var lai leaf area index - returned
      integer, intent(in) :: pnum !@var pnum cover type
      integer, intent(in) :: jday !@var jday julian day
      integer, intent(in) :: hemi !@var hemi =1 in N. hemisphere, =-1 S.hemi
      !-----Local variables------
      real*8 dphi

      dphi = 0.d0
      if ( hemi < 0 ) dphi = 2d0*pi*.5d0

      !* Return lai *!
      lai =  .5d0 * (alamax(pnum) + alamin(pnum))
     $     + .5d0 * (alamax(pnum) - alamin(pnum))
     $     * cos( 2d0*pi*(laday(pnum)-jday)/dble(EDPERY) + dphi )

      end function prescr_calc_lai


!**************************************************************************
      
      real*8 function prescr_calc_shc(vdata) RESULT(shc)
!@sum Returns GISS GCM specific heat capacity for canopy
      use ent_const
      real*8, intent(in) :: vdata(:) !@var pnum cover type
      !-----Local variables------
      real*8 lai, fsum
      integer pnum

      lai = 0.d0
      fsum = 0.d0

      do pnum=COVEROFFSET+1,COVEROFFSET+N_PFT
        lai = lai + .5d0 * (alamax(pnum) + alamin(pnum))*vdata(pnum)
        fsum = fsum + vdata(pnum)
      enddo

      if ( fsum > EPS ) lai = lai/fsum

      shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow

      end function prescr_calc_shc

!**************************************************************************

      real*8 function prescr_calc_shoot(pnum,hdata,dbhdata)
     &     Result(Bshoot)
!@sum Returns GISS GCM veg shoot kg-C per plant for given vegetation type.
!@+   From Moorcroft, et al. (2001), who takes allometry data from
!@+   Saldarriaga et al. (1998).
      integer,intent(in) :: pnum !@var pnum vegetation type
      real*8,intent(in) :: hdata(N_COVERTYPES), dbhdata(N_COVERTYPES)
      !-----Local-------
      real*8 :: wooddens
      integer :: n !covertypes index

      n = pnum + COVEROFFSET
      wooddens = wooddensity_gcm3(pnum)
      Bshoot = 0.069 * (hdata(n))**0.572
     &     * (dbhdata(n))**1.94 * (wooddens**0.931)

      end function prescr_calc_shoot

!**************************************************************************

      subroutine prescr_veg_albedo(hemi, pft, jday, albedo)
!@sum returns albedo for vegetation of type pft 
!@+   as it is computed in GISS modelE
      integer, intent(in) :: hemi !@hemi hemisphere (-1 south, +1 north)
      integer, intent(in) :: pft !@var pftlike iv, plant functional type
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
      real*8, parameter :: ALBVND(NV,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.067,.089,.089,.078,.100,.067,.061,.089,.000,.200,.089,
     2 .500,.062,.100,.100,.073,.055,.067,.061,.100,.000,.200,.100,
     3 .500,.085,.091,.139,.085,.058,.083,.061,.091,.000,.200,.091,
     4 .500,.080,.090,.111,.064,.055,.061,.061,.090,.000,.200,.090,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.080,.107,.107,.093,.120,.080,.073,.107,.000,.200,.107,
     2 .500,.082,.140,.120,.096,.083,.080,.073,.140,.000,.200,.140,
     3 .500,.119,.145,.167,.119,.115,.100,.073,.145,.000,.200,.145,
     4 .500,.102,.126,.132,.081,.087,.073,.073,.126,.000,.200,.126,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.100,.133,.133,.116,.150,.100,.091,.133,.000,.200,.133,
     2 .500,.103,.175,.150,.120,.109,.100,.091,.175,.000,.200,.175,
     3 .500,.148,.182,.208,.148,.144,.125,.091,.182,.000,.200,.182,
     4 .500,.127,.157,.166,.102,.109,.091,.091,.157,.000,.200,.157,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.020,.027,.027,.023,.030,.020,.018,.027,.000,.200,.027,
     2 .500,.021,.035,.030,.024,.022,.020,.018,.035,.000,.200,.035,
     3 .500,.030,.036,.042,.030,.029,.025,.018,.036,.000,.200,.036,
     4 .500,.026,.032,.033,.020,.022,.018,.018,.032,.000,.200,.032
     *     /),(/NV,4,6/) )
C
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
        albedo(l)=wt1*ALBVND(pft,kh1,l)+wt2*ALBVND(pft,kh2,l)
      enddo

      end subroutine prescr_veg_albedo

!**************************************************************************

      subroutine prescr_calc_rootprof(rootprof, pnum)
      !Return array rootprof of fractions of roots in soil layer
      !Cohort/patch level.
      real*8 :: rootprof(:)
      integer :: pnum !plant functional type
      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
 !--- tundr  grass  shrub  trees  decid evrgr  rainf crops bdirt algae  c4grass
      real*8, parameter :: aroot(N_COVERTYPES) = 
     $     (/ 0.d0,12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: broot(N_COVERTYPES) = 
     $     (/ 0.d0, 1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
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
        frdn=aroot(pnum)*z**broot(pnum) !cumulative root distrib.
        !frdn=min(frdn,one)
        frdn=min(frdn,1d0)
        if(l.eq.n)frdn=1.
        rootprof(l) = frdn-frup
        frup=frdn
      end do
      !Return rootprof(:)
      end subroutine prescr_calc_rootprof

!**************************************************************************

      subroutine prescr_get_rootprof(rootprofdata)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH) 
      !---Local--------
      integer :: pnum !plant functional type      

      do pnum=1,N_COVERTYPES
        call prescr_calc_rootprof(rootprofdata(pnum,:), pnum)
        !Return array rootprof of fractions of roots in soil layer
        !by vegetation type.
      end do
      end subroutine prescr_get_rootprof

!**************************************************************************

      subroutine prescr_get_hdata(hdata)
      !* Return array parameter of GISS vegetation heights.
      real*8 :: hdata(N_COVERTYPES) 
      !------
      real*8, parameter :: vhght(N_COVERTYPES) =
      !* bsand tundrv  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
      !GISS ORIGINAL
!     $     (/0.d0, 0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0
!     &     ,0.d0, 0.d0, 1.5d0 /)
      !NYK ADJUSTED VALUES
     $     (/0.d0, 0.1d0, 1.5d0,  5d0,  7.1d0,  20d0,  30d0, 25d0,1.75d0
     &     ,0.d0, 0.d0, 1.5d0 /)

      !* Copy prescr code for calculating seasonal canopy height here.
      ! For prescr Model E replication, don't need to fill in an
      ! i,j array of vegetation height, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file entdata
      ! containing vegetation heights.

      !* Return hdata heights for all vegetation types
      hdata = vhght
      end subroutine prescr_get_hdata

!**************************************************************************

      subroutine prescr_get_initnm(nmdata)
!@sum  Mean canopy nitrogen (nmv; g/m2[leaf])
      real*8 :: nmdata(N_COVERTYPES)
      !-------
      real*8, parameter :: nmv(N_COVERTYPES) =
     $     (/0.d0,1.6d0,3.27d0,2.38d0,1.03d0,1.25d0,2.9d0,2.7d0,2.50d0
     &     ,0.d0, 0.d0, 0.82d0 /)

      !* Return intial nm for all vegetation and cover types
      nmdata = nmv
      end subroutine prescr_get_initnm

!*************************************************************************

      subroutine prescr_get_pop(dbhdata,popdata)
      !* Return array of GISS vegetation population density (#/m2)
      !* Derived from Moorcroft, et al. (2001)
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      !---Local-----------
      integer :: n,pft

      popdata(:) = 0.0 !Zero initialize, and zero bare soil.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        popdata(n) = popdensity(pft,dbhdata(n))
      enddo
      end subroutine prescr_get_pop

!*************************************************************************
      real*8 function popdensity(pft,dbh) Result(popdens)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      !---Local-----------
      real*8 :: Blmax, wooddens

      if (.not.pfpar(pft)%woody) then
        popdens = 10.d0       !Grass ##HACK See Stampfli et al 2008 (~25 seedlings/m2 for cover %1-10, but big range)
      else
        wooddens = wooddensity_gcm3(pft)
        Blmax = 0.0419d0 * (dbh**1.56d0) * (wooddens**0.55d0)
        popdens = (alamax(pft+COVEROFFSET)/pfpar(pft)%sla)/Blmax
      endif

      end function popdensity

!*************************************************************************

      subroutine prescr_get_woodydiameter(hdata, wddata)
      !* Return array of woody plant diameters at breast height (dbh, cm)
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: wddata(N_COVERTYPES)
      !----Local---------
      integer :: n,pft

      wddata(:) = 0.0 !Zero initialize.
      do pft = 1,N_PFT
        n = pft + COVEROFFSET
        if (pfpar(pft)%woody) then !Woody
          if (pft.eq.TUNDRA) then
            wddata(n) = ED_woodydiameter(pft,hdata(n)) * 20 !FUDGE UNTIL HAVE MIXED CANOPIES
          else                  !Most trees
            wddata(n) = ED_woodydiameter(pft,hdata(n))
          end if
        endif
      enddo
      end subroutine prescr_get_woodydiameter

!*************************************************************************

      real*8 function ED_woodydiameter(pft,h) Result(dbh)
      !* Return woody plant diameter (m).
      !* From Moorcroft, et al. (2001)
      integer,intent(in) :: pft !plant functional type
      real*8,intent(in) ::  h !height (m)
      !real*8,intent(out) :: dbh !(cm)

      if (pft.eq.SAVANNA) then
        dbh = 30.0 !Estimate from Tonzi Ranch, NYK
      else
        dbh = ((1/2.34)*h)**(1/0.64)
      endif
      
      end function ED_woodydiameter

!*************************************************************************

      subroutine prescr_get_crownrad(popdata,craddata)
      real*8,intent(in) :: popdata(N_COVERTYPES)
      real*8,intent(out) :: craddata(N_COVERTYPES)
      !---Local----
      integer :: n, pft

      craddata(:) = 0.0 !Zero initialize.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        craddata(n) = 0.5*sqrt(1/popdata(n))
      end do

      end subroutine prescr_get_crownrad

!*************************************************************************

      subroutine prescr_plant_cpools(pft, lai, h, dbh, popdens, cpool )
      !* Calculate plant carbon pools for single plant (g-C/plant)
      !* After Moorcroft, et al. (2001). No assignment of LABILE pool here.
      !* Coarse root fraction is estimated from Zerihun (2007) for Pinus radiata
      !*  at h=20 m for evergr, 50% wood is carbon, and their relation:
      !*  CR(kg-C/tree) = 0.5*CR(kg/tree) = 0.5*exp(-4.4835)*(dbh**2.5064).
      !*  The ratio of CR(kg-C/tree)/HW(kg-C/tree) at h=20 m is ~0.153.
      !*  This is about the same as the ratio of their CR biomass/AG biomass.
      integer,intent(in) :: pft !plant functional type
      real*8, intent(in) :: lai,h,dbh,popdens  !lai, h(m), dbh(cm),popd(#/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      !----Local------
      real*8 :: max_cpoolFOL

      ! just in case, set to 0 to avoid possible NaNs
      cpool(:) = 0.d0 !g-C/individual plant

      max_cpoolFOL = alamax(pft+COVEROFFSET)/pfpar(pft)%sla/popdens*1d3
      cpool(FOL) = lai/pfpar(pft)%sla/popdens *1d3!Bl
      cpool(FR) = cpool(FOL)   !Br
      !cpool(LABILE) = 0.d0      !dummy.  For prescribed growth, labile storage is not needed.
!      if (pft.ne.GRASSC3) then  !Woody
      if (pfpar(pft)%woody ) then !Woody
!        cpool(SW) = 0.00128d0 * pfpar(pft)%sla * cpool(FR) * h  !Bsw 0.00128 is an error in ED paper.
        cpool(SW) = 0.128d0 * pfpar(pft)%sla * max_cpoolFOL * h  !Bsw CONSTANT
        cpool(HW) = 0.069d0*(h**0.572d0)*(dbh**1.94d0) * 
     &       (wooddensity_gcm3(pft)**0.931d0) *1d3
 !       cpoolHW  = 0.069d0*(h**0.572d0)*(dbh**1.94d0) * 
!     &       (wooddensity_gcm3(pft)**0.931d0) *1d3
        cpool(CR) = pfpar(pft)%croot_ratio*cpool(HW) !Estimated from Zerihun (2007)
      else
        cpool(SW) = 0.d0
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif

      end subroutine prescr_plant_cpools

!*************************************************************************
      subroutine prescr_init_Clab(pft,n,cpool)
!@sum prescr_init_Clab - Initializes labile carbon pool to 4x mass of 
!@sum (alamax - alamin) for woody and perennial plants.
!@sum 4x requirement is from Bill Parton (personal communication).
!@sum For herbaceous annuals, assume seed provides 0.5 of alamax mass (guess).

      implicit none
      integer, intent(in) :: pft
      real*8, intent(in) :: n !Density (#/m^2)
      real*8, intent(inout) :: cpool(N_BPOOLS) !g-C/pool/plant
      
!      cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/pfpar(pft)%sla/n*1d3 !g-C/individ.

      if (pfpar(pft)%phenotype.ne.ANNUAL) then
        !Enough to grow peak foliage and fine roots.
        cpool(LABILE) = (alamax(pft+COVEROFFSET)-alamin(pft+COVEROFFSET)
     &       )*4.d0/pfpar(pft)%sla/n*1d3 !g-C/individ.
      else
        cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/
     &       pfpar(pft)%sla/n*1d3 !g-C/individ.
      endif

      end subroutine prescr_init_Clab
!*************************************************************************

      real*8 function wooddensity_gcm3(pft) Result(wooddens)
      integer,intent(in) :: pft
      !* Wood density (g cm-3). Moorcroft et al. (2001).

      wooddens = max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0))

      end function wooddensity_gcm3
!*************************************************************************
      subroutine prescr_get_soilcolor(soil_color)
      !* Return arrays of GISS soil color and texture.
      integer, intent(out) :: soil_color(N_COVERTYPES)
      !------

      !* bsand tundr  grass shrub trees  decid evrgr  rainf crops bdirt algae  c4grass
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
     $     (/1, 2, 2,  2, 2, 2, 2, 2
     &     ,2, 2, 2, 2 /)

      soil_color(:) = soil_color_prescribed(:)
      
      end subroutine prescr_get_soilcolor
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

