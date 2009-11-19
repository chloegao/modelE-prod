      module ent_prescr_veg
!KIM - temporary values for ENT pfts
      use ent_const
      use ent_pfts
      !use GCM_module, only:  GCMi, GCMj !Fix to names from GCM

      implicit none
      private
      save

      public 
     &     prescr_veg_albedo,prescr_calc_rootprof,
     &     prescr_calcconst, prescr_calc_lai
     &     ,alamax !For temporary phenology.
      public prescr_calc_shc,prescr_plant_cpools, prescr_init_Clab
      public prescr_get_hdata,prescr_get_initnm,prescr_get_rootprof,
     &     prescr_get_woodydiameter,prescr_get_pop,prescr_get_crownrad
     &     ,prescr_get_soilcolor,ED_woodydiameter,popdensity
#ifdef ENT_STANDALONE_DIAG
      public print_ent_pfts, alamin
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

!KIM - temp. values, NK-updated
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-update, KIM - modified to have a realistic numbers
c$$$      real*8, parameter :: alamax(N_COVERTYPES) =
c$$$     $     (/ 8.0d0, 8.0d0, 10.0d0, 10.0d0, 6.0d0 ,6.0d0, 4.0d0
c$$$     &     ,10.d0, 1.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
c$$$     &     , 4.5d0, 6.0d0, 0.d0, 0.d0/)

      real*8, parameter :: alamax(N_COVERTYPES) =
      !* Matthews LAI *!
!     $     (/ 8.0d0, 8.0d0, 10.0d0, 10.0d0, 6.0d0 ,6.0d0, 4.0d0
!     &     ,10.d0, 1.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
!     &     , 4.5d0, 6.0d0, 0.d0, 0.d0/)
#ifndef FLUXNETINIT
!      !* Revised Matthews LAI *!
     $     (/ 6.0d0, 6.0d0, 8.0d0, 8.0d0, 6.0d0 ,6.0d0, 4.0d0
     &     ,6.d0, 1.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
     &     , 4.5d0, 6.0d0, 0.d0, 0.d0/)
#else
!      !* FLUXNET LAI *!
     $     (/ 6.0d0, 5.0d0, 8.0d0, 8.0d0, 6.0d0 ,6.0d0, 4.0d0
     &     ,6.d0, 1.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
     &     , 4.5d0, 6.0d0, 0.d0, 0.d0/)
#endif

      real*8, parameter :: alamin(N_COVERTYPES) =
      !* Matthews LAI *!
c$$$     $     (/ 6.0d0, 6.0d0, 8.0d0, 8.0d0, 1.0d0, 1.0d0,	1.0d0
c$$$     &     ,8.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.1d0, 1.0d0
c$$$     &     ,1.0d0, 1.0d0, 0.d0, 0.d0 /)
#ifndef FLUXNETINIT
      !* Revised Matthews LAI *!
     $     (/ 5.0d0, 5.0d0, 6.0d0, 6.0d0, 1.0d0, 1.0d0,	1.0d0
     &     ,1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.1d0, 1.0d0
     &     ,1.0d0, 1.0d0, 0.d0, 0.d0 /)
#else
      !* FLUXNET LAI *!
     $     (/ 5.0d0, 4.0d0, 6.0d0, 6.0d0, 1.0d0, 1.0d0,	1.0d0
     &     ,1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.1d0, 1.0d0
     &     ,1.0d0, 1.0d0, 0.d0, 0.d0 /)
#endif

      integer, parameter :: laday(N_COVERTYPES) =
     $     (/ 196, 196, 196, 196, 196, 196, 196
     &     ,196, 196, 196, 196, 196, 196, 196
     &     , 196, 196, 0, 0 /)

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
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22 !lignin:nitrogen scalar        
        solubfract(n) = 0.85 - (0.018 * lnscl)
        structurallignin(n) = (pfpar(n)%lignin * 0.65 * 2.22) 
     &                      / (1. - solubfract(n))
        lignineffect(n) = exp(-3.0 * structuralLignin(n))
      end do
      print *,
     &     'from prescr_calcconst (ent_prescrveg): lignineffect[n_pft]='
     &       ,lignineffect

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

      dphi = 0
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
C
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated: 
! decid_needle=decid fall-winter & needle spring-summer, 
! c3grass_ann =c3grass opposite seasons
!, cropstree = cold_broad
C
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
     1 .061,.061,.067,.067,.100,.100,.078, .100,.067,.089,.089,.089,
     &     .091, .089, .089, .100, .500,.000,
     2 .061,.061,.067,.067,.055,.055,.073, .067,.100,.100,.100,.100,
     &     .089, .100, .100, .055, .500,.000,
     3 .061,.061,.083,.083,.058,.058,.085, .083,.139,.139,.091,.091,
     &     .100, .091, .091, .058, .500,.000,
     4 .061,.061,.061,.061,.055,.055,.064, .061, .111,.111,.090,.090,
     &     .100,.090, .090, .055, .500,.000,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
     1 .183,.183,.200,.200,.300,.300,.233,.300,.200,.267,.267,.267,
     &     .350,.267,.267,.300,.500,.000,
     2 .183,.183,.200,.200,.218,.218,.241,.200,.300,.300,.350,.350,
     &     .364,.350,.350,.218,.500,.000,
     3 .183,.183,.250,.250,.288,.288,.297,.250,.417,.417,.364,.364,
     &     .267,.364,.364,.288,.500,.000,
     4 .183,.183,.183,.183,.218,.218,.204,.183,.333,.333,.315,.315,
     &     .315,.315,.315,.218,.500,.000,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
     1 .183,.183,.200,.200,.300,.300,.233,.300,.200,.267,.267,.267,
     &     .350,.267,.267,.300,.500,.000,
     2 .183,.183,.200,.200,.218,.218,.241,.200,.300,.300,.350,.350,
     &     .364,.350,.350,.218,.500,.000,
     3 .183,.183,.250,.250,.288,.288,.297,.250,.417,.417,.364,.364,
     &     .267,.364,.364,.288,.500,.000,
     4 .183,.183,.183,.183,.218,.218,.204,.218,.333,.333,.315,.315,
     &     .315,.315,.315,.218,.500,.000,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
     1 .073,.073,.080,.080,.120,.120,.093,.120,.080,.107,.107,.107,
     &     .140,.107,.107,.120,.500,.000,
     2 .073,.073,.080,.080,.083,.083,.096,.080,.120,.120,.140,.140,
     &     .145,.140,.140,.083,.500,.000,
     3 .073,.073,.100,.100,.115,.115,.119,.100,.167,.167,.145,.145,
     &     .107,.145,.145,.115,.500,.000,
     4 .073,.073,.073,.073,.087,.087,.081,.087,.132,.132,.126,.126,
     &     .126,.126,.126,.087,.500,.000,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
     1 .091,.091,.100,.100,.150,.150,.116,.150,.100,.133,.133,.133,
     &     .175,.133,.133,.150,.500,.000,
     2 .091,.091,.100,.100,.109,.109,.120,.100,.150,.150,.175,.175,
     &     .182,.175,.175,.109,.500,.000,
     3 .091,.091,.125,.125,.144,.144,.148,.125,.208,.208,.182,.182,
     &     .133,.182,.182,.144,.500,.000,
     4 .091,.091,.091,.091,.109,.109,.102,.109,.166,.166,.157,.157,
     &     .157,.157,.157,.109,.500,.000,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)

     1 .018,.018,.020,.020,.030,.030,.023,.030,.020,.027,.027,.027,
     &     .035,.027,.027,.030,.500,.000,
     2 .018,.018,.020,.020,.022,.022,.024,.020,.030,.030,.035,.035,
     &     .036,.035,.035,.022,.500,.000,
     3 .018,.018,.025,.025,.029,.029,.030,.025,.042,.042,.036,.036,
     &     .027,.036,.036,.029,.500,.000,
     4 .018,.018,.018,.018,.022,.022,.020,.022,.033,.033,.032,.032,
     &     .032,.032,.032,.022,.500,.000
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

!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
      real*8, parameter :: aroot(N_COVERTYPES) = 
     $     (/ 1.1d0, 1.1d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0
     &     ,0.25d0, 0.8d0, 0.8d0, 0.9d0, 0.9d0, 0.9d0, 0.9d0
     &     ,0.9d0, 0.25d0, 0.d0, 0.d0 /)
      real*8, parameter :: broot(N_COVERTYPES) = 
     $     (/ 0.4d0, 0.4d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
     &     ,2.0d0, 0.4d0, 0.4d0, 0.9d0, 0.9d0, 0.9d0, 0.9d0
     &     , 0.9d0, 2.0d0, 0.0d0, 0.0d0 /)

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
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
     $     (/25d0, 25d0, 30d0, 30d0, 7.1d0, 7.1d0, 5d0
     &     , 30d0, 0.1d0, 5d0, 1.5d0, 1.5d0, 1.5d0, 1.5d0
     &     ,1.75d0, 7.1d0, 0.d0, 0.d0 /)

C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4



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
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
     $     (/2.7d0, 2.7d0, 2.9d0, 2.9d0, 1.25d0, 1.25d0, 1.03d0
     &     , 2.9d0, 1.6d0, 2.38d0, 3.27d0, 0.82d0, 3.27d0
     &     , 3.27d0, 2.50d0, 1.25d0, 0.d0, 0.d0 /)

      !* Return intial nm for all vegetation and cover types
      nmdata = nmv
      end subroutine prescr_get_initnm

!*************************************************************************

      subroutine prescr_get_pop(dbhdata,popdata)
      !* Return array of GISS-derived vegetation population density (#/m2)
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
          wddata(n) = ED_woodydiameter(pft,hdata(n))
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

      dbh = ((1/2.34)*h)**(1/0.64)
        
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

      !* Initialize
      cpool(:) = 0.d0
      
      max_cpoolFOL = alamax(pft+COVEROFFSET)/pfpar(pft)%sla/popdens*1d3
      cpool(FOL) = lai/pfpar(pft)%sla/popdens *1d3!Bl
      cpool(FR) = cpool(FOL)   !Br
      !cpool(LABILE) = 0.d0      !dummy.  For prescribed growth, labile storage is not needed.
      if (pfpar(pft)%woody) then !Woody
!        cpool(SW) = 0.128d0 * pfpar(pft)%sla * cpool(FOL) * h  !Bsw
        cpool(SW) = 0.128d0 * pfpar(pft)%sla * max_cpoolFOL * h  !Bsw
        cpool(HW) = 0.069d0*(h**0.572d0)*(dbh**1.94d0) * 
     &       (wooddensity_gcm3(pft)**0.931d0) *1d3
        cpool(CR) =  pfpar(pft)%croot_ratio*cpool(HW) !Estimated from Zerihun (2007)
      else
        cpool(SW) = 0.d0
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif
      !write(97,*) pft, lai, h, dbh, popdens, cpool

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
       if( n > 0.d0 ) then
        cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/
     &       pfpar(pft)%sla/n*1d3 !g-C/individ.
       else
        cpool(LABILE) = 0.d0
       endif
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
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
     $     (/ 2, 2,  2, 2, 2, 2, 2
     &     , 2, 2, 2, 2, 2 ,2, 2, 2, 2, 1, 2 /)

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

