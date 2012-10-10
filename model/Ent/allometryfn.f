      module allometryfn
!@sum Physics-level functions to calculate intrinsic plant and canopy properties.  
!@sum NO knowledge of entcells, patches, or 
!@sum Needed to allow initialization of vegetation structure to use 
!@sum physics-level functions without causing circular dependencies.
!@sum Particularly functions for allometry, heat capacity.
!@+auth NYK  July 2012

!@sum 7/2/2012 NYK - Only moved the affected functions.  
!@sum     Eventually may want to move all allometry functions from phenology.f

      use ent_const
      use ent_pfts

      implicit none

      public do_geo
      public GRASS, HERB, SHRUB, TREE, BARE
      public acr, bcr, bR, DBHBAmax_cm
      public GISS_shc, wooddensity_gcm3, init_rootdistr
      public dbh2Cdead, dbh2Cfol,dbh2height, height2dbh,height2Cfol
      public dDBHdCdead, Cdead2dbh
      public Cfol_fn, Csw_fn
      public Crown_rad_allom
      public update_plant_cpools, init_Clab, nplant

      logical, save :: do_geo = .false.!Can change on initialization

      !*q: ratio of root to leaf biomass (unitless) (value from ED)
      real*8, parameter :: q=1.0d0 
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282 !NOTE: This value corrects an error in the coefficient in Moorcroft et al. (2001) Appendix D, which had the value too small by a factor of 100, at 0.00128.
!      real*8, parameter :: iqsw=1000.0d0/3900.0d0/2.0d0

      !*hw_fract: ratio of above ground stem to total stem (stem plus structural roots) (value from ED)
      real*8, parameter :: hw_fract = 0.70d0 
      real*8, parameter :: B2C = 0.5d0  !Inverse, to avoid having to divide.

      !* Entpar:  NEED TO BE ADDED TO ent_pfts.f.  DEFINE HERE TEMPORARILY.
      !wdens_g_cm3 is not calculated from wooddensity_gcm3 but is from data.
      !Moved arrays to ent_pfts_ENT.f and ent_pfts_ENT_FLUXNET.f

      contains
!*************************************************************************

      subroutine init_do_geo_flag
      !Flag to do calculations for geographic version instead of Matthews.
      !NYK - created this for modules that call GISS_shc and that need to
      !     calculate mean annual LAI, which is currently really clunky.

      do_geo = .true.
      end subroutine init_do_geo_flag


      integer function gform(PFT)
      integer :: PFT

      if ((PFT.ge.1).and.(PFT.le.8)) then
         gform = TREE
      elseif ((PFT.eq.9).or.(PFT.eq.10)) then
         gform = SHRUB
      elseif ((PFT.ge.11).and.(PFT.le.14)) then
         gform = GRASS
      elseif ((PFT.eq.15)) then
         gform = HERB
      elseif (PFT.eq.16) then
         gform = TREE
      else
         gform = BARE
      endif
      end function gform

!*************************************************************************

      real*8 function GISS_shc(meanLAI) Result(shc)
!@sum Returns GISS GCM specific heat capacity for cohort.
!     meanLAI = (sum over pfts) 
!           {.5d0*(alamax(anum) + alamin(anum)) * vfraction }/sum(vfraction)
!     I.e. GISS GCM computes shc_entcell = shc(mean entcell LAI*vfaction)
!     instead of shc_entcell = mean(shc(patch LAI)*vfraction)
      use ent_const
      real*8,intent(in) :: meanlai    

      !Seems like this ought to use actual LAI, too? - NK
      !-----Local----

      shc = (.010d0+.002d0*meanLAI+.001d0*meanLAI**2)*shw*rhow

      end function GISS_shc

!*************************************************************************
      subroutine init_rootdistr(fracroot, pft)
      !use ent_GISSveg, only : GISS_calc_fracroot
      real*8 :: fracroot(N_DEPTH)
      integer :: pft

      !GISS TEMPORARY - prescribed roots should be passed in
!      call GISS_calc_fracroot(fracroot, pft)

      ! at least set it to zero, since it is called in zero_cohort...
      fracroot(:) = 0.d0

      !* Prognostic roots calculated here *!.

      end subroutine init_rootdistr

!*********************************************************************

      real*8 function wooddensity_gcm3(pft) Result(wooddens)
      !This function might be replaced by array wdens_g_cm3.
      use ent_pfts, only : pfpar
      integer,intent(in) :: pft
      !* Wood density (g cm-3). Moorcroft et al. (2001).

      wooddens = max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0))

      end function wooddensity_gcm3

!*************************************************************************
      real*8 function sla(pft,llspan)
      integer, intent(in) :: pft
      real*8, intent(in) :: llspan
      
      if (pfpar(pft)%phenotype .eq. EVERGREEN .and. 
     &    pfpar(pft)%leaftype .eq. BROADLEAF .and.
     &    llspan .gt. 0.d0) then 
         sla = 10.0**(1.6923-0.3305*log10(llspan))
      else 
         sla = pfpar(pft)%sla
      endif
      end function sla
!*************************************************************************

      real*8 function dbh2Cdead(pft,dbh)
      !gC/plant  Cdead as a function of dbh
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dbh2Cdead = B2C * pfpar(pft)%b1Cd * 
     &     dbh**pfpar(pft)%b2Cd * 1000.0d0 !kgC to gC

      end function  dbh2Cdead

!*************************************************************************
      real*8 function Cdead2dbh(pft,Cdead)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

!      Cdead2dbh = (Cdead*B2C/1000.0d0/pfpar(pft)%b1Cd)
!     &            **(1.0d0/pfpar(pft)%b2Cd)

      Cdead2dbh = (Cdead/(B2C*pfpar(pft)%b1Cd*1000.d0))
     &     **(1.0d0/pfpar(pft)%b2Cd)

      end function Cdead2dbh

!*************************************************************************
      real*8 function dDBHdCdead(pft,Cdead)
      !cm/kg  Derivative w.r.t. Cdead of Inverse of dbh2Cdead
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

!      dDBHdCdead=(C2B/1000.0d0/pfpar(pft)%b1Cd)**(1.0d0/pfpar(pft)%b2Cd)
!     &          *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
!     &          /pfpar(pft)%b2Cd

      dDBHdCdead=(1/(B2C*1000.0d0*pfpar(pft)%b1Cd))**
     &     (1.0d0/pfpar(pft)%b2Cd)
     &     *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
     &     /pfpar(pft)%b2Cd

      end function dDBHdCdead

!*************************************************************************
      real*8 function dbh2Cfol(pft,dbh)
      !gC/plant  Maximum Cfol as a function of dbh
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8 :: maxdbh

      maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+a0h(pft))-a0h(pft))  
     &     /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht

c$$$      if (.not.pfpar(pft)%woody) then !herbaceous
c$$$         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      else !woody
c$$$         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      end if

      dbh2Cfol=B2C *pfpar(pft)%b1Cf *1000.0d0 !kgC to gC
     &        * min(dbh, maxdbh)**pfpar(pft)%b2Cf

      end function dbh2Cfol

!*************************************************************************

      real*8 function iqsw_fn(pft)
      !*iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m) (value from ED)
      !3900.0: leaf area per sapwood area (m2/m2) 
      !1000.0: sapwood density (kg/m3)
      !2.0:  biomass per carbon (kg/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282 !NOTE: This value corrects an error in the coefficient in Moorcroft et al. (2001) Appendix D, which had the value too small by a factor of 100, at 0.00128.
      !use ent_prescr_veg, only : wooddensity_gcm3
      integer, intent(in) :: pft
      !--- Local ---
      real*8, parameter :: iqsw0=1000.0d0/3900.0d0/2.0d0
      real*8 :: rho_wood        !kgC/m3
      real*8 :: Csw

      !Yeonjoo's
      !iqsw = iqsw0*wooddensity_gcm3(pft)

      rho_wood = wdens_g_cm3(pft) * 1000. !dry mass, convert gB/cm3 to kgB/m3
      iqsw_fn = B2C*rho_wood/3900. !kgCwood/m3wood * m2sw/m2fol = kgCsw/m2fol/msw

      end function iqsw_fn
!*************************************************************************
      real*8 function Csw_fn(pft, DBH_cm, h_m)
      !C_sapwood (gC/plant)
      !Csw = iqsw (kgCsw/m2fol/msw) * SLA (m2fol/kgCfol) * h (msw) * Cfol (kgC)
      !    = kgCsw or 1000.*gCsw
      !0.5 converts density from dry biomass to C, both rho_wood and Cfol
      !From data of Cleary et al (2008) for sagebrush, it appears that
      ! the piperatio 3900 = h*factor.  The factor varies with h
      ! to keep h*factor constant.  So I will look at cutting h out of eqn.-nk
      integer :: pft
      real*8 :: DBH_cm, h_m
      !-------
      real*8 :: rho_wood  !kgC/m3
      real*8 :: Csw
      real*8 :: iqsw

      iqsw = iqsw_fn(pft)
      !rho_wood = wdens_g_cm3(pft) * 1000. !dry mass, convert gB/cm3 to kgB/m3
      !iqsw = B2C*rho_wood/3900. !kgCwood/m3wood * m2sw/m2fol = kgCsw/m2fol/msw
	
      if (pfpar(pft)%woody) then
!         if (DBH_cm.lt.10.) then

!         else
            Csw = iqsw * pfpar(pft)%sla * h_m
     &           * B2C *pfpar(pft)%b1Cf*(DBH_cm**pfpar(pft)%b2Cf)
     &        *1.d3 !kgC to gC
!         endif
      else 
         Csw = (0.0)
      endif

      Csw_fn = Csw
      end function Csw_fn


!*************************************************************************
      real*8 function Cfol_fn(pft,DBH_cm, h_m) 
      !Maximum C_foliage (gC/plant) for woody or herbaceous plants.
      integer :: pft
      real*8 :: DBH_cm, h_m
      !-------
      integer :: vform
      real*8 :: rho_wood
      real*8  :: Cfol !gC/plant

      vform = gform(pft)

      if (pfpar(pft)%woody) then
         rho_wood = wdens_g_cm3(pft) * 1000. !convert from g/cm3 to kg/m3
         !Csw = Csw_fn(pft,DBH_cm,h_m)
         if ((h_m.gt.0.0).and.(DBH_cm.gt.0.0)) then
            !* Based on foliage dependence on sapwood.- NYK
            !Cfol = (Csw*3900./(0.5*rho_wood*sla_m2_kgC(pft)*h_m))
            ! Yeonjoo's
            Cfol = dbh2Cfol(pft,DBH_cm)
         else
            Cfol = 0.0
         endif
      elseif ((vform.eq.HERB).or.(vform.eq.GRASS)) then
         if (h_m.gt.0.0) then
            !This departs from Yeonjoo's height2Cfol, instead uses
            !Nancy's formulation.
            !Cfol = 0.5*(ag(pft)*h_m**bg(pft))
            Cfol = B2C*(pfpar(pft)%b1Cf*h_m**pfpar(pft)%b2Cf)*1000.d0
         else
            Cfol = 0
         endif
      else
         Cfol = (0.0)
      endif
      if (isNaN(Cfol)) then
         write(*,*) "Cfol=NaN",pft,vform, rho_wood,DBH_cm,h_m
     &        ,pfpar(pft)%sla,pfpar(pft)%woody
     &        ,pfpar(pft)%b1Cf,pfpar(pft)%b2Cf!,ag(pft),bg(pft)
         !STOP
      endif
      Cfol_fn = Cfol
      end function Cfol_fn

!*************************************************************************
      real*8 function dbh2height(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      !----
      real*8 :: hcrit !max height

      hcrit=a0h(pft)+pfpar(pft)%b1Ht-EPS
      dbh2height = min(a0h(pft) + pfpar(pft)%b1Ht * 
     &             (1.0-exp(pfpar(pft)%b2Ht*dbh)), hcrit)

      end function dbh2height
!*************************************************************************

      real*8 function height2dbh(pft,h)
      !NK - fixed to avoid negative heights for h=0.
      !   !!!a0h>0 is problematic, since get non-zero height for dbh=0.
      !   !!!Need to fix allometry for pfts with a0h>0!!!
      integer,intent(in) :: pft
      real*8, intent(in) :: h
      real*8 :: hcrit, hin

      hcrit=a0h(pft)+pfpar(pft)%b1Ht-EPS
      hin=min(h,hcrit)

      height2dbh = max(0.d0, 
     &     log(1.d0-(hin-a0h(pft))/pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht)
      

      end function height2dbh
!*************************************************************************

      real*8 function height2Cfol(pft,height)
      !Yeonjoo's.  No documentation or units.
      !?? Is this only for annual herbaceous recruitment?
      integer,intent(in) :: pft
      real*8, intent(in) :: height
      real*8,parameter :: h1Cf = 1.66d0 
      real*8,parameter :: h2Cf = 1.50d0 
      real*8,parameter :: nplant = 2500.d0 
      height2Cfol= B2C *h1Cf 
     &        *((height*100.0d0)**h2Cf)/nplant

      end function height2Cfol
!*************************************************************************

      real*8 function Crown_rad_allom(pft, h_m)
      !Return crown radius (m) calculate from allometry.
      !+auth NYK
      integer :: pft
      real*8 :: h_m
      Crown_rad_allom = acr(pft)*h_m

      end function Crown_rad_allom

!*************************************************************************
      real*8 function Crown_rad_max_from_density(ndens) 
      !Return crown horizontal radius (m) based on maximum packing:
      !n = pi/sqrt(12)/(2*(a+da))^2  #Lagrange hexagonal maximum packing density
      ! where n = density, a=crown radius, da = distance between crowns
      !a = 0.5*sqrt(pi/sqrt(12)/n) - da
      real*8 :: ndens !Plant density in #/m^2

      if (ndens.gt.0.d0) then
         Crown_rad_max_from_density = 0.5d0*sqrt(pi/sqrt(12.d0)/ndens)
      else
         Crown_rad_max_from_density = 0.d0
      endif

      end function Crown_rad_max_from_density

!*************************************************************************
      real*8 function crown_radius_horiz(pft,dbh,popdensity)
      !* Return horizontal crown radius (m)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !cm
      real*8, intent(in) :: popdensity !#/m^2

      crown_radius_horiz = min(Crown_rad_max_from_density(popdensity)
     &     ,crown_radius_pft(pft,dbh))

      end function crown_radius_horiz
!*************************************************************************
      real*8 function crown_radius_pft(pft,dbh) Result(cradm)
      !* From Harvard Forest late successional hardward allometry.
      !* with mean conifer dbh_max limit.
      !* Coefficient 0.107 for late-succ hw is approx. mean for all types.
      use ent_pfts, only : is_conifer,is_hw
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      !------
      !dbh_max parameter values from Harvard Forest allometry
      real*8, parameter :: dbh_max_conifer = 42.d0 !cm.  Mean of for pine and late successional conifer.  
      real*8, parameter :: dbh_max_hw = 150.d0 !cm.  
      real*8 :: dbh_max

      if (is_conifer(pft)) then
         dbh_max = dbh_max_conifer
      elseif (is_hw(pft)) then
         dbh_max = dbh_max_hw
      else
         dbh_max = dbh
      endif

      cradm = .107d0 * min(dbh, dbh_max) 

      end function crown_radius_pft
!*************************************************************************
!      real*8 function crown_radius_closed(popdensity) Result(cradm)
!      !* Return plant crown radius (m).
!      !* Assumes closed canopy packing given popdensity, square regular distrib.
!      real*8, intent(in) :: popdensity !#/m^2
!      
!      cradm = 0.5d0*sqrt(1/popdensity)
!      end function crown_radius_closed
!*************************************************************************

      real*8 function crown_radius_vert(h,crx)
      real*8 :: h, crx !Tree height, crown horizontal radius
      !crown_radius_vert = min(0.45*h,crx*2.7d0)
      crown_radius_vert = max(0.45*h,crx)  !##
      end function crown_radius_vert

!*************************************************************************

      subroutine update_plant_cpools(pft,lai,h,dbh,popdens,cpool)
      ! Does NOT update LABILE.
      integer,intent(in) :: pft !plant functional type
      !real*8, intent(in) :: laimax !max lai should not be needed here
      real*8, intent(in) ::lai  !lai
      real*8, intent(in) :: h,dbh,popdens !h(m), dbh(cm),popd(#/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      !----Local-----
      real*8 :: qsw
      !real*8 :: Cfolmax

      !* Initialize
      cpool(:) = 0.d0 

      qsw = pfpar(pft)%sla*iqsw_fn(pft)
      if (.not.pfpar(pft)%woody) qsw = 0.0d0

      if (popdens.eq.0.d0) then
         cpool(FOL) = 0.d0
      else
         cpool(FOL) = lai/pfpar(pft)%sla/popdens *1d3 !Bl
      endif
      !Cfolmax = laimax/pfpar(pft)%sla/popdens *1d3 !Bl - bad old way to calc
      cpool(FR) = q * cpool(FOL)   !Br
      !cpool(SW) = h * qsw* Cfolmax !bad old way to calculate
      cpool(SW) = Csw_fn(pft,dbh,h)
      if (pfpar(pft)%woody) then !Woody
        cpool(HW) = dbh2Cdead(pft,dbh) * hw_fract
        cpool(CR) = cpool(HW) * (1-hw_fract)/hw_fract !=dbh2Cdead*(1-hw_fract)
      else
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif

      end subroutine update_plant_cpools

!**************************************************************************
      subroutine init_Clab(pft,dbh,h,Clabile)
!@sum init_Clab - renamed from prescr_init_Clab (only depends on allometry)
!@sum - Initializes labile carbon pool
!@sum Deciduous woody: Clab = 4 x max Cfol of plant.
!@sum Evergreen woody: Clab = 0.5 x max Cfol of plant.
!@sum 4x requirement is from Bill Parton (personal communication).
!@sum 7/1/2012 - Replaced with different subroutine - NK
      use ent_pfts, only : COVEROFFSET, pfpar !, alamax, alamin
      implicit none
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh !dbh (cm) !0 for herbs
      real*8, intent(in) :: h !height (m)
      !real*8, intent(in) :: n   !Density (#/m^2)
      !real*8, intent(in) :: laimax
      real*8, intent(inout) :: Clabile !g-C/pool/plant
      !----- Local ------
      real*8 :: Cfolmax

      !Cfolmax = laimax/n/pfpar(pft)%sla !kgC/plant - old way based on laimax
      Cfolmax = Cfol_fn(pft,dbh,h) !g-C/individ -- CHECK Cfol_fn!!!

      if (pfpar(pft)%phenotype.eq.EVERGREEN) then
        !Enough to grow peak foliage and fine roots, and some reproduction.
        Clabile = Cfolmax*2.d0
      else
        Clabile = Cfolmax*4.d0
      endif

      end subroutine init_Clab
!*************************************************************************

      real*8 function nplant(pft,dbh,h,lai)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: h
      real*8, intent(in) :: lai

      if (.not.pfpar(pft)%woody) then !grasses/crops/non-woody
         nplant = lai/pfpar(pft)%sla/(height2Cfol(pft,h)/1000.0d0) 
      else
         nplant = lai/pfpar(pft)%sla/(dbh2Cfol(pft,dbh)/1000.0d0)
      end if

      end function nplant


!*************************************************************************

      end module allometryfn

