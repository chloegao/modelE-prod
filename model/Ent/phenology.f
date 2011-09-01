      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence
!@auth Y. Kim
#ifdef ENT_STANDALONE_DIAG
#define PHENOLOGY_DIAG
#define DEBUG
#endif

      use ent_types
      use ent_const
      use ent_pfts
 
      implicit none
!      public veg_init
      public clim_stats
      public pheno_update
!      public frost_hardiness  !DEPENDENCY ISSUES WITH BIOPHYSICS
      public veg_update !may change the name into veg_update 
      public litter_cohort, litter_patch   !Now called from veg_update
      public update_plant_cpools

      private pheno_update_coldwoody
      private pheno_update_coldherbaceous
      private pheno_update_drought
      private growth_cpools_active
      private growth_cpools_structural
!      private senesce_cpools
      private recruit_annual
      private photosyn_acclim
      private phenology_diag

      !************************************************************************
      !** GROWTH MODEL CONSTANTS - phenology & carbon allocation 
      !*l_fract: fraction of leaves retained after leaf fall (unitless) (value from ED) 
      real*8, parameter :: l_fract = 0.50d0 
      !*q: ratio of root to leaf biomass (unitless) (value from ED)
      real*8, parameter :: q=1.0d0 
      !*iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m) (value from ED)
      !3900.0: leaf area per sapwood area (m2/m2) 
      !1000.0: sapwood density (kg/m3)
      !2.0:  biomass per carbon (kg/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282
      real*8, parameter :: iqsw=1000.0d0/3900.0d0/2.0d0
      !*hw_fract: ratio of above ground stem to total stem (stem plus structural roots) (value from ED)
      real*8, parameter :: hw_fract = 0.70d0 
      !*C2B: ratio of biomass to carbon (kg-Biomass/kg-Carbon) 
      real*8, parameter :: C2B = 2.0d0 
      !*temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      !*airtemp_par !base temperature to calculate the growing degree days (gdd)
      !*gdd_par1/2/3: paramters to estimate the threshold for gdd
      !*gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: airtemp_par = 5.d0 
      real*8, parameter :: gdd_par1 = -68.d0 
      real*8, parameter :: gdd_par2 = 638.d0
      real*8, parameter :: gdd_par3 = -0.01d0 
      !*gdd_length - tunning parameters (tuned for HF & MMSF)
      real*8, parameter :: gdd_length = 200.d0 
      !*airt_threshold - tunning parameters (tuned for HF & MMSF)
      real*8, parameter :: airt_max_w = 15.d0
      real*8, parameter :: airt_min_w = 5.d0
      !*soilt_threshold - tunning parameters (tunned for Barrow)
      real*8, parameter :: soilt_max = 10.d0
      real*8, parameter :: soilt_min = 0.d0
      !*sgdd_threshold & length - tunning parameters (tunned for Barrow)
      real*8, parameter :: soiltemp_par = 0.d0
      real*8, parameter :: sgdd_threshold = 100.d0
      real*8, parameter :: sgdd_length=50.d0   
      !*ld_threshold (minute): light length constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: ld_threshold = 655.d0
      real*8, parameter :: ld_min =540.d0
      real*8, parameter :: ld_max =550.d0
      !*tsoil_threshold1, tsoil_threshold2 : soil temperature constraint for cold-deciduous woody PFTs (White et al. 1997)
!      real*8, parameter :: tsoil_threshold1 = 11.15d0
!      real*8, parameter :: tsoil_threshold2 = 2.d0
      !*ddfacu: the rate of leaf fall (1/day) (value from IBIS)
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !*betad : water_stress3  - tunning parameters (sstar/swilt determines betad, then check those first and tune these)
      !_w for woody & _h  (tunned for MMSF) - max/min/resistance parameters for woody
      real*8, parameter :: betad_max_w = 0.1d0
      real*8, parameter :: betad_min_w = 0.d0 
      real*8, parameter :: betad_res_w = 0.25d0
      ! _h for herbaceous  (tunned for Vaira/Tonzi) - max/min/resistance prameters for woody
      real*8, parameter :: betad_max_h = 0.9d0 
      real*8, parameter :: betad_min_h = 0.4d0
      real*8, parameter :: betad_res_h = 1.0d0      
      !*light-controll phenology model (Kim et al; originally for ED2)
      !*different from the original implementation, as it cannot be directly implemented due to model difference.
      !*e.g.) PAR instead of Rshort & other differences in parameterization requires the model to be tunned for Ent.
      !*par_turnover_int & par_turnover_slope  - tunning parameters (tunned for TNF; not finalized)
      real*8, parameter :: par_turnover_int = -12.d0 
      real*8, parameter :: par_turnover_slope = 0.18d0  
      !*r_fract: fraction of excess c going to seed reproduction (value from ED)
      real*8, parameter :: r_fract = 0.3d0
      !*c_fract: fraction of excess c going to clonal reproduction - only for herbaceous (value from ED)
      real*8, parameter :: c_fract = 0.7d0
      !*mort_seedling: mortality rate for seedling (value from ED)
      real*8, parameter :: mort_seedling = 0.90d0 

      contains


      !*********************************************************************
      subroutine clim_stats(dtsec, ecp, config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running mean   
      use soilbgc, only : Soillayer_convert_Ent 
      real*8,intent(in) :: dtsec           !dt in seconds
      type(entcelltype) :: ecp      
      type(ent_config) :: config
      logical, intent(in) :: dailyupdate  
      !-----local--------
      type(patch), pointer :: pp
      type(cohort), pointer :: cop 
      !*local variables for entcell-level envrionment variable
      real*8 :: airtemp        !air temperature degC 
      real*8 :: soiltemp       !soil temperature degC
      real*8 :: par            !photosynthetic active radiation (PAR)
      !*local variables for entcell-level variables, 
      !*updated in this subroutine
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature
      real*8 :: par_10d        !10 day running mean of PAR
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature
      !*PAR-limited phenology parameters
      real*8 :: par_crit       !PAR threshold
      logical :: par_limit     !logical whether PAR-limited phenology parameterization is applied or not for certain PFTs
      real*8 :: turnover0      !turnover amplitude, calculated with the phenology parameterization
      real*8 ::  llspan0       !leaf life span, calculated with the phenology parameterization
      !*soil temperature for CASA layers
      real*8 :: Soiltemp2layer(N_CASA_LAYERS)  
     
  
      airtemp = ecp%TairC
      call Soillayer_convert_Ent(ecp%Soiltemp(:), SOILDEPTH_m, 
     &     Soiltemp2layer) 
      soiltemp = Soiltemp2layer(1)

      soiltemp_10d = ecp%soiltemp_10d
      airtemp_10d = ecp%airtemp_10d
      par_10d = ecp%par_10d  
      gdd = ecp%gdd
      ncd = ecp%ncd
      sgdd = ecp%sgdd

      !*10-day running mean of Air Temperature
      airtemp_10d = running_mean(dtsec, 10.d0, airtemp, airtemp_10d)

      !*10-day running mean of Soil Temperature
      soiltemp_10d = running_mean(dtsec, 10.d0, soiltemp, soiltemp_10d)

      !*10-day running mean of PAR
      par = ecp%IPARdif + ecp%IPARdir !total PAR is the sum of diffused and direct PARs
      par_10d = running_mean(dtsec, 10.d0, par, par_10d)

      !*daylength
      if (ecp%CosZen > 0.d0) then
         ecp%daylength(2) = ecp%daylength(2) + dtsec/60.d0
      end if

      !*GDD & NCD - Update Once a day 
      if (dailyupdate) then
         !*Calculate Growing degree days
         if (airtemp_10d .ge. airtemp_par) then
            gdd = gdd + ( airtemp_10d - airtemp_par )
         end if
         !*Calculate Growing degree days for soil temperature 
         if (soiltemp_10d .ge. soiltemp_par) then
            sgdd = sgdd + ( soiltemp_10d -soiltemp_par )
         end if
         !*Number of chilling days 
         !number of days below airtemp_par - chilling requirements
         if (airtemp_10d .lt. airtemp_par) then
            ncd = ncd +  1.d0
         end if
         !*If the season is fall or not (if fall, it's 1; else, it's 0)
         !1) it is to control the phenological status      
         !2) it is determined according to whether the daylength is decreasing (i.e., fall) or not. 
         if (NInt(ecp%daylength(2)) .lt. NInt(ecp%daylength(1)) ) then
            ecp%fall = 1
         else if (NInt(ecp%daylength(2)).gt.NInt(ecp%daylength(1))) then
            ecp%fall = 0
         end if 
       end if

      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 

        cop => pp%tallest
        do while(ASSOCIATED(cop))
          
          !*10-day running mean of stressH2O (betad)
          cop%betad_10d = running_mean(dtsec, 10.d0, 
     &                    cop%stressH2O, cop%betad_10d)
     &                     
          !*Daily carbon balance
          !it is used for the carbon allocation
          cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n*1000.d0

          !*********************************************
          !* evergreen broadleaf - PAR limited - not finalized yet!
          !********************************************* 

          !if it is evergreen & broadleaf, radiation-limited phenology is working    
          par_limit = ((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &                (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
          !par_limit = .false. !temp. suppress
          !raidation-limited phenology model 
          if (par_limit) then
             par_crit = - par_turnover_int/par_turnover_slope 

             !calculate the turnover amplitude 
             !(relative ratio of turnover compared to its intrinsic turnover rate)
             !based on PAR
             turnover0 = min(100.d0, max(0.01d0, 
     &          par_turnover_slope*par_10d + par_turnover_int))

             if (par_10d .lt. par_crit) turnover0 = 0.01d0

             !calculate 10 day running mean of turnover amplitude
             cop%turnover_amp = running_mean(dtsec,10.d0, 
     &                          turnover0, cop%turnover_amp)

             !calculate the leaf life span based on turnover amplitude
             !lrage is in year, llspan is in month, and then 12 is used to convert the units
             llspan0 = pfpar(cop%pft)%lrage*12.d0/cop%turnover_amp 

             !calculate 90 day running mean of llspan
             cop%llspan = running_mean(dtsec, 90.d0,llspan0,cop%llspan)

          else

             cop%turnover_amp = 1.d0 
             cop%llspan = -999.d0

          endif
          

          !**************************************************************
          !* Update photosynthetic acclimation factor for evergreen veg
          !**************************************************************
          if (config%do_frost_hardiness) then 
             if (((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &           (pfpar(cop%pft)%leaftype.eq.NEEDLELEAF)).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDECID).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDROUGHTDECID)) then
                call photosyn_acclim(dtsec,airtemp_10d,cop%Sacclim) 
             else
                cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
             endif
	  else
              cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
          endif

          cop => cop%shorter  
        end do        
        pp => pp%younger 
      end do 

      ecp%soiltemp_10d = soiltemp_10d
      ecp%airtemp_10d = airtemp_10d
      ecp%par_10d = par_10d
      ecp%gdd = gdd
      ecp%ncd = ncd
      ecp%sgdd = sgdd

      end subroutine clim_stats
      !*********************************************************************   
      subroutine pheno_update(pp)
!@sum Update statstics for phneology_update    
!@sum PHYSICAL time step.
      use ent_const

      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
      integer :: pft
      integer :: phenotype
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature    
      real*8::  betad_10d      !10 day running mean of betad (calculated with stressH2O)  
      real*8 :: ld             !day length in minutes
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature 
      real*8 :: airt_adj       !adjustment for air temperature threshold (airt_max & airt_min)
      real*8 :: soilt_adj      !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !*phenofactor : phenological elongation factor, ranging 0 (no leaf) to 1 (full leaf)
      !_c for the cold-deciduous & _d for the drought-deciduous
      real*8 :: phenofactor
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      !*phenostatus :  define phenological status
      !1 - no leaf (phenofactor, equal to 0; after leaf-off in the  fall, until leaf green-up in the next spring)
      !2 - growing leaf (phenofactor, increasing from 0 to 1 in the spring) 
      !3 - leaf in full growth (phenofactor, equal to 1)
      !4 - leaf senescence (phenofactor, decreasing from 1 to 0 in the fall)
      real*8 :: phenostatus
      real*8 :: phenostatus_c
      real*8 :: phenostatus_d
      !*phenogy type 
      logical :: cold_limit    !.true.=cold deciduous
      logical :: drought_limit  !.true.=drought deciduous
      !*GDD treshold (White et al.)      
      real*8 :: gdd_threshold
      !*whther the season is fall or not (determined in clim_stats)
      logical :: fall         !.true. for fall
      !*whether  PFT is wood or not
      logical :: woody
      !*X day to mature: mature=1+X/1000 
      !1.1 means 100 days to mature.
      !it is devised to prevent abrupt grass green-up 
      !in the middle of winter due to a couple of warm days. 
      real*8 :: mature = 1.1d0 
 
      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      if ( pp%cellptr%fall == 1 ) then
        fall = .true.
      else
        fall = .false.
      endif
      sgdd = pp%cellptr%sgdd
      ld = pp%cellptr%daylength(2)

      cop => pp%tallest
      do while(ASSOCIATED(cop))                
         phenofactor_c=cop%phenofactor_c
         phenofactor_d=cop%phenofactor_d
         phenofactor=cop%phenofactor
         phenostatus_c=cop%phenostatus
         phenostatus_d=cop%phenostatus
         betad_10d=cop%betad_10d
         pft=cop%pft
         phenotype=pfpar(pft)%phenotype
         woody = pfpar(pft)%woody

         !***********************************************
         !*Determine whther PFT cold or drought deciduous
         !***********************************************
         if (phenotype .eq. COLDDECID) then 
            cold_limit = .true.
            drought_limit = .false.
         else if (phenotype .eq. DROUGHTDECID) then 
            cold_limit = .true.
            drought_limit = .true.
         else if (phenotype .eq. EVERGREEN) then
            cold_limit = .false.
            drought_limit = .false.
         else !any of cold and drought deciduous
            cold_limit = .true.
            drought_limit = .true.            
         end if

         !*Set the air/temperature adjustment for specific PFTs
         airt_adj=0.d0
         if (pft.eq.DROUGHTDECIDBROAD)airt_adj=10.d0
         soilt_adj=0.d0
         if (pft .eq. GRASSC3ARCTIC)soilt_adj=-5.d0      


         !*******************************************
         !*Update the phenology for Cold-deciduous
         !*******************************************
         if (cold_limit)then

           !*Cold-deciduous Woody
           if (woody) then  
             !*Update phenofactor and phenostatus
             call pheno_update_coldwoody(cop%phenostatus, 
     i            fall, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, 
     o            phenofactor_c, phenostatus_c)

          !*Cold-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_coldherbaceous(cop%phenostatus,
     i            fall, soiltemp_10d, soilt_adj, 
     o            sgdd, 
     o            phenofactor_c, phenostatus_c)
           end if

         else !cold_limit=.false.
           phenofactor_c = 1.d0
           phenostatus_c  = 3.d0
         end if          
                 
         !*******************************************
         !*Update the phenology for Drought-deciduous
         !*******************************************
         if (drought_limit)then

           !*Drought-deciduous Woody
           if (woody) then 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_w, betad_max_w, betad_res_w, 
     o            phenofactor_d, phenostatus_d)

         !*Drought-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_h, betad_max_h, betad_res_h, 
     o            phenofactor_d, phenostatus_d)
           end if   
    
         else !drought_limit=.false.
           phenofactor_d = 1.d0
           phenostatus_d = 3.d0
         end if  

         !*******************************************
         !*Update the phenology according to PFTs 
         !*******************************************
         if (phenotype .eq. COLDDECID) then 
            phenofactor = phenofactor_c
            phenostatus = phenostatus_c
         else if (phenotype .eq. EVERGREEN) then !leaf in full growth
            phenofactor = 1.d0 
            phenostatus = 3.d0
         else if (phenotype .eq. DROUGHTDECIDBROAD .and. 
     &            .not.phenostatus_c.lt.3.d0) then
             phenofactor = phenofactor_c
             phenostatus = phenostatus_c
         else !any of cold and drought deciduous
            phenofactor = phenofactor_c * phenofactor_d   
           if((phenostatus_c.ge.4.d0.and.phenostatus_d.ge.2.d0).or.
     &       (phenostatus_d.ge.4.d0.and.phenostatus_c.ge.2.d0))then
             phenostatus = max(phenostatus_c, phenostatus_d)
           else
             phenostatus = min(phenostatus_c, phenostatus_d)
           end if
          end if
         
         !*increment phenostatus by 0.001 
         !to track how many days after phenostatus has been changed
         if (aint(cop%phenostatus).eq.aint(phenostatus))then
            phenostatus = cop%phenostatus + 1.d0/1000.d0  
         end if
    
#ifdef DEBUG
             write(202,'(3(i5),100(1pe16.8))') pp%cellptr%fall
     &      ,phenofactor
     &      ,phenofactor_c,phenofactor_d
     &      ,phenostatus, cop%phenostatus
     &      ,phenostatus_c, phenostatus_d
     &      ,betad_10d,soiltemp_10d, mature
     &      ,gdd, ncd
#endif
         
         cop%phenofactor_c=phenofactor_c
         cop%phenofactor_d=phenofactor_d
         cop%phenofactor=phenofactor
         cop%phenostatus=phenostatus
   
         cop => cop%shorter 
      
      end do   

      pp%cellptr%gdd = gdd
      pp%cellptr%ncd = ncd
      pp%cellptr%sgdd = sgdd
      

      end subroutine pheno_update
      !*********************************************************************  
      subroutine pheno_update_coldwoody(phenostatus, 
     i            fall, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, 
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous woody PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall        !.true. if the season is fall
      real*8, intent(in) :: airtemp_10d !10 day running mean of air temperature
      real*8, intent(in) :: airt_adj    !adjustment for air temperature threshold (airt_max & airt_min)
      real*8, intent(in) :: ld          !day length in minutes
      !in/output variables
      real*8, intent(inout) :: gdd      !growing degree days, based on air temperature
      real*8, intent(inout) :: ncd      !number of chilling days, based on air temperature  
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)
      !local variables
      real*8 :: gdd_threshold   !GDD treshold (White et al.)      

      !*GDD threshold for leaf green-up 
      gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)  

      !*Leaf-on in the spring, triggered by thermal sum
      !if gdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if ((.not. fall) .and.
     &   (phenostatus.lt.3.d0).and.(gdd.gt.gdd_threshold))then  
         !determine phenofactor by scaling gdd with gdd_threshold and gdd_length
         phenofactor_c = min (1.d0,(gdd-gdd_threshold)/gdd_length)
         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0 !growing leaf
         else 
            phenostatus_c = 3.d0 !leaf in full grwoth
         end if
      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by air temperature
      !if air temperature is falling below its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and. 
     &   (phenostatus.ge.3.d0).and.
     &   (airtemp_10d.lt.airt_max_w+airt_adj)) then
         !determine phenofactor by scaling air temperature
         !with its minimum(min+adj) and maximum(max+adj).
         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      (airtemp_10d-airt_min_w-airt_adj)/
     &      (airt_max_w-airt_min_w)))
         if (phenofactor_c .eq. 0.d0) then
            phenostatus_c = 1.d0   !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0   !leaf senescence
         end if 
      end if
      !*Leaf-off triggered by day-length
      !if day length is falling shorter than its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and.
     &   (phenostatus.ge.3.d0).and.(ld.lt.ld_max)) then
         !dtermine phenofactor by scaling the day length with its min and max.
         phenofactor_c = min(phenofactor_c, max(0.d0,
     &       (ld - ld_min)/(ld_max-ld_min)))
         if (phenofactor_c .eq. 0.0d0) then
            phenostatus_c =1.d0    !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else
            phenostatus_c = 4.d0   !leaf senescence
         end if
      end if    
       
      end subroutine pheno_update_coldwoody
      !*********************************************************************  
      subroutine pheno_update_coldherbaceous(phenostatus, 
     i            fall, soiltemp_10d, soilt_adj,
     o            sgdd,  
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous herbaceous PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall         !.true. if the season is fall
      real*8, intent(in) :: soiltemp_10d !10 day running mean of soil temperature
      real*8, intent(in) :: soilt_adj    !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !in/output variables
      real*8, intent(inout) :: sgdd      !growing degree days, based on soil temperature
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)
        
      !*Leaf-on in the spring, triggered by thermal sum, based on soil temperature
      !if sgdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if ((.not. fall) .and.
     &   (phenostatus.lt.3.d0) .and.
     &   (sgdd.gt.sgdd_threshold)) then
        !determine phenofactor by scaling sgdd with gsdd_threshold and sgdd_length
         phenofactor_c  
     &      = min (1.d0,(sgdd-sgdd_threshold)/sgdd_length)
         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0    !growing leaf
         else 
            phenostatus_c = 3.d0    !leaf in full growth
         end if
      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by soil temperature
      !if soil temperature is falling below its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and.
     &   (phenostatus.ge.3.d0).and.
     &   (soiltemp_10d.lt.soilt_max+soilt_adj)) then
         !determine phenofactor by scaling soil temperature
         !with its minimum(min+adj) and maximum(max+adj).
         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      (soiltemp_10d-soilt_min-soilt_adj)/(soilt_max-soilt_min)))
         if (phenofactor_c .eq. 0.d0) then
            phenostatus_c = 1.d0    !no leaf
            sgdd = 0.d0             !zero-out sgdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0    !leaf senescence
         end if 
      end if
    

      end subroutine pheno_update_coldherbaceous
      !*********************************************************************  
      subroutine pheno_update_drought(phenostatus, 
     i            mature, betad_10d, betad_min, betad_max, betad_res,
     o            phenofactor_d, phenostatus_d)
!@sum Update phenology for drought-deciduous PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      real*8, intent(in) :: mature      !X day to mature: mature=1+X/1000  (refer pheno_update for details) 
      real*8, intent(in) :: betad_10d   !10 day running mean of betad (calculated with stressH2O)  
      real*8, intent(in) :: betad_min   !betad minimum
      real*8, intent(in) :: betad_max   !betad maximum
      real*8, intent(in) :: betad_res   !betad resistance
      !output variables
      real*8, intent(out) :: phenofactor_d !phenological factor for drought deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_d !phenological status for drought deciduous (refer pheno_update for details)

      !*Leaf-on in the spring, triggered by water stress 
      !if betad is larger than its minimum 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X), 
      !and after long enough after the leaf-off), 
      !determine the phenofactor and corresponding phenostatus.
      if ((phenostatus .gt. mature) .and.
     &   (phenostatus.lt.3.d0).and. (betad_10d.gt.betad_min))then
         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = min(1.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)
         if (phenofactor_d .ge. 0.95d0) then
            phenostatus_d = 3.d0    !leaf in full growth
         else
            phenostatus_d = 2.d0    !growing leaf
         end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by water stress
      !if betad is falling smaller than its maximum 
      !in the fall (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      else if ((phenostatus.ge.3.d0).and. (betad_10d.lt.betad_max))then
         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = max(0.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)
         if (phenofactor_d .le. EPS) then
            phenostatus_d = 1.d0    !no leaf
         else
            phenostatus_d = 4.d0    !leaf senescence
         end if 
      end if    
      
      end subroutine pheno_update_drought
      !*********************************************************************   
      subroutine veg_update(pp,config)
!@sum Update the vegetation state and carbon pools:
!@sum DAILY call.
!@sum LAI, senescefrac, DBH, height 
!@sum carbon pools of foliage, sapwood, fineroot, hardwood, coarseroot
!@sum AND growth respiration from growth and tissue turnnover
      use ent_const
      use ent_prescr_veg
      use cohorts, only : cohort_carbon
      use patches, only : patch_carbon
      implicit none
      type(ent_config) :: config 
      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      integer :: pft
      logical :: woody
      logical :: is_annual 
      logical :: par_limit
      real*8 :: C_fol_old,C_froot_old,C_sw_old,C_hw_old,C_croot_old
      real*8 :: C_fol, C_froot, C_croot, C_sw, C_hw
      real*8 :: C_lab
      ! phenofactor phenological elongation factor [0,1] (unitless)
      real*8 :: phenofactor
      ! Cactive active carbon pool: foliage, sapwood, fine root (gC/pool/individual)
      real*8 :: Cactive
      ! Cactive_max maximum active carbon pool allowed by the allometric constraint 
      real*8 :: Cactive_max
      ! Cdead dead carbon pool, including hardwood and coarse root (gC/pool/individual)
      real*8 :: Cdead
      ! qsw - allometric factor for sapwood as a function of 
      real*8 :: qsw
      ! dbh diameter at the breast height (cm)
      real*8 :: dbh
      ! h plant height (m)
      real*8 :: h
      ! plant population (#-individual/m2-ground)
      real*8 :: nplant
      real*8 :: alloc,ialloc
      real*8 :: senescefrac !This is now net fraction of foliage that is litter.
      ! CB_d daily carbon balance (gC/individual)
      real*8 :: CB_d
      real*8 :: laipatch
      real*8 :: qf
      real*8 :: dCrepro, dC_lab
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      real*8 :: resp_auto_patch, resp_root_patch !kg-C/m/s
      integer :: cohortnum
      real*8 :: C_lab_old, Cactive_old
      real*8 :: turn_leaf,resp_growth1, resp_growth2
      logical :: dormant
      real*8 :: dC_litter_hw,dC_litter_croot
      real*8 :: phenofactor_old, alloc_adj
      real*8 :: Cfol_half
      real*8 :: cpool(N_BPOOLS) 
      logical, parameter :: alloc_new=.false.
      integer, parameter :: irecruit=1
      real*8 :: dummy

!!! debug

      real*8 :: tot_c_old, tot_c
      real*8 :: patch_tot_c_old, patch_tot_c
      real*8 :: d_tot_c(20), tot_closs(20)
      real*8 :: tot_closs_acc, tot_closs_acc_old
      real*8 :: cop_n_old, cop_n

      !Initialize
      laipatch = 0.d0
      Clossacc(:,:,:) = 0.d0 
      resp_auto_patch = 0.d0
      resp_root_patch = 0.d0
      cohortnum = 0
      resp_growth1=0.d0
      resp_growth2=0.d0
      cpool(:) = 0.d0

      patch_tot_c_old = patch_carbon(pp)

      tot_closs_acc_old = 0.d0
      tot_closs_acc = 0.d0

      cop => pp%tallest

      do while(ASSOCIATED(cop))

        tot_c_old = cop%n*cohort_carbon(cop)
        cop_n_old = cop%n
               
cddd         C_fol_old = cop%C_fol
cddd         C_froot_old = cop%C_froot
cddd         C_croot_old = cop%C_croot
cddd         C_sw_old = cop%C_sw
cddd         C_hw_old = cop%C_hw

         pft = cop%pft
         phenofactor = cop%phenofactor        

         cohortnum = cohortnum + 1

         is_annual = .false.

cddd         is_annual = (pfpar(pft)%phenotype .eq. ANNUAL)
cddd
cddd        if (is_annual) then
cddd            if (phenofactor .gt. 0.d0 .AND. cop%C_fol .eq. 0.d0) then
cddd               cop%h = 0.05d0 !min. height = 0.05m
cddd	       select case (irecruit)
cddd	       case(1)
cddd               call recruit_annual(cop%pptr%Reproduction(pft),cop%h,
cddd     o              cop%C_fol, cop%C_froot,cop%n, cop%LAI)    
cddd               cop%pptr%Reproduction(pft) = 0.d0
cddd               case(2)
cddd               cop%C_fol = height2Cfol(pft,cop%h)
cddd!               cop%LAI=cop%n*pfpar(pft)%sla*
cddd               cop%LAI=cop%n*sla(pft,cop%llspan)*
cddd     &              (height2Cfol(pft,cop%h)/1000.0d0) 
cddd               cop%C_froot = q*cop%C_fol
cddd               end select 
cddd           end if
cddd         end if

 
         dbh = cop%dbh
         h = cop%h
         nplant = cop%n
         woody = pfpar(pft)%woody
         C_lab = cop%C_lab 
         C_fol = cop%C_fol
         C_froot = cop%C_froot
         C_sw = cop%C_sw
         C_hw = cop%C_hw
         C_croot = cop%C_croot 


        !------------------------------------------------
        !*calculate allometric relation - qsw, qf, ialloc
        !------------------------------------------------ 
         !qsw  = pfpar(pft)%sla*iqsw
         qsw  = sla(pft,cop%llspan)*iqsw !qsw*h: ratio of sapwood to leaf biomass
         if (.not.woody) qsw = 0.0d0 !for herbaceous, no allocation to the wood        
         
         qf = q !q: ratio of root to leaf biomass
         if (is_annual .and. pfpar(pft)%leaftype.eq.MONOCOT)
     &      qf=q*phenofactor !for annual grasses, fine roots are prop. to the foliage
         
         alloc = phenofactor+qf+h*qsw
         if (alloc_new) alloc=1.d0+qf+h*qsw
         if (alloc .ne. 0.0d0) then
           ialloc = 1.d0/alloc
         else
           ialloc = 0.d0
         end if

        !---------------------------------------- 
        !*determine plant carbon pools 
        !---------------------------------------- 
         Cactive = C_froot + C_fol + C_sw
         if (alloc_new) then
            phenofactor_old=C_fol/C_froot*qf
            alloc_adj=alloc-1.d0+phenofactor_old
            if (alloc_adj.ne.0.d0) then
               Cactive=(C_froot+C_sw+C_fol)*alloc/alloc_adj
            else
               Cactive=C_froot+C_sw+C_fol
            end if
         end if
         Cdead = C_hw + C_croot  
         C_fol_old = C_fol
         C_froot_old = C_froot
         C_croot_old = C_croot
         C_sw_old = C_sw
         C_hw_old = C_hw
         Cactive_old =Cactive
         

         !*calculate the litter from turnover
         call litter_turnover_cohort(SDAY,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc,
     &        turn_leaf,resp_growth1)

         C_lab_old =C_lab
         C_lab = cop%C_lab 
         CB_d = cop%CB_d - (C_lab_old-C_lab)

         if (woody) then  
            Cactive_max=dbh2Cfol(pft,dbh)*(alloc+(1.d0-phenofactor))
         else
            Cactive_max=height2Cfol(pft,5.d0)*(alloc+(1.d0-phenofactor))
           !no allometric constraints in the carbon allocation, 
           !then 5.0d0 is arbitrary number (must be tall enough, then 5m)
         end if
         if (alloc_new) then
            if (woody) then  
               Cactive_max=dbh2Cfol(pft,dbh)*alloc
            else
               Cactive_max=height2Cfol(pft,5.d0)*alloc
            end if
         end if

         call prescr_init_Clab(pft,nplant,cpool)
         Cfol_half =cpool(LABILE)

         !----------------------------------------------------
         !*active growth: increment Cactive and decrease C_lab
         !----------------------------------------------------
         call  growth_cpools_active(pft,phenofactor,ialloc,Cactive_max,
     &        Cfol_half,CB_d,Cactive,C_lab,C_fol)
      
         !----------------------------------------------------  
         !*update the active pools
         !----------------------------------------------------
         cop%C_fol = phenofactor * Cactive *ialloc
         cop%C_froot = Cactive * qf * ialloc
         cop%C_sw = Cactive * h *qsw * ialloc
 
         !--------------------------------------------------------
         !*structural (corresponding active, reproductive) growth
         !-------------------------------------------------------
!          print*,pft,pfpar(pft)%phenotype, COLDDECID,pfpar(pft)%woody
!          print*,cop%phenostatus
         dormant = .false.
         dormant =
     &      (pfpar(pft)%phenotype .eq. COLDDECID .and. 
     &      pfpar(pft)%woody .and. 
     &      cop%phenostatus .lt. 2.d0 .and. cop%phenostatus .ge. 3.d0 )
     
         dCrepro = 0.d0
!#ifdef COMMENT_OUT
         if (.not.dormant)           
     &       call growth_cpools_structural(pft,dbh,h,qsw,qf,phenofactor,
     &       C_sw,Cactive_max,C_fol,CB_d,Cactive_old, 
     &       Cactive,C_lab,Cdead,dCrepro) 
!#endif
          cop%pptr%Reproduction(cop%pft) = 
     &        cop%pptr%Reproduction(cop%pft)+ dCrepro*cop%n

         
         !----------------------------------------------------  
         !*update the active and structural pools
         !----------------------------------------------------
         cop%C_fol = phenofactor * Cactive *ialloc
         cop%C_froot = Cactive * qf * ialloc
         cop%C_sw = Cactive * h *qsw * ialloc
             
         if (.not.config%do_structuralgrowth) then
             dC_litter_hw = max(0.d0,cop%C_hw - C_hw_old)
             dC_litter_croot = max(0.d0,cop%C_croot - C_croot_old)
             cop%C_hw=C_hw_old
             cop%C_croot=C_croot_old
             Cdead = cop%C_hw+cop%C_croot
         else
            dC_litter_hw =0.d0
            dC_litter_croot = 0.d0
            cop%C_hw = Cdead * hw_fract
            cop%C_croot = Cdead * (1-hw_fract)     
         endif

         !----------------------------------------------------   
         !*senesce and accumulate litter
         !---------------------------------------------------- 
         !phenology + turnover + C_lab change + growth respiration
         !senescefrac returned is fraction of foliage that is litter.

!!! HACK !!!
! just to keep the things going reset negavive pools to zero
         cop%C_fol = max( 0.d0, cop%C_fol)
         cop%C_froot = max( 0.d0, cop%C_froot)
         cop%C_croot = max( 0.d0, cop%C_croot)
         cop%C_sw = max( 0.d0, cop%C_sw)
         cop%C_hw = max( 0.d0, cop%C_hw)

         !*calculate the litter from growth
         call litter_growth_cohort(dCrepro, !SDAY,dCrepro,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        dC_litter_hw,dC_litter_croot,cop,Clossacc)

         if (C_fol_old.eq.0.d0) then
           cop%senescefrac = 0.d0
         else
           cop%senescefrac = l_fract *
     &          (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
         endif

         !* Tissue growth respiration is subtracted at physical time step in canopyspitters.f.
         !  -Update of C_growth and C_growth_flux is now done inside litter_growth_cohort
         !   as in litter_cohort.
         !cop%C_growth = (resp_growth1+resp_growth2)*cop%n*1.d-3  
         !cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux

         !Put senesced amount into litterfall into the soil.


         !----------------------------------------------------  
         !*update the plant size, LAI & nitrogen
         !----------------------------------------------------  
         if (woody) then
            cop%dbh = Cdead2dbh(pft,Cdead)
            cop%h = dbh2height(pft,cop%dbh)
         else
            cop%dbh = 0.0d0
            cop%h = Cfol2height(pft,cop%C_fol)
         end if
         
         !cop%LAI=cop%C_fol/1000.0d0*pfpar(pft)%sla*cop%n
         cop%LAI=cop%C_fol/1000.0d0*sla(pft,cop%llspan)*cop%n
         if (cop%LAI .lt. EPS) cop%LAI=EPS
         laipatch = laipatch + cop%lai  

         cop%Ntot = cop%nm * cop%LAI 

         !* Summarize for patch level *!
         !Total respiration flux including growth increment.
         resp_auto_patch = resp_auto_patch  + cop%R_auto 
         resp_root_patch = resp_root_patch + cop%R_root

#ifdef PHENOLOGY_DIAG
         call phenology_diag(cohortnum,cop)  
#endif
            
         !zero-out the daily accumulated carbon 
         cop%CB_d = 0.d0   

         d_tot_c(cohortnum) = cop%n*cohort_carbon(cop) - tot_c_old
         tot_closs_acc = tot_closs_acc + Clossacc(CARBON,LEAF,1)
     &       +Clossacc(CARBON,FROOT,1)
     &       + Clossacc(CARBON,WOOD,1)
         tot_closs(cohortnum) = tot_closs_acc - tot_closs_acc_old
         tot_closs_acc_old = tot_closs_acc
         cop_n = cop%n

cddd         write(901,*)  
cddd         write(901,*) "pft ", cop%pft
cddd         write(901,*) "deltaC ", tot_c - tot_c_old, tot_c_old
cddd         write(901,*) "deltaC*n ", (tot_c - tot_c_old)*cop%n
         !write(901,*) "Clossacc ", Clossacc


         cop => cop%shorter 
      end do !looping through cohorts
  
      !*Update Tpool from all litter.
      call litter_patch(pp, Clossacc) 

      !*Update patch-level LAI
      pp%LAI = laipatch  

      !* Update patch fluxes with growth respiration. *!
      !* The daily respiration fluxes are accumulated in C_growth to
      !* be distributed over the course of a day to avoid pulses at night.
!WRONG      pp%R_auto = resp_auto_patch !Total flux including growth increment.
      pp%R_root = pp%R_root + resp_root_patch !Total flux including growth increment.
!      pp%NPP = pp%GPP - resp_auto_patch ##Distribute with C_growth


      patch_tot_c = patch_carbon(pp)

      if( abs(patch_tot_c - patch_tot_c_old) > 1d-10 ) then
        write(903,*) "P ",patch_tot_c - patch_tot_c_old, patch_tot_c_old
     &      ,"C ",d_tot_c(1:cohortnum)*1000, "S ",tot_closs(1:cohortnum)
     &       ,"N ", cop_n_old, cop_n
      endif

      end subroutine veg_update

      !*********************************************************************
      subroutine growth_cpools_active(pft,phenofactor,ialloc, 
     &     Cactive_max,Cfol_half,CB_d,Cactive,C_lab,C_fol)
     
      integer, intent(in) :: pft
      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: ialloc
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: Cfol_half
      real*8, intent(in) :: CB_d
      real*8, intent(inout) :: Cactive
      real*8, intent(inout) :: C_lab
      real*8, intent(inout) :: C_fol
      real*8 :: Cactive_pot
      real*8 ::dC_lab  !g-C/individual. Negative for reduction of C_lab for growth.
      real*8 :: dC_remainder
      real*8 :: dCactive
      real*8 :: C_labavail
      real*8 :: dCavail
      !option for active growth
      !1 grass growth with no storage; 2 grass growth with storage; 
      !3 grass growth with storage after the certain size; 4 grass/tree growth with storage
      integer, parameter :: AGrowthModel= 3 

      !-------------------------------
      !*calculate the change in C_lab 
      !-------------------------------    
                   
      if (C_lab .gt.0.d0 .and. CB_d .gt. 0.d0) then
         if (phenofactor .eq. 0.d0) then
            dC_lab = 0.d0 !store the carbon in the labile
            dCactive = 0.d0
         else 
            !Cactive_max (max. allowed pool size according to the DBH)
            !Cactive_pot (current size + daily accumulated carbon)  
            !Cactive (cuurent size)
            Cactive_pot = Cactive + min(C_lab, CB_d)!only new carbon is used for growth.
            dCavail = min(Cactive_max, Cactive_pot) - Cactive
            select case (AGrowthModel)
            case(1) !no storage - default
               dCactive = dCavail
            case(2) !grass storage
               if (.not.pfpar(pft)%woody) then !herbaceous
                  dCactive = (1.d0-r_fract) * dCavail 
               else !woody
                  dCactive = dCavail
               end if
            case(3) !grass storage
               if (.not.pfpar(pft)%woody) then !herbaceous
                  if (C_fol .gt. Cfol_half) then
                     dCactive = (1.d0-r_fract) * dCavail 
                  else
                     dCactive = dCavail
                  end if
               else !woody
                  dCactive = dCavail
               end if
            case(4) !grass/tree storage
               if (.not.pfpar(pft)%woody) then !herbaceous
                  dCactive = (1.d0-r_fract) * dCavail 
               else !woody
                  dCactive = (1.d0-r_fract) * dCavail
               end if
            end select
            if (dCactive .lt. 0.d0)then
                dC_lab = - dCactive * l_fract
            else
                dC_lab = - dCactive
            end if
         end if
      else if (C_lab .lt. 0.d0 ) then
         dCactive = C_lab / l_fract
         dC_lab = - C_lab 
      else  
         dCactive = 0.d0
         dC_lab = 0.d0
      end if

#ifdef DEBUG
      write(200,'(100(1pe16.8))') CB_d,C_lab,dC_lab,dCactive,Cactive, 
     &                            Cactive_max,Cactive_pot
#endif

      !------------------------
      !*update the carbon pools
      !------------------------
      C_lab = C_lab + dC_lab
      Cactive = Cactive +dCactive
c$$$      C_fol = phenofactor * Cactive * ialloc
c$$$      dC_remainder = (1.d0-phenofactor )*Cactive*ialloc 
c$$$      Cactive = Cactive + dC_remainder

      end subroutine growth_cpools_active

      !*********************************************************************
      subroutine growth_cpools_structural(pft,dbh,h,qsw,qf,phenofactor,
     &      C_sw,Cactive_max,C_fol,CB_d,Cactive_old,Cactive,
     &      C_lab, Cdead, dCrepro)

      use ent_prescr_veg

      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: h
      real*8, intent(in) :: qsw
      real*8, intent(in) :: qf
      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: C_fol
      real*8, intent(in) :: CB_d
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: C_sw
      real*8, intent(inout) :: Cactive,Cactive_old
      real*8, intent(inout) :: C_lab
      real*8, intent(inout) :: Cdead
      real*8, intent(out) :: dCrepro
      real*8 :: dCdead 
      real*8 :: dCactive
      !gr_fract: fraction of excess c going to structural growth
      real*8 :: gr_fract  
      real*8 :: rp_fract
      real*8 :: qsprime,qs
      real*8 :: dCfoldCdead
      real*8 :: dCfrootdCdead
      real*8 :: dHdCdead
      real*8 :: dCswdCdead
      !option for structural growth
      !1 - based on ED1; 2 - based on ED2; 3 - to reserve Clab (not yet implemented) 
      integer, parameter :: SGrowthModel=1
      real*8 :: Cavail

      !--------------------------------------------------
      !*calculate the growth fraction for different pools
      !--------------------------------------------------        
      if (.not.pfpar(pft)%woody) then !herbaceous
         Cavail = C_lab
         if (C_lab .gt. 0.d0 )then
         if (phenofactor .eq. 0.d0) then
            qs = 0.d0  !no structural pools
            rp_fract = r_fract + c_fract
            gr_fract = 1.d0 - rp_fract
         else
            qs = 0.d0
            rp_fract = 0.d0
            gr_fract = 0.d0
         end if
         else
            qs = 0.d0
            rp_fract = 0.d0
            gr_fract = 1.d0 / l_fract
         end if
      else !woody
         Cavail = min(C_lab,C_sw) !C used for growth is limited by both size of avaiable labile storage 
                                  !& size of sapwood pool.
         if (SGrowthModel.eq.1) then !based on ED1
            if (C_fol .gt. 0.d0 .and. 
     &           Cactive .ge. Cactive_max .and. C_lab .gt. 0.d0) then
               if (dbh .le. maxdbh(pft))then
                  dCfoldCdead = dDBHdCdead(pft,Cdead)
     &                          /dDBHdCfol(pft,C_fol)
                  dCfrootdCdead = qf *dCfoldCdead
                  dHdCdead = dHdDBH(pft, dbh)  * dDBHdCdead(pft,Cdead) 
                  dCswdCdead = qsw*
     &                 (h*dCfoldCdead + C_fol*dHdCdead)
                  qsprime
     &               = 1.d0 / (dCfoldCdead + dCfrootdCdead +
     &                 dCswdCdead)
                  qs=qsprime/(1.d0+qsprime)
                  rp_fract = r_fract
                  gr_fract = 1.d0 - rp_fract
               else
                  qs = 1.d0
                  rp_fract = r_fract
                  gr_fract = 1.d0 - rp_fract 
               end if
            else if (C_lab .le. 0.d0)then
               qs = 0.d0
               rp_fract = 0.d0
               gr_fract = 1.d0 / l_fract
            else
               qs = 0.d0
               rp_fract = 0.d0
               gr_fract = 1.d0
            end if
         else if (SGrowthModel.eq.2) then !based on ED2
            if (C_lab .gt. 0.d0 .and. CB_d .gt.0.d0 )then
               qs = 1.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - rp_fract 
            else
               qs = 0.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - rp_fract
            end if
         else if (SGrowthModel.eq.3) then !option, reserving Clab 
                                          !not yet implemented
            if (C_lab .gt. 0.d0 .and. CB_d .gt.0.d0 )then
               qs = 1.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - rp_fract 
            else
               qs = 0.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - r_fract
            end if
         end if
      end if
       
      dCdead = gr_fract * qs  * Cavail 
      dCactive = gr_fract *(1.d0 - qs) * Cavail
      dCrepro =  rp_fract  * Cavail

!Notes, misc - NK:      
!      cop%pptr%Reproduction(cop%pft) = 
!     &     cop%pptr%Reproduction(cop%pft) + 0.2d0*cop%NPP*dtsec !(kg/m2-patch) !Reprod. fraction in ED is 0.3, in CLM-DGVM 0.1, so take avg=0.2.


#ifdef DEBUG
      write(201,'(100(1pe16.8))') C_lab, dCactive, dCdead,dCrepro
#endif

      !------------------------
      !*update the carbon pools
      !------------------------
      C_lab = C_lab - (dCdead + dCactive+dCrepro)
      Cdead = Cdead + dCdead
      Cactive = Cactive + dCactive

      end subroutine growth_cpools_structural
!*************************************************************************
      subroutine recruit_annual(reproduction,height,
     &           C_fol,C_froot,nplant,lai)
      real*8, intent(in) :: reproduction
      real*8, intent(in) :: height
      real*8, intent(out) :: C_fol
      real*8, intent(out) :: C_froot
      real*8, intent(out) :: nplant
      real*8, intent(out) :: lai
      real*8 :: recruit
      integer :: pft
      
      !this subroutine is sepcifically written for C3 annual grass temporarily
      pft=GRASSC3
      
      !amount of carbon, used for recruit
      recruit = reproduction * (1.d0 - mort_seedling) ! per patch-area

      !properties for new seedling
      C_fol=height2Cfol(pft,height)
      C_froot = q *C_fol
      nplant = recruit / (C_fol + C_froot)
      LAI = nplant*pfpar(pft)%sla*(C_fol/1000.d0)
      
      end subroutine recruit_annual
!*************************************************************************

      subroutine update_plant_cpools(pft, lai,h,dbh,popdens,cpool )
      integer,intent(in) :: pft !plant functional type
      real*8, intent(in) :: lai,h,dbh,popdens  !lai, h(m), dbh(cm),popd(#/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      real*8 :: qsw

      !* Initialize
      cpool(:) = 0.d0 

      qsw = pfpar(pft)%sla*iqsw
      if (.not.pfpar(pft)%woody) qsw = 0.0d0

      cpool(FOL) = lai/pfpar(pft)%sla/popdens *1d3!Bl
      cpool(FR) = q * cpool(FOL)   !Br
      cpool(SW) = h * qsw* cpool(FOL)
      if (pfpar(pft)%woody) then !Woody
        cpool(HW) = dbh2Cdead(pft,dbh) * hw_fract
        cpool(CR) = cpool(HW) * (1-hw_fract)/hw_fract !=dbh2Cdead*(1-hw_fract)
      else
        cpool(HW) = 0.d0
        cpool(CR) = 0.d0
      endif
      
      end subroutine update_plant_cpools

 
      !*********************************************************************
      subroutine accumulate_Clossacc(pft,Closs, Clossacc)
      integer, intent(in) :: pft
      real*8,intent(in) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort by depth.
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !---Local-----
      integer :: i

      !loop through CASA layers-->cumul litter per pool per layer -PK
      do i=1,N_CASA_LAYERS

        !* Accumulate *!
        Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
     &       + Closs(CARBON,LEAF,i)
        Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i) 
     &       + Closs(CARBON,FROOT,i)
        Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i) 
     &       + Closs(CARBON,WOOD,i)
      
        !* NDEAD POOLS *!
        Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i) 
     &       + Closs(CARBON,LEAF,i) * solubfract(pft)
        Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i) 
     &       + Closs(CARBON,FROOT,i) * solubfract(pft)
        Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
     &       + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
        Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i) 
     &       + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
        Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i) 
     &       + Closs(CARBON,WOOD,i)
      end do   

      !* Return Clossacc *!
      end subroutine accumulate_Clossacc

      !*********************************************************************
      subroutine litter_turnover_cohort(dt,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc,turn_leaf,resp_growth)
!@sum litter_turnover_cohort. 
!@sum CALLED BY phenology veg_update.
!@sum DAILY TIME STEP.
!@sum Calculates at daily time step litterfall from cohort 
!@sum     to soil, tissue growth,growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use cohorts, only : calc_CASArootfrac 
      use biophysics, only: Resp_can_growth
      real*8,intent(in) :: dt !seconds, time since last call
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_hw,turn_croot, turn_live !g-C/individual
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep loss less than C_lab
      real*8 :: resp_growth,resp_growth_root !g-C/individ/ms/s
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
!      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: facclim !Frost hardiness parameter - affects turnover rates in winter.

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop,fracrootCASA)

      !* NLIVE POOLS *! 
      facclim = frost_hardiness(cop%Sacclim)
      turnoverdtleaf = facclim*cop%turnover_amp*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
!      turnoverdtleaf = facclim*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
      turnoverdtfroot = facclim*annK(pft,FROOT)*SDAY
      turnoverdtwood = 0.3d0/pfpar(pft)%lrage*
     &     (1.d0-exp(-annK(pft,WOOD)*SDAY)) !Sapwood not hardwood. 0.08d0 is a tuning factor.

      !* Turnover draws down C_lab. *!
      !* Calculate adjustment factor if loss amount is too large for C_lab.
      turn_leaf = C_fol_old * turnoverdtleaf 
      turn_froot =  C_froot_old * turnoverdtfroot

      !Wood losses:  
      turn_hw = C_hw_old * turnoverdtwood 
      turn_croot = C_croot_old * turnoverdtwood
      
      ! no turnover during the winter
      if (C_fol_old .eq. 0.d0) then
      	turn_leaf = 0.d0
      	turn_froot= 0.d0
      	turn_hw = 0.d0
      	turn_croot = 0.d0
      end if

      turn_live = turn_leaf + turn_froot
 
	
      !* Distinguish respiration from turnover vs. from new growth.
      !* ### With constant prescribed structural tissue, dC_sw=0.d0,but
      !* ### there must still be regrowth of sapwood to replace that converted
      !* ### to dead heartwood.  For a hack, turn_hw is regrown as sapwood to
      !* ### maintain a carbon balance. 
      resp_turnover = 0.16d0*turn_froot + 0.014d0*turn_leaf !Coefficients from Amthor (2000) Table 3
      resp_newgrowth = 0.16d0*(max(0.d0,turn_hw)+max(0.d0,turn_croot)) 

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dClab_dbiomass =
     &     max(0.d0,turn_hw) + max(0.d0,turn_croot)
      dC_lab = - dClab_dbiomass
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover


      !* Limit turnover litter if losses and respiration exceed C_lab.*!
      if (cop%C_lab+dC_lab-resp_turnover-resp_newgrowth.lt.0.d0) then
        if ((0.5d0*cop%C_lab -resp_newgrowth).lt.0.d0)
     &       then
          adj = 0.d0            !No turnover litter to preserve C_lab for growth.
                                !C_lab will probably go negative here, but only a short while.
!        else                    !Reduce rate of turnover litter.
        else if ((turn_leaf + turn_froot).ne.0.d0)then
          adj = (0.5d0*cop%C_lab - dClab_dbiomass - resp_newgrowth)/
     &         ((1-l_fract)*(turn_leaf + turn_froot)
     &         + resp_turnover)
        else
          adj = 1.d0
        endif
      else
        adj = 1.d0
      endif

      !* Adjust turnover losses to accommodate low C_lab. *!
      turn_leaf = adj*turn_leaf
      turn_froot = adj*turn_froot
      turn_hw = adj*turn_hw
      turn_croot = adj*turn_croot


      if (adj.lt.1.d0) then
         

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     - (max(0.d0,turn_hw)+max(0.d0,turn_croot))

      end if

      resp_growth_root = 0.16d0 * turn_froot + 0.16d0*turn_croot  
      resp_growth = resp_growth_root + 0.14d0*turn_leaf+0.16d0*turn_hw  

      !* Calculate litter from turnover and from senescence*!
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      ! Senescefrac factor can be calculated by either prescribed or prognostic phenology: ****** NYK!
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%n * (1.d0-l_fract) * turn_leaf
          Closs(CARBON,WOOD,i) = cop%n * (max(0.d0,turn_hw) +
     &         fracrootCASA(i) *max(0.d0,turn_croot))
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i) *max(0.d0,turn_croot))
        end if
        ! both layers have fine root litter 
        Closs(CARBON,FROOT,i) = cop%n * (1.d0-l_fract)
     &       * fracrootCASA(i) * turn_froot
      enddo

      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
        dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &       Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3  ! convert it to kg

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      if ( dC_total < 0.d0 .and. dC_total + cop%C_total < 0.d0 ) then
        Closs(CARBON,:,:) = Closs(CARBON,:,:)
     &       *max( 0.d0, -cop%C_total/dC_total )
        cop%C_total = min(0.d0, cop%C_total)
      else
        cop%C_total = cop%C_total + dC_total
      endif
#else
      cop%C_total = cop%C_total + dC_total
#endif

      !* Update C_lab *!
      cop%C_lab = cop%C_lab + dC_lab

      !* Return Clossacc *!
      call accumulate_Clossacc(pft, Closs, Clossacc)
!      Csum = 0.d0
!      do i=1,NPOOLS
!        Csum = Csum + Clossacc(CARBON,i,1)
!      enddo

      !* Return growth respiration/day *!
      !################ ###################################################
      !#### DUE TO TIMING OF LAI UPDATE IN GISS GCM AT THE DAILY TIME STEP,
      !#### GROWTH RESPIRATION FROM CHANGE IN LAI NEEDS TO BE SAVED AS 
      !#### A RESTART VARIABLE IN ORDER TO SEND THAT FLUX TO THE ATMOSPHERE.
      !#### Igor has put in code to distribute C_growth over the day.
      !####################################################################
      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.
      cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3  !kg-C m-2 day-1
      cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, kg-C m-2 s-1

      end subroutine litter_turnover_cohort
      !*********************************************************************
      subroutine litter_growth_cohort(dCrepro,!dt,dCrepro,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        dC_litter_hw,dC_litter_croot,cop,Clossacc)
!@sum litter_cohort for prognostic growth.
!@sum CALLED BY phenology veg_update.
!@sum DAILY TIME STEP.
!@sum Calculates litterfall from cohort 
!@sum     to soil, tissue growth,growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use cohorts, only : calc_CASArootfrac 
      use biophysics, only: Resp_can_growth
      !real*8,intent(in) :: dt   !seconds, time since last call
      real*8,intent(in) ::dCrepro
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      real*8, intent(in) :: dC_litter_hw, dC_litter_croot
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_hw,turn_croot, turn_live !g-C/individual
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep loss less than C_lab
      real*8 :: resp_growth,resp_growth_root !g-C/individ (mass total per day)
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: facclim !Frost hardiness parameter - affects turnover rates in winter.
      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop,fracrootCASA)
!      print *, 'from litter(pheno*.f): fracrootCASA(:) =', fracrootCASA !***test*** -PK 11/27/06  

      !* NLIVE POOLS *! 
      facclim = frost_hardiness(cop%Sacclim)


      !* Change in plant tissue pools. *!
      dC_fol = cop%C_fol-C_fol_old
      dC_froot = cop%C_froot - C_froot_old
      dC_hw = cop%C_hw - C_hw_old
      dC_sw = cop%C_sw - C_sw_old
      dC_croot = cop%C_croot - C_croot_old

      !* Distinguish respiration from turnover vs. from new growth.
      !* ### With constant prescribed structural tissue, dC_sw=0.d0,but
      !* ### there must still be regrowth of sapwood to replace that converted
      !* ### to dead heartwood.  For a hack, turn_hw is regrown as sapwood to
      !* ### maintain a carbon balance. 

      !* C_lab required for biomass growth or senescence (not turnover)
      dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot)!Growth of new tissue
     &     + max(0.d0,dC_sw)    
     &     + max(0.d0,dC_hw) + max(0.d0,dC_croot)
     &     - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &     + max(0.d0,-dC_sw))

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &     - dClab_dbiomass 
!!!     removed - dCrepro since it doesn't conserve carbon       !Growth (new growth or senescence)

      !* Recalculate respiration.  Distinguish below- vs. above-ground autotrophic respiration.

      resp_growth_root = 0.16d0*(max(0.d0,dC_froot)+max(0.d0,dC_croot))!New biomass growth
      resp_growth = resp_growth_root + 0.14d0 * 
     &      (max(0.d0,dC_fol)+max(0.d0,dC_sw)) + 0.16d0* max(0.d0,dC_hw) 
     
      !* Calculate litter from turnover and from senescence*!
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      ! Senescefrac factor can be calculated by either prescribed or prognostic phenology: ****** NYK!
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%n * (1.d0-l_fract) * 
     &         max(0.d0,-dC_fol)
          Closs(CARBON,WOOD,i) = cop%n * (
     &        max(0.d0,-dC_hw)  + dC_litter_hw +
     &        fracrootCASA(i) * (max(0.d0,-dC_croot) + dC_litter_croot))
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i) * (max(0.d0,-dC_croot) + dC_litter_croot))
        end if
        ! both layers have fine root litter 
        Closs(CARBON,FROOT,i) = cop%n * (1.d0-l_fract)
     &       * fracrootCASA(i) * ( max(0.d0,-dC_froot)
     &       + max(0.d0,-dC_sw) )
      enddo

      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
        dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &       Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3  ! convert it to kg

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      if ( dC_total < 0.d0 .and. dC_total + cop%C_total < 0.d0 ) then
        Closs(CARBON,:,:) = Closs(CARBON,:,:)
     &       *max( 0.d0, -cop%C_total/dC_total )
        cop%C_total = min(0.d0, cop%C_total)
      else
        cop%C_total = cop%C_total + dC_total
      endif
#else
      cop%C_total = cop%C_total + dC_total
#endif

      call accumulate_Clossacc(pft, Closs, Clossacc)

      !Error check
      if ( abs( (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
     &     +dC_lab)*cop%n + Closs(CARBON,LEAF,1)+Closs(CARBON,FROOT,1)
     &     + Closs(CARBON,WOOD,1) ) > 1d-10 ) then
        
         write(901,*) "Closs ", Closs(CARBON,LEAF,1)
     &       +Closs(CARBON,FROOT,1)
     &       + Closs(CARBON,WOOD,1)
     &       ,"dC ", dC_fol, dC_froot, dC_hw, dC_sw, dC_croot
     &       ,dC_lab
     &       ,"dC*n ", (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
     &       +dC_lab)*cop%n
     &       ,"dCrepro ", dCrepro, dCrepro*cop%n
     &       ,"dC_litter_hw ", dC_litter_hw, dC_litter_hw*cop%n
     &       ,"dC_litter_croot ",
     &       dC_litter_croot, dC_litter_croot*cop%n

       endif

cddd      write(901,*) "Closs ", Closs(CARBON,LEAF,1)+Closs(CARBON,FROOT,1)
cddd     &     + Closs(CARBON,WOOD,1)
cddd      write(901,*) "dC ", dC_fol, dC_froot, dC_hw, dC_sw, dC_croot
cddd     &     ,dC_lab
cddd      write(901,*) "dC*n ", (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
cddd     &     +dC_lab)*cop%n
cddd      write(901,*) "dCrepro ", dCrepro, dCrepro*cop%n
cddd      write(901,*) "dC_litter_hw ", dC_litter_hw, dC_litter_hw*cop%n
cddd      write(901,*) "dC_litter_croot ",
cddd     &     dC_litter_croot, dC_litter_croot*cop%n

!      write(992,*) C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
!     &     cop%C_lab,cop%C_fol,cop%C_froot,cop%C_hw,cop%C_sw,
!     &     cop%C_croot, cop%dbh,turn_leaf,turn_froot,turn_hw,turn_croot,
!     &     dC_fol,dC_froot,dC_hw,dC_sw,dC_croot,
!     &     Closs(CARBON,:,:), Clossacc(CARBON,:,:),adj,cop%turnover_amp,
!     &     facclim,turnoverdtleaf,turnoverdtfroot, turnoverdtwood

       !* Update C_lab *!
      cop%C_lab = cop%C_lab + dC_lab
      !at this point, C_lab<0 comes from the rounding errors...
#ifdef COMMENT_OUT
      if (cop%C_lab < 0.d0) cop%C_lab = 0.d0
#endif
      if (cop%C_lab < -1.d-8) then
        write(902,*) "WARNING: Clab ", cop%C_lab, dC_lab, cop%pft
     &       ,"dC ", dC_fol, dC_froot, dC_hw, dC_sw, dC_croot
     &       ,dC_lab
     &       ,"dC*n ", (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
     &       +dC_lab)*cop%n
     &       ,"dCrepro ", dCrepro, dCrepro*cop%n
     &       ,"dC_litter_hw ", dC_litter_hw, dC_litter_hw*cop%n
     &       ,"dC_litter_croot ",
     &       dC_litter_croot, dC_litter_croot*cop%n
      endif
!      if (cop%C_lab < 0.d0) then
!         print*,dC_fol,cop%C_fol,dC_sw,cop%C_sw
!         print*,dC_lab,cop%C_lab, dC_froot, cop%C_froot
!         print*,dC_hw,cop%C_hw,dC_croot,cop%C_croot
!         stop
!      endif

      !* Return Clossacc and resp_growth *!
!      Csum = 0.d0
!      do i=1,NPOOLS
!        Csum = Csum + Clossacc(CARBON,i,1)
!      enddo

      !* Return resp_growth in cop%C_growth *!
      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.
      cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3  !kg-C m-2 day-1
      cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, kg-C m-2 s-1

      end subroutine litter_growth_cohort
      !*********************************************************************
      subroutine litter_cohort(dt,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc)
!@sum litter_cohort for static woody structure 
!@sum CALLED BY ent_prescribed_updates.
!@sum DAILY TIME STEP.
!@sum Calculates litterfall from cohort to soil, tissue growth,
!@sum growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use cohorts, only : calc_CASArootfrac 
      use biophysics, only: Resp_can_growth
      real*8,intent(in) :: dt !seconds, time since last call
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_hw,turn_croot, turn_live !g-C/individual
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep turnover less than C_lab
      real*8 :: resp_growth,resp_growth_root !g-C/individ/ms/s
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: facclim !Frost hardiness parameter - affects turnover rates in winter.

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop,fracrootCASA)

      !* NLIVE POOLS *! 
      facclim = frost_hardiness(cop%Sacclim)
      turnoverdtleaf = facclim*cop%turnover_amp*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
!      turnoverdtleaf = facclim*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
      turnoverdtfroot = facclim*annK(pft,FROOT)*SDAY
      turnoverdtwood = 0.32d0/pfpar(pft)%lrage*
     &     (1.d0-exp(-annK(pft,WOOD)*SDAY)) !Sapwood not hardwood.  0.08d0 is a tuning factor.

      !* Turnover draws down C_lab. *!
      !* Calculate adjustment factor if loss amount is too large for C_lab.
      turn_leaf = C_fol_old * turnoverdtleaf 
      turn_froot =  C_froot_old * turnoverdtfroot
      turn_live = turn_leaf + turn_froot

      !Wood losses:  
      turn_hw = C_hw_old * turnoverdtwood 
      turn_croot = C_croot_old * turnoverdtwood
      !## Need to add hw turn corresponding to woody allocation that goes to litter for static woody structure.

      !* Change in plant tissue pools. *!
      dC_fol = cop%C_fol-C_fol_old
      dC_froot = cop%C_froot - C_froot_old
      dC_hw = cop%C_hw - C_hw_old
      dC_sw = cop%C_sw - C_sw_old
      dC_croot = cop%C_croot - C_croot_old

      !* Distinguish respiration from turnover vs. from new growth.
      !* ### With constant prescribed structural tissue, dC_sw=0.d0,but
      !* ### there must still be regrowth of sapwood to replace that converted
      !* ### to dead heartwood.  For a hack, turn_hw is regrown as sapwood to
      !* ### maintain a carbon balance. 
      resp_turnover = 0.16d0*turn_froot + 0.014d0*turn_leaf !Coefficients from Amthor (2000) Table 3
      resp_newgrowth = 0.16d0*max(0.d0,dC_froot) + 
     &     0.14d0*(max(0.d0,dC_fol)+max(0.d0,dC_sw))
     &     +0.16d0*(max(0.d0,turn_hw)+max(0.d0,turn_croot)) !##THIS IS RESPIRATION FOR REGROWTH OF SAPWOOD TO ACCOUNT FOR CONVERSION TO HEARTWOOD WITH CONSTANT PLANT STRUCTURE.

      !* C_lab required for biomass growth or senescence (not turnover)
      dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot) !Growth of new tissue
     &     + max(0.d0,dC_sw)    !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &     + max(0.d0,turn_hw)+max(0.d0,turn_croot) !### Sapwood growth to replace that converted to heartwood.
     &     - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &     + max(0.d0,-dC_sw))

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     - dClab_dbiomass       !Growth (new growth or senescence)
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

      !* Limit turnover litter if losses and respiration exceed C_lab.*!
      if (cop%C_lab+dC_lab-resp_turnover-resp_newgrowth.lt.0.d0) then
        if ((0.5d0*cop%C_lab - dClab_dbiomass-resp_newgrowth).lt.0.d0)
     &       then
          adj = 0.d0            !No turnover litter to preserve C_lab for growth.
                                !C_lab will probably go negative here, but only a short while.
        else                    !Reduce rate of turnover litter.
          adj = (0.5d0*cop%C_lab - dClab_dbiomass - resp_newgrowth)/
     &         ((1-l_fract)*(turn_leaf + turn_froot)
     &         + resp_turnover)
        endif
      else
        adj = 1.d0
      endif

      !* Adjust turnover losses to accommodate low C_lab. *!
      if (adj < 1.d0) then
        turn_leaf = adj*turn_leaf
        turn_froot = adj*turn_froot
        turn_hw = adj*turn_hw
        turn_croot = adj*turn_croot
        
        !* Recalculate dClab_dbiomass *!
        dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot) !Growth of new tissue
     &       + max(0.d0,dC_sw)  !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &       + max(0.d0,turn_hw)+max(0.d0,turn_croot) !### This is sapwood growth to replace that converted to heartwood.
     &       - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &       + max(0.d0,-dC_sw))
      endif

      !* Recalculate respiration. 
      !  Distinguish below- vs. above-ground autotrophic respiration.
        resp_growth_root = 0.16d0 * ( !Coefficient from Amthor (2000) Table 3
     &       turn_froot            !Turnover growth
     &       + max(0.d0,dC_froot)) !New biomass growth
     &       + 0.16d0*turn_croot   !# Hack for regrowth of sapwood converted to replace senesced coarse root.
        resp_growth = resp_growth_root + 
     &       0.14d0 *              !Coefficient from Amthor (2000) Table 3
     &       ( turn_leaf           !Turnover growth    
     &       +max(0.d0,dC_fol)+max(0.d0,dC_sw)) !New biomass growth
     &       + 0.16d0*turn_hw      !# Hack for regrowth of sapwood converted to replace senesced hw.

!      write(991,*)  facclim,turn_froot,turn_croot,max(0.d0,dC_froot),
!     &     max(0.d0,dC_croot), turn_leaf,turn_hw,
!     &     max(0.d0,dC_fol), max(0.d0,dC_sw) !New biomass growth

      !* Recalculate dC_lab in case adj < 1.0.
      dC_lab = 
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     - dClab_dbiomass       !Growth (new growth or senescence)
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

      !* Calculate litter to soil from turnover and from senescence*!
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      ! Senescefrac factor diagnostic also calculated.
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%n * (1.d0-l_fract) * (turn_leaf +
     &         max(0.d0,-dC_fol))
          Closs(CARBON,WOOD,i) = cop%n * (turn_hw + 
     &         max(0.d0,-dC_hw) +
     &         fracrootCASA(i)*
     &         (turn_croot+max(0.d0,-dC_croot)))
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i)
     &         *(turn_croot+max(0.d0,-dC_croot)))
        end if
        ! both layers have fine root litter 
        Closs(CARBON,FROOT,i) = cop%n * (1.d0-l_fract)
     &       * fracrootCASA(i) 
     &       * (turn_froot + max(0.d0,-dC_froot))
      enddo

      !* Diagnostic
      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
        dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &       Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3  ! convert it to kg

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      if ( dC_total < 0.d0 .and. dC_total + cop%C_total < 0.d0 ) then
        Closs(CARBON,:,:) = Closs(CARBON,:,:)
     &       *max( 0.d0, -cop%C_total/dC_total )
        cop%C_total = min(0.d0, cop%C_total)
      else
        cop%C_total = cop%C_total + dC_total
      endif
#else
      cop%C_total = cop%C_total + dC_total
#endif

      call accumulate_Clossacc(pft,Closs, Clossacc)

!      write(992,*) C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
!     &     cop%C_lab,cop%C_fol,cop%C_froot,cop%C_hw,cop%C_sw,
!     &     cop%C_croot, cop%dbh,turn_leaf,turn_froot,turn_hw,turn_croot,
!     &     dC_fol,dC_froot,dC_hw,dC_sw,dC_croot,
!     &     Closs(CARBON,:,:), Clossacc(CARBON,:,:),adj,cop%turnover_amp,
!     &     facclim,turnoverdtleaf,turnoverdtfroot, turnoverdtwood

      !################ ###################################################
      !#### DUE TO TIMING OF LAI UPDATE IN GISS GCM AT THE DAILY TIME STEP,
      !#### GROWTH RESPIRATION FROM CHANGE IN LAI NEEDS TO BE SAVED AS 
      !#### A RESTART VARIABLE IN ORDER TO SEND THAT FLUX TO THE ATMOSPHERE.
      !#### Igor has put in code to distribute C_growth over the day.
      !####################################################################

      cop%C_lab = cop%C_lab + dC_lab

      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.
      cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3
      cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, C s-1

      !Cactive = Cactive - turn_leaf - turn_froot !No change in active

      !* Diagnostic
      if ( C_fol_old > 0.d0 ) then
        cop%senescefrac = l_fract *
     &     (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
        !cop%senescefrac = max(0.d0,-dC_fol/C_fol_old)
      else
        cop%senescefrac = 0.d0
      endif

      !* Return Clossacc *!
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo

      end subroutine litter_cohort

!**********************************************************************
      subroutine litter_patch(pp, Clossacc)
!@sum litter_dead.  Update soil Tpools following litterfall.
      type(patch),pointer :: pp 
      real*8,intent(in) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !----Local------
      integer :: i

      !* NDEAD POOLS *!
       do i=1,N_CASA_LAYERS
        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
     &     + Clossacc(CARBON,SURFMET,i)
        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
     &     + Clossacc(CARBON,SOILMET,i)
        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
     &     + Clossacc(CARBON,SURFSTR,i)
        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
     &     + Clossacc(CARBON,SOILSTR,i)
        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
     &     + Clossacc(CARBON,CWD,i)
       end do   !loop through CASA layers-->total C per pool per layer -PK

       end subroutine litter_patch

!*********************************************************************
c      subroutine litter_old( pp)
c      !* Determine litter from live carbon pools and update Tpool.
c      !* After CASA, but called at daily time step. - NYK 7/27/06
c      
c      use cohorts, only : calc_CASArootfrac 
c
c      !real*8 :: dtsec           !dt in seconds
c      !type(timestruct) :: tt    !Greenwich Mean Time
c      type(patch),pointer :: pp
c      !--Local-----------------
c      type(cohort),pointer :: cop
c      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
c      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
c      integer :: pft,i
c      real*8 :: fracrootCASA(N_CASA_LAYERS)
c      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdttotal, adj !Total, adjustment factor if larger than C_lab.
c
c      Closs(:,:,:) = 0.d0
c      Clossacc(:,:,:) = 0.d0
c
c      !* Calculate fresh litter from each cohort *!
c      cop => pp%tallest  !changed to => (!) -PK 7/11/06
c      do while(ASSOCIATED(cop)) 
c        pft = cop%pft
c        
c      !assign root fractions for CASA layers -PK
c      call calc_CASArootfrac(cop,fracrootCASA)
c!      print *, 'from litter(pheno*.f): fracrootCASA(:) =', fracrootCASA !***test*** -PK 11/27/06  
c
c       do i=1,N_CASA_LAYERS  !do this over all CASA layers -PK
c        !* NLIVE POOLS *! 
c        turnoverdtleaf = annK(pft,LEAF)*SDAY
c        turnoverdtfroot = annK(pft,FROOT)*SDAY
c        turnoverdtwood = 1.d0-exp(-annK(pft,WOOD)*SDAY) !Sapwood not hardwood
c
c        !* UPDATE C_LAB: Turnover should draw down C_lab. *!
c        ! Check that amount not too large 
c        turnoverdttotal = turnoverdtleaf+turnoverdtfroot+turnoverdtwood
c        !if (turnoverdttotal.lt.cop%C_lab) then
c          adj = 1.d0  !No adjustment, enough C_lab
c        !else
c        !  adj = (cop%C_lab - EPS)/turnoverdttotal !Turnover can reduce C_lab only to EPS.
c        !endif
c        !cop%C_lab = cop%C_lab - adj*turnoverdttotal
c        !* NEED TO PUT RETRANSLOCATION IN HERE, TOO *!
c
c        !* Calculate litter *!
c        ! Senescefrac factor can be calculated by either prescribed or prognostic phenology: ****** NYK!
c        if (i.eq.1) then  !only top CASA layer has leaf and wood litter -PK   
c         Closs(CARBON,LEAF,i) = 
c     &       cop%C_fol * cop%n *
c     &       (adj*turnoverdtleaf + cop%senescefrac) !* x tune factor
c         Closs(CARBON,WOOD,i) = 
c     &       (cop%C_hw + cop%C_croot) * cop%n
c     &       *(adj*turnoverdtwood) !* expr is kdt; x tune factor
c        else    
c         Closs(CARBON,LEAF,i) = 0.d0 
c         Closs(CARBON,WOOD,i) = 0.d0
c        end if
c         Closs(CARBON,FROOT,i) =  !both layers have root litter -PK 
c     &       fracrootCASA(i)*cop%C_froot * cop%n * 
c     &       (adj*turnoverdtfroot + cop%senescefrac) !* x tune factor
c
!        write(98,*) 'In litter: ',dtsec
!        write(98,*) cop%pft, cop%C_fol, cop%n, annK(pft,LEAF),pp%betad
!        write(98,*) cop%pft, cop%C_froot, cop%n, annK(pft,FROOT)
!        write(98,*) cop%pft, cop%C_hw, cop%n, annK(pft,WOOD)
!        write(98,*) 'solubfract(pft)', solubfract(pft)
!        write(98,*) 'Closs(CARBON,LEAF)',Closs(CARBON,LEAF)
!        write(98,*) 'Closs(CARBON,FROOT)',Closs(CARBON,FROOT)
c         Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
c     &        + Closs(CARBON,LEAF,i)
c         Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i) 
c     &        + Closs(CARBON,FROOT,i)
c         Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i) 
c     &        + Closs(CARBON,WOOD,i)
c
c        !* NDEAD POOLS *!
c         Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i) 
c     &        + Closs(CARBON,LEAF,i) * solubfract(pft)
c         Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i) 
c     &        + Closs(CARBON,FROOT,i) * solubfract(pft)
c         Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
c     &        + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
c         Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i) 
c     &        + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
c         Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i) 
c     &        + Closs(CARBON,WOOD,i)
c        end do  !loop through CASA layers-->cumul litter per pool per layer -PK
c     
c        cop => cop%shorter  !added -PK 7/12/06
c      end do  !loop through cohorts

      !* NDEAD POOLS *!
c       do i=1,N_CASA_LAYERS
c        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
c     &     + Clossacc(CARBON,SURFMET,i)
c        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
c     &     + Clossacc(CARBON,SOILMET,i)
c        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
c     &     + Clossacc(CARBON,SURFSTR,i)
c        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
c     &     + Clossacc(CARBON,SOILSTR,i)
c        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
c     &     + Clossacc(CARBON,CWD,i)
c       end do   !loop through CASA layers-->total C per pool per layer -PK
c!       print *, __FILE__,__LINE__,'pp%Tpool=',pp%Tpool(CARBON,:,:) !***test*** -PK 7/24/07  
c
c      end subroutine litter_old

!*********************************************************************

      subroutine photosyn_acclim(dtsec,Ta,Sacc)
!@sum Model for state of acclimation/frost hardiness based 
!@sum for boreal coniferous forests based on Repo et al (1990), 
!@sum Hanninen & Kramer (2007),and  Makela et al (2006)
      implicit none
      real*8,intent(in) :: dtsec ! time step size [sec]
      real*8,intent(in) :: Ta ! air temperature [deg C]
      real*8,intent(inout) :: Sacc ! state of acclimation [deg C]

      !----Local-----
      real*8,parameter :: tau_inv = 2.22222e-6 
                          ! inverse of time constant of delayed
                          ! response to ambient temperature [sec] = 125 hr 
                          ! Makela et al (2004) for Scots pine

!      Use a first-order Euler scheme
       Sacc = Sacc + dtsec*(tau_inv)*(Ta - Sacc) 

!!     Predictor-corrector method requires temperature from next timestep
!       Sacc_old = Sacc
!       Sacc = Sacc_old + dtsec*(1/tau_acclim)*(Ta - Sacc ) 
!       Sacc = Sacc_old + ((1/tau_acclim)*(Ta - Sacc_old)+
!     &                     (1/tau_acclim)*(Ta_next - Sacc))*0.5d*dtsec

      end subroutine photosyn_acclim

!*************************************************************************
      real*8 function running_mean(dtsec,numd,var,var_mean) 
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: numd !number of days for running mean
      real*8, intent(in) :: var
      real*8, intent(in) :: var_mean
      real*8 :: zweight

      zweight=exp(-1.d0/(numd*86400.d0/dtsec))
      running_mean=zweight*var_mean+(1.d0-zweight)*var  
      
      end function running_mean
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
      real*8 function dbh2Cfol(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8 :: maxdbh

      maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
     &     /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht

c$$$      if (.not.pfpar(pft)%woody) then !herbaceous
c$$$         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      else !woody
c$$$         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      end if

      dbh2Cfol=1000.0d0*(1.0/C2B) *pfpar(pft)%b1Cf 
     &        * min(dbh, maxdbh)**pfpar(pft)%b2Cf

      end function dbh2Cfol
!*************************************************************************
      real*8 function maxdbh(pft)
      integer,intent(in) :: pft       

      maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
     &     /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht

c$$$      if (.not.pfpar(pft)%woody) then !grasses/crops 
c$$$         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      else !woody   
c$$$         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      end if
 
      end function maxdbh
!*************************************************************************
      real*8 function dbh2Cdead(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dbh2Cdead = (1.0/C2B) * pfpar(pft)%b1Cd * 
     &             dbh**pfpar(pft)%b2Cd * 1000.0d0

      end function  dbh2Cdead
!*************************************************************************
      real*8 function dbh2height(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dbh2height = 1.3 + pfpar(pft)%b1Ht * 
     &             (1.0-exp(pfpar(pft)%b2Ht*dbh))

      end function dbh2height
!*************************************************************************
      real*8 function height2dbh(pft,h)
      integer,intent(in) :: pft
      real*8, intent(in) :: h
      real*8 :: hcrit, hin

      hcrit=1.3d0+pfpar(pft)%b1Ht-EPS
      hin=min(h,hcrit)

      height2dbh = log(1.d0-(hin-1.3d0)/pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      

      end function height2dbh
!*************************************************************************
      real*8 function Cdead2dbh(pft,Cdead)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

      Cdead2dbh = (Cdead/1000.0d0*C2B/pfpar(pft)%b1Cd)
     &            **(1.0d0/pfpar(pft)%b2Cd)

      end function Cdead2dbh
!*************************************************************************
      real*8 function height2Cfol(pft,height)
      integer,intent(in) :: pft
      real*8, intent(in) :: height
      real*8,parameter :: h1Cf = 1.66d0 
      real*8,parameter :: h2Cf = 1.50d0 
      real*8,parameter :: nplant = 2500.d0 
      height2Cfol=(1.0d0/C2B) *h1Cf 
     &        *((height*100.0d0)**h2Cf)/nplant

      end function height2Cfol
!*************************************************************************
      real*8 function Cfol2height(pft,Cfol)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !gC/pool/plant
      real*8,parameter :: h1Cf = 1.66d0 
      real*8,parameter :: h2Cf = 1.5d0
      real*8,parameter :: nplant = 2500.d0

      if (Cfol.gt.0.0d0) then
         Cfol2height=exp(log(Cfol*C2B/h1Cf*nplant)/h2Cf)/100.0d0
      else
         Cfol2height=0.d0
      end if

      end function Cfol2height
!*************************************************************************
      real*8 function dDBHdCdead(pft,Cdead)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

      dDBHdCdead=(C2B/1000.0d0/pfpar(pft)%b1Cd)**(1.0d0/pfpar(pft)%b2Cd)
     &          *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
     &          /pfpar(pft)%b2Cd

      end function dDBHdCdead
!*************************************************************************
      real*8 function dDBHdCfol(pft,Cfol)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol

      dDBHdCfol=(C2B/1000.0d0/pfpar(pft)%b1Cf)**(1.0d0/pfpar(pft)%b2Cf)
     &          *Cfol**((1.0d0/pfpar(pft)%b2Cf)-1.0d0)
     &          /pfpar(pft)%b2Cf

      end function dDBHdCfol
!*************************************************************************
      real*8 function dHdDBH(pft,dbh)
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh

      dHdDBH = - pfpar(pft)%b1Ht*pfpar(pft)%b2Ht
     &         * exp(pfpar(pft)%b2Ht * dbh)
      end function dHdDBH

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
      real*8 function frost_hardiness(Sacclim) Result(facclim)
!@sum frost_hardiness.  Calculate factor for adjusting photosynthetic capacity
!@sum  due to frost hardiness phenology.
      real*8,intent(in) :: Sacclim 
      !----Local-----
      real*8,parameter :: Tacclim=-5.93d0 ! threshold temperature for photosynthesis [deg C]
                        ! Site specific thres. temp.: state of photosyn.acclim
                        ! Hyytiala Scots Pine, -5.93 deg C Makela et al (2006)
      real*8,parameter :: a_const=0.0595 ! factor to convert from Sacclim [degC] to facclim [-]
                        ! Site specific; conversion (1/Sacclim_max)=1/16.8115
                        ! estimated by using the max S from Hyytiala 1998
!      real*8 :: facclim ! acclimation/frost hardiness factor [-]

      if (Sacclim > Tacclim) then ! photosynthesis occurs 
         facclim = a_const * (Sacclim-Tacclim) 
         if (facclim > 1.d0) facclim = 1.d0
!      elseif (Sacclim < -1E10)then !UNDEFINED
      elseif (Sacclim.eq.UNDEF)then !UNDEFINED
         facclim = 1.d0   ! no acclimation for this pft and/or simualtion
      else
         facclim = 0.01d0 ! arbitrary min value so that photosyn /= zero
      endif

      end function frost_hardiness
    
!*************************************************************************
      subroutine phenology_diag(cohortnum, cop)
      implicit none
      type(cohort), pointer ::cop
      integer :: cohortnum
      real*8 :: fall_real
      if (cop%pptr%cellptr%fall==1) then
        fall_real=1.d0
      else
         fall_real=0.d0
      end if
         write(990,'(2(i5),27(1pe16.8))')
     &        cohortnum, !1
     &        cop%pft,  
     &        cop%phenofactor_c,
     &        cop%phenofactor_d,
     &        cop%phenostatus,
     &        cop%LAI,
     &        cop%C_fol,
     &        cop%C_lab,
     &        cop%C_sw,
     &        cop%C_hw,
     &        cop%C_froot,!11
     &        cop%C_croot,
     &        cop%NPP,
     &        cop%dbh,
     &        cop%h,
     &        cop%CB_d,
     &        cop%senescefrac,
     &        cop%llspan,
     &        cop%turnover_amp,
     &        cop%pptr%cellptr%airtemp_10d,
     &        cop%pptr%cellptr%soiltemp_10d, !21
     &        cop%betad_10d,
     &        cop%pptr%cellptr%par_10d,
     &        cop%pptr%cellptr%gdd,
     &        cop%pptr%cellptr%ncd,
     &        cop%pptr%cellptr%CosZen, 
     &        cop%pptr%cellptr%daylength(1),
     &        cop%pptr%cellptr%daylength(2),
     &        fall_real
!     &        cop%pptr%cellptr%fall

      end subroutine phenology_diag
!*************************************************************************
      end module phenology
