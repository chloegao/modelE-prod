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

      private growth_cpools_active
      private growth_cpools_structural
!      private senesce_cpools
      private recruit_annual
      private photosyn_acclim
      private phenology_diag

      !************************************************************************
      !* GROWTH MODEL CONSTANTS - phenology & carbon allocation 
      !l_fract: fraction of leaves retained after leaf fall (unitless)
      real*8, parameter :: l_fract = 0.50d0 
      !growth_r:  fraction of biomass pool required for growth respration to grow that biomass.
!      real*8, parameter :: growth_r = 0.30d0 !Check same as canopyspitters.f Respauto_NPP_Clabile
      !q: ratio of root to leaf biomass (unitless)
      real*8, parameter :: q=1.0d0 
      !iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m)
      !3900.0: leaf area per sapwood area (m2/m2) 
      !1000.0: sapwood density (kg/m3)
      !2.0:  biomass per carbon (kg/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282
      real*8, parameter :: iqsw=1000.0d0/3900.0d0/2.0d0
      !hw_fract: ratio of above ground stem to total stem (stem plus structural roots)
      real*8, parameter :: hw_fract = 0.70d0 
      !C2B: ratio of biomass to carbon (kg-Biomass/kg-Carbon)
      real*8, parameter :: C2B = 2.0d0 
      !airtemp_par
      real*8, parameter :: airtemp_par = 5.d0 
      !temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      ! gdd_par1/2/3: paramters to estimate the threshold for gdd
      ! gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: gdd_par1 = -68.d0 
      real*8, parameter :: gdd_par2 = 638.d0
      real*8, parameter :: gdd_par3 = -0.01d0 
      !gdd_length
      real*8, parameter :: gdd_length = 200.d0 !50.d0
      !airt_threshold
      real*8, parameter :: airt_max_w = 15.d0
      real*8, parameter :: airt_min_w = 5.d0
      !soilt_threshold
      real*8, parameter :: soilt_max = 10.d0
      real*8, parameter :: soilt_min = 0.d0
      !ld_threshold (minute): light length constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: ld_threshold = 655.d0
      real*8, parameter :: ld_min =540.d0
      real*8, parameter :: ld_max =550.d0
      !tsoil_threshold1, tsoil_threshold2 : soil temperature constraint for cold-deciduous woody PFTs (White et al. 1997)
!      real*8, parameter :: tsoil_threshold1 = 11.15d0
!      real*8, parameter :: tsoil_threshold2 = 2.d0
      !ddfacu: the rate of leaf fall (1/day) (IBIS)
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !paw_: water threshold (unitless)  
      !_w for woody & _h for herbaceous      
      real*8, parameter :: paw_max_w = 0.3d0
      real*8, parameter :: paw_min_w = -0.1d0
      real*8, parameter :: paw_res_w = 0.25d0
!SOILMOIST_OLD
!      real*8, parameter :: paw_max_h = 0.4d0 
!      real*8, parameter :: paw_min_h = 0.1d0 
      real*8, parameter :: paw_max_h = 0.5d0 
      real*8, parameter :: paw_min_h = 0.3d0 
      real*8, parameter :: paw_res_h = 1.0d0
      !betad : water_stress3
      real*8, parameter :: betad_max_w = 0.2d0 !0.5d0 !1.0d0 
      real*8, parameter :: betad_min_w = 0.d0 !-0.5d0 
      real*8, parameter :: betad_res_w = 0.25d0
      real*8, parameter :: betad_max_h = 1.0d0 
      real*8, parameter :: betad_min_h = 0.9d0!0.95d0 !0.8d0 
      real*8, parameter :: betad_res_h = 1.0d0      
      !light-controll: par_turnover_int & par_turnover_slope
      real*8, parameter :: par_turnover_int = -12.d0 !-10.6d0 
      real*8, parameter :: par_turnover_slope = 0.18d0  
      !r_fract: fraction of excess c going to seed reproduction 
      real*8, parameter :: r_fract = 0.3d0
      !c_fract: fraction of excess c going to clonal reproduction - only for herbaceous
      real*8, parameter :: c_fract = 0.7d0
      !mort_seedling: mortality rate for seedling
      real*8, parameter :: mort_seedling = 0.90d0 !0.80d0 !0.95d0

      contains


      !*********************************************************************
      subroutine clim_stats(dtsec, ecp, config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running average   
      use soilbgc, only : Soillayer_convert_Ent 
      real*8,intent(in) :: dtsec           !dt in seconds
      type(entcelltype) :: ecp      
      type(ent_config) :: config
      logical, intent(in) :: dailyupdate
      !-----local--------
      type(patch), pointer :: pp
      type(cohort), pointer :: cop 
      !temperature constrain for cold-deciduous PFTs (Botta et al. 1997)
      !airtemp_par: base temperature to estimate growing degree day (gdd) (degree C)
!      real*8 :: gdd_threshold
      real*8 :: sand, clay
      real*8 :: smpsat, bch
      real*8 :: watsat, watdry 
      real*8 :: paw
      real*8 :: par
      real*8 :: airtemp
      real*8 :: wat   
      real*8 :: soiltemp     
      real*8 :: soiltemp_10d      
      real*8 :: airtemp_10d
      real*8 :: paw_10d
      real*8 :: par_10d
      real*8 :: gdd
      real*8 :: ncd
      real*8 :: sgdd
      real*8 :: ld
      real*8 :: par_crit
      logical :: par_limit
      real*8 :: turnover0, llspan0
      real*8 :: zweight, zweight30, zweight90
      real*8 :: Soilmoist2layer(N_CASA_LAYERS)
      real*8 :: Soiltemp2layer(N_CASA_LAYERS)
      real*8, dimension(2) :: daylength

      sand = ecp%soil_texture(1)*100.d0
      clay = ecp%soil_texture(3)*100.d0
      smpsat = -10.d0 * ( 10.d0**(1.88d0-0.0131d0*sand) )
      bch = 2.91d0 + 0.159d0*clay
      watsat = 0.489d0 - 0.00126d0*sand 
      watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch) 
  
      airtemp = ecp%TairC
      call Soillayer_convert_Ent(ecp%Soilmoist(:), SOILDEPTH_m, 
     &     Soilmoist2layer) 
      wat = Soilmoist2layer(1)  !Top 30 cm average
      call Soillayer_convert_Ent(ecp%Soiltemp(:), SOILDEPTH_m, 
     &     Soiltemp2layer) 
      soiltemp = Soiltemp2layer(1)

      soiltemp_10d = ecp%soiltemp_10d
      airtemp_10d = ecp%airtemp_10d
      paw_10d = ecp%paw_10d
      par_10d = ecp%par_10d  
      gdd = ecp%gdd
      ncd = ecp%ncd
      sgdd = ecp%sgdd

      zweight=exp(-1.d0/(10.d0*86400.d0/dtsec))  !for 10-day running average
      zweight30=exp(-1.d0/(30.d0*86400.d0/dtsec))  !for 30-day running average
      zweight90=exp(-1.d0/(90.d0*86400.d0/dtsec)) 

      !10-day running average of Soil Temperature
      soiltemp_10d=zweight*soiltemp_10d+(1.0d0-zweight)*soiltemp
      
      !10-day running average of Air Temperature
      airtemp_10d=zweight*airtemp_10d+(1.0d0-zweight)*airtemp
      
      !10-day running average of Plant Available Water  
!SOILMOIST_OLD
!      paw = min( max(wat-watdry,0.d0) 
!     &     /(watsat-watdry), 1.d0)
      paw = wat !now Soilmoist is saturated fraction!!!
      paw_10d=zweight*paw_10d+(1.0d0-zweight)*paw

      !10-day running average of PAR
      par = ecp%IPARdif + ecp%IPARdir
      par_10d=zweight*par_10d+(1.0d0-zweight)*par

      !daylength
      if (ecp%CosZen > 0.d0) then
         ecp%daylength(2) = ecp%daylength(2) + dtsec/60.d0
      end if

      !GDD & NCD - Once a day
      if (dailyupdate) then
         !Growing degree days
         if (airtemp_10d .ge. airtemp_par) then
            gdd = gdd + ( airtemp_10d - airtemp_par )
         end if
         !temp.  GDD using soil temperature
         if (soiltemp_10d .ge. 0.d0) then
            sgdd = sgdd + ( soiltemp_10d - 0.d0 )
         end if
         !Number of chilling days
         if (airtemp_10d .lt. airtemp_par) then
            ncd = ncd +  1.d0
         end if
         !is fall?      
         if (NInt(ecp%daylength(2)) .lt. NInt(ecp%daylength(1)) ) then
            ecp%fall = 1
         else if (NInt(ecp%daylength(2)).gt.NInt(ecp%daylength(1))) then
            ecp%fall = 0
         end if 
         !print*,'fall',ecp%daylength(2),ecp%daylength(1),ecp%fall
       end if

      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 

        cop => pp%tallest
        do while(ASSOCIATED(cop))
          
          !10-day running average of stressH2O (betad)
          cop%betad_10d=zweight*cop%betad_10d
     &                 +(1.d0-zweight)*cop%stressH2O
      
          !daily carbon balance
          cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n*1000.d0

          !*********************************************
          !* evergreen broadleaf - PAR limited    
          !*********************************************     
          par_limit = ((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &                (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
          par_limit = .false. !temp. suppress
          if (par_limit) then
             par_crit = - par_turnover_int/par_turnover_slope 
             turnover0 = min(100.d0, max(0.01d0, 
     &          par_turnover_slope*par_10d + par_turnover_int))
             if (par_10d .lt. par_crit) turnover0 = 0.01d0
             cop%turnover_amp=zweight*cop%turnover_amp 
     &          +(1.d0-zweight)*turnover0
     
             llspan0 = pfpar(cop%pft)%lrage*12.d0/cop%turnover_amp
             cop%llspan=zweight90*cop%llspan
     &          +(1.d0-zweight90)*llspan0     
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
          !write(993,*) cop%Sacclim

          cop => cop%shorter  
        end do        
        pp => pp%younger 
      end do 

      ecp%soiltemp_10d = soiltemp_10d
      ecp%airtemp_10d = airtemp_10d
      ecp%paw_10d = paw_10d
      ecp%par_10d = par_10d
      ecp%gdd = gdd
      ecp%ncd = ncd
      ecp%sgdd = sgdd

      end subroutine clim_stats
      !*********************************************************************   
      subroutine pheno_update(dtsec, pp)
!@sum Update statstics for phneology_update    
!@sum PHYSICAL time step.
      use ent_const

      real*8,intent(in) :: dtsec           !dt in seconds
      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
  
      real*8 :: gdd_threshold
      integer :: pft
      integer :: phenotype
      real*8 :: soiltemp_10d      
      real*8 :: airtemp_10d
      real*8 :: paw_10d
      real*8 :: gdd
      real*8 :: ncd
      real*8 :: sgdd
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      real*8 :: phenofactor
      real*8 :: airt_adj
      real*8 :: soilt_adj
      real*8 :: phenostatus
      logical :: temp_limit, water_limit 
      logical :: fall
      logical :: woody
      integer, parameter :: iwater_limit = 3 !SOILMOIST_OLD 0
      integer, parameter :: itemp_limit = 1
      real*8 :: mature = 1.2d0 !X day to mature: mature=1+X/1000
      real*8::  betad_10d
      real*8 :: ld

      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      paw_10d = pp%cellptr%paw_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      if ( pp%cellptr%fall == 1 ) then
        fall = .true.
      else
        fall = .false.
      endif
      sgdd = pp%cellptr%sgdd
      ld = pp%cellptr%daylength(2)

      gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)  

      cop => pp%tallest
      do while(ASSOCIATED(cop))                
         phenofactor_c=cop%phenofactor_c
         phenofactor_d=cop%phenofactor_d
         phenofactor=cop%phenofactor
         phenostatus=cop%phenostatus
         betad_10d=cop%betad_10d
         pft=cop%pft
         phenotype=pfpar(pft)%phenotype


         if (phenotype .eq. COLDDECID) then 
            temp_limit = .true.
            water_limit = .false.
         else if (phenotype .eq. DROUGHTDECID) then 
            temp_limit = .true.
            water_limit = .true.
         else if (phenotype .eq. EVERGREEN) then
            temp_limit = .false.
            water_limit = .false.
         else !any of cold and drought deciduous
            temp_limit = .true.
            water_limit = .true.            
         end if

         airt_adj=0.d0
         if (pft.eq.DROUGHTDECIDBROAD)airt_adj=10.d0

         soilt_adj=0.d0
!!!fix it         if (pft .eq. GRASSC3ARCTIC)soilt_adj=-5.d0

         woody = pfpar(pft)%woody

         !temperature-controlled woody
         if (temp_limit.and.woody) then
             if ((.not. fall) .and.
     &         (phenostatus.lt.3.d0).and.(gdd.gt.gdd_threshold)) then
               phenofactor_c = min (1.d0,(gdd-gdd_threshold)/gdd_length)
               if (phenofactor_c .lt. 1.d0) then
                  phenostatus = 2.d0
               else
                  phenostatus = 3.d0
               end if
            end if
            if (fall .and. 
     &         (phenostatus.ge.3.d0).and.
     &         (airtemp_10d.lt.airt_max_w+airt_adj)) then
               phenofactor_c = min(phenofactor_c,max(0.d0,
     &            (airtemp_10d-airt_min_w-airt_adj)/
     &            (airt_max_w-airt_min_w)))
               if (phenofactor_c .eq. 0.d0) then
                  phenostatus = 1.d0
                  ncd = 0.d0             !zero-out
                  gdd = 0.d0
               else  
                  phenostatus = 4.d0
               end if 
            end if
            if (fall .and.
     &         (phenostatus.ge.3.d0).and.(ld.lt.ld_max)) then
               phenofactor_c = min(phenofactor_c, max(0.d0,
     &          (ld - ld_min)/(ld_max-ld_min)))
               if (phenofactor_c .eq. 0.0d0) then
                 phenostatus =1.d0
                 ncd =0.d0
                 gdd =0.d0
               else
                 phenostatus=4.d0
               end if
            end if 
         end if

         !temperature-controlled herbaceous 
         if (temp_limit .and. (.not.woody) )then
!            phenofactor_c = 1.d0 !not yet implemented
            select case(itemp_limit)
            case(0)
               phenofactor_c = max(0.d0,min(1.d0,
     &              (soiltemp_10d-soilt_min-soilt_adj)
     &              /(soilt_max-soilt_min)))
               if (phenofactor_c .ge. 0.999d0) then
                  phenostatus = 3.d0
               else if (phenofactor_c .le. EPS) then
                  phenostatus = 1.d0
               end if
            case(1)
               if ((.not. fall) .and.
     &           (phenostatus.lt.3.d0).and.(sgdd.gt.100.d0)) then
                 phenofactor_c = min (1.d0,(sgdd-100.d0)/50.d0)
                 if (phenofactor_c .lt. 1.d0) then
                    phenostatus = 2.d0
                 else
                    phenostatus = 3.d0
                 end if
               end if
               if (fall .and.
     &            (phenostatus.ge.3.d0).and.
     &            (soiltemp_10d.lt.soilt_max+soilt_adj)) then
                  phenofactor_c = min(phenofactor_c,max(0.d0,
     &               (soiltemp_10d-soilt_min-soilt_adj)/
     &               (soilt_max-soilt_min)))
                  if (phenofactor_c .eq. 0.d0) then
                     phenostatus = 1.d0
                     sgdd = 0.d0
                  else  
                     phenostatus = 4.d0
                  end if 
               end if
            end select
         end if            
                 
         !phenofactor_c = 1.d0

         !water-controoled woody
         if (water_limit .and. woody .and. (phenostatus.ge.2.d0)) then
            select case (iwater_limit)
            case(0) !default with the 10-day paw
               phenofactor_d = min(1.d0,max(0.d0,
     &          ((paw_10d-paw_min_w)/(paw_max_w-paw_min_w))**paw_res_w))
            case(1) !water_stress 
               phenofactor_d = water_stress(N_DEPTH  
     i         ,pp%cellptr%Soilmp(:)
     i         ,cop%fracroot(:)
     i         ,pp%cellptr%fice(:), pfpar(pft)%hwilt
     o         , cop%stressH2Ol(:)) 
            case(2) !water_stress3 & betad_10d
               if ((phenostatus .gt. mature) .and.
     &           (phenostatus.lt.3.d0).and.
     &           (betad_10d.gt.betad_min_w))then
                  phenofactor_d = min(1.d0,
     &                 ((betad_10d-betad_min_w)
     &                 /(betad_max_w-betad_min_w))**betad_res_w)
                  if (phenofactor_d .ge. 0.999d0) then
                     phenostatus = 3.d0
                  else
                     phenostatus = 2.d0
                  end if
               else if ((phenostatus.ge.3.d0).and.
     &            (betad_10d.lt.betad_max_w))then
                  phenofactor_d = max(0.d0,
     &                 ((betad_10d-betad_min_w)
     &                 /(betad_max_w-betad_min_w))**betad_res_w)
                  if (phenofactor_d .le. EPS) then
                     phenostatus = 1.d0
                  else
                     phenostatus = 4.d0
                  end if 
               end if    
            case(3) !water_stress3 & betad_10d & no controll in phenostatus
               phenofactor_d = max(0.d0,min(1.d0,
     &              ((betad_10d-betad_min_w)
     &              /(betad_max_w-betad_min_w))**betad_res_w))
            end select
         end if
  
 
         !water-controlled herbaceous
         if (water_limit .and. (.not. woody)) then
            select case(iwater_limit)
            case(0) !default with the 10-day paw 
               if ((phenostatus.lt.3.d0).and.(paw_10d.gt.paw_min_h))then
                  phenofactor_d = min(1.d0,
     &                 ((paw_10d-paw_min_h)
     &                 /(paw_max_h-paw_min_h))**paw_res_h)
                  if (phenofactor_d .lt. 1.d0) then
                     phenostatus = 2.d0
                  else
                     phenostatus = 3.d0
                  end if
               else if ((phenostatus.ge.3.d0).and.
     &                 (paw_10d.lt.paw_max_h))then
                  phenofactor_d = max(0.d0,
     &                 ((paw_10d-paw_min_h)
     &                 /(paw_max_h-paw_min_h))**paw_res_h)
                  if (phenofactor_d .eq. 0.d0) then
                     phenostatus = 1.d0
                  else
                     phenostatus = 4.d0
                  end if 
               end if
  
            case(1) !water_stress 
               phenofactor_d= water_stress(N_DEPTH  
     i         ,pp%cellptr%Soilmp(:)
     i         ,cop%fracroot(:)
     i         ,pp%cellptr%fice(:), pfpar(pft)%hwilt
     o         , cop%stressH2Ol(:))
            case(2) !water_stress3 & betad_10d
               if ((phenostatus .gt. mature) .and.
     &           (phenostatus.lt.3.d0).and.
     &           (betad_10d.gt.betad_min_h))then
                  phenofactor_d = min(1.d0,
     &                 ((betad_10d-betad_min_h)
     &                 /(betad_max_h-betad_min_h))**betad_res_h)
                  if (phenofactor_d .ge. 0.999d0) then
                     phenostatus = 3.d0
                  else
                     phenostatus = 2.d0
                  end if
               else if ((phenostatus.ge.3.d0).and.
     &            (betad_10d.lt.betad_max_h))then
                  phenofactor_d = max(0.d0,
     &                 ((betad_10d-betad_min_h)
     &                 /(betad_max_h-betad_min_h))**betad_res_h)
                  if (phenofactor_d .le. EPS) then
                     phenostatus = 1.d0
                  else
                     phenostatus = 4.d0
                  end if 
               end if    
            case(3) !water_stress3 & betad_10d & no controll in phenostatus
               phenofactor_d = max(0.d0,min(1.d0,
     &              ((betad_10d-betad_min_h)
     &              /(betad_max_h-betad_min_h))**betad_res_h))
               if (phenostatus.le.mature)then
                  phenofactor_d=0.d0
               end if
               if (phenofactor_d .ge. 0.999d0) then
                  phenostatus = 3.d0
               else if (phenofactor_d .le. EPS) then
                  phenostatus = 1.d0
               end if
            end select
          
            if (aint(cop%phenostatus).eq.aint(phenostatus))then
              phenostatus = cop%phenostatus + 1.d0/1000.d0  
	    end if
        
         end if  
         
#ifdef DEBUG
            write(202,'(100e16.6)') phenofactor_d,pp%cellptr%Soilmp(:)
     &      ,pp%cellptr%Soilmoist(:)
     &      ,betad_10d, paw_10d
#endif

         if (.not.temp_limit) phenofactor_c = 1.d0
         if (.not.water_limit) phenofactor_d = 1.d0
      
         phenofactor = phenofactor_c * phenofactor_d
         
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
      subroutine veg_update(dtsec,pp,config)
!@sum Update the vegetation state and carbon pools:
!@sum DAILY call.
!@sum LAI, senescefrac, DBH, height 
!@sum carbon pools of foliage, sapwood, fineroot, hardwood, coarseroot
!@sum AND growth respiration from growth and tissue turnnover
      use ent_const
      use ent_prescr_veg
      implicit none
      type(ent_config) :: config 
      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      real*8 :: dtsec
      integer :: pft
      logical :: woody
      logical :: is_annual 
      logical :: par_limit
      ! phenofactor phenological elongation factor [0,1] (unitless)
      real*8 :: phenofactor
      real*8 :: C_fol_old,C_froot_old,C_sw_old,C_hw_old,C_croot_old
      real*8 :: C_fol, C_froot, C_croot, C_sw, C_hw
      real*8 :: C_lab
      ! Cactive active carbon pool, including foliage, sapwood and fine root (gC/pool/individual)
      real*8 :: Cactive
      ! Cactive_max maximum active carbon pool allowed by the allometric constraint 
      real*8 :: Cactive_max
      ! Cdead dead carbon pool, including hardwood and coarse root (gC/pool/individual)
      real*8 :: Cdead
      ! qsw 
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
 
     !Initialize
      laipatch = 0.d0
      Clossacc(:,:,:) = 0.d0 
      resp_auto_patch = 0.d0
      resp_root_patch = 0.d0
      cohortnum = 0
      resp_growth1=0.d0
      resp_growth2=0.d0
      cpool(:) = 0.d0

      cop => pp%tallest

      do while(ASSOCIATED(cop))
               
         pft = cop%pft
         phenofactor = cop%phenofactor        

         cohortnum = cohortnum + 1


         is_annual = (pfpar(pft)%phenotype .eq. ANNUAL)

        if (is_annual) then
            if (phenofactor .gt. 0.d0 .AND. cop%C_fol .eq. 0.d0) then
               cop%h = 0.05d0 !min. height = 0.05m
	       select case (irecruit)
	       case(1)
               call recruit_annual(cop%pptr%Reproduction(pft),cop%h,
     o              cop%C_fol, cop%C_froot,cop%n, cop%LAI)    
               cop%pptr%Reproduction(pft) = 0.d0
               case(2)
               cop%C_fol = height2Cfol(pft,cop%h)
!               cop%LAI=cop%n*pfpar(pft)%sla*
               cop%LAI=cop%n*sla(pft,cop%llspan)*
     &              (height2Cfol(pft,cop%h)/1000.0d0) 
               cop%C_froot = q*cop%C_fol
               end select 
           end if
         end if
 
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
 
         !----------------------------------------------------
         !*structural/active/reproductive
         !----------------------------------------------------
!          print*,pft,pfpar(pft)%phenotype, COLDDECID,pfpar(pft)%woody
!          print*,cop%phenostatus
         dormant = .false.
         dormant =
     &      (pfpar(pft)%phenotype .eq. COLDDECID .and. 
     &      pfpar(pft)%woody .and. 
     &      cop%phenostatus .lt. 2.d0 .and. cop%phenostatus .ge. 3.d0 )
     
         dCrepro = 0.d0
         if (.not.dormant)           
     &       call growth_cpools_structural(pft,dbh,h,qsw,qf,phenofactor,
     &       C_sw,Cactive_max,C_fol,CB_d,Cactive_old, 
     &       Cactive,C_lab,Cdead,dCrepro) 


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

         call litter_growth_cohort(SDAY,dCrepro,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        dC_litter_hw,dC_litter_croot,cop,Clossacc,resp_growth2)

         if (C_fol_old.eq.0.d0) then
           cop%senescefrac = 0.d0
         else
           cop%senescefrac = l_fract *
     &          (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
         endif

         !* Tissue growth respiration is subtracted at physical time step in canopyspitters.f.
         cop%C_growth = (resp_growth1+resp_growth2)*cop%n*1.d-3  

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


         cop => cop%shorter 
      end do !looping through cohorts
 
      call litter_patch(pp, Clossacc) !Update Tpool from all litter.

      pp%LAI = laipatch  !Update

      !* Update patch fluxes with growth respiration. *!
      !* The daily respiration fluxes are accumulated in C_growth to
      !* be distributed over the course of a day to avoid pulses at night.
!WRONG      pp%R_auto = resp_auto_patch !Total flux including growth increment.
      pp%R_root = pp%R_root + resp_root_patch !Total flux including growth increment.
!      pp%NPP = pp%GPP - resp_auto_patch ##Distribute with C_growth

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
      integer, parameter :: CPotModel = 2
      integer, parameter :: AGrowthModel= 3 !1 grass growth with no storage; 2 grass growth with storage; 3 grass growth with storage after the certain size; 4 grass/tree growth with storage
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
            select case (CPotModel)
            case(1)
               Cactive_pot = Cactive + C_lab
            case(2)
               Cactive_pot = Cactive + min(C_lab, CB_d)!only new carbon is used for growth.
            end select
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
c$$$      subroutine senesce_cpools(is_annual, C_fol_old,C_fol,Cactive,
c$$$     &                          senescefrac, dC_lab)
c$$$      logical, intent(in) :: is_annual  !NOT USED
c$$$      real*8, intent(in) :: C_fol_old
c$$$      real*8, intent(in) :: C_fol
c$$$      !real*8, intent(in) :: C_lab
c$$$      real*8, intent(inout) :: Cactive
c$$$      real*8, intent(out) :: senescefrac
c$$$      real*8, intent(out) :: dC_lab
c$$$      !---- Local ------
c$$$      real*8 :: dC_fol
c$$$      real*8 :: dCactive
c$$$
c$$$      senescefrac = 0.d0
c$$$      if (C_fol_old .gt. C_fol) then 
c$$$         !with senescening, the part of foliage carbon goes to the litter pool
c$$$         ! and the remained goes to the labile.
c$$$         dC_fol = C_fol - C_fol_old !negative
c$$$         dCactive = dC_fol 
c$$$         dC_lab = - dCactive * l_fract
c$$$         if (C_fol_old .ne. 0.d0 ) then 
c$$$            senescefrac = -dC_fol * (1.d0 - l_fract) 
c$$$     &           /C_fol_old
c$$$         endif    
c$$$      else
c$$$         dCactive = 0.d0
c$$$         dC_lab = 0.d0         
c$$$      endif
c$$$
c$$$      !C_lab = C_lab + dC_lab  
c$$$      !Instead, return dC_lab
c$$$      Cactive= Cactive + dCactive
c$$$      
c$$$      end subroutine senesce_cpools



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
          Closs(CARBON,WOOD,i) = cop%n * (turn_hw +
     &         fracrootCASA(i) *turn_croot)
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i) *turn_croot)
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
      end do   !loop through CASA layers-->cumul litter per pool per layer -PK

      !################ ###################################################
      !#### DUE TO TIMING OF LAI UPDATE IN GISS GCM AT THE DAILY TIME STEP,
      !#### GROWTH RESPIRATION FROM CHANGE IN LAI NEEDS TO BE SAVED AS 
      !#### A RESTART VARIABLE IN ORDER TO SEND THAT FLUX TO THE ATMOSPHERE.
      !#### Igor has put in code to distribute C_growth over the day.
      !####################################################################

      cop%C_lab = cop%C_lab + dC_lab

      !* Return Clossacc *!
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo
      !## ? Should Csum be added to C_growth ? ## NK

      !* Return Clossacc

      end subroutine litter_turnover_cohort
      !*********************************************************************
      subroutine litter_growth_cohort(dt,dCrepro,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        dC_litter_hw,dC_litter_croot,cop,Clossacc,resp_growth)
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
      real*8,intent(in) :: dt !seconds, time since last call
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
     &     - dClab_dbiomass - dCrepro        !Growth (new growth or senescence)

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
     &       * fracrootCASA(i) * max(0.d0,-dC_froot)
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
      end do                    !cumul litter per pool per layer 

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
      !at this point, C_lab<0 comes from the rounding errors...
      if (cop%C_lab < 0.d0) cop%C_lab = 0.d0
!      if (cop%C_lab < 0.d0) then
!         print*,dC_fol,cop%C_fol,dC_sw,cop%C_sw
!         print*,dC_lab,cop%C_lab, dC_froot, cop%C_froot
!         print*,dC_hw,cop%C_hw,dC_croot,cop%C_croot
!         stop
!      endif

      !* Return Clossacc *!
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo
      !## ? Should Csum be accumulated in C_growth ? ##NK
      
      !* Return Clossacc

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
      end do  !loop through CASA layers-->cumul litter per pool per layer -PK

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
      cop%C_growth = resp_growth*cop%n*1.d-3  

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
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0
      real*8,parameter :: nplant = 3000.0d0
      height2Cfol=(1.0d0/C2B) *h1Cf 
     &        *((height*100.0d0)**h2Cf)/nplant

      end function height2Cfol
!*************************************************************************
      real*8 function Cfol2height(pft,Cfol)
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !gC/pool/plant
      real*8,parameter :: h1Cf = 5.0d0  !1.660d0
      real*8,parameter :: h2Cf = 1.20d0 !1.500d0
      real*8,parameter :: nplant = 3000.0d0

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
!---------------------------------------------------------------------!
      function water_stress(nlayers, soilmp, fracroot, fice,
     &     hwilt, betadl) Result(betad)
      !1. Rosensweig & Abramopoulos water stress fn.

      implicit none
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  soilmp(:) !Soil matric potential (m)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(in) :: hwilt  !Wilting point of pft, matric pot. (m)
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress
      !---Local-----------
      integer :: k
      
      betad = 0.d0
      do k = 1,nlayers
        betadl(k) = (1.d0-fice(k))*fracroot(k)
!     &       *max((hwilt-soilmp(k))/hwilt,0.d0) !R&A original
     &       *min(1.d0,max((hwilt-soilmp(k))/(hwilt + 25.d0),0.d0))  !With unstressed range to h=-25 m.
        betad = betad + betadl(k) 
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress

!----------------------------------------------------------------------!
    
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
     &        cop%betad_10d, !cop%pptr%cellptr%paw_10d,
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
      !*********************************************************************
c$$$      subroutine veg_init(pp)
c$$$
c$$$      use ent_pfts
c$$$      use entcells
c$$$
c$$$      implicit none
c$$$      save
c$$$
c$$$      type(patch) :: pp
c$$$      type(cohort), pointer :: cop 
c$$$      integer :: pft
c$$$      real*8 :: dbh
c$$$      real*8 :: h
c$$$      real*8 :: qsw
c$$$      real*8 :: ialloc
c$$$
c$$$   
c$$$      pp%cellptr%soiltemp_10d = 0.0d0
c$$$      pp%cellptr%airtemp_10d = 0.0d0
c$$$      pp%cellptr%paw_10d = 0.50d0
c$$$      pp%cellptr%par_10d = 100.d0
c$$$      pp%cellptr%light_prev = pp%cellptr%CosZen
c$$$      !call entcell_print( 6, pp%cellptr )
c$$$
c$$$      cop => pp%tallest   
c$$$
c$$$      do while(ASSOCIATED(cop))
c$$$         pft=cop%pft
c$$$         h=cop%h
c$$$         if (.not.pfpar(pft)%woody) then !grasses/crops/non-woody
c$$$            dbh = 0.0d0  
c$$$            cop%n=cop%LAI/pfpar(pft)%sla/(height2Cfol(pft,h)/1000.0d0) 
c$$$            !write(992,*) "In phenology veg_init herb,cop%n=",cop%n
c$$$            qsw = 0.d0
c$$$            cop%C_hw = 0.0d0
c$$$            cop%C_croot =  0.0d0
c$$$         else
c$$$            dbh=height2dbh(pft,h)
c$$$            cop%n=cop%LAI/pfpar(pft)%sla/(dbh2Cfol(pft,dbh)/1000.0d0)
c$$$            !write(992,*) "In phenology veg_init woody,cop%n=",cop%n
c$$$            qsw = pfpar(pft)%sla*iqsw
c$$$            cop%C_hw = dbh2Cdead(pft,dbh) * hw_fract
c$$$            cop%C_croot = dbh2Cdead(pft,dbh) * (1.0d0 - hw_fract)
c$$$         end if
c$$$         cop%dbh=dbh
c$$$         cop%C_fol = 1000.0d0*cop%LAI/pfpar(pft)%sla/cop%n
c$$$         cop%C_froot =  q*1000.0d0*cop%LAI/pfpar(pft)%sla/cop%n
c$$$         ialloc = (1.0d0+q+h*qsw)
c$$$         cop%C_sw = cop%C_fol * h * qsw  
c$$$         cop%llspan = 36.d0  !late successional tropical
c$$$         cop%phenostatus = 1
c$$$         cop%C_lab = 0.5 * cop%C_fol
c$$$
c$$$         cop => cop%shorter  
c$$$
c$$$      end do
c$$$      end subroutine veg_init
      end module phenology
