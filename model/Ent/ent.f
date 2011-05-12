      module ent
      
!@sum Contains main routines to perform single time step Ent model simulation
!@sum on a single grid cell, entcell, whose parameters are set previously by
!@sum ent_driver routines and the subroutine ent_model called by a 
!@sum main program.

!@auth N. Kiang

!#define ENT_STANDALONE_DIAG TRUE
!#define DEBUG

      !Ent MODULES TO USE
      use ent_const
      use ent_types

      implicit none
      private
      save

      public ent_ecosystem_dynamics,  ent_biophysics
      !public ent_integrate_GISS,
      public ent_integrate !Added by KIM

      contains
      !*********************************************************************
      subroutine ent_ecosystem_dynamics(dtsec,tt,ecp,ALBEDO_FLAG)
!@sum Ent ecosystem dynamics.
      use phenology
      use disturbance
      use canopyrad
      use disturbance
      use reproduction
      use cohorts, only : reorganize_cohorts
      use patches, only : reorganize_patches, summarize_patch
      use entcells, only : summarize_entcell

      real*8,intent(in) :: dtsec
      type(timestruct),pointer :: tt
      type(entcelltype) :: ecp
      logical, optional, intent(in) :: ALBEDO_FLAG
      !---
      type(patch), pointer :: pp
!      write(*,*) 'Ecosystem dynamics for (long,lat)=(',

      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
      !#### THIS LOOP: NEED TO REPLACE ALL CALLS WITH ecp TO CALLS WITH pp ##

        if (ALBEDO_FLAG) then
          call get_patchalbedo(pp)
        end if

!        call ent_integrate(dtsec,ecp,0.0) !Biophysics, respiration

        if (STRUCT_FLAG(tt,ecp)) then
          call reproduction_calc(dtsec, tt, pp)
          call reorganize_cohorts(pp)
!          call phenology_update (dtsec,tt, pp) !UPDATE LAI
          call recalc_radpar (pp) !UPDATE canopy radiative transfer
        end if
        call summarize_patch(pp)
      
        !Flag when it's time to update disturbance.
        !May be at set time intervals, or function of biomass accumulation, etc.
        !For now, monthly update as a place holder.
        if (STRUCT_FLAG_MONTH(tt,ecp)) then
        !* Update phenology and disturbance
        !call phenology_update (dtsec,tt, pp) !UPDATE LAI - put in ent_integrate
          call fire_frequency_cell (dtsec,tt, ecp) !DUMMY
          call recalc_radpar_cell (ecp) !
          call reorganize_patches(ecp)
          call calc_cell_disturbance_rates(dtsec,tt,ecp)
        else
          call calc_cell_disturbance_rates(dtsec,tt,ecp)
        end if
        pp => pp%younger
      end do

      call summarize_entcell(ecp)

        end subroutine ent_ecosystem_dynamics

      !*********************************************************************
!      subroutine ent_integrate_GISS(ecp, dtsec)
!@sum Ent biophysics/biogeochemistry - THIS SUBROUTINE WILL NOT BE NEEDED.
!      use reproduction
!      use cohorts
!      use patches
!      use biophysics, only : photosynth_cond
!      use growthallometry, only : uptake_N
!      use phenology, only : litter
!      use soilbgc, only : soil_bgc
!      use phenology, only : phenology_update
!      use canopyrad, only : recalc_radpar
!      use entcells, only : summarize_entcell, entcell_print
!
!      implicit none
!      type(entcelltype) :: ecp
!      real*8 :: dtsec  !dt in seconds
!!      type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
!      !-----local--------
!      type(patch),pointer :: pp
!
!      pp => ecp%oldest  !changed to => (!) -PK 7/11/06
!      do while (ASSOCIATED(pp)) 
!        call photosynth_cond(dtsec, pp)
!        call uptake_N(dtsec, pp) !Dummy
!        call litter(pp)  !Update litter pools
!        call soil_bgc(dtsec, pp)
!        pp%CO2flux = -pp%NPP + pp%Soil_resp
!        pp%age = pp%age + dtsec
!        call summarize_patch(pp)
!        pp => pp%younger  !changed to => (!) -PK 7/11/06
!      end do
!      call summarize_entcell(ecp)
!
!      end subroutine ent_integrate_GISS

      !*********************************************************************

      subroutine ent_integrate(dtsec, ecp, update_day, config)
!@sum Ent biophysics/biogeochemistry/patch dynamics
      use reproduction
      use cohorts
      use patches
      use biophysics, only : photosynth_cond
      use growthallometry, only : uptake_N
      use soilbgc, only : soil_bgc
!      use phenology, only : litter
      use phenology, only : clim_stats, pheno_update, veg_update
      use canopyrad, only : recalc_radpar
      use entcells, only : summarize_entcell, entcell_print

      implicit none
      real*8 :: dtsec  !dt in seconds
      !type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: ecp
      logical :: update_day
      type(ent_config) :: config 

      !-----local--------
      integer :: patchnum
      type(patch),pointer :: pp


      call clim_stats(dtsec,ecp,config,update_day)

      !* Loop through patches
      patchnum = 0
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        patchnum = patchnum + 1
        !call patch_print(771,pp," ff ")
        call photosynth_cond(dtsec, pp)

        if (config%do_phenology_activegrowth) then
          !call uptake_N(dtsec, pp) !?
          !call growth(...)
          !call patch_print(771,pp," bb ")
          if (update_day) then
            call pheno_update(dtsec,pp)
            !call patch_print(771,pp," bb ")
            call veg_update(dtsec,pp,config)
            !call litter(pp) !Litter is now called within veg_update

          end if
          !call patch_print(771,pp," aa ")
        endif

        if (config%do_soilresp) then
          call soil_bgc(dtsec, pp)
        endif
        pp%CO2flux = -pp%NPP + pp%Soil_resp
        
        pp%age = pp%age + dtsec

          !*********** DIAGNOSTICS FOR PLOTTING ********************!
#ifdef  ENT_STANDALONE_DIAG         
        call summarize_patch(pp)
        call ent_diagnostics(patchnum, pp)
#endif
          !*********************************************************!
        pp => pp%younger 
      end do

      call summarize_entcell(ecp)

      if (config%do_patchdynamics) then
      !  call patch_dynamics(pp,monthlyupdate)
      ! call summarize_entcell(ecp)
      endif !do_patchdynamics
      
      if (update_day) then
        ecp%daylength(1) = ecp%daylength(2)
        ecp%daylength(2) = 0.d0
      end if

#ifdef DEBUG
      print *,"End of ent_biophysics"
      call entcell_print(6, ecp)
#endif

#ifdef ENT_STANDALONE_DIAG
      call ent_diagnostics_entcell(ecp)
!      call entcell_print(6, ecp)
#endif

      end subroutine ent_integrate


      !*********************************************************************
      subroutine ent_biophysics(dtsec, ecp, config)
!@sum  Photosynthesis CO2 uptake.
!@+    If do_soilresp, then also  soil respiration for net CO2 fluxes.
      use biophysics, only : photosynth_cond
      use soilbgc, only : soil_bgc
      use patches, only : summarize_patch
      use entcells, only : summarize_entcell, entcell_print
!      use phenology,only : litter !### Igor won't like this here.
      implicit none
      real*8 :: dtsec  !dt in seconds
      type(entcelltype) :: ecp
      type(ent_config) :: config
      !---Local--------
      type(patch),pointer :: pp
      integer :: patchnum
       
      patchnum = 0
      pp => ecp%oldest
      do while(ASSOCIATED(pp))
        patchnum = patchnum + 1
        !print*,'NEXT PATCH'
        !print*,'Calling photosynth_cond'
        call photosynth_cond(dtsec, pp)
        if (config%do_soilresp) then 
          !print*,'Calling soil_bgc'
          call soil_bgc(dtsec,pp)
          pp%CO2flux = -pp%NPP + pp%Soil_resp
        ! Litter is updated daily in ent_prescribe_vegupdates.
        ! is do_soilresp flag ok or different flag is needed ?
        ! if ( dailyupdate ) call litter(pp) 

          !*********** DIAGNOSTICS FOR PLOTTING ********************!
#ifdef ENT_STANDALONE_DIAG         
          call summarize_patch(pp)
          call ent_diagnostics(patchnum,pp)
#endif
          !*********************************************************!
        else 
          pp%CO2flux = UNDEF
        endif
        pp%age = pp%age + dtsec
        !call summarize_patch(pp)

        pp => pp%younger

      end do
      call summarize_entcell(ecp)

#ifdef DEBUG
      print *,"End of ent_biophysics"
      call entcell_print(6, ecp)
      print *,"*"
#endif
      end subroutine ent_biophysics
      !*********************************************************************

      function STRUCT_FLAG(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !------local------
        update_struct = STRUCT_FLAG_DAY(tt,ecp)

      end function STRUCT_FLAG
      !*********************************************************************

      function STRUCT_FLAG_DAY(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple end-of-day flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !-----local----------
        real*8 :: hourfrac
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
!        if (hourfrac.le.dtsec) then !Midnight
        if (hourfrac.eq.0.0) then  !Midnight
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_DAY

      !*********************************************************************

      function STRUCT_FLAG_MONTH(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple beginning-of-the-month flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        
        real*8 :: hourfrac
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
        if ((tt%day.eq.1).and.
     &       (hourfrac.eq.0.0)) then
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_MONTH


!*****************************************************************************
      subroutine ent_diagnostics(patchnum, pp)
      !*********** DIAGNOSTICS FOR PLOTTING ********************!
      !use patches, only : patch_print
      implicit none
      integer :: patchnum
      type(patch),pointer :: pp
      !---Local------
      integer :: tmp_pft
      real*8 :: tmp_n,tmp_senescefrac,tmp_Sacclim
      type(cohort),pointer :: cop
!#ifdef DEBUG

      tmp_pft = -1
      tmp_n = -1
      tmp_senescefrac = 0.d0
      tmp_Sacclim = -1

      !call patch_print(6,pp)
      if ( ASSOCIATED(pp%tallest) ) then
        tmp_pft = pp%tallest%pft
        tmp_n = pp%tallest%n
        tmp_senescefrac = pp%tallest%senescefrac
        tmp_Sacclim = pp%tallest%Sacclim
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          cop => cop%shorter
        end do
      endif

      write(995,'(i5,3(1pe16.8),i5,100(1pe16.8))') !Fluxes are positive up.
     &     patchnum,pp%cellptr%IPARdir,pp%cellptr%IPARdif, 
     &     pp%cellptr%coszen,
     &     tmp_pft,tmp_n,pp%lai, pp%h, pp%Tpool(CARBON,:,:), 
     &     pp%C_fol, pp%C_w, pp%C_froot, pp%C_root, pp%C_lab,
     &     pp%Reproduction(tmp_pft),
     &     pp%TRANS_SW,
     &     pp%Ci, pp%GPP,pp%R_auto,pp%Soil_resp,
     &     pp%NPP,pp%CO2flux,pp%GCANOPY,pp%IPP,
     &     tmp_senescefrac,tmp_Sacclim,pp%c_total,
!### HACK: c_growth is Igor's hack to store daily growth respiration somewhere
!### HACK: N_up is temporarily litterfall, using unused variable -NK
     &     pp%c_growth,pp%N_up,pp%betad

      if (pp%GPP.lt.0.d0) then
        print *,"ent.f: BAD GPP:",pp%lai, pp%GPP
      endif
!      write(999,*) pp%cellptr%Soilmp, pp%Soilmoist
!     &     ,pp%cellptr%betad, pp%cellptr%betadl
!      write(994,*) pp%cellptr%GCANOPY
!#endif
      end subroutine ent_diagnostics
          !*********************************************************!

      subroutine ent_diagnostics_entcell(ecp)
      implicit none
      type(entcelltype) :: ecp

      write(996,*) ecp%area,ecp%IPARdir,ecp%IPARdif,ecp%LAI,ecp%fv
     &     ,ecp%Tpool(CARBON,:,:)
     &     ,ecp%C_fol, ecp%C_w, ecp%C_froot, ecp%C_root, ecp%C_lab
     &     ,ecp%TRANS_SW
     &     ,ecp%Ci, ecp%GPP,ecp%R_auto,ecp%Soil_resp
     &     ,ecp%NPP,ecp%CO2flux,ecp%GCANOPY

      end subroutine ent_diagnostics_entcell

      end module ent
