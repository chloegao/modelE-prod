#include "rundeck_opts.h"

      module irrigate_crop

!@sum Contains routines to read in irrigation data and to compute the
!@sum potential irrigation rate
!@auth M.J. Puma, R. Ruedy

      implicit none
      private

      public irrigate_flux
      public init_irrigate
      public irrigate_extract

      SAVE

!@var Monthly irrigation arrays read in to compute daily irrigation
      REAL*8,ALLOCATABLE,DIMENSION(:,:) :: irr_month_0,irr_month_1
!@var iu_irrigate is unit number for irrigation climatologies
      integer :: iu_irrigate
!@dbparam irrig_cycl determines whether prescribed irrigation data
!@+       vary year to year:  0=no; 1=yes
      integer :: irrig_cycl
!@var Potential irrigation rate [m/s] 
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: irrig_water_pot

      contains

!-----------------------------------------------------------------------

      subroutine init_irrigate()
      use filemanager, only : openunit
      use Dictionary_mod, only : sync_param
      USE DOMAIN_DECOMP_ATM, ONLY : grid, getDomainBounds
      implicit none

      integer :: I_0,I_1,J_0,J_1
      INTEGER :: I_0H, I_1H, J_1H, J_0H
      INTEGER :: IER, i,j

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE( irr_month_0   ( I_0H:I_1H , J_0H:J_1H ),
     &          irr_month_1   ( I_0H:I_1H , J_0H:J_1H ),
     &          STAT=IER)
      ALLOCATE( irrig_water_pot ( I_0H:I_1H , J_0H:J_1H ), STAT=IER)

      call openunit("IRRIG",iu_irrigate,.true.,.true.)
      call sync_param("irrig_cycl",irrig_cycl)
      call getDomainBounds(grid,I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

!**** Compute irrigation rate
      call irrigate_flux(.false.)

      end subroutine init_irrigate

!-----------------------------------------------------------------------

      subroutine irrigate_flux(endofday)
!@sum  Calculates daily irrigation from monthly irrigation data.  Routine
!@sum  sets imon, irr_month_0, and irr_month_1 depending on jday and
!@sum  jyear (cyclical case does not need jyear).
      USE CONSTANT, only : rhow,teeny,shw,sday
      use model_com, only : Itime,JDmidOfM,itimei,JMperY
      USE DOMAIN_DECOMP_ATM, ONLY : grid, getDomainBounds, am_i_root
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR,
     &                             INT_MONTHS_PER_YEAR
     &                          ,READT_PARALLEL,REWIND_PARALLEL
     &                          ,READ_PARALLEL, MREAD_PARALLEL
     &                          ,BACKSPACE_PARALLEL
      USE FILEMANAGER, only : NAMEUNIT
      USE DIAG_COM, only : aij=>aij_loc, ij_irrW_tot
      implicit none

!@var endofday flag to determine whether this is an initial or daily call
      logical, INTENT (IN) :: endofday ! true if called from daily
!@var m_yr_mnth: title of the monthly irrigation data (YYYYMM)
      integer :: m_yr_mnth
!@var JDLAST julian day that irrigate_flux was last called
      integer, SAVE :: jdlast=0
!@var IMON counter to identify data for later of two months
!@var needed for interpolation (irr_month_1) ->
!@var IMON-1 is for irr_month_0
      integer, SAVE :: imon = 0
      integer :: i,j,m,lstmon
      integer :: I_0,I_1,J_0,J_1
      real*8 :: FRAC

      call getDomainBounds(grid,I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)

      if (irrig_cycl == 0) then
!**** Non-cyclical case: data vary from year to year
         if(jdlast == 0)then ! first call
!           Position within the irrigation file based on jyear and jmon
            lstmon = jyear*100 + jmon
            imon = jmon + 1
            if (jday < JDmidOFM(imon))then
               imon = jmon
               lstmon=lstmon-1
               if(jmon==1) lstmon = lstmon-100+12
            endif

            call READ_PARALLEL(m_yr_mnth, iu_irrigate)

            do while (m_yr_mnth < lstmon)
               call READ_PARALLEL(m_yr_mnth, iu_irrigate)
            enddo

            if(m_yr_mnth /= lstmon ) then
              call stop_model(
     &        'Read error: irrig. data not available for start date'
     &         ,255)
            endif

!           Read in irr_month_0
            call BACKSPACE_PARALLEL( iu_irrigate )
            call MREAD_PARALLEL(grid,iu_irrigate,
     &                 NAMEUNIT(iu_irrigate),m,irr_month_0)

            if(AM_I_ROOT())then
               write(6,*) 'Read in irrigate for month,',jmon,m
            endif

         else

           if (jday /= jdlast+1) then !have a new year: reset imon 13->1
             imon=imon-INT_MONTHS_PER_YEAR
             go to 330
           end if

           if (jday <= JDmidOFM(imon)) go to 330

           imon=imon+1          ! read in new month of climatological data
           irr_month_0 = irr_month_1

         endif

!        Read in irr_month_1
         CALL MREAD_PARALLEL(grid,iu_irrigate,
     &          NAMEUNIT(iu_irrigate),m,irr_month_1)
         if(AM_I_ROOT())then
            WRITE(6,*) 'Read in irrigate for next month,',jmon,m
         endif

 330     jdlast=jday

!**** End of non-cyclical case

      else

!**** Cyclical case: data does not from year to year
         if (JDLAST == 0) then ! need to read in first month climatology
           IMON=1          ! IMON=January
           if (JDAY < 16)  then ! JDAY in Jan 1-15, first month is Dec
             print *, "iu_irrigate, IPOS =", iu_irrigate, 12
             CALL READT_PARALLEL(grid,iu_irrigate,NAMEUNIT(iu_irrigate),
     &                        irr_month_0,12)
             CALL REWIND_PARALLEL( iu_irrigate )
           else            ! JDAY is in Jan 16 to Dec 16, get first month
  620        IMON=IMON+1
             IF(JDAY > JDmidOFM(IMON) .and. IMON <= INT_MONTHS_PER_YEAR)
     &          GO TO 620
             CALL READT_PARALLEL
     &       (grid,iu_irrigate,NAMEUNIT(iu_irrigate),irr_month_0,IMON-1)
             if (IMON == 13)  CALL REWIND_PARALLEL( iu_irrigate )
           endif
         ELSE                      ! Do we need to read in second month?
           IF (JDAY /= JDLAST+1) THEN ! Check that data is read in daily
             IF (JDAY /= 1 .or. JDLAST /= INT_DAYS_PER_YEAR) THEN
             WRITE(6,*) 'Incorrect values in irrigate,JDAY,JDLAST=',JDAY
     &           ,JDLAST
               call stop_model(
     &           'ERROR READING IN SETTING IRRIGATION CLIMATOLOGY',255)
             END IF
             IMON=IMON-INT_MONTHS_PER_YEAR        ! New year
             GO TO 630
           END IF
           IF (JDAY <= JDmidOFM(IMON)) GO TO 630
           IMON=IMON+1          ! read in new month of climatological data
           irr_month_0 = irr_month_1
           IF (IMON == 13) CALL REWIND_PARALLEL( iu_irrigate )
         END IF
         CALL READT_PARALLEL(grid,iu_irrigate,NAMEUNIT(iu_irrigate),
     &                    irr_month_1,1)
 630     JDLAST=JDAY
      endif

!**** Interpolate daily irrigation depth to the current day
      FRAC = REAL(JDmidOFM(IMON)-JDAY,KIND=8)/
     &           (JDmidOFM(IMON)-JDmidOFM(IMON-1))
      if(am_I_root())
     &write(6,*) 'interp. irrigate to it,day,fr0',itime,jday,frac
!**** Compute irrigation rate
      do j=J_0,J_1
         do i=I_0,I_1
!****      Input irrigation values in m/s
           irrig_water_pot(i,j)= ( FRAC*irr_month_0(i,j)
     &                          + (1.d0-FRAC)*irr_month_1(i,j) )
!****      Check to make sure no negative irrigation rates
           if( irrig_water_pot(i,j) < 0.d0 ) then
              irrig_water_pot(i,j) = 0.d0
           endif
C**** diagnostic
           if (endofday) aij(i,j,ij_irrW_tot)=aij(i,j,ij_irrW_tot)+
     &          irrig_water_pot(i,j)*SECONDS_PER_DAY*1000.d0

         enddo
      enddo

      end subroutine irrigate_flux

!-----------------------------------------------------------------------

      subroutine irrigate_extract(i,j,MWL,GML,MLDLK,tlake,flake
     *     ,hlake_min,MWL_to_irrig,GML_to_irrig,irrig_gw,irrig_gw_energy
     *     ,irrig_water_act,irrig_energy_act
#ifdef TRACERS_WATER
     *     ,TRML,TRML_to_irrig,irrig_tracer_act,irrig_gw_tracer
#endif
     *       )
!@sum irrigate_extract calculate water used for irrigation
!@auth Michael Puma
      USE CONSTANT, only : rhow,teeny,shw
      USE MODEL_COM, only : dtsrc
!      USE sle001, only : tp ! tp is not saved for each gridpoint
      USE GHY_COM, only : tearth
#ifdef TRACERS_WATER
      USE TRACER_COM, only : ntm
      USE FLUXES, only : atmlnd
#endif
! fixed i,j arrays - feed in from call?
      USE GEOM, only : axyp
      USE GHY_COM, only : fearth

      implicit none
      integer,intent(in):: i,j

C**** Inputs
!@var MWL,GML,MLDLK,tlake,flake local versions of lake variables
      REAL*8, INTENT(IN):: MWL,GML,MLDLK,tlake,flake,hlake_min
#ifdef TRACERS_WATER
     *     ,TRML(NTM,2)
#endif
C**** Outputs
!@var MWL_to_irrig (kg), GML_to_irrig (J) mass/energy changes
      REAL*8, INTENT(OUT) :: MWL_to_irrig,GML_to_irrig
#ifdef TRACERS_WATER
!@var TRML_to_irrig (kg) tracer changes
     *     ,TRML_to_irrig(NTM,2)
#endif
!@var irrig_gw (m/s), irrig_gw_energy (W/m2) ground water diagnostics 
!@var irrig_water_act (m/s), irrig_energy_act (W/m2) actual irrigation diagnostics 
      REAL*8, INTENT(OUT) :: irrig_gw,irrig_gw_energy,irrig_water_act
     *     ,irrig_energy_act
#ifdef TRACERS_WATER
!@var irrig_tracer_act (kg/s) actual irrigation tracer diagnostics 
!@var irrig_gw_tracer_act (kg/s) implied ground water tracer diagnostics 
     *     ,irrig_tracer_act(NTM),irrig_gw_tracer(NTM)
#endif

!@var Mass of irrigation water withdrawn from rivers/groundwater [kg]
      real*8 :: m_irr_pot
!@var Mass of water available for irrigation [kg]
      real*8 :: m_avail
!@var Temperature of abstracted irrigation water [deg C]
      real*8 :: T_irr, T_irr_g, T_irr2
!@var Vegetated and bare fractions of the grid cell
      real*8 :: fv, fb
!@param Conversion factor from m^3 water/(m^2 grid cell) to kg water
      real*8 :: m_to_kg
!@param Flag for external irrigation source (from deep aquifers - not modeled) 
      integer,parameter :: flag_irrig_grndwat = 1  ! =1 use ground water

C**** set default output
      irrig_water_act = 0.d0 ; irrig_energy_act = 0.d0
      irrig_gw = 0.d0 ; irrig_gw_energy = 0.d0
      MWL_to_irrig = 0 ; GML_to_irrig = 0
#ifdef TRACERS_WATER
      irrig_tracer_act = 0 
      irrig_gw_tracer = 0 
      TRML_to_irrig = 0
#endif

      m_to_kg = rhow*axyp(i,j)

      call get_fb_fv( fb, fv, i, j )

!***  Irrigation based on potential values in m/s per area grid cell 
      if ( (irrig_water_pot(i,j) > teeny)  .and. 
     &     (fearth(i,j)          > teeny)  .and.
     &     (fv                   > teeny)) then

         m_irr_pot = irrig_water_pot(i,j)* m_to_kg *dtsrc

         if (flake > 0.d0) then
            m_avail = mwl - hlake_min*flake*m_to_kg
            m_avail = max(m_avail, 0.d0)
            T_irr = tlake
            T_irr2 =0.
            if (mwl.gt.flake*mldlk*m_to_kg+teeny) T_irr2 = 
     *           (gml-mldlk*m_to_kg*flake*tlake*shw)/
     *           (mwl-mldlk*m_to_kg*flake+teeny)/shw 
         else
            m_avail = mwl
            T_irr = gml/(mwl*shw+teeny)
            T_irr2 = T_irr
         endif
!        Check these limits !!!!
         T_irr = max(T_irr, 0.d0)
         T_irr2 = max(T_irr2, 0.d0)
! need to reconstuct local tp(1,2) using ground hydrology code
!         T_irr_g = max(tp(1,2),0.d0)
         T_irr_g = max(tearth(i,j),0.d0)

!***     Set actual irrigation rates and update mwl and gml (if necessary)
         if (m_avail <= teeny) then

            if(flag_irrig_grndwat .ne. 0) then
               irrig_water_act = irrig_water_pot(i,j)
               irrig_energy_act= irrig_water_act*shw*
     &                                 T_irr_g*rhow
               irrig_gw        = irrig_water_act
               irrig_gw_energy = irrig_energy_act
#ifdef TRACERS_WATER
               irrig_tracer_act = irrig_water_act*atmlnd%gtracer(:,i,j)
               irrig_gw_tracer  = irrig_gw*atmlnd%gtracer(:,i,j)
#endif
            endif

         elseif (m_avail >= m_irr_pot) then

            mwl_to_irrig = irrig_water_pot(i,j)*m_to_kg*dtsrc
            if (flake.gt.0 .and. mwl_to_irrig .gt. mldlk*m_to_kg*flake) 
     *           then           ! need layer 2 water
              gml_to_irrig = mldlk*m_to_kg*flake*shw*T_irr + 
     *             (mwl_to_irrig-mldlk*m_to_kg*flake)*shw*T_irr2
#ifdef TRACERS_WATER
              trml_to_irrig(:,1)=trml(:,1)
              trml_to_irrig(:,2)=trml(:,2)/(mwl-mldlk*m_to_kg*flake)
     *             *(mwl_to_irrig-mldlk*m_to_kg*flake)
#endif
            else
              gml_to_irrig = mwl_to_irrig*shw*T_irr
#ifdef TRACERS_WATER
              if (flake.gt.0) then
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/
     *                (mldlk*m_to_kg*flake)
              else
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/mwl
              endif                 
              trml_to_irrig(:,2)=0.
#endif
            end if
            irrig_water_act = irrig_water_pot(i,j)
            irrig_energy_act= gml_to_irrig*rhow / (m_to_kg*dtsrc) 
#ifdef TRACERS_WATER
            irrig_tracer_act = (trml_to_irrig(:,1)+trml_to_irrig(:,2))/
     *           (m_to_kg*dtsrc)
#endif
            
         else !!! (m_avail < m_irr_pot)

            mwl_to_irrig = m_avail
            if (flake.gt.0 .and. mwl_to_irrig .gt. mldlk*m_to_kg*flake) 
     *           then           ! need layer 2 water 
              gml_to_irrig = mldlk*m_to_kg*shw*T_irr*flake + 
     *             (mwl_to_irrig-mldlk*m_to_kg*flake)*shw*T_irr2
#ifdef TRACERS_WATER
              trml_to_irrig(:,1)=trml(:,1)
              trml_to_irrig(:,2)=trml(:,2)/(mwl-mldlk*m_to_kg*flake)
     *             *(mwl_to_irrig-mldlk*m_to_kg*flake)
#endif
            else
              gml_to_irrig = m_avail*shw*T_irr
#ifdef TRACERS_WATER
              if (flake.gt.0) then
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/
     *                (mldlk*m_to_kg*flake)
              else
                 trml_to_irrig(:,1)=trml(:,1)*mwl_to_irrig/mwl
              end if
              trml_to_irrig(:,2)=0.
#endif
            end if

            if(flag_irrig_grndwat == 0) then
               irrig_water_act  = m_avail / (m_to_kg*dtsrc) 
               irrig_energy_act = gml_to_irrig*rhow / (m_to_kg*dtsrc) 
#ifdef TRACERS_WATER
               irrig_tracer_act=(trml_to_irrig(:,1)+trml_to_irrig(:,2))/
     *              (m_to_kg*dtsrc)
#endif
            else
               irrig_water_act  = irrig_water_pot(i,j)
               irrig_energy_act = irrig_water_act*shw*T_irr*rhow
               irrig_gw = irrig_water_act - m_avail / (m_to_kg*dtsrc)
               irrig_gw_energy  = irrig_energy_act - gml_to_irrig*rhow /
     *              (m_to_kg*dtsrc)  
#ifdef TRACERS_WATER
               irrig_tracer_act = irrig_water_act*(trml_to_irrig(:,1)
     *              +trml_to_irrig(:,2)) / mwl_to_irrig
               irrig_gw_tracer = irrig_tracer_act - (trml_to_irrig(:,1)
     *              +trml_to_irrig(:,2)) / (m_to_kg*dtsrc)
#endif
            endif

         endif

      endif

!!!!! HACK HACK test test
      if( (irrig_water_act-irrig_water_pot(i,j)) > 1.d-30 ) then
         write(6,*) 'error in irrigation: act> pot'
         write(6,*) 'i,j,irrig_water_act, irrig_water_pot(i,j):',
     &              i,j,irrig_water_act, irrig_water_pot(i,j)
         write(6,*) 'm_avail, m_irr_pot:', m_avail, m_irr_pot
         irrig_water_act=irrig_water_pot(i,j)
      endif
!!!!!
      end subroutine irrigate_extract


!-----------------------------------------------------------------------

      end module irrigate_crop
