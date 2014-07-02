#include "rundeck_opts.h" 
      MODULE AEROSOL_SOURCES
!@sum repository for Koch aerosol sources, features, etc.
!@auth Dorothy Koch
!@ subroutines in this file include:
!@ alloc_aerosol_sources
!@ get_O3_offline
!@ read_mon3Dsources
!@ READ_OFFHNO3
!@ READ_OFFSS
!@ read_DMS_sources
!@ read_seasalt_sources
!@ aerosol_gas_chem
!@ SCALERAD
!@ GET_SULFATE
!@ GET_BC_DALBEDO
!@ GRAINS
!@ read_mon_3D
!@ read_seawifs_chla
      use OldTracer_mod, only: om2oc
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: ndmssrc  = 1
!@var DMSinput           DMS ocean source (kg/s/m2)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSinput ! DMSinput(im,jm,12)
c!@var DMS_AER           DMS prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: DMS_AER  !(im,jm,366)
c!@var SS1_AER        SALT bin 1 prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: SS1_AER  !(im,jm,366)
c!@var SS2_AER        SALT bin 2 prescribed by AERONET (kg S/day/box)
      real*4, ALLOCATABLE, DIMENSION(:,:,:) :: SS2_AER  !(im,jm,366)
#ifdef TRACERS_AEROSOLS_OCEAN
!@var OC_SS_enrich_fact OCocean enrichment factor of seasalt1
      real*8, ALLOCATABLE, DIMENSION(:,:) :: OC_SS_enrich_fact !(im,jm)
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifndef TRACERS_AEROSOLS_SOA
!@var OCT_src    OC Terpene source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: OCT_src !(im,jm,12)
#endif  /* TRACERS_AEROSOLS_SOA */
!@var ss_src  Seasalt sources in 2 bins (kg/s/m2)
      INTEGER, PARAMETER :: nsssrc = 2
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: ss_src !(im,jm,nsssrc)
      INTEGER, PARAMETER :: nso2src_3d  = 1
!@var SO2_src_3D SO2 volcanic sources (and biomass) (kg/s)
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: SO2_src_3D !(im,jm,lm,nso2src_3d)
!@var PBLH boundary layer height
!@var MDF is the mass of the downdraft flux
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: 
     *   oh,dho2,perj,tno3,o3_offline  !im,jm,lm
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: ohr,dho2r,perjr,
     *   tno3r,ohsr  !im,jm,lm,12   DMK jmon
      integer :: JmonthCache = -1
      real*8, allocatable, dimension(:,:,:) :: 
     &     ohrCache, dho2rCache, perjrCache, tno3rCache
      
      real*8, ALLOCATABLE, DIMENSION(:,:) :: snosiz
#ifdef TRACERS_RADON
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: rn_src
#endif
!var off_HNO3 off-line HNO3 field, used for nitrate and AMP when gas phase chemistry turned off
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)     ::  off_HNO3, off_SS
!@dbparam tune_ss1, tune_ss2 factors to tune seasalt sources
      real*8 :: tune_ss1=1.d0, tune_ss2=1.d0
!@var BBinc enhancement factor of BB carbonaceous aerosol emissions (Kostas: should this be applied to all BB emitted tracers?)
      real*8:: BBinc=1.0d0
#ifdef TRACERS_AEROSOLS_VBS
!@var VBSemifact factor that distributes organic aerosols in volatility bins
      real*8, allocatable, dimension(:) :: VBSemifact
#endif /* TRACERS_AEROSOLS_VBS */

      END MODULE AEROSOL_SOURCES

      SUBROUTINE alloc_aerosol_sources(grid)
!@auth D. Koch
      use domain_decomp_atm, only: dist_grid, getDomainBounds
      use TRACER_COM, only: NTM, n_OCII
      use AEROSOL_SOURCES, only: DMSinput,DMS_AER,SS1_AER,SS2_AER,
#ifndef TRACERS_AEROSOLS_SOA
     * OCT_src,
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
     * OC_SS_enrich_fact,
#endif  /* TRACERS_AEROSOLS_OCEAN */
     * nsssrc,ss_src,nso2src_3d,SO2_src_3D,
     * ohr,dho2r,perjr, tno3r, 
     * ohrCache, dho2rCache, perjrCache, tno3rCache,
     * oh,dho2,perj,tno3,ohsr
     * ,o3_offline
     * ,snosiz
     * ,off_HNO3,off_SS
#ifdef TRACERS_RADON
     * ,rn_src
#endif
#ifdef TRACERS_AEROSOLS_VBS
     * ,VBSemifact
      use TRACERS_VBS, only: vbs_tr
#endif

      use RESOLUTION, only: im,lm
      
      IMPLICIT NONE
      type (dist_grid), intent(in) :: grid
      integer ::  J_1H, J_0H, I_0H, I_1H
      integer :: IER
      logical :: init = .false.

      if(init)return
      init=.true.

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      allocate( DMSinput(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
      allocate( DMS_AER(I_0H:I_1H,J_0H:J_1H,366) ,STAT=IER)
      allocate( SS1_AER(I_0H:I_1H,J_0H:J_1H,366) ,STAT=IER)
      allocate( SS2_AER(I_0H:I_1H,J_0H:J_1H,366) ,STAT=IER) 
#ifdef TRACERS_AEROSOLS_OCEAN
      allocate( OC_SS_enrich_fact(I_0H:I_1H,J_0H:J_1H) ,STAT=IER) 
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifndef TRACERS_AEROSOLS_SOA
      allocate( OCT_src(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
#endif  /* TRACERS_AEROSOLS_SOA */
      allocate( ss_src(I_0H:I_1H,J_0H:J_1H,nsssrc) ,STAT=IER)
      allocate( SO2_src_3D(I_0H:I_1H,J_0H:J_1H,lm,nso2src_3d),STAT=IER )
      allocate( oh(I_0H:I_1H,J_0H:J_1H,lm),dho2(I_0H:I_1H,J_0H:J_1H,lm),
     * perj(I_0H:I_1H,J_0H:J_1H,lm),tno3(I_0H:I_1H,J_0H:J_1H,lm)
     * ,o3_offline(I_0H:I_1H,J_0H:J_1H,lm),STAT=IER )
      allocate( ohr(I_0H:I_1H,J_0H:J_1H,lm),
     * dho2r(I_0H:I_1H,J_0H:J_1H,lm),
     * perjr(I_0H:I_1H,J_0H:J_1H,lm),tno3r(I_0H:I_1H,J_0H:J_1H,lm),
     * ohsr(I_0H:I_1H,J_0H:J_1H,lm),STAT=IER )
      allocate( ohrCache(I_0H:I_1H,J_0H:J_1H,lm),
     * dho2rCache(I_0H:I_1H,J_0H:J_1H,lm),
     * perjrCache(I_0H:I_1H,J_0H:J_1H,lm),
     *     tno3rCache(I_0H:I_1H,J_0H:J_1H,lm))
      allocate( snosiz(I_0H:I_1H,J_0H:J_1H) ,STAT=IER)
#ifdef TRACERS_RADON
      allocate( rn_src(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
#endif
c off line 
      allocate(  off_HNO3(I_0H:I_1H,J_0H:J_1H,LM)     )
      allocate(  off_SS(I_0H:I_1H,J_0H:J_1H,LM)     )
#ifdef TRACERS_AEROSOLS_VBS
      allocate(VBSemifact(vbs_tr%nbins))
#endif

      return
      end SUBROUTINE alloc_aerosol_sources      
      SUBROUTINE get_O3_offline
!@sum read in ozone fields for aqueous oxidation
c
C**** GLOBAL parameters and variables:
C
      use resolution, only: im,jm,lm
      use model_com, only: modelEclock
      use filemanager, only: openunit,closeunit
      use aerosol_sources, only: o3_offline
      use domain_decomp_atm, only: grid, getDomainBounds, write_parallel
C     
      implicit none

C**** Local parameters and variables and arguments:
c
!@var nmons: number of monthly input files
      integer, parameter :: nmons=1,levo3=23
      integer, dimension(nmons) :: mon_units,imon
      integer :: i,j,k,l
      integer :: jdlast=0
      logical :: ifirst=.true.
      character*80 title
      character(len=300) :: out_line
      character*10 :: mon_files(nmons) = (/'O3_FIELD'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 :: frac
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,levo3,1):: src
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb
      save jdlast,mon_units,imon,ifirst,tlca,tlcb
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0H, I_1H

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      if (ifirst) then
        call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        I_0H = grid%I_STRT_HALO
        I_1H = grid%I_STOP_HALO

        allocate(tlca(i_0H:i_1H,j_0H:j_1H,levo3,nmons),
     &           tlcb(i_0H:i_1H,j_0H:j_1H,levo3,nmons))
        ifirst=.false.
      endif
      k=1
      call openunit(mon_files(k),mon_units(k),mon_bins(k),.true.)
      call read_mon3Dsources(levo3,mon_units(k),jdlast,
     & tlca(:,:,:,k),tlcb(:,:,:,k),src(:,:,:,k),frac,imon(k))     
      call closeunit(mon_units(k))
      
      jdlast = modelEclock%dayOfYear()

      if(levo3 /= LM) then ! might be ok, but you should check
        write(out_line,*)
     &  'make sure levels are correct in get_O3_offline'
        call write_parallel(trim(out_line))
        call stop_model('check on get_O3_offline',255)
      endif
      
      o3_offline(:,J_0:J_1,:)=src(:,J_0:J_1,:,k)

      write(out_line,*)
     &'offline ozone interpolated to current day',frac
      call write_parallel(trim(out_line))
      
      return
      END SUBROUTINE get_O3_offline
      
      
      SUBROUTINE read_mon3Dsources(Ldim,iu,jdlast,tlca,tlcb,data1,
     & frac,imon)
! we need to combine this with the one that is used in TRACERS_
! SPECIAL_Shindell, (called read_monthly_3Dsources) and then move it
! into more general tracer code...
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others / Greg Faluvegi
      use resolution, only: im,jm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only: idofm=>JDmidOfM
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR,INT_MONTHS_PER_YEAR
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL,
     &     REWIND_PARALLEL,write_parallel
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,jdlast
      character(len=300) :: out_line
      real*8 :: frac
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::A2D,B2D,dummy
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::tlca,tlcb,data1

      integer :: J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
C
      if (jdlast == 0) then   ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1                ! imon=January
        if (modelEclock%dayOfYear() <= 16)  then ! DAYOFYEAR in Jan 1-15, first month is Dec
          do L=1,Ldim*11
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          CALL REWIND_PARALLEL( iu )
        else              ! DAYOFYEAR is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (modelEclock%dayOfYear() > idofm(imon) .AND. 
     *      imon <= INT_MONTHS_PER_YEAR) go to 120
          do L=1,Ldim*(imon-2)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          if (imon == 13)  CALL REWIND_PARALLEL( iu )
        end if
      else                         ! Do we need to read in second month?
        if (modelEclock%dayOfYear() /= jdlast+1) then ! Check that data is read in daily
          if (modelEclock%dayOfYear() /= 1 .OR. 
     &      jdlast /= INT_DAYS_PER_YEAR) then
            write(out_line,*)'Bad day values in read_monthly_3Dsources'
     &      //': JDAY,JDLAST=',modelEclock%dayOfYear(),JDLAST
            call write_parallel(trim(out_line),crit=.true.)
            call stop_model('Bad values in read_monthly_3Dsources',255)
          end if
          imon=imon-INT_MONTHS_PER_YEAR             ! New year
          go to 130
        end if
        if (modelEclock%dayOfYear() <= idofm(imon)) go to 130
        imon=imon+1                ! read in new month of data
        if (imon == 13) then
          CALL REWIND_PARALLEL( iu  )
        else
          do L=1,Ldim*(imon-1)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
        endif
      end if
      DO L=1,Ldim
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      END DO
 130  continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-modelEclock%dayOfYear())
     &     / (idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     &     tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      return
      end SUBROUTINE read_mon3Dsources

      SUBROUTINE READ_OFFHNO3(OUT)
      use resolution, only: im,jm,lm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only : JDendOFM
      USE DOMAIN_DECOMP_ATM, only : grid,am_i_root
      IMPLICIT NONE
      include 'netcdf.inc'
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),intent(out) :: OUT
!@param  nlevnc vertical levels of off-line data  - 4x5 model=23, 2x2.5 model = 40
      INTEGER, PARAMETER :: nlevnc =40
      REAL*4, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,nlevnc) ::
     &     IN1_nohalo, IN2_nohalo
      REAL*8, DIMENSION(:,:,:), pointer, save :: IN1, IN2
!@var netcdf integer
      INTEGER :: ncid,id
      INTEGER, save :: step_rea=0, first_call=1
!@var time interpoltation
      REAL*8 :: tau
      integer start(4),count(4),status,l
      integer :: i_0,i_1,j_0,j_1
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c ----------------------------------------------------------------     

      if (first_call==1) then
        first_call=0
        allocate( IN1(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
        allocate( IN2(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
      endif
      if (step_rea.ne.modelEclock%month()) then 
        step_rea = modelEclock%month()
        if ( am_i_root() ) then
          print*,'READING HNO3 OFFLINE ',modelEclock%month(), step_rea
        endif
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
        status=NF_OPEN('OFFLINE_HNO3.nc',NCNOWRIT,ncid)
        status=NF_INQ_VARID(ncid,'FELD',id)
C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read
c   this is still latlon-specific.
c   will call read_dist_data for cubed sphere compatibility
c -----------------------------------------------------------------
        i_0 = grid%i_strt
        i_1 = grid%i_stop
        j_0 = grid%j_strt
        j_1 = grid%j_stop
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        start(4)=step_rea
        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=nlevnc
        count(4)=1

        status=NF_GET_VARA_REAL(ncid,id,start,count,IN1_nohalo)
        start(4)=step_rea+1
        if (start(4).gt.12) start(4)=1
        status=NF_GET_VARA_REAL(ncid,id,start,count,IN2_nohalo)

        status=NF_CLOSE(ncid)

        IN1(I_0:I_1,J_0:J_1,:) = IN1_nohalo(I_0:I_1,J_0:J_1,:)
        IN2(I_0:I_1,J_0:J_1,:) = IN2_nohalo(I_0:I_1,J_0:J_1,:)

      endif
C-----------------------------------------------------------------
      tau = (modelEclock%date()-.5)/(JDendOFM(modelEclock%month()) - 
     &     JDendOFM(modelEclock%month()-1))
         do l=1,lm
         OUT(:,:,l) = (1.-tau)*IN1(:,:,l)+tau*IN2(:,:,l)  
         enddo
c -----------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_OFFHNO3
c -----------------------------------------------------------------

      SUBROUTINE READ_OFFSS(OUT)
      use resolution, only: im,jm,lm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only : JDendOFM
      USE DOMAIN_DECOMP_ATM, only : grid,am_i_root
      IMPLICIT NONE
      include 'netcdf.inc'
!@param  nlevnc vertical levels of off-line data  
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),intent(out) :: OUT
      INTEGER, PARAMETER :: nlevnc =23
      REAL*4, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,nlevnc) ::
     &     IN1_nohalo, IN2_nohalo
      REAL*8, DIMENSION(:,:,:), pointer, save :: IN1_ss, IN2_ss
!@var netcdf integer
      INTEGER :: ncid,id
      INTEGER, save :: step_rea_ss=0, first_call_ss=1
!@var time interpoltation
      REAL*8 :: tau
      integer start(4),count(4),status,l,i,j
      integer :: i_0,i_1,j_0,j_1
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c ----------------------------------------------------------------     

      if (first_call_ss==1) then
        first_call_ss=0
        allocate( IN1_ss(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc))
        allocate( IN2_ss(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc))
      endif
      if (step_rea_ss.ne.modelEclock%month()) then 
        step_rea_ss = modelEclock%month()
        if ( am_i_root() ) then
          print*,'READING SEAS OFFLINE ',modelEclock%month(),step_rea_ss
        endif
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
        status=NF_OPEN('OFFLINE_SEAS.nc',NCNOWRIT,ncid)
        status=NF_INQ_VARID(ncid,'FELD',id)
C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read
c   this is still latlon-specific.
c   will call read_dist_data for cubed sphere compatibility
c -----------------------------------------------------------------
        i_0 = grid%i_strt
        i_1 = grid%i_stop
        j_0 = grid%j_strt
        j_1 = grid%j_stop
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        start(4)=step_rea_ss

        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=nlevnc
        count(4)=1

        status=NF_GET_VARA_REAL(ncid,id,start,count,IN1_nohalo)
        start(4)=step_rea_ss+1
        if (start(4).gt.12) start(4)=1
        status=NF_GET_VARA_REAL(ncid,id,start,count,IN2_nohalo)

        status=NF_CLOSE(ncid)

        IN1_ss(I_0:I_1,J_0:J_1,:) = IN1_nohalo(I_0:I_1,J_0:J_1,:)
        IN2_ss(I_0:I_1,J_0:J_1,:) = IN2_nohalo(I_0:I_1,J_0:J_1,:)

      endif

C-----------------------------------------------------------------------
      tau = (modelEclock%date()-.5)/(JDendOFM(modelEclock%month())
     *     - JDendOFM(modelEclock%month()-1))
         do l=1,lm
         OUT(:,:,l) = (1.-tau)*IN1_ss(:,:,l)+tau*IN2_ss(:,:,l)  
         enddo
c -----------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_OFFSS
c -----------------------------------------------------------------

      SUBROUTINE read_DMS_sources(swind,itype,i,j,DMS_flux) !!! T
!@sum generates DMS ocean source
!@auth Koch
c Monthly DMS ocean concentration sources are read in and combined
c  with wind and ocean temperature functions to get DMS air surface
c  concentrations
c want kg DMS/m2/s
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE GEOM, only: axyp
      use OldTracer_mod, only: tr_mm
      USE TRACER_COM, only: n_DMS,OFFLINE_DMS_SS
      use resolution, only: lm
      use model_com, only: modelEclock
      USE AEROSOL_SOURCES, only: DMSinput,DMS_AER
#ifdef old_DMS_emis
      USE FLUXES, only: GTEMP
#endif
      implicit none
      integer jread
      REAL*8 akw,erate,SCH,SCHR !!! T,Tc
#ifdef old_DMS_emis
      real*8 Tc ! YHL - FOR another DMS source
#endif
      real*8, PARAMETER :: E1=0.17d0
      real*8, PARAMETER :: E2=2.85d0
      real*8, PARAMETER :: E3=0.612d0
      real*8, PARAMETER :: E4=5.9d0
      real*8, PARAMETER :: E5=26.79d0
      real*8, PARAMETER :: E6=0.612d0
      real*8, PARAMETER :: SCHT=600.d0
      real*8, INTENT(OUT) :: DMS_flux
      real*8, INTENT(IN) :: swind
      integer, INTENT(IN) :: itype,i,j

      DMS_flux=0.d0
        erate=0.d0
        !!!Tc=T-273.d0
        if (OFFLINE_DMS_SS.ne.1) then
        if (itype.eq.1) then
c       if (lm.lt.40) then 
#ifndef old_DMS_emis
c Nightingale et al
        akw = 0.23d0*swind*swind + 0.1d0 * swind
        akw = akw * 0.24d0
        erate=akw*DMSinput(i,j,modelEclock%month())*1.d-9*62.d0 !*tr_mm(nt)
     *       /SECONDS_PER_DAY
#endif

#ifdef old_DMS_emis
c YUNHA - Liss and Merlivat (1986) code is from GISS GCM II'. 
c Liss and Merlivat (1986), use for > lm=40 to moderate DMS flux

       Tc=GTEMP(1,1,I,J) ! YUNHA GTEMP is already Celcius. 

       SCH=2674.d0-147.12d0*Tc+3.726d0*Tc*Tc-0.038d0*Tc*Tc*Tc
       IF(Tc.gt.47.) print*,'BAD_TEMPERATURE_DMS',i,j,Tc,SCH
       SCHR=SCHT/SCH
       if (swind.lt.3.6) then
        akw=0.041*(SCHR)**(2.d0/3.d0)*swind
       else if (swind.lt.13.) then
         akw=(0.68*SWIND - 2.31)*DSQRT(SCHR)
       else
       akw=(1.42*SWIND - 11.8)*DSQRT(SCHR)
       endif  !swind
       erate=akw*DMSinput(i,j,modelEclock%month())*1.d-9*62.d0/
     *      SECONDS_PER_DAY     !not sure of units

#endif
c       if (lm.ge.40) erate=erate/5.d0   !I think there was an error in input files
c       else
c Liss and Merlivat (1986), use for > lm=40 to moderate DMS flux
c       SCH=2674.d0-147.12d0*Tc+3.726d0*Tc*Tc-0.038d0*Tc*Tc*Tc
c       SCHR=SCHT/SCH
c       if (swind.lt.3.6) then
c       akw=E1*(SCHR)**(2.d0/3.d0)*swind
c       else if (swind.lt.13.) then
c       akw=E2*DSQRT(SCHR)*(swind-3.6d0)
c    *     +E3*(SCHR)**(2.d0/3.d0)
c       else
c       akw=E4*(swind-13.d0)*DSQRT(SCHR)+E5*(swind-3.6d0)*
c    *      DSQRT(SCHR)+E6*(SCHR)**(2.d0/3.d0)
c       endif  !swind
c       erate=akw*DMSinput(i,j,jmon)*1.d-9/sday !not sure of units
c       endif ! lm
        endif !itype
        else !AEROCOM run, prescribed flux
c if after Feb 28 skip the leapyear day
         jread=modelEclock%dayOfYear()
         if (modelEclock%dayOfYear().gt.59) 
     *        jread=modelEclock%dayOfYear()+1
c         if (j.eq.1.or.j.eq.46) DMS_AER(i,j,jread)
c     *      =DMS_AER(i,j,jread)*72.d0
         erate=DMS_AER(i,j,jread)/SECONDS_PER_DAY/axyp(i,j)*
     &         tr_mm(n_DMS)/32.d0
        endif
        DMS_flux=erate          ! units are kg/m2/s
c
      return
      end SUBROUTINE read_DMS_sources

      SUBROUTINE read_seasalt_sources(swind,itype,ibin,i,j,ss,tr)
!@sum determines wind-speed dependent oceanic seasalt source
!@auth Koch
c want kg seasalt/m2/s, for now in 2 size bins
      USE TRACER_COM, only: OFFLINE_DMS_SS,OFFLINE_SS
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE GEOM, only: axyp
      use model_com, only: modelEclock
      USE AEROSOL_SOURCES, only: SS1_AER,SS2_AER,tune_ss1,tune_ss2
#ifdef TRACERS_AEROSOLS_OCEAN
     &                          ,OC_SS_enrich_fact
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_TOMAS
      USE TOMAS_EMIS, only : scalesizeSalt
#endif
      !USE FLUXES, only: gtemp !Jaegle
      use Dictionary_mod, only: sync_param
      implicit none
      REAL*8 erate,swind_cap
      integer jread,nb
      integer, INTENT(IN)::itype,ibin,i,j
      REAL*8, INTENT(IN)::swind
      REAL*8, INTENT(OUT)::ss
      character*8, intent(in) :: tr
c
      ss=0.
      erate=0.d0
#ifndef TRACERS_TOMAS
      if (OFFLINE_DMS_SS.ne.1.and.OFFLINE_SS.ne.1) then
        if (itype.eq.1) then
c Monahan 1971, bubble source, important for small (<10um) particles
!          swind_cap=swind !modelE
!          if (swind.gt.10.d0) swind_cap=10.d0 !modelE
!          erate= 1.373d0 * swind_cap**(3.41d0) !modelE
!          erate=1.373d0*swind**3.41d0 !Monahan
!          erate=(swind/10.d0)**2.5d0 !Lewis and Schwartz
!          erate=(0.3d0+0.1d0*gtemp(1,1,i,j)-0.0076d0*gtemp(1,1,i,j)**2+
!     *           0.00021d0*gtemp(1,1,i,j)**3)*1.373d0*swind**3.41d0 !Jaegle
          erate=1.373d0*swind**3.41d0 !Gong
          if (ibin.eq.1) then ! submicron (0.1 < r_d < 1.)
!            ss=tune_ss1*erate*2.11d-14 !modelE
!            ss=tune_ss1*erate*1.750298d-14 !Monahan
!            ss=tune_ss1*erate*3.431418d-11 !Lewis and Schwartz
!            ss=tune_ss1*erate*1.336825d-14 !Jaegle
            ss=tune_ss1*erate*1.336825d-14 !Gong
#ifdef TRACERS_AEROSOLS_OCEAN
            if (trim(tr).eq.'OCocean') then
              ss=ss*OC_SS_enrich_fact(i,j)
            else
              ss=ss*(1.d0-OC_SS_enrich_fact(i,j))
            endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
          else ! supermicron (1. < r_d < 4.)
!            ss=tune_ss2*erate*7.78d-14 !modelE
!            ss=tune_ss2*erate*7.097854d-14 !Monahan
!            ss=tune_ss2*erate*1.650674d-9 !Lewis and Schwartz
!            ss=tune_ss2*erate*8.676763d-14 !Jaegle
            ss=tune_ss2*erate*8.676763d-14 !Gong
          endif
c     units are kg salt/m2/s
        endif
      else
c if after Feb 28 skip the leapyear day
        jread=modelEclock%dayOfYear()
        if (modelEclock%dayOfYear().gt.59) 
     &       jread=modelEclock%dayOfYear()+1
        if (ibin.eq.1) then
          ss=SS1_AER(i,j,jread)/(SECONDS_PER_DAY*axyp(i,j))
#ifdef TRACERS_AEROSOLS_OCEAN
          if (trim(tr).eq.'OCocean') then
            ss=ss*OC_SS_enrich_fact(i,j)
          else
            ss=ss*(1.d0-OC_SS_enrich_fact(i,j))
          endif
#endif  /* TRACERS_AEROSOLS_OCEAN */
        else 
          ss=SS2_AER(i,j,jread)/(SECONDS_PER_DAY*axyp(i,j))
        endif
      endif
#endif /* TRACERS_TOMAS */

#ifdef TRACERS_TOMAS 
      if (OFFLINE_DMS_SS.ne.1.and.OFFLINE_SS.ne.1) then
        if (itype.eq.1) then
          erate=1.373d0*swind**3.41d0 !Gong
          ss=tune_ss1*erate*scalesizeSalt(ibin)
!#ifdef TRACERS_AEROSOLS_OCEAN
!            if (trim(tr).eq.'OCocean') then
!              ss=ss*OC_SS_enrich_fact(i,j)
!            else
!              ss=ss*(1.d0-OC_SS_enrich_fact(i,j))
!            endif
!#endif  /* TRACERS_AEROSOLS_OCEAN */
         endif
       else
C TOMAS don't allow offline SS! 
        call stop_model('TOMAS: NO offline Sea-salt emission',255)
       endif

#endif 


      return
      end SUBROUTINE read_seasalt_sources


      SUBROUTINE aerosol_gas_chem
!@sum aerosol gas phase chemistry
!@vers 2013/03/27
!@auth Dorothy Koch
      use OldTracer_mod, only: trname, tr_mm
      use TRACER_COM, only: ntm, oh_live, no3_live, trm
      use TRACER_COM, only: coupled_chem, n_BCIA, n_BCII, n_DMS,n_H2O2_s
      use TRACER_COM, only: rsulf1, rsulf2, rsulf3, rsulf4
      use TRACER_COM, only: n_MSA, N_OCII, n_OX, n_SO2, n_OCIA
      use TRACER_COM, only: n_SO4, n_SO4_d1, n_SO4_d2, n_SO4_d3
      use TRACER_COM, only: nChemistry, nChemLoss
#if (defined TRACERS_HETCHEM) || (defined TRACERS_NITRATE)
      use TRACER_COM, only: rxts1, rxts2, rxts3
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AECOB, n_AECIL, n_AOCOB, n_AOCIL
      use TRACER_COM, only: n_AECIL, n_H2SO4, nbins
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only: n_H2SO4
#endif
      USE TRDIAG_COM, only : 
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     *     jls_OHconk,jls_HO2con,jls_NO3,jls_phot
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ,jls_OHcon
#endif
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT, getDomainBounds 
      USE DOMAIN_DECOMP_ATM, only: DREAD8_PARALLEL,DREAD_PARALLEL
      USE DOMAIN_DECOMP_ATM, only : GRID, write_parallel
      USE RESOLUTION, only : ls1
      use resolution, only: im,jm,lm
      use atm_com, only : t,q
      use model_com, only: modelEclock
      USE MODEL_COM, only: dtsrc
      USE ATM_COM, only: pmid,MA,pk,LTROPO,byMA
      USE PBLCOM, only : dclev
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      USE AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,
     & dho2,perj,tno3,ohsr,o3_offline, JmonthCache,
     &      ohrCache, dho2rCache, perjrCache, tno3rCache
      USE CONSTANT, only : mair
      use TimeConstants_mod, only: SECONDS_PER_DAY
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: which_trop
#endif
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only : h2so4_chem
#endif
#ifdef TRACERS_AEROSOLS_VBS
      use CONSTANT, only : gasc
      use TRACERS_VBS, only: vbs_tracers, vbs_conditions, 
     &                       vbs_calc, vbs_tr
#endif /* TRACERS_AEROSOLS_VBS */
c Aerosol chemistry
      implicit none
      logical :: ifirst=.true.
      real*8 ppres,te,tt,mm,dmm,ohmc,r1,d1,r2,d2,ttno3,r3,d3,
     * ddno3,dddms,ddno3a,fmom,dtt
      real*8 rk4,ek4,r4,d4
      real*8 r6,d6,ek9,ek9t,ch2o,eh2o,dho2mc,dho2kg,eeee,xk9,
     * r5,d5,dmssink,bdy
#ifdef TRACERS_HETCHEM
     *       ,d41,d42,d43,o3mc,rsulfo3
#endif
      real*8 bciage,ociage
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: ohsr_in
      integer i,j,l,n,iuc,iun,itau,ixx1,ixx2,ichemi,itt,
     * ittime,isp,iix,jjx,llx,ii,jj,ll,iuc2,it,najl,j_0,j_1,
     * j_0s,j_1s,mmm,J_0H,J_1H,I_0,I_1
#ifdef TRACERS_SPECIAL_Shindell
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
#endif
#ifdef TRACERS_AEROSOLS_VBS
      type(vbs_tracers) :: vbs_tr_old ! concentrations, ug m-3
      type(vbs_conditions) :: vbs_cond ! current box conditions (meteo+chem)
!@var kg2ugm3 factor to convert kilograms gridbox-1 to ug m-3
      real*8 :: kg2ugm3
#endif /* TRACERS_AEROSOLS_VBS */
      integer maxl,nrecs_skip
      logical :: newMonth
      save ifirst
#ifdef TRACERS_TOMAS
      REAL*8 TAU_hydro
      integer k
#endif

      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1,
     *    J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,J_STRT_SKP=J_0S,
     * J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** initialise source arrays
        tr3Dsource(:,j_0:j_1,:,1,n_DMS)=0. ! DMS chem sink
#ifndef TRACERS_AMP
#ifndef TRACERS_TOMAS
        tr3Dsource(:,j_0:j_1,:,1,n_MSA)=0. ! MSA chem sink
        tr3Dsource(:,j_0:j_1,:,1,n_SO4)=0. ! SO4 chem source
#endif
#endif
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO2)=0. ! SO2 chem source
        tr3Dsource(:,j_0:j_1,:,nChemloss,n_SO2)=0. ! SO2 chem sink
        if(n_H2O2_s>0) tr3Dsource(:,j_0:j_1,:,1,n_H2O2_s)=0. ! H2O2 chem source
        if(n_H2O2_s>0) tr3Dsource(:,j_0:j_1,:,2,n_H2O2_s)=0. ! H2O2 chem sink
#ifdef TRACERS_AMP
        tr3Dsource(:,j_0:j_1,:,2,n_H2SO4)=0. ! H2O2 chem sink
#endif
#ifdef TRACERS_TOMAS
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_H2SO4)=0. ! H2O2 chem sink
        H2SO4_chem(:,j_0:j_1,:)=0.0
        do k=1,nbins
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AECOB(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AECIL(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AOCOB(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AOCIL(k))=0.
        enddo
#endif
#ifdef TRACERS_HETCHEM
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d1) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d2) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,1,n_SO4_d3) =0. ! SO4 on dust
#endif
        if (n_BCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_BCII)=0. ! BCII sink
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_BCIA)=0. ! BCIA source
        end if
        if (n_OCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_OCII)=0. ! OCII sink
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_OCIA)=0. ! OCIA source
        end if
#ifdef TRACERS_AEROSOLS_VBS
        tr3Dsource(:,j_0:j_1,:,nChemistry,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nChemloss,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nOther,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nChemistry,vbs_tr%iaer)=0.
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
C Coupled mode: use on-line radical concentrations
      if (coupled_chem.eq.1) then
        oh(:,j_0:j_1,:)=oh_live(:,j_0:j_1,:)
        tno3(:,j_0:j_1,:)=no3_live(:,j_0:j_1,:)
c Set h2o2_s =0 and use on-line h2o2 from chemistry
        if(n_H2O2_s>0) trm(:,j_0:j_1,:,n_h2o2_s)=0.0
      endif

      if (coupled_chem.eq.0) then
c Use this for chem inputs from B4360C0M23, from Drew
c      if (ifirst) then
        newMonth = jMonthCache /= modelEclock%month()
        if (newMonth) then
          jMonthCache = modelEclock%month()
          call openunit('AER_CHEM',iuc,.true.)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),ohrCache,
     &         recs_to_skip=5*(modelEclock%month()-1)+1)    ! 5 recs/month + ichemi for this month
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),dho2rCache)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),perjrCache)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),tno3rCache)
          call closeunit(iuc)
        end if
        ohr   = ohrCache  
        dho2r = dho2rCache
        perjr = perjrCache
        tno3r = tno3rCache
        if (im.eq.72) then
        call openunit('AER_OH_STRAT',iuc2,.true.)
        nrecs_skip=lm*(modelEclock%month()-1) ! skip all the preceding months
        do ll=1,lm
          call DREAD_PARALLEL(grid,iuc2,nameunit(iuc2),ohsr_in,
     &       recs_to_skip=nrecs_skip)
          ohsr(:,:,ll)=ohsr_in(:,:)*1.D5
          nrecs_skip=0 ! do not skip any more records
        enddo
        call closeunit(iuc2)

c skip poles because there was a bug in the input file over the pole
        do j=j_0s,j_1s   
        do i=i_0,i_1
          maxl=ltropo(i,j)
        do l=maxl,lm
          ohr(i,j,l)=ohsr(i,j,l)
        end do
        end do
        end do
cdmk turning these off because I don't have 2x2.5
c       ifirst=.false.
        else    !need to scale inputs (10^5 mol/cm3)
        ohr(:,:,:)=ohr(:,:,:)*1.D5
        tno3r(:,:,:)=tno3r(:,:,:)*1.D5
        dho2r(:,:,:)=dho2r(:,:,:)*1.D7
        perjr(:,:,:)=perjr(:,:,:)*1.D2
        endif   !im.eq.72
c I have to read in every timestep unless I can find a better way
c
c impose diurnal variability
        CALL SCALERAD
c       write(6,*) ' RRR OXID2 ',ohr(10,45,1),
c    *   oh(10,45,1),dho2r(3,45,1),dho2(3,45,1)
       endif   !coupled_chem.eq.0

C Calculation of gas phase reaction rates
C Now called from tracer_3Dsource
c#ifndef  TRACERS_SPECIAL_Shindell  
c        CALL GET_SULF_GAS_RATES
c#endif

#ifdef TRACERS_HETCHEM
c calculation of heterogeneous reaction rates: SO2 on dust 
      CALL SULFDUST
c calculation of heterogeneous reaction rates: SO2 on seasalt
c      CALL SULFSEAS 
c     if (COUPLED_CHEM.ne.1) then
c     CALL GET_O3_OFFLINE
c     endif
#endif
#endif
      dtt=dtsrc
      !efold time of 2.7 days
      bciage=(1.d0-exp(-dtsrc/(2.7d0*SECONDS_PER_DAY)))/dtsrc 
      !efold time of 1.6 days
      ociage=(1.d0-exp(-dtsrc/(1.6d0*SECONDS_PER_DAY)))/dtsrc
C**** THIS LOOP SHOULD BE PARALLELISED
      do 20 l=1,lm
      do 21 j=j_0,j_1
      do 22 i=i_0,imaxj(j)
C Initialise       
c       maxl = ltropo(i,j)
        bdy = dclev(i,j) 
#ifdef TRACERS_SPECIAL_Shindell
        if(which_trop.eq.1)maxl=ls1-1
#endif
c     if(l.le.maxl) then
      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm = MA(l,i,j)*axyp(i,j)
      tt = 1.d0/te

c DMM is number density of air in molecules/cm3
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3

      do 23 n=1,NTM ! ===== THIS IS CHEMISTRY OF Koch AEROSOLS =====

        select case (trname(n))
c    Aging of industrial carbonaceous aerosols 
        case ('BCII')
          tr3Dsource(i,j,l,nChemistry,n)=-bciage*trm(i,j,l,n)
          tr3Dsource(i,j,l,nChemistry,n_BCIA)=bciage*trm(i,j,l,n)

#ifdef TRACERS_AEROSOLS_VBS
        case ('vbsAm2') ! This handles all VBS tracers
          kg2ugm3=1.d9*(1.d2*pmid(l,i,j))*mair/
     &            (MA(l,i,j)*axyp(i,j)*gasc*te)
          vbs_cond%dt=dtsrc
          vbs_cond%OH=ohmc
          vbs_cond%temp=te
          vbs_cond%nvoa=(trm(i,j,l,n_BCII)
     &                  +trm(i,j,l,n_BCIA)
     &                  +trm(i,j,l,n_BCB)
#ifdef TRACERS_AEROSOLS_SOA
     &                  +trm(i,j,l,n_isopp1a)
     &                  +trm(i,j,l,n_isopp2a)
     &                  +trm(i,j,l,n_apinp1a)
     &                  +trm(i,j,l,n_apinp2a)
#endif /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
     &                  +trm(i,j,l,n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
     &                  +trm(i,j,l,n_msa)
     &                  +trm(i,j,l,n_so4)
#ifdef TRACERS_NITRATE
     &                  +trm(i,j,l,n_nh4)
     &                  +trm(i,j,l,n_no3p)
#endif
     &                  )*kg2ugm3
          vbs_tr_old%gas=trm(i,j,l,vbs_tr%igas)*kg2ugm3
          vbs_tr_old%aer=trm(i,j,l,vbs_tr%iaer)*kg2ugm3

          call vbs_calc(vbs_tr_old,vbs_cond)

          tr3Dsource(i,j,l,nChemistry,vbs_tr%igas)=
     &      vbs_tr%chem_prod/kg2ugm3/vbs_cond%dt
          tr3Dsource(i,j,l,nChemloss,vbs_tr%igas)=
     &      vbs_tr%chem_loss/kg2ugm3/vbs_cond%dt
          tr3Dsource(i,j,l,nOther,vbs_tr%igas)=
     &      -vbs_tr%partition/kg2ugm3/vbs_cond%dt ! partitioning
          tr3Dsource(i,j,l,nChemistry,vbs_tr%iaer)=
     &      vbs_tr%partition/kg2ugm3/vbs_cond%dt
!     &      (vbs_tr%gas-vbs_tr_old%gas)/kg2ugm3/vbs_cond%dt
!      if (sum(vbs_tr_old%gas)+sum(vbs_tr_old%aer) /= 0.) then
!        print '(a,3e)','KOSTAS gas',
!     &                 sum(vbs_tr_old%gas),
!     &                 sum(vbs_tr%gas),
!     &                 sum(vbs_tr_old%gas)+sum(vbs_tr_old%aer)
!        print '(a,3e)','KOSTAS aer',
!     &                 sum(vbs_tr_old%aer),
!     &                 sum(vbs_tr%aer),
!     &                 sum(vbs_tr%gas)+sum(vbs_tr%aer)
!        print '(a,3e)','KOSTAS bud',
!     &                 sum(vbs_tr%chem_prod),
!     &                 sum(vbs_tr%chem_loss),
!     &                 sum(vbs_tr%partition)
!      endif
#else
        case ('OCII')
          tr3Dsource(i,j,l,nChemistry,n)=-ociage*trm(i,j,l,n)
          tr3Dsource(i,j,l,nChemistry,n_OCIA)=ociage*trm(i,j,l,n)
#endif /* TRACERS_AEROSOLS_VBS */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
        case ('DMS')
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2

          r1=rsulf1(i,j,l)*ohmc 
          d1 = exp(-r1*dtsrc)
          r2=rsulf2(i,j,l)*ohmc
          d2 = exp(-r2*dtsrc)

c     NO3 is in mixing ratio: convert to molecules/cm3
c - not necessary for Shindell source
          if (l.gt.bdy) then
            ttno3=0.d0
            go to 87
          endif
          ttno3 = tno3(i,j,l)   !*6.02d20*ppres/(.082056d0*te)
 87       r3=rsulf3(i,j,l)*ttno3
          d3= exp(-r3*dtsrc)
          ddno3=r3*trm(i,j,l,n)/tr_mm(n)*1000.d0*dtsrc
          dddms=trm(i,j,l,n)/tr_mm(n)*1000.d0
          if (ddno3.gt.dddms) ddno3=dddms

          ddno3=ddno3*0.9
C DMS losses: eqns 1, 2 ,3

          tr3Dsource(i,j,l,1,n) = trm(i,j,l,n)*(d1*d2-1.)/dtsrc

          dmssink=ddno3*tr_mm(n)/1000.d0

          if (dmssink.gt.trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc)
     *         dmssink=trm(i,j,l,n)+tr3Dsource(i,j,l,1,n)*dtsrc
          tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) - dmssink/dtsrc
          
        case ('MSA')
C MSA gain: eqn 1

          tr3Dsource(i,j,l,1,n) = 0.25d0*Tr_mm(n)/Tr_mm(n_dms)*trm(i,j
     *         ,l,n_dms)*(1.d0 -D1)*SQRT(D2)/dtsrc
          
        case ('SO2')
c SO2 production from DMS
          tr3Dsource(i,j,l,nChemistry,n) = (0.75*tr_mm(n)/tr_mm(n_dms)
     *         *trm(i,j,l
     *         ,n_dms)*(1.d0 - d1)*sqrt(d2)+ tr_mm(n)/tr_mm(n_dms)*trm(i
     *         ,j,l,n_dms)*(1.d0 - d2)*sqrt(d1)+dmssink*tr_mm(n)
     *         /tr_mm(n_dms))/dtsrc
#ifdef TRACERS_TOMAS 
! EC/OC aging 
        case ('AECIL_01')
           TAU_hydro=1.5D0*SECONDS_PER_DAY !24.D0*3600.D0 !1.5 day 

           DO K=1,nbins
              tr3Dsource(i,j,l,nChemistry,n_AECIL(K))=
     &             trm(i,j,l,n_AECOB(K))*
     &             (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 

              tr3Dsource(i,j,l,nChemistry,n_AECOB(K))=
     &             -trm(i,j,l,n_AECOB(K))*
     &            (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 

!              IF(am_i_root())
!         print*,'ECOB aging',k,n_AECOB(K),(1.D0-EXP(-dtsrc/TAU_hydro))
!     &  , trm(i,j,l,n_AECOB(K))
           ENDDO
           
        case ('AOCIL_01')
           TAU_hydro=1.5D0*SECONDS_PER_DAY !24.D0*3600.D0 !1.5 day 

           DO K=1,nbins
              tr3Dsource(i,j,l,nChemistry,n_AOCIL(K))
     &             =trm(i,j,l,n_AOCOB(K))*
     &             (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc
!     &             4.3D-6
              tr3Dsource(i,j,l,nChemistry,n_AOCOB(K))
     &             =-trm(i,j,l,n_AOCOB(K))*
     &            (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 
!     &             4.3D-6
!              IF(am_i_root())
!         print*,'OCOB aging',k,n_AOCOB(K),(1.D0-EXP(-dtsrc/TAU_hydro))
!     &  , trm(i,j,l,n_AOCOB(K))

           ENDDO   

#endif
#ifndef TRACERS_TOMAS                   
          najl = jls_NO3
          call inc_tajls2(i,j,l,najl,ttno3)
#endif
#endif
        end select
        
 23   CONTINUE ! ===== END OF CHEMISTRY OF Koch AEROSOLS ====
c       endif
 22   CONTINUE
 21   CONTINUE
 20   CONTINUE

      do 30 l=1,lm
      do 31 j=j_0,j_1
      do 32 i=i_0,imaxj(j)

c     maxl = ltropo(i,j)
#ifdef TRACERS_SPECIAL_Shindell
      if(which_trop.eq.1)maxl=ls1-1
#endif
c     if(l.le.maxl) then

      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm = MA(l,i,j)*axyp(i,j)
      tt = 1.d0/te
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3
#ifdef TRACERS_HETCHEM
      if (COUPLED_CHEM.ne.1) then
      o3mc=o3_offline(i,j,l)*dmm*(28.0D0/48.0D0)*BYAXYP(I,J)*byMA(L,I,J)
      else
        o3mc=trm(i,j,l,n_Ox)*dmm*(28.0D0/48.0D0)*BYAXYP(I,J)*byMA(L,I,J)
      endif
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      do 33 n=1,NTM
        select case (trname(n))
        case ('SO2')
c oxidation of SO2 to make SO4: SO2 + OH -> H2SO4

          r4=rsulf4(i,j,l)*ohmc
          d4 = exp(-r4*dtsrc)

          IF (d4.GE.1.) d4=0.99999d0
#ifdef TRACERS_HETCHEM
      rsulfo3 = 4.39d11*exp(-4131/te)+( 2.56d3*exp(-966/te)) * 10.d5 !assuming pH=5
      rsulfo3 = exp(-rsulfo3*o3mc *dtsrc) !O3 oxidation Maahs '83
       d41 = exp(-rxts1(i,j,l)*dtsrc)     
       d42 = exp(-rxts2(i,j,l)*dtsrc)     
       d43 = exp(-rxts3(i,j,l)*dtsrc)     
       tr3Dsource(i,j,l,nChemloss,n) = (-trm(i,j,l,n)*(1.d0-d41)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d4)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d42)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d43)/dtsrc)
#else
       tr3Dsource(i,j,l,nChemloss,n) = -trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#ifdef TRACERS_AMP
       tr3Dsource(i,j,l,2,n_H2SO4)=trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#endif  
#ifdef TRACERS_TOMAS
       H2SO4_chem(i,j,l)=trm(i,j,l,n)*(1.d0-d4)/dtsrc 
     &      *tr_mm(n_H2SO4)/tr_mm(n) 
#endif      
#endif        
c diagnostics to save oxidant fields
c No need to accumulate Shindell version here because it
c   is done elsewhere
c#ifdef TRACERS_SPECIAL_Shindell
c         najl = jls_OHcon
c#else
          najl = jls_OHconk
c#endif
          if (najl > 0) call inc_tajls2(i,j,l,najl,oh(i,j,l))
          najl = jls_HO2con
          if (najl > 0) call inc_tajls2(i,j,l,najl,dho2(i,j,l))

#ifdef TRACERS_HETCHEM
       case ('SO4_d1')
c sulfate production from SO2 on mineral dust aerosol due to O3 oxidation

       tr3Dsource(i,j,l,1,n)=tr3Dsource(i,j,l,1,n)+tr_mm(n)/tr_mm(n_so2)
     *         *(1.d0-d41)*trm(i,j,l,n_so2)            !  SO2
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d2')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d42)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d3')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,1,n) = tr3Dsource(i,j,l,1,n) +  tr_mm(n)/
     *             tr_mm(n_so2)*(1.d0-d43)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc

#endif
        case('SO4')
C SO4 production
          tr3Dsource(i,j,l,nChemistry,n) = 
     *         tr3Dsource(i,j,l,nChemistry,n)+tr_mm(n)
     *         /tr_mm(n_so2)*trm(i,j,l,n_so2)*(1.d0 -d4)/dtsrc
        case('H2O2_s')

          if (coupled_chem.ne.1) then

c hydrogen peroxide formation and destruction:
C***5.H2O2 +hv -> 2OH
C***6.H2O2 + OH -> H2O + HO2
C***7.H2O2 + SO2 -> H2O + SO3 (in-cloud, in CB)
C***9.HO2 + HO2 ->H2O2 + O2
C     HO2 + HO2 + M ->
C     HO2 + HO2 + H2O ->
C     HO2 + HO2 + H2O + M ->

          r6 = 2.9d-12 * exp(-160.d0*tt)*ohmc
          d6 = exp(-r6*dtsrc)
          ek9 = 2.2d-13*exp(600.d0*tt)
          ek9t = 1.9d-20*dmm*0.78d0*exp(980.d0*tt)*1.d-13
          ch2o = q(i,j,l)*6.02d20*28.97d0/18.d0*ppres/(.082d0*te)
          eh2o = 1.+1.4d-21*exp(2200.d0*tt)*ch2o
          dho2mc = dho2(i,j,l)  !/mm*1.292/.033*6.02e17
          dho2kg = dho2(i,j,l)*mm*te*.082056d0/(ppres*28.97d0*6.02d20)
          eeee = eh2o*(ek9+ek9t)*dtt*dho2mc
          xk9 = dho2kg*eeee
c         if (i.eq.2.and.l.eq.1.and.j.eq.46) write(6,*) 
c    *    'RRR CHEM DEBUG ',i,j,xk9,dho2kg,eeee,dho2mc
c    *    ,eh2o,ek9,ek9t,dtt
c         if (i.eq.72.and.l.eq.1.and.j.le.46) write(6,*) 
c    *    'RRR CHEM DEBUG ',i,j,xk9,dho2kg,eeee,dho2mc
c H2O2 production: eqn 9
         
          tr3Dsource(i,j,l,1,n) = tr_mm(n)*xk9/dtsrc
c        if (i.eq.10.and.j.eq.45.and.l.eq.1) then
c        write(6,*) 'RRR OXID H2O2',xk9,dho2kg,eeee
c         endif
c H2O2 losses:5 and 6
          r5 = perj(i,j,l)
          d5 = exp(-r5*dtsrc)

          tr3Dsource(i,j,l,2,n)=(trm(i,j,l,n))*(d5*d6-1.d0)
     *         /dtsrc
          
          najl = jls_phot
          if (najl > 0) call inc_tajls(i,j,l,najl,perj(i,j,l))
          endif
        end select

 33   CONTINUE
#endif
c     endif
 32   CONTINUE
 31   CONTINUE
 30   CONTINUE

      RETURN
      END SUBROUTINE aerosol_gas_chem

      SUBROUTINE SCALERAD
      use constant, only : pi
      use resolution, only: im,jm,lm
      use AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,dho2,perj,tno3
      USE DOMAIN_DECOMP_ATM, only:GRID, getDomainBounds
      use RAD_COM, only: cosz1,cosz_day,sunset
      implicit none
      real*8, parameter ::
     &     night_frac_min=.01d0 ! minimum night_frac for tno3 scaling
      real*8 stfac,night_frac
      integer i,j,l,i_0,i_1,j_0,j_1
c      real*8, DIMENSION(JM) :: tczen
c      integer, DIMENSION(JM) :: nradn

      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1, J_STRT=J_0,J_STOP=J_1)

c      nradn(:)=0
c      tczen(:)=0.d0
c      do 100 j = j_0,j_1
c      do 100 i = 1, im
c      if (cosz1(i,j).gt.0.) then
c      tczen(j)=tczen(j)+cosz1(i,j)
c      else
c      nradn(j)=nradn(j)+1
c      endif
c 100  continue

      do j = j_0,j_1
      do i = i_0,i_1
c Get NO3 only if dark, weighted by number of dark hours
        night_frac = 1.-sunset(i,j)/pi
c        night_frac = real(nradn(j))/real(im)
        if (cosz1(i,j).le.0.and.night_frac.gt.night_frac_min) then
          tno3(i,j,:)=tno3r(i,j,:)/night_frac !DMK jmon
        else
          tno3(i,j,:)=0.d0
        endif
        if (cosz1(i,j).gt.0.) then
c          stfac=cosz1(i,j)/tczen(j)*real(IM)
          stfac=cosz1(i,j)/cosz_day(i,j)
          oh(i,j,:)=ohr(i,j,:)*stfac
          perj(i,j,:)=perjr(i,j,:)*stfac
          dho2(i,j,:)=dho2r(i,j,:)*stfac
        end if
      end do
      end do

c        if (I.EQ.1.AND.L.EQ.1) write(6,*)'NO3R',TAU,J,NRADN(I,J)
c            if (l.eq.1.and.j.eq.23.and.i.eq.10) write(6,*)
c    *     'RRR SCALE ',stfac,cosz1(i,j),tczen(j),oh(i,j,l),ohr(i,j,l)
      RETURN
      END SUBROUTINE SCALERAD



      SUBROUTINE GET_SULFATE(L,temp_in,fcloud,
     *  wa_vol,wmxtr,sulfin,sulfinc,sulfout,tr_left,
     *  tm,tmcl,airm,LHX,dt_sulf,fcld0,no_plume)

!@sum  GET_SULFATE calculates formation of sulfate from SO2 and H2O2
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch
c
C**** GLOBAL parameters and variables:
      USE CONSTANT, only: BYGASC, MAIR,teeny,mb2kg,gasc,LHE
      use OldTracer_mod, only: trname, mass2vol, tr_mm
      use OldTracer_mod, only: tr_RKD, tr_DHD
      USE TRACER_COM, only: n_H2O2_s,n_SO2
     *     ,NTM
     *     ,lm,n_SO4,n_H2O2,coupled_chem
      USE CLOUDS, only: PL,NTIX,NTX,DXYPIJ
      USE MODEL_COM, only: dtsrc
c
      IMPLICIT NONE
c
C**** Local parameters and variables and arguments:
c
!@param BY298K unknown meaning for now (assumed= 1./298K)
!@var Ppas pressure at current altitude (in Pascal=kg/s2/m)
!@var TFAC exponential coeffiecient of tracer condensation temperature
!@+   dependence (mole/joule)
!@var FCLOUD fraction of cloud available for tracer condensation
!@var SSFAC dummy variable (assumed units= kg water?)
!@var L index for altitude loop
!@var N index for tracer number loop
!@var RKD dummy variable (= tr_RKD*EXP[ ])
!@var WA_VOL is the cloud water volume in L
!@var WMX_INC is the change in cloud water ratio
!@var CLWC is the cloud liquid water content: L water/L air
!@var sulfin is the amount of SO2 and H2O2 used to make sulfate
!@var sulfout is the amount of sulfate generated
!@var tr_left is the amount of SO2 and H2O2 left after sulfate is made
!@+    and is now available to condense
!@var PPH is the partial pressure of the gas in M/kg
!@var amass is airmass in kg airm(l)*mb2kg
!@var trd is dissolved portion in moles/L
!@var sulfinc is the change in dissolved SO2 and H2O2 as we use those
!@+  to form sulfate (in addition to what is dissolved)
      REAL*8, PARAMETER :: BY298K=3.3557D-3
      REAL*8 Ppas, tfac, ssfac, RKD, Henry_const(ntm)
      real*8 clwc,rk1f,rkdm(ntm),amass,trd(ntm),trdr(ntm),
     * dso4g,dso4d,pph(ntm),trmol(ntm),trdmol(ntm),dso4gt,dso4dt
      integer n,ih,is,is4,ihx,isx,is4x
      real*8, parameter :: rk1=1.3d-2 !M
      real*8, parameter :: dh1=-1.6736d4 !J/mol
      real*8, parameter :: rk=6.357d14    !1/(M*M*s)
      real*8, parameter :: ea=3.95d4 !J/mol
      REAL*8,  INTENT(IN) :: fcloud,temp_in,wa_vol,wmxtr,LHX
      real*8, dimension(lm,ntm) :: tm
      real*8, dimension(ntm) :: tmcl,sulfin,sulfout,tr_left
     *  ,sulfinc
      real*8, dimension(lm) :: airm
c     REAL*8,  INTENT(OUT)::
!@var dt_sulf accumulated diagnostic of sulfate chemistry changes
      real*8, dimension(ntm), intent(inout) :: dt_sulf
      real*8 finc,fcld0,temp
      INTEGER, INTENT(IN) :: L
!@var NO_PLUME : false ==> TM is already TMP (TM in plume), so fcloud shouldn't 
!@+   be applied tm (added by Yunha Lee, 2012)
      LOGICAL NO_PLUME 
      do n=1,ntx
        sulfin(N)=0.
        sulfinc(N)=0.
        sulfout(N)=0.
        tr_left(N)=1.
      end do
c is this needed?
c#if (defined TRACERS_COSMO) || (defined TRACERS_
c      do n=1,ntx
c       select case (trname(ntix(n)))
c       case ('Pb210   ','Be7     ','Be10    ','Rn222')
c       go to 333
c       end select
c      end do
C#endif
c
c
C**** CALCULATE the fraction of tracer mass that becomes condensate:
c
      if (LHX.NE.LHE.or.fcloud.lt.teeny) go to 333
      finc=(fcloud-fcld0)/fcloud
      if (finc.lt.0d0) finc=0.d0
c First allow for formation of sulfate from SO2 and H2O2. Then remaining
c  gases may be allowed to dissolve (amount given by tr_left)
C H2O2 + SO2 -> H2O + SO3 -> H2SO4
      amass=airm(l)*mb2kg*DXYPIJ
      Ppas = PL(L)*1.D2
      ! calls to this subroutine are sometimes made at stages of
      ! the cloud scheme at which some but not all tendencies have
      ! been applied to temp_in, so we impose a lower limit
      ! (liquid water is very unlikely to exist below 230 K)
      temp = max(temp_in, 230d0) ! K
      tfac = (1./temp - by298k)*bygasc  !mol/J
c  cloud liquid water content
      clwc=wmxtr*mair*ppas/temp*bygasc/1.D6/fcloud
      rk1f=rk1*exp(-dh1*tfac)

      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO2')
       is=ntix(n)
       isx=n
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(is)=tr_rkd(is)*(1.+ rk1f/3.2d-5)
c mole of tracer, used to limit so4 production
      if(no_plume)then
        trmol(is)=1000.*tm(l,isx)/tr_mm(is)*fcloud
      else
        trmol(is)=1000.*tm(l,isx)/tr_mm(is)
      endif
c partial pressure of gas x henry's law coefficient
      pph(is)=mass2vol(is)*1.d-3*ppas/amass*
     *   tr_rkd(is)*exp(-tr_dhd(is)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
      henry_const(is)=rkdm(is)*exp(-tr_dhd(is)*tfac)
      pph(is)=pph(is)/(1+(henry_const(is)*clwc*gasc*temp))
c again all except tmcl(n)
      if(CLWC==0.)then
        trdr(is)=0.d0
      else
        trdr(is)=mass2vol(is)*ppas/amass*bygasc
     &  /temp*1.D-3/CLWC  !"CLWC" is newly added by Yunha Lee. This makes "mole/liter of cloud water/kg of tracer"  
      end if
c dissolved moles
      trdmol(is)=tmcl(isx)*1000./tr_mm(is)  !trdr(is) is replaced with tmcl(isx) by Yunha Lee. This makes "trdmol" to be the amount of tracer already dissolved in the cloud water. 

       case('H2O2','H2O2_s')

         if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 400
         if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto
     *        400

       ih=ntix(n)
       ihx=n
c modified Henry's Law coefficient assuming pH of 4.5
      rkdm(ih)=tr_rkd(ih)
c mole of tracer, used to limit so4 production
      if(no_plume)then
      trmol(ih)=1000.*tm(l,ihx)/tr_mm(ih)*fcloud
      else
      trmol(ih)=1000.*tm(l,ihx)/tr_mm(ih)
      endif
c partial pressure of gas x henry's law coefficient
      pph(ih)=mass2vol(ih)*1.D-3*ppas/amass*
     *   tr_rkd(ih)*exp(-tr_dhd(ih)*tfac)
c the following is from Phil:
c      reduction in partial pressure as species dissolves
c yhl - revisit the pph calculation in future - is "*clwc*gasc*temp" needed?
      henry_const(ih)=rkdm(ih)*exp(-tr_dhd(ih)*tfac)
      pph(ih)=pph(ih)/(1+(henry_const(ih)*clwc*gasc*temp))
c all except tmcl(n)
      if(CLWC==0.) then
        trdr(ih)=0.d0
      else
        trdr(ih)=mass2vol(ih)*ppas/amass*bygasc/temp*1.D-3/CLWC  !"CLWC" is newly added by Yunha Lee. This makes "mole/liter of cloud water/kg of tracer"   !M/kg
      end if
c dissolved moles
      trdmol(ih)=tmcl(ihx)*1000./tr_mm(ih)  !trdr(ih) is replaced with tmcl(ih) by Yunha Lee. This makes "trdmol" to be the amount of tracer already dissolved in the cloud water.

 400   CONTINUE

      end select
      end do
      if (tm(l,ihx).lt.teeny.or.tm(l,isx).lt.teeny) then
      dso4g=0.
      go to 21
      endif
c this part from gas phase:moles/kg/kg
      dso4g=rk*exp(-ea/(gasc*temp))*rk1f
     *    *pph(ih)*pph(is)*dtsrc*wa_vol
c dmk  incremental water volume
      dso4g=dso4g*finc
c should probably be (finc+tr_lef) but then tr_lef has to be saved   
c check to make sure no overreaction: moles of production:
      dso4gt=dso4g*tm(l,ihx)*tm(l,isx)
c can't be more than moles going in:
      if (dso4gt.gt.trmol(is)) then
        dso4g=trmol(is)/(tm(l,ihx)*tm(l,isx))
      endif
      dso4gt=dso4g*tm(l,ihx)*tm(l,isx)
      if (dso4gt.gt.trmol(ih)) then
        dso4g=trmol(ih)/(tm(l,ihx)*tm(l,isx))
      endif
c this part from dissolved gases
 21    dso4d=rk*exp(-ea/(gasc*temp))*rk1f
     *    *trdr(ih)*trdr(is)*dtsrc*wa_vol  !mole/kg/kg

      if (tmcl(ihx).lt.teeny.or.tmcl(isx).lt.teeny) then   !"trdr" is replaced by "tmcl" by YUNHA LEE.
      dso4d=0.
      go to 22
      endif

c check to make sure no overreaction: moles of production:
      dso4dt=dso4d*tmcl(ihx)*tmcl(isx)  !"trdr [M/kg]" is replaced by "tmcl [kg]" by YUNHA LEE. This makes "dso4dt" to be "mole of SO4 formed". 
c can't be more than moles going in:
      if (dso4dt.gt.trdmol(is)) then
        if(tmcl(ihx)==0. .or. tmcl(isx)==0.)then
          dso4d=0.d0
        else
          dso4d=trdmol(is)/(tmcl(ihx)*tmcl(isx))  !"trdr" is replaced by "tmcl" by YUNHA LEE.
        end if
      endif
      dso4dt=dso4d*tmcl(ihx)*tmcl(isx) !"trdr" is replaced by "tmcl" by YUNHA LEE.
      if (dso4dt.gt.trdmol(ih)) then
        if(tmcl(ihx)==0. .or. tmcl(isx)==0.)then
          dso4d=0.d0
        else
          dso4d=trdmol(ih)/(tmcl(ihx)*tmcl(isx))   !"trdr" is replaced by "tmcl" by YUNHA LEE.
        end if
      endif

 22   continue
      do n=1,ntx
       select case (trname(ntix(n)))
       case('SO4','M_ACC_SU','ASO4__01')
       is4=ntix(n)

       sulfout(is4)=tr_mm(is4)/1000.*(dso4g*tm(l,isx)*tm(l,ihx)
     *  +dso4d*tmcl(isx)*tmcl(ihx)) !kg

       dt_sulf(is4) = dt_sulf(is4) + sulfout(is4)

       case('SO2')
       is=ntix(n)
       isx=n

! is ih/ihx set here, then why isn't is/isx?
       sulfin(is)=-dso4g*tm(l,ihx)*tr_mm(is)/1000. !dimnless
       sulfinc(is)=-dso4d*tmcl(ihx)*tr_mm(is)/1000.
       sulfinc(is)=max(-1d0,sulfinc(is))
       sulfin(is)=max(-1d0,sulfin(is))
       tr_left(isx)=0.
       if (fcloud.gt.abs(sulfin(is))) then
         tr_left(isx)=(fcloud+sulfin(is))
       endif

       dt_sulf(is)=
     & dt_sulf(is)+sulfin(is)*tm(l,isx)+sulfinc(is)*tmcl(isx) !trdr(is) is replaced by tmcl(isx) by YUNHA LEE. 

       case('H2O2','H2O2_s')

         if (trname(ntix(n)).eq."H2O2" .and. coupled_chem.eq.0) goto 401
         if (trname(ntix(n)).eq."H2O2_s" .and. coupled_chem.eq.1) goto
     *        401

       ih=ntix(n)
       ihx=n
       sulfin(ih)=-dso4g*tm(l,isx)*tr_mm(ih)/1000.
       sulfinc(ih)=-dso4d*tmcl(isx)*tr_mm(ih)/1000.
       sulfinc(ih)=max(-1d0,sulfinc(ih))
       sulfin(ih)=max(-1d0,sulfin(ih))
       tr_left(ihx)=0.
       if (fcloud.gt.abs(sulfin(ih))) then
         tr_left(ihx)=fcloud+sulfin(ih)
       endif


 401   CONTINUE

       dt_sulf(ih)=
     & dt_sulf(ih)+sulfin(ih)*tm(l,ihx)+sulfinc(ih)*tmcl(ihx) !trdr(ih) is replaced by tmcl(isx) by YUNHA LEE.

      end select
      END DO

 333  RETURN
      END SUBROUTINE GET_SULFATE

      SUBROUTINE GET_BC_DALBEDO(i,j,bc_dalb)
!@sum Calculates change to albedo of snow on ice and snow on land due
!@+     to BC within the snow.
!@+     Parameterization based on Warren and Wiscombe (1980) (21 inputs)
!@+     or actually on Flanner et al. fig 2 r_e=500 (14 inputs)
!@+     or Warren and Wiscombe (1985) (18 input, old vs new, then I
!@+     continue linearly from 19-29 off the plot)
!@+auth Dorothy Koch
c
      USE CONSTANT, only: rhow
      USE GHY_COM, only: tr_wsn_ij, wsn_ij
      USE SEAICE_COM, only : si_atm
      use TRACER_COM, only:
#ifdef TRACERS_AEROSOLS_Koch
      use TRACER_COM, only: n_BCB,n_BCII,n_BCIA
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only:n_M_BC1_BC,n_M_BC2_BC,n_M_BC3_BC,n_M_DBC_BC
     *  ,n_M_BOC_BC,n_M_BCS_BC,n_M_MXX_BC
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AECIL,n_AECOB,nbins
#endif
      !USE VEG_COM, only: afb
      USE RADPAR, only: agesn
      USE FLUXES, only: atmice
      IMPLICIT NONE
c Warren and Wiscombe 1985 includes age dependence
      real*8, parameter :: bc(29)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
     * 6.d0,7.d0,8.d0,9.d0,10.d0,20.d0,30.d0,40.d0,50.d0,60.d0,
     * 70.d0,80.d0,90.d0,100.d0,110.d0,120.d0,130.d0,140.d0,
     * 150.d0,160.d0,170.d0,180.d0,190.d0,200.d0/)
      real*8, parameter :: daln(29)=(/0.d0,0.1d0,0.1d0,0.2d0,
     * 0.2d0,0.2d0,0.2d0,0.3d0,0.3d0,0.4d0,0.7d0,0.9d0,1.1d0,
     * 1.3d0,1.5d0,1.6d0,1.8d0,2.d0,2.2d0,2.4d0,2.6d0,2.8d0,
     * 3.d0,3.2d0,3.4d0,3.6d0,3.8d0,4.d0,4.2d0/)
      real*8, parameter :: dalo(29)=(/0.1d0,0.2d0,0.4d0,0.5d0,
     * 0.6d0,0.7d0,0.8d0,0.9d0,1.d0,1.d0,2.d0,2.6d0,3.2d0,
     * 3.8d0,4.3d0,4.8d0,5.2d0,5.5d0,5.9d0,6.3d0,6.7d0,7.1d0,
     * 7.5d0,7.9d0,8.3d0,8.7d0,9.1d0,9.5d0,9.9d0/)
c Flanner et al
c     real*8, parameter :: bc(14)=(/25.d0,50.d0,100.d0,150.d0,
c    * 200.d0,
c    *250.d0,300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,
c    * 1000.d0/)
c     real*8, parameter :: dal(14)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
c    * 6.d0,7.d0,8.d0,9.d0,10.d0,11.d0,11.5d0,12.d0,12.5d0/)
c Warren and Wiscomb 1980
c     real*8, parameter :: bc(21)=(/0.05d0,0.075d0,0.1d0,0.2d0,
c    * 0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.d0,2.d0,
c    * 3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0/)
c     real*8, parameter :: dal(21)=(/2.d0,3.d0,4.d0,6.d0,8.d0,
c    * 9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,16.d0,18.d0,20.d0,
c    * 22.d0,24.d0,26.d0,28.d0,30.d0,32.d0,34.d0/)
      REAL*8  scon,icon,bcsnowb,bcsnowv,bcice,fv,fb,
     * sconb,sconv,bcc,rads
      INTEGER n,ic,ib
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: bc_dalb
#ifdef TRACERS_AEROSOLS_Koch
      integer, parameter :: nspBC=3
#endif
#ifdef TRACERS_AMP
      integer, parameter :: nspBC=7
#endif
#ifdef TRACERS_TOMAS
      integer, parameter :: nspBC=nbins+nbins
#endif
      integer, dimension(nspBC) :: spBC
c
c tr_wsn_ij(n,nsl,2,i,j) tracer in snow layer l multiplied by fraction snow, kg/m2
c wsn_ij(nsl,2,i,j)
c trsi(n,nsi,i,j) tracer in sea ice in layer l, kg/m2
c snowi(i,j) snow amount on sea ice, kg/m2
c afb(i,j)=fb, fraction that is bare soil
c fv=1-fb fraction that is vegetated
c rads is the snow grain size determined in GRAINS
c gtracer(n,i,j) is tracer concentration in snow on sea ice?
c Maybe I need tracer in snow on sea ice, or mass of sea ice...?
c Does trsi accumulate for ALL tracers?
c fractions??
c
#ifdef TRACERS_AEROSOLS_Koch
      spBC(1)=n_BCII
      spBC(2)=n_BCIA
      spBC(3)=n_BCB
#endif
#ifdef TRACERS_AMP
      spBC(1)=n_M_BC1_BC
      spBC(2)=n_M_BC2_BC
      spBC(3)=n_M_BC3_BC
      spBC(4)=n_M_DBC_BC
      spBC(5)=n_M_BOC_BC
      spBC(6)=n_M_BCS_BC
      spBC(7)=n_M_MXX_BC
#endif
#ifdef TRACERS_TOMAS
      do n=1,nbins
         spBC(n)=n_AECOB(n)
         spBC(n+nbins)=n_AECIL(n)
      enddo
#endif
      bc_dalb=0.
      scon=0.
      sconb=0.
      sconv=0.
      icon=0.
      bcsnowb=0.
      bcsnowv=0.
      bcice=0.
      !fb=afb(i,j)
      !fv=1.-fb
      sconb=0.d0
      sconv=0.d0
      call get_fb_fv( fb, fv, i, j )
      if (wsn_ij(1,1,i,j).gt.0.) then
        do n=1,nspBC
          bcsnowb=bcsnowb+tr_wsn_ij(spBC(n),1,1,i,j)
        enddo
        sconb=bcsnowb/wsn_ij(1,1,i,j)/rhow
      endif
      if (wsn_ij(1,2,i,j).gt.0.) then
        do n=1,nspBC
          bcsnowv=bcsnowv+tr_wsn_ij(spBC(n),1,2,i,j)
        enddo
        sconv=bcsnowv/wsn_ij(1,2,i,j)/rhow
      endif
      scon=(fb*sconb+fv*sconv)*1.D9   !kg/kg to ppmw
      if (si_atm%snowi(i,j).gt.0.) then
        do n=1,nspBC
          icon=icon+atmice%gtracer(spBC(n),i,j)*1.d9
        enddo
      endif
      bcc=DMAX1(icon,scon)
      call GRAINS(i,j,rads)

      do ib=1,28
        if (bcc.gt.bc(ib).and.bcc.lt.bc(ib+1)) then
          bc_dalb=-(daln(ib)
     *       +(rads-100.d0)/900.d0*(dalo(ib)-daln(ib)))/100.
          go to 33
        endif
      end do
 33   continue
      if (bcc.ge.bc(29)) bc_dalb=-(daln(29)
     *    + (rads-100.d0)/900.d0*(dalo(29)-daln(29)))/100.
c     if (bc_dalb.ne.0.) write(6,*) 'alb_write',i,j,bc_dalb,bcc,rads
      RETURN
      END SUBROUTINE GET_BC_DALBEDO

      SUBROUTINE GRAINS(i,j,rads)
!@sum Estimates snow grain size (microns) based on air temperature
!@+     and snow age. From Susan Marshall's PhD thesis
!@+auth Dorothy Koch
c
      USE CONSTANT, only: pi,gasc
      USE FLUXES, only: atmsrf
      use TimeConstants_mod, only: DAYS_PER_YEAR
      USE GHY_COM, only: snoage
      USE AEROSOL_SOURCES, only: snosiz
      IMPLICIT none
      REAL*8 E,A,age,r0,radmm,ert,
     * tfac,area,delrad
      INTEGER n,ic,ib
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: rads
      DATA E,A /26020.d0, 29100.d0/
c tsavg(i,j) surface air temperature
c snoage(k,i,j) k=1 ocean ice, k=2 land ice, k=3 land
c   (do these differ within a gridbox?) in days
c RN1 radius from previous timestep
c RADS snow grain radius
c
c Find the age of snow, I assume the age does not
c  vary within the gridbox so just take the max?
       age=DMAX1(snoage(1,i,j),snoage(2,i,j),snoage(3,i,j))
c Use Temperature to check if melting or non-melting snow
       IF (atmsrf%tsavg(i,j).le.273.15) then
c Non-melting snow; distinguish between initial or
c  secondary growth rate
        IF (age.lt.13.5) then
c initial growth
         r0=50.
         radmm=r0+ (0.008d0+age)
         rads=radmm*1000.d0
         snosiz(i,j)=rads+r0
        ELSE
c secondary growth
         r0=150.d0
         ert = dexp(-E/(GASC*atmsrf%TSAVG(I,J)))
         tfac = a*ert
         area = (TFAC/DAYS_PER_YEAR) * (AGE-12.5d0)
         radmm=dsqrt(area/pi)
         delrad = radmm*1000.d0
         rads=delrad + r0
         snosiz(i,j)=rads
        ENDIF
       ELSE
c melting snow
        radmm=dsqrt(0.05d0)
        rads=snosiz(i,j)+(radmm*100.d0)
        snosiz(i,j)=rads
       ENDIF
       rads=DMIN1(rads,1000.d0)
       rads=DMAX1(rads,100.d0)
      RETURN
      END SUBROUTINE GRAINS
      
      SUBROUTINE read_mon_3D
     & (Ldim,iu,data1,trans_emis,yr1,yr2)
!@sum Read in monthly sources and interpolate to current day
!@auth Jean Lerner and others / Greg Faluvegi
! taken from TRACERS_SPECIAL_Shindell, in case we
!  we run aerosols independent of gases
      use resolution, only: im,jm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only: idofm=>JDmidOfM
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL
     &     ,REWIND_PARALLEL,write_parallel,backspace_parallel,am_i_root
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,ipos,k,nn
      character(len=300) :: out_line
      real*8 :: frac, alpha
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::A2D,B2D,dummy
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::tlca,tlcb,data1
     *     ,sfc_a,sfc_b
      logical, intent(in):: trans_emis
      integer, intent(in):: yr1,yr2
      integer :: kstep=10 !<< please note this hard-code

      integer :: J_0, J_1
      integer :: year, dayOfYear

      call modelEclock%getDate(year=year, dayOfYear=dayOfYear)
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

C No doubt this code can be combined/compressed, but I am going to
C do the transient and non-transient cases separately for the moment:

! -------------- non-transient emissions ----------------------------!
      if(.not.trans_emis) then
C
      imon=1                ! imon=January
      if (dayOfYear <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        call rewind_parallel(iu)
      else              ! DAYOFYEAR is in Jan 16 to Dec 16, get first month
        do while(dayOfYear > idofm(imon) .AND. 
     *    imon <= INT_MONTHS_PER_YEAR)
          imon=imon+1
        enddo
        if(imon/=2)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:'
          call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13) call rewind_parallel(iu)
      end if
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-dayOfYear)/(idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      write(out_line,*) '3D source monthly factor=',frac
      call write_parallel(trim(out_line))

! --------------- transient emissions -------------------------------!
      else
        ipos=1
        alpha=0.d0 ! before start year, use start year value
        if(year>yr2.or.(year==yr2.and.dayOfYear>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(yr2-yr1)/kstep
        endif
        do k=yr1,yr2-kstep,kstep
          if(year>k .or. (year==k.and.dayOfYear>=183)) then
            if(year<k+kstep .or. (year==k+kstep.and.dayOfYear<183))then
              ipos=1+(k-yr1)/kstep ! (integer artithmatic)
              alpha=real(year-k)/real(kstep)
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
      imon=1                ! imon=January
      if (dayOfYear <= 16)  then ! JDAY in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! JDAY is in Jan 16 to Dec 16, get first month
        do while(dayOfYear > idofm(imon) .AND. 
     &    imon <= INT_MONTHS_PER_YEAR)
          imon=imon+1
        enddo
        if(imon/=2 .or. ipos/=1)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:'
          call readt_parallel
     &    (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCC write(6,*) 'Not using this first record:'
CCCCC call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
      frac = float(idofm(imon)-dayOfYear)/(idofm(imon)-idofm(imon-1))
      sfc_a(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      call rewind_parallel( iu )

      ipos=ipos+1
      imon=1                ! imon=January
      if (dayOfYear <= 16)  then ! DAYOFYEAR in Jan 1-15, first month is Dec
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! JDAY is in Jan 16 to Dec 16, get first month
        do while(dayOfYear > idofm(imon) .AND. 
     &    imon <= INT_MONTHS_PER_YEAR)
          imon=imon+1
        enddo
        write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCCCwrite(6,*) 'Not using this first record:'
CCCCCCcall readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      enddo
      frac = float(idofm(imon)-dayOfYear)/(idofm(imon)-idofm(imon-1))
      sfc_b(:,J_0:J_1,:) =
     & tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)

! now interpolate between the two time periods:

      data1(:,J_0:J_1,:) =
     & sfc_a(:,J_0:J_1,:)*(1.d0-alpha) + sfc_b(:,J_0:J_1,:)*alpha

      write(out_line,*) '3D source at',
     &100.d0*alpha,' % of period this day ',k,' to this day ',k+kstep,
     &' and monthly fraction= ',frac
      call write_parallel(trim(out_line))

      endif ! transient or not

      return
      end SUBROUTINE read_mon_3D
      

#ifdef TRACERS_AEROSOLS_OCEAN
      subroutine read_seawifs_chla(imon)
!@sum read_seawifs_chla Reads in SeaWiFS chlorophyll-a concentration
!@auth Kostas Tsigaridis

      use filemanager, only: openunit,closeunit,nameunit
      use domain_decomp_atm, only: grid,get,readt_parallel
      use aerosol_sources, only: OC_SS_enrich_fact
      implicit none

      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &  seawifs_chla
      integer :: iu, imon
      integer :: j,i,j_0,j_1,i_0,i_1

      call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
      call getDomainBounds(grid, i_strt=i_0, i_stop=i_1)

! read chla SeaWiFS data
      call openunit('SeaWiFS_chla',iu,.true.,.true.)
      call readt_parallel(grid,iu,nameunit(iu),seawifs_chla,imon)
      call closeunit(iu)

! SS enrichment factor of OC (Vignati et al., AE, 2009, 670)
! if chla is more than 1.43 mg m-3 the enrichment factor remains
! constant, 76%.
      do j=j_0,j_1; do i=i_0,i_1
        OC_SS_enrich_fact(i,j)=min(seawifs_chla(i,j),1.43d0)*
     &                         0.435d0+0.13805d0
      enddo ; enddo

      end subroutine read_seawifs_chla
#endif  /* TRACERS_AEROSOLS_OCEAN */
