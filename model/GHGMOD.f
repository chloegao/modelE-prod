#include "rundeck_opts.h"

      module O3mod
!@sum O3mod administers reading of ozone files
!@auth M. Kelley and original development team
      use timestream_mod, only : timestream
      implicit none
      save
!@var O3stream interface for reading and time-interpolating O3 files
!@+   See usage notes in timestream_mod
      type(timestream) :: O3stream,delta_O3stream
#ifdef HIGH_FREQUENCY_O3_INPUT
      type(timestream) :: OxHFstream,PSFforO3stream
#endif

!@dbparam use_sol_Ox_cycle if =1, a cycle of ozone is appled to
!@+ o3year, as a function of the solar constant cycle.
      integer :: use_sol_Ox_cycle = 0
      real*8 :: S0min, S0max

!@dbparam ozone_use_ppm_interp = 1 uses ppm interpolation in the
!@+ timestream. Otherwise uses linm2m.
      integer :: ozone_use_ppm_interp = 1

!@var have_o3_file whether an O3file was specified in the rundeck
      logical :: have_o3_file

!@param NLO3_traditional assumed number of layers in ozone data files.
      integer, parameter :: NLO3_traditional = 49
!@var NLO3 number of layers in ozone data files, as read from file.
      integer :: NLO3 = 0
!@var PLBO3_traditional assumed edge pressures in O3 input file.
      real*8 :: PLBO3_traditional(NLO3_traditional+1) = (/
     *       984d0, 934d0, 854d0, 720d0, 550d0, 390d0, 285d0, 210d0,
     *       150d0, 125d0, 100d0,  80d0,  60d0,  55d0,  50d0,
     *        45d0,  40d0,  35d0,  30d0,  25d0,  20d0,  15d0,
     *       10.d0,  7.d0,  5.d0,  4.d0,  3.d0,  2.d0,  1.5d0,
     *        1.d0,  7d-1,  5d-1,  4d-1,  3d-1,  2d-1,  1.5d-1,
     *        1d-1,  7d-2,  5d-2,  4d-2,  3d-2,  2d-2,  1.5d-2,
     *        1d-2,  7d-3,  5d-3,  4d-3,  3d-3,  1d-3,  1d-7/)
!@var PLBO3 edge pressures in O3 input file, as read from file
      real*8, allocatable :: PLBO3(:)

      REAL*8, dimension(:,:,:), pointer :: o3jday,o3jref
#ifdef HIGH_FREQUENCY_O3_INPUT
      REAL*8, dimension(:,:,:), pointer :: o3jday_HF_modelLevels
#endif
      INTEGER :: use_o3_ref=0

      contains

      subroutine UPDO3D(JYEARO,JJDAYO,O3JDAY,O3JREF)
      use dictionary_mod
      use resolution, only : psf
      use domain_decomp_atm, only: grid, getdomainbounds
      use timestream_mod, only : init_stream,read_stream,
     &     getname_firstfile
      use pario, only : par_open,par_close,read_dist_data,
     & variable_exists,get_dimlen,read_data
      use filemanager, only : file_exists
      implicit none
      integer, intent(in) :: JYEARO,JJDAYO
      real*8, dimension(:,:,:), pointer :: o3jday,o3jref

      integer :: i,j,l,jyearx,fid
      logical, save :: init = .false.
      logical :: cyclic,exists
      real*8, allocatable :: o3arr(:,:,:)
      character(len=6) :: method
      character(len=32) :: fname1st

      integer :: j_0, j_1, i_0, i_1

      call getdomainbounds(grid, j_strt=j_0,j_stop=j_1,
     &                           i_strt=i_0,i_stop=i_1)

      jyearx = abs(jyearo)

      if (.not. init) then
        init = .true.

        have_o3_file = file_exists('O3file')

        if(have_o3_file) then

          ! Initialize the timestream for the O3 data file:
          cyclic = jyearo < 0

          call sync_param("ozone_use_ppm_interp",ozone_use_ppm_interp)
          if(ozone_use_ppm_interp==1)then
            method = 'ppm'
          else
            method = 'linm2m'
          endif
          call init_stream(grid,O3stream,'O3file','O3',
     &         0d0,1d30,trim(method),jyearx,jjdayo,cyclic=cyclic)
          ! query the layering
          call getname_firstfile(O3stream,fname1st)
          fid = par_open(grid,trim(fname1st),'read')
          if(variable_exists(grid,fid,'ple'))then
            nlo3=get_dimlen(grid,fid,'ple') - 1 ! coord var but one less
            if(nlo3.ne.get_dimlen(grid,fid,'plm'))call
     &        stop_model('ple/plm dim problem in O3file',255)
            allocate(plbo3(nlo3+1))
            call read_data(grid,fid,'ple',plbo3,bcast_all=.true.)
          else
            call stop_model('missing ple info in o3file',255)
          endif
          call par_close(grid,fid)
        else
          nlo3=nlo3_traditional
          allocate(plbo3(nlo3+1))
          plbo3(:)=plbo3_traditional(:)
        endif

        allocate(o3jday(nlo3,grid%i_strt:grid%i_stop,
     &                       grid%j_strt:grid%j_stop))
        o3jday = 0.

        allocate(o3jref(nlo3_traditional,grid%i_strt:grid%i_stop,
     &                                   grid%j_strt:grid%j_stop))
        o3jref = 0.

! The next line is brought over from the original UPDO3D. I think
! it is to prevent "losing" some ozone in the REPART interpolation
! if the (fixed) lowest O3 level pressure is at lower pressure than
! the the (fixed) lowest nominal model pressure:
        if(plbo3(1) < psf) plbo3(1) = psf
        if(plbo3_traditional(1) < psf) plbo3_traditional(1) = psf

! Read the 3D field for O3 RCOMPX reference calls.
! (There is no need to allow for this one on flexible # of levels)
        inquire(file='Ox_ref',exist=exists)
        if(exists) then
          allocate(o3arr(grid%i_strt_halo:grid%i_stop_halo,
     &                   grid%j_strt_halo:grid%j_stop_halo,
     &                   nlo3_traditional))
          fid = par_open(grid,'Ox_ref','read')
          call read_dist_data(grid,fid,'O3',o3arr)
          call par_close(grid,fid)
          do j=j_0,j_1
          do i=i_0,i_1
            O3JREF(:,I,J)=O3ARR(I,J,:)
          enddo
          enddo
          deallocate(o3arr) ! note quick deallocation as will be resized below
        endif

      endif  ! end init

      if(have_o3_file)then
        allocate(o3arr(grid%i_strt_halo:grid%i_stop_halo, ! resizing
     &                 grid%j_strt_halo:grid%j_stop_halo,nlo3))

        call read_stream(grid,O3stream,jyearx,jjdayo,o3arr)
        do j=j_0,j_1
        do i=i_0,i_1
          O3JDAY(:,I,J)=O3ARR(I,J,:)
        enddo
        enddo
        deallocate(o3arr)
      endif ! there is an o3file

      return
      end subroutine UPDO3D


#ifdef HIGH_FREQUENCY_O3_INPUT
      subroutine UPDO3D_highFrequency(JYEARO,JJDAYO,
     & o3jday_HF_modelLevels)
      use resolution, only : LM, mtop,mfix,mfrac,mfixs
      use domain_decomp_atm, only: grid, getdomainbounds
      use timestream_mod, only : init_stream,read_stream
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      use atm_com, only: pedn
      use constant, only : bygrav,tf,rgas,mb2kg,kg2mb
      implicit none
      integer, intent(in) :: JYEARO,JJDAYO
      real*8, dimension(:,:,:), pointer :: o3jday_HF_modelLevels

      integer :: i,j,l,jyearx,fid
      logical, save :: init = .false.
      logical :: cyclic,exists
      real*8, allocatable :: OxHFarr(:,:,:)
      real*8, allocatable :: psf4o3arr(:,:)
      real*8, dimension(LM):: OxHFarr_Interpolated, OxHFarr_Converted,
     &                        airmass
      real*8, dimension(LM+1)::modelPressureBottoms, filePressureBottoms
      real*8 :: numerator, denominator, mvar

      integer :: j_0, j_1, i_0, i_1

      call getdomainbounds(grid, j_strt=j_0,j_stop=j_1,
     &                           i_strt=i_0,i_stop=i_1)

      allocate(OxHFarr(  grid%i_strt_halo:grid%i_stop_halo,
     &                   grid%j_strt_halo:grid%j_stop_halo,LM))
      allocate(psf4o3arr(grid%i_strt_halo:grid%i_stop_halo,
     &                   grid%j_strt_halo:grid%j_stop_halo))

      jyearx = abs(jyearo)

      if (.not. init) then
        init = .true.
        allocate(o3jday_HF_modelLevels(LM,grid%i_strt:grid%i_stop,
     &                                    grid%j_strt:grid%j_stop))
        o3jday_HF_modelLevels = 0.

        cyclic = jyearo < 0

        call init_stream(grid,OxHFstream,'OxHFfile','Ox',
     &        0d0,1d30,'none',jyearx,jjdayo,cyclic=cyclic)
        call init_stream(grid,PSFforO3stream,'OxHFfile','p_surf',
     &        0d0,1d30,'none',jyearx,jjdayo,cyclic=cyclic)

      endif  ! end init

      call read_stream(grid,OxHFstream,jyearx,jjdayo,OxHFarr)
      call read_stream(grid,PSFforO3stream,jyearx,jjdayo,psf4o3arr)

      do j=j_0,j_1
      do i=i_0,i_1
        modelPressureBottoms(:)=pedn(:,i,j)
        ! approximate the air mass concurrent with ozone input:
        mvar = psf4o3arr(i,j)*mb2kg - mfixs - mtop
        filePressureBottoms(LM+1) = mtop
        do L=LM,1,-1
           airmass(L) = mfix(L) + mvar*mfrac(L)
           filePressureBottoms(L) = filePressureBottoms(L+1) +
                                    airmass(L)*kg2mb  ;  EndDo

        ! to avoid potentially losing some of the column ozone, adjust
        ! bottom level edge (similar to how routine UPDO3D does:
        ! if(plbo3(1) < psf) plbo3(1) = psf ) Though we are not altering
        ! the airmass at the same time. Should we?
        filePressureBottoms(1)=
     &   max(filePressureBottoms(1),modelPressureBottoms(1))

        ! Convert units from pppv (volume mixing ratio input) to atm-cm.
        ! For now using the approx. input file air mass, but this could
        ! be changed such that the vmr is the fundamental quantity and
        ! the mass is based on current model air mass...
        !
        ! Explanation of conversion:
        ! Get "numerator" which is the amount of ozone we are inputting
        ! in a given layer in units of kg(O3)/m2, where m2 is horizontal
        ! surface area. Then get a "denominator" that is the density of
        ! ozone at standard atmospheric pressure and 0 deg C in units of
        ! kg(O3)/m3. This is a constant.
        ! Then the num/den ratio has units of {kg/m2} / {kg/m3} = m and
        ! represents "how much" ozone you would have if it were at those
        ! pressure and temperature conditions, expressed as a thickness.
        ! Call that an "atm-m". Then in the end conversion to atm-cm or
        ! Dobson Unit are just powers of 10.
        !
        ! The numerator is obtained starting with read-in mole fraction
        ! and denote a mole as n, our starting units are: n(O3)/n(air).
        ! n(O3)/n(air) * [ratio of molecular weights, ozone to air] is:
        ! n(O3)/n(air) * [48. g(O3)/n(O3)  /  28.9655d g(air)/n(air)]
        !  --> g(O3)/g(air) = kg(O3)/kg(air). So now we have a mass
        ! mixing ratio. Then multiply by the air mass in kg/m2:
        ! kg(O3)/kg(air) * [kg(air)/m2] --> kg(O3)/m2.
        !
        ! The demoninator is obtained starting with the density of ozone
        ! at 1 atmosphere and 0 deg C. At those conditions, air density
        ! is p/RT. I.e. p, R, T are constants here and reference Earth's
        ! atmosphere: p=101325 Pa, R=rgas in J kg-1 K-1, T=tf in K,
        ! so units work out to kg(air)/m3. Convert from air to ozone
        ! again using the ratio of molecular weights, e.g. on Earth:
        ! 1.2922 kg(air)/m3 * [48. g(O3)/n(O3) / 28.9655d g(air)/n(air)]
        ! --> 2.1415 kg(O3)/m3. Note that this ratio of molecular weights
        ! appears in the numerator and denominator so is skipped below.
        !
        ! Now do numerator/denominator and obtain atm-m units, and
        ! multiply by 100 to get the desired atm-cm units. This is the
        ! cm thickness of O3 one would have under those those specific
        ! atmoserpheric conditions. All that results in just:

        do L=1,LM
          numerator=OxHFarr(i,j,L)*airmass(L) ! kg O3 / m2 we have
          denominator=101325.d0/(rgas*tf)     ! kg O3 / m3 @ 1 atm and 0 deg C
          OxHFarr_Converted(L)=1.d2*numerator/denominator
        enddo

        ! Now, interpolate vertically, but this interpolation is not onto
        ! the rad code O3 levels, it is just an adjustment over the same
        ! LM levels but allowing for different surface pressure than was
        ! concurrent when this model input was saved from a previous run:
        call repart(OxHFarr_Converted, filePressureBottoms,  LM+1,   ! IN
     &           OxHFarr_Interpolated, modelPressureBottoms, LM+1)   ! OUT
        ! save for use in rad code proper:
        o3jday_HF_modelLevels(:,i,j)=OxHFarr_Interpolated(:)
      enddo
      enddo

      deallocate(OxHFarr, psf4o3arr)

      return
      end subroutine UPDO3D_highFrequency
#endif /* HIGH_FREQUENCY_O3_INPUT */


      SUBROUTINE UPDO3D_solar(jjdayo,S0,o3jday)
!@sum UPDO3D_solar adds solar cycle variability to O3JDAY
      use dictionary_mod
      use domain_decomp_atm, only: grid, getdomainbounds, am_i_root
      use timestream_mod, only : init_stream,read_stream
      use pario, only : par_open,par_close,read_data
      implicit none
      integer :: jjdayo
      real*8 :: S0
      real*8, dimension(:,:,:), pointer :: o3jday
!@var delta_o3_now the difference in O3 between solar max and solar min,
!@+   interpolated to the current day
      real*8, allocatable :: delta_o3_now(:,:,:)
!@var add_sol is [S00WM2(now)-1/2(S00WM2min+S00WM2max)]/
!@+ [S00WM2max-S00WM2min] so that O3(altered) = O3(default) +
!@+ add_sol*delta_O3_now
      real*8 :: add_sol
      logical, save :: init = .false.
      integer :: i,j,l,fid,jyearx

      integer :: j_0, j_1, i_0, i_1

      jyearx = 2000 ! nominal year

      if (.not. init) then
        init = .true.

        call sync_param("use_sol_Ox_cycle",use_sol_Ox_cycle)

        if(use_sol_Ox_cycle /= 1) return

        fid = par_open(grid,'delta_O3','read')
        call read_data(grid,fid,'S0min',S0min,bcast_all=.true.)
        call read_data(grid,fid,'S0max',S0max,bcast_all=.true.)
        call par_close(grid,fid)

        call init_stream(grid,delta_O3stream,'delta_O3','O3',
     &       -1d30,1d30,'linm2m',jyearx,jjdayo,cyclic=.true.)

      endif

      if(use_sol_Ox_cycle /= 1) return

      call getdomainbounds(grid, j_strt=j_0,j_stop=j_1,
     &                           i_strt=i_0,i_stop=i_1)

      add_sol = (S0-0.5d0*(S0min+S0max))/(S0max-S0min)
      if(am_i_root()) then
        write(6,661)JJDAYO,S0,S0min,S0max,add_sol
      endif

      allocate(delta_o3_now(grid%i_strt_halo:grid%i_stop_halo,
     &                      grid%j_strt_halo:grid%j_stop_halo,nlo3))

      call read_stream(grid,delta_O3stream,jyearx,jjdayo,delta_o3_now)
      do j=j_0,j_1
      do i=i_0,i_1
        O3JDAY(:,I,J) = O3JDAY(:,I,J) + add_sol*delta_O3_now(i,j,:)
      enddo
      enddo

      deallocate(delta_o3_now)

  661 format('JJDAYO,S0,S0min,S0max,frac=',I4,3F9.2,F7.3)
      RETURN
      END SUBROUTINE UPDO3D_solar

      subroutine get_o3col(igcm,jgcm,lm_gcm,nl,plb0,fulgas3,o3col)
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
#endif
      implicit none
      integer, intent(in) :: igcm,jgcm,lm_gcm,nl
      real*8, dimension(nl+1) :: plb0
      real*8 :: fulgas3
      real*8, dimension(nl) :: o3col

!!!                   CALL GETO3D(ILON,JLAT) ! may have to be changed ??
      if(use_o3_ref > 0 )then
        CALL REPART (O3JREF(1,IGCM,JGCM),
     *          PLBO3_traditional,NLO3_traditional+1, ! in
     *                        O3col,PLB0, NL+1)  ! out, ok if L1>1 ?
        ! next block may seem weird but it is here to allow RCOMPX calls with
        ! reference ozone in part of the atmosphere and tracer below:
        FULGAS3=1.d0
      else
        CALL REPART (O3JDAY(1,IGCM,JGCM),PLBO3,NLO3+1, ! in
     *                        O3col,PLB0, NL+1)   ! out, ok if L1>1 ?
#ifdef HIGH_FREQUENCY_O3_INPUT
        ! Overwrite the lm_gcm levels with higher frequency ozone, leaving
        ! climatology above those levels:
        O3col(1:lm_gcm)=O3JDAY_HF_modelLevels(1:lm_gcm,IGCM,JGCM)
        FULGAS3=1.d0
#endif
#ifdef SCM
        if(SCMopt%ozone)then
        ! Overwrite specified SCM levels (indicated by non-zero values),
        ! leaving climatology above those levels:
          do k = 1,lm_gcm
            if(SCMin%O3(k) > 0.) O3col(k)=SCMin%O3(k)
          enddo
          FULGAS3=1.d0
        endif
#endif
      endif
      end subroutine get_o3col

      end module O3mod

      module ghgmod
!@sum GHGmod administers setting of radiatively active gas amounts.
!@+   This coding has been extracted from module RADPAR and
!@+   transformed as necessary.   The primary interface is
!@+   subroutine getgas().
!@auth A. Lacis, R. Ruedy, original development team
!@extractor M. Kelley
      use constant, only: pO2
      use constant, only: avog, mair, grav, loschmidt_constant
      use atm_com, only : lm_req
      use resolution, only : lm_gcm=>lm
      implicit none

!@var ppmv_to_cm_at_stp Conversion factor for conversion from PPMV to cm at
!                       STP. Also needs an additional factor dP for the
!                       conversion.
      REAL*8, PARAMETER :: ppmv_to_cm_at_stp = 1.0D-05*avog/
     *      (grav*mair*loschmidt_constant)

!@var h2o_mmr_to_cm_at_stp Conversion factor for conversion from mass
!                          mixing ratio to cm at STP for water vapor.
!                          Also needs an additional factor dP for the
!                          conversion.
      REAL*8, PARAMETER :: h2o_mmr_to_cm_at_stp = ppmv_to_cm_at_stp*
     *      1.0D+06*mair/18.0153D0


C     ------------------------------------------------------------------
C          NO2 Trace Gas Vertical Distribution and Concentration Profile
C     ------------------------------------------------------------------
      REAL*8, PARAMETER ::
     *     CMANO2(42)=(/            ! every 2 km starting at 0km
     1  8.66E-06,5.15E-06,2.85E-06,1.50E-06,9.89E-07,6.91E-07,7.17E-07,
     2  8.96E-07,3.67E-06,4.85E-06,5.82E-06,6.72E-06,7.77E-06,8.63E-06,
     3  8.77E-06,8.14E-06,6.91E-06,5.45E-06,4.00E-06,2.67E-06,1.60E-06,
     4  8.36E-07,3.81E-07,1.58E-07,6.35E-08,2.57E-08,1.03E-08,4.18E-09,
     5  1.66E-09,6.57E-10,2.58E-10,1.02E-10,4.11E-11,1.71E-11,7.73E-12,
     6  9.07E-12,4.63E-12,2.66E-12,1.73E-12,1.28E-12,1.02E-12,1.00E-30/)

C-----------------------------------------------------------------------
C     GHG 1980 Reference Concentrations and Vertical Profile Definitions
C-----------------------------------------------------------------------

!@var KTREND if > 0 table GHG concentrations (Trend G) are used for
!@+             yr/day KYEARG/KJDAYG; if KTREND=0, GHG are set to PPMVK0
      INTEGER :: KTREND=1

!@var PPMV80  reference GHG concentrations (ppm)
      REAL*8, dimension(13) ::
C     GAS NUMBER    1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
     *   PPMV80=(/0d0, 337.90d0, 0d0,pO2*1.d6,0d0,  .3012d0,   1.5470d0
     *     ,.1666d-03,.3003d-03, 0d0,   .978D-04,  .0010D-10,  .0420d0/)
C              CCL3F1    CCL2F2   N2     CFC-Y       CFC-Z         SO2
C     GAS NUMBER    8         9   10        11          12          13

!@var PPMVK0  user set  GHG concentrations (ppm), used if KTREND=0
      REAL*8, dimension(12) ::
C     GAS  NUMBER   1         2    3      4    5         6           7
C                 H2O       CO2   O3     O2  NO2       N2O         CH4
     *   PPMVK0=(/0d0, 337.90d0, 0d0, 21.d4, 0d0,  .3012d0,   1.5470d0
     *               ,.1666d-03,  .3003d-03, 0d0, .978D-04, 0.0010D-10/)
C                        CCL3F1      CCL2F2   N2     CFC-Y       CFC-Z
C     GAS  NUMBER             8           9   10        11          12

C     Makiko GHG Trend Compilation  GHG.1850-2050.Dec1999 in GTREND
C     ---------------------------------------------------------------
!@var nghg nr. of well-mixed GHgases: CO2 N2O CH4 CFC-11 CFC-12 others
!@var nyrsghg max.number of years of prescr. greenhouse gas history
      INTEGER, PARAMETER :: nghg=6

!@var ghgyr1,ghgyr2 first and last year of GHG history
      INTEGER ghgyr1,ghgyr2
!@var ghgam,xref,xnow     GHG-mixing ratios in ppm,ppm,ppm,ppb,ppb,ppb
      REAL*8 XREF(nghg+1),XNOW(nghg+1)
      real*8, allocatable :: ghgam(:,:)

C     GTREND:  1980.,  337.9,  .3012,  1.547,  .1666,  .3003,  .0978,
C     ---------------------------------------------------------------


!@var KPGRAD,KLATZ0 control parameters for vertical GHG profiles
!@+   -----------------------------------------------------------------
!@+   Minschwaner et al JGR (1998) CH4, N2O, CFC-12 Vertical profiles
!@+   Exponential decrease by EXP(-(Z-Z0)/H) is superimposed on this.
!@+   IF(KLATZ0 > 0) Then: Z0 depends on latitude, KGGVDF not used
!@+   KPGRAD>0: Pole-to-Pole lat. gradient (PPGRAD) is also superimposed
!@+   ------------------------------------------------------------------
!@var PPGRAD  Pole-to-Pole latitud.gradient for GHG (used if KPGRAD > 0)
      INTEGER :: KPGRAD=1, KLATZ0=1

!@var Z0,ZH   scale heights used for vertical profile (km)
      REAL*8, dimension(12) ::
C     NUMBER   1    2    3    4  5    6    7    8     9   10   11  12
C             H2O  CO2  O3   O2 NO2  N2O  CH4 CFC11 CFC12 N2 CF-Y  CF-Z
     *   Z0=(/0.0, 0.0,0.0, 0.0,0.0, 16., 16., 16., 16., 0.0, 16., 16./)
     *  ,ZH=(/8.0, 8.0,8.0, 8.0,8.0, 30., 50., 30., 30., 0.0, 30., 30./)


C     GAS  NUMBER   1     2    3    4    5         6         7
C                 H2O   CO2   O3   O2  NO2       N2O       CH4
     *  ,PPGRAD=(/0.0,  0.0, 0.0, 0.0, 0.0,   0.0100,   0.0900,
     *               0.0600,   0.0600, 0.0,   0.0600,   0.0600/)
C                    CCL3F1    CCL2F2   N2     CFC-Y     CFC-Z
C     GAS  NUMBER         8         9   10        11        12

!@var FULGAS scales the various atmospheric constituents:
!@+         H2O CO2 O3 O2 NO2 N2O CH4 F11 F12 N2C CFC11 CFC12 SO2
!@+   Note: FULGAS(1) only acts in the stratosphere (unless LS1_loc=1)
      REAL*8 :: FULGAS(13) = (/    ! scales ULGAS

C      H2O CO2  O3  O2 NO2 N2O CH4 F11 F12 N2C CFC11+ CFC12+ SO2
C        1   2   3   4   5   6   7   8   9  10    11     12   13
     +   1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,   1.,    1.,  0./)
#ifdef ALTER_RADF_BY_LAT
!@var FULGAS_orig saves initial FULGAS values
      REAL*8, dimension(13) :: FULGAS_orig
#endif

!@var chem_IN column variable for importing ozone(1) and methane(2)
!@+   fields from rest of model
!@var use_tracer_chem:set U0GAS(L, )=chem_IN( ,L), L=L1,use_tracer_chem( )
      REAL*8 :: chem_IN(2,lm_gcm+lm_req) ! would 1:lm be ok?
      INTEGER :: use_tracer_chem(2)

      integer, parameter :: lxghg = lm_gcm + lm_req
      real*8, dimension(lxghg+1) :: plb0ghg,hlb0ghg

      real*8, dimension(lxghg,13), private :: u0gas

      contains

      subroutine init_ghgmod(plb_in,hlb_in)
      real*8, dimension(lxghg+1) :: plb_in,hlb_in
      plb0ghg = plb_in
      hlb0ghg = hlb_in
      end subroutine init_ghgmod

      SUBROUTINE SETGHG(JYEARG,JJDAYG)
      IMPLICIT NONE
C
C
C     ---------------------------------------------------------------
C     SETGHG  Sets Default Greenhouse Gas Reference Year (for FULGAS)
C
C     Control Parameter:
C                     KTREND (specified in RADPAR) activates GH Trend
C               Default
C     KTREND  =   1
C     Selects   GTREND
C     ---------------------------------------------------------------
      INTEGER, INTENT(IN) :: JYEARG,JJDAYG
      REAL*8 TREF,TNOW
      INTEGER I
C
      TREF=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND==0) THEN
        XREF(1)=PPMV80(2)
        XREF(2)=PPMV80(6)
        XREF(3)=PPMV80(7)
        XREF(4)=PPMV80(8)*1000.D0
        XREF(5)=PPMV80(9)*1000.D0
        XREF(6)=PPMV80(11)*1000.D0  ! YREF11=PPMV80(11)*1000.D0
        XREF(7)=PPMV80(12)*1000.D0  ! ZREF12=PPMV80(12)*1000.D0
        RETURN
      END IF

      CALL GTREND(XREF,TREF)     ! finds xref 1-6 (yref11=xx6=xref(6))
      XREF(7)=1.D-13             ! ZREF12=1.D-13
      DO 120 I=1,NGHG
      IF(XREF(I) < 1.D-06) XREF(I)=1.D-06
  120 CONTINUE
      PPMV80(2)=XREF(1)
      PPMV80(6)=XREF(2)
      PPMV80(7)=XREF(3)
      PPMV80(8)=XREF(4)/1000.D0
      PPMV80(9)=XREF(5)/1000.D0
      PPMV80(11)=XREF(6)/1000.d0   ! YREF11/1000.D0
      PPMV80(12)=XREF(7)/1000.D0   ! ZREF12/1000.D0
      RETURN
      end SUBROUTINE SETGHG
C
C--------------------------------
!      ENTRY UPDGHG(JYEARG,JJDAYG)
C--------------------------------
      subroutine UPDGHG(JYEARG,JJDAYG)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JYEARG,JJDAYG
      REAL*8 TREF,TNOW
      INTEGER I
C
      TNOW=JYEARG+(JJDAYG-0.999D0)/366.D0
C
      IF(KTREND==0) THEN
        FULGAS(2)=PPMVK0(2)/XREF(1)
        FULGAS(6)=PPMVK0(6)/XREF(2)
        FULGAS(7)=PPMVK0(7)/XREF(3)
        FULGAS(8)=PPMVK0(8)/XREF(4)
        FULGAS(9)=PPMVK0(9)/XREF(5)
        FULGAS(11)=PPMVK0(11)/XREF(6) ! YREF11
        FULGAS(12)=PPMVK0(12)/XREF(7) ! .../ZREF12
        RETURN
      END IF

      CALL GTREND(XNOW,TNOW) ! finds xnow 1-6 (ynow11=xx6=xnow(6))
      XNOW(7)=1.D-20         ! ZNOW12=1.D-20
      FULGAS(2)=XNOW(1)/XREF(1)
      FULGAS(6)=XNOW(2)/XREF(2)
      FULGAS(7)=XNOW(3)/XREF(3)
      FULGAS(8)=XNOW(4)/XREF(4)
      FULGAS(9)=XNOW(5)/XREF(5)
      FULGAS(11)=XNOW(6)/XREF(6) ! YNOW11/YREF11
      FULGAS(12)=XNOW(7)/XREF(7) ! ZNOW12/ZREF12
C
      RETURN
      end subroutine UPDGHG

      subroutine getgas(
     &     igcm,jgcm  ! gcm-grid i,j indices (currently only for retrieving O3)
     &     ,jlat46    ! lat index on reference 72x46 lon-lat grid
     &     ,plb       ! layer edge pressures (mb)
     &     ,ulgas     ! per-layer gas amounts (cm-atm)
     &     ,u0gas_out ! optional: ulgas without multiplication by fulgas
     &     )
! fills elements 2-13 of ulgas
! Only side effect on module variables: setting of fulgas(3); not sure how to manage that
!      use O3mod, only : plbo3,nlo3,plbo3_traditional,nlo3_traditional
!      use O3mod, only : use_o3_ref,o3jref,o3jday
      use O3mod, only : get_o3col
      IMPLICIT NONE
      integer, intent(in) :: igcm,jgcm,jlat46
      real*8, dimension(lxghg+1), intent(in) :: plb
      real*8, dimension(lxghg,13) :: ulgas
      real*8, dimension(lxghg,13), optional :: u0gas_out
!
      integer, parameter :: l1=1
      REAL*8, PARAMETER :: PI=3.141592653589793D0  ! yapi: yet another pi
      INTEGER I,L,J,K,N,NL0,NL
      REAL*8 DP
     &     ,ZT,ZB,EXPZT,EXPZB,PARTTG
     &     ,PTRO,DL,DLS,DLN
     &     ,Z0LAT
     &     ,ACM,HI,FI,HL,HJ,FJ,DH,FF,SINLAT
      real*8, dimension(lxghg) :: pl,dpl,o3col
      real*8, dimension(lxghg+1) :: plb0,hlb0
      real*8 :: lat46 ! center latitude on 72x46 grid

      nl = lxghg
      nl0 = nl

      plb0 = plb0ghg
      hlb0 = hlb0ghg

      lat46 = -94d0 + 4d0*real(jlat46,kind=8)

      do l=1,nl
        dpl(l)=plb(l)-plb(l+1)
        pl(l)=(plb(l)+plb(l+1))*0.5d0
      enddo



C --------------------------------------------------------------------
C Uniformly Mixed Gas Distributions (can be overwritten/modified below)
C --------------------------------------------------------------------
      DO K=2,13
        U0GAS(1:NL0,K)=PPMV80(K)*ppmv_to_cm_at_stp*DPL(1:NL0)
      ENDDO

C                -----------------------------------------------------
C                N20,CH4,F11,F12 Specified Latitudinal Z0 Distribution
C                -----------------------------------------------------

      IF(KLATZ0 > 0) THEN
        PTRO=100.D0
        DL=LAT46
        DLS=-40.D0
        DLN= 40.D0
        IF(DL < DLS) PTRO=189.D0-(DL+40.D0)*2.22D0
        IF(DL > DLN) PTRO=189.D0+(DL-40.D0)*2.22D0
        DO L=1,NL0
          IF(PLB0(L) >= PTRO) Z0LAT=HLB0(L)  ! orig. hlb not hlb0
        END DO
        DO K=6,12
          IF(K==10) CYCLE
          DO L=1,NL0
            IF(PLB0(1) >= PTRO) THEN ! safety check until P,H hard-coding removed
              ZT=(HLB0(L+1)-Z0LAT)/ZH(K) ! orig. hlb not hlb0
              IF(ZT <= 0.D0) CYCLE
              ZB=(HLB0(L)-Z0LAT)/ZH(K) ! orig. hlb not hlb0
              EXPZT=EXP(-ZT)
              EXPZB=EXP(-ZB)
              IF(ZB < 0.D0) EXPZB=1.D0-ZB
              U0GAS(L,K)=U0GAS(L,K)*(EXPZB-EXPZT)/max(ZT-ZB,1d-6)
            ENDIF               ! safety check
          ENDDO
        ENDDO
      ENDIF


C                                         ----------------------------
C     IGAS=5                              Global Mean NO2 Distribution
C                                         ----------------------------
      u0gas(:,5)=0.
      ACM=0.D0
      HI=0.D0
      HJ=HI+2.D0
      FI=CMANO2(1)
      FJ=CMANO2(2)
      HL=HLB0(2)
      L=1
      J=2
      DO
        DH=HJ-HI
        IF(HJ <= HL) THEN
          ACM=ACM+(FI+FJ)*DH*0.5D0
          HI=HJ
          FI=FJ
          J=J+1
          IF(J > 42) EXIT
          HJ=HI+2.D0
          FJ=CMANO2(J)
        ELSE
          FF=FI+(FJ-FI)*(HL-HI)/DH
          DH=HL-HI
          ACM=ACM+(FI+FJ)*DH*0.5D0
          U0GAS(L,5)=ACM
          ACM=0.D0
          HI=HL
          FI=FF
          IF(L==NL0) EXIT
          L=L+1
          HL=HLB0(L+1)
        ENDIF
      ENDDO
      u0gas(l,5) = acm
      ! layer thickness amount scaling
      do l=l1,nl0
        u0gas(l,5) = u0gas(l,5) * (plb(l)-plb(l+1))/(plb0(l)-plb0(l+1))
      enddo

C****
C**** O3 section
C****
      call get_o3col(igcm,jgcm,lm_gcm,nl,plb0,fulgas(3),o3col)
      u0gas(1:nl,3) = o3col(1:nl)
      if(use_tracer_chem(1) > 0) then
        FULGAS(3)=1.d0
      endif

C****
C**** Scaling section
C****

      IF(KPGRAD > 0) THEN ! "pole-to-pole gradient"
        SINLAT = SIN(LAT46*PI/180.D0)
        DO K=2,12
          IF(K==3) CYCLE
          PARTTG=(1.D0+0.5D0*PPGRAD(K)*SINLAT)
          DO L=L1,NL0           ! =L1,NL for GCM use, =1,NL0 for offline use
            U0GAS(L,K)=U0GAS(L,K)*PARTTG
          ENDDO
        ENDDO
      ENDIF

      ! fulgas scaling
      do k=2,13
        ulgas(1:nl0,k)=u0gas(1:nl0,k)*fulgas(k)
      enddo

      if(present(u0gas_out)) then
        u0gas_out(1:nl0,2:13)=u0gas(1:nl0,2:13)
      endif

      RETURN
      END SUBROUTINE GETGAS

      end module ghgmod

      SUBROUTINE GTREND(XNOW,TNOW)
C
      use ghgmod, only: nghg,ghgyr1,ghgyr2,ghgam
      IMPLICIT NONE
      REAL*8 xnow(nghg),tnow,year,dy,frac
      INTEGER iy,n
C
C-------------------------------------------------------------
C        Makiko GHG Trend Compilation  GHG.1850-2050.Dec1999
C
C        Annual-Mean      Greenhouse Gas Mixing Ratios
C-------------------------------------------------------------
C                 CO2     N2O     CH4   CFC-11  CFC-12  others
C        Year     ppm     ppm     ppm     ppb     ppb     ppb
C-------------------------------------------------------------
C     Read from external file - outside table: use value from
C                                      years ghgyr1 or ghgyr2
      YEAR=TNOW
      IF(TNOW <= ghgyr1+.5D0) YEAR=ghgyr1+.5D0
      IF(TNOW >= ghgyr2+.49999D0) YEAR=ghgyr2+.49999D0
      DY=YEAR-(ghgyr1+.5D0)
      IY=DY
      frac=DY-IY
      IY=IY+1
C
C     CO2 N2O CH4 CFC-11 CFC-12 other_GHG  SCENARIO
C--------------------------------------------------
C
      do n=1,nghg
        XNOW(N)=GHGAM(N,IY)+frac*(GHGAM(N,IY+1)-GHGAM(N,IY))
      end do
C
      RETURN
      END SUBROUTINE GTREND

      SUBROUTINE GHGHST(iu,ghg_yr)
!@sum  reads history for nghg well-mixed greenhouse gases
!@auth R. Ruedy

      use domain_decomp_atm, only : write_parallel
      USE GHGMOD, only : nghg,ghgyr1,ghgyr2,ghgam
      IMPLICIT NONE
      INTEGER :: iu,ghg_yr
      integer :: n,k,nhead=4,iyr
      CHARACTER*80 title
      character(len=300) :: out_line

      write(out_line,*)  ! print header lines and first data line
      call write_parallel(trim(out_line),unit=6)
      do n=1,nhead+1
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
      end do
      if(title(1:2).eq.'--') then                 ! older format
        read(iu,'(a)') title
        write(out_line,'(1x,a80)') title
        call write_parallel(trim(out_line),unit=6)
        nhead=5
      end if

!**** find range of table: ghgyr1 - ghgyr2
      read(title,*) ghgyr1
      do ; read(iu,'(a)',end=20) title ; end do
   20 read(title,*) ghgyr2
      rewind iu  !   position to data lines
      do n=1,nhead ; read(iu,'(a)') ; end do

      allocate (ghgam(nghg,ghgyr2-ghgyr1+1))
      do n=1,ghgyr2-ghgyr1+1
        read(iu,*) iyr,(ghgam(k,n),k=1,nghg)
        do k=1,nghg ! replace -999. by reasonable numbers
          if(ghgam(k,n).lt.0.) ghgam(k,n)=ghgam(k,n-1)
        end do
        if(ghg_yr>0 .and. abs(ghg_yr-iyr).le.1) then
          write(out_line,'(i5,6f10.4)') iyr,(ghgam(k,n),k=1,nghg)
          call write_parallel(trim(out_line),unit=6)
        endif
      end do
      write(out_line,*) 'read GHG table for years',ghgyr1,' - ',ghgyr2
      call write_parallel(trim(out_line),unit=6)
      return
      end SUBROUTINE GHGHST
