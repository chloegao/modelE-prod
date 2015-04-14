#include "rundeck_opts.h"
      module trdust_drv
!@sum  trdust_drv routines with resolution dependent variables for
!@+               soil dust aerosols
!@auth Jan Perlwitz
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

      use filemanager,only: nameunit,openunit,closeunit
      use constant, only: rgas
      use resolution, only: im,jm,lm
      use Dictionary_mod, only : sync_param
      use domain_decomp_atm, only: am_i_root,grid,dread_parallel
     &     ,broadcast,write_parallel,getDomainBounds
      use model_com, only: ioread,iowrite,irsfic,irsficno
     &     ,irerun,itime
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR,INT_MONTHS_PER_YEAR
      use fluxes, only: dust_flux_glob
#if (defined TRACERS_DUST) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     &     ,dust_flux2_glob
#endif
#ifdef TRACERS_DRYDEP
     &     ,depo_turb_glob,depo_grav_glob
#endif
#ifdef TRACERS_WATER
     &     ,trprec
#else
     &     ,trprec_dust
#endif
      use OldTracer_mod, only: trname
      use tracer_com,only: n_clay,n_clayilli,n_sil1quhe,ntm_dust,trm
     &     ,coupled_chem
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: dodrydep
#endif
#ifdef TRACERS_WATER
      use OldTracer_mod, only: dowetdep
#endif
      use trdiag_com, only: trcsurf,trcSurfByVol,to_conc,set_to_conc
      use tracers_dust
#ifdef NEW_IO
      use pario, only: defvar,read_dist_data,write_dist_data
#endif
      use PolynomialInterpolator_mod, only: interpolator3D

      implicit none

      include 'netcdf.inc'

      integer :: i_0,i_1,j_0,j_1
      type(interpolator3D), save :: wsgInterp

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP || TRACERS_TOMAS*/

      contains

c init_soildust
      subroutine init_soildust
!@sum init_soildust  initializiations for soil dust/mineral dust aerosols
!@+    at startup
!@auth Jan Perlwitz
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

      implicit none

      logical, save :: qfirst = .true.

      integer :: n, n1

      if ( .not. qfirst ) return
      qfirst = .false.

#ifdef TRACERS_DUST
      n_soilDust = n_clay
#else
#ifdef TRACERS_MINERALS
      n_soilDust = n_clayilli
#else
#ifdef TRACERS_QUARZHEM
      n_soilDust = n_sil1quhe
#endif
#endif
#endif

#ifndef TRACERS_AMP
#ifndef TRACERS_TOMAS
c**** initialize dust names
      do n=1,Ntm_dust
        n1=n_soilDust+n-1
        dust_names(n) = trname(n1)
      end do

c**** insert to_conc_soildust into to_conc

      if ( .not. 
     &     any([(to_conc(n),n=n_soilDust,n_soilDust+ntm_dust-1)] > 0 )
     &     ) then
        call sync_param( 'to_conc_soildust', to_conc_soildust )
        do n = n_soilDust, n_soilDust+ntm_dust-1
          call set_to_conc(n, to_conc_soildust)
        end do
      end if
#endif
#endif

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
      call sync_param( 'brittleFactor', brittleFactor, Mtrac )
      call sync_param( 'calcMineralAggrProb', calcMineralAggrProb )
#endif

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP || TRACERS_TOMAS*/
      return
      end subroutine init_soildust

c tracer_ic_soildust
      subroutine tracer_ic_soildust
!@sum tracer_ic_soildust  reads in source and parameter files for
!@+    dust/mineral tracer at itime_tr0
!@auth Jan Perlwitz
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)

      IMPLICIT NONE

      integer :: i,ierr,j,io_data,k,ires
      INTEGER startd(3),countd(3),statusd
      INTEGER idd1,idd2,idd3,idd4,ncidd1,ncidd2,ncidd3,ncidd4
      REAL*8 :: zsum
c**** temporary array to read in data
      REAL*4,DIMENSION(grid%i_strt:grid%i_stop,
     &                 grid%j_strt:grid%j_stop) :: work ! no halo

      LOGICAL,SAVE :: qfirst=.TRUE.
      CHARACTER :: cierr*3,name*256
      CHARACTER(80) :: OptModelVers='No Entry'

      IF (.NOT. qfirst) RETURN
      qfirst=.FALSE.

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1, 
     &     I_STRT=I_0, I_STOP=I_1)

c**** prescribed AEROCOM dust emissions
      IF (imDust == 1) THEN

c these netcdf reads are still latlon-specific.
c will call read_dist_data for cubed sphere compatibility

        statusd=NF_OPEN('dust_bin1',NCNOWRIT,ncidd1)
        statusd=NF_OPEN('dust_bin2',NCNOWRIT,ncidd2)
        statusd=NF_OPEN('dust_bin3',NCNOWRIT,ncidd3)
        statusd=NF_OPEN('dust_bin4',NCNOWRIT,ncidd4)

        statusd=NF_INQ_VARID(ncidd1,'dust',idd1)
        statusd=NF_INQ_VARID(ncidd2,'dust',idd2)
        statusd=NF_INQ_VARID(ncidd3,'dust',idd3)
        statusd=NF_INQ_VARID(ncidd4,'dust',idd4)

        startd(1)=i_0
        startd(2)=j_0

        countd(1)=1+(i_1-i_0)
        countd(2)=1+(j_1-j_0)
        countd(3)=1

        DO k=1,INT_DAYS_PER_YEAR

          IF (k > 59) THEN
            startd(3)=k+1
          ELSE
            startd(3)=k
          END IF

          statusd=NF_GET_VARA_REAL(ncidd1,idd1,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,1,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd2,idd2,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,2,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd3,idd3,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,3,k)=DBLE(work(i_0:i_1,j_0:j_1))
          statusd=NF_GET_VARA_REAL(ncidd4,idd4,startd,countd,work)
          d_dust(i_0:i_1,j_0:j_1,4,k)=DBLE(work(i_0:i_1,j_0:j_1))

        END DO

        statusd=NF_CLOSE(ncidd1)
        statusd=NF_CLOSE(ncidd2)
        statusd=NF_CLOSE(ncidd3)
        statusd=NF_CLOSE(ncidd4)

        CALL write_parallel(' Read from file dust_bin[1-4]',unit=6)

      ELSE IF (imDust == 2) THEN
c**** legacy emission scheme
        IF (Im /= 72 .OR. Jm /= 46) CALL stop_model
     & ('Stopped in init_dust: imDust=2 works only for Im=72 and Jm=46'
     &       ,255)

c**** Read input: threshold speed
        CALL openunit('VTRSH',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),vtrsh)
        CALL closeunit(io_data)

c**** Read input: fraction clay
        CALL openunit('FRCLAY',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),frclay)
        CALL closeunit(io_data)

c**** Read input: fraction silt
        CALL openunit('FRSILT',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),frsilt)
        CALL closeunit(io_data)

c**** Read input: prec-evap data
        CALL openunit('DRYHR',io_data,.TRUE.,.TRUE.)
        CALL dread_parallel(grid,io_data,nameunit(io_data),dryhr)
        CALL closeunit(io_data)

      ELSE IF (imDust == 0) THEN
c**** Probability density function scheme for dust emission

c**** Read input: ERS data
        CALL openunit('ERS',io_data,.TRUE.,.TRUE.)
        DO k=1,INT_MONTHS_PER_YEAR
          CALL dread_parallel(grid,io_data,nameunit(io_data),
     &         ers_data(:,:,k))
        END DO
        CALL closeunit(io_data)

c**** Read input: source function data
        call openunit('DSRC',io_data,.true.,.true.)
        call dread_parallel(grid,io_data,nameunit(io_data)
     &       ,dustSourceFunction)
        CALL closeunit(io_data)

        ires = 1 ! default
c**** set parameters depending on the preferred sources chosen
        ires = 0
        select case(prefDustSources)
        case(0)                 ! Ginoux 2001 sources w/ vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.036636533D0
              fracSiltPDFscheme = 0.023985135D0
#else
              fracClayPDFscheme = 0.091393298D0
              fracSiltPDFscheme = 0.10973922D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.039218312D0
              fracSiltPDFscheme = 0.021643261D0
#else
              fracClayPDFscheme = 0.10309873D0
              fracSiltPDFscheme = 0.080987688D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.091393298D0  ! not optimized yet
              fracSiltPDFscheme = 0.10973922D0   ! not optimized yet
            else
              fracClayPDFscheme = 0.10309873D0   ! not optimized yet
              fracSiltPDFscheme = 0.080987688D0  ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.091393298D0  ! not optimized yet
              fracSiltPDFscheme = 0.10973922D0   ! not optimized yet
            else
              fracClayPDFscheme = 0.10309873D0   ! not optimized yet
              fracSiltPDFscheme = 0.080987688D0  ! not optimized yet
            end if
            ires=4
          case(90)              ! same as for prefDustSources = 1
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1' //
     &             ' and coupled chemistry = 0'
            else
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1'
            end if
            ires=5
          end select
        case(1)                 ! Ginoux 2009 sources w/ vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.034805282D0
              fracSiltPDFscheme = 0.0307729D0
#else
              fracClayPDFscheme = 0.086874723D0
              fracSiltPDFscheme = 0.097082074D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.053818764d0
              fracSiltPDFscheme = 0.063865154d0
              OptModelVers='AR 5 branch, 10/18/2011, 10:00 AM EDT'
#else
              fracClayPDFscheme = 0.091387274D0
              fracSiltPDFscheme = 0.15582714D0
              OptModelVers='02/20/2010, 23:30 EST'
#endif
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.086874723D0 ! not optimized yet
              fracSiltPDFscheme = 0.097082074D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.091387274D0 ! not optimized yet
              fracSiltPDFscheme = 0.15582714D0  ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.086874723D0 ! not optimized yet
              fracSiltPDFscheme = 0.097082074D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.091387274D0 ! not optimized yet
              fracSiltPDFscheme = 0.15582714D0  ! not optimized yet
            end if
            ires=4
          case(90)
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.059819419D0 ! same as for dust only case
              fracSiltPDFscheme = 0.13299809D0  ! same as for dust only case
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for coupled chemistry = 0'
            else
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT'
            end if
            ires=5
          end select
        case(2)                 ! Ginoux 2009 sources w/o vegetation mask
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.028216341D0
              fracSiltPDFscheme = 0.010625339D0
#else
              fracClayPDFscheme = 0.051027254D0
              fracSiltPDFscheme = 0.047415049D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.03028031D0
              fracSiltPDFscheme = 0.097915021D-1
#else
              fracClayPDFscheme = 0.059176377D0
              fracSiltPDFscheme = 0.056770932D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.051027254D0 ! not optimized yet
              fracSiltPDFscheme = 0.047415049D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.059176377D0 ! not optimized yet
              fracSiltPDFscheme = 0.056770932D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.051027254D0 ! not optimized yet
              fracSiltPDFscheme = 0.047415049D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.059176377D0 ! not optimized yet
              fracSiltPDFscheme = 0.056770932D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! same as for prefDustSources = 1
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1' //
     &             ' and coupled chemistry = 0'
            else
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1'
            end if
            ires=5
          end select
        case(3)                 ! Grini/Zender sources
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.015809495D0
              fracSiltPDFscheme = 0.059532525D-1
#else
              fracClayPDFscheme = 0.03586315D0
              fracSiltPDFscheme = 0.028140008D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.016681537D0
              fracSiltPDFscheme = 0.054484947D-1
#else
              fracClayPDFscheme = 0.036893354D0
              fracSiltPDFscheme = 0.039649002D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.03586315D0  ! not optimized yet
              fracSiltPDFscheme = 0.028140008D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.036893354D0 ! not optimized yet
              fracSiltPDFscheme = 0.039649002D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.03586315D0  ! not optimized yet
              fracSiltPDFscheme = 0.028140008D0 ! not optimized yet
            else
              fracClayPDFscheme = 0.036893354D0 ! not optimized yet
              fracSiltPDFscheme = 0.039649002D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! same as for prefDustSources = 1
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1' //
     &             ' and coupled chemistry = 0'
            else
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1'
            end if
            ires=5
          end select
        case(4)                 ! Tegen sources
          select case(im)
          case(72)              ! uses old values for Ginoux 2001 source file
            fracClayPDFscheme = 0.092335D0 ! not optimized yet
            fracSiltPDFscheme = 0.226916D0 ! not optimized yet
            ires=1
          case(144)
            if (coupled_chem == 1) then
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.046381759D0
              fracSiltPDFscheme = 0.040315694D-1
#else
              fracClayPDFscheme = 0.083940461D0
              fracSiltPDFscheme = 0.06811547D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            else
#ifdef NUDGE_ON
              fracClayPDFscheme = 0.048860416D0
              fracSiltPDFscheme = 0.029009232D-1
#else
              fracClayPDFscheme = 0.09620785D0
              fracSiltPDFscheme = 0.058139336D0
#endif
              OptModelVers='02/20/2010, 23:30 EST'
            end if
            ires=2
          case(288)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.083940461D0 ! not optimized yet
              fracSiltPDFscheme = 0.06811547D0  ! not optimized yet
            else
              fracClayPDFscheme = 0.09620785D0  ! not optimized yet
              fracSiltPDFscheme = 0.058139336D0 ! not optimized yet
            end if
            ires=3
          case(360)             ! uses values for im=144 for now
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.083940461D0 ! not optimized yet
              fracSiltPDFscheme = 0.06811547D0  ! not optimized yet
            else
              fracClayPDFscheme = 0.09620785D0  ! not optimized yet
              fracSiltPDFscheme = 0.058139336D0 ! not optimized yet
            end if
            ires=4
          case(90)              ! same as for prefDustSources = 1
            if (coupled_chem == 1) then
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1' //
     &             ' and coupled chemistry = 0'
            else
              fracClayPDFscheme = 0.059819419D0
              fracSiltPDFscheme = 0.13299809D0
              OptModelVers = 'AR5 branch, 07/04/2011, 11:55 PM EDT' //
     &             ', same as for prefDustSources = 1'
            end if
            ires=5
          end select
        end select

        if (am_i_root()) then
          write(6,*) ' Actually used parameters for soil dust emission:'
          write(6,'(1x,a28,f12.9)') '  Clay: fracClayPDFscheme = '
     &         ,fracClayPDFscheme
          write(6,'(1x,a28,f12.9)') '  Silt: fracSiltPDFscheme = '
     &         ,fracSiltPDFscheme
          if (prefDustSources <= numDustSourceOpt-1 .and. ires > 0) then
            write(6,*) '  For prefDustSources =',prefDustSources,','
            write(6,*) '  these parameters are the optimized values for'
            write(6,*) '  following file with preferred dust sources:'
            write(6,*) '  >> '
     &           ,trim(dustSourceFile(prefDustSources,ires)),' <<'
            write(6,*) '  for model version: ',trim(OptModelVers)
            if (coupled_chem == 1) then
              write(6,*) '  for all aerosols and chemistry'
            else
              write(6,*) '  for dust aerosols w/o other aerosols and'
     &             ,' chemistry'
            end if
            write(6,*) '  For a free choice of emission parameters set'
            write(6,*) '  prefDustSources >',numDustSourceOpt-1
            write(6,*) '  and set fracClayPDFscheme and'
     &           ,' fracSiltPDFscheme in rundeck'
          end if
        end if

c**** Read input: EMISSION LOOKUP TABLE data
        IF (am_i_root()) THEN
          CALL openunit('LKTAB',io_data,.TRUE.,.TRUE.)
          DO k=1,Lkm
            READ(io_data,IOSTAT=ierr) ((table(i,j,k),i=1,Lim),j=1,Ljm)
          END DO
          name=nameunit(io_data)
          CALL closeunit(io_data)
          IF (ierr == 0) THEN
            write(6,*) ' Read from file '//TRIM(name)
          ELSE
            WRITE(cierr,'(I2)') ierr
            write(6,*) ' READ ERROR ON FILE '//TRIM(name)//' rc='//cierr
          END IF
        END IF
        CALL broadcast(grid,ierr)
        if(ierr.ne.0) CALL stop_model('init_dust: READ ERROR',255)
        CALL broadcast(grid,table)

c**** index of table for threshold velocity from 6.5 to 17 m/s
        DO k=1,Lkm
          x3(k)=6.d0+0.5d0*k
        END DO

c**** index of table for sub grid scale velocity (sigma) from .0001 to 30 m/s
        zsum=0.d0
        DO j=1,Ljm
          IF (j <= 30) THEN
            zsum=zsum+0.0001d0+FLOAT(j-1)*0.00008d0
            x2(j)=zsum
          ELSE IF (j > 30) THEN
            zsum=zsum-0.055254d0+0.005471d0*FLOAT(j)-
     &           1.938365d-4*FLOAT(j)**2.d0+
     &           3.109634d-6*FLOAT(j)**3.d0-
     &           2.126684d-8*FLOAT(j)**4.d0+
     &           5.128648d-11*FLOAT(j)**5.d0
            x2(j)=zsum
          END IF
        END DO
c**** index of table for GCM surface wind speed from 0.0001 to 30 m/s
        x1(:)=x2(:)
      ELSE
        CALL stop_model
     &     ('Stopped in init_dust: parameter imDUST must be <= 2',255)
      END IF
      wsgInterp = interpolator3D(x1,x2,x3,table)

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
!     read mineral fractions from input file
      call read_minfr_claquin1999_netcdf

!     apply mineral dependent brittle factor to mineral fractions
      do m = 1,mtrac
        minfr( i_0:i_1, j_0:j_1, m ) = brittleFactor( m ) * minfr(
     &       i_0:i_1, j_0:j_1, m )
      end do
      call fraction_of_mineral_claquin1999
#endif /* TRACERS_MINERAL || TRACERS_QUARZHEM */

#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM || TRACERS_AMP || TRACERS_TOMAS*/
      RETURN

      end subroutine tracer_ic_soildust

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c io_trDust
      subroutine io_trDust(kunit,iaction)
!@sum  io_trDust I/O of specific dust aerosol diagnostics (not netcdf)
!@auth Jan Perlwitz
      use domain_decomp_1d, only: pack_data,unpack_data
      implicit none

!@var kunit unit number of read/write
!@var iaction flag for reading or writing to file
      integer,intent(in) :: iaction,kunit

!@var iostat I/O status
      integer :: iostat
!@var header text with information on variables
      character(len=80) :: header

!@var dustDiagSubdd_glob global structured variable of dust aerosol diagnostics
      type(dustDiagSubdd) :: dustDiagSubdd_glob

      allocate(dustDiagSubdd_glob%dustEmission(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustEmission2(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustDepoTurb(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustDepoGrav(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustMassInPrec(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustSurfMixR(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustSurfConc(im,jm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustMass(im,jm,lm,Ntm_dust))
      allocate(dustDiagSubdd_glob%dustConc(im,jm,lm,Ntm_dust))

      select case(iaction)

      case(:iowrite)            ! write to restart file
        call pack_data(grid,dustDiagSubdd_acc%dustEmission
     &       ,dustDiagSubdd_glob%dustEmission)
        call pack_data(grid,dustDiagSubdd_acc%dustEmission2
     &       ,dustDiagSubdd_glob%dustEmission2)
        call pack_data(grid,dustDiagSubdd_acc%dustDepoTurb
     &       ,dustDiagSubdd_glob%dustDepoTurb)
        call pack_data(grid,dustDiagSubdd_acc%dustDepoGrav
     &       ,dustDiagSubdd_glob%dustDepoGrav)
        call pack_data(grid,dustDiagSubdd_acc%dustMassInPrec
     &       ,dustDiagSubdd_glob%dustMassInPrec)
        call pack_data(grid,dustDiagSubdd_acc%dustSurfMixR
     &       ,dustDiagSubdd_glob%dustSurfMixR)
        call pack_data(grid,dustDiagSubdd_acc%dustSurfConc
     &       ,dustDiagSubdd_glob%dustSurfConc)
        call pack_data(grid,dustDiagSubdd_acc%dustMass
     &       ,dustDiagSubdd_glob%dustMass)
        call pack_data(grid,dustDiagSubdd_acc%dustConc
     &       ,dustDiagSubdd_glob%dustConc)
        header='For subdaily dust tracers diagnostics: dustDiagSubdd'
        if (am_i_root()) then
          write(kunit,iostat=iostat) header
     &         ,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_glob%dustConc
          if (iostat > 0) call stop_model
     &         ('In io_trdust_drv: Restart file write error',255)
        end if

      case(ioread:)
        select case(iaction)    ! read from restart file
        case(ioread,irerun,irsfic,irsficno) ! restarts
          if (am_i_root()) then
            read(kunit,iostat=iostat) header
     &         ,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_glob%dustConc
            if (iostat > 0) call stop_model
     &           ('In io trdust_drv: Restart file read error',255)
          end if
          call unpack_data(grid,dustDiagSubdd_glob%dustEmission
     &         ,dustDiagSubdd_acc%dustEmission)
          call unpack_data(grid,dustDiagSubdd_glob%dustEmission2
     &         ,dustDiagSubdd_acc%dustEmission2)
          call unpack_data(grid,dustDiagSubdd_glob%dustDepoTurb
     &         ,dustDiagSubdd_acc%dustDepoTurb)
          call unpack_data(grid,dustDiagSubdd_glob%dustDepoGrav
     &         ,dustDiagSubdd_acc%dustDepoGrav)
          call unpack_data(grid,dustDiagSubdd_glob%dustMassInPrec
     &         ,dustDiagSubdd_acc%dustMassInPrec)
          call unpack_data(grid,dustDiagSubdd_glob%dustSurfMixR
     &         ,dustDiagSubdd_acc%dustSurfMixR)
          call unpack_data(grid,dustDiagSubdd_glob%dustSurfConc
     &         ,dustDiagSubdd_acc%dustSurfConc)
          call unpack_data(grid,dustDiagSubdd_glob%dustMass
     &         ,dustDiagSubdd_acc%dustMass)
          call unpack_data(grid,dustDiagSubdd_glob%dustConc
     &         ,dustDiagSubdd_acc%dustConc)
        end select

      end select

      deallocate(dustDiagSubdd_glob%dustEmission)
      deallocate(dustDiagSubdd_glob%dustEmission2)
      deallocate(dustDiagSubdd_glob%dustDepoTurb)
      deallocate(dustDiagSubdd_glob%dustDepoGrav)
      deallocate(dustDiagSubdd_glob%dustMassInPrec)
      deallocate(dustDiagSubdd_glob%dustSurfMixR)
      deallocate(dustDiagSubdd_glob%dustSurfConc)
      deallocate(dustDiagSubdd_glob%dustMass)
      deallocate(dustDiagSubdd_glob%dustConc)

      return
      end subroutine io_trDust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/

#ifdef NEW_IO
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c def_rsf_trdust
      subroutine def_rsf_trdust(fid)
!@sum def_rsf_trdust defines control info in restart files specifically for
!@+                  soil dust aerosols
!@auth Jan Perlwitz

      implicit none

!@var fid file id
      integer,intent(in) :: fid

      call defvar(grid,fid,dustDiagSubdd_acc%dustEmission(:,:,:)
     &     ,'dustEmission'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustEmission2(:,:,:)
     &     ,'dustEmission2'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustDepoTurb(:,:,:)
     &     ,'dustDepoTurb'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustDepoGrav(:,:,:)
     &     ,'dustDepoGrav'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustMassInPrec(:,:,:)
     &     ,' dustMassInPrec'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustSurfMixR(:,:,:)
     &     ,'dustSurfMixR'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustSurfConc(:,:,:),'
     &     dustSurfConc'//'(dist_im,dist_jm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustMass(:,:,:,:)
     &     ,'dustMass'//'(dist_im,dist_jm,lm,Ntm_dust)')
      call defvar(grid,fid,dustDiagSubdd_acc%dustConc(:,:,:,:)
     &     ,'dustConc'//'(dist_im,dist_jm,lm,Ntm_dust)')

      return
      end subroutine def_rsf_trdust

c new_io_trdust
      subroutine new_io_trdust(fid,iaction)
!@sum  new_io_trdust netcdf I/O specifically for soil dust aerosols
!@auth Jan Perlwitz

      implicit none

!@var fid file id
!@var iaction flag for reading or writing to file
      integer,intent(in) :: fid,iaction

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'dustEmission'
     &       ,dustDiagSubdd_acc%dustEmission(:,:,:))
        call write_dist_data(grid,fid,'dustEmission2'
     &       ,dustDiagSubdd_acc%dustEmission2(:,:,:))
        call write_dist_data(grid,fid,'dustDepoTurb'
     &       ,dustDiagSubdd_acc%dustDepoTurb(:,:,:))
        call write_dist_data(grid,fid,'dustDepoGrav'
     &       ,dustDiagSubdd_acc%dustDepoGrav(:,:,:))
        call write_dist_data(grid,fid,'dustMassInPrec'
     &       ,dustDiagSubdd_acc%dustMassInPrec(:,:,:))
        call write_dist_data(grid,fid,'dustSurfMixR'
     &       ,dustDiagSubdd_acc%dustSurfMixR(:,:,:))
        call write_dist_data(grid,fid,'dustSurfConc'
     &       ,dustDiagSubdd_acc%dustSurfConc(:,:,:))
        call write_dist_data(grid,fid,'dustMass'
     &       ,dustDiagSubdd_acc%dustMass(:,:,:,:))
        call write_dist_data(grid,fid,'dustConc'
     &       ,dustDiagSubdd_acc%dustConc(:,:,:,:))
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'dustEmission'
     &       ,dustDiagSubdd_acc%dustEmission(:,:,:))
        call read_dist_data(grid,fid,'dustEmission2'
     &       ,dustDiagSubdd_acc%dustEmission2(:,:,:))
        call read_dist_data(grid,fid,'dustDepoTurb'
     &       ,dustDiagSubdd_acc%dustDepoTurb(:,:,:))
        call read_dist_data(grid,fid,'dustDepoGrav'
     &       ,dustDiagSubdd_acc%dustDepoGrav(:,:,:))
        call read_dist_data(grid,fid,'dustMassInPrec'
     &       ,dustDiagSubdd_acc%dustMassInPrec(:,:,:))
        call read_dist_data(grid,fid,'dustSurfMixR'
     &       ,dustDiagSubdd_acc%dustSurfMixR(:,:,:))
        call read_dist_data(grid,fid,'dustSurfConc'
     &       ,dustDiagSubdd_acc%dustSurfConc(:,:,:))
        call read_dist_data(grid,fid,'dustMass'
     &       ,dustDiagSubdd_acc%dustMass(:,:,:,:))
        call read_dist_data(grid,fid,'dustConc'
     &       ,dustDiagSubdd_acc%dustConc(:,:,:,:))
      end select

      return
      end subroutine new_io_trdust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/
#endif /*NEW_IO*/

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_QUARZHEM)
c accSubddDust
      subroutine accSubddDust(dustDiagSubdd_acc)
!@sum  accSubddDust accumulates specific soil dust aerosol variables for
!@+                 subdaily diagnostics
!@auth Jan Perlwitz
      use atm_com, only : t,byMA,pk,pmid
      implicit none

      type(dustDiagSubdd),intent(inout) :: dustDiagSubdd_acc

      integer :: l,n,n1

      do n=1,Ntm_dust
        n1=n_soilDust+n-1
        dustDiagSubdd_acc%dustEmission(:,:,n)
     &       =dustDiagSubdd_acc%dustEmission(:,:,n)+dust_flux_glob(:,:,n
     &       )
#if (defined TRACERS_DUST) || (defined TRACERS_AMP)
        dustDiagSubdd_acc%dustEmission2(:,:,n)
     &       =dustDiagSubdd_acc%dustEmission(:,:,n)+dust_flux2_glob(:,:
     &       ,n)
#endif
#ifdef TRACERS_DRYDEP
        if (dodrydep(n1)) then
          dustDiagSubdd_acc%dustDepoTurb(:,:,n)
     &         =dustDiagSubdd_acc%dustDepoTurb(:,:,n)
     &         +depo_turb_glob(:,:,n1)
          dustDiagSubdd_acc%dustDepoGrav(:,:,n)
     &         =dustDiagSubdd_acc%dustDepoGrav(:,:,n)
     &         +depo_grav_glob(:,:,n1)
        end if
#endif
#ifdef TRACERS_WATER
        if (dowetdep(n1)) then
          dustDiagSubdd_acc%dustMassInPrec(:,:,n)
     &         =dustDiagSubdd_acc%dustMassInPrec(:,:,n)+trprec(n1,:,:)
        end if
#else
        dustDiagSubdd_acc%dustMassInPrec(:,:,n)
     &       =dustDiagSubdd_acc%dustMassInPrec(:,:,n)+trprec_dust(n,:,:)
#endif
        dustDiagSubdd_acc%dustSurfMixR(:,:,n)
     &       =dustDiagSubdd_acc%dustSurfMixR(:,:,n)+trcSurf(:,:,n1)
        dustDiagSubdd_acc%dustSurfConc(:,:,n)
     &       =dustDiagSubdd_acc%dustSurfConc(:,:,n)+trcSurfByVol(:,:,n1)
        dustDiagSubdd_acc%dustMass(:,:,:,n)=dustDiagSubdd_acc%dustMass(:
     &       ,:,:,n)+trm(:,:,:,n1)
        do l=1,lm
          dustDiagSubdd_acc%dustConc(:,:,l,n) =
     &         dustDiagSubdd_acc%dustConc(:,:,l,n) + trm(:,:,l,n1)
     &         *byMA(l,:,:)*1d2*pmid(l,:,:)/(rgas*t(:,:,l)*pk(l,:,:))
        end do
      end do

      return
      end subroutine accSubddDust
#endif /*TRACERS_DUST || TRACERS_MINERALS || TRACERS_QUARZHEM*/

#if (defined TRACERS_MINERALS) || (defined TRACERS_QUARZHEM)
c read_minfr_claquin1999_netcdf
      subroutine read_minfr_claquin1999_netcdf
!@sum read_minfr_claquin1999_netcdf  reads mineral fractions in soil from file
!@+                                  according to Claquin et al. JGR, 1999 data
!@auth jan perlwitz

      implicit none

      logical, dimension(Mtrac) :: qminfr = .false.

      integer :: m, n, ncid
      integer :: count(2), start(2)
      integer, dimension(mtrac) :: min_index, varid
      real(kind=8), dimension( grid%i_strt:grid%i_stop,
     &                 grid%j_strt:grid%j_stop ) :: work ! no halo

      call getDomainBounds( grid, j_strt=j_0, j_stop=j_1, i_strt=i_0, i_stop=i_1 )

      call write_parallel( 'Read from file MINFR', unit=6 )

      call check_netcdf( nf_open( 'MINFR', ncnowrit, ncid ),
     &     'nf_open called from read_minfr_claquin1999_netcdf ' //
     &     'in TRDUST_DRV.f' )

      start( 1 ) = i_0
      start( 2 ) = j_0

      count( 1 ) = 1 + ( i_1 - i_0 )
      count( 2 ) = 1 + ( j_1 - j_0 )

      m = 0
      do n = 1,Ntm_dust

        select case ( dust_names( n ) )

#ifdef TRACERS_MINERALS
        case( "ClayIlli" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralClayIllite',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 1
          qminfr( min_index( m ) ) = .true.

        case( "ClayKaol" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralClayKaolinite'
     &         , varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 2
          qminfr( min_index( m ) ) = .true.

        case( "ClaySmec" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralClaySmectite'
     &         , varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 3
          qminfr( min_index( m ) ) = .true.

        case( "ClayCalc" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralClayCalcite',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 4
          qminfr( min_index( m ) ) = .true.

        case( "ClayQuar" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralClayQuartz',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 5
          qminfr( min_index( m ) ) = .true.

        case( "Sil1Quar" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltQuartz',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 6
          qminfr( min_index( m ) ) = .true.

        case( "Sil1Feld" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltFeldspar',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 7
          qminfr( min_index( m ) ) = .true.

        case( "Sil1Calc" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltCalcite',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 8
          qminfr( min_index( m ) ) = .true.

        case( "Sil1Hema" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltHematite',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 9
          qminfr( min_index( m ) ) = .true.

        case( "Sil1Gyps" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltGypsum',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 10
          qminfr( min_index( m ) ) = .true.

#else /* TRACERS_MINERALS */
#ifdef TRACERS_QUARZHEM
        case( "Sil1QuHe" )
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltQuartz',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 6
          qminfr( min_index( m ) ) = .true.
          m = m + 1
          call check_netcdf( nf_inq_varid ( ncid, 'mineralSiltHematite',
     &         varid( m ) ), 'nf_inq_varid called from ' //
     &         'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
          min_index( m ) = 9
          qminfr( min_index( m ) ) = .true.
#endif /* TRACERS_QUARZHEM */
#endif

        end select

      end do                    ! Ntm_dust


      do m = 1,mtrac

        if ( .not. qminfr( min_index ( m ) ) ) cycle

        call check_netcdf( nf_get_vara_double ( ncid, varid ( m ), start
     &       , count, work ), 'nf_get_vara_double called from ' //
     &       'read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )
        minfr( i_0:i_1, j_0:j_1, min_index( m ) ) = max( work( i_0:i_1,
     &       j_0:j_1 ), 0.d0 )

      end do

      call check_netcdf( nf_close( ncid ), 'nf_close called ' //
     &     'from read_minfr_claquin1999_netcdf in TRDUST_DRV.f' )

      return
      end subroutine read_minfr_claquin1999_netcdf

c fraction_of_mineral_claquin1999
      subroutine fraction_of_mineral_claquin1999
!@sum fraction of mineral  determines fraction of each mineral for each bin
!@+                        from input mineral fractions of clay and silt
!@+                        based on data provided by Claquin et al., JGR 1999.
!@auth jan perlwitz

      implicit none

      integer :: i, j, m, n, n1
      real(kind=8) :: fac2, fac3, summinfr, y
      real(kind=8), dimension( grid%i_strt:grid%i_stop,
     &     grid%j_strt:grid%j_stop, ntm_dust ) :: work
      real(kind=8), dimension( grid%i_strt:grid%i_stop,
     &     grid%j_strt:grid%j_stop ) :: fac1
!@var hematiteAggrProb probability of iron oxide in aggregate with other mineral
      real(kind=8), dimension( grid%i_strt:grid%i_stop,
     &     grid%j_strt:grid%j_stop, Mtrac ) :: hematiteAggrProb

c**** Treatment of Quartz, iron oxide mineral (Hematite), and internally
c**** mixed Quartz-Hematite aggregates in grid box, when aggregates are
c**** switched on:
c**** 1. a) pureByTotalHematite - prescribed initial fraction of
c****       Hematite that stays a pure mineral relative to all Hematite
c****       (pure plus aggregating with Quartz).
c****    b) pureByTotalHematite - calculated as function of aggregation
c****       probability depending on mineral fractions in soil;
c****       probability of Hematite-other mineral aggegrate = (1 -
c****       fraction of Hematite) * fraction of other mineral.
c**** 2. frHemaInQuarAggr - prescribes the ratio of Hematite mass in
c****    internally mixed aggregate to total Quartz-Hematite aggregate
c****    mass.
c**** 3. pure Quartz fraction equals initial Quartz fraction minus
c****    fraction of the available Quartz fraction that goes into
c****    aggregates. The Quartz fraction needed for aggregate formation
c****    increases with decreasing parameter frHemaInQuarAggr,
c****    decreasing the pure Quartz fraction.
c**** 4. pure Hematite fraction equals prescribed initial pure Hematite
c****    fraction plus Hematite fraction that is not forming aggregates
c****    anymore, after the available Quartz fraction has been
c****    exhausted.
c**** 5. the fraction of internally mixed Quartz-Hematite aggregates
c****    depends on the fraction of Hematite that is allowed to form
c****    aggregates and the available Quartz fraction. It can't be
c****    larger than the initial Quartz fraction plus the available
c****    Hematite fraction mixed in according to parameter
c****    frHemaInQuarAggr.

      if (  calcMineralAggrProb == 1 ) then
        hematiteAggrProb = 0.d0
! loop over silt mineral fractions
        m = 5
        do while ( m < Mtrac )
          m = m + 1
          if ( m == 9 ) cycle ! skip Hematite
          do j = j_0,j_1
            do i = i_0,i_1
! Although the original fractions of the minerals in soil add up to 1,
! this is not the case anymore, if the fractions are weighted, e.g., by
! multiplying with a brittle factor, to account for different
! fractionating of the minerals during emission. Therefore the fractions
! are renormalized to 1 using summinfr for calculating the aggregation
! probabilities.
              summinfr = sum( minfr( i, j, 6:mtrac ) ) + tiny( sum(
     &             minfr( i, j, 6:mtrac ) ) )
              hematiteAggrProb( i, j, m ) = (1.d0 - minfr( i, j, 9 ) /
     &             summinfr) * minfr( i, j, m ) / summinfr
            end do
          end do
        end do
! Since Quartz-iron oxide is currently the only aggregate, all
! aggregation probabilities are lumped together as probability of
! Quartz-iron oxide formation.
        fac1( i_0:i_1, j_0:j_1 ) = sum( hematiteAggrProb( i_0:i_1,
     &       j_0:j_1, 6:mtrac), dim=3 )
      else
        fac1 = 1.d0 - pureByTotalHematite
      end if

      fac2 = 1.d0 / (frHemaInQuarAggr + tiny( frHemaInQuarAggr ))
      y = 1.d0 - frHemaInQuarAggr
      fac3 = 1.d0 / (y + tiny( y ))

      do n = 1,ntm_dust

        select case( dust_names( n ) )

        case( 'ClayIlli', 'ClayKaol', 'ClaySmec', 'ClayCalc', 'ClayQuar'
     &       , 'Sil1Quar', 'Sil1Feld', 'Sil1Calc', 'Sil1Hema',
     &       'Sil1Gyps')
          n1 = n
        case( 'Sil2Quar', 'Sil2Feld', 'Sil2Calc', 'Sil2Hema',
     &         'Sil2Gyps')
          n1 = n - 5
        case( 'Sil3Quar', 'Sil3Feld', 'Sil3Calc', 'Sil3Hema',
     &         'Sil3Gyps')
          n1 = n - 10

        end select

        select case( dust_names( n ) )

#ifdef TRACERS_MINERALS
        case ( 'Sil1Quar', 'Sil2Quar', 'Sil3Quar' )
          do j = j_0,j_1
            do i = i_0,i_1
! initial Quartz fraction
              work( i, j, n ) = minfr( i, j, n1 )
#ifdef TRACERS_QUARZHEM
! calculate the pure Quartz fraction from initial Quartz fraction minus
! the fraction used for Quartz-Hematite aggregates up to the exhaustion
! of available Quartz
     &             - min( fac1( i, j ) * minfr ( i, j, 9 ) * (fac2 -
     &             1.d0), minfr( i, j, n1 ) )
#endif
            end do
          end do
#ifdef TRACERS_QUARZHEM
        case ( 'Sil1Hema', 'Sil2Hema', 'Sil3Hema' )
! calculate the pure Hematite fraction as sum of non-aggregated Hematite
! plus leftover Hematite after exhaustion of Quartz available for
! aggregation with Hematite
          do j = j_0,j_1
            do i = i_0,i_1
              work( i, j, n ) = (1.d0 - fac1( i, j )) * minfr( i, j, n1
     &             ) + max( fac1( i, j ) * minfr( i, j, n1 ) - minfr( i,
     &             j, 6 ) * fac3 / fac2, 0.d0 )
            end do
          end do
#endif
#endif
#ifdef TRACERS_QUARZHEM
        case ( 'Sil1QuHe', 'Sil2QuHe', 'Sil3QuHe' )
! calculate the fraction of Quartz-Hematite aggregates from the Hematite
! fraction up to the exhaustion of available Quartz
          do j = j_0,j_1
            do i = i_0,i_1
              work( i, j, n ) = min( fac1( i, j ) * minfr( i, j, 9 ) *
     &             fac2, minfr( i, j, 6 ) * fac3 )
            end do
          end do
#endif
#ifdef TRACERS_MINERALS
        case default
! non-aggegrated minerals
          work( i_0:i_1, j_0:j_1, n ) = minfr( i_0:i_1, j_0:j_1, n1 )
#endif

        end select

        mineralFractions( i_0:i_1, j_0:j_1, n ) = work( i_0:i_1, j_0:j_1
     &       , n )

      end do

      return
      end subroutine fraction_of_mineral_claquin1999
#endif /*  TRACERS_MINERALS || TRACERS_QUARZHEM */

      subroutine check_netcdf( status, infostring )
!@sum check_netcdf  checks for netcdf error messages
!@auth jan perlwitz

      implicit none
      include 'netcdf.inc'

      integer, intent(in) :: status
      character(len=*), intent(in) :: infostring

      character(len=256) :: outstring

      if( status /= nf_noerr ) then
        outstring = trim( nf_strerror( status ) ) // ' in ' // trim(
     &       infostring )
        call stop_model( trim( outstring ), 255 )
      end if

      return
      end subroutine check_netcdf

      end module trdust_drv
