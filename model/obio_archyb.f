#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      subroutine obio_archyb(nn,dpav,temav,salav,th3av,
     .                       dpmxav,oiceav)
c
c --- write archive file for time level n to flnm ( b i n a r y  hycom fmt)
c
      USE JulianCalendar_mod, only : jdendofm
      USE MODEL_COM, only : modelEclock,
     *  itime,iyear1,nday,aMON,xlabel,lrunid
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE HYCOM_SCALARS, only : nstep,time,lp,theta,onem
     &     ,thref,baclin
#if (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
     || (defined TRACERS_ZEBRA)
     .     ,diag_counter
#endif
      USE HYCOM_DIM_GLOB, only : jj,jdm,kk,ntrcr,ii,idm,kdm
      USE HYCOM_ARRAYS_GLOB, only: tracer,temp,saln,p,dpmixl
#if (defined TRACERS_AGE_OCEAN) || (defined TRACERS_OCEAN_WATER_MASSES) \
     || (defined TRACERS_ZEBRA)
     .  ,plevav,tracav
#endif
#ifdef TRACERS_OceanBiology    
      USE obio_com, only : pCO2_glob,ao_co2flux_glob,pCO2av
     .                    ,ao_co2fluxav,diag_counter,tracav,plevav
     .                    ,cexp_glob,pp2tot_day_glob
     .                    ,cexpav,pp2tot_dayav
#ifdef TRACERS_Alkalinity
     .                    ,caexp_glob,caexpav
#endif
#endif
c
      implicit none
      real*8, dimension(idm,jdm) :: dpmxav,oiceav
      real*8, dimension(idm,jdm,kdm) :: dpav,temav,salav,th3av
      real factor
      integer i,j,k,l,kn,nn
c
      integer no,nop,nt
      character flnm*40,intvl*3,title*80
      integer :: year, month, dayOfYear, date, hour
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)

      call getdte(Itime,Nday,Iyear1,year,month,dayOfYear,date,hour,amon)
c --- check if ogcm date matches agcm date
      if (nstep.eq.1) then
        write(flnm,'(a3,i4.4,2a)') amon,0,'.obio',xlabel(1:lrunid)
!     elseif (abs((itime+1.)/nday-time).gt.1.e-5) then
!       write(0,*) 'obio: mismatching archive date in agcm/ogcm=',
!    .     (itime+1.)/nday,time
!       stop 'obio_archyb: mismatching archive date'
      else
        write(flnm,'(a3,i4.4,2a)') amon,year,'.obio',xlabel(1:lrunid)
      endif
c
      if (date.lt.100) then
        write (intvl,'(i3.3)') date
      else
        stop 'date >100'
      endif
c
      nop=12
      write (lp,'(a/9x,a)') 'storing history data in',flnm

      factor=baclin/(date*SECONDS_PER_DAY)
c
      open (unit=nop,file=flnm,status='unknown',
     .      form='unformatted')

      ! ocean layer depth
      do k=1,kk
        write(title,'(a,i4)')'depth, k=',k
        call write2giss(nop,p(:,:,k)/onem,title)
      enddo

      ! ocean temperature
      do k=1,kk
      kn=k+nn
        write(title,'(a,i4)')'temp, k=',k
        call write2giss(nop,temp(:,:,kn),title)
      enddo

      ! ocean salinity
      do k=1,kk
      kn=k+nn
        write(title,'(a,i4)')'saln, k=',k
        call write2giss(nop,saln(:,:,kn),title)
      enddo

      !mld
        write(title,'(a)')'MLD (m)'
        call write2giss(nop,dpmixl(:,:,2)/onem,title)

      !tracer
      do nt=1,ntrcr
      do k=1,kk
        write(title,'(a,i4,a,i4)')'tracer, nt=',nt,', k=',k
        call write2giss(nop,tracer(:,:,k,nt),title)
      enddo
      enddo

#ifdef TRACERS_OceanBiology    
      !pco2
        write(title,'(a)')'pCO2, water'
        call write2giss(nop,pCO2_glob,title)

      !ao co2 flux
        write(title,'(a)')'AO CO2 flux'
        call write2giss(nop,ao_co2flux_glob,title)

      !cexp
        write(title,'(a)')'cexp, water'
        call write2giss(nop,cexp_glob,title)

#ifdef TRACERS_Alkalinity
      !caexp
        write(title,'(a)')'caexp, water'
        call write2giss(nop,caexp_glob,title)
#endif

      !pp2tot_day
        write(title,'(a)')'pp2tot_day, water'
        call write2giss(nop,pp2tot_day_glob,title)
#endif


! time averaged arrays

      !dpav
      do k=1,kk
        write(title,'(a,i4)')'dpav, k=',k
        call write2giss(nop,dpav(:,:,k)/onem,title)
      enddo
   
      !tempav
      do k=1,kk
        write(title,'(a,i4)')'tempav, k=',k
        call write2giss(nop,temav(:,:,k),title)
      enddo
   
      !salav
      do k=1,kk
        write(title,'(a,i4)')'salav, k=',k
        call write2giss(nop,salav(:,:,k),title)
      enddo
   
      !th3av
      do k=1,kk
        write(title,'(a,i4)')'th3av, k=',k
        call write2giss(nop,th3av(:,:,k),title)
      enddo
   
      !dpmxav
        write(title,'(a)')'dpmxav'
        call write2giss(nop,dpmxav(:,:)/onem,title)
   
      !oiceav
        write(title,'(a)')'oiceav'
        call write2giss(nop,oiceav(:,:),title)
   
      !tracav
      !no need to divide by diag_counter because pressumably
      !this is already done when divide by plevav
      do nt=1,ntrcr
      where(plevav(:,:,:).ne.0.d0)
       tracav(:,:,:,nt)=tracav(:,:,:,nt)/plevav(:,:,:)
      endwhere
      do k=1,kk
        write(title,'(a,i4,a,i4)')'tracav, nt=',nt,', k=',k
        call write2giss(nop,tracav(:,:,k,nt),title)
      enddo
      enddo
   
#ifdef TRACERS_OceanBiology    
      !pco2av
      if (diag_counter .ne. 0.d0) then
        pco2av=pco2av/diag_counter
        pp2tot_dayav=pp2tot_dayav/diag_counter
        cexpav=cexpav/diag_counter
#ifdef TRACERS_Alkalinity
        caexpav=caexpav/diag_counter
#endif
      endif
        write(title,'(a)')'pCO2av'
        call write2giss(nop,pCO2av,title)
        write(title,'(a)')'pp2tot_dayav'
        call write2giss(nop,pp2tot_dayav,title)
        write(title,'(a)')'cexpav'
        call write2giss(nop,cexpav,title)
#ifdef TRACERS_Alkalinity
        write(title,'(a)')'caexpav'
        call write2giss(nop,caexpav,title)
#endif
      !ao_co2fluxav
      if (diag_counter .ne. 0.d0) then
        ao_co2fluxav=ao_co2fluxav/diag_counter
      endif
        write(title,'(a)')'ao_co2fluxav'
        call write2giss(nop,ao_co2fluxav,title)
#endif

!zero out for next diagnostic period
      diag_counter= 0.
#ifdef TRACERS_OceanBiology    
      pco2av = 0.
#endif

      return
      end subroutine obio_archyb

      module write2giss_mod
      use resolution, only : im,jm
      use fluxes, only : focean
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : pack_data
      implicit none
      real*8, dimension(im,jm) :: focean_glob
      contains
      subroutine write2giss_init
      call pack_data(grid, focean, focean_glob)
      end subroutine write2giss_init
      end module write2giss_mod

      subroutine write2giss(nop,array_o,title)

      USE HYCOM_DIM_GLOB, only : iia,jja,iio,jjo
      USE HYCOM_CPLER, only: ssto2a_global
      USE HYCOM_SCALARS, only : onem
      USE write2giss_mod, only: focean=>focean_glob

      integer nop
      real*4 array4(iia,jja)
      real*8 array8(iia,jja),array_o(iio,jjo)
      character title*80

        !convert to atmosgrid
        call ssto2a_global(array_o,array8)
        array4=array8
        do j=1,jja
        do i=1,iia
         if (focean(i,j).le.0.) array4(i,j)=1.e33
        enddo
        enddo
        write (nop)title,array4 

      end subroutine write2giss
