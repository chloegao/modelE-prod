#include "rundeck_opts.h"

#ifdef OBIO_ON_GARYocean
      subroutine obio_model(atm)
#else
      subroutine obio_model(mm,atm)
#endif

!@sum  OBIO_MODEL is the main ocean bio-geo-chem routine 
!@auth Natassa Romanou/Watson Gregg

      USE Dictionary_mod
      USE obio_dim
      USE obio_incom
      USE obio_forc, only: solz,tirrq,Ed,Es
     .                    ,rmud,atmFe,avgq,ihra,sunz
     .                    ,wind
     .                    ,alk
     .                    ,tirrq3d
     .                    ,surfN
      USE obio_com,  only: dobio,gcmax,day_of_month,hour_of_day
     .                    ,temp1d,dp1d,obio_P,det,car,avgq1d
     .                    ,ihra_ij,gcmax1d,atmFe_ij,covice_ij
     .                    ,P_tend,D_tend,C_tend,saln1d
     .                    ,pCO2_ij,p1d,wsdet
     .                    ,rhs,alk1d
     .                    ,tzoo,tfac,rmuplsr,rikd,wshc,Fescav
     .                    ,tzoo2d,tfac3d,rmuplsr3d,rikd3d
     .                    ,wshc3d,Fescav3d 
     .                    ,acdom,pp2_1d,pp2tot_day,pp2diat_day
     .                    ,pp2chlo_day,pp2cyan_day,pp2cocc_day
     .                    ,tot_chlo,acdom3d
     .                    ,itest,jtest
     .                    ,obio_ws
     .                    ,cexp,flimit,kzc
     .                    ,rhs_obio,chng_by,Kpar,Edz,Esz,Euz
     .                    ,delta_temp1d
#ifdef TOPAZ_params
     .                    ,ca_det_calc1d
#endif
      use obio_com, only: caexp
      use obio_com, only: build_ze

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      use obio_com, only: rnitrmflo_loc
      use obio_com, only: rnitrconc_loc
#endif
#ifdef DIC_RUNOFF
      use obio_com, only: rdicconc_loc
#endif
#ifdef DOC_RUNOFF
      use obio_com, only: rdocconc_loc
#endif
#ifdef SILI_RUNOFF
      use obio_com, only: rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
      use obio_com, only: rironconc_loc
#endif
#ifdef POC_RUNOFF
      use obio_com, only: rpocconc_loc
#endif
#ifdef ALK_RUNOFF
      use obio_com, only: ralkconc_loc
#endif
#endif

      use runtimecontrols_mod, only: tracers_alkalinity
      use obio_com, only: tracer,nstep0
#ifdef OBIO_ON_GARYocean
      use obio_com, only: obio_deltat
#endif
      USE obio_diag, only : ij_pCO2,ij_dic,ij_nitr,ij_diat
     .                 ,ij_amm,ij_sil,ij_chlo,ij_cyan,ij_cocc,ij_herb
     .                 ,ij_doc,ij_iron,ij_alk,ij_Ed,ij_Es,ij_pp,ij_dayl
     .                 ,ij_cexp,ij_lim,ij_wsd,ij_ndet,ij_xchl
     .                 ,ij_sunz,ij_solz
     .                 ,ij_pp1,ij_pp2,ij_pp3,ij_pp4
     .                 ,ij_rhs,ij_flux,ij_fca

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
      USE obio_diag, only: ij_rnitrconc
!     .                 ,ij_rnitrmflo
#endif
#ifdef DIC_RUNOFF
      USE obio_diag, only: ij_rdicconc
#endif
#ifdef DOC_RUNOFF
      USE obio_diag, only:  ij_rdocconc
#endif
#ifdef SILI_RUNOFF
      USE obio_diag, only:  ij_rsiliconc
#endif
#ifdef IRON_RUNOFF
      USE obio_diag, only:  ij_rironconc
#endif
#ifdef POC_RUNOFF
      USE obio_diag, only:  ij_rpocconc
#endif
#ifdef ALK_RUNOFF
      USE obio_diag, only:  ij_ralkconc
#endif
#endif
      USE obio_diag, only : oijl=>obio_ijl,ijl_avgq,ijl_kpar,ijl_dtemp
      use obio_diag, only: oij=>obio_ij
      use ocalbedo_mod, only: ocalbedo

      USE MODEL_COM, only: modelEclock
     . ,itime,iyear1,aMON,dtsrc
     . ,xlabel,lrunid
      use JulianCalendar_mod, only: jdendofm
      use TimeConstants_mod, only: HOURS_PER_DAY

      USE FILEMANAGER, only: openunit,closeunit,file_exists

      USE obio_com, only : co2flux


#ifdef OBIO_ON_GARYocean
      USE MODEL_COM,  only : nstep=>itime,itimei
      USE OCEAN,       only : oLON_DG,oLAT_DG
      USE CONSTANT,   only : grav
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : kdm=>lmo,dzo
      USE OFLUXES,    only : oice=>oRSI,oAPRESS
      USE OCEAN,      only : g0m,s0m,mo,dxypo,ip=>focean,lmm
     .                      ,trmo,txmo,tymo,tzmo
#else
      USE hycom_dim
      USE hycom_arrays, only: tracer_h=>tracer,dpinit,temp,saln,oice
     .                            ,p,dpmixl,latij,lonij,scp2
      USE  hycom_arrays_glob, only: latij_glob=>latij,lonij_glob=>lonij
      USE hycom_scalars, only: nstep,onem,time,lp,baclin
      USE obio_com, only: ao_co2fluxav_loc,
     .     pCO2av_loc,pp2tot_dayav_loc,cexpav_loc,caexpav_loc,
     .     pp2diat_dayav_loc,pp2chlo_dayav_loc,pp2cyan_dayav_loc,
     .     pp2cocc_dayav_loc,pHav,pHav_loc

      USE obio_com, only: diag_counter
#endif
      use obio_com, only: ze

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      use TimerPackage_mod
      use obio_diag, only: reset_obio_diag

      use exchange_types, only : atmocn_xchng_vars
      implicit none
      type(atmocn_xchng_vars) :: atm

      REAL*4, parameter  :: obio_tr_mm(16)= (/ 14., 14., 28.055,
     &     55.845, 1., 1., 1., 1., 1., 14., 14., 28.055, 55.845,
     &     12., 12., 1. /)
      integer i,j,k,l,km,mm

      integer ihr,ichan,iyear,nt,ihr0,lgth,kmax
      integer ll,ilim
      real    tot,dummy(6),dummy1
      real    rod(nlt),ros(nlt)
#ifdef OBIO_ON_GARYocean
      Real*8,External   :: VOLGSP
      real*8 temgs,g,s,temgsp,pres
      real*8 time,dtr,ftr,rho_water
      integer i_0,i_1,j_0,j_1
#endif
      real*8 trmo_unit_factor(kdm,ntrac)
      integer :: idx_co2

      logical vrbos,noon,errcon
      integer :: year, month, dayOfYear, date, hour
      real, allocatable, dimension(:), save :: eda_frac, esa_frac
      integer :: iu_bio
      logical, save :: initialized=.false.

      if(.not.dobio) return

      call start(' obio_model')
!--------------------------------------------------------
      diagno_bio=.false.
c
      call modelEclock%get(year=year, month=month, date=date,
     .  hour=hour, dayOfYear=dayOfYear)

      if (JDendOfM(month).eq.dayOfYear.and.hour.eq.12) then
          if (mod(nstep,2).eq.0)
     .    diagno_bio=.true. ! end of month,mid-day
      endif

!Cold initialization

      call start(' obio_init')
      !nstep0=0 : cold initialization
      !nstep0>0 : warm initialization, is the timestep of current restart run
      !nstep    : current timestep   

      print*, 'nstep,nstep0 =',
     .         nstep,nstep0

      call build_ze
      if (nstep0==0) then
        print*, 'COLD INITIALIZATION....'

        if (AM_I_ROOT()) write(*,'(a)')'BIO:Ocean Biology starts ....'


        call obio_init
      !tracer array initialization.
      !note: we do not initialize obio_P,det and car

#ifndef OBIO_ON_GARYocean
        ao_co2fluxav_loc  = 0.
        pCO2av_loc = 0
        pp2tot_dayav_loc = 0
        pp2diat_dayav_loc = 0
        pp2chlo_dayav_loc = 0
        pp2cyan_dayav_loc = 0
        pp2cocc_dayav_loc = 0
        cexpav_loc = 0
        caexpav_loc = 0
        pHav_loc = 0
#endif
        call obio_bioinit
      endif   !cold restart



!Warm initialization

#ifdef OBIO_ON_GARYocean
      time = float(nstep)
#endif
      if ((nstep0>0).and..not.initialized) then
         write(*,'(a)')'For restart runs.....'
         write(*,'(a,2i9,f10.3)')
     .            'nstep0,nstep,time=',nstep0,nstep,time
         call obio_init
         print*,'WARM INITIALIZATION'

      endif !for restart only
      call stop(' obio_init')

#ifdef OBIO_ON_GARYocean
      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP
#endif

      if (AM_I_ROOT()) then
#ifdef OBIO_ON_GARYocean
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      oLAT_DG(jtest,2),oLON_DG(itest,2)
#else
      write(*,'(a,2i5,2e12.4)')'TEST POINT at: ',itest,jtest,
     .      latij_glob(itest,jtest,3),lonij_glob(itest,jtest,3)
#endif
       endif

       call sync_param( "solFe", solFe)
       print*, 'solfe=',solFe

!--------------------------------------------------------

       day_of_month=date
       hour_of_day =hour

       !if (time.kmod(nstep,24*30).eq.0) then
       !   day_of_month=1
       !else
       !   day_of_month=day_of_month+1
       !endif

       !if (mod(nstep,24).eq.0) then
       !  hour_of_day=1
       !else
       !  hour_of_day=hour_of_day+1
       !endif
 
      if (AM_I_ROOT()) then
       write(*,'(a,i15,1x,f9.3,2x,3i5)')
     .    'BIO: nstep,time,day_of_month,hour_of_day,dayOfYear=',
     .    nstep,time,day_of_month,hour_of_day,dayOfYear
      endif

         ihr0 = int(hour_of_day/2)

      if (diagno_bio) then
      endif  !diagno_bio

#ifdef OBIO_ON_GARYocean
      write(*,'(/,a,2i5,2e12.4)')'obio_model, test point=',
     .      itest,jtest,oLON_DG(itest,1),oLAT_DG(jtest,1)
#endif

      !print out tracer integrals just before main loop

#ifndef OBIO_ON_GARYocean     /* HYCOM only */
      diag_counter=diag_counter+1
#endif

      if ((dayofyear==1+jdendofm(month-1)).and.
     &       modeleclock%isbeginningofday()) call reset_obio_diag

      call start('  obio main loop')

      atm%chl_defined=.true.
       do j=j_0,j_1
       do i=i_0,i_1
       dp1d(:) = 0.
       IF(ip(I,J)==0) cycle

cdiag  if (nstep.eq.1)
cdiag. write(*,'(a,3i5,2e12.4)')'obio_model, step,i,j=',nstep,i,j,
cdiag.          olon_dg(i,1),olat_dg(j,1)

       vrbos=.false.
       if (i.eq.itest.and.j.eq.jtest) vrbos=.true.

       !fill in reduced rank arrays
       ihra_ij=ihra(i,j)
       !!covice_ij=covice(i,j)  !for standalone hycom
       covice_ij=oice(i,j)      !for modelE-hycom
       pCO2_ij=atm%gtracer(atm%n_co2n,i,j)
     
#ifdef OBIO_ON_GARYocean
       pres = oAPRESS(i,j)    !surface atm. pressure
       do k=1,lmm(i,j)
         pres=pres+MO(I,J,k)*GRAV*.5
         g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
!!!!     temp1d(k)=TEMGS(g,s)           !potential temperature
         temp1d(k)=TEMGSP(g,s,pres)     !in situ   temperature
          saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
          !dp1d(k)=dzo(k)               !thickenss of each layer in meters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! array tracer has units [mole]/m3. convert to/from trmo with units kg as follows:
! trmo = tracer * MB*1e-3/rho_water     * mo * dxypo
! note that occassionally tracer is in milimoles and othertimes in micromoles
! trmo_unit_factor below has units m^3 and depends on which 
! layer you are: trmo = tracer * factor,  tracer=trmo/factor
! txmo,tzmo,tymo are all computed based on trmo
! this should be done at every timestep except when COLD INITIALIZATION
! because then trmo=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         rho_water = 1d0/VOLGSP(g,s,pres)
!!!!!!   rho_water = 1035.
         dp1d(k)=MO(I,J,K)/rho_water   !local thickenss of each layer in meters
         if(vrbos.and.k.eq.1)write(*,'(a,4e12.4)')
     .             'obio_model,t,s,p,rho= '
     .             ,temp1d(k),saln1d(k),dp1d(k),rho_water

         !add missing part of density to get to the bottom of the layer
         !now pres is at the bottom of the layer
         pres=pres+MO(I,J,k)*GRAV*.5

         do nt=1,ntrac
         
                trmo_unit_factor(k,nt) =  1d-3*1d-3*obio_tr_mm(nt)        ! milimoles/m3=> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg

         if (nt.eq.4.or.nt.eq.13) 
     .          trmo_unit_factor(k,nt) =  1d-6*1d-3*obio_tr_mm(nt)        ! nanomoles/lt=> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg

         if (nt.eq.11)
     .          trmo_unit_factor(k,nt) =  1d-6 *1d-3/1d-3                 ! micro-grC/lt -> kg/m3
     .                      *  MO(I,J,k)*DXYPO(J)/rho_water               ! kg/m3 => kg

         if (nt.ge.5.and.nt.le.9) 
     .          trmo_unit_factor(k,nt) =  
     .                          cchlratio * 1d-6                          ! miligr,chl/m3=> kg,C/m3
     .                       *  MO(I,J,k)*DXYPO(J)/rho_water              ! kg/m3 => kg

#ifdef TRACERS_Alkalinity
           if (nt.eq.16)    !factor for alkalinity
     .         trmo_unit_factor(k,nt) = 1d-6*1d-3*obio_tr_mm(nt)      ! umol/kg=micro-mol/kg=> kg,trac/kg,air
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,air=> kg,trac

#ifdef TOPAZ_params
           if (nt.eq.17)    !factor for ca_det
     .         trmo_unit_factor(k,nt) = 1d-6*1d-3*obio_tr_mm(nt)      ! umol/kg=micro-mol/kg=> kg,trac/kg,air
     .                                *  MO(I,J,k)*DXYPO(J)           ! kg,trac/kg,air=> kg,trac
#endif
#endif

           if (nstep0>0) then
              tracer(i,j,k,nt) = trmo(i,j,k,nt) / trmo_unit_factor(k,nt)
           endif
         enddo


#ifdef IncreaseN
       tracer(i,j,1,1) = 100.d0 * tracer(i,j,1,1)
#endif
#ifdef Relax2SurfN
#ifdef Relax2maxSurfN
       tracer(i,j,1,1) = max(tracer(i,j,1,1),surfN(i,j))
#else
       tracer(i,j,1,1) = surfN(i,j)
#endif
#endif

#else
       trmo_unit_factor=1
       if (nstep0>0) tracer(i,j,:,:)=tracer_h(i,j,:,:)/trmo_unit_factor
       do k=1,kdm
         km=k+mm
         temp1d(k)=temp(i,j,km)
         saln1d(k)=saln(i,j,km)
         dp1d(k)=dpinit(i,j,k)/onem
#endif
         avgq1d(k)=avgq(i,j,k)
         gcmax1d(k)=gcmax(i,j,k)
         tirrq(k)=tirrq3d(i,j,k)
#ifdef TRACERS_Alkalinity
         nt=1
         alk1d(k)=tracer(i,j,k,ntyp+n_inert+ndet+ncar+1)
#ifdef TOPAZ_params
         nt=2
         ca_det_calc1d(k)=tracer(i,j,k,ntyp+n_inert+ndet+ncar+nt)
#endif
#else
              !NOT for INTERACTIVE alk
         alk1d(k)=alk(i,j,k)
#endif
              !----daysetbio/daysetrad arrays----!
         tzoo=tzoo2d(i,j)
         tfac(k)=tfac3d(i,j,k)
         do nt=1,nchl
           rmuplsr(k,nt)=rmuplsr3d(i,j,k,nt)
           rikd(k,nt)=rikd3d(i,j,k,nt)
         enddo
         wshc(k)=wshc3d(i,j,k)
         Fescav(k)=Fescav3d(i,j,k)
         do nt=1,nlt
           acdom(k,nt)=acdom3d(i,j,k,nt)
         enddo
              !----daysetbio arrays----!
         do nt=1,ntyp+n_inert
           obio_P(k,nt)=tracer(i,j,k,nt)
         enddo
         do nt=1,ndet
           det(k,nt)=tracer(i,j,k,ntyp+n_inert+nt)
         enddo
         do nt=1,ncar
           car(k,nt)=tracer(i,j,k,ntyp+n_inert+ndet+nt)
         enddo
       enddo  !k=1,kdm or lmm

       p1d(1)=0.
       do k=2,kdm+1
          p1d(k)=p1d(k-1)+dp1d(k-1)    !in meters
       enddo

      !if(vrbos) write(*,'(a,15e12.4)')'obio_model, strac conc:',
      if(vrbos) write(*,*)'obio_model, strac conc:',
     .    obio_P(1,:),det(1,:),car(1,:)

#ifdef TOPAZ_params
      if(vrbos) write(*,'(a,2e12.4)')'obio_model, sALK  conc:',
     . alk1d(1),ca_det_calc1d(1)
#endif

#ifdef OBIO_ON_GARYocean
      kmax = lmm(i,j)
cdiag if(vrbos)write(*,'(a,4i5)')'nstep,i,j,kmax= ',
cdiag.         nstep,i,j,kmax

#else

! --- find number of non-massless layers in column
      do kmax=1,kdm
c       if (dp1d(kmax).lt.1.e-2) exit
        if (dp1d(kmax).lt.1.) exit
      end do
      kmax=kmax-1
cdiag write(*,'(a,4i5)')'nstep,i,j,kmax= ',nstep,i,j,kmax

      !ensure that all massless layers have same concentration
      !as last mass one
      do k=kmax+1,kdm
       do nt=1,ntyp+n_inert
        if(obio_P(kmax,nt).le.0.)obio_P(kmax,nt)=1.e-8
        obio_P(k,nt)=obio_P(kmax,nt)
       enddo
       do nt=1,ndet
        det(k,nt)=   det(kmax,nt)
       enddo
       do nt=1,ncar
        car(k,nt)=   car(kmax,nt)
       enddo
      enddo
#endif

      !find layer index for zc, max depth of sinking phytoplankton
      kzc = 1
      do k=kmax+1,1,-1
           if (p1d(k).gt.zc) kzc = k
      enddo
      if (kzc.lt.1) kzc=1
      if (kzc.gt.kmax) kzc=kmax

      if (vrbos) write(*,*) 'compensation depth, kzc = ',kzc


       !solz is read inside hycom.f and forfun.f
!      do ihr=1,12
!       solz2(ihr)=solz_all(i,j,ihr,l0)*w0
!    .            +solz_all(i,j,ihr,l1)*w1
!    .            +solz_all(i,j,ihr,l2)*w2
!    .            +solz_all(i,j,ihr,l3)*w3

!       sunz2(ihr)=acos(solz2(ihr))*rad   !in degs
!      enddo

!        solz=solz2(ihr0)     !because of bi-hourly values
!        sunz=sunz2(ihr0)     !in degs
         solz=atm%cosz1(i,j) !osolz(i,j)      !read instead from modelE
         sunz=acos(solz)*rad  !in degs


!      wind=wndspd(i,j,l0)*w0+wndspd(i,j,l1)*w1
!    .     +wndspd(i,j,l2)*w2+wndspd(i,j,l3)*w3

       wind=atm%wsavg(i,j) !owind(i,j)

!      atmFe_ij=atmFe_all(i,j,l0)*w0 + atmFe_all(i,j,l1)*w1
!    .         +atmFe_all(i,j,l2)*w2 + atmFe_all(i,j,l3)*w3

       !atmospheric deposition iron 
       atmFe_ij=atmFe(i,j,month)
#ifdef zero_ironflux
        atmFe_ij=0.d0
#endif
       if (vrbos) then
         write(*,'(/,a,3i5,4e12.4)')'obio_model, forcing: ',
     .   nstep,i,j,solz,sunz,wind,atmFe_ij
       endif

#ifdef Relax2SurfN
       if (vrbos) then
         write(*,'(/,a,3i5,e12.4)')'obio_model, surfN: ',
     .   nstep,i,j,surfN(i,j)
       endif
#endif

       idx_co2=atm%gasex_index%getindex(atm%n_co2n)
       if (idx_co2>0) co2flux=atm%trgasex(idx_co2, i, j)

       !------------------------------------------------------------
       !at the beginning of each day only
       if ((hour_of_day.le.1).or..not.initialized) then

          if (day_of_month.eq.1)ihra_ij=1
          call obio_daysetrad(vrbos,i,j)
          ihra_ij = 0
          call obio_daysetbio(vrbos,i,j)

             !fill in the 3d arrays to keep in memory for rest of the day
             !when daysetbio is not called again
             tzoo2d(i,j)=tzoo
             do  k=1,kdm
               tfac3d(i,j,k)=tfac(k)
             do nt=1,nchl
               rmuplsr3d(i,j,k,nt)=rmuplsr(k,nt)
               rikd3d(i,j,k,nt)=rikd(k,nt)
             enddo
               wshc3d(i,j,k)=wshc(k)
               Fescav3d(i,j,k)=Fescav(k)
             do nt=1,nlt
               acdom3d(i,j,k,nt)=acdom(k,nt)
             enddo
             enddo

         if (day_of_month.eq.1)ihra_ij=1


cdiag    if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       nitr    ',
cdiag.'   ammo     sili     iron',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                        obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrsetrad   dpth     dp       diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p1d(k+1),dp1d(k),obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                        obio_P(k,8),obio_P(k,9),k=1,kdm)
cdiag    endif
 103     format(i9,2i5,a,a/(25x,i3,6(1x,es9.2)))
 104     format(i9,2i5,a,a/(25x,i3,7(1x,es9.2)))

       endif   !end of calculations for the beginning of day


       do ichan = 1,nlt
         Ed(ichan) = 0.0
         Es(ichan) = 0.0
       enddo
       rmud = 0.0
       iyear=2001


         !compute rod and ros only here. not ocean albedo.
         !ocean albedo is computed in ALBEDO.f
         !have to have hygr =  .true. 
       call ocalbedo(wind,solz,dummy,dummy,dummy1,
     .                      rod,ros,.true.,i,j)


cdiag    if (vrbos)
cdiag.     write(*,105)nstep,
cdiag.  ' channel, surf refl dir, surf refl diff',
cdiag.                        (k,rod(k),ros(k),k=1,nlt)
 105  format(i9,a/(18x,i3,2(1x,es9.2)))

cdiag    do k=1,nlt
cdiag    write(*,'(a,4i5,2e12.4)')'obio_model,dir-diff rad:',
cdiag.    nstep,i,j,k,rod(k),ros(k)
cdiag    enddo

         !only call obio_sfcirr for points in light
       tot = 0.0
       if (.not.allocated(eda_frac)) then
         allocate(eda_frac(nlt), esa_frac(nlt))
         call openunit('eda_esa_ratios',iu_bio,.false.,.false.)
         do ichan=1,nlt
           read(iu_bio,'(3f13.8)')dummy1,eda_frac(ichan),esa_frac(ichan)
         enddo
         close(iu_bio)
       endif
       if (allocated(atm%dirvis)) then
         do ichan = 1,nlt
           if (ichan .le. 18) then     !visible + uv
             Ed(ichan) = atm%dirvis(i,j) * eda_frac(ichan)
             Es(ichan) = atm%difvis(i,j) * esa_frac(ichan)
           else               !nir
             Ed(ichan) = atm%dirnir(i,j) * eda_frac(ichan)
             Es(ichan) = atm%difnir(i,j) * esa_frac(ichan)
           endif

           tot = tot + Ed(ichan)+Es(ichan)

       !integrate over all ichan
           OIJ(I,J,IJ_ed) = OIJ(I,J,IJ_ed) + Ed(ichan) ! direct sunlight   
           OIJ(I,J,IJ_es) = OIJ(I,J,IJ_es) + Es(ichan) ! diffuse sunlight   
         enddo  !ichan
       endif
         noon=.false.
         if (hour_of_day.eq.12)then
          if (i.eq.itest.and.j.eq.jtest)noon=.true.
         endif
         if (tot .ge. 0.1) call obio_sfcirr(noon,rod,ros,vrbos)

      if (vrbos) then
        write(*,'(a,3i9)')
     .       'obio_model: counter days,   i,j,nstep=',i,j,nstep
        write(*,'(3(a,i9))')'hour of day=',hour_of_day,
     .                       ', day of month=',day_of_month,
     .                       ', ihr0=',ihr0

cdiag   write(*,'(a)')
cdiag.'    k     dp          u            v         temp         saln'
cdiag   do k=1,kdm
cdiag   write(*,'(i5,5e12.4)')
cdiag   write(*,*)
cdiag.  k,dp1d(k),u(i,j,k+mm),v(i,j,k+mm),
cdiag.    temp(i,j,k+mm),saln(i,j,k+mm)
cdiag   enddo

        write(*,'(a)')
     .'    k     P(1)      P(2)         P(3)       P(4)         P(5) '
        do k=1,kdm
        write(*,'(i5,7e12.4)')
!       write(*,*)
     .   k,obio_P(k,1),obio_P(k,2),obio_P(k,3),obio_P(k,4),
     .   obio_P(k,5),obio_P(k,6),obio_P(k,7)
        enddo

        write(*,'(a)')
     .'    k     P(8)      P(9)         P(11)      P(12)        P(13)'
        do k=1,kdm
        write(*,'(i5,7e12.4)')
!       write(*,*)
     .   k,obio_P(k,8),obio_P(k,9),det(k,1),det(k,2),
     .   det(k,3),car(k,1),car(k,2)
        enddo

        write(*,'(2a)')
cdiag.'    Ed          Es          solz         sunz',
cdiag.'       atmFe       wind'
!       write(*,'(6e12.4)')
cdiag   write(*,*)
cdiag.   Ed(ichan),Es(ichan),solz,sunz,atmFe_ij,wind
      endif

cdiag    if (vrbos)
cdiag.     write(*,106)nstep,
cdiag.    '      channel, dir dwn irr, diff dwn irr,  tot',
cdiag.                  (ichan,Ed(ichan),Es(ichan),
cdiag.                  tot,ichan=1,nlt)
 106  format(i9,a/(18x,i3,3(1x,es9.2)))


         !this part decomposes light into gmao bands-- dont need it
         !         m = indext2   !use the indicator of "past"
         !         call gmaoirr(ihr,nwater,Edg,Esg)
         !         call oasimbands(idmo,ihr,nwater,Edg,Esg)

         do k=1,kdm
          tirrq(k) = 0.0
         enddo

         kpar(:) = 0.d0
         delta_temp1d(:) = 0.d0
         if (tot .ge. 0.1) call obio_edeu(kmax,vrbos,i,j)

#ifdef OBIO_ON_GARYocean
       do k=1,kdm
       OIJL(I,J,k,IJL_kpar) = OIJL(I,J,k,IJL_kpar) + Kpar(k) !  kpar
     .                      * MO(i,j,k) * dxypo(j)    !in order to get landmask-have to set denom_ijl in odiag_com
       OIJL(I,J,k,IJL_dtemp) = OIJL(I,J,k,IJL_dtemp) + delta_temp1d(k) !  change in T due to kpar
     .                      * MO(i,j,k) * dxypo(j)    !in order to get landmask-have to set denom_ijl in odiag_com
       enddo
#endif

         if (vrbos) then
cdiag      write(*,107)nstep,
cdiag.       '           k   avgq    tirrq',
cdiag.                 (k,avgq1d(k),tirrq(k),k=1,kdm)
 107  format(i9,a/(18x,i3,2(1x,es9.2)))

cdiag    do k=1,kdm
cdiag    write(*,'(a,4i5,2e12.5)')'obio_model,k,avgq,tirrq:',
cdiag.       nstep,i,j,k,avgq1d(k),tirrq(k)
cdiag    enddo
         endif

         if (tot .ge. 0.1) ihra_ij = ihra_ij + 1


       !------------------------------------------------------------
       !compute tendency terms on the m level
cdiag     if (vrbos) then
cdiag     write (*,103) nstep,i,j,
cdiag.' bfreptend  dpth       dp       nitr    ',
cdiag.'     ammo     sili     iron',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                       obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     diat       chlo    ',
cdiag.' cyan       cocc     herb',
cdiag.    (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' bfreptend  dpth     nitr_det  sili_det iron_det   doc   ',
cdiag.'    dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif
       !------------------------------------------------------------

cdiag  if (vrbos)write(*,*)'bfre obio_ptend: ',
cdiag.     nstep,(k,tirrq(k),k=1,kmax)

       call obio_ptend(vrbos,kmax,i,j)

       !------------------------------------------------------------
cdiag  if (vrbos)then
cdiag   write(*,108)nstep,' aftrptend dpth      dp     P_tend(1:9)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),P_tend(k,1),P_tend(k,2),
cdiag.     P_tend(k,3),P_tend(k,4),P_tend(k,5),P_tend(k,6),
cdiag.     P_tend(k,7),P_tend(k,8),P_tend(k,9),k=1,kdm)
cdiag   write(*,109)nstep,
cdiag.    ' aftrptend dpth      dp     D_tend(1:3) and C_tend(1:2)',
cdiag.    (k,p(i,j,k+1)/onem,dp1d(k),D_tend(k,1),D_tend(k,2),
cdiag.     D_tend(k,3),C_tend(k,1),C_tend(k,2),k=1,kdm)
cdiag  endif
 108  format(i6,a/(12x,i3,11(1x,es9.2)))
 109  format(i6,a/(12x,i3,7(1x,es9.2)))


       !------------------------------------------------------------

#ifdef OBIO_ON_GARYocean
       !update biology to new time level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl BEFORE update
       call obio_sinksettl(vrbos,kmax,errcon,i,j)
#ifdef noBIO
       do k=1,kmax
       P_tend(k,:)= 0.d0
       C_tend(k,1)= 0.d0     !C_tend(k,2) not set to zero only flux term
       D_tend(k,:)= 0.d0
       obio_ws(k,:) = 0.d0
       wsdet(k,:) = 0.d0
       enddo
#endif
       call obio_update(vrbos,kmax,i,j)
#else
       !update biology from m to n level
       !also do phyto sinking and detrital settling here
       !MUST CALL sinksettl AFTER update
#ifdef noBIO
       do k=1,kmax
       P_tend(k,:)= 0.d0
       C_tend(k,1)= 0.d0     !C_tend(k,2) not set to zero only flux term
       D_tend(k,:)= 0.d0
       obio_ws(k,:) = 0.d0
       wsdet(k,:) = 0.d0
       enddo
#endif
       call obio_update(vrbos,kmax,i,j)
       call obio_sinksettl(vrbos,kmax,errcon,i,j)
       if (errcon) then
          write (*,'(a,3i5)') 'error update at i,j =',i,j,kmax
          do k=1,kdm
          write (*,'(i5,3e12.4)')
     .           k,dp1d(k),p1d(k),wsdet(k,1)
          enddo
          stop
       endif
#endif

cdiag     if (vrbos) then
cdiag     write (*,*)'     '
cdiag     write (*,103) nstep,i,j,
cdiag.' aftrupdate dpth     dp        nitr      ammo      sili',
cdiag.'      iron',
cdiag.  (k,p(i,j,k+1)/onem,dp1d(k),obio_P(k,1),obio_P(k,2),
cdiag.                     obio_P(k,3),obio_P(k,4),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     diat      chlo        cyan',
cdiag.'      cocc      herb',
cdiag.  (k,p(i,j,k+1)/onem,obio_P(k,5),obio_P(k,6),obio_P(k,7),
cdiag.                       obio_P(k,8),obio_P(k,9),k=1,kdm)

cdiag     write (*,104) nstep,i,j,
cdiag.' aftrupdate dpth     nitr_det  sili_det    iron_det',
cdiag.'      doc       dic ',
cdiag.    (k,p(i,j,k+1)/onem,det(k,1),det(k,2),det(k,3),
cdiag.                       car(k,1),car(k,2),k=1,kdm)
cdiag     endif


      !------------------------------------------------------------
      !integrate rhs vertically 
      do ll= 1, 17
      do nt= 1, ntrac-1
      rhs_obio(i,j,nt,ll) = 0.d0
      do k = 1, kdm
          rhs_obio(i,j,nt,ll) = rhs_obio(i,j,nt,ll) +
     .                 rhs(k,nt,ll)*dp1d(k)    
      enddo  !k

      if (vrbos) then
      write(*,'(a,5i5,1x,e20.13)')'rhs_obio (mass,trac/m2/hr):',
     .   nstep,i,j,nt,ll,rhs_obio(i,j,nt,ll)
      endif
      enddo  !ntrac

      !convert all to mili-mol,C/m2
      rhs_obio(i,j,5,ll) = rhs_obio(i,j,5,ll) * mgchltouMC
      rhs_obio(i,j,6,ll) = rhs_obio(i,j,6,ll) * mgchltouMC
      rhs_obio(i,j,7,ll) = rhs_obio(i,j,7,ll) * mgchltouMC
      rhs_obio(i,j,8,ll) = rhs_obio(i,j,8,ll) * mgchltouMC
      rhs_obio(i,j,9,ll) = rhs_obio(i,j,9,ll) * mgchltouMC
      rhs_obio(i,j,10,ll) = rhs_obio(i,j,10,ll) / uMtomgm3

      chng_by(i,j,1) = (rhs_obio(i,j,5,9)+rhs_obio(i,j,6,9)
     .                  +rhs_obio(i,j,7,9)+rhs_obio(i,j,8,9))
     .               +  rhs_obio(i,j,9,9)
      chng_by(i,j,2) = (rhs_obio(i,j,5,5)+rhs_obio(i,j,6,6)
     .               +  rhs_obio(i,j,7,7)+rhs_obio(i,j,8,8))
     .               +  rhs_obio(i,j,10,5)
      chng_by(i,j,3) = (rhs_obio(i,j,5,13)+rhs_obio(i,j,6,13)
     .               +  rhs_obio(i,j,7,13)+rhs_obio(i,j,8,13))
     .               + rhs_obio(i,j,14,6)
      chng_by(i,j,4) = (rhs_obio(i,j,5,14)+rhs_obio(i,j,6,14)
     .               +  rhs_obio(i,j,7,14)+rhs_obio(i,j,8,14))
     .               +  rhs_obio(i,j,13,5)
      chng_by(i,j,5) = (rhs_obio(i,j,5,15)+rhs_obio(i,j,6,15)
     .               +  rhs_obio(i,j,7,15)+rhs_obio(i,j,8,15))
     .               +  rhs_obio(i,j,14,5)
      chng_by(i,j,6) =  rhs_obio(i,j,9,10) 
     .                            + rhs_obio(i,j,13,9)
      chng_by(i,j,7) =rhs_obio(i,j,9,11)
     .                            + rhs_obio(i,j,10,7) 
      chng_by(i,j,8) = rhs_obio(i,j,9,12)
     .                            + rhs_obio(i,j,10,8)
      chng_by(i,j,9) = rhs_obio(i,j,9,14)
     .                            + rhs_obio(i,j,13,12) 
      chng_by(i,j,10)= rhs_obio(i,j,9,15)
     .                            + rhs_obio(i,j,14,15)
      chng_by(i,j,11)= rhs_obio(i,j,10,14)
     .                            + rhs_obio(i,j,13,13)
      chng_by(i,j,12)= rhs_obio(i,j,10,9)
     .                            + rhs_obio(i,j,13,10)
      chng_by(i,j,13)= rhs_obio(i,j,10,10)
     .                            + rhs_obio(i,j,14,10)
      chng_by(i,j,14)= rhs_obio(i,j,13,14)
     .                            + rhs_obio(i,j,14,14)
 
      enddo  !ll

      call obio_chkbalances(vrbos,nstep,i,j)

      do nt=1,ntrac-1   ! don't include unused inert tracer
      do ll=1,17
      OIJ(I,J,IJ_rhs(nt,ll)) = OIJ(I,J,IJ_rhs(nt,ll))
     .                                    + rhs_obio(i,j,nt,ll)  ! all terms in rhs
      enddo
      enddo

      if (vrbos) then
       print*, 'OBIO TENDENCIES, 1-17, 1,7'
       do k=1,1
        write(*,'(17(e9.2,1x))')((rhs(k,nt,ll),ll=1,17),nt=1,7)
         print*, 'OBIO TENDENCIES, 1-17, 8,14'
        write(*,'(17(e9.2,1x))')((rhs(k,nt,ll),ll=1,17),nt=8,ntrac-1)
       enddo
      endif


!      if (i.eq.169 .and. j.eq.59) then  ! test loc at nile mouth
!        write(*,'(a,17(e12.4,1x))') 'dic rhs = ',(rhs(1,14,ll),ll=1,17)
!        write(*,'(a,17(e12.4,1x))') 'dic rhs_obio='
!     .                   ,(rhs_obio(i,j,14,ll),ll=1,17)
!        write(*,*) 'OIJ(I,J,IJ_rhs(14,17))=',OIJ(I,J,IJ_rhs(14,17))
!      endif



!      write(*,'(a,3i5,16(e12.4,1x))')'obio_model, ironrhs:',
!     .  nstep,i,j,(rhs_obio(i,j,4,ll),ll=1,16)

       !------------------------------------------------------------
      !!if(vrbos) write(*,'(a,15e12.4)')'obio_model, strac conc2:',
!      if(vrbos) write(*,*)'obio_model, strac conc2:',
!     .    obio_P(1,:),det(1,:),car(1,:)

       !update 3d tracer array
       do k=1,kmax
        do nt=1,ntyp+n_inert
         tracer(i,j,k,nt)=obio_P(k,nt)
        enddo
        do nt=1,ndet
         tracer(i,j,k,ntyp+n_inert+nt)=det(k,nt)
        enddo
        do nt=1,ncar
         tracer(i,j,k,ntyp+n_inert+ndet+nt)=car(k,nt)
        enddo
#ifdef TRACERS_Alkalinity
         nt=1
         tracer(i,j,k,ntyp+n_inert+ndet+ncar+nt)=alk1d(k)
#ifdef TOPAZ_params
         nt=2
         tracer(i,j,k,ntyp+n_inert+ndet+ncar+nt)=ca_det_calc1d(k)
#endif
#endif
        !update avgq and gcmax arrays
        avgq(i,j,k)=avgq1d(k)
        OIJL(I,J,k,IJL_avgq)= OIJL(I,J,k,IJL_avgq) + avgq1d(k)
        gcmax(i,j,k)=gcmax1d(k)
        tirrq3d(i,j,k)=tirrq(k)
       enddo !k

#ifdef OBIO_ON_GARYocean
      !update trmo etc arrays
       do k=1,kmax
       do nt=1,ntrac
          dtr = tracer(i,j,k,nt) * trmo_unit_factor(k,nt) 
     .        - trmo(i,j,k,nt)

        if (dtr.lt.0) then
          ftr = -dtr/trmo(i,j,k,nt)
          txmo(i,j,k,nt)=txmo(i,j,k,nt)*(1.-ftr)
          tymo(i,j,k,nt)=tymo(i,j,k,nt)*(1.-ftr)
          tzmo(i,j,k,nt)=tzmo(i,j,k,nt)*(1.-ftr)       
        endif

        trmo(i,j,k,nt)=trmo(i,j,k,nt)+dtr

       enddo
       enddo
#else
       tracer_h(i,j,:,:)=tracer(i,j,:,:)*trmo_unit_factor
#endif

       ihra(i,j)=ihra_ij

       !compute total chlorophyl at surface layer
       tot_chlo(i,j)=0.
       do nt=1,nchl
          tot_chlo(i,j)=tot_chlo(i,j)+obio_P(1,nnut+nt)
       enddo
       atm%chl(i,j) = tot_chlo(i,j)
       if (vrbos) then
          !!!write(*,'(/,a,3i5,e12.4)')
          write(*,*)
     .       'obio_model, tot_chlo= ',nstep,i,j,tot_chlo(i,j)
       endif

       !compute total primary production per day
       !sum pp over species and depth, NOT over day
       pp2tot_day(i,j)=0.
       pp2diat_day(i,j)=0.
       pp2chlo_day(i,j)=0.
       pp2cyan_day(i,j)=0.
       pp2cocc_day(i,j)=0.
!      if (p1d(kdm+1).gt.200.) then    !if total depth > 200m
          do nt=1,nchl
          do k=1,kdm
          if (nt.eq.1)pp2diat_day(i,j)
     .       =pp2diat_day(i,j)+pp2_1d(k,1)*24.d0  !mg,C/m2/day
          if (nt.eq.2)pp2chlo_day(i,j)
     .       =pp2chlo_day(i,j)+pp2_1d(k,2)*24.d0  !mg,C/m2/day
          if (nt.eq.3)pp2cyan_day(i,j)
     .       =pp2cyan_day(i,j)+pp2_1d(k,3)*24.d0  !mg,C/m2/day
          if (nt.eq.4)pp2cocc_day(i,j)
     .       =pp2cocc_day(i,j)+pp2_1d(k,4)*24.d0  !mg,C/m2/day

              pp2tot_day(i,j)=pp2tot_day(i,j)+pp2_1d(k,nt)     !mg,C/m2/hr
     &                                       * HOURS_PER_DAY   !->mg,C/m2/day
          enddo
          enddo
!       endif
       if (vrbos) then
       write(*,'(a,i8,2i5,7e12.4)')'obio_model, pp:',
     .    nstep,i,j,tot_chlo(i,j),tot,pp2diat_day(i,j),
     .    pp2chlo_day(i,j),pp2cyan_day(i,j),pp2cocc_day(i,j),
     .    pp2tot_day(i,j)
       endif

       atm%gtracer(atm%n_co2n, i,j)=pCO2_ij

!diagnostics
       if (solz.gt.0) then
           OIJ(I,J,IJ_dayl) = OIJ(I,J,IJ_dayl) + 1.d0   !number of timesteps daylight   
       endif
       OIJ(I,J,IJ_solz) = OIJ(I,J,IJ_solz) + solz  !cos zenith angle
       OIJ(I,J,IJ_sunz) = OIJ(I,J,IJ_sunz) + sunz  !solar zenith angle in degrees  

       OIJ(I,J,IJ_nitr) = OIJ(I,J,IJ_nitr) + tracer(i,j,1,1) ! surf ocean nitrates
       OIJ(I,J,IJ_amm) = OIJ(I,J,IJ_amm) + tracer(i,j,1,2) ! surf ocean nitrates
       OIJ(I,J,IJ_sil) = OIJ(I,J,IJ_sil) + tracer(i,j,1,3) ! surf ocean nitrates
       OIJ(I,J,IJ_iron) = OIJ(I,J,IJ_iron) + tracer(i,j,1,4) ! surf ocean nitrates

       OIJ(I,J,IJ_diat) = OIJ(I,J,IJ_diat) + tracer(i,j,1,5) ! surf ocean diatoms
       OIJ(I,J,IJ_chlo) = OIJ(I,J,IJ_chlo) + tracer(i,j,1,6) ! surf ocean chlorophytes
       OIJ(I,J,IJ_cyan) = OIJ(I,J,IJ_cyan) + tracer(i,j,1,7) ! surf ocean cyanobacteria
       OIJ(I,J,IJ_cocc) = OIJ(I,J,IJ_cocc) + tracer(i,j,1,8) ! surf ocean coccolithophores
       OIJ(I,J,IJ_herb) = OIJ(I,J,IJ_herb) + tracer(i,j,1,9) ! surf ocean herbivores
       OIJ(I,J,IJ_pp) = OIJ(I,J,IJ_pp) + pp2tot_day(i,j)     ! surf ocean pp

       OIJ(I,J,IJ_pp1) = OIJ(I,J,IJ_pp1) + pp2diat_day(i,j)     ! ocean pp from diatoms
       OIJ(I,J,IJ_pp2) = OIJ(I,J,IJ_pp2) + pp2chlo_day(i,j)     ! ocean pp from chlorophytes
       OIJ(I,J,IJ_pp3) = OIJ(I,J,IJ_pp3) + pp2cyan_day(i,j)     ! ocean pp from cyanobacteria
       OIJ(I,J,IJ_pp4) = OIJ(I,J,IJ_pp4) + pp2cocc_day(i,j)     ! ocean pp from coccolithophores

       OIJ(I,J,IJ_doc) = OIJ(I,J,IJ_doc) + tracer(i,j,1,14) ! surf ocean doc
       OIJ(I,J,IJ_dic) = OIJ(I,J,IJ_dic) + tracer(i,j,1,15) ! surf ocean dic
       OIJ(I,J,IJ_pCO2) = OIJ(I,J,IJ_pCO2) + pCO2_ij*(1.-oice(i,j)) ! surf ocean pco2

       OIJ(I,J,IJ_cexp) = OIJ(I,J,IJ_cexp) + cexp             ! export production
       OIJ(I,J,IJ_ndet) = OIJ(I,J,IJ_ndet) + tracer(i,j,4,11) ! ndet at 74m
       OIJ(I,J,IJ_wsd)= OIJ(I,J,IJ_wsd)+ wsdet(4,1)           ! sink. vel. n/cdet at 74m
       if (4<=kmax) then
         OIJ(I,J,IJ_xchl) = OIJ(I,J,IJ_xchl)
     .              + tracer(i,j,4,5)*trmo_unit_factor(4,5)*obio_ws(4,1)
     .              + tracer(i,j,4,6)*trmo_unit_factor(4,6)*obio_ws(4,1)
     .              + tracer(i,j,4,7)*trmo_unit_factor(4,7)*obio_ws(4,1)
     .              + tracer(i,j,4,8)*trmo_unit_factor(4,8)*obio_ws(4,1) !total phyto cexp at 74m
       else
         oij(i,j,ij_xchl)=0
       endif

       !limitation diags surface only (for now)
       k = 1
       do nt=1,nchl
       do ilim=1,5
       OIJ(I,J,IJ_lim(nt,ilim)) = OIJ(I,J,IJ_lim(nt,ilim)) 
     .                          + flimit(k,nt,ilim) 
       enddo
       enddo

       OIJ(I,J,IJ_flux) = OIJ(I,J,IJ_flux) + co2flux      !air-sea CO2 flux(watson)

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!       OIJ(I,J,IJ_rnitrmflo) = OIJ(I,J,IJ_rnitrmflo)+rnitrmflo_loc(i,j)  ! riverine nitr mass flow from dC (kg/s)
       OIJ(I,J,IJ_rnitrconc) = OIJ(I,J,IJ_rnitrconc)+rnitrconc_loc(i,j)  ! riverine nitr conc from dC (kg/kg)
#endif
#ifdef DIC_RUNOFF
       OIJ(I,J,IJ_rdicconc) = OIJ(I,J,IJ_rdicconc)+rdicconc_loc(i,j)     ! riverine dic conc from dC (kg/kg)
#endif
#ifdef DOC_RUNOFF
       OIJ(I,J,IJ_rdocconc) = OIJ(I,J,IJ_rdocconc)+rdocconc_loc(i,j)     ! riverine doc conc from dC (kg/kg)
#endif
#ifdef SILI_RUNOFF
       OIJ(I,J,IJ_rsiliconc) = OIJ(I,J,IJ_rsiliconc)+rsiliconc_loc(i,j)  ! riverine silica conc from dC (kg/kg)
#endif
#ifdef IRON_RUNOFF
       OIJ(I,J,IJ_rironconc) = OIJ(I,J,IJ_rironconc)+rironconc_loc(i,j)  ! riverine iron conc from dC (kg/kg)
#endif
#ifdef POC_RUNOFF
       OIJ(I,J,IJ_rpocconc) = OIJ(I,J,IJ_rpocconc)+rpocconc_loc(i,j)     ! riverine poc conc from dC (kg/kg)
#endif
#ifdef ALK_RUNOFF
       OIJ(I,J,IJ_ralkconc) = OIJ(I,J,IJ_ralkconc)+ralkconc_loc(i,j)     ! riverine alkalinity conc from A-S (mol/kg)
#endif
#endif

       if (tracers_alkalinity) then
         OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) + tracer(i,j,1,16)    ! surf ocean alkalinity
         OIJ(I,J,IJ_fca) = OIJ(I,J,IJ_fca) + caexp               ! carbonate export
       else
         OIJ(I,J,IJ_alk) = OIJ(I,J,IJ_alk) + alk(i,j,1)          ! surf ocean alkalinity
       endif

#ifndef OBIO_ON_GARYocean    /* HYCOM ACCUMULATED DIAGNOSTICS */
      ao_co2fluxav_loc(i,j)=ao_co2fluxav_loc(i,j) + co2flux
      pCO2av_loc(i,j)=pCO2av_loc(i,j)+pCO2_ij
      pp2tot_dayav_loc(i,j) = pp2tot_dayav_loc(i,j) 
     .                      + pp2tot_day(i,j)
      pp2diat_dayav_loc(i,j) = pp2diat_dayav_loc(i,j)
     .                      + pp2diat_day(i,j)
      pp2chlo_dayav_loc(i,j) = pp2chlo_dayav_loc(i,j)
     .                      + pp2chlo_day(i,j)
      pp2cyan_dayav_loc(i,j) = pp2cyan_dayav_loc(i,j)
     .                      + pp2cyan_day(i,j)
      pp2cocc_dayav_loc(i,j) = pp2cocc_dayav_loc(i,j)
     .                      + pp2cocc_day(i,j)
      cexpav_loc(i,j)=cexpav_loc(i,j)+cexp
#ifdef TRACERS_Alkalinity
            caexpav_loc(i,j)=caexpav_loc(i,j)+caexp
#endif
#endif

      end do
      end do
      call stop('  obio main loop')

      call stop(' obio_model')
      initialized=.true.
      nstep0=nstep0+1
      return

      end subroutine obio_model
