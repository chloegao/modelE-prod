!     SCMDATA_TWPICE.f
!@sum program to read in SCM forcing data and get it in format for passing
!     to model
!@auth  Audrey Wolf
!@ver  Version to read TWP ICE  IOP                          
!
C--------------------------------------------------------------------------------

      SUBROUTINE init_scmdata

!     read data from file store time history and set up initial conditions
      USE SCMCOM , only : SCM_SURFACE_FLAG,NARM,TAUARM,NRINIT,IKT, 
     &                  AMEANPS, SG_T, SG_Q,
     &                  SG_U,SG_V,ASWINDSPD,AQS,AVS,AUS,ATSAIR,ATSKIN,
     &                  iu_scm_prt,NSTEPSCM
      USE RESOLUTION , only : LM,LS1,PTOP,PSF
      USE ATM_COM , only : P,PK,T,Q,U,V
      USE FLUXES, only : FLAND,FOCEAN,FLICE,FLAKE0,FEARTH0
      USE GHY_COM, only : FEARTH
      USE LAKES_COM, only : FLAKE
      USE CONSTANT , only : KAPA,TF   
      USE FLUXES, only : atmocn,atmlnd,atmsrf
      USE FILEMANAGER
      USE Dictionary_mod
      IMPLICIT NONE

      INTEGER L,I,J       

      call openunit("scm.prt",iu_scm_prt,.false.,.false.)
      write(0,*) 'I/J Targets set ',1,1
      write(iu_scm_prt,*) 'I/J Targets set ',1,1

      if (SCM_SURFACE_FLAG.eq.0) then
          write(0,*) 
     &        'SCM set to run with GCM calculated surface fluxes'
          write(iu_scm_prt,*) 
     &        'SCM set to run with GCM calculated surface fluxes'
      else if (SCM_SURFACE_FLAG.eq.1) then
          write(0,*) 
     &        'SCM set to run with ARM prescribed surface fluxes'
          write(iu_scm_prt,*) 
     &        'SCM set to run with ARM prescribed surface fluxes'
      else if (SCM_SURFACE_FLAG.eq.2) then
          write(0,*) 
     &        'SMC set to run with ARM srf tmps + GCM calc srf fluxes'
          write(iu_scm_prt,*)
     &        'SCM set to run with ARM srf tmps + GCM calc srf fluxes'
      endif
      NARM = 6                        ! for 30 min time steps from 3-hourly input 
      TAUARM = 0                      ! not used for now 
      NRINIT = 0                     ! when to reinitialize t,q profiles to data
      IKT = 1                         ! index to data
 

      write(iu_scm_prt,25) FLAND(1,1),
     &   FOCEAN(1,1),FLICE(1,1),
     &   FLAKE0(1,1),
     &   FEARTH0(1,1),FEARTH(1,1)
 25   format(1x,'init flags  land ocean lice lake earth0 earth ',
     &   6(f8.3))

      call init_read_surface
      call init_read_layers

c     initialize variables from SCM data
      call pass_scm_surface
      call pass_scm_layers


c     get I_TARG and J_TARG from MODEL_TARGET_COM.f  as included in MODEL_COM 
c     module   
      P(1,1) = AMEANPS - PTOP     

      call CALC_AMPK(LM)
      
      do L = 1,LM
c        fill T with potential temperature as ARM provided
         Q(1,1,L) = SG_Q(L)
         T(1,1,L) = SG_T(L)/PK(L,1,1)
      enddo 

c     also note for U,V initialize boxes around target I,J for mean wind
c     calculations in ATURB 
c     check at some point into passing mean surface winds.

      
      do L=1,LM
c        write(iu_scm_prt,'(a10,i3,i5,i5,2(f9.3))') 
c    &       'L J I u v ',l,1,1,sg_u(L),sg_v(l)
         U(1,1,L) = SG_U(L)
         V(1,1,L) = SG_V(L)
      enddo 
c
cccc
c     set surface variables - - - find new variable names
c 
      atmsrf%WSAVG(1,1)  = ASWINDSPD        !BLDATA(1)  
      atmsrf%USAVG(1,1) = AUS
      atmsrf%VSAVG(1,1) = AVS

      atmlnd%GTEMP(1,1) = ATSKIN        !GDATA(4)    
      atmocn%GTEMP(1,1) = ATSKIN 
      atmocn%GTEMP2(1,1) = ATSKIN 
      atmocn%GTEMPR(1,1) = ATSKIN + TF
      atmlnd%GTEMPR(1,1) = ATSKIN + TF

      write(iu_scm_prt,120) atmocn%GTEMP(1,1),
     &      atmocn%GTEMP2(1,1),atmlnd%GTEMP(1,1),
     &      atmocn%GTEMPR(1,1),atmlnd%GTEMPR(1,1)
 120  format(1x,
     &    'initial temps  gtemp11 gtemp21 gtemp14 gtempr1 gtempr4 ',
     &     5(f10.3))

      end subroutine init_scmdata



      SUBROUTINE init_read_surface
C    
      USE FILEMANAGER, only : openunit,closeunit
      USE SCMCOM , only : MCT,IKT,NARM,stmstep,AMPS,ATSA,ATSK,
     &                    APRCHR,ALHHR,ASHHR,AWSHR,AVSHR,AUSHR,
     &                    AQSHR,ALWPHR,ARHHR,iu_scm_prt  

      IMPLICIT NONE
      INTEGER NTARM,NVARSRF 
      parameter(NTARM=199,NVARSRF=43)  !TWPICE IOP
      REAL*4 d(NVARSRF,NTARM), rdat(NTARM,NVARSRF)
      character*50 var_name(NVARSRF)    


      INTEGER II,K,KT,ITP,ITIM
      INTEGER iu,iv,NTOT     
      REAL*4 DT,DP,RARM
      REAL*4 PSMEAN
      character*12 scmdummy

      integer itest
  
 
c      ARM VARIABLEs
c      data var_name/
c     +     'Calday',          'Year',                  'Month',  
c     +     'Day',             'Hour',                  'Minute', 
c     +     'Prec[Surface Precipitation] (mm/hour)',
c     +     'LH[Surf Latent Heat Flux, upward positive] (W/m2)',
c     +     'SH[Surf Sensible Heat Flux,upward positive] (W/m2)', 
c     +     'PSA[Surf pressure averaged over the domain](mb)',
c     +     'PSI[Surf Pressure at center of domain(mb)', 
c     +     'Ts_Air[2m air temperature(C)', 
c     +     'Tskin [Surf skin temperature (C)',
c     +     'RHair[2m air relative humidity (%)',  
c     +     'WSPD[10m wind speed](m/s)', 
c     +     'u_wind[10m u component(m/s)',
c     +     'v_wind[10m v component(m/s)',
c     +     'Srf_Net_Dn_Rad[surf net radiatin downward positive(W/m2)', 
c     +     'FLUT [TOA_LW_Up_positive](W/m2)', 
c     +     'FSNT[TOA_net_SW_Down_positive] (W/m2)', 
c     +     'SOLIN[TOA_Insolation(W/m2)',       
c     +     'SAT_Lowcld(%)',  
c     +     'SAT_Midcld(%) ',
c     +     'SAT_Hghcld(%)',     
c     +     'SAT_Totcld(%)',  
c     +     'Cld_Thickness(km)', 
c     +     'Cld_Top_ht(km)',   
c     +     'LWP[MWR_Cld_liquid water path](cm)', 
c     +     'CDH2ODT[column integrated dH2O/dt](mm/hour)',       
c     +     'CH2OADV[column integrated H2O_Advection](mm/hour)',
c     +     'Evap[Srf_Evaporation](mm/hour)',
c     +     'CDSDT[column d(dry static energy)/dt](W/m2) ' 
c     +     'CSADV[Column_Dry_Static_Energy_Advection](W/m2)', 
c     +     'CRAD[Column_Radiative_Heating](W/m2)',        
c     +     'CLH[Column_Latent_heating](W/m2)' 
c     +     'OMEGAS[omega_surface](mb/hr)',
c     +     'qs[wm water vapor mixing ratio](kg/kg)',
c     +     'S[2m dry static energy](K)',
c     +     'PW[MWR column precip_water](cm)',
c     +     'FLUS[surface upwelling LW](W/m2)',
c     +     'FLDS[surface downwelling LW](W/m2)',
c     +     'FSUS[surface upwelling SW](W/m2)',
c     +      FSDS[surface downwelling SW](W/m2)' 

      itest = 0

      KT = IKT
      call openunit('SCMSRF',iu,.false.,.true.)    

      read(iu,*) scmdummy 
      read(iu,*) scmdummy 
      read(iu,*) scmdummy 
      read(iu,*) scmdummy 


      do iv=1,NVARSRF
         read(iu,111)var_name(iv)
111      format(a50)
         read(iu,112)(d(iv,k),k=1,NTARM)
112      format(5e15.7)
      enddo
      call closeunit(iu)  

c     do itim=1,NTARM
c        write(iu_scm_prt,115) d(2,itim),d(3,itim),d(4,itim),d(5,itim),
c    &               d(10,itim),d(13,itim)    
c115     format(1x,f6.0,3(f3.0),f8.2,f8.3)
c     enddo

C
C     fill array with surface pressure by hour
c 
      RARM = NARM
      itim=1
      stmstep(itim) = d(1,1)
      do itp = 1,NTARM-1
         dt = d(1,itp+1)-d(1,itp)
         dt = dt/RARM
         do ii=1,NARM
            itim = itim+1
            stmstep(itim) = stmstep(itim-1)+dt
         enddo
      enddo

      itim = 1
      AMPS(itim)=d(10,1)
      do itp = 1,NTARM-1
         dp = d(10,itp+1)-d(10,itp)
         dp = dp/RARM
         do ii = 1,NARM
            itim = itim+1
            AMPS(itim) = AMPS(itim-1)+dp
         enddo
      enddo 

      PSMEAN=0.0
      NTOT=itim
      do itim=1,NTOT
         PSMEAN = PSMEAN + AMPS(itim)
      enddo
      PSMEAN=PSMEAN/NTOT
c     write(iu_scm_prt,120) NTOT,PSMEAN
c120  format(1x,'NTOT   PSMEAN ',i5,f10.2)

      itim = 1
      ATSK(itim)=d(13,1)
      do itp = 1,NTARM-1
         dt = d(13,itp+1)-d(13,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            ATSK(itim) = ATSK(itim-1)+dt
         enddo
      enddo 

      itim = 1
      ATSA(itim)=d(12,1)
      do itp = 1,NTARM-1
         dt = d(12,itp+1)-d(12,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            ATSA(itim) = ATSA(itim-1)+dt
         enddo
      enddo 


      itim = 1
      APRCHR(itim)=d(7,1)
      do itp = 1,NTARM-1
         dt = d(7,itp+1)-d(7,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            APRCHR(itim) = APRCHR(itim-1)+dt
         enddo
      enddo 

      itim = 1
      ALHHR(itim)=d(8,1)
      do itp = 1,NTARM-1
         dt = d(8,itp+1)-d(8,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            ALHHR(itim) = ALHHR(itim-1)+dt
         enddo
      enddo 

      itim = 1
      ARHHR(itim) = d(14,1)
      do itp = 1,NTARM-1 
         dt = d(14,itp+1)-d(14,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim + 1
            ARHHR(itim) = ARHHR(itim-1) + dt
         enddo
      enddo

      itim = 1
      ASHHR(itim)=d(9,1)
      do itp = 1,NTARM-1
         dt = d(9,itp+1)-d(9,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            ASHHR(itim) = ASHHR(itim-1)+dt
         enddo
      enddo

      itim = 1
      AWSHR(itim)=d(15,1)
      do itp = 1,NTARM-1
         dt = d(15,itp+1)-d(15,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            AWSHR(itim) = AWSHR(itim-1)+dt
         enddo
      enddo 

      itim = 1
      AUSHR(itim)=d(16,1)
      do itp = 1,NTARM-1
         dt = d(16,itp+1)-d(16,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            AUSHR(itim) = AUSHR(itim-1)+dt
         enddo
      enddo 

      itim = 1
      AVSHR(itim)=d(17,1)
      do itp = 1,NTARM-1
         dt = d(17,itp+1)-d(17,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            AVSHR(itim) = AVSHR(itim-1)+dt
         enddo
      enddo 

      itim = 1
      AQSHR(itim) = d(37,1)      
      do itp = 1,NTARM-1
         dt = d(37,itp+1)-d(37,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim + 1
            AQSHR(itim) = AQSHR(itim-1) + dt
         enddo
      enddo

      itim = 1
      ALWPHR(itim) = d(28,1)
      do itp = 1,NTARM-1
         dt = d(28,itp+1)-d(28,itp)
         dt = dt/RARM
         do ii = 1,NARM
            itim = itim+1
            ALWPHR(itim) = ALWPHR(itim-1)+dt
         enddo
      enddo

c     itim = 1
c     ADWDTHR(itim) = d(29,1)
c     do itp = 1,NTARM-1
c        dt = d(29,itp+1)-d(29,itp)
c        dt = dt/RARM
c        do ii = 1,NARM
c           itim = itim+1
c           ADWDTHR(itim) = ADWDTHR(itim-1)+dt
c        enddo
c     enddo

c     itim = 1
c     ADWADVHR(itim) = d(30,1)
c     do itp = 1,NTARM-1
c        dt = d(30,itp+1)-d(30,itp)
c        dt = dt/RARM
c        do ii = 1,NARM
c           itim = itim+1
c           ADWADVHR(itim) = ADWADVHR(itim-1)+dt
c        enddo
c     enddo

c     itim = 1
c     ATLWUPHR(itim) = d(19,1)
c     do itp = 1,NTARM-1
c        dt = d(19,itp+1)-d(19,itp)
c        dt = dt/RARM
c        do ii = 1,NARM
c           itim = itim+1
c           ATLWUPHR(itim) = ATLWUPHR(itim-1)+dt
c        enddo
c     enddo

c     itim = 1
c     ATSWDNHR(itim) = d(20,1)
c     do itp = 1,NTARM-1
c        dt = d(20,itp+1)-d(20,itp)
c        dt = dt/RARM
c        do ii = 1,NARM
c           itim = itim+1
c           ATSWDNHR(itim) = ATSWDNHR(itim-1)+dt
c        enddo
c     enddo

c     itim = 1
c     ATSWINHR(itim) = d(21,1)
c     do itp = 1,NTARM-1
c        dt = d(21,itp+1)-d(21,itp)
c        dt = dt/RARM
c        do ii = 1,NARM
c           itim = itim+1
c           ATSWINHR(itim) = ATSWINHR(itim-1)+dt
c        enddo
c     enddo

      return

      end subroutine init_read_surface



      SUBROUTINE pass_scm_surface
        
      USE SCMCOM
      IMPLICIT NONE

      INTEGER KT


      KT = NSTEPSCM + IKT 


      ASTIME = stmstep(KT)
c
c     in this test version apply an averaged LH and SH over this time
c     step by getting an average of time now and time next
      ALH = ALHHR(KT)
      ASH = ASHHR(KT)
      AMEANPS = AMPS(KT)
      APREC = APRCHR(KT)
      ATSAIR = ATSA(KT)
      ATSKIN = ATSK(KT)
      ASWINDSPD = AWSHR(KT)
      AUS = AUSHR(KT)
      AVS = AVSHR(KT)
      AQS = AQSHR(KT)
      ARHSAIR = ARHHR(KT) 
      ALWP = ALWPHR(KT)
      ADWDT = ADWDTHR(KT)
      ADWADV = ADWADVHR(KT)
      ATLWUP = ATLWUPHR(KT)
      ATSWDN = ATSWDNHR(KT)
      ATSWIN = ATSWINHR(KT)

c     write(iu_scm_prt,100) nstepscm,kt,AMEANPS,ATSAIR,ATSKIN     
c100  format(1x,'pass scm variables  nstepscm kt ',
c    &       2(i6),f10.2,f10.4,f10.4)

      return
      end subroutine pass_scm_surface



      SUBROUTINE init_read_layers


      USE FILEMANAGER , only : openunit,closeunit   
      USE SCMCOM

      IMPLICIT NONE
      INTEGER NTARM,NVARLAY,NPARM,NPM
      parameter(NTARM=199,NVARLAY=19, NPARM=40)  !TWPICE IOP
      parameter(NPM=NPARM-1)

      real*4 time(NTARM),      !Calenday day
     +     yy(NTARM),        !Year
     +     mo(NTARM),        !Month
     +     dd(NTARM),        !Day
     +     hh(NTARM),        !Hour
     +     mm(NTARM)         !Minutes

      real*4 parm(NPARM)          !Pressure (mb)
      real*4 d(NVARLAY,NPARM,NTARM) !Multi-Layer fields as described below
c                          -9999.0 is assigned to levels below 
c                          the ground level

      COMMON /CTHREE/ t3hr(NPARM,NTARM),q3hr(NPARM,NTARM),
     &   u3hr(NPARM,NTARM),v3hr(NPARM,NTARM),om3hr(NPARM,NTARM),
     &   wd3hr(NPARM,NTARM),hta3hr(NPARM,NTARM),vsa3hr(NPARM,NTARM),
     &   hqa3hr(NPARM,NTARM),vqa3hr(NPARM,NTARM),acld3hr(NPARM,NTARM),
     &   tm1hr(MCT)
      REAL*4 t3hr,q3hr,u3hr,v3hr,om3hr,wd3hr,hta3hr,vsa3hr,hqa3hr,
     &   vqa3hr,acld3hr,tm1hr

      REAL*4 RARM 
      INTEGER I,J,K,ii,itp,IP   
      INTEGER iu,nv     
      character*50 var_name2(nvarlay)
      real*8 QSC
      parameter (QSC=1000000.0)
      real*8 hrsec,dt  
      integer itest
      character*12 scmdummy

c     data var_name2/
c    +   'Temp[temperature](K)',
c    +   'w[Water Vapor Mixing ratio](g/kg)',  
c    +   'u[horizontal u-wind component ](m/s)',
c    +   'v[horizontal v-wind component ](m/s)',             
c    +   'omega[verticle velocity ]_(mb/hour)',
c    +   'DIV[Horizontal Wind Divergence] (1/s)',   
c    +   'TadvH[Horizontal_Temp__Advec_](K/hour)',
c    +   'TadvV[Vertical_T_Advec](K/hour)',
c    +   'QadvH[Horizontal_q_Advec_](g/kg/hour)', 
c    +   'QadvV[Vertical_q_Advec](g/kg/hour)',
c    +   's(Dry_Static_Energy)(K)',
c    +   'SadvH[Horizontal_s_Advec_](K/hour)',
c    +   'SadvV[Vertical_s_Advec](K/hour)',
c    +   'ds/dt[d(dry static energy)/dt](K/hour)',       
c    +   'DT/dt[d(temperature)/dt](K/hour)',
c    +   'dq/dt[d(water vapor mixing ratio)/dt](g/kg/hour)',   
c    +   'Q1[Apparent Heat Sources Yanai (1973)] (k/hour)',
c    +   'Q2[Apparent Moisture sinks Yanai (1973)] (K/hour)', 
c    +   'CLD[Cloud Fraction](%)'
c
c
      hrsec = 1./3600.
      itest = 0

c   set up input files to read SCM layer data 

      call openunit('SCMLAY',iu,.false.,.true.)    

      read(iu,*,ERR=290,END=295) scmdummy
      read(iu,*,ERR=290,END=295) scmdummy 
      read(iu,*,ERR=290,END=295) scmdummy 
      read(iu,*,ERR=290,END=295) scmdummy 
      if (itest.eq.1) stop 299      
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295) parm
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)time
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)yy
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)mo
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)dd
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)hh
      read(iu,*,ERR=290,END=295)
      read(iu,112,ERR=290,END=295)mm
      read(iu,*,ERR=290,END=295)
      read(iu,('(i8)'),ERR=290,END=295) nv

c     do i=1,NTARM
c        write(iu_scm_prt,90) i,yy(i),mo(i),dd(i),hh(i),mm(i)
c90      format(1x,i5,5(f8.2))
c     enddo
         
      do i=1,nv
         read(iu,111,ERR=290,END=295)var_name2(i)
         do j=1,NPARM
            read(iu,112,ERR=290,END=295)(d(i,j,k),k=1,NTARM)
c           write(iu_scm_prt,95) i,j,d(i,j,1)
c95         format(1x,'read var(layers) var# layer value ',
c    &             i5,i5,f15.5) 
         enddo
      enddo
111   format(a50)
112   format(5e15.7)


      call closeunit(iu)
      

c
C
c     fill 3 hour arrays of full time line and send to interpolation routine
c 
      do itp = 1,NTARM
         do ip = 1,NPARM
            t3hr(ip,itp) = d(1,ip,itp)
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            q3hr(ip,itp) = d(2,ip,itp)/1000.0
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            u3hr(ip,itp) = d(3,ip,itp)
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            v3hr(ip,itp) = d(4,ip,itp)
         enddo
      enddo

      do itp = 1,NTARM 
         do ip = 1,NPARM 
            om3hr(ip,itp) = d(5,ip,itp)
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            wd3hr(ip,itp) = d(6,ip,itp)
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            hta3hr(ip,itp) = d(7,ip,itp)*hrsec
         enddo
      enddo
C     use vertical S advection instead of vertical T advection - see CAUTION
C                   in CASE 3 Intercomparison Instructions
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            vsa3hr(ip,itp) = d(13,ip,itp)*hrsec
         enddo
      enddo
      do itp = 1,NTARM
         do ip = 1,NPARM 
            hqa3hr(ip,itp) = (d(9,ip,itp)/1000.0)*hrsec
         enddo
      enddo
      do itp = 1,NTARM 
         do ip = 1,NPARM 
            vqa3hr(ip,itp) = (d(10,ip,itp)/1000.0)*hrsec
         enddo
      enddo

      do itp = 1,NTARM 
         do ip = 1,NPARM 
            acld3hr(ip,itp) = d(19,ip,itp) 
c           write(iu_scm_prt,101)  itp,ip,acld3hr(ip,itp)   
c101        format(1x,'itp ip  clds ',i5,i5,f10.3)
         enddo
      enddo

      RARM = NARM
      i = 1 
      STMSTEPL(i) = time(1)
      do itp = 1,NTARM-1
         dt = time(itp+1)-time(itp)
         dt = dt/RARM
         do ii = 1,NARM
            i=i+1
            STMSTEPL(i) = STMSTEPL(i-1)+dt
         enddo
      enddo

C     call interpolation routine to interpolate to hourly and then
C     to SCM sigma levels
      call arm_to_sig(parm)

      return

 290  write(iu_scm_prt,*) 'error reading scm layers data ' 
      return
 295  write(iu_scm_prt,*) 'unexpected end of SCM layers data '
      return
     
      end subroutine init_read_layers 




      SUBROUTINE pass_scm_layers 
c         
      USE RESOLUTION , only : LM,LS1,PTOP,PSF 
      USE DYNAMICS, only : SIG
      USE SCMCOM


      IMPLICIT NONE
      INTEGER L,KT


      KT = NSTEPSCM + IKT 
C     fill variables already interpolated to SCM sig levels

      ALTIME = STMSTEPL(KT)
      do L = 1,LS1-1       
         SG_P(L) = SIG(L)*(AMPS(KT)-PTOP) + PTOP
      enddo
      do L=LS1,LM   
         SG_P(L) = SIG(L)*(PSF-PTOP) + PTOP
      enddo
      do L = 1,LM
         SG_T(L) = THR(L,KT)
      enddo
      do L = 1,LM
         SG_Q(L) = QHR(L,KT)
      enddo
      do L = 1,LM
         SG_U(L) = UHR(L,KT)
      enddo
      do L = 1,LM
         SG_V(L) = VHR(L,KT)
      enddo
      do L = 1,LM
         SG_OMEGA(L) = OMGHR(L,KT)
      enddo
      do L = 1,LM
         SG_WINDIV(L) = WDHR(L,KT)
      enddo 
c         
      do L = 1,LM
         SG_HOR_TMP_ADV(L) = HTA_HR(L,KT)
      enddo
      do L = 1,LM
         SG_VER_S_ADV(L) = VSA_HR(L,KT)
      enddo
      do L = 1,LM
         SG_HOR_Q_ADV(L) = HQA_HR(L,KT)
      enddo
      do L = 1,LM
         SG_VER_Q_ADV(L) = VQA_HR(L,KT)
      enddo
c     do l=1,LM
c        write(iu_scm_prt,100) L,SG_HOR_TMP_ADV(L),SG_VER_S_ADV(L),
c    &               SG_HOR_Q_ADV(L),SG_VER_Q_ADV(L)
c100     format(1x,'pass_scm_layers l tadv qadv ',i5,4(E15.5))
c     enddo 


      do L=1,LM
         SG_ARSCL(L) = ACLDHR(L,KT)
      enddo
 
      return

      end subroutine pass_scm_layers 

 
      subroutine pass_SCMDATA

      USE RESOLUTION , only : LM,LS1,PTOP,PSF
      USE ATM_COM , only  : P,PK,T,Q,U,V
      USE FLUXES, only : FLAND,FOCEAN,FLICE,FLAKE0,FEARTH0
      USE GHY_COM, only : FEARTH
      USE LAKES_COM, only : FLAKE
      USE FLUXES, only : atmocn,atmlnd,atmsrf
      USE CONSTANT, only : tf,KAPA 
      USE SCMCOM
C     
C
      INTEGER MODINT
      INTEGER I,J,L 
C                
      MODINT = 9999
      if (NRINIT.gt.0) MODINT = MOD(NSTEPSCM,NRINIT)
      
      write(iu_scm_prt,25) FLAND(1,1),
     &   FOCEAN(1,1),FLICE(1,1),
     &   FLAKE0(1,1),
     &   FEARTH0(1,1),FEARTH(1,1)
 25   format(1x,'pass flags  land ocean lice lake earth0 earth ',
     &   6(f8.3))
c
c     if you want to change the land/water flags --- this is the place to do it
c
c
      call pass_scm_surface 
      call pass_scm_layers 

c * * * * indices
      P(1,1) = AMEANPS - PTOP   
      call CALC_AMPK(LM)
 
 
      if ((ALTIME-ASTIME).gt.0.005) stop 4000 


C     if we are updating profiles of T and Q to ARM data then check
C     if it is time to reinitialize them

      if (MODINT.eq.0) then
c         do L = 1,LM
c            write(iu_scm_prt,300) L,Q(1,1,L),
c    &                  T(1,1,L)
c300         format(1x,'pass_SCMDATA OLD L Q T ',i5,E10.4,f8.2)
c         enddo
          do L = 1,LM
             Q(1,1,L) = SG_Q(L)
C            get potential temperature 
C* * * * check how to do this now
             T(1,1,L) = SG_T(L) / PK(L,1,1)    
c            write(iu_scm_prt,310) L,Q(1,1,L),SG_T(L),
c    &                     T(1,1,L)
c310         format(1x,'NEW ICS  L Q SGT T ',i5,E10.4,f9.2,f8.2)
          enddo 
          atmlnd%GTEMP(1,1) = ATSKIN !GDATA(4)    
          atmocn%GTEMP(1,1) = ATSKIN 
          atmocn%GTEMP2(1,1) = ATSKIN 
          atmocn%GTEMPR(1,1) = ATSKIN + TF
          atmlnd%GTEMPR(1,1) = ATSKIN + TF
          write(iu_scm_prt,340) atmlnd%GTEMP(1,1),
     &       atmocn%GTEMP(1,1),atmocn%GTEMP2(1,1),
     &       atmocn%GTEMPR(1,1)
 340      format(1x,'SCM GTEMP14 GTEMP11 GTEMP21 GTEMPR1 ',4(f10.3))
      endif

      do L=1,LM
         U(1,1,L) = SG_U(L)
         V(1,1,L) = SG_V(L)
      enddo 
c
c     set surface variables 
c    

      atmsrf%WSAVG(1,1)  = ASWINDSPD        !BLDATA(1)
      atmsrf%USAVG(1,1) = AUS
      atmsrf%VSAVG(1,1) = AVS
      if (SCM_SURFACE_FLAG.eq.1) then
          atmlnd%GTEMP(1,1) = ATSKIN !GDATA(4)    
          atmocn%GTEMP(1,1) = ATSKIN 
          atmocn%GTEMP2(1,1) = ATSKIN 
          atmocn%GTEMPR(1,1) = ATSKIN + TF
          atmlnd%GTEMPR(1,1) = ATSKIN + TF
          write(iu_scm_prt,360) atmlnd%GTEMP(1,1),
     &       atmocn%GTEMP(1,1),atmocn%GTEMP2(1,1),
     &       atmocn%GTEMPR(1,1),atmlnd%GTEMPR(1,1)
 360      format(1x,
     &      'pass SCM GTEMP14 GTEMP11 GTEMP21 GTEMPR1 GTEMPR4',
     &        5(f10.3))

      endif
 
      return 

      end subroutine pass_SCMDATA 


C--------------------------------------------------------------------
C
C     routine to interpolate arm data to GCM sigma levels, weighted by
C     pressure
C
      SUBROUTINE arm_to_sig(parm)

      USE RESOLUTION , only : LM,LS1,PTOP,PSF
      USE DYNAMICS, only : SIG,SIGE
      USE SCMCOM  
      USE CONSTANT , only : grav
C
      IMPLICIT NONE
      INTEGER NTARM,NPARM 
      parameter(NTARM=199,NPARM=40)  
      
      COMMON /CTHREE/ t3hr(NPARM,NTARM),q3hr(NPARM,NTARM),
     &     u3hr(NPARM,NTARM),v3hr(NPARM,NTARM),om3hr(NPARM,NTARM),
     &     wd3hr(NPARM,NTARM),hta3hr(NPARM,NTARM),vsa3hr(NPARM,NTARM),
     +     hqa3hr(NPARM,NTARM),vqa3hr(NPARM,NTARM),acld3hr(NPARM,NTARM),
     &     tm1hr(MCT)
      real*4 t3hr,q3hr,u3hr,v3hr,om3hr,wd3hr,hta3hr,vsa3hr,hqa3hr,
     +       vqa3hr,acld3hr,tm1hr
C     pass three hour data and interpolate to 1 hr data


      real*4 parm(NPARM)

      real*4 t1hr(NPARM,MCT),q1hr(NPARM,MCT),u1hr(NPARM,MCT),
     &       v1hr(NPARM,MCT),om1hr(NPARM,MCT),wd1hr(NPARM,MCT),
     &       hta1hr(NPARM,MCT),vsa1hr(NPARM,MCT),
     &       hqa1hr(NPARM,MCT),vqa1hr(NPARM,MCT),acld1hr(NPARM,MCT)
   
      real*8 APE(NPARM+1),AP(NPARM)
      real*8 DELTAP,DELAG,DELAG1,DELAG2 
      real*8 QSC
      parameter (QSC=1000000.0)
      real*8 totpa,totpg
      real*8 sumaqh,sumaqv,sumath,sumasv,
     +       sumat,sumaq,sumau,sumav,sumawd,sumaom,sumac1
      real*8 sumgqh,sumgqv,sumgth,sumgsv,sumgt,sumgq,sumgu,sumgv,
     +       sumgwd,sumgom,sumgc1


      real*4 RARM,dx       
      INTEGER L,n,ni,ip,itp,itt,ii,i,ihr      
      INTEGER IASTART,NB,NE,ITOP      
 


c     first get temperature and humidity on an hourly basis
      RARM = NARM
      ni = 1
      do ip = 1,NPARM 
         t1hr(ip,ni) = t3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM
            dx = t3hr(ip,itp+1)-t3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               t1hr(ip,i) = t1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM  
         q1hr(ip,ni) = q3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM 
            dx = q3hr(ip,itp+1)-q3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               q1hr(ip,i) = q1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         u1hr(ip,ni) = u3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM   
            dx = u3hr(ip,itp+1)-u3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               u1hr(ip,i) = u1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         v1hr(ip,ni) = v3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM
            dx = v3hr(ip,itp+1)-v3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               v1hr(ip,i) = v1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         om1hr(ip,ni) = om3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM 
            dx = om3hr(ip,itp+1)-om3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               om1hr(ip,i) = om1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         wd1hr(ip,ni) = wd3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM 
            dx = wd3hr(ip,itp+1)-wd3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               wd1hr(ip,i) = wd1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         hta1hr(ip,ni) = hta3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM 
            dx = hta3hr(ip,itp+1)-hta3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               hta1hr(ip,i) = hta1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         vsa1hr(ip,ni) = vsa3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM    
            dx = vsa3hr(ip,itp+1)-vsa3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               vsa1hr(ip,i) = vsa1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         hqa1hr(ip,ni) = hqa3hr(ip,1)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM  
            dx = hqa3hr(ip,itp+1)-hqa3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               hqa1hr(ip,i) = hqa1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM 
         vqa1hr(ip,ni) = vqa3hr(ip,ni)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM   
            dx = vqa3hr(ip,itp+1)-vqa3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i = ni+ii
               vqa1hr(ip,i) = vqa1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo 

      ni = 1
      do ip = 1,NPARM  
         acld1hr(ip,ni) = acld3hr(ip,ni)
      enddo
      do itp = 1,NTARM-1
         do ip = 1,NPARM 
            dx = acld3hr(ip,itp+1)-acld3hr(ip,itp)
            dx = dx/RARM
            do ii = 1,NARM
               i=ni+ii
               acld1hr(ip,i) = acld1hr(ip,i-1)+dx
            enddo
         enddo
         ni = ni + NARM
      enddo


c     write(iu_scm_prt,*) 'to time step interpolation done'

c     do l=1,LM
c        write(iu_scm_prt,*) 'l sig sige ',l,sig(l),sige(l)
c     enddo



c     for all the time steps
      do ihr = 1,MCT
         if (AMPS(ihr).gt.0.0) then
C            fill pressure levels
ccccc check how to fill pressure levels
             do L = 1,LS1-1 
c               write(iu_scm_prt,*) ihr,L,SIG(L),PTOP,AMPS(ihr)
                SG_P(L) = SIG(L)*(AMPS(ihr)-PTOP) + PTOP
                SGE_P(L) = SIGE(L)*(AMPS(ihr)-PTOP) + PTOP
             enddo
             do L=LS1,LM
c               write(iu_scm_prt,*) ihr,L,SIG(L),PTOP,AMPS(ihr)
                SG_P(L) = SIG(L)*(PSF-PTOP)+PTOP
                SGE_P(L) = SIGE(L)*(PSF-PTOP)+PTOP 
             enddo
             SGE_P(LM+1) = SIGE(LM+1)*(PSF-PTOP) + PTOP
       
c            do L=1,LM
c               write(iu_scm_prt,*) 'l sge_p sg_p ',l,SGE_P(l),SG_P(l)   
c            enddo
c            write(iu_scm_prt,*) 'l sge_p      ',lm+1,sge_p(lm+1)


             do i=1,NPARM 
                AP(i) = parm(i)
             enddo
C
c            check TWP ARM pressure levels
C            create array of ARM pressure level endpoints
ccccc        APE(1) = AMPS(ihr)
             do n=1,NPARM 
                APE(n) = AP(n)+12.5
             enddo
             APE(NPARM+1) = AP(NPARM)-12.5
             
c            do n=1,NPARM   
c               write(iu_scm_prt,*) 'n APE AP ',n,APE(n),AP(n)        
c            enddo 
c            write(iu_scm_prt,*) 'n APE    ',nparm+1,APE(nparm+1)


             IASTART=1
             if (AMPS(ihr).lt.APE(1)) IASTART=2
             if (AMPS(ihr).lt.APE(2)) IASTART=3
c
             sumaqh = 0.0
             sumaqv = 0.0
             sumath = 0.0
             sumasv = 0.0
             sumat = 0.0
             sumaq = 0.0
             sumau = 0.0
             sumav = 0.0
             sumawd = 0.0
             sumaom = 0.0
             totpa = 0.0
c
             do n=IASTART,NPARM+1   
                if (AMPS(ihr).gt.APE(n)) then
                    if (AMPS(ihr).lt.APE(n-1)) then
                        deltap = AMPS(ihr)-APE(n)
                    else
                        deltap = APE(n-1)-APE(n)
                    endif 
c                   deltap = APE(n-1)-APE(n)
                    sumat = sumat + deltap*t1hr(n-1,ihr)
                    sumaq = sumaq + deltap*q1hr(n-1,ihr)*100./grav
                    sumaqh = sumaqh + deltap*hqa1hr(n-1,ihr)
                    sumaqv = sumaqv + deltap*vqa1hr(n-1,ihr)
                    sumath = sumath + deltap*hta1hr(n-1,ihr)
                    sumasv = sumasv + deltap*vsa1hr(n-1,ihr)
                    sumau = sumau + deltap*u1hr(n-1,ihr)
                    sumav = sumav + deltap*v1hr(n-1,ihr)
                    sumawd = sumawd + deltap*wd1hr(n-1,ihr)
                    sumaom = sumaom + deltap*om1hr(n-1,ihr)
                    totpa = totpa + deltap
                endif  
             enddo
c            write(iu_scm_prt,*) 'sums done'
c            write(iu_scm_prt,414) totpa,sumat
c414          format(1x,'totpa  sumat ',f10.2,f12.2)
c            write(iu_scm_prt,415) totpa,sumaq
c415          format(1x,'totpa  sumaq ',f10.2,f12.5)
c            write(iu_scm_prt,416) totpa,sumaqh*QSC
c416          format(1x,'totpa  sumaqh ',f10.2,f12.5)
c            write(iu_scm_prt,417) totpa,sumaqv*QSC
c417          format(1x,'totpa  sumaqv ',f10.2,f12.5)
c            write(iu_scm_prt,418) totpa,sumath
c418          format(1x,'totpa  sumath ',f10.2,f12.5)
c            write(iu_scm_prt,419) totpa,sumasv
c419          format(1x,'totpa  sumasv ',f10.2,f12.5)
            
C            redistribute ARM data from arm layers to gcm sigma layers
             totpg = 0.0
             sumgqh = 0.0
             sumgqv = 0.0
             sumgth = 0.0
             sumgsv = 0.0
             sumgt = 0.0
             sumgq = 0.0
             sumgu = 0.0
             sumgv = 0.0
             sumgwd = 0.0
             sumgom = 0.0

             do L=1,LM
C               for temperature and humidity advective terms get a weighted
c               average over the layer
C               Find beginning and ending indices for ARM data corresponding to
C               GCM layer
                NB=0
                do n = IASTART,NPARM    
                   if ((APE(n).le.SGE_P(L)).and.
     &                     (APE(n).gt.SGE_P(L+1))) NB=N
                   if (NB.gt.0) go to 500
                enddo 
500             continue
                NE=0
                do N=IASTART,NPARM    
                   if (APE(n).le.SGE_P(L).and.
     &                     APE(n).gt.SGE_P(L+1)) NE=N
c                  if (NE.gt.0) go to 510 
                enddo 
510             continue
 
                DELTAP = SGE_P(L)-SGE_P(L+1)

c               4 cases
c                 1: GCM pressure < ARM TOP - have reached top of
c                    ARM profile 
c                 2: NB=0    arm layer completely contains gcm layer
c                 3: NB=NE   gcm layer overlaps 2 different arm layers
c                 4: NE>NB   gcm layer overlaps more than 2 arm layers
                if (NB.eq.1.and.NE.eq.2) then
                    write(iu_scm_prt,513) IHR,NB,NE
 513                format(1x,'513-SPECIAL CASE IHR NB NE ',i5,i5,i5)
C                   set NB = NE in this case so we go through different case code
                    NB = NE
                endif
                ITOP = 0
                if (SGE_P(L+1).lt.APE(NPARM+1)) then
                    if (SGE_P(L).gt.APE(NPARM+1)) then
c                       part of layer is within top layer of arm
                        if (SGE_P(L).gt.APE(NPARM)) then
c                           model layer overlaps the 2 top ARM layers
                            DELAG1=SGE_P(L)-APE(NPARM)
                            QHR(L,ihr) = DELAG1*q1hr(NPARM-1,ihr)
                            THR(L,ihr) = DELAG1*t1hr(NPARM-1,ihr)
                            UHR(L,ihr) = DELAG1*u1hr(NPARM-1,ihr)
                            VHR(L,ihr) = DELAG1*v1hr(NPARM-1,ihr)
                            OMGHR(L,ihr) = DELAG1*om1hr(NPARM-1,ihr)
                            WDHR(L,ihr) = DELAG1*wd1hr(NPARM-1,ihr)
                            HTA_HR(L,ihr) = DELAG1*hta1hr(NPARM-1,ihr)
                            VSA_HR(L,ihr) = DELAG1*vsa1hr(NPARM-1,ihr)
                            HQA_HR(L,ihr) = DELAG1*hqa1hr(NPARM-1,ihr)
                            VQA_HR(L,ihr) = DELAG1*vqa1hr(NPARM-1,ihr)
                            ACLDHR(L,ihr) = DELAG1*acld1hr(NPARM-1,ihr)
                            DELAG2 = APE(NPARM)-APE(NPARM+1)
                            QHR(L,ihr)=QHR(L,ihr)+DELAG2*q1hr(NPARM,ihr)
                            THR(L,ihr)=THR(L,ihr)+DELAG2*t1hr(NPARM,ihr)
                            UHR(L,ihr)=UHR(L,ihr)+DELAG2*u1hr(NPARM,ihr)
                            VHR(L,ihr)=VHR(L,ihr)+DELAG2*v1hr(NPARM,ihr)
                            OMGHR(L,ihr)=
     +                          OMGHR(L,ihr)+DELAG2*om1hr(NPARM,ihr)
                            WDHR(L,ihr) =
     +                          WDHR(L,ihr)+DELAG2*wd1hr(NPARM,ihr)
                            HTA_HR(L,ihr) = HTA_HR(L,ihr) +
     +                                      DELAG2*hta1hr(NPARM,ihr)
                            VSA_HR(L,ihr) = VSA_HR(L,ihr) +
     +                                      DELAG2*vsa1hr(NPARM,ihr)
                            HQA_HR(L,ihr) = HQA_HR(L,ihr) +
     +                                      DELAG2*hqa1hr(NPARM,ihr)
                            VQA_HR(L,ihr) = VQA_HR(L,ihr) +
     +                                      DELAG2*vqa1hr(NPARM,ihr)
                            ACLDHR(L,ihr) = ACLDHR(L,ihr) + 
     +                                      DELAG2*acld1hr(NPARM,ihr)
                            QHR(L,ihr) = QHR(L,ihr)/DELTAP
                            THR(L,ihr) = THR(L,ihr)/DELTAP
                            UHR(L,ihr) = UHR(L,ihr)/DELTAP
                            VHR(L,ihr) = VHR(L,ihr)/DELTAP
                            OMGHR(L,ihr) = OMGHR(L,ihr)/DELTAP
                            WDHR(L,ihr) = WDHR(L,ihr)/DELTAP
                            HTA_HR(L,ihr) = HTA_HR(L,ihr)/DELTAP
                            VSA_HR(L,ihr) = VSA_HR(L,ihr)/DELTAP
                            HQA_HR(L,ihr) = HQA_HR(L,ihr)/DELTAP
                            VQA_HR(L,ihr) = VQA_HR(L,ihr)/DELTAP
                            ACLDHR(L,ihr) = ACLDHR(L,ihr)/DELTAP
                        else
                           DELAG = SGE_P(L)-APE(NPARM+1)
                           QHR(L,ihr) = q1hr(NPARM,ihr)
                           THR(L,ihr) = t1hr(NPARM,ihr)
                           UHR(L,ihr) = u1hr(NPARM,ihr)
                           VHR(L,ihr) = v1hr(NPARM,ihr)
                           OMGHR(L,ihr) = om1hr(NPARM,ihr)
                           WDHR(L,ihr) = wd1hr(NPARM,ihr)
                           HTA_HR(L,ihr)=DELAG*hta1hr(NPARM,ihr)/DELTAP
                           VSA_HR(L,ihr)=DELAG*vsa1hr(NPARM,ihr)/DELTAP
                           HQA_HR(L,ihr)=DELAG*hqa1hr(NPARM,ihr)/DELTAP
                           VQA_HR(L,ihr)=DELAG*vqa1hr(NPARM,ihr)/DELTAP
                           ACLDHR(L,ihr) = acld1hr(NPARM,ihr)
                        endif
                    else
                        ITOP = 1
                        HTA_HR(L,ihr) = 0.0
                        VSA_HR(L,ihr) = 0.0
                        HQA_HR(L,ihr) = 0.0
                        VQA_HR(L,ihr) = 0.0
c                       what to do at TOP above ARM layers for the
C                       following variables
                        QHR(L,ihr) = QHR(L-1,ihr)
                        THR(L,ihr) = THR(L-1,ihr)
                        UHR(L,ihr) = UHR(L-1,ihr)
                        VHR(L,ihr) = VHR(L-1,ihr)
                        OMGHR(L,ihr) = OMGHR(L-1,ihr)
                        ACLDHR(L,ihr) = 0.0
                        WDHR(L,ihr) = WDHR(L-1,ihr) 
                    endif
                else if (NB.eq.0 .or. (NB.eq.NE.and.NB.eq.1)) then
c                  gcm layer completely contained within the arm layer 
                   do n=1,NPARM    
                     if (SGE_P(L).le.APE(n).and.
     +                  SGE_P(L+1).gt.APE(n+1)) then
                        DELAG = SGE_P(L)-SGE_P(L+1)
                        IF (NB.eq.1) then
                            DELAG = (APE(NB)-APE(NB+1))
     +                           -(SGE_P(L+1)-APE(NB+1))
                        endif
                        QHR(L,ihr) = (q1hr(n,ihr)*DELAG)/DELTAP
                        THR(L,ihr) = (t1hr(n,ihr)*DELAG)/DELTAP
                        UHR(L,ihr) = (u1hr(n,ihr)*DELAG)/DELTAP
                        VHR(L,ihr) = (v1hr(n,ihr)*DELAG)/DELTAP
                        OMGHR(L,ihr) = (om1hr(n,ihr)*DELAG)/DELTAP
                        WDHR(L,ihr) = (wd1hr(n,ihr)*DELAG)/DELTAP
                        HTA_HR(L,ihr) = (hta1hr(n,ihr)*DELAG)/DELTAP
                        VSA_HR(L,ihr) = (vsa1hr(n,ihr)*DELAG)/DELTAP
                        HQA_HR(L,ihr) = (hqa1hr(n,ihr)*DELAG)/DELTAP
                        VQA_HR(L,ihr) = (vqa1hr(n,ihr)*DELAG)/DELTAP
                        ACLDHR(L,ihr) = (acld1hr(n,ihr)*DELAG)/DELTAP
                     endif
                   enddo
                elseif (NB.eq.NE) then
c                   gcm layer overlaps 2 different arm layers 
                    DELAG1 = SGE_P(L)-APE(NB)
                    QHR(L,ihr) = DELAG1*q1hr(NB-1,ihr)
                    THR(L,ihr) = DELAG1*t1hr(NB-1,ihr)
                    UHR(L,ihr) = DELAG1*u1hr(NB-1,ihr)
                    VHR(L,ihr) = DELAG1*v1hr(NB-1,ihr)
                    OMGHR(L,ihr) = DELAG1*om1hr(NB-1,ihr)
                    WDHR(L,ihr) = DELAG1*wd1hr(NB-1,ihr)
                    HTA_HR(L,ihr) = DELAG1*hta1hr(NB-1,ihr)
                    VSA_HR(L,ihr) = DELAG1*vsa1hr(NB-1,ihr)
                    HQA_HR(L,ihr) = DELAG1*hqa1hr(NB-1,ihr)
                    VQA_HR(L,ihr) = DELAG1*vqa1hr(NB-1,ihr)
                    ACLDHR(L,ihr) = DELAG1*acld1hr(NB-1,ihr)
 
                    DELAG2 = APE(NB)-SGE_P(L+1)             
                    QHR(L,ihr) = QHR(L,ihr) + DELAG2*q1hr(NB,ihr)
                    THR(L,ihr) = THR(L,ihr) + DELAG2*t1hr(NB,ihr)
                    UHR(L,ihr) = UHR(L,ihr) + DELAG2*u1hr(NB,ihr)
                    VHR(L,ihr) = VHR(L,ihr) + DELAG2*v1hr(NB,ihr)
                    OMGHR(L,ihr) = OMGHR(L,ihr) + DELAG2*om1hr(NB,ihr)
                    WDHR(L,ihr) = WDHR(L,ihr) + DELAG2*wd1hr(NB,ihr)
                    HTA_HR(L,ihr)=HTA_HR(L,ihr)+DELAG2*hta1hr(NB,ihr)
                    VSA_HR(L,ihr)=VSA_HR(L,ihr)+DELAG2*vsa1hr(NB,ihr)
                    HQA_HR(L,ihr)=HQA_HR(L,ihr)+DELAG2*hqa1hr(NB,ihr)
                    VQA_HR(L,ihr)=VQA_HR(L,ihr)+DELAG2*vqa1hr(NB,ihr)
                    ACLDHR(L,ihr)=ACLDHR(L,ihr)+DELAG2*acld1hr(NB,ihr)
c
                    QHR(L,ihr) = QHR(L,ihr)/DELTAP
                    THR(L,ihr) = THR(L,ihr)/DELTAP
                    UHR(L,ihr) = UHR(L,ihr)/DELTAP
                    VHR(L,ihr) = VHR(L,ihr)/DELTAP
                    OMGHR(L,ihr) = OMGHR(L,ihr)/DELTAP
                    WDHR(L,ihr) = WDHR(L,ihr)/DELTAP
                    HTA_HR(L,ihr) = HTA_HR(L,ihr)/DELTAP
                    VSA_HR(L,ihr) = VSA_HR(L,ihr)/DELTAP
                    HQA_HR(L,ihr) = HQA_HR(L,ihr)/DELTAP
                    VQA_HR(L,ihr) = VQA_HR(L,ihr)/DELTAP
                    ACLDHR(L,ihr) = ACLDHR(L,ihr)/DELTAP

                elseif (NE.gt.NB) then 
                    QHR(L,ihr) = 0.0
                    THR(L,ihr) = 0.0
                    UHR(L,ihr) = 0.0
                    VHR(L,ihr) = 0.0
                    OMGHR(L,ihr) = 0.0
                    WDHR(L,ihr) = 0.0
                    HTA_HR(L,ihr) = 0.0
                    VSA_HR(L,ihr) = 0.0
                    HQA_HR(L,ihr) = 0.0
                    VQA_HR(L,ihr) = 0.0
                    ACLDHR(L,ihr) = 0.0
                    do N=NB,NE
                     if (N.eq.NB) then
                         DELAG = SGE_P(L)-APE(N) 
                     else
                         DELAG = APE(N-1)-APE(N)
                     endif
                     QHR(L,ihr) = QHR(L,ihr) +DELAG*q1hr(N-1,ihr)
                     THR(L,ihr) = THR(L,ihr) +DELAG*t1hr(N-1,ihr)
                     UHR(L,ihr) = UHR(L,ihr) +DELAG*u1hr(N-1,ihr)
                     VHR(L,ihr) = VHR(L,ihr) +DELAG*v1hr(N-1,ihr)
                     OMGHR(L,ihr) = OMGHR(L,ihr) +DELAG*om1hr(N-1,ihr)
                     WDHR(L,ihr) = WDHR(L,ihr) +DELAG*wd1hr(N-1,ihr)
                     HTA_HR(L,ihr)=HTA_HR(L,ihr)+DELAG*hta1hr(N-1,ihr)
                     VSA_HR(L,ihr)=VSA_HR(L,ihr)+DELAG*vsa1hr(N-1,ihr)
                     HQA_HR(L,ihr)=HQA_HR(L,ihr)+DELAG*hqa1hr(N-1,ihr)
                     VQA_HR(L,ihr)=VQA_HR(L,ihr)+DELAG*vqa1hr(N-1,ihr)
                     ACLDHR(L,ihr)=ACLDHR(L,ihr)
     &                                  +DELAG*acld1hr(N-1,ihr)
                    end do
                    DELAG = APE(NE)-SGE_P(L+1)
                    QHR(L,ihr) = QHR(L,ihr) + DELAG*q1hr(NE,ihr)
                    THR(L,ihr) = THR(L,ihr) + DELAG*t1hr(NE,ihr)
                    UHR(L,ihr) = UHR(L,ihr) + DELAG*u1hr(NE,ihr)
                    VHR(L,ihr) = VHR(L,ihr) + DELAG*v1hr(NE,ihr)
                    OMGHR(L,ihr) = OMGHR(L,ihr) +DELAG*om1hr(NE,ihr)
                    WDHR(L,ihr) = WDHR(L,ihr) +DELAG*wd1hr(NE,ihr)
                    HTA_HR(L,ihr) = HTA_HR(L,ihr) +DELAG*hta1hr(NE,ihr)
                    VSA_HR(L,ihr) = VSA_HR(L,ihr) +DELAG*vsa1hr(NE,ihr)
                    HQA_HR(L,ihr) = HQA_HR(L,ihr) +DELAG*hqa1hr(NE,ihr)
                    VQA_HR(L,ihr) = VQA_HR(L,ihr) +DELAG*vqa1hr(NE,ihr) 
                    ACLDHR(L,ihr) = ACLDHR(L,ihr)+DELAG*acld1hr(NE,ihr)
                    QHR(L,ihr) = QHR(L,ihr)/DELTAP
                    THR(L,ihr) = THR(L,ihr)/DELTAP
                    UHR(L,ihr) = UHR(L,ihr)/DELTAP
                    VHR(L,ihr) = VHR(L,ihr)/DELTAP
                    OMGHR(L,ihr) = OMGHR(L,ihr)/DELTAP
                    WDHR(L,ihr) = WDHR(L,ihr)/DELTAP
                    HTA_HR(L,ihr) = HTA_HR(L,ihr)/DELTAP
                    VSA_HR(L,ihr) = VSA_HR(L,ihr)/DELTAP
                    HQA_HR(L,ihr) = HQA_HR(L,ihr)/DELTAP
                    VQA_HR(L,ihr) = VQA_HR(L,ihr)/DELTAP
                    ACLDHR(L,ihr) = ACLDHR(L,ihr)/DELTAP
                else
                    write(iu_scm_prt,561) NB,NE
561                 format(1x,'561-ERROR  NB NE ',i5,i5,'  which case')
                endif                 
                if (ITOP.eq.0) then
                    totpg = totpg + DELTAP
                    sumgt = sumgt + DELTAP*THR(L,ihr)
                    sumgq = sumgq + DELTAP*QHR(L,ihr)*100./GRAV 
                    sumgth = sumgth + DELTAP*HTA_HR(L,ihr)
                    sumgsv = sumgsv + DELTAP*VSA_HR(L,ihr)
                    sumgqh = sumgqh + DELTAP*HQA_HR(L,ihr)
                    sumgqv = sumgqv + DELTAP*VQA_HR(L,ihr) 
                endif
             enddo
             QHR(LM,ihr) = QHR(LM-1,ihr) 
             THR(LM,ihr) = THR(LM-1,ihr) 
             UHR(LM,ihr) = UHR(LM-1,ihr)
             VHR(LM,ihr) = VHR(LM-1,ihr)
             OMGHR(LM,ihr) = OMGHR(LM-1,ihr)
             WDHR(LM,ihr) = WDHR(LM-1,ihr)
             ACLDHR(LM,ihr) = 0.0 
             HTA_HR(LM,ihr) = 0.0
             VSA_HR(LM,ihr) = 0.0
             HQA_HR(LM,ihr) = 0.0
             VQA_HR(LM,ihr) = 0.0
             DELTAP = SGE_P(LM)-SGE_P(LM+1)
             totpg = totpg + DELTAP
             sumgt = sumgt + DELTAP*THR(LM,ihr)
             sumgq = sumgq + DELTAP*QHR(LM,ihr)
c            write(iu_scm_prt,572) LM ,DELTAP,totpg 
 572         format(1x,'LM DELTAP totpg',i4,f10.2,f10.2)
c            write(iu_scm_prt,734) totpg,sumgt
734          format(1x,'totpg sumgt ',f10.2,f12.2)
c            write(iu_scm_prt,735) totpg,sumgq
735          format(1x,'totpg sumgq ',f10.2,f12.5)
c            write(iu_scm_prt,738) totpg,sumgqh*QSC
738          format(1x,'totpg sumgqh ',f10.2,f12.5)
c            write(iu_scm_prt,739) totpg,sumgqv*QSC
739          format(1x,'totpg sumgqv ',f10.2,f12.5)
c            write(iu_scm_prt,736) totpg,sumgth
736          format(1x,'totpg sumgth ',f10.2,f12.5)
c            write(iu_scm_prt,737) totpg,sumgsv
737          format(1x,'totpg sumgsv ',f10.2,f12.5)
c                         
c            write(iu_scm_prt,1205) ihr
c1205        format(/1x,'IHR ',i4)
c            do L = 1,LM
c               write(iu_scm_prt,1210) L,SG_P(L),THR(L,ihr),
c    +            QHR(L,ihr)*1000.0,HTA_HR(L,ihr),VSA_HR(L,ihr) 
c1210           format(1x,'SCM P T Q HTA VSA',
c    &                   i3,f8.2,f8.2,f10.5,F12.7,f12.7)
c            enddo
         endif
      enddo        

      return
      end subroutine arm_to_sig  
