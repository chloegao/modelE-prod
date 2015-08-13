#include "rundeck_opts.h"
c ----------------------------------------------------------------
      subroutine obio_init
c --- biological/light setup
c ----------------------------------------------------------------
c 
      USE FILEMANAGER, only: openunit,closeunit,file_exists
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, pack_data

      USE obio_dim
      USE obio_incom
      USE obio_forc, only : ihra,atmFe_glob,atmFe,alk
      USE obio_com, only : npst,npnd,WtoQ,obio_ws,P_tend,D_tend
     .                    ,C_tend,wsdet,gro,obio_deltath,obio_deltat 

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!     .                    ,rnitrmflo_loc
     .                    ,rnitrconc_loc
#endif
#ifdef DIC_RUNOFF
     .                    ,rdicconc_loc
#endif
#ifdef DOC_RUNOFF
     .                    ,rdocconc_loc
#endif
#ifdef SILI_RUNOFF
     .                    ,rsiliconc_loc
#endif
#ifdef IRON_RUNOFF
     .                    ,rironconc_loc
#endif
#ifdef POC_RUNOFF
     .                    ,rpocconc_loc
#endif
#ifdef ALK_RUNOFF
     .                    ,ralkconc_loc
#endif
#endif
#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : idm=>imo,jdm=>jmo,kdm=>lmo
      USE OCEAN, only : ZOE=>ZE,focean,lmm
      USE MODEL_COM, only: dtsrc
      USE OCEANR_DIM, only : ogrid
      use obio_com, only : focean_glob,lmom_glob
#else
      USE hycom_dim_glob, only : idm,jdm,kdm
      USE hycom_scalars, only : nstep,baclin
      USE hycom_dim, only : ogrid
#endif
      USE pario

      use ocalbedo_mod, only: lam, ocalbedo_init=>init
      use RunTimeControls_mod, only: chl_from_seawifs

      implicit none  

      integer i,j,k
      integer iu_bio
      integer nt,nl
      integer imon,ihr,nrec,ichan
      integer lambda,ic
      integer icd,ntr,ich,ih,iu_fac
      integer fid
      real saw,sbw,sac,sbc
      real*4  facirr4(nh,nch,5,ncd)

      real planck,c,hc,oavo,rlamm,rlam450,Sdom,rlam,hcoavo
     .    ,rnn,rbot,pi
     .    ,dummy

      real fldo2(idm,jdm,kdm)
      real fldoz(idm,jdm,kdm)

      character*50 title
!     character*50 cfle
      character cacbc*11,cabw*10
      character*80 filename,fn

      data cacbc,cabw /'acbc25b.dat','abw25b.dat'/

c 
      if (AM_I_ROOT()) print*, 'Ocean Biology setup starts'

#ifdef OBIO_ON_GARYocean
      if(am_i_root()) then
        allocate(focean_glob(idm,jdm),lmom_glob(idm,jdm))
      else
        allocate(focean_glob(1,1),lmom_glob(1,1))
      endif
      call pack_data(ogrid,focean,focean_glob)
      call pack_data(ogrid,lmm,lmom_glob)
#endif

! time steps
#ifdef OBIO_ON_GARYocean
      obio_deltath = dtsrc/3600.d0  !time step in hours
      obio_deltat = obio_deltath    !time step in hrs because all rates are in hrs
#else
      obio_deltath = baclin/3600.d0  !time step in hours
      obio_deltat = obio_deltath    !time step in hrs because all rates are in hrs
#endif

      if (AM_I_ROOT()) 
     . print*, 'Ocean Biology time step(per hour)=',obio_deltath

c  Read in constants, light data
c  Computes constants over entire run of model, reads in required
c  data files, and otherwise obtains one-time-only information
c  necessary for the run.

c  Degrees to radians conversion
      pi = dacos(-1.0D0)
      pi2 = pi*2.0
      rad = 180.0D0/pi

      do nt = 1,nchl
       rkn(nt) = 0.0
       rks(nt) = 0.0
       rkf(nt) = 0.0
      enddo
      pHsfc = 8.0
c
c  Phytoplankton group parameters
      do nt = 1,nchl
       obio_wsd(nt)    = 0.0
       obio_wsh(nt)    = 0.0
      enddo
      do nt = 1,nchl
       rmumax(nt) = 0.0
        rik(1,nt) = 0.0
        rik(2,nt) = 0.0
        rik(3,nt) = 0.0
      enddo
      Pdeep(1) = 32.0    !bottom BC for nitrate
      Pdeep(2) = 0.1     !ammonium
      Pdeep(3) = 60.0    !silica
      Pdeep(4) = 0.6     !iron from Archer and Johnson 2000
      do nt = nnut+1,ntyp
       Pdeep(nt) = 0.0    !chl and herbivores
      enddo
      do nt = 1,ndet
       detdeep(nt) = 0.0 !detritus
      enddo
      cardeep(1) = 0.0   !DOC
      cardeep(2) = 2330.0  !DIC uM(C) from Goyet et al (2000)
c
c  Carbon:chl ratios for different adaptation states
      cchl(1) = 25.0
      cchl(2) = 50.0
      cchl(3) = 80.0
c      cchl(1) = 20.0
c      cchl(2) = 60.0
c      cchl(3) = 100.0
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

!change: March 15, 2010
       bn = cchl(2)/cnratio         !N:chl ratio (uM/ugl)
       bf = cchl(2)/cfratio         !Fe:chl ratio (nM/ugl)
       cchlratio = cchl(2)          !C:chl ratio (ugl/ugl)
       mgchltouMC = cchlratio/uMtomgm3
c
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then
c  Diatoms
      nt = 1
!change: March 15, 2010
!     rmumax(nt) = 1.50       !u max in /day at 20C
      rmumax(nt) = 2.00       !u max in /day at 20C
#ifdef OBIO_ON_GARYocean
#ifdef unlimitDIATOMS
      obio_wsd(nt)    = 0.50  !sinking rate in m/day
#else
      obio_wsd(nt)    = 0.75  !sinking rate in m/day
#endif
#else
      obio_wsd(nt)    = 0.50  !sinking rate in m/day   !!change Oct27,2008
#endif
      rik(1,nt)  = 90.0       !low light-adapted Ik (<50 uE/m2/s)
      rik(2,nt)  = 93.0       !medium light-adapted Ik (50-200 uE/m2/s)
      rik(3,nt)  = 184.0      !high light adapted Ik (>200 uE/m2/s)
      rkn(nt) = 1.0           !M-M half-sat constant for nitrogen
      rks(nt) = 0.2           !M-M half-sat constant for silica
      rkf(nt) = 0.12          !M-M half-sat constant for iron

      endif
!!#endif
!!#if NCHL_DEFINED > 1
      if (nchl > 1) then
c  Chlorophytes
      nt = 2
      rmumax(nt) = rmumax(nt-1)*0.840
      obio_wsd(nt)    = 0.25
      rik(1,nt)  = rik(1,nt-1)*1.077
      rik(2,nt)  = rik(2,nt-1)*0.935
      rik(3,nt)  = rik(3,nt-1)*0.781
      rkn(nt) = rkn(nt-1)*0.75
      rkn(nt) = rkn(nt-1)*0.667   !1/3 distance bet. cocco and dia
      rkf(nt) = rkf(nt-1)*0.835   !midway between cocco's and diatoms
      rkf(nt) = rkf(nt-1)*0.779   !1/3 distance bet. cocco and dia

      endif
!!#endif
!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
c  Cyanobacteria
      nt = 3
      rmumax(nt) = rmumax(nt-2)*0.670
      obio_wsd(nt)    = 0.0085
      rik(1,nt)  = rik(1,nt-2)*0.723
      rik(2,nt)  = rik(2,nt-2)*0.710
      rik(3,nt)  = rik(3,nt-2)*0.256
      rkn(nt) = rkn(nt-2)*0.50
      rkf(nt) = rkf(nt-2)*0.67  !equals cocco

      endif
!!#endif
!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
c  Coccolithophores
      nt = 4
#ifdef RMUMAX_allcocco
       !change 11/5/09
       rmumax(nt) = rmumax(nt-3)*0.663   !all coccos
#else
       !default
       rmumax(nt) = rmumax(nt-3)*0.755   !E. huxleyi only
#endif
c      rmumax(nt) = rmumax(nt-3)*0.781   !E. huxleyi only (no Sunda/Hunts)
      obio_wsd(nt)    = 0.82
      obio_wsd(nt)    = 0.648
      rik(1,nt)  = rik(1,nt-3)*0.623
      rik(2,nt)  = rik(2,nt-3)*0.766
      rik(3,nt)  = rik(3,nt-3)*0.899
      rkn(nt) = rkn(nt-3)*0.5
      rkf(nt) = rkf(nt-3)*0.67

      endif
!!#endif
!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
c  Dinoflagellates
      nt = 5
      rmumax(nt) = rmumax(nt-4)*0.335
      obio_wsd(nt)    = 0.0
      rik(1,nt)  = rik(1,nt-4)*1.321
      rik(2,nt)  = rik(2,nt-4)*1.381
      rik(3,nt)  = rik(3,nt-4)*1.463
      rkn(nt) = rkn(nt-4)*1.0
      rkf(nt) = rkf(nt-4)*0.67

      endif
!!#endif
      do nt = 1,nchl
       obio_wsh(nt) = obio_wsd(nt)/24.0  !convert to m/hr
      enddo

c  Detrital sinking rates m/h
#ifdef limitEXPORT
      wsdeth(1) = 50.0/24.0     !nitrogen
#else
      !default
!change: March 10, 2010
!     wsdeth(1) = 30.0/24.0     !nitrogen
      wsdeth(1) = 20.0/24.0     !nitrogen
#endif
      wsdeth(2) = 50.0/24.0     !silica
!     wsdeth(3) = 20.0/24.0     !iron
!change June 1, 2010
      wsdeth(3) =  5.0/24.0     !iron
!endofchange
c
c  Detrital remineralization rates /hr
!change: March 10, 2010
!     remin(1) = 0.010/24.0            !nitrogen
      remin(1) = 0.020/24.0            !nitrogen
      remin(2) = 0.0001/24.0           !silica
!     remin(3) = 0.020/24.0            !iron
!change June 1, 2010
      remin(3) = 0.50/24.0            !iron
!endofchange
      fescavrate(1) = 2.74E-5/24.0      !low fe scavenging rate/hr
      fescavrate(2) = 50.0*fescavrate(1) !high fe scavenging rate/hr
c
c  (originally done inside lidata subroutine of obio_daysetrad)
c  Reads in radiative transfer data: specifically
c  water data (seawater absorption and total scattering coefficients,
c  and chl-specific absorption and total scattering data for
c  several phytoplankton groups).  PAR (350-700) begins at index 3,
c  and ends at index 17.
c     
c  Water data files
!     cfle = cabw                       
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')


c  Phytoplankton group chl-specific absorption and total scattering
c  data.  Chl-specific absorption data is normalized to 440 nm; convert
c  here to actual ac*(440)
!     cfle = cacbc
!     open(4,file='/explore/nobackup/aromanou/2.0deg/'//cfle
!    .      ,status='old',form='formatted')
      call ocalbedo_init
      call openunit('cfle2',iu_bio)
      do ic = 1,6
       read(iu_bio,'(a50)')title
      enddo
      do nt = 1,nchl
       read(iu_bio,'(a50)')title
       do nl = 1,19
        read(iu_bio,30)lambda,sac,sbc
        ac(nt,nl) = sac
        bc(nt,nl) = sbc
       enddo
       do nl = 20,nlt
        ac(nt,nl) = 0.0
        bc(nt,nl) = 0.0
       enddo
      enddo
      call closeunit(iu_bio)
 30   format(i4,2f10.4)

!ifst part from daysetrad.f
c      h = 6.6256E-34   !Plancks constant J sec
       planck = 6.6256E-34   !Plancks constant J sec
       c = 2.998E8      !speed of light m/sec
c      hc = 1.0/(h*c)
       hc = 1.0/(planck*c)
       oavo = 1.0/6.023E23   ! 1/Avogadros number
       hcoavo = hc*oavo
       do nl = npst,npnd
        rlamm = float(lam(nl))*1.0E-9  !lambda in m
        WtoQ(nl) = rlamm*hcoavo        !Watts to quanta conversion
       enddo
       !CDOM absorption exponent
       rlam450 = 450.0
       Sdom = 0.014
       do nl = 1,nlt
        if (lam(nl) .eq. 450)nl450 = nl
        rlam = float(lam(nl))
        excdom(nl) = exp(-Sdom*(rlam-rlam450))
       enddo
       if (nl450.eq.0) stop 'obio_init: nl450=0'
       !First time thru set ihra to 1 to assure correct value read in
       !from restart file
       !!do j=1,jj
       !!do l=1,isp(j)
       !!do i=ifp(j,l),ilp(j,l)
        !!ihra(i,j) = 1
       !!enddo
       !!enddo
       !!enddo

!ifst part from edeu.f
       bbw = 0.5            !backscattering to forward scattering ratio
       rmus = 1.0/0.83      !avg cosine diffuse down
       Dmax = 500.0         !depth at which Ed = 0

       rnn = 1.341
       rmuu = 1.0/0.4             !avg cosine diffuse up
       rbot = 0.0                 !bottom reflectance
       rd = 1.5   !these are taken from Ackleson, et al. 1994 (JGR)
       ru = 3.0

c  Read in factors to compute average irradiance
! (this part originally done inside obio_edeu)
      if (AM_I_ROOT()) then
      print*, '    '
      print*,'Reading factors for mean irradiance at depth...'
      print*,'nh,nch,ncd=',nh,nch,ncd
      endif

      call openunit('facirr',iu_fac)
      do icd=1,ncd
       do ntr=1,5
        do ich=1,nch
         do  ih=1,nh
          read(iu_fac,*)facirr4(ih,ich,ntr,icd)
           facirr(ih,ich,ntr,icd)=1.D0*facirr4(ih,ich,ntr,icd)
         enddo
        enddo
       enddo
      enddo
      call closeunit(iu_fac)
!     if (AM_I_ROOT()) then
!     print*,'nstep, facirr(10,18,1)=',nstep,facirr4(10,18,1,1)
!     print*,'nstep, facirr(10,18,1)=',nstep,facirr(10,18,1,1)
!     print*, '    '
!     endif

!ifst part from ptend.f
       do k=1,kdm

         do nt=1,nchl
          obio_ws(k,nt) = 0.0
          gro(k,nt) = 0.0
         enddo
 
         do nt=1,ntyp+n_inert
          P_tend(k,nt) = 0.0
         enddo
 
         do nt=1,ndet
          D_tend(k,nt) = 0.0
          wsdet(k,nt) = 0.0
         enddo
 
         do nt=1,ncar
          C_tend(k,nt) = 0.0
         enddo
       enddo
       do nt=1,nchl
        obio_ws(kdm+1,nt)=0.0
       enddo
       do nt=1,ndet
        wsdet(kdm+1,nt) = 0.0
       enddo
 
!read in atmospheric iron deposition (this will also be changed later...)
      if (AM_I_ROOT()) then
      print*, '    '
      print*, 'reading iron data.....'
      print*, '    '
      endif

      if(file_exists('ironflux')) then ! read netcdf format flux
	fid = par_open(ogrid,'ironflux','read')
	call read_dist_data(ogrid,fid,'ironflux',atmFe)
	call par_close(ogrid,fid)
      else
!     open(unit=iu_bio,file='atmFedirect0'
!    . ,form='unformatted',status='old',access='direct' 
!    . ,recl=idm*jdm*8/4)
!     do imon=1,12  !1 year of monthly values
!      nrec=imon
!      read (iu_bio,rec=nrec)((atmFe_glob(i,j,imon),i=1,idm),j=1,jdm)
!     enddo
!     close(iu_bio)
        filename='atmFe_inicond'
#ifdef OBIO_ON_GARYocean
        call bio_inicond2D_g(filename,atmFe(:,:,:),.true.)
#else
        call bio_inicond2D(filename,atmFe(:,:,:),.true.)
#endif
      endif ! netcdf iron or not

#ifdef OBIO_RUNOFF
! read in nutrient concentrations, already regridded to model grid
	if (AM_I_ROOT()) then
	print*, '    '
	print*, 'reading nutrient runoff data.....'
	print*, '    '
	endif
#ifdef NITR_RUNOFF
!        filename='rnitr_mflo'
        filename='rnitr_conc'
	fid=par_open(ogrid,filename,'read')
!	call read_dist_data(ogrid,fid,'din',rnitrmflo_loc)
	call read_dist_data(ogrid,fid,'din',rnitrconc_loc)
	call par_close(ogrid,fid)
#endif
#ifdef DIC_RUNOFF
	filename='rdic_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'dic',rdicconc_loc)
	call par_close(ogrid,fid)
	write(*,*)'reading dic from',filename
#endif
#ifdef DOC_RUNOFF
	filename='rdoc_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'doc',rdocconc_loc)
	call par_close(ogrid,fid)
#endif
#ifdef SILI_RUNOFF
	filename='rsili_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'sil',rsiliconc_loc)
	call par_close(ogrid,fid)
#endif
#ifdef IRON_RUNOFF
	filename='riron_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'fe',rironconc_loc)
	call par_close(ogrid,fid)
#endif
#ifdef POC_RUNOFF
	filename='rpoc_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'poc',rpocconc_loc)
	call par_close(ogrid,fid)
#endif
#ifdef ALK_RUNOFF
	filename='ralk_conc'
	fid=par_open(ogrid,filename,'read')
	call read_dist_data(ogrid,fid,'alk',ralkconc_loc)
	call par_close(ogrid,fid)
#endif
#endif

#ifdef TRACERS_Alkalinity
! Alkalinity will be read in from obio_bioinit
! don't do anything here
#else
!read in alkalinity annual mean file
      if (ALK_CLIM.eq.1) then      !read from climatology
        filename='alk_inicond'
#ifdef OBIO_ON_GARYocean
        call bio_inicond_g(filename,fldo2,fldoz)
        alk(:,:,:)=fldo2

        !remove negative values
        !negs are over land or under ice due to GLODAP missing values in the Arctic Ocean
        !for under ice missing values, use climatological minimums for sets of layers based
        !on GLODAP, rather than setting to the same global min.
        do j=1,jdm
        do i=1,idm
        do k=1,kdm
         if (alk(i,j,k).lt.0. .and. focean_glob(i,j).gt.0) then
            if (zoe(k).le.150.) alk(i,j,k)=2172.      !init neg might be under ice,
            if (zoe(k).gt.150. .and. zoe(k).lt.1200.) alk(i,j,k)=2200.
            if (zoe(k).ge.1200.) alk(i,j,k)=2300.
         endif
         if (focean_glob(i,j).le.0) then
           alk(i,j,k)=0.
         endif
        enddo
        enddo
        enddo
#else
        call bio_inicond(filename,alk(:,:,:))
#endif
      else      !set to zero, obio_carbon sets alk=tabar*sal/sal_mean
        alk = 0.
      endif
#endif

! printout some key information
      if (AM_I_ROOT()) then
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'           INITIALIZATION                         '

      write(*,'(a,i5)') 'OBIO - NUMBER OF TRACERS=',ntrac

      write(*,*)'ALK_CLIM = ', ALK_CLIM
      if (ALK_CLIM.eq.0) write(*,*) 'ALKALINITY, from SALINITY'
      if (ALK_CLIM.eq.1) write(*,*) 'ALKLNTY, GLODAP annmean'
      if (ALK_CLIM.eq.2) write(*,*) 'ALKALINITY prognostic'

#ifdef pCO2_ONLINE
      print*, 'PCO2 is computed online and not through lookup table'
#else
      print*, 'PCO2 is computed through lookup table'
#endif
      write(*,'(a,4e12.4)')'obio_init, sinking rates for chlorophyll: ',
     .    obio_wsd(1),obio_wsd(2),obio_wsd(3),obio_wsd(4)
      write(*,'(a,3e12.4)')'obio_init, settling rates for detritus: ',
     .      wsdeth(1),  wsdeth(2),  wsdeth(3)

#ifdef limitDIC1
      print*,'limit DIC to +0.5%'
#endif
#ifdef limitDIC2
      print*,'limit DIC to +0.2%'
#endif

      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      endif

      return
      end
c------------------------------------------------------------------------------
#ifdef OBIO_ON_GARYocean

      subroutine bio_inicond2D_g(filename,fldo,dateline)

!read in a field at 1x1 resolution
!convert to atmospheric grid
!convert to ocean grid 
!this routine only for (i,j,monthly) arrays

      USE FILEMANAGER, only: openunit,closeunit
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only : idm=>imo,jdm=>jmo
      USE OCEAN, only : oDLATM=>DLATM

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,unpack_data

      implicit none

      integer i,j,k,l,n
      integer, parameter :: igrd=360,jgrd=180,kgrd=12
      integer iu_file,lgth

      real data(igrd,jgrd,kgrd)
      real data_min(kgrd),data_max(kgrd)
      real data_mask(igrd,jgrd)
      real fldo(ogrid%I_STRT_HALO:ogrid%I_STOP_HALO,
     .          ogrid%J_STRT_HALO:ogrid%J_STOP_HALO,kgrd)
      real fldo_glob(idm,jdm,kgrd)
      real idl_n

      logical vrbos,dateline

      character*80 filename
      real*8 dlata,offib,datmis

      if ( AM_I_ROOT() ) then

!--------------------------------------------------------------
      lgth=len_trim(filename)
      print*, 'obio_init: reading from file...',filename(1:lgth)
      call openunit(filename,iu_file,.false.,.true.)

!iron gocart data start from dateline
!      missing values are -9999
       do k=1,kgrd
       data_min(k)=1.e10
       data_max(k)=-1.e10
         do i=1,igrd
           do j=1,jgrd
             read(iu_file,'(e12.4)')data(i,j,k)
             !preserve the mean for later
              if (data(i,j,k)>0.) then
                data_min(k)=min(data_min(k),data(i,j,k))
                data_max(k)=max(data_max(k),data(i,j,k))
             endif
           enddo
         enddo
      enddo
      call closeunit(iu_file)

!--------------------------------------------------------------
      do k=1,kgrd

      do i=1,igrd
      do j=1,jgrd

        !mask
        data_mask(i,j)=0.d0
        if (data(i,j,k)>=0.d0) data_mask(i,j)=1.d0

cdiag   if (k.eq.1)
cdiag.    write(*,'(a,3i5,2e12.4)')'before hntr80 ',
cdiag.    i,j,1,data(i,j,k),data_mask(i,j)


      enddo    ! i-loop
      enddo    ! j-loop

      !iron gocart data start from dateline
      if (dateline) idl_n=0.d0   !no of poits away from dateline
      dlata = 60d0
      offib = 0d0
      datmis = -9999d0
      call HNTR80(igrd,jgrd,idl_n,dlata,
     .             idm,jdm,offib,oDLATM,datmis)

      !use hntr8p in order to get correct polar value: 
      !i.e. average longitudinal value everywhere

      call HNTR8P (data_mask,data(:,:,k),fldo_glob(:,:,k))

cdiag if (k.eq.1) then
cdiag do i=1,idm
cdiag do j=1,jdm
cdiag  write(*,'(a,3i5,e12.4)')'after hntr8p ',
cdiag.    i,j,k,fldo_glob(i,j,k)
cdiag enddo
cdiag enddo
cdiag endif

      enddo    ! k-loop

!--------------------------------------------------------------

      endif   !AM_I_ROOT

      call unpack_data(ogrid, fldo_glob, fldo)

      return
  
      end subroutine bio_inicond2D_g

#endif

c------------------------------------------------------------------------------
#ifndef OBIO_ON_GARYocean
      subroutine bio_inicond2D(filename,fldo,dateline)

!read in a field at 1x1 resolution
!convert to atmospheric grid
!convert to ocean grid (using Shana's routine) 
!this routine only for (i,j,monthly) arrays

c --- mapping flux-like field from agcm to ogcm
c     input: flda (W/m*m), output: fldo (W/m*m)
c

      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only : DLATM      !here okay to use dlatm because interpolate from atmos

      USE hycom_dim_glob, only : jj,isp,ifp,ilp,iia,jja,iio,jjo
      USE hycom_dim, only : ogrid,i_0h,i_1h,j_0h,j_1h
      USE hycom_cpler, only: wlista2o,ilista2o,jlista2o,nlista2o
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,unpack_data
      implicit none

      integer i,j,k,l,n
      integer, parameter :: igrd=360,jgrd=180,kgrd=12
      integer iu_file,lgth
      integer i1,j1,iii,jjj,isum,kmax

      real data_mask(igrd,jgrd)
      real data(igrd,jgrd,kgrd)
      real data2(iia,jja,kgrd)
      real data_min(kgrd),data_max(kgrd)
      real sum1
      real dummy1(36,jja,kgrd),dummy2(36,jja,kgrd)
      real fldo(i_0h:i_1h,j_0h:j_1h,kgrd)
      real fldo_glob(iio,jjo,kgrd)
      real idl_n

      logical vrbos,dateline

      character*80 filename

      if ( AM_I_ROOT() ) then

!--------------------------------------------------------------
      lgth=len_trim(filename)
      if (AM_I_ROOT())
     .print*, 'obio_init: reading from file...',filename(1:lgth)
      call openunit(filename,iu_file,.false.,.true.)

!NOTE: data starts from Greenwich
!      missing values are -9999
!iron gocart data starts from dateline
       do k=1,kgrd
       data_min(k)=1.e10
       data_max(k)=-1.e10
         do i=1,igrd
           do j=1,jgrd
             read(iu_file,'(e12.4)')data(i,j,k)
             !preserve the mean for later
              if (data(i,j,k)>0.) then
                data_min(k)=min(data_min(k),data(i,j,k))
                data_max(k)=max(data_max(k),data(i,j,k))
             endif
           enddo
         enddo
      enddo
      call closeunit(iu_file)

!--------------------------------------------------------------
! convert to the atmospheric grid
! this code is also present in obio_bioinit_g
      do k=1,kgrd
      do i=1,igrd
      do j=1,jgrd
        !mask
        data_mask(i,j)=0.d0
        if (data(i,j,k)>=0.d0) data_mask(i,j)=1.d0
cdiag   if (k.eq.1)
cdiag.    write(*,'(a,3i5,2e12.4)')'before hntr80 ',
cdiag.    i,j,1,data(i,j,k),data_mask(i,j)
      enddo    ! i-loop
      enddo    ! j-loop

      !compute glb average and replace missing data
      if (dateline) idl_n=0.d0   !no of poits away from dateline
      call HNTR80(igrd,jgrd,idl_n,60.d0,
     .             iia,jja,0.d0,DLATM,-9999.d0)
      call HNTR8P (data_mask,data(:,:,k),data2(:,:,k))
      enddo    ! k-loop

!     !this is needed for dic
!     do k=1,kgrd
!     do i=1,iia
!     do j=1,jja
!     if (data2(i,j,k).gt.0) then
!       if (data2(i,j,k)<data_min(k)) data2(i,j,k)=data_min(k)
!       if (data2(i,j,k)>data_max(k)) data2(i,j,k)=data_max(k)
!     endif
!     enddo
!     enddo
!     enddo

      !--------------------------------------------------------
      !***************** important! start from dateline
!     if (.not.dateline) then
!     !move to dateline
!     dummy1=data2(1:36,:,:);
!     dummy2=data2(37:72,:,:);
!     data2(1:36,:,:)=dummy2;
!     data2(37:72,:,:)=dummy1;
!     endif

      !--------------------------------------------------------
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)

      do 9 k=1,kgrd
      fldo_glob(i,j,k)=0.
c
      do 91 n=1,nlista2o(i,j)
      fldo_glob(i,j,k)=fldo_glob(i,j,k)
     .           +data2(ilista2o(i,j,n),jlista2o(i,j,n),k)
     .                       *wlista2o(i,j,n)
 91   continue
 9    continue
 8    continue

!     !this is needed for dic
!     do k=1,kgrd
!     do i=1,iio
!     do j=1,jjo
!       if (fldo(i,j,k)<data_min(k)) fldo(i,j,k)=data_min(k)
!       if (fldo(i,j,k)>data_max(k)) fldo(i,j,k)=data_max(k)
!     enddo
!     enddo
!     enddo

      !--------------------------------------------------------
      endif   !if am-i-root

      call unpack_data(ogrid, fldo_glob, fldo)
      return
  
      end subroutine bio_inicond2D
#endif

      subroutine setup_obio
#ifdef OBIO_ON_GARYocean
      use ocn_tracer_com, only: add_ocn_tracer
      use runtimecontrols_mod, only: tracers_alkalinity
#endif
      use exchange_types, only: rad_coupling
      implicit none

      rad_coupling=.true.
#ifdef OBIO_ON_GARYocean
      call add_ocn_tracer('Nitr      ', i_ntrocn=-4)
      call add_ocn_tracer('Ammo      ', i_ntrocn=-6)
      call add_ocn_tracer('Sili      ', i_ntrocn=-4)
      call add_ocn_tracer('Iron      ', i_ntrocn=-8)
      call add_ocn_tracer('Diat      ', i_ntrocn=-8)
      call add_ocn_tracer('Chlo      ', i_ntrocn=-8)
      call add_ocn_tracer('Cyan      ', i_ntrocn=-8)
      call add_ocn_tracer('Cocc      ', i_ntrocn=-8)
      call add_ocn_tracer('Herb      ', i_ntrocn=-8)
      call add_ocn_tracer('Inert     ', i_ntrocn=-4)
      call add_ocn_tracer('N_det     ', i_ntrocn=-6)
      call add_ocn_tracer('S_det     ', i_ntrocn=-6)
      call add_ocn_tracer('I_det     ', i_ntrocn=-10)
      call add_ocn_tracer('DOC       ', i_ntrocn=-6)
      call add_ocn_tracer('DIC       ', i_ntrocn=-3)
      if (tracers_alkalinity)
     &             call add_ocn_tracer('Alk       ', i_ntrocn=-6)
#endif   /* #ifdef OBIO_ON_GARYocean */

      return
      end subroutine setup_obio
