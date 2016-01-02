#include "rundeck_opts.h"
      module bio_inicond_mod
      implicit none

      contains

      subroutine bio_inicond_read(filename, array, depth)
      use pario, only : par_open,par_close
     &     ,read_data,read_dist_data,get_dimlens
#ifdef OBIO_ON_GARYocean
      use oceanr_dim, only : ogrid
#else
      USE hycom_dim, only : ogrid
#endif
      implicit none
      character(len=*), intent(in) :: filename
      real*8, dimension(:,:,:), allocatable, intent(out) :: array
      real*8, dimension(:), allocatable, intent(out), optional :: depth
      integer :: fid, dlens(7), ndims

      fid=par_open(ogrid, filename, 'read')
      call get_dimlens(ogrid, fid, 'array', ndims, dlens)
      if ((dlens(1)/=ogrid%im_world).or.(dlens(2)/=ogrid%jm_world))
     &   call stop_model('dimension mismatch: '//trim(filename), 255)
      allocate(array(ogrid%i_strt:ogrid%i_stop,
     &                  ogrid%j_strt:ogrid%j_stop, dlens(3)))
      if (present(depth)) then
        allocate(depth(dlens(3)))
        call read_data(ogrid, fid, 'depth', depth, bcast_all=.true.)
      endif
      call read_dist_data(ogrid, fid, 'array', array)
      call par_close(ogrid, fid)
      end subroutine bio_inicond_read

      subroutine bio_inicond(filename, fldo)
      use obio_com, only: ze
#ifdef OBIO_ON_GARYocean
      use oceanres, only : kdm=>lmo
      use oceanr_dim, only : ogrid
      use ocean, only : ip=>focean, lmm
#else
      USE hycom_dim, only : kdm, ogrid, ip
#endif
      implicit none
      character(len=*), intent(in) :: filename
      real, dimension(ogrid%i_strt:ogrid%i_stop,
     &    ogrid%j_strt:ogrid%j_stop, kdm), intent(out) ::  fldo
      real*8, dimension(:,:,:), allocatable :: array
      real*8, dimension(:), allocatable :: depth, nodc_d
      real*8, dimension(:), allocatable :: dummy
      logical :: regrid
      integer :: i, j, k, kmax, nodc_kmax
      interface
        Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ, missing)
        Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM), RZ(LM), missing
        end Subroutine VLKtoLZ
      end interface

      call bio_inicond_read(filename, array, depth)
      regrid=kdm/=size(depth)
      if (.not.regrid) regrid=all(abs(depth-
     &                      ze(ogrid%i_strt, ogrid%j_strt, :))<1d0)
      if (regrid) then
#ifdef OBIO_ON_GARYocean
        allocate(dummy(size(fldo, 3)))
        do i=ogrid%i_strt,ogrid%i_stop
          do j=ogrid%j_strt,ogrid%j_stop
            if (ip(i, j)==0) cycle
            call vlktolz(size(depth), lmm(i, j), depth, ze(i, j, :),
     &           array(i, j, :), fldo(i, j, :), dummy, -999999.)
          end do
        end do
#else
        fldo=-9999.d0
        allocate(nodc_d(size(depth)+1))
        do j=ogrid%j_strt,ogrid%j_stop
        do i=ogrid%i_strt,ogrid%i_stop
          if (ip(i,j)==0) cycle
          kmax=1
          do k=1,kdm
            if (ze(i, j, k) .gt. ze(i, j, k-1)) kmax=k+1
          enddo
          do k=1,size(depth)
            if (depth(k) .le. ze(i, j, min(20,kmax))) then
              nodc_d(k)=depth(k)
              nodc_kmax=k
            endif
          enddo
          nodc_d(nodc_kmax+1)=ze(i, j, min(20,kmax))
          call remap1d_plm(array(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
     .             fldo(i,j,1:kdm),ze(i, j, :),kdm,.false.,i,j)
        enddo
        enddo
#endif
      else
        fldo=array
      endif
      end subroutine bio_inicond

      end module bio_inicond_mod


c ----------------------------------------------------------------
      subroutine obio_init
c --- biological/light setup
c ----------------------------------------------------------------
c 
      USE FILEMANAGER, only: openunit,closeunit,file_exists
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT

      USE obio_dim
      USE obio_incom
      use bio_inicond_mod, only: bio_inicond_read
      USE obio_forc, only : ihra,atmFe,alk,surfN
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
      USE OCEANRES, only : kdm=>lmo
      USE MODEL_COM, only: dtsrc
      USE OCEANR_DIM, only : ogrid
#else
      USE hycom_dim_glob, only : kdm
      USE hycom_scalars, only : baclin
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

      character*50 title
!     character*50 cfle
      character cacbc*11,cabw*10
      character*80 filename,fn

      data cacbc,cabw /'acbc25b.dat','abw25b.dat'/

c 
      if (AM_I_ROOT()) print*, 'Ocean Biology setup starts'

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
#ifdef increaseNremin
      remin(1) = 0.5/24.0            !nitrogen
#endif
#ifdef increaseNremin2
      remin(1) = 0.1/24.0            !nitrogen
#endif
#ifdef increaseNremin3
      remin(1) = 0.3/24.0            !nitrogen
#endif
      remin(2) = 0.0001/24.0           !silica
#ifdef increaseSremin
      remin(2) = 0.002/24.0           !silica
#endif
!     remin(3) = 0.020/24.0            !iron
!change June 1, 2010
      remin(3) = 0.50/24.0            !iron
!endofchange
#ifdef increaseIremin
      remin(3) = 0.70/24.0            !iron
#endif

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
        call bio_inicond_read('atmFe_inicond', atmFe)
#ifdef Relax2SurfN
        allocate(surfn(ogrid%i_strt:ogrid%i_stop,
     &                       ogrid%j_strt:ogrid%j_stop))
        call bio_surfN('nitrates_inicond',surfn)
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
        call init_alk(alk)
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

       write(*,'(a,3(f8.6,1x))'), 'OBIO remin rates (per day)=',
     . remin(1)*24.,remin(2)*24.,remin(3)*24.

      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      write(*,*)'**************************************************'
      endif

      return
      end subroutine obio_init

c------------------------------------------------------------------------------

      subroutine bio_surfn(filename, fldo)
      use bio_inicond_mod, only: bio_inicond_read
      use dictionary_mod, only: sync_param
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE OCEAN, only : dlatm
#else
      USE hycom_dim, only : ogrid
      USE GEOM, only : dlatm
#endif
      implicit none
      character(len=*), intent(in) :: filename
      real, intent(out) :: fldo(ogrid%I_STRT:ogrid%I_STOP,
     .          ogrid%J_STRT:ogrid%J_STOP)
      real*8, dimension(:, :, :), allocatable :: array

      call bio_inicond_read(filename, array)
      fldo=array(:, :, 1)
      end subroutine bio_surfn


c------------------------------------------------------------------------------

      subroutine setup_obio
#ifdef OBIO_ON_GARYocean
      use ocn_tracer_com, only: add_ocn_tracer
#endif
      use runtimecontrols_mod, only: tracers_alkalinity
      use exchange_types, only: rad_coupling
      use obio_dim, only: ntrac
      use obio_diag
      implicit none
      integer, dimension(1) :: con_idx
      character(len=10), dimension(1) :: con_str
      integer :: nt, ilim, ll
      character(len=5) :: str1
      character(len=9) :: str2
      character(len=1), parameter :: lim_sym(4)=(/'d', 'h', 'b', 'c'/)  ! diatoms, chloroph, cyanobact, coccoliths
      character(len=4), parameter :: rhs_sym(15)=(/ 'nitr', 'ammo',
     &     'sili', 'iron', 'diat', 'chlo', 'cyan', 'cocc', 'herb',
     &     'ndet', 'sdet', 'idet', 'doc_', 'dic_', 'alk_' /)

      con_idx=[12]
      con_str=['OCN BIOL']
      rad_coupling=.true.
#ifdef OBIO_ON_GARYocean
      call add_ocn_tracer('Nitr      ', i_ntrocn=-4, i_ntrocn_delta=-12,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Ammo      ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Sili      ', i_ntrocn=-4, i_ntrocn_delta=-12,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Iron      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Diat      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Chlo      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Cyan      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Cocc      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Herb      ', i_ntrocn=-8, i_ntrocn_delta=-16,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('Inert     ', i_ntrocn=-4, i_ntrocn_delta=-12,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('N_det     ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('S_det     ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('I_det     ', i_ntrocn=-10,i_ntrocn_delta=-18,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('DOC       ', i_ntrocn=-6, i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      call add_ocn_tracer('DIC       ', i_ntrocn=-3, i_ntrocn_delta=-11,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
      if (tracers_alkalinity)
     &  call add_ocn_tracer('Alk       ',i_ntrocn=-6,i_ntrocn_delta=-14,
     &                 i_con_point_idx=con_idx, i_con_point_str=con_str)
#endif   /* #ifdef OBIO_ON_GARYocean */
#ifdef TOPAZ_params
      call add_diag("co3 ", "oij_co3",
     &              "????????", .false., IJ_co3)
#endif
      call add_diag("ocean surface pH", "oij_pH",
     &              "pH units", .false., IJ_pH)
      call add_diag("Cos Solar Zenith Angle", "oij_solz",
     &              "xxxx", .false., IJ_solz)
      call add_diag("Solar Zenith Angle", "oij_sunz",
     &              "degrees", .false., IJ_sunz)
      call add_diag("Daylight length", "oij_dayl",
     &              "timesteps", .false., IJ_dayl)
      call add_diag("Surface Ocean Direct Sunlight",
     &              "oij_Ed", "quanta", .false., IJ_Ed)
      call add_diag("Surface Ocean Diffuse Sunlight",
     &              "oij_Es", "quanta", .false., IJ_Es)
      call add_diag("Surface ocean Nitrates", "oij_nitr",
     &              "uM", .false., IJ_nitr)
      call add_diag("Surface ocean Ammonium", "oij_amm",
     &              "uM", .false., IJ_amm)
      call add_diag("Surface ocean Silicate", "oij_sil",
     &              "uM", .false., IJ_sil)
      call add_diag("Surface ocean Iron", "oij_iron",
     &              "nM", .false., IJ_iron)
      call add_diag("Surface ocean Diatoms", "oij_diat",
     &              "mg/m3", .false., IJ_diat)
      call add_diag("Surface ocean Chlorophytes", "oij_chlo",
     &              "mg/m3", .false., IJ_chlo)
      call add_diag("Surface ocean Cyanobacteria", "oij_cyan",
     &              "mg/m3", .false., IJ_cyan)
      call add_diag("Surface ocean Coccolithophores", "oij_cocc",
     &              "mg/m3", .false., IJ_cocc)
      call add_diag("Surface ocean Herbivores", "oij_herb",
     &              "mg/m3", .false., IJ_herb)
      call add_diag("Surface ocean DOC", "oij_doc",
     &              "uM", .false., IJ_doc)
      call add_diag("Surface ocean DIC", "oij_dic",
     &              "uM", .false., IJ_dic)
      call add_diag("Surface ocean partial CO2 pressure",
     &              "oij_pCO2", "uatm", .false., IJ_pCO2)
      call add_diag("Surface ocean alkalinity", "oij_alk",
     &              "umol/kg", .false., IJ_alk)
      call add_diag("AO Flux CO2 (ogrid,grC/m2/yr)", "oij_flux",
     &              "grC/m2/yr", .false., IJ_flux)
      call add_diag("C export flux at compensation depth", "oij_cexp",
     &              "mili-grC/m2/hr", .false., IJ_cexp)
      call add_diag("N/C detritus at 74m", "oij_ndet",
     &              "ugC/l", .false., IJ_ndet)
      call add_diag("sink vel n/cdet at 74m", "oij_wsd",
     &              "m/hr", .false., IJ_wsd)
      call add_diag("C export due to chloroph", "oij_xchl",
     &              "kg,C*m/hr", .false., IJ_xchl)
      if (tracers_alkalinity) then
        call add_diag("CaCO3 export flux at compensation depth",
     &                "oij_fca", "mili-g,C/m2/hr", .false., IJ_fca)
      endif

#ifdef OBIO_RUNOFF
#ifdef NITR_RUNOFF
!      call add_diag("Nitrate mass flow from rivers", "oij_rnitrmflo",
!     &               "kg/s", IJ_rnitrmflo)
      call add_diag("Nitrate conc in runoff", "oij_rnitrconc",
     &              "kg/kg", .false., IJ_rnitrconc)
#endif
#ifdef DIC_RUNOFF
      call add_diag("DIC conc in runoff", "oij_rdicconc",
     &              "kg/kg", .false., IJ_rdicconc)
#endif
#ifdef DOC_RUNOFF
      call add_diag("DOC conc in runoff", "oij_rdocconc",
     &              "kg/kg", .false., IJ_rdocconc)
#endif
#ifdef SILI_RUNOFF
      call add_diag("silica conc in runoff", "oij_rsiliconc",
     &              "kg/kg", .false., IJ_rsiliconc)
#endif
#ifdef IRON_RUNOFF
      call add_diag("iron conc in runoff", "oij_rironconc",
     &              "kg/kg", .false., IJ_rironconc)
#endif
#ifdef POC_RUNOFF
      call add_diag("poc conc in runoff", "oij_rpocconc",
     &              "kg/kg", .false., IJ_rpocconc)
#endif
#ifdef ALK_RUNOFF
      call add_diag("alkalinity conc in runoff", "oij_ralkconc",
     &              "mol/kg", .false., IJ_ralkconc)
#endif
#endif
      call add_diag("Depth integrated PP", "oij_pp",
     &              "mg,C/m2/day", .false., IJ_pp)
      call add_diag("PP-diat", "oij_pp1",
     &              "mg,C/m2/day", .false., IJ_pp1)
      call add_diag("PP-chlor", "oij_pp2",
     &              "mg,C/m2/day", .false., IJ_pp2)
      call add_diag("PP-cyan", "oij_pp3",
     &              "mg,C/m2/day", .false., IJ_pp3)
      call add_diag("PP-cocc", "oij_pp4",
     &              "mg,C/m2/day", .false., IJ_pp4)
      do nt=1, 4
        do ilim=1, 5
          write(str1, '(A1,A3,I1)') lim_sym(nt), 'lim', ilim
          call add_diag(str1, str1, "?", .false., ij_lim(nt, ilim))
        end do
      end do
      do nt=1, ntrac-1      ! don't include unused inert tracer
        do ll=1, 17
          write(str2, '(A4,A3,I2.2)') rhs_sym(nt), 'rhs', ll
          call add_diag(str2, str2, "?", .false., ij_rhs(nt, ll))
        end do
      end do

      call add_diag("Mean daily irradiance", "avgq",
     &              "quanta", .true., IJL_avgq)
      call add_diag("KPAR", "kpar",
     &              "??", .true., IJL_kpar)
      call add_diag("dtemp due to kpar", "dtemp_par",
     &              "C", .true., IJL_dtemp)

      return
      end subroutine setup_obio
