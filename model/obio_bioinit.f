#include "rundeck_opts.h"

      subroutine obio_bioinit
 
!note: 
!obio_bioinit is called only for a cold start and reads in INITIAL conditions and interpolates them
!to the ocean grid. Such fields are nitrates,silicate,dic
!obio_init  is called for every start of the run and reads in BOUNDARY conditions and interpolates 
!them to ocean grid. such fields are iron,alkalinity,chlorophyl

c based on /g6/aromanou/Watson_new/BioInit/rstbio.F

c  Makes initialization data files for biological variables.
c  This is for the global model.  To subset, use the routines in
c  /u2/gregg/bio/biodat/subreg.
c  Uses NOAA 2001 atlas for NO3 and SiO2 distributions.
c  Includes initial iron distributions.
c
c  Particle type 1  = nitrate
c  Particle type 2  = ammonium
c  Particle type 3  = silicate
c  Particle type 4  = iron
c  Particle type 5  = diatoms
c  Particle type 6  = chlorophytes
c  Particle type 7  = cyanobacteria
c  Particle type 8  = coccolithophores
c  Particle type 9  = dinoflagellates
c  Particle type 10 = zooplankton
c
c  Detritus type 1  = carbon/nitrogen
c  Detritus type 2  = silica
c  Detritus type 3  = iron
c
c  Carbon type 1    = semi-labile DOC
c  Carbon type 2    = DIC
 
      USE FILEMANAGER, only: openunit,closeunit

      USE obio_dim
      USE obio_incom
      USE obio_forc, only: avgq
      USE obio_com, only: gcmax, tracer
 
#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo,dzo
      USE OCEAN, only : ZOE=>ZE, ip=>focean
      USE OCEANR_DIM, only : ogrid
#else
      USE hycom_dim, only : ip,kdm,ogrid
      USE hycom_arrays, only : dpinit
      USE hycom_scalars, only: onem
#endif

      implicit none

      real, parameter, dimension(0:13) :: fer_values=
     &  (/ 0E0,
     &     3.5E-3,      !Antarctic
     &     20.0E-3,     !South Indian
     &     4.5E-3,      !South Pacific
     &     20.0E-3,     !South Atlantic
     &     22.5E-3,     !Equatorial Indian
     &     4.0E-3,      !Equatorial Pacific
     &     20.0E-3,     !Equatorial Atlantic
     &     25.0E-3,     !North Indian
     &     6.0E-3,      !North Central Pacific
     &     20.0E-3,     !North Central Atlantic
     &     6.0E-3,      !North Pacific
     &     10.0E-3,     !North Atlantic
     &     20.0E-3 /)   !Mediterranean

      integer i,j,k,l


      integer nir(nrg), nt

      INTEGER :: j_0,j_1,i_0,i_1

      real zz

      integer, ALLOCATABLE, DIMENSION(:,:)   :: ir
      real,  ALLOCATABLE, DIMENSION(:,:,:) :: fer,dic

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP

      ALLOCATE(ir(i_0:i_1,j_0:j_1))
      ALLOCATE(fer(i_0:i_1,j_0:j_1,kdm))
      allocate(dic(i_0:i_1,j_0:j_1,kdm))

      tracer(:,:,:,1:ntyp)=0.d0
      Fer(:,:,:) = 0.d0

      call bio_inicond('nitrates_inicond',tracer(:,:,:,1))
      call bio_inicond('silicate_inicond',tracer(:,:,:,3))


#ifdef obio_TRANSIENTRUNS
! in the transient runs keep the dic read in from RSF/AIC file
      dic = tracer(:,:,:,15)
#else
! otherwise take from rundeck
      call bio_inicond('dic_inicond',dic)
#endif

#ifdef TRACERS_Alkalinity
      call init_alk(tracer(:,:,:,ntrac))
#endif

!     /archive/u/aromanou/Watson_new/BioInit/iron_ron_4x5.asc
!     /archive/u/aromanou/Watson_new/BioInit/CHL_WG_4x5


!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the hycom grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model

      do k=1,kdm
       do j=j_0,j_1
        do i=i_0,i_1
          if(tracer(i,j,k,1).le.0.)tracer(i,j,k,1)=0.085d0
          if(tracer(i,j,k,3).le.0.)tracer(i,j,k,3)=0.297d0
          if (dic(i,j,k).le.0.) dic(i,j,k)=1837.d0
          dic(i,j,k)=max(dic(i,j,k),1837.0)   !set minimum =1837
        enddo
       enddo
      enddo

      write(*,'(a,2e12.4)')'BIO: bioinit: dic min-max=',
     .       minval(dic),maxval(dic)

c  Obtain region indicators
      write(6,*)'calling fndreg...'
      call fndreg(ir)
 
c  Define Fe:NO3 ratios by region, according to Fung et al. (2000)
c  GBC.  Conversion produces nM Fe, since NO3 is as uM

      do j=j_0,j_1
        do i=i_0,i_1
          fer(i,j,:)=fer_values(ir(i,j))
        enddo  !i-loop
      enddo  !j-loop
 
c  Create arrays 
      write(6,*)'Creating bio restart data for ',ntyp,' arrays and'
     . ,kdm,'  layers...'

      do j=j_0,j_1
      do i=i_0,i_1
        if (ip(i,j)==0) cycle

        zz=0.d0
        do k=1,kdm
          !Nitrate
          !read earlier from file

#ifndef OBIO_ON_GARYocean
          zz=zz+dpinit(i,j,k)/onem   !depth in m
#endif
          !Ammonium
          tracer(i,j,k,2) = 0.5
          !!!if (zz.gt.4000.d0)tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file
          !!!if (zz.gt. 4000.d0)tracer(i,j,k,2) = Pdeep(3)

          !Iron
          tracer(i,j,k,4) = Fer(i,j,k)*tracer(i,j,k,1)  !Fung et al. 2000
c          if (ir(nw) .eq. 3)then
c           P(i,j,k,4) = 0.04*float(k-1) + 0.2
c           P(i,j,k,4) = 0.06*float(k-1) + 0.2
c           P(i,j,k,4) = 0.08*float(k-1) + 0.2
c           P(i,j,k,4) = min(P(i,j,k,4),0.65)
c           P(i,j,k,4) = min(P(i,j,k,4),0.75)
c          endif
          if (ir(i,j) .eq. 1)then
           tracer(i,j,k,4) = Fer(i,j,k)*0.5*tracer(i,j,k,1)
          endif
          tracer(i,j,k,4) = max(tracer(i,j,k,4),0.01)
          !!!if (zz.gt. 4000.)tracer(i,j,k,4) = Pdeep(4)

          !Herbivores
          do nt = nnut+1,ntyp-nzoo
           tracer(i,j,k,nt) = 0.05
          enddo
          do nt = ntyp-nzoo+1,ntyp
           tracer(i,j,k,nt) = 0.05  !in chl units mg/m3
c          tracer(i,j,k,nt) = 0.05*50.0  !in C units mg/m3
          enddo

          !inert tracer
          do nt = ntyp+1,ntyp+n_inert
           tracer(i,j,k,nt) = tracer(i,j,k,1)
          enddo

          !DIC
          !read earlier from file
#ifdef limitDIC1
!!!       dic(i,j,k)=dmax1(1837d0,0.99*dic(i,j,k))  !!! g6hh
          dic(i,j,k)=dmax1(1837d0,1.005*dic(i,j,k))  !!! g6hh2
#endif
#ifdef limitDIC2
          dic(i,j,k)=dmax1(1837d0,1.002*dic(i,j,k))  !!! g6hh3
#endif

         enddo
      end do
      end do


c  Detritus (set to 0 for start up)
      write(6,*)'Detritus...'
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

      do j=j_0,j_1
       do i=i_0,i_1
         do k=1,kdm
#ifdef OBIO_ON_GARYocean
          if (.true.) then
#else
          if (dpinit(i,j,k)/onem .gt. 0.0)then
#endif
           !only detritus components
           tracer(i,j,k,ntyp+n_inert+1) = tracer(i,j,k,1)*0.25*cnratio !as carbon
           tracer(i,j,k,ntyp+n_inert+2) = tracer(i,j,k,3)*0.1
           tracer(i,j,k,ntyp+n_inert+3) = tracer(i,j,k,4)*0.25
           tracer(i,j,k,ntyp+n_inert+1) = 0.0
           tracer(i,j,k,ntyp+n_inert+2) = 0.0
           tracer(i,j,k,ntyp+n_inert+3) = 0.0
          endif
         enddo
        enddo
      enddo

c  Carbon (set to 0 for start up)
c   DIC is derived from GLODAP.  Using mean H from exp601,
c   mean DIC for these values is computed.  Surface DIC is taken
c   as the mean for 020m deeper than the mixed layer, converted from
c   uM/kg to uM
      write(6,*)'Carbon...'
c    conversion from uM to mg/m3
      do j=j_0,j_1
       do i=i_0,i_1
         do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+1) = 0.0
          tracer(i,j,k,ntyp+n_inert+ndet+2) = 0.0
         enddo
       enddo
      enddo

      !only carbon components
      do j=j_0,j_1
       do i=i_0,i_1
         if (ip(i,j)==0) cycle
         do k = 1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+2) = dic(i,j,k)
         enddo
c         car(i,j,k,1) = 3.0  !from Bissett et al 1999 (uM(C))
c         car(i,j,k,1) = 0.0  !from Walsh et al 1999
        enddo
      enddo

c  Light saturation data
      avgq = 0.0

      do j=j_0,j_1
        do i=i_0,i_1
          if (ip(i,j)==0) cycle
          do k=1,kdm
            avgq(i,j,k) = 25.0
          enddo
        enddo
      enddo

c  Coccolithophore max growth rate
      do j=j_0,j_1
       do k=1,kdm
        do i=i_0,i_1
          gcmax(i,j,k) = 0.0
        enddo
       enddo
      enddo
 
      return
      end subroutine obio_bioinit

      subroutine init_alk(alk)
#ifdef OBIO_ON_GARYocean
      use oceanr_dim, only: ogrid
      use oceanres, only: kdm=>lmo
      use ocean, only : ze
#else
      use hycom_dim, only: ogrid, kdm
#endif
      implicit none
      real, dimension(ogrid%i_strt:ogrid%i_stop,
     &      ogrid%j_strt:ogrid%j_stop,kdm), intent(out) :: alk
      integer :: i, j, k

      call bio_inicond('alk_inicond',alk)
      do k=1,kdm
      do j=ogrid%j_strt,ogrid%j_stop
      do i=ogrid%i_strt,ogrid%i_stop
#ifdef OBIO_ON_GARYocean
        if (alk(i,j,k).lt.0.) then
          if (ze(k).le.150.) alk(i,j,k)=2172.      !init neg might be under ice,
          if (ze(k).gt.150. .and. ze(k).lt.1200.)
     &                           alk(i,j,k)=2200.
          if (ze(k).ge.1200.) alk(i,j,k)=2300.
        endif
#else
        if (alk(i,j,k).lt.0.) then
           alk(i,j,k)=0.
         else
           alk(i,j,k)=dmax1(alk(i,j,k),2000.d0)   !set minimum =2000
        endif
#endif
      enddo
      enddo
      enddo
      return
      end subroutine init_alk

c------------------------------------------------------------------------------
      subroutine fndreg(ir)
 
c  Finds nwater indices corresponding to significant regions,
c  and defines arrays.  Variables representative of the regions will 
c  later be kept in these array indicators.
c  Regions are defined as follows:
c        1 -- Antarctic
c        2 -- South Indian
c        3 -- South Pacific
c        4 -- South Atlantic
c        5 -- Equatorial Indian
c        6 -- Equatorial Pacific
c        7 -- Equatorial Atlantic
c        8 -- North Indian
c        9 -- North Central Pacific
c       10 -- North Central Atlantic
c       11 -- North Pacific
c       12 -- North Atlantic
c       13 -- Mediterranean/Black Seas
 

      USE obio_dim
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid   
      Use OCEAN,      only : ZOE=>ZE, oLON_DG,oLAT_DG,ip=>focean
#else
      USE hycom_dim, only : ip,ogrid
      USE hycom_arrays, only : lonij,latij,dpinit
      USE hycom_scalars, only: onem
#endif

      implicit none

      integer, DIMENSION(ogrid%i_strt:ogrid%i_stop,
     &             ogrid%j_strt:ogrid%j_stop), intent(out)   :: ir

      integer i,j,l
      integer iant,isin,ispc,isat,iein,iepc,ieat,incp
     .       ,inca,inat,imed,inin,inpc,nr,ntot

      real rlat,rlon

      integer nir(nrg)

      real antlat,rnpolat
      data antlat,rnpolat /-40.0, 40.0/
 

c  Set up indicators
      iant = 1   !antarcic region
      isin = 2   !south indian ocean
      ispc = 3   !south pacific
      isat = 4   !south atlantic
      iein = 5   !equatorial indian ocean
      iepc = 6   !equatorial pacific ocean
      ieat = 7   !equatorial atlantic ocean
      inin = 8   !north indian ocean
      incp = 9   !north-central pacific
      inca = 10  !north-central atlantic
      inpc = 11  !north pacific
      inat = 12  !north atlantic
      imed = 13  !mediterranean/black sea
 
c  Initialize region indicator array
      ir = 0
      nir = 0
 
c  Find nwater values corresponding to regions
       do j=ogrid%j_strt,ogrid%j_stop
       do i=ogrid%i_strt,ogrid%i_stop

#ifdef OBIO_ON_GARYocean
        rlon=oLON_DG(i,1)
        rlat=oLAT_DG(j,1)
#else
        rlon=lonij(i,j,3)
        rlat=latij(i,j,3)
#endif

        if (rlon .gt. 180)rlon = rlon-360.0

        if (ip(i,j)==0) cycle
#ifndef OBIO_ON_GARYocean
        if (dpinit(i,j,1)/onem .le. 0.0) cycle
#endif

c   Antarctic region
        if (rlat .le. antlat)then
         ir(i,j) = iant
         nir(ir(i,j)) = nir(ir(i,j))+1
        endif

c   South Indian region
        if (rlat .le. -30.0 .and. rlat .gt. antlat)then
         if (rlon .gt. 20.0 .and. rlon .lt. 150.0)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 142.5)then
          ir(i,j) = isin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   South Pacific region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .le. -70.0 .and. rlon .ge. -180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .ge. 150.0 .and. rlon .lt. 180.0)then
          ir(i,j) = ispc
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -10.0 .and. rlat .gt. -30.0)then
         if (rlon .ge. 142.5 .and. rlon .lt. 180.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ispc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   South Atlantic region
        if (rlat .le. -10.0 .and. rlat .gt. antlat)then
         if (rlon .gt. -70.0 .and. rlon .le. 20.0)then
          ir(i,j) = isat
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Indian Ocean region
        if (rlat .le. -8.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 142.5)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -6.0 .and. rlat .gt. -8.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 108.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iein
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif
        if (rlat .le. -4.0 .and. rlat .gt. -6.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 105.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. -2.0 .and. rlat .gt. -4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 2.0 .and. rlat .gt. -2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 104.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 4.0 .and. rlat .gt. 2.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 103.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .le. 7.0 .and. rlat .gt. 4.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 102.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .lt. 10.0 .and. rlat .gt. 7.0)then
         if (rlon .gt. 20.0 .and. rlon .le. 99.0)then
          ir(i,j) = iein
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Equatorial Pacific region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .gt. 90.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
         if (rlon .gt. -180.0 .and. rlon .lt. -70.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = iepc
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   Equatorial Atlantic region
        if (rlat .lt. 10.0 .and. rlat .gt. -10.0)then
         if (rlon .lt. 20.0 .and. rlon .ge. -74.0)then
          if (ir(i,j) .eq. 0)then
           ir(i,j) = ieat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

c   North Indian Ocean
        if (rlat .le. 30.0 .and. rlat .ge. 10.0)then
         if (rlon .gt. 20.0 .and. rlon .lt. 99.0)then
          ir(i,j) = inin
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Pacific region
        if (rlat .ge. 10.0 .and. rlat .le. rnpolat)then
         if (rlon .gt. 99.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
         if (rlon .le. -100.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 18.0)then
         if (rlon .gt. -100.0 .and. rlon .le. -90.0)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .ge. 10.0 .and. rlat .lt. 14.0)then
         if (rlon .gt. -90.0 .and. rlon .le. -84.5)then
          ir(i,j) = incp
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   Mediterranean/Black Seas region
        if (rlat .ge. 30.0 .and. rlat .le. 43.0)then
         if (rlon .ge. -5.5 .and. rlon .lt. 60.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif
        if (rlat .gt. 43.0 .and. rlat .le. 48.0)then
         if (rlon .ge. 0.0 .and. rlon .lt. 70.0)then
          ir(i,j) = imed
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c   North Central Atlantic region
        if (rlat .ge. 7.0 .and. rlat .le. 40.0)then !pickup Carib.Sea
         if (ir(i,j) .eq. 0)then
          ir(i,j) = inca
          nir(ir(i,j)) = nir(ir(i,j))+1
         endif
        endif

c  North Pacific and Atlantic
        if (rlat .gt. rnpolat)then
c  North Pacific
         if (rlon .lt. -105.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else if (rlon .gt. 120.0)then
          ir(i,j) = inpc
          nir(ir(i,j)) = nir(ir(i,j))+1
         else
c  North Atlantic
          if (ir(i,j) .ne. imed)then
           ir(i,j) = inat
           nir(ir(i,j)) = nir(ir(i,j))+1
          endif
         endif
        endif

!       write(121,'(2(i4,1x),4(f8.3,1x),e12.4,i4)')
!    .        i,j,lonij(i,j,3),latij(i,j,3),
!    .        rlon,rlat,dpinit(i,j,1)/onem,ir(i,j)
       end do
       end do
 
c  Set nir to minimum 1 value to prevent error in division
      do nr = 1,nrg
       nir(nr) = max(nir(nr),1)
      enddo
 
c  Total up points for check
      ntot = 0
      do nr = 1,nrg
       ntot = ntot + nir(nr)
       write(6,*)'Region, no. points = ',nr,nir(nr)
      enddo
      write(6,*)'Total ocean points = ',ntot
 
      return
      end subroutine fndreg
c------------------------------------------------------------------------------
      subroutine bio_inicond(filename,fldo2)
      use bio_inicond_mod, only: bio_inicond_read
#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : idm=>imo, jdm=>jmo, kdm=>lmo
      USE OCEANR_DIM, only : ogrid
      USE OCEAN, only : DLATM,ZOE=>ZE, focean, lmm
#else
      use hycom_dim, only: aj_0h,aj_1h,kdm
      use hycom_dim, only: ip,ogrid,iia,jja,idm
      USE GEOM, only : DLATM
      use hycom_arrays, only: dpinit
      USE hycom_scalars, only: onem
      USE hycom_cpler, only: flxa2o
#endif
      implicit none

      character(len=*), intent(in) :: filename
      real, intent(out) :: fldo2(ogrid%i_strt:ogrid%i_stop,
     &    ogrid%j_strt:ogrid%j_stop, kdm)
      integer, parameter :: kgrd=33
#ifdef OBIO_ON_GARYocean
      real fldo(idm,jdm,kgrd)
      interface
        Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)
        Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM) 
        Real*8, optional :: RZ(LM)
        end Subroutine VLKtoLZ
      end interface
#else
      real data2(iia,jja,kgrd)
      real fldo(idm,ogrid%j_strt:ogrid%j_stop,kgrd)
      real pinit(kdm+1)
#endif
      real nodc_depths(kgrd),nodc_d(kgrd+1)
      data nodc_depths/0,  10,  20,  30,  50,  75, 100, 125, 150, 200,
     .          250, 300, 400, 500, 600, 700, 800, 900,1000,1100,1200,
     .    1300,1400,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500/
      integer :: i, j, k, kmax, nodc_kmax

#ifdef OBIO_ON_GARYocean
      call bio_inicond_read(filename, dlatm, 180d0, .true., fldo)
      
      fldo2=-9999.d0
      do j=ogrid%j_strt, ogrid%j_stop
      do i=ogrid%i_strt, ogrid%i_stop
        IF (FOCEAN(i,j).gt.0) then
          call VLKtoLZ(kgrd,lmm(i,j),nodc_depths,ZOE,
     &                             fldo(i,j,:),fldo2(i,j,:))
        ENDIF
      enddo
      enddo
#else
      call bio_inicond_read(filename, dlatm, 180d0, .true., data2)

      fldo=0.
      do k=1,kgrd
        call flxa2o(data2(:,aj_0h:aj_1h,k),fldo(:,:,k))
      end do

      !--------------------------------------------------------
      !use dpinit(i,j,k)/onem

       fldo2=-9999.d0
       do j=ogrid%j_strt,ogrid%j_stop
       do i=ogrid%i_strt,ogrid%i_stop
         if (ip(i,j)==0) cycle

         pinit=0
         do k=1,kdm
           pinit(k+1)=pinit(k)+dpinit(i,j,k)/onem
         end do
          !match nodc and model bottom pressure
          !the top match already
          kmax=1
          do k=2,kdm+1
           if (pinit(k) .gt. pinit(k-1)) kmax=k
          enddo

          !model bottom at kmax+1
          do k=1,kgrd
           if (nodc_depths(k) .le. pinit(min(21,kmax+1))) then
               nodc_d(k)=nodc_depths(k)
               nodc_kmax=k
           endif
cdiag      write(*,'(a,3i5,2e12.4,i5)')'bioinit: ',
cdiag.               i,j,k,fldo(i,j,k),nodc_d(k),nodc_kmax
          enddo
          nodc_d(nodc_kmax+1)=pinit(min(21,kmax+1))

!        call remap1d_pcm(fldo(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
!    .             fldo2(i,j,:),pinit(i,j,:),kdm,.false.,i,j)

         call remap1d_plm(fldo(i,j,1:nodc_kmax),nodc_d,nodc_kmax,
     .             fldo2(i,j,1:kdm),pinit(1:kdm+1),kdm,.false.,i,j)

       enddo
       enddo
#endif
      !--------------------------------------------------------

      end subroutine bio_inicond

      subroutine remap1d_plm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PLM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko,n
      real colin,clout,slope,wgta,wgtb,wgtc,yinteg,
     .     ylft(kold),yrgt(kold),xlo,xhi,ra,rb,ya,yb,q,plmslp,
     .     yrka,ylk,yrk,ylkb
      external plmslp
      logical at_top
      real,parameter    :: onemu=1.e-6, acurcy=1.e-6, flag=-999.
      integer,parameter :: iter=2
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.2))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      do 3 k=1,kold
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- replace each flat segment of stairstep curve by
c --- a slanting segment, using PLM-type limiters.
c
      ylft(   1)=yold(   1)
      yrgt(   1)=yold(   1)
      ylft(kold)=yold(kold)
      yrgt(kold)=yold(kold)
      do 6 n=1,iter                        !  iterate to optimize limiters
      do 2 k=2,kold-1
      if (n.eq.1) then
        yrka=yold(k-1)
        ylk=yold(k)
        yrk=yold(k)
        ylkb=yold(k+1)
      else
        yrka=yrgt(k-1)
        ylk=ylft(k)
        yrk=yrgt(k)
        ylkb=ylft(k+1)
      end if
      wgta=max(onemu,xold(k  )-xold(k-1))
      wgtb=max(onemu,xold(k+1)-xold(k  ))
      wgtc=max(onemu,xold(k+2)-xold(k+1))
      if (k.eq.     1) wgta=onemu
      if (k.eq.kold-1) wgtc=onemu
      slope=plmslp((wgtb*yrka+wgta*ylk)/(wgtb+wgta),
     .     yold(k),(wgtb*ylkb+wgtc*yrk)/(wgtb+wgtc))
      ylft(k)=yold(k)-slope
 2    yrgt(k)=yold(k)+slope
      if (vrbos) print '(8x,a,12x,a,5x,a/(i3,f9.2,5x,2f9.2))',
     .  'y','ylft','yrgt',(ko,yold(ko),ylft(ko),yrgt(ko),ko=1,kold)
 6    continue
      if (vrbos) print '(a/(10f8.2))',
     .  '  target values:',(ynew(ko),ko=1,knew)
c
c --- y in k-th interval now varies from ylft at xold(k) to yrgt at xold(k+1).
c --- find ynew(k) by requiring
c --- that the integral over y*dx from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
ccc      if (vrbos) print '(a,2f9.3)','xlo,xhi =',xlo,xhi
      if (xhi.gt.xlo) then
        at_top=.false.
        do 5 ko=1,kold
        if (xold(ko).ge.xhi) go to 1
c --- integrate over sloping portions of y(x) curve:
        ra=max(xlo,min(xhi,xold(ko  )))
        rb=max(xlo,min(xhi,xold(ko+1)))
cnat    ya=ylft(k)
cnat    yb=yrgt(k)
        ya=ylft(ko)
        yb=yrgt(ko)
        wgta=flag
        wgtb=flag
        if (xold(ko+1).ne.xold(ko)) then
          if (ra.ge.xold(ko).and.ra.le.xold(ko+1)) then
            wgta=(xold(ko+1)-ra)/(xold(ko+1)-xold(ko))
            ya=ylft(ko)*wgta+yrgt(ko)*(1.-wgta)
          end if
          if (rb.ge.xold(ko).and.rb.le.xold(ko+1)) then
            wgtb=(xold(ko+1)-rb)/(xold(ko+1)-xold(ko))
            yb=ylft(ko)*wgtb+yrgt(ko)*(1.-wgtb)
          end if
        end if
        yinteg=yinteg+.5*(ya+yb)*(rb-ra)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,ko,ra,rb,wgta,wgtb,ya,yb,yinteg
 5      continue
        yinteg=yinteg+yb*(xhi-rb)
ccc        if (vrbos) print '(2i4,4f9.3,3f11.1)',
ccc     .    k,0,rb,xhi,wgta,wgtb,yb,yb,yinteg
 1      ynew(k)=yinteg/(xhi-xlo)
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
ccc      if (vrbos) print '(a,f11.1)','ynew =',ynew(k)
      clout=clout+ynew(k)*(xnew(k+1)-xnew(k))
 4    continue
c
      if (abs(clout-colin).gt.acurcy*10.*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)

      return
      end
c
c
      function plmslp(ylft,ymid,yrgt)
c
c --- get slope at point 'ymid' for piecewise linear interpolation
      if (ymid.le.min(ylft,yrgt) .or.
     .    ymid.ge.max(ylft,yrgt)) then
        plmslp=0.
      else if ((yrgt-ylft)*(ylft+yrgt-2.*ymid).gt.0.) then
        plmslp=ymid-ylft
      else
        plmslp=yrgt-ymid
      end if
      return
      end

      subroutine remap1d_pcm(yold,xold,kold,ynew,xnew,knew,vrbos,i,j)
c
c --- consider two stepwise constant functions -yold,ynew- whose
c --- discontinuities are at abscissa values -xold,xnew- respectively.
c --- treat -ynew- as unknown. solve for -ynew- under the condition that
c --- the integral over y*dx is preserved (integration based on PCM).
c
      implicit none
      integer,intent(IN) :: kold,knew
      integer,intent(IN) :: i,j                !  current location in horiz.grid
      real,intent(IN)    :: yold(kold),xold(kold+1),xnew(knew+1)
      real,intent(OUT)   :: ynew(knew)
      logical,intent(IN) :: vrbos        !  if true, print diagnostics
      integer k,ko
      real colin,clout,colmx,yinteg,xlo,xhi,xa,xb
      logical at_top
      real,parameter    :: acurcy=1.e-6
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- old profile:     x         y',
     .  (k,xold(k),yold(k),k=1,kold),kold+1,xold(kold+1)
 101  format (2i5,a/(i30,f12.1,f10.3))
c
c --- column integrals (colin/clout) are computed for diagnostic purposes only
      colin=0.
      clout=0.
      colmx=0.
      do 3 k=1,kold
      colmx=max(colmx,abs(yold(k)))
 3    colin=colin+yold(k)*(xold(k+1)-xold(k))
c
c --- find ynew(k) by requiring that the integral over y*dx
c --- from xnew(k) to xnew(k+1) be preserved.
c
      at_top=.true.
      do 4 k=1,knew
      yinteg=0.
      xlo=xnew(k  )
      xhi=xnew(k+1)
      if (xhi.gt.xlo) then
        at_top=.false.
        xb=xlo
        do 5 ko=1,kold
        xa=xb
        xb=min(xhi,max(xlo,xold(ko+1)))
        yinteg=yinteg+yold(ko)*(xb-xa)
        if (xa.ge.xhi) go to 1
 5      continue
        xa=xb
        xb=xhi
        yinteg=yinteg+yold(kold)*(xb-xa)
 1      ynew(k)=yinteg/(xhi-xlo)
        clout=clout+yinteg
      else if (at_top) then
        ynew(k)=yold(   1)
      else                              !  at end
        ynew(k)=yold(kold)
      end if
 4    continue
c
      if (abs(clout-colin).gt.acurcy*colmx*xold(kold+1))
     .  write (*,100) i,j,' remap1d - column intgl.error',
     .   colin,clout,(clout-colin)/colin
 100  format (2i5,a,2es14.6,es9.1)
c
      if (vrbos)
     . write (*,101) i,j,' remap1d -- new profile:     x         y',
     .  (k,xnew(k),ynew(k),k=1,knew),knew+1,xnew(knew+1)
      return
      end

#ifndef OBIO_ON_GARYocean
      Subroutine HNTR80 (IMA,JMA,OFFIA,DLATA,
     *                   IMB,JMB,OFFIB,DLATB, DATMIS)
C****
C**** HNTR80 fills in the common block HNTRCB with coordinate
C**** parameters that will be used by subsequent calls to HNTR8.
C**** The 5 Real input values are expected to be Real*8.
C****
C**** Input: IMA = number of cells in east-west direction of grid A
C****        JMA = number of cells in north-south direction of grid A
C****      OFFIA = number of cells of grid A in east-west direction
C****              from IDL (180) to western edge of cell IA=1
C****      DLATA = minutes of latitude for non-polar cells on grid A
C****        IMB = number of cells in east-west direction of grid B
C****        JMB = number of cells in north-south direction of grid B
C****      OFFIB = number of cells of grid B in east-west direction
C****              from IDL (180) to western edge of cell IB=1
C****      DLATB = minutes of latitude for non-polar cells on grid B
C****     DATMIS = missing data value inserted in output array B when
C****              cell (IB,JB) has integrated value 0 of WTA
C****
C**** Output: common block /HNTRCB/
C**** SINA(JA) = sine of latitude of northern edge of cell JA on grid A
C**** SINB(JB) = sine of latitude of northern edge of cell JB on grid B
C**** FMIN(IB) = fraction of cell IMIN(IB) on grid A west of cell IB
C**** FMAX(IB) = fraction of cell IMAX(IB) on grid A east of cell IB
C**** GMIN(JB) = fraction of cell JMIN(JB) on grid A south of cell JB
C**** GMAX(JB) = fraction of cell JMAX(JB) on grid A north of cell JB
C**** IMIN(IB) = western most cell of grid A that intersects cell IB
C**** IMAX(IB) = eastern most cell of grid A that intersects cell IB
C**** JMIN(JB) = southern most cell of grid A that intersects cell JB
C**** JMAX(JB) = northern most cell of grid A that intersects cell JB
C****
      Implicit Real*8 (A-H,O-Z)
      Parameter (TWOPI=6.283185307179586477d0)
      Real*8 OFFIA,DLATA, OFFIB,DLATB, DATMIS,DATMCB
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMCB, INA,JNA, INB,JNB
C****
      INA = IMA  ;  JNA = JMA
      INB = IMB  ;  JNB = JMB
      DATMCB = DATMIS
      If (IMA<1 .or. IMA>10800 .or. JMA<1 .or. JMA>5401 .or.
     *    IMB<1 .or. IMB>10800 .or. JMB<1 .or. JMB>5401)  GoTo 400
C****
C**** Partitions in east-west (I) direction
C**** Domain, around the globe, is scaled to fit from 0 to IMA*IMB
C****
      DIA = IMB  !  width of single A grid cell in scaled domain
      DIB = IMA  !  width of single B grid cell in scaled domain
      IA  = 1
      RIA = (IA+OFFIA - IMA)*IMB  !  scaled longitude of eastern edge
      IB  = IMB
      Do 150 IBp1=1,IMB
      RIB = (IBp1-1+OFFIB)*IMA    !  scaled longitude of eastern edge
  110 If (RIA-RIB)  120,130,140
  120 IA  = IA  + 1
      RIA = RIA + DIA
      GoTo 110
C**** Eastern edges of cells IA of grid A and IB of grid B coincide
  130 IMAX(IB) = IA
      FMAX(IB) = 0
      IA  = IA  + 1
      RIA = RIA + DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 0
      GoTo 150
C**** Cell IA of grid A contains western edge of cell IB of grid B
  140 IMAX(IB) = IA
      FMAX(IB) = (RIA-RIB)/DIA
      IMIN(IBp1) = IA
      FMIN(IBp1) = 1-FMAX(IB)
  150 IB = IBp1
      IMAX(IMB) = IMAX(IMB) + IMA
C       WRITE (0,915) 'IMIN=',IMIN(1:IMB)
C       WRITE (0,915) 'IMAX=',IMAX(1:IMB)
C       WRITE (0,916) 'FMIN=',FMIN(1:IMB)
C       WRITE (0,916) 'FMAX=',FMAX(1:IMB)
C****
C**** Partitions in the north-south (J) direction
C**** Domain is measured in minutes (1/60-th of a degree)
C****
      FJEQA = .5*(1+JMA)
      Do 210 JA=1,JMA-1
      RJA = (JA+.5-FJEQA)*DLATA  !  latitude in minutes of northern edge
  210 SINA(JA) = Sin (RJA*TWOPI/(360*60))
      SINA(0)  = -1
      SINA(JMA)=  1
C****
      FJEQB = .5*(1+JMB)
      Do 220 JB=1,JMB-1
      RJB = (JB+.5-FJEQB)*DLATB  !  latitude in minutes of northern edge
  220 SINB(JB) = Sin (RJB*TWOPI/(360*60))
      SINB(0)  = -1
      SINB(JMB)=  1
C****
      JMIN(1) = 1
      GMIN(1) = 0
      JA = 1
      Do 350 JB=1,JMB-1
  310 If (SINA(JA)-SINB(JB))  320,330,340
  320 JA = JA + 1
      GoTo 310
C**** Northern edges of cells JA of grid A and JB of grid B coincide
  330 JMAX(JB) = JA
      GMAX(JB) = 0
      JA = JA + 1
      JMIN(JB+1) = JA
      GMIN(JB+1) = 0
      GoTo 350
C**** Cell JA of grid A contains northern edge of cell JB of grid B
  340 JMAX(JB) = JA
      GMAX(JB) = SINA(JA) - SINB(JB)
      JMIN(JB+1) = JA
      GMIN(JB+1) = SINB(JB) - SINA(JA-1)
  350 Continue
      JMAX(JMB) = JMA
      GMAX(JMB) = 0
C       WRITE (0,915) 'JMIN=',JMIN(1:JMB)
C       WRITE (0,915) 'JMAX=',JMAX(1:JMB)
C       WRITE (0,916) 'GMIN=',GMIN(1:JMB)
C       WRITE (0,916) 'GMAX=',GMAX(1:JMB)
      Return
C****
C**** Invalid parameters or dimensions out of range
C****
  400 Write (0,940) IMA,JMA,OFFIA,DLATA, IMB,JMB,OFFIB,DLATB, DATMIS
      Stop 400
C****
C 915 Format (/ 1X,A5 / (20I6))
C 916 Format (/ 1X,A5 / (20F6.2))
  940 Format ('0Arguments received by HNTRP0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATA   = minutes of latitude for interior grid cell'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid cells from',
     *                    ' IDL to western edge of grid cell I=1'/
     *  E24.8,' = DLATB   = minute of latitude for interior grid cell'/
     *  E24.8,' = DATMIS  = missing data value to be put in B array',
     *                    ' when integrated WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      End Subroutine HNTR80

      Subroutine HNTR8P (WTA,A,B)
C****
C**** HNTR8P is similar to HNTR8 but polar values are replaced by
C**** their longitudinal mean.
C**** The 3 Real input values are expected to be Real*8.
C****
      Implicit Real*8 (A-H,O-Z)
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
      Call HNTR8 (WTA,A,B)
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      Do 40 JB=1,JMB,JMB-1
      BMEAN  = DATMIS
      WEIGHT = 0
      VALUE  = 0
      Do 10 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
      If (B(IJB) == DATMIS)  GoTo 20
      WEIGHT = WEIGHT + 1
      VALUE  = VALUE  + B(IJB)
   10 Continue
      If (WEIGHT.ne.0)  BMEAN = VALUE/WEIGHT
   20 Do 30 IB=1,IMB
      IJB  = IB + IMB*(JB-1)
   30 B(IJB) = BMEAN
   40 Continue
      Return
      End Subroutine HNTR8P


      Subroutine HNTR8 (WTA,A,B)
C**** 
C**** HNTR8 performs a horizontal interpolation of per unit area or per
C**** unit mass quantities defined on grid A, calculating the quantity
C**** on grid B.  B grid values that cannot be calculated because the
C**** covering A grid boxes have WTA = 0, are set to the value DATMIS.
C**** The area weighted integral of the quantity is conserved.
C**** The 3 Real input values are expected to be Real*8.
C****
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit area or per unit mass quantity
C**** Output:  B = horizontally interpolated quantity on B grid
C****
      Implicit Real*8 (A-H,O-Z)  
      Real*8 WTA(*), A(*), B(*), DATMIS
      Common /HNTRCB/ SINA(0:5401),SINB(0:5401),
     *       FMIN(10800),FMAX(10800),GMIN(5401),GMAX(5401),
     *       IMIN(10800),IMAX(10800),JMIN(5401),JMAX(5401),
     *       DATMIS, IMA,JMA, IMB,JMB
C****
C**** Interpolate the A grid onto the B grid
C****
      Do 20 JB=1,JMB
      JAMIN = JMIN(JB)
      JAMAX = JMAX(JB)
      Do 20 IB=1,IMB 
      IJB  = IB + IMB*(JB-1)
      WEIGHT= 0
      VALUE = 0 
      IAMIN = IMIN(IB)
      IAMAX = IMAX(IB)
      Do 10 JA=JAMIN,JAMAX
      G = SINA(JA)-SINA(JA-1)
      If (JA==JAMIN)  G = G - GMIN(JB)
      If (JA==JAMAX)  G = G - GMAX(JB)
      Do 10 IAREV=IAMIN,IAMAX
      IA  = 1 + Mod(IAREV-1,IMA)
      IJA = IA + IMA*(JA-1)
      F   = 1
      If (IAREV==IAMIN)  F = F - FMIN(IB)
      If (IAREV==IAMAX)  F = F - FMAX(IB)
      WEIGHT = WEIGHT + F*G*WTA(IJA)
   10 VALUE  = VALUE  + F*G*WTA(IJA)*A(IJA)
      B(IJB) = DATMIS 
      If (WEIGHT.ne.0)  B(IJB) = VALUE/WEIGHT
   20 Continue
      Return
      End Subroutine HNTR8
#endif

      Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)        !  2008/02/05
C****
C**** VLKtoLZ assumes a continuous piecewise linear tracer distribution,
C**** defined by input tracer concentrations RK at KM specific points.
C**** MK in the downward vertical mass coordinate.
C**** R(M) = {RK(K-1)*[MK(K)-M] + RK(K)*[M-MK(K-1)]} / [MK(K)-MK(K-1)]
C****               when MK(K-1) < M < MK(K).
C**** R(M) = RK(1)  when M < MK(1).
C**** R(M) = is undefined when MK(KM) < M.
C****
C**** VLKtoLZ integrates this tracer distribution over the LM output
C**** layers defined by their layer edges ME, calculating the tracer
C**** mass RM of each layer and the vertical gradient RZ.
C**** RNEW(M) = RL(L) + RZ(L)*[M-MC(L)]/dM(L) when ME(L-1) < M < ME(L)
C**** where MC(L) = .5*[ME(L-1)+ME(L)] and dM(L) = ME(L)-ME(L-1)
C**** Mean concentration of output layers is RL(L) = RM(L)/dM(L).
C****
C**** If ME(L-1) < MK(KM) < ME(L), then RL(L) and RZ(L) are calculated
C**** from the input profile up to MK(KM); RL(L+1:LM) and RZ(L+1:LM)
C**** for deeper layers are undefined, set to DATMIS.
C****
C**** Input:  KM = number of input edges
C****         LM = number of output cells
C****         MK = mass coordinates of input points (kg/m^2)
C****         ME = mass coordinates of output layer edges (kg/m^2)
C****         RK = tracer concentration at input points
C****
C**** Output: RL = mean tracer concentration of each output layer
C****         RZ = vertical gradient of tracer mass of each output layer
C****
C**** Internal: RM = integrated tracer mass of output layers (kg/m^2)
C****           RQ = integrated tracer mass times mass (kg^2/m^4)
C****
      Implicit Real*8 (A-H,M-Z)
      Parameter (DATMIS = -999999)
      Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM),RM(1024),RQ(1024)
      Real*8, optional :: RZ(LM)

C     If (LM > 1024)  Stop 'LM exceeds internal dimentions in VLKtoLZ'
C****
      RM(1:LM) = 0
      RQ(1:LM) = 0
      K = 1
      L = 1
      MC = .5*(ME(L)+ME(L-1))
C****
C**** Integrate layers with M < MK(1)
C****
      If (ME(0) < MK(1))  GoTo 20
C**** MK(1) <= ME(0), determine K such that MK(K-1) <= ME(0) < MK(K)
   10 If (K == KM)  GoTo 200  ;  K = K+1
      If (MK(K) <= ME(0))  GoTo 10
      GoTo 130  !  MK(K-1) <= ME(0) < MK(K)
C**** ME(0) < MK(1), determine output cell containing MK(1)
   20 If (MK(1) < ME(L))  GoTo 30
C**** ME(L-1) < ME(L) < MK(1), integrate RM from ME(L-1) to ME(L)
      RM(L) = RK(1)*(ME(L)-ME(L-1))
      RQ(L) = 0
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      GoTo 20
C**** ME(L-1) < MK(1) < ME(L), integrate RM from ME(L-1) to MK(1)
   30 RM(L) = RK(1)*(MK(1)-ME(L-1))
      RQ(L) = RK(1)*(MK(1)-ME(L-1))*(.5*(MK(1)+ME(L-1))-MC)
      If (K == KM)  GoTo 220  ;  K = K+1
C****
C**** Integrate layers with MK(1) < M < MK(KM)
C****
  100 If (ME(L) < MK(K))  GoTo 120
C**** ME(L-1) < MK(K-1) < MK(K) < ME(L), integrate from MK(K-1) to MK(K)
      RM(L) = RM(L) + (RK(K)-RK(K-1))*(MK(K)+MK(K-1))/2 +
     +                RK(K-1)*MK(K)-RK(K)*MK(K-1)
      RQ(L) = RQ(L) +
     +  (RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +  (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+MK(K-1))/2 +
     +  (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C**** ME(L-1) < MK(K-1) < ME(L) < MK(K), integrate from MK(K-1) to ME(L)
  120 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(ME(L)+MK(K-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-MK(K-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+MK(K-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-MK(K-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
  130 If (MK(K) < ME(L))  GoTo 160
C**** MK(K-1) < ME(L-1) < ME(L) < MK(K), integrate from ME(L-1) to ME(L)
  140 RM(L) = ((RK(K)-RK(K-1))*(ME(L)+ME(L-1))/2 +
     +         (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-ME(L-1)) /
     /        (MK(K)-MK(K-1))
      RQ(L) =
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      If (ME(L) < MK(K))  GoTo 140
C**** MK(K-1) < ME(L-1) < MK(K) < ME(L), integrate from ME(L-1) to MK(K)
  160 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(MK(K)+ME(L-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (MK(K)-ME(L-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (MK(K)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C****
C**** Calculate RL and RZ from RM and RQ when MK(KM) < ME(LM)
C****
C**** MK(KM) <= ME(0)
  200 RL(:) = DATMIS
      if (present(rz)) RZ(:) = DATMIS
      Return
C**** ME(L-1) < MK(KM) < ME(L)
  220 RL(1:L-1) =   RM(1:L-1) / (ME(1:L-1)-ME(0:L-2))
      RL(L)  =   RM(L)  / (MK(KM)-ME(L-1))
      RL(L+1:LM) = DATMIS
      if (present(rz)) then
        RZ(1:L-1) = 6*RQ(1:L-1) / (ME(1:L-1)-ME(0:L-2))**2
        RZ(L)= 6*(RQ(L) + .5*(ME(L)-MK(KM))*RM(L)) / (MK(KM)-ME(L-1))**2
C**** Vertical gradient is extrapolated half way to .5*[MK(KM)+ME(L)]
        RZ(L)  = RZ(L) * (.5*(MK(KM)+ME(L))-ME(L-1)) / (MK(KM)-ME(L-1))
        RZ(L+1:LM) = DATMIS
      endif
      Return
C****
C**** Calculate RL and RZ from RM and RQ when ME(LM) < MK(KM)
C****
  300 RL(1:lm) =   RM(1:lm) / (ME(1:lm)-ME(0:lm-1))
      if (present(rz)) RZ(1:lm) = 6*RQ(1:lm) / (ME(1:lm)-ME(0:lm-1))**2
      End
