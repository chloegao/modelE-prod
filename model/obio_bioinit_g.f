#include "rundeck_opts.h"

#ifdef OBIO_ON_GARYocean
      subroutine obio_bioinit_g
!note: 
!obio_bioinit is called only for a cold start and reads in 
!INITIAL conditions and interpolates them
!to the ocean grid. Such fields are nitrates,silicate,dic
!obio_init  is called for every start of the run and reads 
!in BOUNDARY conditions and interpolates 
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
      USE obio_com, only: gcmax,tracer

      USE OCEANRES, only : kdm=>lmo,dzo
      USE OCEAN, only : ZOE=>ZE
      USE OCEANR_DIM, only : ogrid

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
      integer i1,i2

      integer nir(nrg),nt

      INTEGER :: j_0h,j_1h,i_0h,i_1h

      integer, ALLOCATABLE, DIMENSION(:,:)   :: ir
      real,  ALLOCATABLE, DIMENSION(:,:,:) :: Fer,dicmod,dic

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO


      ALLOCATE(ir(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(Fer(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(dicmod(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(dic(i_0h:i_1h,j_0h:j_1h,kdm))

c  Initialize

      tracer(:,:,:,1:ntyp)=0.

      Fer(:,:,:) = 0.d0


      call bio_inicond_g('nitrates_inicond',tracer(:,:,:,1)) !because 1 is nintrate
      call bio_inicond_g('silicate_inicond',tracer(:,:,:,3)) !because 3 is silicate

#ifdef obio_TRANSIENTRUNS
! in the transient runs keep the dic read in from RSF/AIC file
      dic = tracer(:,:,:,15)
#else
! otherwise take from rundeck
      call bio_inicond_g('dic_inicond',dic)
#endif

#ifdef TRACERS_Alkalinity
      call init_alk(tracer(:,:,:,ntrac))
#endif

!!these rno3 and so2_init fields are not correct. There are void points due to
!mismatch of the noaa grid and the ocean grid. To fill in have to do
!the interpolation in matlab (furtuna).
!at the same time use dps and interpolate to layer depths from the model

      do k=1,kdm
      do j=j_0h,j_1h
       do i=i_0h,i_1h
         if(tracer(i,j,k,1).le.0.d0)tracer(i,j,k,1)=0.085d0
         if(tracer(i,j,k,3).le.0.d0)tracer(i,j,k,3)=0.297d0
          if (dic(i,j,k).le.0.d0) dic(i,j,k)=1837d0
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

      do j=j_0h,j_1h
        do i=i_0h,i_1h
          fer(i,j,:)=fer_values(ir(i,j))
        enddo  
      enddo  
 
c  Create arrays 
      write(6,*)'Creating bio restart data for ',ntyp,' arrays and'
     . ,kdm,'  layers...'

      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm


          !Nitrate
          !read earlier from file

          !Ammonium
          tracer(i,j,k,2) = 0.5
          !!!if (ZOE(k) .gt. 4000.d0) tracer(i,j,k,2) = Pdeep(2)

          !Silica
          !read earlier from file

          !Iron
          !initialize iron distribution in the ocean
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
          !!!if (ZOE(k) .gt. 4000.d0) tracer(i,j,k,4) = Pdeep(4)

!         write(*,'(a,3i5,e12.4,i5,2e12.4)')'obio_bioinit, iron:',
!    .       i,j,k,Fer(i,j,k),ir(i,j),tracer(i,j,k,4),tracer(i,j,k,1)

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
          dicmod(i,j,k)=dic(i,j,k)

#ifdef limitDIC1
!!!       dic(i,j,k)=dmax1(1837d0,0.99*dic(i,j,k))  !!! g6hh
          dic(i,j,k)=dmax1(1837d0,1.005*dic(i,j,k))  !!! g6hh2
#endif
#ifdef limitDIC2
          dic(i,j,k)=dmax1(1837d0,1.002*dic(i,j,k))  !!! g6hh3
#endif
          dicmod(i,j,k)=dic(i,j,k)

      end do
      end do
      end do


c  Detritus (set to 0 for start up)
      write(6,*)'Detritus...'
      cnratio = 106.0/16.0*12.0    !C:N ratio (ugl:uM)
      csratio = 106.0/16.0*12.0    !C:Si ratio (ugl:uM)
      cfratio = 150000.0*12.0*1.0E-3    !C:Fe ratio (ugl:nM)

      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
           !only detritus components
           tracer(i,j,k,ntyp+n_inert+1) = tracer(i,j,k,1)*0.25*cnratio !as carbon
           tracer(i,j,k,ntyp+n_inert+2) = tracer(i,j,k,3)*0.1
           tracer(i,j,k,ntyp+n_inert+3) = tracer(i,j,k,4)*0.25
           tracer(i,j,k,ntyp+n_inert+1) = 0.0
           tracer(i,j,k,ntyp+n_inert+2) = 0.0
           tracer(i,j,k,ntyp+n_inert+3) = 0.0
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

      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+1) = 0.0
          tracer(i,j,k,ntyp+n_inert+ndet+2) = 0.0
      enddo
      enddo
      enddo

      !only carbon components
      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          tracer(i,j,k,ntyp+n_inert+ndet+2) = dicmod(i,j,k)
c         car(i,j,k,1) = 3.0  !from Bissett et al 1999 (uM(C))
c         car(i,j,k,1) = 0.0  !from Walsh et al 1999
      enddo
      enddo
      enddo

c  Light saturation data
      write(6,*)'Light saturation data...'
      avgq = 0.0d0

      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          avgq(i,j,k) = 25.0
      enddo
      enddo
      enddo

c  Coccolithophore max growth rate
      write(6,*)'Coccolithophore max growth rate...'
      do j=j_0h,j_1h
      do i=i_0h,i_1h
      do k=1,kdm
          gcmax(i,j,k) = 0.0
      enddo
      enddo
      enddo
 
      call obio_trint(0)

      return
      end subroutine obio_bioinit_g

      subroutine init_alk(alk)
      use oceanr_dim, only: ogrid
      use oceanres, only: lmo
      use ocean, only : ze
      implicit none
      real, dimension(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &      ogrid%j_strt_halo:ogrid%j_stop_halo,lmo), intent(out) :: alk
      integer :: i, j, k

      call bio_inicond_g('alk_inicond',alk)

      !remove negative values
      !negs are over land or under ice due to GLODAP missing values in the Arctic Ocean
      !for under ice missing values, use climatological minimums for sets of layers based
      !on GLODAP, rather than setting to the same global min.
      do j=ogrid%j_strt_halo,ogrid%j_stop_halo
      do i=ogrid%i_strt_halo,ogrid%i_stop_halo
      do k=1,lmo
       if (alk(i,j,k).lt.0.) then
          if (ze(k).le.150.) alk(i,j,k)=2172.      !init neg might be under ice,
          if (ze(k).gt.150. .and. ze(k).lt.1200.)
     &                           alk(i,j,k)=2200. 
          if (ze(k).ge.1200.) alk(i,j,k)=2300.      
       endif
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
 

      USE OCEANR_DIM, only : ogrid   
      Use OCEAN,      only : ZOE=>ZE, oLON_DG,oLAT_DG,focean

      USE obio_dim
      implicit none


      integer, DIMENSION(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &    ogrid%j_strt_halo:ogrid%j_stop_halo), intent(out)   :: ir
      integer i,j,l,k
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
       do 1000 j=ogrid%j_strt_halo,ogrid%j_stop_halo
       do 1000 i=ogrid%i_strt_halo,ogrid%i_stop_halo

        rlon=oLON_DG(i,1)
        rlat=oLAT_DG(j,1)

        if (rlon .gt. 180)rlon = rlon-360.0

        if (focean(i,j).le.0) go to 100

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
 100    continue

 1000 continue
 
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

      subroutine bio_inicond_g(filename,fldo2)
      use bio_inicond_mod, only: bio_inicond_read
      USE OCEANRES, only : idm=>imo, jdm=>jmo, kdm=>lmo
      USE OCEANR_DIM, only : ogrid
      USE OCEAN, only : DLATM,ZOE=>ZE, focean, lmm
      implicit none

      character(len=*), intent(in) :: filename
      real, intent(out) :: fldo2(ogrid%i_strt_halo:ogrid%i_stop_halo,
     &    ogrid%j_strt_halo:ogrid%j_stop_halo, kdm)
      integer, parameter :: kgrd=33
      real fldo(idm,jdm,kgrd)
      real nodc_depths(kgrd)
      data nodc_depths/0,  10,  20,  30,  50,  75, 100, 125, 150, 200,
     .          250, 300, 400, 500, 600, 700, 800, 900,1000,1100,1200,
     .    1300,1400,1500,1750,2000,2500,3000,3500,4000,4500,5000,5500/
      integer :: i, j, k, lm
      interface
        Subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ)
        Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM) 
        Real*8, optional :: RZ(LM)
        end Subroutine VLKtoLZ
      end interface

      call bio_inicond_read(filename, dlatm, 180d0, .true., fldo)
      
      fldo2=-9999.d0
      do j=ogrid%j_strt_halo, ogrid%j_stop_halo
      do i=ogrid%i_strt_halo, ogrid%i_stop_halo
        IF (FOCEAN(i,j).gt.0) then
          lm=lmm(i,j)
          call VLKtoLZ(kgrd,lm,nodc_depths,ZOE,fldo(i,j,:),fldo2(i,j,:))
        ENDIF
      enddo
      enddo

      return
      end subroutine bio_inicond_g

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

#endif /*  OBIO_ON_GARYocean */
