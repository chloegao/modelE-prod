#include "rundeck_opts.h"

#if (defined CONSTANT_MESO_DIFFUSIVITY)
#define USE_1D_MESODIFF
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
#define USE_1D_MESODIFF
#endif

#if defined(SIMPLE_MESODIFF) || defined(OCN_GISS_MESO) || defined(USE_1D_MESODIFF)
#else
#define ORIG_MESODIFF
#endif

      module ocnmeso_com
      implicit none

!@dbparam kvismult multiplier for Visbeck K when using simple_mesodiff
      real*8 :: kvismult=1d0

!@dbparam kbg (m2/s) minimum mesoscale diffusivity
      real*8 :: kbg=100d0


!@var g3d potential enthalpy (J/kg)
!@var s3d salinity (kg/kg)
!@var v3d specific volume (ref to mid point pressure)
!@var r3d density (ref to mid point pressure)
!@var p3d mid point pressure
!@var rhox along-layer x-gradient of potential density (ref. to local pres)
!@var rhoy along-layer y-gradient of potential density (ref. to local pres)
!@var rhomz minus the z-gradient of potential density (ref. to local pres)
!@var byrhoz 1/rhomz
!@var dze3d approx distance between layer midpoints
!@var bydze3d 1/dze3d
      real*8, allocatable, dimension(:,:,:) ::
     &     g3d,s3d,p3d,r3d,v3d,
     &     rhox,rhoy,rhomz,byrhoz,
     &     dze3d,bydze3d

#if (defined CONSTANT_MESO_DIFFUSIVITY)
!@dbparam meso_diffusivity_const constant diffusivity (m2/s)
!@+       for sensitivity studies
      real*8 :: meso_diffusivity_const
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
!@dbparam meso_diffusivity_z0
!@dbparam meso_diffusivity_zscale
!@dbparam meso_diffusivity_deep
      real*8 ::
     & meso_diffusivity_z0,meso_diffusivity_zscale,meso_diffusivity_deep
#elif (defined CONSTANT_MESO_LENSCALE)
!@dbparam meso_lenscale_const a fixed length scale (meters) to use in
!@+       lieu of Rossby radius RD in AINV = AMU * RD**2 * BYTEADY
      real*8 :: meso_lenscale_const
#endif

      end module ocnmeso_com

      subroutine alloc_ocnmeso_com
      use ocean, only : im,lmo
      use ocnmeso_com
      use dictionary_mod, only : sync_param
      use oceanr_dim, only : ogrid
      use domain_decomp_1d, only : getdomainbounds
      use dictionary_mod
      implicit none
      integer :: j_0h,j_1h

      call getdomainbounds(ogrid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      allocate( g3d (lmo,im,j_0h:j_1h) )
      allocate( s3d (lmo,im,j_0h:j_1h) )
      allocate( p3d (lmo,im,j_0h:j_1h) )
      allocate( r3d (lmo,im,j_0h:j_1h) )
      allocate( v3d (lmo,im,j_0h:j_1h) )
      allocate( rhox   (lmo,im,j_0h:j_1h) )
      allocate( rhoy   (lmo,im,j_0h:j_1h) )

      allocate( rhomz  (im,j_0h:j_1h,lmo) )
      allocate( byrhoz (im,j_0h:j_1h,lmo) )
      allocate( dze3d  (im,j_0h:j_1h,lmo) )
      allocate( bydze3d(im,j_0h:j_1h,lmo) )

#ifdef SIMPLE_MESODIFF
      call sync_param("ocean_kvismult",kvismult)
#endif

#if (defined CONSTANT_MESO_DIFFUSIVITY)
       call get_param('meso_diffusivity_const',meso_diffusivity_const)
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
       call get_param('meso_diffusivity_z0',meso_diffusivity_z0)
       call get_param('meso_diffusivity_zscale',meso_diffusivity_zscale)
       call get_param('meso_diffusivity_deep',meso_diffusivity_deep)
#elif (defined CONSTANT_MESO_LENSCALE)
       call get_param( 'meso_lenscale_const', meso_lenscale_const )
#endif

      end subroutine alloc_ocnmeso_com

      subroutine ocnmeso_drv
      use ocean, only : im,jm,lmo,lmm,mo,g0m,s0m,dts
      use ocean, only :
     *     gxmo,gymo,gzmo,
     *     sxmo,symo,szmo
      use ocean, only :
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use domain_decomp_1d, only : getdomainbounds
      use oceanr_dim, only : grid=>ogrid
      use odiag, only : oijl=>oijl_loc,oij=>oij_loc,
     *    ijl_ggmfl,ijl_sgmfl
     &   ,ijl_mfub,ijl_mfvb,ijl_mfwb
      use odiag, only : ij_gmsc!,ij_gmscz
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
      USE OCEAN, only : trmo,
     &     txmo,tymo,tzmo,txxmo,tyymo,tzzmo,txymo,tyzmo,tzxmo
      Use ODIAG, Only: toijl=>toijl_loc,toijl_gmfl
#endif
      implicit none
      integer i,j,l,n

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d,gmscz
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &   k3d

      logical, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     zeroing_mask
#ifdef OCN_GISS_MESO
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &    g0m0,gxmo0,gymo0,s0m0,sxmo0,symo0
#endif

#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: trentry
#endif

c**** Extract domain decomposition info

      integer :: j_0, j_1, j_0h,j_1h, j_0s,j_1s
      logical :: have_south_pole,have_north_pole

      call getdomainbounds(grid, j_strt = j_0, j_stop = j_1,
     &     j_strt_skp = j_0s, j_stop_skp = j_1s,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h,
     &     have_south_pole = have_south_pole,
     &     have_north_pole = have_north_pole)

C**** Calculate horizontal and vertical density gradients.
      call densgrad

#ifdef OCN_GISS_MESO
      G0M0=G0M; GXMO0=GXMO; GYMO0=GYMO
      S0M0=S0M; SXMO0=SXMO; SYMO0=SYMO
#endif

C**** Apply GM + Redi tracer fluxes
      k3d = 0.
C**** Calculate mesoscale diffusivity
#ifdef USE_1D_MESODIFF
      call get_1d_mesodiff(k3d)
#endif

      gmscz = 0.
#ifdef SIMPLE_MESODIFF
      call simple_mesodiff(k2d)
#endif

#if defined(ORIG_MESODIFF) || defined(OCN_GISS_MESO)
      call orig_mesodiff(k2d,zeroing_mask)
#endif

#if defined(ORIG_MESODIFF) || defined(SIMPLE_MESODIFF)
      call make_k3d(k2d,k3d)
#endif

#if defined(OCN_GISS_MESO)
      CALL OCN_mesosc(k3d)
! This zero-out snippet was added so that previous results under
! OCN_GISS_MESO option are preserved identically after separating
! the various calculations of mesoscale diffusivity.  It corresponds
! to the following logical in old routine DENSGRAD (now ORIG_MESODIFF):
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        if(zeroing_mask(i,j)) k3d(i,j,:) = 0.
        k3d(i,j,:) = k2d(i,j)
      enddo
      enddo
      enddo
#endif

      if(have_south_pole) then
        do l=1,lmo
          k3d(2:im,1,l) = k3d(1,1,l)
        enddo
      endif
      if(have_north_pole) then
        do l=1,lmo
          k3d(2:im,jm,l) = k3d(1,jm,l)
        enddo
      endif


      call gmkdif(k3d,1d0)
      call gmfexp(g0m,gxmo,gymo,gzmo,.false.,oijl(1,j_0h,1,ijl_ggmfl))
      call gmfexp(s0m,sxmo,symo,szmo,.true. ,oijl(1,j_0h,1,ijl_sgmfl))
#ifdef TRACERS_OCEAN
      do n = 1,tracerlist%getsize()
        trentry=>tracerlist%at(n)
        call gmfexp(trmo(1,j_0h,1,n),
     &         txmo(1,j_0h,1,n),tymo(1,j_0h,1,n),tzmo(1,j_0h,1,n),
     &         trentry%t_qlimit,toijl(1,j_0h,1,toijl_gmfl,n))
      end do
#endif

#ifdef OCN_GISS_MESO
c     CALL MESO_D(G0M0,GXMO0,GYMO0,G0M,GXMO,GYMO,GZMO)
      CALL MESO_D(G0M0,GXMO0,GYMO0,GXMO,GYMO,GZMO)
      G0M=G0M0
c     CALL MESO_D_TEST(G0M0,G0M,GXMO0,GYMO0,GZMO)
c     CALL MESO_D(S0M0,SXMO0,SYMO0,S0M,SXMO,SYMO,SZMO)

      CALL MESO_A(G0M,GXMO,GYMO,GZMO)
c     CALL MESO_A(S0M,SXMO,SYMO,SZMO)
#endif

      ! set diagnostics
      do j=j_0,j_1
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        oij(i,j,ij_gmsc) = oij(i,j,ij_gmsc) + k3d(i,j,1)
c        oij(i,j,ij_gmscz) = oij(i,j,ij_gmscz) + gmscz(i,j)
      enddo
      enddo

      end subroutine ocnmeso_drv

      subroutine densgrad
!@sum  DENSGRAD calculates all horizontal and vertical density gradients
      use ocean_dyn, only : bydh
      use ocnmeso_com, only : g3d,s3d,p3d,rhox,rhoy,rhomz,byrhoz,
     &     rho=>r3d,vbar=>v3d,bydzv=>bydze3d,dzv=>dze3d
      use constant, only: grav
      !use ocean, only : bydyv,bydxp
      use ocean, only : dyvo,dxpo
      use ocean, only : G0M,GZM=>GZMO, S0M,SZM=>SZMO, OPRESS, FOCEAN, MO
      use ocean, only : im,jm,lmo,lmm,lmu,lmv,dxypo
      use ocean_dyn, only  : dh
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds,halo_update_column
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      implicit none

c
      REAL*8  BYRHO,DZVLM1,ARHO,ARHOX,ARHOY,ARHOZ,DH0,R1,R2,P12
      INTEGER IM1,IMAX
      Real*8,External   :: VOLGSP

      integer :: i,j,k,l,n,il,ir,jl,jr

      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE

      Real*8, Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Real*8, dimension(lmo) :: gup,gdn,sup,sdn,pm
!@var dVBARdZ specific volume vertical difference (ref to lower point pressure)
      Real*8 :: vup,vdn,vupu,vdnu,bym,dvbardz

      real*8, dimension(0:lmo) :: pe

c**** Extract domain decomposition info
      CALL GETDOMAINBOUNDS(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE)


      RHOMZ = -0.0; BYRHOZ = -0.0;
      RHO = -0.0; BYDH = -0.0;
      DZV = -0.0; BYDZV = -0.0;

      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        pe(0) = opress(i,j)
        do l=1,lmm(i,j)
          bym = 1d0/(mo(i,j,l)*dxypo(j))
          pe(l) = pe(l-1) + mo(i,j,l)*grav
          g3d(l,i,j) = g0m(i,j,l)*bym
          s3d(l,i,j) = s0m(i,j,l)*bym
          p3d(l,i,j) = .5*(pe(l)+pe(l-1))

          GUP(L)=(G0M(I,J,L)-2*z12eH*GZM(I,J,L))*BYM
          GDN(L)=(G0M(I,J,L)+2*z12eH*GZM(I,J,L))*BYM
          SUP(L)=(S0M(I,J,L)-2*z12eH*SZM(I,J,L))*BYM
          SDN(L)=(S0M(I,J,L)+2*z12eH*SZM(I,J,L))*BYM
          sup(l) = max(0d0,sup(l))
          sdn(l) = max(0d0,sdn(l))
C**** Calculate potential specific volume (ref to mid-point pr)
          PM(L) = P3D(L,I,J)
          VUP = VOLGSP (GUP(L),SUP(L),PM(L))
          VDN = VOLGSP (GDN(L),SDN(L),PM(L))
          VBAR(L,I,J) = (VUP + VDN)*.5
C**** Vertical gradient calculated using lower box mid-point pr
          IF (L.gt.1) then 
            VUPU = VOLGSP (GUP(L-1),SUP(L-1),PM(L))
            VDNU = VOLGSP (GDN(L-1),SDN(L-1),PM(L))
            dVBARdZ = .5* (VUP + VDN - VUPU - VDNU)
            DZVLM1     = 0.5* (DH(I,J,L) + DH(I,J,L-1))
            DZV(I,J,L-1) = DZVLM1
            BYDZV(I,J,L-1) = 1d0/DZVLM1
            RHOMZ(I,J,L-1)=MAX(0d0,
     *             -dVBARdZ*BYDZV(I,J,L-1)/VBAR(L-1,I,J)**2)
C**** minus vertical gradient
            IF(RHOMZ(I,J,L-1).ne.0.)
     *             BYRHOZ(I,J,L-1)=1./RHOMZ(I,J,L-1)
          end if
C**** RHO(I,J,L)  Density=1/specific volume
          RHO(L,I,J)  = 1d0/VBAR(L,I,J)
          BYDH(I,J,L) = 1d0/DH(I,J,L)

        enddo
      enddo
      enddo
      enddo

C**** Copy to all longitudes at poles
      If(have_north_pole) Then
        Do L=1,LMM(1,JM)
          RHO(L,2:IM,JM) = RHO(L,1,JM)
          VBAR(L,2:IM,JM) = VBAR(L,1,JM)
          G3D(L,2:IM,JM) = G3D(L,1,JM)
          S3D(L,2:IM,JM) = S3D(L,1,JM)
          P3D(L,2:IM,JM) = P3D(L,1,JM)
          DZV(2:IM,JM,L) = DZV(1,JM,L)
          BYDZV(2:IM,JM,L) = BYDZV(1,JM,L)
          BYDH(2:IM,JM,L) = BYDH(1,JM,L)
          RHOMZ(2:IM,JM,L) = RHOMZ(1,JM,L)
          BYRHOZ(2:IM,JM,L) = BYRHOZ(1,JM,L)
        EndDo
      EndIf
      If(have_south_pole) Then
        Do L=1,LMM(1,1)
          RHO(L,2:IM,1) = RHO(L,1,1)
          VBAR(L,2:IM,1) = VBAR(L,1,1)
          G3D(L,2:IM,1) = G3D(L,1,1)
          S3D(L,2:IM,1) = S3D(L,1,1)
          P3D(L,2:IM,1) = P3D(L,1,1)
          DZV(2:IM,1,L) = DZV(1,1,L)
          BYDZV(2:IM,1,L) = BYDZV(1,1,L)
          BYDH(2:IM,1,L) = BYDH(1,1,L)
          RHOMZ(2:IM,1,L) = RHOMZ(1,1,L)
          BYRHOZ(2:IM,1,L) = BYRHOZ(1,1,L)
        EndDo
      EndIf

      call halo_update_column(grid,rho)
      call halo_update_column(grid,g3d)
      call halo_update_column(grid,s3d)
      call halo_update_column(grid,p3d)

C**** Calculate density gradients
      rhox = 0.
      rhoy = 0.
      do j=j_0s,j_1s
      do n=1,nbyzu(j,1)
      do i=i1yzu(n,j,1),i2yzu(n,j,1)
        il = i
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif
        do l=1,lmu(i,j)
          p12 = .5d0*(p3d(l,il,j)+p3d(l,ir,j)) 
          rhox(l,i,j) = (
     &         1d0/volgsp(g3d(l,ir,j),s3d(l,ir,j),p12)
     &        -1d0/volgsp(g3d(l,il,j),s3d(l,il,j),p12)
c     &         )*bydxp(j)!/dxpo(j)
     &         )*(1d0/dxpo(j))
        enddo
      enddo
      enddo
      enddo
      do j=max(2,j_0-1),j_1s
      do n=1,nbyzv(j,1)
      do i=i1yzv(n,j,1),i2yzv(n,j,1)
        jl = j
        jr = j+1
        do l=1,lmv(i,j)
          p12 = .5d0*(p3d(l,i,jl)+p3d(l,i,jr)) 
          rhoy(l,i,j) = (
     &         1d0/volgsp(g3d(l,i,jr),s3d(l,i,jr),p12)
     &        -1d0/volgsp(g3d(l,i,jl),s3d(l,i,jl),p12)
c     &         )*bydyv(j)!/dyvo(j)
     &         )*(1d0/dyvo(j))
        enddo
      enddo
      enddo
      enddo

      end subroutine densgrad

      subroutine simple_mesodiff(ainv)
      use domain_decomp_1d, only : getdomainbounds
      use oceanr_dim, only : grid=>ogrid
      use constant, only : grav
      use ocnmeso_com, only : rho=>r3d,rhomz,rhox,rhoy
      use ocnmeso_com, only : kvismult,kbg
      use ocean, only : im,jm,lmo,lmm,lmu,lmv,sinpo,ze,dzo
      use ocean, only : sinic,cosic
      use ocean, only : nbyzv,i1yzv,i2yzv
      implicit none
! A prescription of mesoscale K a la "Visbeck" that is being used for
! a GISS-MIT MIP (excepting the latitude dependence).
! It is like the default method in that
! K = factor*L*L*[upper-ocean average of 1/t_eady=N*isopycnal_slope]
! but with the following differences:
! 1. The length scale L is a constant rather than the deformation
!    radius R_d=NH/f (or beta-plane R_d=sqrt(NH/beta) at low latitudes)
! 2. Off the equator, the explicit latitude dependence is 1/sin(lat)
!    rather than the 1/sin**2(lat) arising from use of L=R_d.  Near the
!    equator, the 1/sin(lat) dependence is capped, while the use of
!    beta-plane R_d in the default method makes K independent of N.
! 3. As N goes to zero, K goes to zero via imposition of a maximum
!    isopycnal slope rather than via R_d.  As N goes to infinity
!    but 1/t_eady goes to zero, K goes to zero, while the N**2 in
!    R_d in the default method permits large K (off the equator)
! 4. N*isopycnal_slope is computed locally and averaged over the upper
!    ocean rather than computed using averages of density gradients.
!    Where the ocean is shallower than an imposed limit, the
!    upper-ocean average of 1/t_eady is computed as per the
!    comments above for the mindepth parameter.
! 5. Minimum and maximum K values are imposed.

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: ainv

! 
!@var snavg vertical average of the product of isoneutral slope times
!@+         N (sqrt of B-V freq).  snavg = 1/t_eady
      real*8 :: snavg,snsum,delz,ksurf,arhoh,arho,arhox,arhoy,arhoz
      real*8, parameter ::
!@param maxslope maximum slope used in t_eady calculation
     &         maxslope=4d-4
!@param mindepth (m) minimum averaging depth in t_eady calculation.  If
!@+              the local depth is less than mindepth, the calculation
!@+              is performed as if the depth were mindepth and isoneutral
!@+              slopes were zero between the local depth and mindepth.
     &        ,mindepth=400d0
!@param kfac the value of AMU used for SIMPLE_MESODIFF
     &        ,kfac0=12d-3  !*2d0
!@param lscale (m) fixed length scale for calculating mesoscale diffusivity
     &        ,lscale=2.5d5
!@param maxk (m2/s) upper bound for mesoscale diffusivity
     &        ,maxk=6000d0
!@param minsinlat minimum for 1/sin(lat) scaling of mesoscale diffusivity
     &        ,minsinlat=.1d0 ! roughly sin(6 degrees)
!
!@var mink (m2/s) lower bound for mesoscale diffusivity
      real*8 :: mink=100d0

      integer :: i,j,l,n,im1,lav
      integer :: j_0, j_1, j_0s, j_1s
      logical :: have_north_pole
      integer :: l1k
      real*8 :: kfac

      mink = kbg

      kfac = kfac0*kvismult

C**** Calculate level at 1km depth
      do l1k=1,lmo-1
        if(ze(l1k+1) .ge. 1d3) exit
      enddo

      call getdomainbounds(grid, j_strt = j_0, j_stop = j_1,
     &               j_strt_skp  = j_0s,   j_stop_skp  = j_1s,
     &               have_north_pole = have_north_pole)

      do j=j_0s,j_1s
      do i=1,im
        if(lmm(i,j).eq.0) cycle
        im1 = i-1
        if(i.eq.1) im1=im
        arho = 0.
        lav = min(l1k,lmm(i,j))
        snsum = 0.
        do l=1,lav
          arho  = arho  + rho(l,i,j)
          arhoz = .5d0*(rhomz(i,j,max(1,l-1))
     &                + rhomz(i,j,min(l,lmm(i,j)-1)) )
          if(arhoz.le.0.) cycle
          arhox = .5d0*(rhox(l,im1,j  )+rhox(l,i,j))
          arhoy = .5d0*(rhoy(l,i  ,j-1)+rhoy(l,i,j))
          arhoh = sqrt(arhox*arhox+arhoy*arhoy)
          snsum = snsum + dzo(l)*sqrt(arhoz)*min(maxslope, arhoh/arhoz)
        enddo
        arho  = arho / real(lav,kind=8)
        delz = max(mindepth, .5*(ze(lav)+ze(lav-1)-ze(1)))
        snavg = snsum * sqrt(grav/arho) / delz
        ksurf = (kfac * (lscale**2)) *snavg
        ksurf = ksurf/max(abs(sinpo(j)),minsinlat)
        ksurf = max(mink, min(maxk, ksurf))
        ainv(i,j) = ksurf
      enddo
      enddo

      if(have_north_pole) then
        arho = 0.
        lav = min(l1k,lmm(1,jm))
        snsum = 0.
        do l=1,lav
          i = 1
          j = jm
          arho  = arho  + rho(l,i,j)
          arhoz = .5d0*(rhomz(i,j,max(1,l-1))
     &                + rhomz(i,j,min(l,lmm(i,j)-1)) )
          if(arhoz.le.0.) cycle
          arhox = 0.
          arhoy = 0.
          j = jm-1
          do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            arhox = arhox - sinic(i)*rhoy(l,i,j)
            arhoy = arhoy + cosic(i)*rhoy(l,i,j)
          enddo
          enddo
          arhox = arhox*2/im
          arhoy = arhoy*2/im
          arhoh = sqrt(arhox*arhox+arhoy*arhoy)
          snsum = snsum + dzo(l)*sqrt(arhoz)*min(maxslope, arhoh/arhoz)
        enddo
        i = 1
        j = jm
        arho  = arho / real(lav,kind=8)
        delz = max(mindepth, .5*(ze(lav)+ze(lav-1)-ze(1)))
        snavg = snsum * sqrt(grav/arho) / delz
        ksurf = (kfac * (lscale**2)) *snavg
        ksurf = ksurf/max(abs(sinpo(j)),minsinlat)
        ksurf = max(mink, min(maxk, ksurf))
        ainv(i,j) = ksurf
        ainv(2:im,j) = ainv(1,j)
      endif

      end subroutine simple_mesodiff

      subroutine make_k3d(k2d,k3d)
! to be replaced with nontrivial version in next commit
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm,im,lmo
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     k3d
c
      integer :: i,j,l,n

      integer :: j_0,j_1

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1)

      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        k3d(i,j,:) = k2d(i,j)
      enddo
      enddo
      enddo

      end subroutine make_k3d

      subroutine orig_mesodiff(k2d,zeroing_mask)
!@auth Gavin Schmidt/Dan Collins
      use constant, only : omega,radius,grav
      use ocean, only : im,jm,lmo,dzo,ze,sinpo,cospo,dypo,lmm,lmu,lmv
      use ocnmeso_com, only : rhox,rhoy
      use ocnmeso_com, only : g3d,s3d,p3d,vbar=>v3d,rho=>r3d
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k2d
      logical, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     zeroing_mask ! temporarily passed out to preserve results
                        ! for OCN_GISS_MESO option
c
      REAL*8  BYRHO,CORI,BETA,ARHO,ARHOX,ARHOY,ARHOZ,AN,RD
     *     ,BYTEADY,DZSUMX,DZSUMY,R1,R2,P12
      REAL*8 :: HUP
      INTEGER I,J,L,IM1,LAV
      Real*8,External   :: VOLGSP

      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE 
      real*8 :: rhoz1k,vdn,vup
      integer :: lavm
!@var LUP level corresponding to 1km depth
      INTEGER :: LUP

!@var AMU = Visbeck scheme scaling parameter (1)
      REAL*8, PARAMETER :: AMU = 0.13d0

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)



C**** Calculate level at 1km depth
        LUP=0
   10   LUP=LUP + 1
        IF (ZE(LUP+1).lt.1d3) GOTO 10
        HUP = ZE(LUP)

C**** Calculate VMHS diffusion = amu* min(NH/f,equ.rad)^2 /Teady
      K2D = 0.
      zeroing_mask = .false.
      DO J=J_0S,J_1S
        CORI = ABS(2d0*OMEGA*SINPO(J))
        BETA = ABS(2d0*OMEGA*COSPO(J)/RADIUS)
        IM1=IM
        DO I=1,IM
          IF (LMM(I,J).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
            ARHO  = 0.
            ARHOX = 0.
            ARHOY = 0.
            DZSUMX = 0.
            DZSUMY = 0.
            LAV = MIN(LUP,LMM(I,J))
            DO L=1,LAV
              ARHO  = ARHO  + RHO(L,I,J)
              IF(LMU(IM1,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(L,IM1,J)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMU(I  ,J).ge.L) THEN
                ARHOX = ARHOX + RHOX(L,I  ,J)*DZO(L)
                DZSUMX = DZSUMX + DZO(L)
              END IF
              IF(LMV(I,J-1).ge.L) THEN
                ARHOY = ARHOY + RHOY(L,I,J-1)*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
              IF(LMV(I,J  ).ge.L) THEN
                ARHOY = ARHOY + RHOY(L,I,J  )*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
            ARHO  = ARHO / REAL(LAV,KIND=8)
            IF (DZSUMX.gt.0.) ARHOX = ARHOX / DZSUMX
            IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
            IF (LAV.gt.1) THEN
C**** Vertical potential density gradient in top 1km
              LAVM = MAX(LAV/2,1) ! mid depth
              VUP = VOLGSP (G3D(1,I,J),S3D(1,I,J),P3D(LAVM,I,J))
              VDN = VOLGSP (G3D(LAV,I,J),S3D(LAV,I,J),P3D(LAVM,I,J))
              RHOZ1K = (VUP - VDN)/VBAR(LAVM,I,J)**2
              ARHOZ=2*RHOZ1K/(ZE(LAV)+ZE(LAV-1)-ZE(1))
            ELSE
              ARHOZ = 0.
            END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
            IF (ARHOZ.gt.0) THEN
              AN = SQRT(GRAV * ARHOZ / ARHO)
#ifdef CONSTANT_MESO_LENSCALE
              RD = meso_lenscale_const
#else
              RD = AN * HUP / CORI
              IF (RD.gt.ABS(J-.5*(JM+1))*DYPO(J)) RD=SQRT(AN*HUP/BETA)
#endif
              BYTEADY = GRAV * SQRT(ARHOX*ARHOX + ARHOY*ARHOY) / (AN
     *             *ARHO)
              K2D(I,J) = AMU * RD**2 * BYTEADY ! was = AIN
            ELSE
              zeroing_mask(i,j) = .true.
            END IF
          END IF
          IM1=I
        END DO
      END DO
C**** North pole
      if ( HAVE_NORTH_POLE ) then
        IF (LMM(1,JM).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,JM))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(L,1,JM)
            DO I=1,IM
              IF(LMV(I,JM-1).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(L,I,JM-1))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
C**** Vertical potential density gradient in top 1km
            LAVM = MAX(LAV/2,1) ! mid depth
            VUP = VOLGSP (G3D(1,1,JM),S3D(1,1,JM),P3D(LAVM,1,JM))
            VDN = VOLGSP (G3D(LAV,1,JM),S3D(LAV,1,JM),P3D(LAVM,1,JM))
            RHOZ1K = (VUP - VDN)/VBAR(LAVM,1,JM)**2
            ARHOZ=2*RHOZ1K/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            K2D(1,JM) = AMU * RD**2 * BYTEADY ! was = AIN
          ELSE
            zeroing_mask(1,jm) = .true.
          END IF
        END IF
        K2D(2:IM,JM)=K2D(1,JM)
      endif
C**** South pole
      if ( HAVE_SOUTH_POLE ) then
        IF (LMM(1,1).gt.0) THEN
C**** Calculate average density + gradients over [1,LUP]
          ARHO  = 0. ; ARHOY = 0. ;  DZSUMY = 0.
          LAV = MIN(LUP,LMM(1,1))
          DO L=1,LAV
            ARHO  = ARHO  + RHO(L,1,1)
            DO I=1,IM
              IF(LMV(I,2).ge.L) THEN
! take abs to get a non-directional scale
                ARHOY = ARHOY + ABS(RHOY(L,I,2))*DZO(L)
                DZSUMY = DZSUMY + DZO(L)
              END IF
            END DO
          END DO
          ARHO  = ARHO / REAL(LAV,KIND=8)
          IF (DZSUMY.gt.0.) ARHOY = ARHOY / DZSUMY
          IF (LAV.gt.1) THEN
C**** Vertical potential density gradient in top 1km
            LAVM = MAX(LAV/2,1) ! mid depth
            VUP = VOLGSP (G3D(1,1,1),S3D(1,1,1),P3D(LAVM,1,1))
            VDN = VOLGSP (G3D(LAV,1,1),S3D(LAV,1,1),P3D(LAVM,1,1))
            RHOZ1K = (VUP - VDN)/VBAR(LAVM,1,1)**2
            ARHOZ=2*RHOZ1K/(ZE(LAV)+ZE(LAV-1)-ZE(1))
          ELSE
            ARHOZ = 0.
          END IF
C**** avoid occasional inversions. IF ARHOZ<=0 then GM is pure vertical
C**** so keep at zero, and let KPP do the work.
          IF (ARHOZ.gt.0) THEN
            AN = SQRT(GRAV * ARHOZ / ARHO)
            CORI = ABS(2d0*OMEGA*SINPO(JM))
#ifdef CONSTANT_MESO_LENSCALE
            RD = meso_lenscale_const
#else
            RD = AN * HUP / CORI
#endif
            BYTEADY = GRAV * ARHOY / (AN*ARHO)
            K2D(1,1) = AMU * RD**2 * BYTEADY ! was = AIN
          ELSE
            zeroing_mask(1,1) = .true.
          END IF
        END IF
        K2D(2:IM,1)=K2D(1,1)
      endif

      return
      end subroutine orig_mesodiff

#ifdef USE_1D_MESODIFF
      subroutine get_1d_mesodiff(k3d)
      use ocean, only : nbyzm,i1yzm,i2yzm,lmo,im
      use oceanr_dim, only : grid=>ogrid
      use ocnmeso_com
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &   k3d
      integer :: i,j,l,n
      real*8 :: k1d(lmo)
#if (defined CONSTANT_MESO_DIFFUSIVITY)
      k1d(:) = meso_diffusivity_const
#elif (defined EXPDECAY_MESO_DIFFUSIVITY)
      do l=1,lmo
        k1d(l) = meso_diffusivity_z0*
     &       exp(-(.5*(ze(l-1)+ze(l)))/meso_diffusivity_zscale)
     &       + meso_diffusivity_deep
      enddo
#endif

      do l=1,lmo
      do j=grid%j_strt,grid%j_stop
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        k3d(i,j,l) = k1d(l)
      enddo
      enddo
      enddo
      enddo

      end subroutine get_1d_mesodiff
#endif
