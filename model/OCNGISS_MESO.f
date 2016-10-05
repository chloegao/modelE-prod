#include "rundeck_opts.h"

#ifdef OCN_GISS_MESO

      module gissmeso_com
      implicit none

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
     &  omrx3d,omry3d,etx3d,ety3d,etz3d
      real*8, allocatable, dimension(:,:,:) ::
     &  rhoxc,rhoyc,uoc,voc
      integer, allocatable, dimension(:,:) :: lhstar2d

      end module gissmeso_com

      subroutine alloc_gissmeso_com
      use ocean, only : im,lmo
      use gissmeso_com
      use dictionary_mod, only : sync_param
      use oceanr_dim, only : ogrid
      use domain_decomp_1d, only : getdomainbounds
      use dictionary_mod
      implicit none
      integer :: j_0h,j_1h

      call getdomainbounds(ogrid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      allocate(   omrx3d (im,j_0h:j_1h,lmo) )
      allocate(   omry3d (im,j_0h:j_1h,lmo) )
      allocate(    etx3d (im,j_0h:j_1h,lmo) )
      allocate(    ety3d (im,j_0h:j_1h,lmo) )
      allocate(    etz3d (im,j_0h:j_1h,lmo) )
      allocate( lhstar2d (im,j_0h:j_1h) )

      allocate( rhoxc(lmo,im,j_0h:j_1h) )
      allocate( rhoyc(lmo,im,j_0h:j_1h) )
      allocate( uoc(im,j_0h:j_1h,lmo) )
      allocate( voc(im,j_0h:j_1h,lmo) )

      lhstar2d = 2;
      omrx3d = 0.; omry3d = 0.; etx3d = 0.; ety3d = 0.; etz3d=0.
      rhoxc=0.; rhoyc=0.; uoc=0.; voc=0.

      end subroutine alloc_gissmeso_com

      subroutine giss_meso(k3d,zeroing_mask)
!@sum  giss_meso calculates the skew and diffusion coefficients
!@+    using the giss meso model, updating the variables
!@+    k3d,xix3d,xiy3d,etx3d,ety3d and lhstar2d
!@auth ALeboissetier/YCheng/AHoward
!@ver  2016/07/04
!@var  epe eddy potential energy (m2/s2)

      use ocnmeso_com, only : rhox,rhoy,rhomz,byrhoz,
     &     rho=>r3d,vbar=>v3d,bydzv=>bydze3d,dzv=>dze3d,kbg
      use gissmeso_com
      USE OCEAN, only : im,jm,lmo,lmm,lmu,lmv,dts,cospo,sinpo,ze,dxypo
     *     ,mo,dypo,dyvo,dxpo,dzo,uo,vo,uod,vod,focean,sinvo,cosvo
      use ocean_dyn, only  : dh
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds,halo_update_column,
     &                             halo_update,north,south

      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use constant, only : grav,omega,omega2,teeny,pi,radius
      use kpp_com, only : kpl

      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc
     *     ,ijl_k3d,ijl_eke,ijl_sxc,ijl_syc,ijl_n2c,ijl_dopplerspd
     *     ,ijl_k31,ijl_k32,ijl_k33,ijl_fvbm,ijl_taper,ijl_epe
     *     ,ij_rd,ij_ekes,ij_udrift,ij_vdrift,ij_driftspd,ij_hstar
     *     ,ij_var1,ij_var2,ij_ekes_d,ij_sstar,ij_divri

      implicit none
 
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LMO),
     &  intent(out) :: K3D
      logical, dimension(im,grid%j_strt_halo:grid%j_stop_halo),
     &  intent(in) :: zeroing_mask ! temporarily passed out to preserve results

      integer :: i,j,k,l,n

      integer :: j_0, j_1, j_0h,j_1h, j_0s,j_1s
      logical :: have_south_pole,have_north_pole

      real*8, parameter :: n2min=1d-8

      real*8, dimension(lmo) :: zg,zi,n2,byrho
     &    ,uc,vc,sxc,syc
      real*8, dimension(lmo) :: sx,sy,dopplerspd,eke,kappam
     &    ,omrx,omry,etx,ety,etz,k31,k32,k33,fvbm,taper,epe

      real*8 f,dfdy,af,sf,rd
     &    ,byrz_bot,byrz_top,ekes,udr,vdr,veldr,byrhomz,rhomzl
     &    ,var1,var2,ekes_d,sstar,divri
      integer im1,lmij,lhstar,lm1

      ! horizonal c-grid:
      !  
      !               (i,j): v,ud,rhoy,sy,xiy,ety
      !           ------x-----
      !          -            -
      !         -     T(i,j)   -
      ! U(im1,j) x      *       x (i,j): u,vd,rhox,sx,xix,etx
      !       -                   -
      !      -                     -
      !     ------------x------------
      !              V(i,j-1)
   


      ! Vertical grid diagram
      !
      !       grid levels, zg<0                  interface levels, zi<0
      !
      !                   -----------------------------  surf
      !                1    - - - - - - - - - - - - -
      !                   -----------------------------  1
      !                2    - - - - - - - - - - - - -
      !                   -----------------------------  2
      !               l-1   - - - - - - - - - - - - -
      !                   -----------------------------  l-1
      !zg(l)<0,tracer  l    - - - - - - - - - - - - -
      !                   -----------------------------  l, zi(l),byrhoz(l)
      !               l+1   - - - - - - - - - - - - -
      !                   -----------------------------  l+1
      !                     - - - - - - - - - - - - -
      !                   -----------------------------  lmi1-1 (=lme) 
      !             lmij    - - - - - - - - - - - - -
      !                   -----------------------------  lmij

      !             dh(l)=zi(l-1)-zi(l) is the layer thickness


      k3d=0.d0; eke=0.; dopplerspd=0.
      omrx3d = 0.; omry3d = 0.; etx3d = 0.; ety3d = 0.; etz3d=0.
      lhstar2d = kpl;
      rhoxc=0.; rhoyc=0.; uoc=0.; voc=0.

      ! it seems that I must define them at l=lmij ?
      ! zg < 0 at tracer level
      ! zi < 0 at interface level

c**** Extract domain decomposition info
      call getdomainbounds(grid, j_strt = j_0, j_stop = j_1,
     &     j_strt_skp = j_0s, j_stop_skp = j_1s,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h,
     &     have_south_pole = have_south_pole,
     &     have_north_pole = have_north_pole)

      call get_qxc_qyc(mo,rhox,rhoy,0,rhoxc,rhoyc)
      call get_qxc_qyc_ijl(mo,uo,vo,0,uoc,voc)

      do j=j_0,j_1

      f=omega2*sinpo(j)      !evaluated at tracer points
      dfdy=omega2*cospo(j)/radius
      af=abs(f)
      sf=f/af

      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        im1=i-1
        if(i.eq.1) im1=im
        lmij=lmm(i,j)

        ! zg (<0, at layer middle), zi (<0, at interface level)
        zg(1)=-.5d0*dh(i,j,1)
        zi(1)=-dh(i,j,1)
        do l=2,lmij
          zg(l)=zi(l-1)-.5d0*dh(i,j,l)
          zi(l)=zi(l-1)-dh(i,j,l)
        enddo
        !blht=.5d0*dh(i,j,lmij) ! bottom layer half thickness
        !blht=0.d0
        

        if(kpl(i,j).ge.lmij) then
          ! no A-regime
          lhstar2d(i,j)=lmij
          cycle
        endif

        do l=1,lmij
          byrho(l)=1.d0/rho(l,i,j)
          uc(l)=uoc(i,j,l)   
          vc(l)=voc(i,j,l)
        end do

        do l=2,lmij
          ! cell center and layer moddle:
          !   uc,vc,sxc,syc,rhoxc,rhoyc,byrho,n2
          ! cell east edge and layer middle:
          !   ue,ve,sxe,sye,rhoxe,rhoye
          ! cell north edge and layer middle:
          !   vn,un,sxn,syn,rhoyn,rhoxn
              rhomzl=.5d0*(rhomz(i,j,l)+rhomz(i,j,l-1))
              n2(l)=grav*byrho(l)*rhomzl
              if(n2(l).le.n2min) then
                n2(l)=n2min
                rhomzl=n2(l)/(grav*byrho(l)+teeny)
c               sxc(l)=sign(smin,sxc(l))
c               syc(l)=sign(smin,syc(l))
              endif
              sxc(l)=rhoxc(l,i,j)/(rhomzl+teeny)
              syc(l)=rhoyc(l,i,j)/(rhomzl+teeny)
        end do
        sxc(1)=sxc(2)
        syc(1)=syc(2)
        n2(1)=n2(2)

        call giss_meso_1d(
        ! in:
     &     lmij,kpl(i,j),i,j,f,dfdy,dxypo(j),zg
     &    ,uc,vc,sxc,syc,n2,mo(i,j,1:lmij)
        ! out:
     &    ,lhstar,rd,ekes,ekes_d,sstar,udr,vdr,veldr,var1,var2,divri
     &    ,sx,sy,dopplerspd,kappam,eke,omrx,omry,etx,ety,etz,fvbm
     &    ,taper,epe)

        k31=0.; k32=0.; k33=0.;
        do l=1,lmij
          if(l.le.lhstar) then
            k31(l)=etx(l)
            k32(l)=ety(l)
            k33(l)=etz(l)
          else
            k31(l)=sx(l)*(2.d0-omrx(l))
            k32(l)=sy(l)*(2.d0-omry(l))
            k33(l)=sx(l)**2+sy(l)**2
          endif
          k31(l)=k31(l)/(2.d0*sx(lhstar))
          k32(l)=k32(l)/(2.d0*sy(lhstar))
          k33(l)=k33(l)/(sx(lhstar)**2+sy(lhstar)**2)
        end do

         if(lhstar.ge.lmij) then
           ! no A-regime
           cycle
         endif
        lhstar2d(i,j)=lhstar

        oij(i,j,ij_rd)= oij(i,j,ij_rd) + rd 
        oij(i,j,ij_ekes)= oij(i,j,ij_ekes) + ekes
        oij(i,j,ij_ekes_d)= oij(i,j,ij_ekes_d) + ekes_d
        oij(i,j,ij_sstar)= oij(i,j,ij_sstar) + sstar
        oij(i,j,ij_udrift)=oij(i,j,ij_udrift) + udr
        oij(i,j,ij_vdrift)=oij(i,j,ij_vdrift) + vdr
        oij(i,j,ij_driftspd)=oij(i,j,ij_driftspd) + veldr
        oij(i,j,ij_hstar)= oij(i,j,ij_hstar) - zg(lhstar) 
        oij(i,j,ij_var1)= oij(i,j,ij_var1) + var1
        oij(i,j,ij_var2)= oij(i,j,ij_var2) + var2
        oij(i,j,ij_divri)= oij(i,j,ij_divri) + divri
        ! k3d diag. is done in meso_drv
        oijl(i,j,:,ijl_eke)= oijl(i,j,:,ijl_eke) + eke(:) 
        oijl(i,j,:,ijl_sxc)= oijl(i,j,:,ijl_sxc) + sx(:) 
        oijl(i,j,:,ijl_syc)= oijl(i,j,:,ijl_syc) + sy(:) 
        oijl(i,j,:,ijl_n2c)= oijl(i,j,:,ijl_n2c) + n2(:) 
        oijl(i,j,:,ijl_dopplerspd)=
     &    oijl(i,j,:,ijl_dopplerspd) + dopplerspd(:)
        oijl(i,j,:,ijl_k31)= oijl(i,j,:,ijl_k31) + k31(:) 
        oijl(i,j,:,ijl_k32)= oijl(i,j,:,ijl_k32) + k32(:) 
        oijl(i,j,:,ijl_k33)= oijl(i,j,:,ijl_k33) + k33(:) 
        oijl(i,j,:,ijl_fvbm)= oijl(i,j,:,ijl_fvbm) + fvbm(:) 
        oijl(i,j,:,ijl_taper)= oijl(i,j,:,ijl_taper) + taper(:) 
        oijl(i,j,:,ijl_epe)= oijl(i,j,:,ijl_epe) + epe(:) 

        do l=1,lmij
          k3d(i,j,l)=kappam(l)
          omrx3d(i,j,l)=omrx(l)
          omry3d(i,j,l)=omry(l)
          etx3d(i,j,l)=etx(l)
          ety3d(i,j,l)=ety(l)
          etz3d(i,j,l)=etz(l)
        end do

      end do ! end of main loops over i
      end do ! end of main loops over n
      end do ! end of main loops over j

      if(have_south_pole) then
        do l=1,lmo
             k3d(2:im,1,l) =         k3d(1,1,l)
          omrx3d(2:im,1,l) =      omrx3d(1,1,l)
          omry3d(2:im,1,l) =      omry3d(1,1,l)
           etx3d(2:im,1,l) =       etx3d(1,1,l)
           ety3d(2:im,1,l) =       ety3d(1,1,l)
           etz3d(2:im,1,l) =       etz3d(1,1,l)
            oijl(2:im,1,l,ijl_eke)= oijl(1,1,l,ijl_eke) 
            oijl(2:im,1,l,ijl_sxc)= oijl(1,1,l,ijl_sxc) 
            oijl(2:im,1,l,ijl_syc)= oijl(1,1,l,ijl_syc) 
            oijl(2:im,1,l,ijl_n2c)= oijl(1,1,l,ijl_n2c) 
            oijl(2:im,1,l,ijl_dopplerspd)= oijl(1,1,l,ijl_dopplerspd) 
            oijl(2:im,1,l,ijl_k31)= oijl(1,1,l,ijl_k31) 
            oijl(2:im,1,l,ijl_k32)= oijl(1,1,l,ijl_k32) 
            oijl(2:im,1,l,ijl_k33)= oijl(1,1,l,ijl_k33) 
            oijl(2:im,1,l,ijl_fvbm)= oijl(1,1,l,ijl_fvbm) 
            oijl(2:im,1,l,ijl_taper)= oijl(1,1,l,ijl_taper) 
            oijl(2:im,1,l,ijl_epe)= oijl(1,1,l,ijl_epe) 
        enddo
        lhstar2d(2:im,1) =    lhstar2d(1,1)
             oij(2:im,1,ij_rd)=    oij(1,1,ij_rd) 
             oij(2:im,1,ij_ekes)=  oij(1,1,ij_ekes) 
             oij(2:im,1,ij_ekes_d)=  oij(1,1,ij_ekes_d) 
             oij(2:im,1,ij_sstar)=  oij(1,1,ij_sstar) 
             oij(2:im,1,ij_udrift)=oij(1,1,ij_udrift) 
             oij(2:im,1,ij_vdrift)=oij(1,1,ij_vdrift) 
             oij(2:im,1,ij_driftspd)=oij(1,1,ij_driftspd) 
             oij(2:im,1,ij_hstar)= oij(1,1,ij_hstar) 
             oij(2:im,1,ij_var1)= oij(1,1,ij_var1) 
             oij(2:im,1,ij_var2)= oij(1,1,ij_var2) 
             oij(2:im,1,ij_divri)= oij(1,1,ij_divri) 
      endif
      if(have_north_pole) then
        do l=1,lmo
             k3d(2:im,jm,l) =         k3d(1,jm,l)
          omrx3d(2:im,jm,l) =      omrx3d(1,jm,l)
          omry3d(2:im,jm,l) =      omry3d(1,jm,l)
           etx3d(2:im,jm,l) =       etx3d(1,jm,l)
           ety3d(2:im,jm,l) =       ety3d(1,jm,l)
           etz3d(2:im,jm,l) =       etz3d(1,jm,l)
            oijl(2:im,jm,l,ijl_eke)= oijl(1,jm,l,ijl_eke) 
            oijl(2:im,jm,l,ijl_sxc)= oijl(1,jm,l,ijl_sxc) 
            oijl(2:im,jm,l,ijl_syc)= oijl(1,jm,l,ijl_syc) 
            oijl(2:im,jm,l,ijl_n2c)= oijl(1,jm,l,ijl_n2c) 
            oijl(2:im,jm,l,ijl_dopplerspd)= oijl(1,jm,l,ijl_dopplerspd) 
            oijl(2:im,jm,l,ijl_k31)= oijl(1,jm,l,ijl_k31) 
            oijl(2:im,jm,l,ijl_k32)= oijl(1,jm,l,ijl_k32) 
            oijl(2:im,jm,l,ijl_k33)= oijl(1,jm,l,ijl_k33) 
            oijl(2:im,jm,l,ijl_fvbm)= oijl(1,jm,l,ijl_fvbm) 
            oijl(2:im,jm,l,ijl_taper)= oijl(1,jm,l,ijl_taper) 
            oijl(2:im,jm,l,ijl_epe)= oijl(1,jm,l,ijl_epe) 
        enddo
        lhstar2d(2:im,jm) =    lhstar2d(1,jm)
             oij(2:im,jm,ij_rd) =   oij(1,jm,ij_rd) 
             oij(2:im,jm,ij_ekes)=  oij(1,jm,ij_ekes) 
             oij(2:im,jm,ij_ekes_d)=  oij(1,jm,ij_ekes_d) 
             oij(2:im,jm,ij_sstar)=  oij(1,jm,ij_sstar) 
             oij(2:im,jm,ij_udrift)=oij(1,jm,ij_udrift) 
             oij(2:im,jm,ij_vdrift)=oij(1,jm,ij_vdrift) 
             oij(2:im,jm,ij_driftspd)=oij(1,jm,ij_driftspd) 
             oij(2:im,jm,ij_hstar)= oij(1,jm,ij_hstar) 
             oij(2:im,jm,ij_var1)= oij(1,jm,ij_var1) 
             oij(2:im,jm,ij_var2)= oij(1,jm,ij_var2) 
             oij(2:im,jm,ij_divri)= oij(1,jm,ij_divri) 
      endif


      return
      end subroutine giss_meso

      subroutine giss_meso_1d(
        ! in:
     &     lm,kpl,i,j,fc0,dfcdy,dxypoj,z
     &    ,uc,vc,sx0,sy0,n20,mo
        ! out:
     &    ,lhstar,rd,ks,ksd,sstar,udr,vdr,veldr,var1,var2,divri
     &    ,sx,sy,dopplerspd,kappam,eke,omrx,omry,etx,ety,etz
     &    ,fvbm,taper,epe)
!@var epe eddy potental energy (m/s)^2

      use constant, only : pi,teeny
      implicit none
      integer, intent(in) :: lm,kpl,i,j
      real*8, intent(in) :: fc0,dfcdy,dxypoj
      ! z < 0
      real*8, dimension(lm), intent(in) :: z,uc,vc !z < 0
      real*8, dimension(lm), intent(in) :: n20,sx0,sy0,mo

      integer, intent(out) :: lhstar
      real*8, intent(out) :: rd,ks,ksd,sstar,udr,vdr,veldr,var1,var2
     &    ,divri
      real*8, dimension(lm), intent(out) :: sx,sy,dopplerspd,kappam,eke
     &    ,omrx,omry,etx,ety,etz,fvbm,taper,epe

      ! local:
      real*8, parameter :: 
     &   ck=12.d0,delb=40.d0
     &  ,n2max=1d-4,smax=1d-2,smin=1d-6,s2min=smin**2
     &  ,ksmin=1d-8 ! (m/s)^2
     &  ,ksmax=100.d0
     &  ,rmax=2.d0
     &  ,rmin=.5d0
     &  ,etmax=2.*smax
     &  ,kappammax=5d4
c    &  ,tau=86400.d0
     &  ,rdmin=1d2
     &  ,rdmax=5.d5
     &  ,sigma_t=1d0

      real*8, dimension(lm) :: b1,gamh,gamhi,gam,gami,dgamhdz
     &  ,dgamhdzsxi,dgamhdzsyi,gamhui,gamhvi,gamh3i
     &  ,s2,gg,ggi,integ1,integ2,ui,vi,ox,oy,upar,uper,s,phi
     &  ,ta1,ta2,db1dz,n2,p,q
      integer l,iter,flg
      real*8 qty,rd2,oxm,oym,bx,by,bya,bd,afc_eff
      real*8 a0,byhh_star,um,vm,dsxdzm,dsydzm,fac,hstar,nzmax,frd
      real*8 crx,cry,hh,tmp1,tmp2,coef,zet,srks,ome
      real*8 etmax1,etmax2,zero,byalf,factor
      real*8 afc,fc,utild,vtild,dbdx,dbdy,uper_star,rr
      real*8 rdm,k31,k32,k33,pr,ksa

      kappam=0.; dopplerspd=0.; eke=0.
      omrx=0.; omry=0.; etx=0.; ety=0.; etz=0.; fvbm=0.; taper=0.
      
      fc=fc0
      afc=abs(fc)

      ! find lhstar

      do iter=2,lm

        lhstar=iter
        do l=1,lm
          n2(l)=n20(l)
          sx(l)=sign(max(abs(sx0(l)),smin),sx0(l))
          sy(l)=sign(max(abs(sy0(l)),smin),sy0(l))
          s2(l)=sx(l)**2+sy(l)**2
          s(l)=sqrt(s2(l))
        end do

        if(lhstar.ge.lm) then

          ! no A-regime
          lhstar=lm
          b1=1.d0
          gamh=1.d0
          gam=1.d0
          um=0.d0
          vm=0.d0
          frd=rdmin*afc
          rd=rdmin
          rd2=rd**2
          ks=0.d0
          ksd=0.d0
          sstar=0.d0
          udr=0.d0
          vdr=0.d0
          veldr=0.d0
          var1=0.d0
          var2=0.d0
          divri=0.d0
          dopplerspd=0.d0
          kappam=0.d0;eke=0.d0
          omrx=0.d0;omry=0.d0;etx=0.d0;ety=0.d0;etz=0.d0;fvbm=0.d0;
          taper=0.
          return
        endif

          do l=lhstar+1,lm ! A-region except the 2 layers nearer to lhstar
            if(s(l).gt.smax) then
              factor=smax/s(l)
              s(l)=smax
              s2(l)=smax**2
              sx(l)=sx(l)*factor
              sy(l)=sy(l)*factor
            endif
            n2(l)=min(n2(l),n2max)
          end do

          call wkb(
          ! in:
     &      lm,z,n2,lhstar ! z < 0
          ! out:
     &     ,frd,b1,db1dz)
          a0=sqrt(abs(b1(lm)))
          do l=1,lhstar-1
            gamh(l)=1.
            gam(l)=1.
            dgamhdz(l)=0.
          end do
          do l=lhstar,lm
            gamh(l)=abs((a0+b1(l))/(1+a0))
            gam(l)=gamh(l)**2
            if((a0+b1(l)).le.0.) then
              dgamhdz(l)=-dB1dz(l)/(1+a0)
            else
              dgamhdz(l)=dB1dz(l)/(1+a0)
            endif
          end do
          call intgr_botup(z,gamh,gamhi,lhstar,lm)
          ! gamhi(l) is int of gamh fr z(lm) to z
          ! gamhi(lhstar) is int of gamh fr z(lm) to z(lhstar)
          byhh_star=1.d0/(gamhi(lhstar))

          do l=lhstar,lm
            gamhui(l)=gamh(l)*uc(l)
            gamhvi(l)=gamh(l)*vc(l)
          end do
          call intgr_botup(z,gamhui,gamhui,lhstar,lm)
          call intgr_botup(z,gamhvi,gamhvi,lhstar,lm)
          ! gamhui(l) is int of gamh*uc fr z(lm) to z
          ! gamhvi(l) is int of gamh*vc fr z(lm) to z
          ! um=<u>, vm=<v>
          um=gamhui(lhstar)*byhh_star
          vm=gamhvi(lhstar)*byhh_star

          frd=max(frd,rdmin*afc)
          rd=frd/afc
          !rr=((um**2+vm**2)/(dfcdy*dfcdy))**.25d0
          rd=max(rd,rdmin)
          rd=min(rd,rdmax)
          rd2=rd*rd

          do l=lhstar,lm
            dgamhdzsxi(l)=dgamhdz(l)*sx(l)
            dgamhdzsyi(l)=dgamhdz(l)*sy(l)
          end do
          call intgr_botup(z,dgamhdzsxi,dgamhdzsxi,lhstar,lm)
          call intgr_botup(z,dgamhdzsyi,dgamhdzsyi,lhstar,lm)
          ! dgamhdzsxi(l) is int of dgamhdz*sx fr z(lm) to z
          ! dgamhdzsyi(l) is int of dgamhdz*sy fr z(lm) to z
          ! <dsxdz>=dsxdzm, <dsydz>=dsydzm
          dsxdzm=(gamh(lhstar)*sx(lhstar)-gamh(lm)*sx(lm)
     &          -(dgamhdzsxi(lhstar)))*byhh_star
          dsydzm=(gamh(lhstar)*sy(lhstar)-gamh(lm)*sy(lm)
     &          -(dgamhdzsyi(lhstar)))*byhh_star

          ! drift velocity udr,vdr
          ! cr neglected in the manuscript
          ! crx=-len2*dfcdy    
          crx=-rd**2*dfcdy    
          cry=0.d0
          udr=crx+um+fc*rd2*(dsydzm-sy(lhstar)*byhh_star)
          vdr=    vm-fc*rd2*(dsxdzm-sx(lhstar)*byhh_star)
        ! D-regime:

        call dqtfg1(z,uc,ui,1,lhstar)
        call dqtfg1(z,vc,vi,1,lhstar)
        do l=1,lhstar
          ! upar: u parallel; uper: u perpendicular
          ! uarrow=ox/z, varrow=oy/z
          s(l)=sqrt(s2(l))+teeny
          ! make sure ox(1)=z(1)*(u(1)-udr), oy(1)=z(1)*(v(1)-vdr)
          ox(l)=z(l)*(uc(l)-2*udr)+ui(l)
          oy(l)=z(l)*(vc(l)-2*vdr)+vi(l)
          p(l)=(sy(l)*ox(l)-sx(l)*oy(l))/(s(l)*s(l)+teeny)/(fc*rd**2)
          q(l)=(sx(l)*ox(l)+sy(l)*oy(l))/(s(l)*s(l)+teeny)/(fc*rd**2)
        end do
        hstar=-z(lhstar)

        do l=1,lhstar
          phi(l)=(z(l)/hstar)**2*min(n2(l)/n2(lhstar),1.d0)
          !pr=min(p(l)/p(lhstar),1d0)
          pr=min(-p(l),1d0)
          etx(l)=pr*(1+phi(l))*sx(l)+q(l)*(1-phi(l))*sy(l)
          ety(l)=pr*(1+phi(l))*sy(l)-q(l)*(1-phi(l))*sx(l)
          if(abs(etx(l)).gt.etmax) etx(l)=sign(etmax,etx(l))
          if(abs(ety(l)).gt.etmax) ety(l)=sign(etmax,ety(l))
          if(etx(l)*sx(lhstar).lt.0.d0) etx(l)=sign(0.d0,sx(lhstar))
          if(ety(l)*sy(lhstar).lt.0.d0) ety(l)=sign(0.d0,sy(lhstar))
          if(abs(etx(l)).gt.abs(2*sx(lhstar)))
     &      etx(l)=sign(2*sx(lhstar),sx(lhstar))
          if(abs(ety(l)).gt.abs(2*sy(lhstar)))
     &      ety(l)=sign(2*sy(lhstar),sy(lhstar))
          etz(l)=pr*s2(l)*phi(l)
          etz(l)=max(etz(l),0.d0)
          etz(l)=min(etz(l),s2(lhstar))
          !qty=.5d0*(etx(l)*sx(l)+ety(l)*sy(l))/(s2(l)+teeny)
          taper(l)=(etx(l)*sx(l)+ety(l)*sy(l)-etz(l))/(s2(l)+teeny) ! from (3.5)
        end do
        if(taper(lhstar).ge.1.d0) then
          exit
        endif

      end do  ! end of iter loop
      
      !** surface kinetic enery **!
      
      hh=-z(lm)
      hstar=-z(lhstar)
      ! ksmax=(kappammax/len)**2
      ! note gam=1 fr lhstar to 0
      call intgr_botup(z,gam,gami,1,lm)
      do l=1,lm
      !do l=lhstar,lm ! yc
        gamh3i(l)=gamh(l)*gam(l)
      end do
      !call intgr_botup(z,gamh3i,gamh3i,lhstar,lm) ! yc
      call intgr_botup(z,gamh3i,gamh3i,1,lm)
      fac=sqrt(2.d0/pi)*hh/delb
      do l=lhstar,lm
        zet=(z(l)+hh)/delb
        gg(l)=fac*exp(-.5d0*zet*zet)*gam(l)
      end do
      call intgr_botup(z,gg,ggi,lhstar,lm)
      bd=ggi(lhstar)/gami(1)
      byalf=ck*rd/gamh3i(1)
      !byalf=ck*rd/gamh3i(lhstar) ! yc
      !coef=byalf/(1.d0+bd)
      coef=byalf*(1.d0-bd)

      ! ksa is the surface ks fr A-regime 
 
      do l=lhstar,lm
        integ1(l)=gamh(l)*s2(l)*n2(l)
      end do
      call intgr_botup(z,integ1,integ2,lhstar,lm)
      ! integ2(lhstar) is the int of integ1 fr -H to -hstar
      ksa=coef*rd*integ2(lhstar)
      ksa=min(max(ksa,ksmin),ksmax)

      ! ksd is the surface ks fr D-regime 

      ! bx==dbdx=-n2*sx, by==dbdy=-n2*sy, dbdz=n2
      do l=1,lhstar
        bx=-n2(l)*sx(l)
        by=-n2(l)*sy(l)
        ta1(l)=ox(l)*by-oy(l)*bx
      end do
      call dqtfg1(z,ta1,integ1,1,lhstar)
      !byalf=ck*rd/hstar ! yc
      !coef=byalf/(1.d0+bd) ! yc
      ksd=-coef/(rd*fc)*integ1(lhstar)
      ksd=min(max(ksd,ksmin),ksmax)

      ! ks is the total surface ke
 
      ks=ksd+ksa
      !ks=min(max(ks,ksmin),ksmax)
      srks=sqrt(ks)

      sstar=sqrt(sx(lhstar)**2+sy(lhstar)**2)

      !** kappam(z) **!

      do l=1,lm
        ome=1d0/( 1d0 + ((uc(l)-udr)**2+(vc(l)-vdr)**2)/(2*ks*gam(l)) )
        ome=max(ome,0.7d0)
        !ome=1d0
        eke(l)=ks*gam(l)
        kappam(l)=min(rd*srks*gamh(l)*ome,kappammax)
      end do
      ! modified OCNGM.f requires omrx,omry be 0 in D-regime
      do l=lhstar+1,lm
         omrx(l)=-1/(sx(l)*gamh(l)+teeny)*( -gamh(lm)*sx(lm)
     &        -dgamhdzsxi(l)-1/(fc*rd2)*(gamhvi(l)-gamhi(l)*vm)
     &        -gamhi(l)*(dsxdzm-sx(lhstar)*byhh_star) )
         omry(l)=-1/(sy(l)*gamh(l)+teeny)*( -gamh(lm)*sy(lm)
     &        -dgamhdzsyi(l)+1/(fc*rd2)*(gamhui(l)-gamhi(l)*um)
     &        -gamhi(l)*(dsydzm-sy(lhstar)*byhh_star) )
         ! rmin <= r <= rmax
         ! 1-rmax =< 1-r <= 1-rmin
         omrx(l)=min(max(omrx(l),1.d0-rmax),1.d0-rmin)
         omry(l)=min(max(omry(l),1.d0-rmax),1.d0-rmin)
      end do

      veldr=sqrt(udr**2+vdr**2)
      do l=1,lm
        dopplerspd(l)=sqrt((uc(l)-udr)**2+(vc(l)-vdr)**2)
        if(l.le.lhstar) then
          k31=etx(l)
          k32=ety(l)
          k33=etz(l)
        else
          k31=sx(l)*(2.d0-omrx(l))
          k32=sy(l)*(2.d0-omry(l))
          k33=sx(l)**2+sy(l)**2
        endif
        fvbm(l)=-kappam(l)*n2(l)*(-sx(l)*k31-sy(l)*k32+k33)
      end do
CAH160916 Sum kinematic power times gridbox mass to get integrated power.
CAH160916 Calculate integrated power here.
! Psxi stored in var1
      ! Ps=integ1(lhstar)
      ! Ps-Peps=integ2(lhstar)
      ! Peps=integ1(lhstar)-integ2(lhstar)
      var1=0.D0
      var2=0.D0
      do l=lhstar,lm
        var1 = var1 +
     &    kappam(l)*n2(l)*(sx(l)**2+sy(l)**2)*mo(l)*dxypoj
        var2 = var2 +
     &    kappam(l)*n2(l)*(sx(l)**2*omrx(l)+sy(l)**2*omry(l))
     &      *mo(l)*dxypoj
      end do
      ! epe = sigma_t (frd/N)^2 K_s (dB1/dz)^2
      do l=1,lm        
        epe(l) = sigma_t*(frd**2/n2(l))*ks*db1dz(l)**2
        integ1(l) = (eke(l)+epe(l))*(vc(l)-vdr)
      end do
      call intgr_botup(z,integ1,integ2,1,lm)
      divri=fc/(fc*rd)**2*integ2(1)

      return
      end subroutine giss_meso_1d

      subroutine wkb(lm,z,n2,lhstar,frd,B1,dB1dz)
!@sum wkb solves for the eigenvalues and eigenfunctions
!@+   using wkb method in a two-point boundary problem
!@auth A.Leboissetier/Y.Cheng
!@var B1 the first baroclinic mode solution of the eigenvalue problem
!@+   d/dz ((dB1/dz)/N2) + B1/(frd)**2 = 0
!@+   dB1/dz = 0 at z=0 and z(lmij) (z < 0)
!@var frd the eigenvalue

      use constant, only : pi,teeny

      IMPLICIT NONE
      
      integer, intent(in) :: lm    ! total layer number
      real*8, intent(in) :: z(lm) ! depth (<0), z0(1)=surface, z0(m)=bottom
      real*8, intent(in) :: n2(lm)! Brunt Vaisala frequency squared (1/s**2)
      integer, intent(in) :: lhstar ! layer edge number at hstar
      real*8, intent(out) :: frd   ! the eigenvalue, frd=f*rd, f is Coriol param
      real*8, intent(out) :: B1(lm),dB1dz(lm)! the eigenfunction
      
      real*8, dimension(lm) :: n,ni,phi,ta1,integ1
      real*8 :: hh,B
      integer l

      do l=1,lm
        n(l)=sqrt(n2(l))
      end do

      call intgr_botup(z,n,ni,lhstar,lm)

      frd=ni(lhstar)/pi
      frd=max(frd,teeny)
      hh=-z(lm)

      do l=lhstar,lm
        phi(l)=(n(l)/frd)**(-.5d0)*sin(ni(l)/frd)
        ta1(l)=(z(l)+hh)*n2(l)*phi(l)
      end do
      call intgr_botup(z,ta1,integ1,lhstar,lm)
      B=hh/integ1(lhstar)
      do l=lhstar,lm
        dB1dz(l)=B*n2(l)*phi(l)
      end do
      call dqtfg1(z,dB1dz,integ1,lhstar,lm)
      do l=lhstar,lm
        B1(l)=integ1(l)+1.d0
      end do
       ! frd and B1 are the final outputs
       do l=1,lhstar-1
         B1(l)=1.d0
         dB1dz(l)=0.d0
       end do

      RETURN
      end subroutine wkb

      subroutine intgr_botup(x,y,z,n1,n)
!@sum
!@+  integrate y bottom-up 
!@+  from x(n) to x(i) upto x(n1) and store it in z(i)
      implicit none
      integer, intent(in) :: n1,n
      real*8, intent(in) :: x(n),y(n)
      real*8, intent(out) :: z(n)
      integer i
      real*8 sum1,sum2
      if(n.lt.n1) then
        write(*,*) 'n.lt.n1 in intgr_botup, stop.'
        stop
      elseif(n.eq.n1) then
       z(n1)=0.d0
       return
      endif
      sum2=0.d0
      do i=n,n1+1,-1
        sum1=sum2
        sum2=sum2+.5d0*(x(i-1)-x(i))*(y(i-1)+y(i))
        z(i)=sum1
      end do
      z(n1)=sum2
      return
      end subroutine intgr_botup

      subroutine dqtfg1(x,y,z,n1,n)
      !@sum integrate y from x(n1) to x(i) upto x(n) and store it in z(i)
      implicit none
      integer, intent(in) :: n1,n
      real*8, intent(in) :: x(n),y(n)
      real*8, intent(out) :: z(n)
      integer i
      real*8 sum1,sum2
      sum2=0.
c     if(n-1)4,3,1
      if(n.lt.n1) then
        write(*,*) 'n.lt.n1 in dqtfg1, stop.'
        stop
      elseif(n.eq.n1) then
       z(n1)=0.d0
       return
      endif
      do i=n1+1,n
        sum1=sum2
        sum2=sum2+.5d0*(x(i)-x(i-1))*(y(i)+y(i-1))
        z(i-1)=sum1
      end do
      z(n)=sum2
      return
      end subroutine dqtfg1

      subroutine get_qxc_qyc(mokgm2,qxe_in,qyn_in,flag,qxc,qyc)
      use domain_decomp_1d, only :
     &     getDomainBounds,halo_update,south,north
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : dxpo,dyvo,dxypo
      use ocean, only : lmu,lmm,
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use ocean, only : im,jm,lmo,ivnp,sinic,cosic,sinu,cosu
      implicit none
      real*8, dimension(lmo,im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mokgm2,      ! units: kg/m2
     &     qxe_in,qyn_in, ! extensive or intesive units
     &     qxc,qyc        ! outputs have intensive units
      integer flag ! 1: q_in is extensive; 0: q_in is intensive

      real*8, dimension(lmo,im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     qxe,qyn
      integer :: i,j,l,n
      integer :: j_0s,j_1s,j_0,j_1
      logical :: have_north_pole
      real*8 :: unp,vnp

      call getdomainbounds(grid,
     &     j_strt=j_0, j_stop=j_1,
     &     j_strt_skp=j_0s, j_stop_skp=j_1s,
     &     have_north_pole=have_north_pole)

      if(flag.eq.1) then ! q_in is extensive
        ! convert q to intensive units
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          qxe(l,i,j) = qxe_in(l,i,j)/(mokgm2(l,i,j)*dxypo(j))
          qyn(l,i,j) = qyn_in(l,i,j)/(mokgm2(l,i,j)*dxypo(j))
        enddo
        enddo
        enddo
        enddo
      else               ! q_in is intensive
        qxe=qxe_in
        qyn=qyn_in
      endif

      ! average cell-edge gradients to cell centers.
      call halo_update(grid,qyn,from=south)

      qxc = 0.
      qyc = 0.
      do l=1,lmo
      do j=j_0s,j_1s
        do n=1,nbyzm(j,l)
        do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
          qxc(l,i,j) = .5*(qxe(l,i-1,j)+qxe(l,i,j))
          qyc(l,i,j) = .5*(qyn(l,i,j-1)+qyn(l,i,j))
        enddo
        enddo
        i=1
        if(lmm(i,j).ge.l) then
          qxc(l,i,j) = .5*(qxe(l,im,j)+qxe(l,i,j))
          qyc(l,i,j) = .5*(qyn(l,i,j-1)+qyn(l,i,j))
        endif
      enddo
      enddo
 
      if(have_north_pole) then
c at the north pole
        unp = 0.
        vnp = 0.
        j = jm-1
        do l=1,lmo
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              unp = unp - sinic(i)*qyn(l,i,j)
              vnp = vnp + cosic(i)*qyn(l,i,j)
            enddo
          enddo
          unp = unp*2/im
          vnp = vnp*2/im
          do i=1,im
            qxc(l,i,jm) = unp*cosu(i)  + vnp*sinu(i)
            qyc(l,i,jm) = vnp*cosic(i) - unp*sinic(i)
          enddo
c         qx(l,im,jm) = unp   ! as a result of the above loop
c         qx(l,ivnp,jm) = vnp ! as a result of the above loop
        enddo
      endif
      return
      end subroutine get_qxc_qyc

      subroutine get_qxc_qyc_ijl(mokgm2,qxe_in,qyn_in,flag,qxc,qyc)
      use domain_decomp_1d, only :
     &     getDomainBounds,halo_update,south,north
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : dxpo,dyvo,dxypo
      use ocean, only : lmu,lmm,
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use ocean, only : im,jm,lmo,ivnp,sinic,cosic,sinu,cosu
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mokgm2,      ! units: kg/m2
     &     qxe_in,qyn_in, ! extensive or intesive units
     &     qxc,qyc        ! outputs have intensive units
      integer flag ! 1: q_in is extensive; 0: q_in is intensive

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     qxe,qyn
      integer :: i,j,l,n
      integer :: j_0s,j_1s,j_0,j_1
      logical :: have_north_pole
      real*8 :: unp,vnp

      call getdomainbounds(grid,
     &     j_strt=j_0, j_stop=j_1,
     &     j_strt_skp=j_0s, j_stop_skp=j_1s,
     &     have_north_pole=have_north_pole)

      if(flag.eq.1) then ! q_in is extensive
        ! convert q to intensive units
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          qxe(i,j,l) = qxe_in(i,j,l)/(mokgm2(i,j,l)*dxypo(j))
          qyn(i,j,l) = qyn_in(i,j,l)/(mokgm2(i,j,l)*dxypo(j))
        enddo
        enddo
        enddo
        enddo
      else               ! q_in is intensive
        qxe=qxe_in
        qyn=qyn_in
      endif

      ! average cell-edge gradients to cell centers.
      call halo_update(grid,qyn,from=south)

      qxc = 0.
      qyc = 0.
      do l=1,lmo
      do j=j_0s,j_1s
        do n=1,nbyzm(j,l)
        do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
          qxc(i,j,l) = .5*(qxe(i-1,j,l)+qxe(i,j,l))
          qyc(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
        enddo
        enddo
        i=1
        if(lmm(i,j).ge.l) then
          qxc(i,j,l) = .5*(qxe(im,j,l)+qxe(i,j,l))
          qyc(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
        endif
      enddo
      enddo
 
      if(have_north_pole) then
c at the north pole
        unp = 0.
        vnp = 0.
        j = jm-1
        do l=1,lmo
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              unp = unp - sinic(i)*qyn(i,j,l)
              vnp = vnp + cosic(i)*qyn(i,j,l)
            enddo
          enddo
          unp = unp*2/im
          vnp = vnp*2/im
          do i=1,im
            qxc(i,jm,l) = unp*cosu(i)  + vnp*sinu(i)
            qyc(i,jm,l) = vnp*cosic(i) - unp*sinic(i)
          enddo
c         qx(l,im,jm) = unp   ! as a result of the above loop
c         qx(l,ivnp,jm) = vnp ! as a result of the above loop
        enddo
      endif
      return
      end subroutine get_qxc_qyc_ijl

#endif
