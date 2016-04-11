#include "rundeck_opts.h"
#include "hycom_mpi_hacks.h"

#ifndef STANDALONE_HYCOM

#ifndef CUBED_SPHERE

      module hycom_cpler
      USE CONSTANT, only : tf
      USE HYCOM_DIM_GLOB, only : iia,jja,iio,jjo,isp,ifp,ilp,ii,jj,ip
      USE HYCOM_SCALARS, only : flnma2o,flnma2o_s,flnmo2a,flnmo2a_f
     &   ,flnma2o_tau,flnmcoso,lp
      USE DOMAIN_DECOMP_1D, only : dist_grid
      USE HYCOM_DIM, only : ogrid
     &    ,aJ_0, aJ_1, aJ_0H, aJ_1H,
     &      J_0,  J_1,  J_0H,  J_1H
      USE DOMAIN_DECOMP_1D, only : am_i_root,pack_data,unpack_data
      use filemanager, only : findunit
c
      implicit none
      private

      type(dist_grid), pointer, public :: agrid

      public ssta2o,ssto2a,veca2o,flxa2o,flxo2a,veco2a,tempro2a,cpl_wgt
      public ssto2a_global,flxa2o_global

      public nwgta2o,nwgto2a

      public wlista2o, wtaua2o, wlista2o_s, wlisto2a, wlisto2a_f
     .     ,ilista2o_s, jlista2o_s,nlista2o_s, coso, sino
     .     ,ilista2o,  jlista2o,  nlista2o
     .     ,itaua2o,   jtaua2o,   ntaua2o
     .     ,ilisto2a,  jlisto2a,  nlisto2a
     .     ,ilisto2a_f,jlisto2a_f,nlisto2a_f

      integer nwgta2o,nwgto2a
#ifdef ATM4x5
#ifdef HYCOM2deg
      parameter (nwgta2o=18,nwgto2a=39)
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM2deg
      parameter (nwgta2o=36,nwgto2a=19)
#endif
#endif
#ifdef ATM2x2h
#ifdef HYCOM1deg
      parameter (nwgta2o=37,nwgto2a=48)
#endif
#endif
c
      real*8 wlista2o(iio,jjo,nwgta2o),wtaua2o(iio,jjo,nwgta2o)
     .    ,wlista2o_s(iio,jjo,nwgta2o)
     .    ,wlisto2a(iia,jja,nwgto2a), wlisto2a_f(iia,jja,nwgto2a)
     .    ,coso(iio,jjo),sino(iio,jjo)
      integer ilista2o_s(iio,jjo,nwgta2o),jlista2o_s(iio,jjo,nwgta2o)
     .                                 ,nlista2o_s(iio,jjo)
     .       ,ilista2o(iio,jjo,nwgta2o),jlista2o(iio,jjo,nwgta2o)
     .                                 ,nlista2o(iio,jjo)
     .       ,itaua2o(iio,jjo,nwgta2o), jtaua2o(iio,jjo,nwgta2o)
     .                                 , ntaua2o(iio,jjo)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a  (iia,jja,nwgto2a)
     .                                 ,nlisto2a    (iia,jja)
     .       ,ilisto2a_f(iia,jja,nwgto2a),jlisto2a_f(iia,jja,nwgto2a)
     .                                 ,nlisto2a_f  (iia,jja)
      integer iu1,iu2,iu3,iu4,iu5,iu6,iu7,iu8

      contains
      subroutine ssta2o(flda,fldo)
c --- mapping scalar-like field from agcm A grid to ogcm A grid
c     input: flda, output: fldo 
c
      implicit none
      integer i,j,l,n
      real*8, intent(in)  :: flda(iia,jja)
      real*8, intent(out) :: fldo(iio,jjo)
c
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o_s(i,j)
      fldo(i,j)=fldo(i,j)+flda(ilista2o_s(i,j,n),jlista2o_s(i,j,n))
     .                        *wlista2o_s(i,j,n)
 9    continue
 8    continue
c
      return
      end subroutine ssta2o

      subroutine ssto2a(fldo_loc,flda_loc)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo_loc, output: flda_loc
c
      implicit none
      integer n,ia,ja
      real*8 fldo_loc(iio,J_0H:J_1H),flda_loc(iia,aJ_0H:aJ_1H)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo))
      else
        allocate(flda(1,1),fldo(1,1))
      endif
      call pack_data(ogrid,fldo_loc,fldo)
      if(am_i_root()) then
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
 17   flda(ia,ja)=flda(ia,ja)+fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .                                              *wlisto2a(ia,ja,n)
 16   continue
c
      endif ! am_i_root
      call unpack_data(agrid,flda,flda_loc)
      deallocate(flda,fldo)
      return
      end subroutine ssto2a

      subroutine ssto2a_global(fldo,flda)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo, output: flda
c
      implicit none
      integer n,ia,ja
      real*8 flda(iia,jja),fldo(iio,jjo)
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
 17   flda(ia,ja)=flda(ia,ja)+fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .                                              *wlisto2a(ia,ja,n)
 16   continue
c
      return
      end subroutine ssto2a_global
c
c
      subroutine veca2o(tauxa_loc,tauya_loc,tauxo_loc,tauyo_loc)
c --- mapping vector like stress from agcm to ogcm, both on A grid
c --- input  tauxa/tauya: E_/N_ward on agcm A grid
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercador domain)
c
      implicit none
      integer i,j,l,n
      real*8, dimension(iia,aJ_0H:aJ_1H) :: tauxa_loc,tauya_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real*8, dimension(:,:), allocatable :: tauxa,tauya,tauxo,tauyo,
     &     sward,eward
      if(am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       sward(iio,jjo),eward(iio,jjo)
     &       )
      else
        allocate(
     &       tauxa(1,1),tauya(1,1),
     &       tauxo(1,1),tauyo(1,1),
     &       sward(1,1),eward(1,1)
     &       )
      endif
      call pack_data(agrid,tauxa_loc,tauxa)
      call pack_data(agrid,tauya_loc,tauya)
      if(am_i_root()) then
c
c --- mapping tauxa/tauya to ogcm grid
      do 6 j=1,jj               
      do 6 l=1,isp(j)           
      do 6 i=ifp(j,l),ilp(j,l)
      eward(i,j)=0.
      sward(i,j)=0.
c
      do 7 n=1,ntaua2o(i,j)
      eward(i,j)=eward(i,j)+tauxa(itaua2o(i,j,n),jtaua2o(i,j,n))
     .                                          *wtaua2o(i,j,n)
 7    sward(i,j)=sward(i,j)-tauya(itaua2o(i,j,n),jtaua2o(i,j,n))
     .                                          *wtaua2o(i,j,n)
 6    continue
c
c --- rotate sward/eward to fit onto Panam grid
      do 9 j=1,jj
      do 9 l=1,isp(j)           
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue
      endif ! am_i_root
      call unpack_data(ogrid,tauxo,tauxo_loc)
      call unpack_data(ogrid,tauyo,tauyo_loc)
      deallocate(tauxa,tauya,tauxo,tauyo,sward,eward)
      return
      end subroutine veca2o
c
      subroutine flxa2o(flda_loc,fldo_loc)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
      integer i,j,l,n
      real*8 flda_loc(:,:),fldo_loc(:,:)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo))
      else
        allocate(flda(1,1),fldo(1,1))
      endif
      call pack_data(agrid,flda_loc,flda)
c
      if(am_i_root()) then
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o(i,j)
      fldo(i,j)=fldo(i,j)+flda(ilista2o(i,j,n),jlista2o(i,j,n))
     .                        *wlista2o(i,j,n)
 9    continue
 8    continue
      endif ! am_i_root
c
      call unpack_data(ogrid,fldo,fldo_loc)
      deallocate(flda,fldo)
      return
      end subroutine flxa2o

      subroutine flxa2o_global(flda,fldo)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
      integer i,j,l,n
      real*8, intent(in)  :: flda(iia,jja)
      real*8, intent(out) :: fldo(iio,jjo)
c
      do 8 j=1,jj
      do 8 l=1,isp(j)
      do 8 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
c
      do 9 n=1,nlista2o(i,j)
      fldo(i,j)=fldo(i,j)+flda(ilista2o(i,j,n),jlista2o(i,j,n))
     .                        *wlista2o(i,j,n)
 9    continue
 8    continue
c
      return
      end subroutine flxa2o_global
c
      subroutine flxo2a(fldo,flda)
c --- mapping flux-like field from ogcm A grid to agcm A grid
c     input: fldo (W/m*m), output: flda (W/m*m)
c
      implicit none
      real*8, intent(in)  :: fldo(iio,jjo)
      real*8, intent(out) :: flda(iia,jja)
      integer n,ia,ja
c
      do 8 ja=1,jja
      do 8 ia=1,iia
      flda(ia,ja)=0.
c
      do 9 n=1,nlisto2a_f(ia,ja)
      flda(ia,ja)=flda(ia,ja)+ 
     . fldo(ilisto2a_f(ia,ja,n),jlisto2a_f(ia,ja,n))*wlisto2a_f(ia,ja,n)
 9    continue
 8    continue
c
      return
      end subroutine flxo2a
c
      subroutine veco2a(tauxo_loc,tauyo_loc,tauxa_loc,tauya_loc)
c --- mapping vector like velocity from C grid ogcm to A grid agcm
c --- input  tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxa/tauya: E_/N_ward on agcm A grid
c
      implicit none
      integer i,j,l,n,ia,ja,jb
      real*8, dimension(iia,aJ_0H:aJ_1H) :: tauxa_loc,tauya_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real*8, dimension(:,:), allocatable :: tauxa,tauya,tauxo,tauyo,
     &     nward,eward
      real*8 sine
      if(am_i_root()) then
        allocate(
     &       tauxa(iia,jja),tauya(iia,jja),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       nward(iio,jjo),eward(iio,jjo)
     &       )
      else
        allocate(
     &       tauxa(1,1),tauya(1,1),
     &       tauxo(1,1),tauyo(1,1),
     &       nward(1,1),eward(1,1)
     &       )
      endif
      call pack_data(ogrid,tauxo_loc,tauxo)
      call pack_data(ogrid,tauyo_loc,tauyo)
      if(am_i_root()) then
c
      do 10 j=1,jj
      do 10 i=1,ii
      nward(i,j)=0.
 10   eward(i,j)=0.
c
c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)           
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c
c --- weights are for mapping nward/eward from ogcm to agcm, both on A grid
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      tauxa(ia,ja)=0.
      tauya(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
      tauxa(ia,ja)=tauxa(ia,ja)+eward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)
 17   tauya(ia,ja)=tauya(ia,ja)+nward(ilisto2a(ia,ja,n)
     .            ,jlisto2a(ia,ja,n))*wlisto2a(ia,ja,n)
 16   continue
c
      endif ! am_i_root
      call unpack_data(agrid,tauxa,tauxa_loc)
      call unpack_data(agrid,tauya,tauya_loc)
      deallocate(tauxa,tauya,tauxo,tauyo,nward,eward)
      return
      end subroutine veco2a
c
      subroutine tempro2a(fldo_loc,flda_loc)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: fldo in deg C; outout: flda in deg K
c
      implicit none
      integer n,ia,ja
      real*8 fldo_loc(iio,J_0H:J_1H),flda_loc(iia,aJ_0H:aJ_1H)
      real*8, allocatable :: flda(:,:),fldo(:,:)
      if(am_i_root()) then
        allocate(flda(iia,jja),fldo(iio,jjo))
      else
        allocate(flda(1,1),fldo(1,1))
      endif
      call pack_data(ogrid,fldo_loc,fldo)
      if(am_i_root()) then
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
 17   flda(ia,ja)=flda(ia,ja)+(fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .    +tf)**4*wlisto2a(ia,ja,n)
      flda(ia,ja)=sqrt(sqrt(flda(ia,ja)))       ! Kelvin for radiation
 16   continue
c
      endif ! am_i_root
      call unpack_data(agrid,flda,flda_loc)
      deallocate(flda,fldo)
      return
      end subroutine tempro2a
c
      subroutine cpl_wgt
      implicit none
      integer :: iz,jz
c
#ifdef  ATM4x5
#ifdef  HYCOM2deg
       integer, parameter :: nsize1=10249200, nsize2=2079936
#endif     
#endif     
#ifdef ATM2x2h
#ifdef HYCOM2deg
       integer, parameter :: nsize1=20358000, nsize2=3991680
#endif
#endif     
#ifdef ATM2x2h
#ifdef HYCOM1deg
       integer, parameter :: nsize1=83034720, nsize2=10005120
#endif
#endif     

       if(.not. am_i_root()) return

c --- read in all weights
      if (iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8).ne.nsize1 .or.
     .    iia*jja*((nwgto2a*2+1)*4+nwgto2a*8).ne.nsize2) then
        write(lp,'(a,2i12,a,2i12)') 'wrong size in cpler '
     .  ,iio*jjo*((nwgta2o*2+1)*4+nwgta2o*8)
     .  ,iia*jja*((nwgto2a*2+1)*4+nwgto2a*8)
     .  ,' should be ',nsize1,nsize2
        stop ' wrong size in cpler'
      endif
c
      call findunit(iu1)
      open(iu1,file=flnma2o,form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(iu1,rec=1) ilista2o,jlista2o,wlista2o,nlista2o
      close(iu1)
c
      call findunit(iu2)
      open(iu2,file=flnma2o_tau,form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(iu2,rec=1) itaua2o,jtaua2o,wtaua2o,ntaua2o
      close(iu2)
c
      call findunit(iu3)
      open(iu3,file=flnmo2a,form='unformatted',status='old',
     .  access='direct',recl=nsize2)
      read(iu3,rec=1) ilisto2a,jlisto2a,wlisto2a,nlisto2a
      close(iu3)
c
      call findunit(iu4)
      open(iu4,file=flnmcoso,form='unformatted',status='old')
      read(iu4) iz,jz,coso,sino
      close(iu4)
      if (iz.ne.iio .or. jz.ne.jjo) then
        write(lp,*) ' iz,jz=',iz,jz
        stop '(wrong iz/jz in cososino.8bin)'
      endif
c
#ifdef ATM2x2h
#ifdef HYCOM1deg
      call findunit(iu5)
      open(iu5,file=flnmo2a_f,form='unformatted',status='old',     ! TNL
     .  access='direct',recl=nsize2)
      read(iu5,rec=1) ilisto2a_f,jlisto2a_f,wlisto2a_f,nlisto2a_f
      close(iu5)
c
      call findunit(iu6)
      open(iu6,file=flnma2o_s,form='unformatted',status='old',   ! TNL
     .  access='direct',recl=nsize1)
      read(iu6,rec=1) ilista2o_s,jlista2o_s,wlista2o_s,nlista2o_s
      close(iu6)
#endif
#endif
c
      return
      end subroutine cpl_wgt

      end module hycom_cpler

      module hycom_dynsi_cpler  ! latlon atm not using this yet
      contains
      subroutine init_hycom_dynsi_cpler
      end subroutine init_hycom_dynsi_cpler
      end module hycom_dynsi_cpler

#else  /* coupler with cubed-sphere atm is in cplercs.f */

      module hycom_dynsi_cpler
      USE HYCOM_DIM_GLOB, only : iio,jjo,isp,ifp,ilp,ii,jj,ip
      USE DOMAIN_DECOMP_1D, only : am_i_root,get,pack_data,unpack_data
      USE HYCOM_DIM, only : ogrid, J_0H, J_1H
      USE ICEDYN, only : igrid=>grid_ICDYN,iii=>imicdyn,jji=>jmicdyn
      USE HYCOM_CPLER, only : coso=>coso_glob, sino=>sino_glob
      USE HYCOM_ATM, only : dmui_loc,dmvi_loc
      use filemanager, only : findunit
      implicit none
      private
      save
      public :: init_hycom_dynsi_cpler,veci2o,veco2i,scai2o
      public :: reset_dynsi_accum,do_dynsi_accum,
     &     taux_dynsi,tauy_dynsi,ustar_dynsi

      real*8, dimension(:,:,:), allocatable :: wlisto2i,wlisti2o
      integer, dimension(:,:,:), allocatable ::
     &     ilisto2i,jlisto2i, ilisti2o,jlisti2o
      integer, dimension(:,:), allocatable :: nlisto2i,nlisti2o

      real*8, dimension(:,:), allocatable ::
     &     taux_dynsi,tauy_dynsi,ustar_dynsi

      contains

      subroutine reset_dynsi_accum
      integer :: i_j0h,i_j1h
      if(.not.igrid%have_domain) return
      if(.not.allocated(taux_dynsi)) then
        i_j0h = igrid%j_strt_halo
        i_j1h = igrid%j_stop_halo
        allocate(taux_dynsi(iii,i_j0h:i_j1h)
     &          ,tauy_dynsi(iii,i_j0h:i_j1h)
     &          ,ustar_dynsi(iii,i_j0h:i_j1h)
     &       )
      endif
      taux_dynsi = 0.
      tauy_dynsi = 0.
      ustar_dynsi = 0.
      return
      end subroutine reset_dynsi_accum

      subroutine do_dynsi_accum(dt,dtacc)
      use hycom_scalars, only : thref
      real*8 :: dt,dtacc
      integer :: i,j
      if(.not.igrid%have_domain) return
      do j=igrid%j_strt,igrid%j_stop
        do i=1,iii
c taux,tauy are time integrals over the time interval dt
          taux_dynsi(i,j) = taux_dynsi(i,j) + dmui_loc(i,j)/dtacc
          tauy_dynsi(i,j) = tauy_dynsi(i,j) + dmvi_loc(i,j)/dtacc
          ustar_dynsi(i,j) = ustar_dynsi(i,j) + (dt/dtacc)*
     &         sqrt(sqrt(dmui_loc(i,j)**2 +dmvi_loc(i,j)**2)*thref/dt)
        enddo
      enddo
      return
      end subroutine do_dynsi_accum

      subroutine init_hycom_dynsi_cpler
      integer :: iu
      integer :: nwgti2o,nwgto2i,nsize1,nsize2

      if(.not.am_i_root()) return

      if(iio==387 .and. jjo==360 .and. iii==180 .and. jji==180) then
        nwgti2o=46; nwgto2i=24
        nsize1=103096800; nsize2=12571200 ! file sizes on disk
      else
        call stop_model('hycom_dynsi_cpler: unknown resolution',255)
      endif

      if (iio*jjo*((nwgti2o*2+1)*4+nwgti2o*8).ne.nsize1 .or.
     .    iii*jji*((nwgto2i*2+1)*4+nwgto2i*8).ne.nsize2) then
        write(6,'(a,2i12,a,2i12)') 'wrong size in cpler '
     .  ,iio*jjo*((nwgti2o*2+1)*4+nwgti2o*8)
     .  ,iii*jji*((nwgto2i*2+1)*4+nwgto2i*8)
     .  ,' should be ',nsize1,nsize2
        call stop_model(' wrong size in hycom_dynsi_cpler',255)
      endif

      allocate(wlisto2i(iii,jji,nwgto2i)
     &        ,ilisto2i(iii,jji,nwgto2i)
     &        ,jlisto2i(iii,jji,nwgto2i)
     &        ,nlisto2i(iii,jji))
      allocate(wlisti2o(iio,jjo,nwgti2o)
     &        ,ilisti2o(iio,jjo,nwgti2o)
     &        ,jlisti2o(iio,jjo,nwgti2o)
     &        ,nlisti2o(iio,jjo))

      call findunit(iu)
      open(iu,file='taui2o',form='unformatted',status='old',
     .  access='direct',recl=nsize1)
      read(iu,rec=1) ilisti2o,jlisti2o,wlisti2o,nlisti2o
      close(iu)
c
      call findunit(iu)
      open(iu,file='tauo2i',form='unformatted',status='old',
     .  access='direct',recl=nsize2)
      read(iu,rec=1) ilisto2i,jlisto2i,wlisto2i,nlisto2i
      close(iu)

      return
      end subroutine init_hycom_dynsi_cpler

      subroutine veci2o(tauxi_loc,tauyi_loc,tauxo_loc,tauyo_loc)
c --- mapping vector like stress from B grid ice model to A grid ogcm
c --- input tauxi/tauyi: E_/N_ward on ice B grid located at i+1/2 & j+1/2 corner
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercador domain)
c
      implicit none
      integer i,j,l,n
      real*8, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     tauxi_loc,tauyi_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc
      real*8, dimension(:,:), allocatable :: tauxi,tauyi,tauxo,tauyo,
     &     sward,eward
      if(am_i_root()) then
        allocate(
     &       tauxi(iii,jji),tauyi(iii,jji),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       sward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      if(igrid%have_domain) then
        call pack_data(igrid,tauxi_loc,tauxi)
        call pack_data(igrid,tauyi_loc,tauyi)
      endif
      if(am_i_root()) then
c
c --- mapping B-grid tauxi/tauyi to A-grid ogcm
      do 6 j=1,jj               
      do 6 l=1,isp(j)           
      do 6 i=ifp(j,l),ilp(j,l)
      eward(i,j)=0.
      sward(i,j)=0.
c
      do 7 n=1,nlisti2o(i,j)
      eward(i,j)=eward(i,j)+tauxi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                           *wlisti2o(i,j,n)
 7    sward(i,j)=sward(i,j)-tauyi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                           *wlisti2o(i,j,n)
 6    continue
c
c --- rotate sward/eward to fit onto Panam grid
      do 9 j=1,jj
      do 9 l=1,isp(j)           
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue
      endif ! am_i_root
      call unpack_data(ogrid,tauxo,tauxo_loc)
      call unpack_data(ogrid,tauyo,tauyo_loc)
      if(am_i_root()) then
        deallocate(tauxi,tauyi,tauxo,tauyo,sward,eward)
      endif
      return
      end subroutine veci2o
c
      subroutine veco2i(tauxo_loc,tauyo_loc,tauxi_loc,tauyi_loc)
c --- mapping vector like velocity from C grid ocean model to B grid ice model
c --- input tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxi/tauyi: E_/N_ward on ice B grid
c
      implicit none
      integer i,j,l,n,ia,ja,jb
      real*8, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     tauxi_loc,tauyi_loc
      real*8, dimension(iio, J_0H: J_1H) :: tauxo_loc,tauyo_loc

      real*8, dimension(:,:), allocatable :: tauxi,tauyi,tauxo,tauyo,
     &     nward,eward
      real*8 sine
      if(am_i_root()) then
        allocate(
     &       tauxi(iii,jji),tauyi(iii,jji),
     &       tauxo(iio,jjo),tauyo(iio,jjo),
     &       nward(iio,jjo),eward(iio,jjo)
     &       )
      endif
      call pack_data(ogrid,tauxo_loc,tauxo)
      call pack_data(ogrid,tauyo_loc,tauyo)
      if(am_i_root()) then
      nward(:,:)=0.
      eward(:,:)=0.
c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=1,jj
      jb=mod(j,jj)+1
      do 12 l=1,isp(j)           
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c
c --- mapping nward/eward from A-grid ogcm to B-grid ice model
c
      do 16 j=1,jji
      do 16 i=1,iii
      tauxi(i,j)=0.
      tauyi(i,j)=0.
c
      do 17 n=1,nlisto2i(i,j)
      tauxi(i,j)=tauxi(i,j)+eward(ilisto2i(i,j,n)
     .          ,jlisto2i(i,j,n))*wlisto2i(i,j,n)
 17   tauyi(i,j)=tauyi(i,j)+nward(ilisto2i(i,j,n)
     .          ,jlisto2i(i,j,n))*wlisto2i(i,j,n)
 16   continue
c
      endif ! am_i_root
      if(igrid%have_domain) then
        call unpack_data(igrid,tauxi,tauxi_loc)
        call unpack_data(igrid,tauyi,tauyi_loc)
      endif
      if(am_i_root()) then
        deallocate(tauxi,tauyi,tauxo,tauyo,nward,eward)
      endif
      return
      end subroutine veco2i

      subroutine scai2o(fldi_loc,fldo_loc)
c --- mapping scalar from B grid ice model to A grid ogcm
c --- input fldi: quantity on ice B grid located at i+1/2 & j+1/2 corner
c --- output fldo: quantity on ogcm A grid
c
      implicit none
      integer i,j,l,n
      real*8, dimension(iii,igrid%J_STRT_HALO:igrid%J_STOP_HALO) ::
     &     fldi_loc
      real*8, dimension(iio, J_0H: J_1H) :: fldo_loc
      real*8, dimension(:,:), allocatable :: fldi,fldo
      if(am_i_root()) then
        allocate(
     &       fldi(iii,jji),
     &       fldo(iio,jjo)
     &       )
      endif
      if(igrid%have_domain) then
        call pack_data(igrid,fldi_loc,fldi)
      endif
      if(am_i_root()) then
c
c --- mapping B-grid fldi to A-grid ogcm
      do 6 j=1,jj               
      do 6 l=1,isp(j)           
      do 6 i=ifp(j,l),ilp(j,l)
      fldo(i,j)=0.
      do n=1,nlisti2o(i,j)
        fldo(i,j)=fldo(i,j)+fldi(ilisti2o(i,j,n),jlisti2o(i,j,n))
     .                                          *wlisti2o(i,j,n)
      enddo
 6    continue
      endif ! am_i_root
      call unpack_data(ogrid,fldo,fldo_loc)
      if(am_i_root()) then
        deallocate(fldi,fldo)
      endif
      return
      end subroutine scai2o

      end module hycom_dynsi_cpler

#endif /* cubed-sphere versus lat-lon atm */

#else /* code below is a trivial coupler for the standalone ocean */

      module hycom_cpler
      USE HYCOM_DIM_GLOB, only : iio,jjo,jj
      USE DOMAIN_DECOMP_1D, only : dist_grid
      USE HYCOM_DIM, only : ogrid, J_0,J_1, J_0H,J_1H
     &    ,isp,ifp,ilp,ip
c
      implicit none
      private

      type(dist_grid), pointer, public :: agrid

      public ssto2a,veca2o,flxa2o,veco2a,tempro2a,cpl_wgt

      public coso_glob, sino_glob

      real*8, dimension(:,:), allocatable :: coso,sino
      real*8, dimension(:,:), allocatable :: coso_glob,sino_glob

      contains

      subroutine ssto2a(fldo,flda)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo, output: flda
c
      implicit none
      real*8 fldo(iio,J_0H:J_1H),flda(iio,J_0H:J_1H)
      flda(:,:) = fldo(:,:)
      return
      end subroutine ssto2a

      subroutine veca2o(tauxa,tauya,tauxo,tauyo)
c --- mapping vector like stress from agcm to ogcm, both on A grid
c --- input  tauxa/tauya: E_/N_ward on agcm A grid
c --- output tauxo/tauyo: +i_/+j_ward on ogcm A grid (S_/E_ward in Mercador domain)
c
      implicit none
      real*8, dimension(iio,J_0H:J_1H) :: tauxa,tauya,tauxo,tauyo
      real*8, dimension(iio,J_0H:J_1H) :: sward,eward
      integer i,j,l,n

      eward(:,:) = +tauxa(:,:)
      sward(:,:) = -tauya(:,:)

c
c --- rotate sward/eward to fit onto Panam grid
      do 9 j=j_0,j_1
      do 9 l=1,isp(j)           
      do 9 i=ifp(j,l),ilp(j,l)
      tauxo(i,j)= sward(i,j)*coso(i,j)+eward(i,j)*sino(i,j)
      tauyo(i,j)= eward(i,j)*coso(i,j)-sward(i,j)*sino(i,j)
 9    continue

      end subroutine veca2o
c
      subroutine flxa2o(flda,fldo)
c --- mapping flux-like field from agcm A grid to ogcm A grid
c     input: flda (W/m*m), output: fldo (W/m*m)
c
      implicit none
      real*8 flda(iio,J_0H:J_1H),fldo(iio,J_0H:J_1H)
      fldo(:,:) = flda(:,:)
      return
      end subroutine flxa2o

      subroutine veco2a(tauxo,tauyo,tauxa,tauya)
c --- mapping vector like velocity from C grid ogcm to A grid agcm
c --- input  tauxo/tauyo: +i_/+j_ward (S_/E_ward in Mercador domain) on ocean C grid (@ i-1/2 & j-1/2)
c --- output tauxa/tauya: E_/N_ward on agcm A grid
c
      use domain_decomp_1d, only : halo_update
      implicit none
      integer i,j,l,n,ia,ja,jb
      real*8, dimension(iio,J_0H:J_1H) :: tauxo,tauyo,tauxa,tauya
      real*8 sine

      call halo_update(ogrid,tauyo)

      tauxa(:,:) = 0.
      tauya(:,:) = 0.

c --- average tauxo/tauyo from C to A grid & rotate to n/e orientation at A grid
c --- check velocity bounds
      do 12 j=j_0,j_1
      jb = PERIODIC_INDEX(j+1, jj) !mod(j,jj)+1
      do 12 l=1,isp(j)           
      do 12 i=ifp(j,l),ilp(j,l)
      if (ip(i,j).eq.1) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      tauya(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      tauxa(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue

      return
      end subroutine veco2a
c
      subroutine tempro2a(fldo,flda)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: fldo in deg C; outout: flda in deg K
c
      implicit none
      real*8 fldo(iio,J_0H:J_1H),flda(iio,J_0H:J_1H)
      flda(:,:) = fldo(:,:)+tf
      return
      end subroutine tempro2a
c
      subroutine cpl_wgt
      USE HYCOM_SCALARS, only : flnmcoso,lp
      use filemanager, only : findunit
      USE DOMAIN_DECOMP_1D, only : am_i_root,unpack_data
      implicit none
      integer :: iu1,iz,jz

c read rotation coeffs between hycom gridlines and geographic north
      allocate(coso(iio,j_0h:j_1h),sino(iio,j_0h:j_1h))
      if(am_i_root()) then
        allocate(coso_glob(iio,jjo),sino_glob(iio,jjo))
        call findunit(iu1)
        open(iu1,file=flnmcoso,form='unformatted',status='old')
        read(iu1) iz,jz,coso_glob,sino_glob
        close(iu1)
        if (iz.ne.iio .or. jz.ne.jjo) then
          write(6,*) ' iz,jz=',iz,jz
          stop '(wrong iz/jz in cososino.8bin)'
        endif
      endif
      call unpack_data(ogrid,coso_glob,coso)
      call unpack_data(ogrid,sino_glob,sino)

      return
      end subroutine cpl_wgt

      end module hycom_cpler

      module hycom_dynsi_cpler
      contains
      subroutine init_hycom_dynsi_cpler
      end subroutine init_hycom_dynsi_cpler
      end module hycom_dynsi_cpler


#endif /* standalone ocean */
