#include "rundeck_opts.h"
  
      module lightning
!@sum  lightning variables for lightning parameterization
!@auth Colin Price (modelEification by Greg Faluvegi)
      use resolution, only : LM

      implicit none
      save
      
!@param JNlight NH latitude boundary where lightning changes  
!@param JSlight SH latitude boundary where lightning changes  
      real*8, parameter :: JNlight=30.d0, JSlight=-30.d0 
!@dbparam tune_lt_land multiplier of flash rate over land to match obs
!@dbparam tune_lt_sea multiplier of flash rate over ocean to match obs
! Set these in the rundeck and note that they are a fuction of the 
! horizontal resolution of the GCM: See Price and Rind (1994) Mon.Wea.Rev.
! Values here default for 4x5 model:
      real*8:: tune_lt_land=4.774d0 ! =2.2d0*2.17d0 for 4x5 model
     &        ,tune_lt_sea= 8.463d0 ! =3.9d0*2.17d0 for 4x5 model
#ifdef TRACERS_SPECIAL_Shindell
!@var RNOx_lgt NOx production rate from lightning
!@var srclight NOx source from lightning column array
!@var HGT_lgt Pickering vertical lightning distributions (1998)
      real*8, allocatable, dimension(:,:) :: RNOx_lgt
      real*8, dimension(LM) :: srclight
      real*8, dimension(2,2,16) :: HGT_lgt
      integer :: i
      
      data (HGT_lgt(1,1,i),i=1,16)/8.2,1.9,2.1,1.6,1.1,1.6,3.0,5.8,
     &  7.6,9.6,10.5,12.3,11.8,12.5,8.1,2.3/
      data (HGT_lgt(1,2,i),i=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/
      data (HGT_lgt(2,1,i),i=1,16)/20.1,2.3,0.8,1.5,3.4,5.3,3.6,3.8,
     &    5.4,6.6,8.3,9.6,12.8,10.0,6.2,0.3/
      data (HGT_lgt(2,2,i),i=1,16)/5.8,2.9,2.6,2.4,2.2,2.1,2.3,6.1,
     & 16.5,14.1,13.7,12.8,12.5,2.8,0.9,0.3/ 
#endif
!@var saveLightning array to save flash rate for SUBDD (instantaneous)
!@var saveC2gLightning to save CtoG flash rate for SUBDD (instantaneous)
      real*8,allocatable,dimension(:,:)::saveLightning,saveC2gLightning

      end module lightning 
      
      
      subroutine alloc_lightning(grid)
!@SUM  alllocate lightning arrays for current grid
!@auth G.Faluvegi
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use LIGHTNING, only : saveC2gLightning,saveLightning
#ifdef TRACERS_SPECIAL_Shindell
      use LIGHTNING, only     : RNOx_lgt
#endif
      implicit none

      type (dist_grid), intent(in) :: grid
      integer :: J_1H, J_0H, I_0H, I_1H

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &                 I_STRT_HALO=I_0H, I_STOP_HALO=I_1H )
 
#ifdef TRACERS_SPECIAL_Shindell
      allocate( RNOx_lgt(I_0H:I_1H,J_0H:J_1H) )
#endif
      allocate( saveC2gLightning(I_0H:I_1H,J_0H:J_1H) )
      allocate(    saveLightning(I_0H:I_1H,J_0H:J_1H) )
  
      return
      end subroutine alloc_lightning
     
 
      subroutine calc_lightning(i,j,lmax,lfrz)
!@sum calc_lightning calculates lightning flash amount and cloud-
!@+   to-ground amount, based on cloud top height. WARNING: this 
!@+   routine is apparently resolution dependant.  See the comments.
!@auth Colin Price (modelEifications by Greg Faluvegi)

      use TimeConstants_mod, only: SECONDS_PER_MINUTE
      use lightning, only : JNlight,JSlight,tune_lt_land,tune_lt_sea
     & ,saveC2gLightning,saveLightning
#ifdef TRACERS_SPECIAL_Shindell
     & ,RNOx_lgt
#endif
      use model_com, only : DTsrc
      use geom,      only : lat2d_dg,axyp,byaxyp
      use constant,  only : bygrav
      use atm_com,   only : gz
      use diag_com,  only : ij_CtoG,ij_flash,aij=>aij_loc
      use fluxes, only : fland
      implicit none

!@var lmax highest layer of current convective cloud event
!@var lmax_temp local copy of LMAX (alterable)
!@var lfrz freezing level
!@var htcon height of convection?
!@var htfrz height of freezing?
!@var flash lightning flashes per minute
!@var th,th2,th3,th4 thickness of cold sector, squared, cubed, etc.
!@var cg fraction of lightning that is cloud-to-ground
!@var zlt ?
      integer, intent(in) :: lmax,lfrz,i,j
      integer:: lmax_temp
      real*8 :: htcon,htfrz,flash,th,th2,th3,th4,zlt,cg,area_ref
!@param tune_NOx multiplier of NOx production rate from lightning
      real*8, parameter :: tune_NOx=0.500d0
      real*8 :: minPerTimeStep

      minPerTimeStep=DTsrc/SECONDS_PER_MINUTE
    
! The folowing simple algorithm calculates the lightning
! frequency in each gridbox using the moist convective cloud
! scheme.  The algorithm is based on the Price and Rind (1992)
! JGR paper, and differentiates between oceanic and continental
! convective clouds.  Only the non-entraining plume is used for
! the lightning.  However, I [Greg] have removed the IC variable
! from this routine because Gavin says: "IC distinguishes between
! entraining and non-entraining, but since non-entraining always
! go higher, those are the ones represented by LMCMAX."
!
! The lightning calculation uses the maximum height of the cloud,
! or the maximum depth of the mass flux penetration, LMAX, at
! every timestep.  In each timestep there may be a number of
! different values of LMAX, associated with clouds originating
! from different levels.  However, there are always less values
! of LMAX then the number of vertical layers. gz is the
! geopotential height, used from DYNAMICS module.
 
      lmax_temp=lmax
      if(lat2d_dg(i,j)<JSlight.or.lat2d_dg(i,j)>JNlight)
     & lmax_temp=lmax+1 ! extra-tropics             
      htcon=gz(i,j,lmax_temp)*bygrav*1.d-3
      htfrz=gz(i,j,lfrz)*bygrav*1.d-3

! If the gridbox is over land or coastline use the continental
! parameterization. We use a threshold of 1% for continental
! gridboxes. The units are flashes per minute.

      If (fland(i,j) <= 0.01) then
        flash= tune_lt_sea*6.2d-4*(htcon**1.73d0)  ! ocean
      else
        flash=tune_lt_land*3.44d-5*(htcon**4.92d0) ! continent
      end if

#ifdef CUBED_SPHERE
! rescale flash rate by gridbox area
      area_ref = 6d10 ! avg 2x2.5 tropical area for reference
      flash = flash * axyp(i,j)/area_ref
#endif

! The formulation by Price and Rind (1993) in Geophysical Res.
! Letters is used to calculate the fraction of total lightning
! that is cloud-to-ground lightning.  This uses the thickness of
! cold sector of the cloud (Hmax-Hzero) as the determining
! parameter.  The algorithm is only valid for thicknesses greater
! than 5.5km and less than 14km.

      th=(htcon-htfrz)
      th=min(max(th,5.5d0),14.d0)
      th2=th*th
      th3=th2*th
      th4=th3*th
      zlt = 0.021d0*th4 - 0.648d0*th3 + 7.493d0*th2 
     &    - 36.544d0*th + 63.088d0
      cg=flash/(1.+zlt)

! *If* flash is indeed in flashes/min, accumulate it (for aij)
! in flashes/m2. It will later be divided by DTsrc and become
! flashes/m2/s:
      aij(i,j,ij_flash)= aij(i,j,ij_flash) + 
     & flash*minPerTimeStep*byaxyp(i,j)
      aij(i,j,ij_CtoG) = aij(i,j,ij_CtoG)  + 
     &    cg*minPerTimeStep*byaxyp(i,j)
! Save these for SUBDDdiags instantaneous output (in flashes/m2/s),
! and (in the case of saveC2gLightning) for igniting fires:
      saveLightning(i,j)    =  flash*minPerTimeStep*byaxyp(i,j)/DTsrc
      saveC2gLightning(i,j) =     cg*minPerTimeStep*byaxyp(i,j)/DTsrc

#ifdef TRACERS_SPECIAL_Shindell
! Given the number of cloud-to-ground flashes, we can now calculate
! the NOx production rate based on the Price et al. (1997) JGR paper.
! The units are grams of Nitrogen per minute:

      RNOx_lgt(i,j)=15611.d0*(cg + 0.1d0*(flash-cg))*tune_NOx
#endif
      end subroutine calc_lightning

             

#ifdef TRACERS_SPECIAL_Shindell
      subroutine get_lightning_NOx
!@sum  get_lightning_NOx to define the 3D source of NOx from lightning
!@auth Colin Price / Greg Faluvegi
 
      use TimeConstants_mod, only: SECONDS_PER_MINUTE
      use geom, only       : lat2d_dg,byaxyp,imaxj
      use fluxes, only     : tr3Dsource,fland
      use tracer_com, only : n_NOx,nOther
      use lightning, only  : HGT_lgt,JSlight,JNlight,srclight,RNOx_lgt
      use constant, only   : bygrav
      use resolution, only  : LM
      use atm_com, only    : ltropo, phi
      use domain_decomp_atm, only : GRID, getDomainBounds
#ifdef ACCMIP_LIKE_DIAGS
      use trdiag_com, only : taijls=>taijls_loc,ijlt_NOxLgt
#endif
 
      implicit none
 
!@var latindx Pickering latitude index
!@var landindx Pickering surface type index
!@var ih Pickering altitude index
!@var height Pickering HGT_lgt variable specific to 3D point
!@var levtrop local variable to hold the tropopause level
!@var alttrop altitude of tropopause
!@var alttop altitude at the top of ?
!@param pmin2psec to convert from min-1 to sec-1
      real*8, parameter :: pmin2psec = 1.d0/SECONDS_PER_MINUTE
      integer:: latindx,landindx,ih,levtrop,i,j,L
      real*8, dimension(16):: height
      real*8:: alttrop,alttop
      
      integer :: J_1, J_0, J_0H, J_1H, I_0, I_1

      I_0 = grid%I_STRT; J_0 = grid%J_STRT
      I_1 = grid%I_STOP; J_1 = grid%J_STOP
 
      do j=J_0,J_1
      do i=I_0,imaxj(j)
! Lightning source function altitude dependence:
! Determine if latitude is tropical:
         latindx=2
         if(lat2d_dg(i,j)>JSlight.or.lat2d_dg(i,j)<=JNlight)latindx=1
! Determine if over land or ocean:
         landindx=2
         if(fland(i,j) >= 0.5) landindx=1
! Choose appropriate height file:
         height(:)=HGT_lgt(latindx,landindx,:)*0.01d0
! Store local tropopause
         levtrop=ltropo(i,j)
         alttrop=phi(i,j,levtrop)*bygrav*1.d-3
         if(alttrop==0.)call stop_model('lightn 1',255)! ?
         srclight(:) = 0.d0
! Determine which GCM level the height belongs, and sum NOx source:
         ih_loop: do ih=1,16
           L=1
           do
             alttop=(phi(i,j,l)+(phi(i,j,l+1)-
     &       phi(i,j,l))*0.5d0)*bygrav*1.d-3
             if(ih <= alttop)then
               srclight(L)=srclight(L)+RNOx_lgt(i,j)*height(ih)
               cycle ih_loop
             elseif(ih-1 < alttop)then
               srclight(L)=srclight(L)+RNOx_lgt(I,J)*height(ih)*
     &         (alttop-ih+1.d0)
               srclight(L+1)=srclight(L+1)+RNOx_lgt(i,j)*height(ih)*
     &         (ih-alttop)
               cycle ih_loop
             endif
             L=L+1
           enddo
         enddo ih_loop

! Save tracer 3D source. Convert from gN/min to kgN/s:
         do L=1,levtrop
           tr3Dsource(i,j,L,nOther,n_NOx) = 
     &     srclight(L)*pmin2psec*1.d-3
         enddo 
         do L=levtrop+1,LM
           tr3Dsource(i,j,levtrop,nOther,n_NOx) = 
     &     tr3Dsource(i,j,levtrop,nOther,n_NOx) + 
     &     srclight(L)*pmin2psec*1.d-3
         enddo  
#ifdef ACCMIP_LIKE_DIAGS
         do L=1,LM
           taijls(i,j,L,ijlt_NOxLgt)=taijls(i,j,L,ijlt_NOxLgt) +
     &     tr3Dsource(i,j,L,nOther,n_NOx)*byaxyp(i,j)
         enddo
#endif
      enddo  ! I
      enddo  ! J
         
      end subroutine get_lightning_NOx
#endif
