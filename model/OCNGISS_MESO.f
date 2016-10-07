
#include "rundeck_opts.h"

      MODULE MESOML

      USE MODEL_COM,  only : nstep=>itime
      USE CONSTANT,   only : grav,omega
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : idm=>imo,jdm=>jmo,kdm=>lmo,dzo
      USE OFLUXES,    only : oRSI,oAPRESS
      USE OCEAN,      only : ZOE=>ZE,g0m,s0m,mo,dxypo,focean,lmm,cospo
     .                      ,oLON_DG,oLAT_DG,uo,vo,sinpo,im,dxpo,dypo
     .                      ,kpl
      USE OCEAN_DYN, Only : DH
      USE OCNMESO_COM, only: RHOX, RHOY
#ifdef OCN_GISS_MESO
      USE ODIAG, only: oij=>oij_loc,oijl=>oijl_loc
     .            ,ij_eke,ij_rd,ijl_ueddy,ijl_veddy,ijl_n2
      USE OCEAN,      only : auvel,avvel,kappam3d_sm
#endif
      USE ODIAG, only: zoc

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
     ., HALO_UPDATE, HALO_UPDATE_COLUMN, NORTH, SOUTH
      use TimerPackage_mod

      IMPLICIT NONE

      SAVE

      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ktap
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: kappam3d0,
c     REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: 
     *  OMEGAx3d,OMEGAy3d,Uplus,Vplus,Wplus
     * ,Unew,Vnew,Wnew,LX3d,LY3d,kappamsigmsx3d,kappamsigmsy3d
     * ,kappamsigx3d,kappamsigy3d

      contains

      SUBROUTINE ALLOC_MESOML
        implicit none
        allocate(ktap (idm,ogrid%j_strt_halo:ogrid%j_stop_halo) )
        allocate(OMEGAx3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(OMEGAy3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Uplus (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Vplus (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Wplus (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Unew (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Vnew (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(Wnew (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(LX3d (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(LY3d (idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(kappam3d0(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm))
        allocate(kappamsigmsx3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,
     *           kdm))
        allocate(kappamsigmsy3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,
     *           kdm))
        allocate(kappamsigx3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,
     *           kdm))
        allocate(kappamsigy3d(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,
     *           kdm))
      END SUBROUTINE ALLOC_MESOML

      END MODULE MESOML

      subroutine OCN_mesosc(kappam3d)
c     subroutine OCN_mesosc

cc    USE MODEL_COM,  only : nstep=>itime,itimei
cc   .                    ,JMON,jhour,nday,jdate,jday
cc   . ,iyear1,jdendofm,jyear,aMON,dtsrc
cc   . ,xlabel,lrunid
c     USE MODEL_COM,  only : nstep=>itime
c     USE CONSTANT,   only : grav,omega,sday
c     USE OCEANR_DIM, only : ogrid
c     USE OCEANRES,   only : idm=>imo,jdm=>jmo,kdm=>lmo,dzo
c     USE OFLUXES,    only : oRSI,oAPRESS
c     USE OCEAN,      only : ZOE=>ZE,g0m,s0m,mo,dxypo,focean,lmm,cospo
c    .                      ,oLON_DG,oLAT_DG,uo,vo,sinpo,im,dxpo,dypo
c     USE KPP_COM,    only : kpl
c     USE OCEAN_DYN, Only : DH
c     USE GM_COM, only: RHOX, RHOY
c#ifdef OCN_GISS_MESO
c     USE ODIAG, only: oij=>oij_loc,oijl=>oijl_loc
c    .            ,ij_eke,ij_rd,ijl_ueddy,ijl_veddy,ijl_n2
c     USE OCEAN,      only : auvel,avvel,kappam3d_sm
c#endif
c     USE ODIAG, only: zoc

c     USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
c    ., HALO_UPDATE, NORTH, SOUTH
c     use TimerPackage_mod

      USE MESOML


      implicit none

      integer, parameter :: itest=97,jtest=128
c     integer, parameter :: it=16,jt=66
      integer it,jt
      integer i,j,k,l,ndepi
      integer i_0,i_1,j_0,j_1

      Real*8,External   :: VOLGSP,TEMGSP,TEMGS
      Real*8,External   :: VOLGS
      real*8  g,s,pres
      real*8   p1d(kdm+1),dp1d(kdm),temp1d(kdm),saln1d(kdm)
      real*8   rho_water,amld_cgs,Rd
      real*8   K0,Ustar(kdm),Vstar(kdm)
      real*8   z_cm(kdm+1),dens_cgs(kdm)
     .        ,uvel_cgs(kdm),vvel_cgs(kdm)
     .        ,drhodz_cgs(kdm),coriol,n2(kdm)
     .        ,drhodx_cgs(kdm),drhody_cgs(kdm)
      real*8  drhopdz_cgs(kdm),densulave_cgs(kdm),n20(kdm)
      REAL*8,ALLOCATABLE :: presi(:),densu(:),densl(:)
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     .   kappam3d,zma,zh,zt
c     REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
c    .   OMEGAx3d,OMEGAy3d,dzOMEGAx3d,dzOMEGAy3d,
c    .   dxOMEGAx3d,dyOMEGAy3d,Uplus,Vplus,Wplus 
c    .  ,temp3d,dTdx,dTdy,dTdz,Fdiffx,Fdiffy
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     .   dzOMEGAx3d,dzOMEGAy3d,
     .   dxOMEGAx3d,dyOMEGAy3d,
     .   temp3d,dTdx,dTdy,dTdz,Fdiffx,Fdiffy,
     .   OMEGAxT3d,OMEGAyT3d,OMEGAx1mT3d,OMEGAy1mT3d
c    .  ,dxOMEGAy3di,dyOMEGAy3di,dxOMEGAy3de,dyOMEGAy3de
     .  ,dyOMEGAx3d,dxOMEGAy3d
     .  ,dxkappamsigx3d,dykappamsigx3d,
     .   dxkappamsigy3d,dykappamsigy3d,dzkappamsigx3d,dzkappamsigy3d
c     INTEGER,DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo) :: ktap
      REAL*8 wta,wta1m,wta3m,beta
      integer ip1,im1

      logical vrbos

      INTEGER, SAVE :: IFIRST = 1

      if ( IFIRST.ne.0 ) then
        IFIRST = 0
        call ALLOC_MESOML
      endif

      it=42; jt=127


      kappam3d=0. 
      kappam3d0=0.

      ktap=0
      OMEGAx3d=0.
      OMEGAy3d=0.
      dzOMEGAx3d=0.
      dzOMEGAy3d=0.
      dxOMEGAx3d=0.
      dyOMEGAy3d=0.
      Uplus=0.
      Vplus=0.
      Wplus=0.
      temp3d=0.
      n20=0.
      OMEGAxT3d=0.; OMEGAyT3d=0.; OMEGAx1mT3d=0.; OMEGAy1mT3d=0.

      kappamsigx3d=0.; kappamsigy3d=0.
      dxkappamsigx3d=0.; dykappamsigx3d=0.
      dxkappamsigy3d=0.; dykappamsigy3d=0.
      dzkappamsigx3d=0.; dzkappamsigy3d=0.
      Unew=0.; Vnew=0.; Wnew=0.
      LX3d=0.; LY3d=0.
      kappamsigmsx3d=0.; kappamsigmsy3d=0.

      if (nstep.eq.0) return


      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP

C-- z at tracer level
      DO j=j_0,j_1
      DO i=i_0,i_1
        zh(i,j,1)=DH(i,j,1)/2.
        zma(i,j,1)=DH(i,j,1)
        DO k=2,kdm-1
          zh(i,j,k)=zma(i,j,k-1)+DH(i,j,k)/2.
          zma(i,j,k)=zma(i,j,k-1)+DH(i,j,k)
        ENDDO
      ENDDO
      ENDDO
      zh(:,:,kdm)=zoc(kdm)

      zt=-zh*100.
C--

      CALL HALO_UPDATE_COLUMN(ogrid,rhox)
      CALL HALO_UPDATE_COLUMN(ogrid,rhoy)
      CALL HALO_UPDATE(ogrid,
     *                 vo(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)

      do 1000 j=j_0,j_1
      do 1000 i=i_0,i_1
      IF(FOCEAN(I,J).gt.0.) THEN

      im1=i-1
      IF(i.eq.1) im1=idm

      vrbos=.false.
c     vrbos=.true.
      if (i.eq.itest.and.j.eq.jtest) vrbos=.false.

! max depth cell
      ndepi=lmm(i,j)

C-- TONY - 03/07/11
C-- Calculate an interface pressure array
      ALLOCATE(presi(0:lmm(i,j)),densu(lmm(i,j)),densl(lmm(i,j)))
      presi=0.; densu=0.; densl=0.
      presi(0) = 0.
      DO k=1,lmm(i,j)
        presi(k) = presi(k-1) + MO(i,j,k)*GRAV
      ENDDO
C-- Calculate a pot. density referenced to the upper and lower interfaces
      DO k=1,lmm(i,j)
        g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        densu(k) = 1d0/VOLGSP(g,s,presi(k-1)) 
        densl(k) = 1d0/VOLGSP(g,s,presi(k))
c       densulave_cgs(k) = 0.5*(densu(k)+densl(k))*1.d-3
c       densu(k) =1d0/volgs(g,s)
      ENDDO
C--

!change units
       pres = oAPRESS(i,j)    !surface atm. pressure
       do k=1,lmm(i,j)
         pres=pres+MO(I,J,k)*GRAV*.5
         g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         temp1d(k)=TEMGSP(g,s,pres)    !in situ   temperature
         temp3d(I,J,k)=TEMGSP(g,s,pres)
         saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
         rho_water = 1d0/VOLGSP(g,s,pres)
         dp1d(k)=MO(I,J,K)/rho_water   !local thickenss of each layer in meters
         !add missing part of density to get to the bottom of the layer
         !now pres is at the bottom of the layer
         pres=pres+MO(I,J,k)*GRAV*.5
      enddo

C-- convert in CGS
      p1d(1)=0.
      do k=2,kdm+1
      p1d(k)=p1d(k-1)+dp1d(k-1)
      enddo
      do k=1,kdm
       z_cm(k) = p1d(k)*100.d0                !depth in cm
       dens_cgs(k) = MO(i,j,k)/max(1.d0,dp1d(k)) *1.d-3 
       uvel_cgs(k) = uo(i,j,k) * 100.d0
       vvel_cgs(k) = vo(i,j,k) * 100.d0
       IF(uo(i,j,k).ne.0..AND.uo(im1,j,k).ne.0.)
     *   uvel_cgs(k) = 0.5*(uo(i,j,k)+uo(im1,j,k)) * 100.d0
       IF(vo(i,j,k).ne.0..AND.vo(i,j-1,k).ne.0.)
     *   vvel_cgs(k) = 0.5*(vo(i,j,k)+vo(i,j-1,k)) * 100.d0
      enddo

c     IF(j.eq.90) WRITE(*,*)'uvel:',i,uvel_cgs(1)

      z_cm(kdm+1) = p1d(kdm+1)*100.d0                !depth in cm
      coriol = 2d0*omega*sinpo(j)      !evaluated at tracer points
      beta = 2.*omega*cospo(j)/6371000.
      beta = beta/100. !in cgs
!drhodz_cgs
      do k=1,kdm-1
         if (dens_cgs(k+1).ne.0.d0.and.dens_cgs(k).ne.0.d0) then
           drhodz_cgs(k) = (dens_cgs(k+1) - dens_cgs(k))          ! at (i,j,k+1/2)   interface
c    .                 /(max(1.d0,z_cm(k+1)+z_cm(k)))/2.
     .                 /(max(1.d0,z_cm(k+1)-z_cm(k)))
         elseif (dens_cgs(k+1).eq.0.d0.
     *       and.dens_cgs(k).ne.0.d0.
     *       and.dens_cgs(max(1,k-1)).ne.0.d0) then
           drhodz_cgs(k)=drhodz_cgs(max(1,k-1))
         else
           drhodz_cgs(k) = 1.d0
         endif
      enddo
      do k=1,kdm
         drhodx_cgs(k)=rhox(k,i,j)*1.d-3/100.d0
         drhody_cgs(k)=rhoy(k,i,j)*1.d-3/100.d0
         IF(rhox(k,i,j).ne.0..AND.rhox(k,im1,j).ne.0.)
     *     drhodx_cgs(k)=0.5*(rhox(k,i,j)+rhox(k,im1,j))*1.d-3/100.d0
         IF(rhoy(k,i,j).ne.0..AND.rhoy(k,i,j-1).ne.0.)
     *     drhody_cgs(k)=0.5*(rhoy(k,i,j)+rhoy(k,i,j-1))*1.d-3/100.d0
      enddo
C-- TONY - 03/07/11
C-- Calculate a new drhodz using the potential densities
      drhopdz_cgs = 1.d30
      DO k=1,lmm(i,j)-1
        drhopdz_cgs(k)=(1.d-3)*(densu(k+1)-densl(k))/
     *                 ((-zt(i,j,k+1))-(-zt(i,j,k)))
      ENDDO
      drhopdz_cgs(lmm(i,j)) = drhopdz_cgs(max(1,lmm(i,j)-1))
C
! compute Brunt-Vaisala
      do k=1,kdm-1
c       if (drhodz_cgs(k).ne.1d30) then
c       n2(k) = GRAV* drhodz_cgs(k)*100.d0     !no minus sign here, because drhodz defined as rho(k+1)-rho(k)
c    .            /(max(1.d0,dens_cgs(k+1)+dens_cgs(k))/2.)   !at (i,j,k+1/2)   interface
C-- TONY - 03/07/11 - Updated n2 using the pot density instead of the in situ one
        if (drhopdz_cgs(k+1).ne.1d30) then
        n2(k) = GRAV*100.d0 * drhopdz_cgs(k)
     *         /(max(1.d0,(1.d-3)*(densu(k+1)+densl(k)))/2.)
        else
        n2(k) = 1.d30
        endif
        IF(k.eq.lmm(i,j)) THEN
        IF(lmm(i,j).eq.2) THEN
          n2(k)=n2(k-1)
        ELSE
          n2(k)=n2(k-1)+((n2(k-1)-n2(k-2))/(z_cm(k-1)-z_cm(k-2)))
     *          *(z_cm(k)-z_cm(k-1))
        ENDIF
        ENDIF
      enddo
      DEALLOCATE(presi,densu,densl)
      n20=n2
! mixed layer depth
      amld_cgs = 0.d0
      do k=1,kpl(i,j)
      if (dp1d(k).ne.1.d30) then
          amld_cgs = amld_cgs + dp1d(k)*100.d0      ! at (i,j,k) the middle of last layer in MLD
      endif
      enddo
      if (amld_cgs.gt.p1d(ndepi)*100.d0) amld_cgs=p1d(ndepi)*100.d0

#ifdef OCN_GISS_MESO
C-- time-averaging of u,v
      wta=exp(-1./240.)
      auvel(i,j,:) = wta*auvel(i,j,:) + (1.-wta)*uvel_cgs(:)
      avvel(i,j,:) = wta*avvel(i,j,:) + (1.-wta)*vvel_cgs(:)

      uvel_cgs=auvel(i,j,:)
      vvel_cgs=avvel(i,j,:)
#endif

      call mesoscales1d(kdm,ndepi,
     .      dens_cgs,uvel_cgs,vvel_cgs,
     .      n2,drhodx_cgs,drhody_cgs,drhopdz_cgs,coriol,amld_cgs
     .     ,Rd,K0,Ustar,Vstar,i,j,kappam3d(i,j,:),beta,
     .      OMEGAxT3d(i,j,:),OMEGAyT3d(i,j,:),
     .      OMEGAx1mT3d(i,j,:),OMEGAy1mT3d(i,j,:),ktap(i,j),
     .      kappamsigx3d(i,j,:),kappamsigy3d(i,j,:),
     .      LX3d(i,j,:),LY3d(i,j,:),
     .      kappamsigmsx3d(i,j,:),kappamsigmsy3d(i,j,:))

c     IF(i.eq.257.AND.j.eq.91) THEN
c       WRITE(*,*)'DEBUG-OMEGA:',i,j,ktap(i,j)
c       DO k=1,kdm
c         WRITE(*,*)k,OMEGAx3d(i,j,k),OMEGAy3d(i,j,k)
c       ENDDO
c     ENDIF

C--   Convert kappam3d from CGS to MKS
      kappam3d(i,j,:)=kappam3d(i,j,:)*1.d-4

#ifdef OCN_GISS_MESO
      wta1m=exp(-1./1440.)
      wta3m=exp(-1./(3.*1440.))

c     IF(i.eq.42.AND.j.eq.127) WRITE(*,*)'DEBUG-k3dsm before:',
c    *  nstep,i,j,kappam3d(i,j,1),kappam3d_sm(i,j,1)

      kappam3d_sm(i,j,:)=wta3m*kappam3d_sm(i,j,:)
     *                  +(1.-wta3m)*kappam3d(i,j,:)
      kappam3d(i,j,:)=kappam3d_sm(i,j,:)

c     IF(i.eq.42.AND.j.eq.127) WRITE(*,*)'DEBUG-k3dsm after:',
c    *  nstep,i,j,kappam3d(i,j,1),kappam3d_sm(i,j,1)
#endif


#ifdef OCN_GISS_MESO
      if (vrbos) then
      write(*,'(a,2i6)')'MESOSCALES:',nstep
      endif


       OIJ(I,J,IJ_eke)  = OIJ(I,J,IJ_eke) + K0      ! eddy kinetic energy, ocean
       OIJ(I,J,IJ_rd )  = OIJ(I,J,IJ_rd ) + Rd      ! Rossby radius of deformation
C
       DO k=1,kdm
          OIJL(I,J,k,IJL_n2   )= OIJL(I,J,k,IJL_n2   ) + n20(k)    ! brunt vaisala squared
          OIJL(I,J,k,IJL_ueddy)= OIJL(I,J,k,IJL_ueddy) + ustar(k) ! ustar, eddy induced velocity (Canuto)
          OIJL(I,J,k,IJL_veddy)= OIJL(I,J,k,IJL_veddy) + vstar(k) ! vstar, eddy induced velocity (Canuto)
       ENDDO
#endif

      endif    !focean

 1000 continue


c      DO j=j_0,j_1
c     IF(j.eq.127) WRITE(*,*)'DEBUG-k3dsm final:',
c    *  nstep,42,j,kappam3d(42,j,1),i_0,i_1,j_0,j_1
c      ENDDO


C--   Convert kappam3d from CGS to MKS
c     kappam3d=kappam3d*1.d-4
C--   Maximum of 15000 m2/s
C--   Minimum value of 1 m2/s
      DO k=1,kdm
       DO j=j_0,j_1
         DO i=i_0,i_1
         IF(FOCEAN(I,J).gt.0.) THEN
c     IF(i.eq.42.AND.j.eq.127) WRITE(*,*)'DEBUG-k3dsm final:',
c    *  nstep,i,j,kappam3d(i,j,k)
           kappam3d(i,j,k)=min(kappam3d(i,j,k),15000.d0)
           kappam3d(i,j,k)=max(kappam3d(i,j,k),1.d0)
         ENDIF
         ENDDO
       ENDDO
      ENDDO
      kappam3d0=kappam3d

c     DO j=ogrid%J_STRT,ogrid%J_STOP
c     IF(j.eq.127) THEN
c       WRITE(*,*)'DEBUG-kappam3d-MESO:',kappam3d0(42,j,1),
c    *    kappam3d(42,j,1),dxpo(j),dypo(j)
c     ENDIF
c     ENDDO
C-- Convert OMEGA from cgs to MKS
      OMEGAx3d = OMEGAxT3d + OMEGAx1mT3d*(kappam3d0*1.d4)
      OMEGAy3d = OMEGAyT3d + OMEGAy1mT3d*(kappam3d0*1.d4)
      OMEGAx3d=1.d-4*OMEGAx3d
      OMEGAy3d=1.d-4*OMEGAy3d



C--
C-- Calculate the new eddy induced velocity Uplus
C--- Calculate dzOMEGAx,dzOMEGAy
       DO j=j_0,j_1
         DO i=i_0,i_1
         IF(FOCEAN(I,J).gt.0..AND.ktap(i,j).gt.1) THEN
C?? z_cm??
           CALL d1sym(ktap(i,j),zt(i,j,1:ktap(i,j))/100.,
     *                OMEGAx3d(i,j,1:ktap(i,j)),
     *                dzOMEGAx3d(i,j,1:ktap(i,j)),1)
           CALL d1sym(ktap(i,j),zt(i,j,1:ktap(i,j))/100.,
     *                OMEGAy3d(i,j,1:ktap(i,j)),
     *                dzOMEGAy3d(i,j,1:ktap(i,j)),2)
cc       IF(i.eq.42.AND.j.eq.127)  THEN
c        IF(i.eq.it.AND.j.eq.jt)  THEN
cc       IF(i.eq.it.AND.j.eq.jt.OR.
cc   *      i.eq.it+1.AND.j.eq.jt.OR.
cc   *      i.eq.it.AND.j.eq.jt+1) THEN
c          WRITE(*,*)'DEBUG-dzOMEGA:',nstep,i,j,ktap(i,j)
c          DO k=1,ktap(i,j)
c            WRITE(*,*)k,OMEGAx3d(i,j,k),dzOMEGAx3d(i,j,k),
c    *         OMEGAxT3d(i,j,k),
c    *         OMEGAx1mT3d(i,j,k)*(kappam3d0(i,j,k)*1.d4),
c    *         OMEGAx1mT3d(i,j,k),kappam3d0(i,j,k),
c    *         OMEGAyT3d(i,j,k),
c    *         OMEGAy1mT3d(i,j,k)*(kappam3d0(i,j,k)*1.d4),
c    *         OMEGAy1mT3d(i,j,k)
c          ENDDO
c        ENDIF
         ENDIF
         ENDDO
       ENDDO
C--- Calculate dxOMEGAx,dyOMEGAy
      CALL HALO_UPDATE(ogrid,
     *                OMEGAy3d(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *              dyOMEGAy3d(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)

c     do 102 j=j_0,j_1
c     do 102 i=i_0,i_1
c       ip1=i+1
c       im1=i-1
c       IF(i.eq.1) im1=288
c       IF(i.eq.288) ip1=1
c       do k=1,ktap(i,j)
C-- dxOMEGAx
c       IF (OMEGAx3d(i,j,k).ne.0.) THEN
c         IF(OMEGAx3d(ip1,j,k).ne.0.) THEN
c          dxOMEGAx3d(i,j,k)=(OMEGAx3d(ip1,j,k)-OMEGAx3d(i,j,k))/dxpo(j)
c         ELSEIF(OMEGAx3d(ip1,j,k).eq.0..AND.OMEGAx3d(im1,j,k).ne.0.
c    *           .AND.i.ne.288) THEN
c          dxOMEGAx3d(i,j,k)=dxOMEGAx3d(im1,j,k)
c         ENDIF
c       ENDIF
C-- dyOMEGAy
c       IF (OMEGAy3d(i,j,k).ne.0.) THEN
c         IF(OMEGAy3d(i,j+1,k).ne.0.) THEN
c          dyOMEGAy3d(i,j,k)=(OMEGAy3d(i,j+1,k)-OMEGAy3d(i,j,k))/dypo(j)
c         ELSEIF(OMEGAy3d(i,j+1,k).eq.0..
c    *           AND.OMEGAy3d(i,j-1,k).ne.0.) THEN
c          dyOMEGAy3d(i,j,k)=dyOMEGAy3d(i,j-1,k)
c         ENDIF
c       ENDIF
c       enddo
c102  CONTINUE
c     dxOMEGAx3d=dxOMEGAx3d/100.
c     dyOMEGAy3d=dyOMEGAy3d/100.
C??? dxpo,dypo n cm???

      CALL get_gradients(mo,OMEGAx3d,0,dxOMEGAx3d,dyOMEGAx3d)
      CALL get_gradients(mo,OMEGAy3d,0,dxOMEGAy3d,dyOMEGAy3d)
c     CALL get_gradients(mo,OMEGAy3d,1,dxOMEGAy3de,dyOMEGAy3de)

c     DO j=j_0,j_1
c       DO i=i_0,i_1
c       IF(i.eq.42.AND.j.eq.127)  THEN
c         WRITE(*,*)'DEBUG-GRADs:',nstep,i,j,ktap(i,j),lmm(i,j)
c         DO k=1,ktap(i,j)
c           WRITE(*,123)k,OMEGAy3d(i,j,k),OMEGAy3d(i,j+1,k),
c    *        dyOMEGAy3d(i,j,k),dyOMEGAy3di(i,j,k),
c    *        dyOMEGAy3de(i,j,k)
c         ENDDO
c       ENDIF
c       ENDDO
c     ENDDO
c123  FORMAT(I,3E)

      Uplus = -dzOMEGAx3d
      Vplus = -dzOMEGAy3d
      Wplus = dxOMEGAx3d + dyOMEGAy3d
c     DO j=j_0,j_1
c       DO i=i_0,i_1
c       IF(i.eq.it.AND.j.eq.jt)  THEN
c         WRITE(*,*)'DEBUG-Wplus:',nstep,i,j,ktap(i,j)
c         DO k=1,ktap(i,j)
c           WRITE(*,*)k,Wplus(i,j,k),dxOMEGAx3d(i,j,k),dyOMEGAy3d(i,j,k)
c         ENDDO
c       ENDIF
c       ENDDO
c     ENDDO
c     it=16; jt=66+1
c     DO j=j_0,j_1
c       DO i=i_0,i_1
c       IF(i.eq.it.AND.j.eq.jt)  THEN
c         WRITE(*,*)'DEBUG-Wplus:',nstep,i,j,ktap(i,j)
c         DO k=1,ktap(i,j)
c           WRITE(*,*)k,Wplus(i,j,k),dxOMEGAx3d(i,j,k),dyOMEGAy3d(i,j,k)
c         ENDDO
c       ENDIF
c       ENDDO
c     ENDDO

C-- Calculate Fdiff
C-- Calculate the temperature gradients
c     dTdx=0.
c     dTdy=0.
c     CALL HALO_UPDATE(ogrid,
c    *                 temp3d(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                 FROM=SOUTH+NORTH)
cc    WRITE(*,*)'DEBUG-BORNES:',ogrid%gid,i_0,i_1,j_0,j_1
c     do 103 j=j_0,j_1
c     do 103 i=i_0,i_1
c       ip1=i+1
c       im1=i-1
c       IF(i.eq.1) im1=288
c       IF(i.eq.288) ip1=1
c       do k=1,ktap(i,j)
c       IF (temp3d(i,j,k).ne.0.) THEN
C-- dTdx
c         IF(temp3d(ip1,j,k).ne.0.) THEN
c           dTdx(i,j,k)=(temp3d(ip1,j,k)-temp3d(i,j,k))/dxpo(j)
c         ELSEIF(temp3d(ip1,j,k).eq.0..AND.temp3d(im1,j,k).ne.0.
c    *           .AND.i.ne.288) THEN
c           dTdx(i,j,k)=dTdx(im1,j,k)
c         ENDIF
C-- dTdy
c         IF(temp3d(i,j+1,k).ne.0.) THEN
c           dTdy(i,j,k)=(temp3d(i,j+1,k)-temp3d(i,j,k))/dypo(j)
c         ELSEIF(temp3d(i,j+1,k).eq.0..AND.temp3d(i,j-1,k).ne.0.) THEN
c           dTdy(i,j,k)=dTdy(i,j-1,k)
c         ENDIF
c       ENDIF
c       enddo
c103  CONTINUE
c     dTdx=dTdx/100.
c     dTdy=dTdy/100.
C--- Calculate dzT
c     dTdz=0.
c     DO j=j_0,j_1
c       DO i=i_0,i_1
c         CALL d1sym(ktap(i,j),zt(i,j,1:ktap(i,j)),
c    *               temp3d(i,j,1:ktap(i,j)),
c    *               dTdz(i,j,1:ktap(i,j)))
c       ENDDO
c     ENDDO
C-- Fdiff = -(kappam*grad(tau) + OMEGA*dtaudz)
c     Fdiffx = 0.
c     Fdiffy = 0.
c     Fdiffx = -(kappam3d*dTdx + OMEGAx3d*dTdz)
c     Fdiffy = -(kappam3d*dTdy + OMEGAy3d*dTdz)

c     do 104 j=j_0,j_1
c     do 104 i=i_0,i_1

c     IF(i.eq.42.AND.j.eq.127) THEN
c       WRITE(*,*)'DEBUG-Fdiff:',i,j,ktap(i,j)
c       DO k=1,ktap(i,j)
c         WRITE(*,*)k,Fdiffx(i,j,k),Fdiffy(i,j,k)
c       ENDDO
c     ENDIF

c104  CONTINUE

C--
C-- A-Regime:Calculate the new eddy induced velocity Unew
C--- Calculate dzkappamsigx,dzkappamsigy
C-- first: convert to MKS:
       kappamsigx3d=1.d-4*kappamsigx3d
       kappamsigy3d=1.d-4*kappamsigy3d
       kappamsigmsx3d=1.d-4*kappamsigmsx3d
       kappamsigmsy3d=1.d-4*kappamsigmsy3d
C--
       dzkappamsigx3d=0.; dzkappamsigy3d=0.
       DO j=j_0,j_1
         DO i=i_0,i_1
         IF(FOCEAN(I,J).gt.0..AND.lmm(i,j).gt.1) THEN
           CALL d1sym(lmm(i,j),zt(i,j,1:lmm(i,j))/100.,
     *      kappamsigx3d(i,j,1:lmm(i,j)),dzkappamsigx3d(i,j,1:lmm(i,j))
     &      ,3)
           CALL d1sym(lmm(i,j),zt(i,j,1:lmm(i,j))/100.,
     *      kappamsigy3d(i,j,1:lmm(i,j)),dzkappamsigy3d(i,j,1:lmm(i,j))
     &      ,4)
         ENDIF
         ENDDO
       ENDDO

C--- Calculate dxkappamsigx3d,dykappamsigy3d
      CALL HALO_UPDATE(ogrid,
     *            kappamsigy3d(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *            FROM=SOUTH+NORTH)
      CALL get_gradients(mo,kappamsigx3d,0,
     *                   dxkappamsigx3d,dykappamsigx3d)
      CALL get_gradients(mo,kappamsigy3d,0,
     *                   dxkappamsigy3d,dykappamsigy3d)


      Unew = -dzkappamsigx3d
      Vnew = -dzkappamsigy3d
      Wnew = +(dxkappamsigx3d + dykappamsigy3d)

c     do 104 j=j_0,j_1
c     do 104 i=i_0,i_1
c     IF(i.eq.42.AND.j.eq.127) THEN
c       WRITE(*,*)'DEBUG-Unew-in meso',nstep,i,j
c       DO k=1,lmm(i,j)
c         WRITE(*,*)k,Unew(i,j,k),Vnew(i,j,k),Wnew(i,j,k)
c       ENDDO
c       WRITE(*,*)'DEBUG-Unew-in meso',nstep,i,j
c       DO k=1,lmm(i,j)
c         WRITE(*,*)k,kappamsigx3d(i,j,k),dxkappamsigx3d(i,j,k),
c    *     kappamsigy3d(i,j,k),dykappamsigy3d(i,j,k)
c       ENDDO
c     ENDIF
c104  CONTINUE
     


      RETURN
      end subroutine OCN_mesosc

      subroutine mesoscales1d(km,ndepi,
     .      rhoi,UI,VI,n2,dxrho,dyrho,dzrho,f,ml
     .     ,rdm,K02,Ustar,Vstar,
     .     igrid,jgrid,KAPPAMZa,beta,OMEGAxT,
     .     OMEGAyT,OMEGAx1mT,OMEGAy1mT,ktap,kappamsigx,kappamsigy,
     .     LX3d,LY3d,kappamsigmsx,kappamsigmsy)

        USE MODEL_COM,  only : nstep=>itime
        USE ODIAG, only : zoc
        USE OCEAN_DYN, Only : DH
        USE OCEAN,    only : kpl

        IMPLICIT NONE

        real*8, parameter :: a02=0.03
        INTEGER k,km,kmli,igrid,jgrid,ndepi
        REAL*8 z_cm(km),zma(km),zh(km)
        REAL*8 Ustar(km),Vstar(km)
        REAL*8 KAPPAM,KAPPAMZa(km)
        REAL*8 SIGMAT,f,ml,rd,rdm,pi,SIGMAS
        PARAMETER(SIGMAT=1.)
        INTEGER kkpmax2,kkpint,ktap,kkpint2,ktap2
        REAL*8 DZLXB1AVE,DZLYB1AVE,kpmax2
        REAL*8 UI(km),VI(km)
        REAL*8 UB1AVE,VB1AVE,UL,VL
        REAL*8 rhoi(km)
        REAL*8 frd,n2(km),lr,rd0
        REAL*8 K02,CK,K02A,K02D,CA,CD,INTGAMMAA,INTGAMMAD,INTGAMMA,
     *    INTGAMMA1,INTEXP,CD1,INTGAMMAS
        REAL*8 dxrho(km),dyrho(km),dzrho(km)
        REAL*8, ALLOCATABLE :: n2t(:),n2t1(:),z(:),zm(:),b1(:),b1t(:),
     *   U(:),V(:),rho(:),n2tr(:),
     *   LX(:),LY(:),DZLX(:),DZLY(:),
     *   UINT(:),VINT(:),UT(:),VT(:),
     *   dxb(:),dyb(:),KP(:),KINT(:),
     *   zint(:),UREV(:),VREV(:),KPINT(:),
     *   b1rev(:),DZLXREV(:),DZLYREV(:),
     *   UTINT(:),VTINT(:),DZU(:),DZV(:),
     *   F1X(:),F2X(:),F1Y(:),F2Y(:),KP2(:),KINT2(:),KP2INT(:),
     *   KAPPAMZ(:),GAMMA(:),n2sm(:),Tap(:),
     *   Gsx(:),Gsy(:),DZGAMMA(:),Gsxint(:),Gsyint(:),
     *   DZKAPPAMZ(:),Iintx(:),Iinty(:),Ix(:),Iy(:),alpha(:),
     *   Ixd(:),Iyd(:),alphad(:),ung(:),vng(:),alphax(:),alphay(:)
        INTEGER im,kmi
        REAL*8 Nm
        REAL*8 zw(km),zwt(km),zt(km)
        REAL*8 k02count,Tapcount,tappos,hstaroh,n20(km)
        REAL*8, ALLOCATABLE :: ud(:),vd(:),MK(:),GAMMAREV(:)
        REAL*8 a02b,b0,beta,Tapten(21),stap,ztap,Ro(km),Isx,Isy
        INTEGER it,jt
c       REAL*8 OMEGAx(km),OMEGAy(km)
        REAL*8 OMEGAxT(km),OMEGAyT(km),OMEGAx1mT(km),OMEGAy1mT(km)
        REAL*8 kappamsigx(km),kappamsigy(km),
     *   kappamsigmsx(km),kappamsigmsy(km),LX3d(km),LY3d(km)
        REAL*8, ALLOCATABLE :: KAPPAX(:),KAPPAY(:),FvM(:),
     *   KAPPAlx(:),KAPPAly(:),n2wkb(:)

      it=42; jt=127
 
      pi=4.*atan(1.)
      rd=0.
      Ro=0.
      rd0=0.
      kkpint2=1
      ktap=1

      zh(1)=DH(igrid,jgrid,1)/2.
      zma(1)=DH(igrid,jgrid,1)
      DO k=2,km-1
        zh(k)=zma(k-1)+DH(igrid,jgrid,k)/2.
        zma(k)=zma(k-1)+DH(igrid,jgrid,k)
      ENDDO
      zh(km)=zoc(km)

      z_cm=zh*100.

      if (ndepi.eq.0) then
        K02=0.
        KAPPAMZa=0.
        goto 10
      endif

      do k=1,km
         zw(k) = z_cm(k)
         zt(k) =-z_cm(k)
      enddo
      zwt(1)=0.
      do k=2,km
        zwt(k) = zw(k-1)
      enddo

      kmi = ndepi

      kmli=km+1   !outside bounds
c     DO k=1,kmi-1
c       IF(zw(k).le.ml.and.ml.lt.zw(k+1)) kmli=k+1
c     ENDDO
C
c     kmli=min(kmli,kmi)
      kmli=kpl(igrid,jgrid)
C
      ALLOCATE(n2t(kmi),n2t1(kmi),z(kmi),zm(kmi),b1(kmi),b1t(kmi),
     * U(kmi),V(kmi),rho(kmi),n2tr(kmi),
     * LX(kmi),LY(kmi),DZLX(kmi),DZLY(kmi),
     * UINT(kmi),VINT(kmi),UT(kmi),VT(kmi),
     * dxb(kmi),dyb(kmi),KP(kmi),
     * zint(kmi),UREV(kmi),VREV(kmi),KPINT(kmi),
     * b1rev(kmi),DZLXREV(kmi),DZLYREV(kmi),
     * UTINT(kmi),VTINT(kmi),DZU(kmi),DZV(kmi),
     * F1X(kmi),F1Y(kmi),F2X(kmi),F2Y(kmi),KP2(kmi),KP2INT(kmi))

      ALLOCATE(Tap(kmi),GAMMA(kmi),n2sm(kmi))
      Tap=0.; GAMMA=0.; n2sm=0.
      z=0.
C-- TONY - 01/28/2013 - Test on n2
      IF(maxval(n2(1:min(kmli+1,kmi))).lt.0..OR.kmli.EQ.0) THEN
c       IF(jgrid.ne.180) THEN
c       WRITE(*,*)'DEBUG-N2!!!',igrid,jgrid,kmli,kmi
c       DO k=1,kmi
c         WRITE(*,*)k,n2(k)
c       ENDDO
c       ENDIF
        K02=0.
        GAMMA=1.
        ktap=1
        rdm=0.
        GO TO 101
      ENDIF
C--
      n20=n2
C-- M1: Calculate B1, frd and rd
      n2t(1)=n2(1)
C     n2t(1)=max(n2t(1),1.d-8)
      DO k=2,kmi
        n2t(k)=n2(k-1)
C       n2t(k)=max(n2t(k),1.d-8)
      ENDDO
      DO k=1,kmi-1
        n2tr(k)=0.5*(n2(k)+n2(k+1))
      ENDDO
      n2tr(kmi)=n2(kmi-1)
C
      z=-zw(1:kmi)
      zm=-z
C     n2t=((7.d-3)**2)*exp(2.*z/100000.)
c     n2t1=n2t
      n2t1=n2(:)
C
c     CALL baroclin1st(zm,n2t1,kmi,b1,im,Nm,frd)
c     CALL baroclin1st(zm,n2(:),kmi,b1,im,Nm,frd,igrid,jgrid)
      ALLOCATE(n2wkb(kmi))
      n2wkb=n2(1:kmi)
      CALL WKB(zm,n2wkb,kmi,frd,b1)
      DEALLOCATE(n2wkb)
      rd=dabs(frd/f)
      rd0=rd

C-- Tony - 01/28/2013
      IF(rd.EQ.0.) THEN
        K02=0.
        GAMMA=1.
        ktap=1
        rdm=0.
        GO TO 101
      ENDIF
C--

c     rdm=rd
      DO k=1,kmi-1
        b1t(k)=0.5*(b1(k)+b1(k+1))
      ENDDO
      b1t(kmi)=b1(kmi)
      b1=b1t
      U=UI(1:kmi)
      V=VI(1:kmi)
      rho=rhoi(1:kmi)
c
      DZU=0.; DZV=0.
      IF(kmi.gt.1) CALL d1sym(kmi,z,U,DZU,5)
      IF(kmi.gt.1) CALL d1sym(kmi,z,V,DZV,6)

c
c     ALLOCATE(n2sm(kmi))
      n2sm = n2(1:kmi)
      DO k=1,50
      CALL z121_mesosc(n2sm,kmi,kmi)
      ENDDO
      do k=1,kmi
c        if (n2tr(k).ne.0.d0)then 
c               LX(k)=-(f/n2tr(k))*DZV(k)
c               LY(k)=+(f/n2tr(k))*DZU(k)
         if (n2sm(k).ne.0.d0)then
                LX(k)=-(f/n2sm(k))*DZV(k)
                LY(k)=+(f/n2sm(k))*DZU(k)
c               LX(k)=-dxrho(k)/(-n2tr(k)*rho(k)/981.)
c               LY(k)=-dyrho(k)/(-n2tr(k)*rho(k)/981.)
         else
                LX(k)=0.d0
                LY(k)=0.d0
         endif
C-- TONY 03/08/11 -- limiter for isopycnal slopes
c        LX(k)=min(1.d-3,LX(k))
c        LY(k)=min(1.d-3,LY(k))
      enddo

c     IF(igrid.eq.272.AND.jgrid.EQ.45) THEN
c       WRITE(*,*)'DEBUG-LX in meso:',igrid,jgrid,kmi
c       DO k=1,kmi
c         WRITE(*,*)k,LX(k)
c       ENDDO
c     ENDIF
C
      UREV=0.
      VREV=0.
      DO k=1,kmi
        zint(k)=zt(kmi-k+1)
        UREV(k)=U(kmi-k+1)
        VREV(k)=V(kmi-k+1)
        b1rev(k)=b1(kmi-k+1)
      ENDDO
C
C-- Calculate GAMMA
c     ALLOCATE(GAMMA(kmi))
c     GAMMA=0.
c     GAMMA=(a02+(dabs(b1))**2)/(1.+a02)
C-- New Gamma
      DO k=kmi,1,-1
        IF(B1(k).ne.0.) THEN
          a02b=dabs(B1(k))
          GO TO 112
        ENDIF
      ENDDO
 112  CONTINUE
c     b0=dsqrt(a02b)/1.5
      b0=dsqrt(a02b)/dsqrt(2.d0)
      GAMMA=(a02b+(dabs(b1))**2+2.*B1*b0)/(1.+a02b+2.*b0)
      DO k=1,kmi
        IF(GAMMA(k).lt.0.) THEN
          WRITE(*,*)'DEBUG-GAMMA<0 i,j,k,Gamma:',igrid,jgrid,k,GAMMA(k)
          CALL STOP_MODEL ('GAMMA in meso',255)
        ENDIF
      ENDDO


      z=zt(1:kmi)
c     CALL B1AVERAGE(UREV,zint,b1rev,kmi,UB1AVE,kmli+0)
c     CALL B1AVERAGE(VREV,zint,b1rev,kmi,VB1AVE,kmli+0)
      CALL B1AVERAGEG(UREV,zint,b1rev,kmi,UB1AVE,kmli+0,GAMMA)
      CALL B1AVERAGEG(VREV,zint,b1rev,kmi,VB1AVE,kmli+0,GAMMA)
C
      DZLX=0; DZLY=0.
      IF(kmi.gt.1) CALL d1sym(kmi,z,LX,DZLX,7)
      IF(kmi.gt.1) CALL d1sym(kmi,z,LY,DZLY,8)

      DO k=1,kmi
        DZLXREV(k)=DZLX(kmi-k+1)
        DZLYREV(k)=DZLY(kmi-k+1)
      ENDDO

c     CALL B1AVERAGE(DZLXREV,zint,b1rev,kmi,DZLXB1AVE,kmli+0)
c     CALL B1AVERAGE(DZLYREV,zint,b1rev,kmi,DZLYB1AVE,kmli+0)
      CALL B1AVERAGEG(DZLXREV,zint,b1rev,kmi,DZLXB1AVE,kmli+0,GAMMA)
      CALL B1AVERAGEG(DZLYREV,zint,b1rev,kmi,DZLYB1AVE,kmli+0,GAMMA)
C--
C-- Calculate UT and VT
      DO k=1,kmi
        CALL dqtfg(zint(k:kmi),UREV(k:kmi),UINT(k:kmi),kmi-k+1)
        CALL dqtfg(zint(k:kmi),VREV(k:kmi),VINT(k:kmi),kmi-k+1)
        if (zint(k).ne.0.d0) then
          UTINT(k)=-(1./zint(k))*UINT(kmi)
          VTINT(k)=-(1./zint(k))*VINT(kmi)
        else
          UTINT(k)=0.d0
          VTINT(k)=0.d0
        endif
      ENDDO
      UTINT(kmi)=0.
      VTINT(kmi)=0.
      DO k=1,kmi
        UT(k)=UTINT(kmi-k+1)
        VT(k)=VTINT(kmi-k+1)
      ENDDO
c
C-- Calulate the Rhines Scale and replace rd by lr between +/-30 degrees
      lr=dsqrt(dsqrt(UT(kmli)**2+VT(kmli)**2)/(2.*2.1*1.d-13))
C-- Tony - 01/28/2013
      IF(kmi.EQ.1.OR.kmli.EQ.1) THEN
        K02=0.
        GAMMA=1.
        ktap=1
        rdm=0.
        GO TO 101
      ENDIF
      rd=min(rd,lr)
c     IF(jgrid.ge.70.AND.jgrid.le.110.AND.rd.eq.0.) rd=lr
c     IF(jgrid.ge.70.AND.jgrid.le.110) rd=lr
      rdm=rd
 
C-- Calculate F1
      if (rd.ne.0.d0) then                !Natassa
c        F1X = DZLXB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(-VB1AVE)
         F1X = DZLXB1AVE + (1.)*(1./f)*(1./rd**2)*(-VB1AVE)
      else
         F1X = 0.d0
      endif

      if (rd.ne.0.d0) then                !Natassa
c         F1Y = DZLYB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(+UB1AVE)
          F1Y = DZLYB1AVE + (1.)*(1./f)*(1./rd**2)*(+UB1AVE)
      else
          F1Y = 0.d0
      endif

C-- Calculate F2
      if (rd.ne.0.d0) then                !Natassa
c     F2X = (1./f)*(1./rd**2)*((1./SIGMAT)*(+VT)+V)
c     F2Y = (1./f)*(1./rd**2)*((1./SIGMAT)*(-UT)-U)
      F2X = 0.
      F2Y = 0.
      else
      F2X = 0.d0
      F2Y = 0.d0
      endif

C-- Calculate DXB and DYB
      dxb=-981.*dxrho(1:kmi)/rho(1:kmi)
      dyb=-981.*dyrho(1:kmi)/rho(1:kmi)

C-- Calculate Tapering function T(z)
      ktap=1
c     ALLOCATE(Tap(kmi))
c     Tap=0.
      DO k=1,kmi
        IF(((LX(k)**2+LY(k)**2)*n2tr(k)).ne.0.)
     *  Tap(k)=-z(k)*((F1X(k)+F2X(k))*dxb(k)+(F1Y(k)+F2Y(k))*dyb(k))
     *      /((LX(k)**2+LY(k)**2)*n2tr(k))
        IF(Tap(k).gt.1.) THEN
          ktap=k
          GO TO 116
        ENDIF
      ENDDO
 116  CONTINUE
c     IF(igrid.EQ.it.AND.jgrid.EQ.jt) THEN
cc       IF(igrid.eq.it.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it+1.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it.AND.jgrid.eq.jt+1) THEN
c       WRITE(*,*)'DEBUG-Tap:',nstep,igrid,jgrid,ktap,kmi
c       DO k=1,ktap
c         WRITE(*,118)k,Tap(k),z(k),(F1X(k)+F2X(k))*dxb(k),
c    *      (F1Y(k)+F2Y(k))*dyb(k),(LX(k)**2+LY(k)**2),
c    *      n2tr(k),n2(k),n20(k)
c       ENDDO
c     ENDIF
c118  FORMAT(I,8E)
C-- 2013-04-12: ktap before Tap=1??? or just for Uplus/Fdiff?
      ktap=max(ktap-1,1)
C--

c     Ustar=0.
c     Vstar=0.
c     Ustar(1:kmi)=Tap
c     DEALLOCATE(Tap)

C-- Calculate GAMMA
c     ALLOCATE(GAMMA(kmi))
c     GAMMA=(a02+(dabs(b1))**2)/(1.+a02)

C-- Calculate new K
      CK=27.
      KP2=sqrt(GAMMA)*z*((F1X+F2X)*dxb+(F1Y+F2Y)*dyb)
c     KP2=-z*((F1X+F2X)*sx+(F1Y+F2Y)*sy)*n2tr
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
c     kkpmax2=1
c     kpmax2=0.
c     DO k=1,kmi
c       IF(dabs(KP2(k)).gt.kpmax2) THEN
c         kpmax2=dabs(KP2(k))
c         kkpmax2=k
c       ENDIF
c     ENDDO
C
c     kkpint=max(1,kkpmax2)
c     kkpint=max(1,kmli)
c     kkpint=max(1,max(kmli,kkpmax2))

c     ktap2=max(ktap,kmli)
c     kkpint=max(1,ktap2)
      kkpint=max(1,ktap)
c     IF(ktap.le.kmli) kkpint=1

      K02D=0.
      ALLOCATE(KINT2(kkpint))
      CALL dqtfg(zint(kmi-(kkpint-1):kmi),KP2INT(kmi-(kkpint-1):kmi),
     *           KINT2,kkpint)
      if (rd.ne.0d0.and.z(kkpint).ne.0.d0) then
        K02D=-(rd**2)*KINT2(kkpint)
      else
        K02D=0.d0
      endif
      DEALLOCATE(KINT2)

      IF(K02D.lt.0.) K02D=0.
c     IF(ktap.eq.kmi.OR.kmli.eq.kmi) K02D=0.

      KAPPAMZa=0.

C-- TONY -  1/20/2012 - additional term
      K02A=0.
      KP2=0.
      KP2INT=0.
      KP2=sqrt(GAMMA)*(LX**2+LY**2)*n2tr
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      kkpint=max(1,ktap)
      kkpint2=max(1,kmi-kkpint)
      IF(ktap.eq.1) kkpint2=max(1,kmi-kmli)
      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      IF(ktap.ne.kmi.and.kmi.ne.1)
     *  K02A = (rd**2)*KINT2(kkpint2)
      DEALLOCATE(KINT2)

      IF(ktap.eq.kmi) K02A=0.
c     IF(K02D.eq.0.) K02A=0.
      IF(K02A.lt.0.) K02A=0.

C-- K02=K02/int(GAMMA**3/2) between -H and 0
C--- CA=((3/2)*Ko)^3/2 with Ko=4-> CA=14.7
      CA=14.7
C--- CD=((3/2)*Ko)^3/2 with Ko=4 -> CD=14.7
      CD=14.7
 
      INTGAMMAA=0.
      INTGAMMAD=0.

      KP2=0.
      KP2INT=0.
      KP2=GAMMA**(3./2.)
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO

      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      INTGAMMAA=KINT2(kkpint2)
      DEALLOCATE(KINT2)

      ALLOCATE(KINT2(kkpint))
      KINT2=0.
      CALL dqtfg(zint(kmi-(kkpint-1):kmi),KP2INT(kmi-(kkpint-1):kmi),
     *           KINT2,kkpint)
      INTGAMMAD=KINT2(kkpint)
      DEALLOCATE(KINT2)

c?    IF(K02A.EQ.0.) INTGAMMAA=0.
c?    IF(K02D.EQ.0.) INTGAMMAD=0.

      INTGAMMA=INTGAMMAA/CA+INTGAMMAD/CD
      K02=0.

      INTGAMMA1=0.
      KP2=0.
      KP2INT=0.
      KP2=GAMMA
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      ALLOCATE(KINT2(kmi))
      KINT2=0.
      CALL dqtfg(zint,KP2INT,KINT2,kmi)
      INTGAMMA1=KINT2(kmi)
      DEALLOCATE(KINT2)

      INTEXP=0.
      CD1=0.
      KP2=0.
      KP2INT=0.
      KP2=GAMMA*exp(-((z+abs(z(kmi)))**2)/(2.*4000.**2))
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      INTEXP=KINT2(kkpint2)
      DEALLOCATE(KINT2)
      CD1=INTEXP*sqrt(2./pi)*(abs(z(kmi)/4000.))

      IF(INTGAMMA1.ne.0.) INTGAMMAS=CD1*INTGAMMA/INTGAMMA1+INTGAMMAD/CD
      IF(INTGAMMAS.ne.0.) K02=(K02A+K02D)/INTGAMMAS

      IF(GAMMA(kmi).eq.1.) K02=0.

      Ro=0.
      IF(rd*f.NE.0.) Ro(1:kmi)=sqrt(K02*GAMMA)/(rd*f)

C-- Tony - 05/31/2013
c     IF(dabs(Ro(1)).GT.1..OR.Tap(1).LT.0..OR.Tap(1).GT.1.) THEN
      IF(dabs(Ro(1)).GT.1..OR.Tap(1).LT.0..OR.Tap(1).GT.1.
     *  .OR.minval(Tap(1:ktap)).LT.0.) THEN
        K02=0.
        GAMMA=1.
        ktap=1
        rdm=0.
        GO TO 101
      ENDIF

 101  CONTINUE

      IF (K02.le.0.) THEN
c       K02=0.d0
        K02=1.
        k02count=1.
      ELSE
        k02count=0.
      ENDIF
 
C---  Calculate udrift (ud,vd)
C---- Calculate ud0,vd0
      ALLOCATE(ud(kmi),vd(kmi))
      ud=0.
      vd=0.
c     ud = UB1AVE - 0.5 * f * (rd**2) * (-DZLYB1AVE)
c     vd = VB1AVE - 0.5 * f * (rd**2) * (+DZLXB1AVE)
C-- New ud,vd
      ALLOCATE(KINT2(kkpint2),GAMMAREV(kmi))
      KINT2=0.
      GAMMAREV=0.
      DO k=1,kmi
        GAMMAREV(k)=GAMMA(kmi-k+1)
      ENDDO
      CALL dqtfg(zint(1:kkpint2),sqrt(GAMMAREV(1:kkpint2)),
     *           KINT2,kkpint2)

      ALLOCATE(Gsx(kmi),Gsy(kmi),DZGAMMA(kmi),
     *         Gsxint(kkpint2),Gsyint(kkpint2))
      Gsx=0; Gsy=0.; DZGAMMA=0.; Gsxint=0; Gsyint=0.
      Isx=0; Isy=0.

      DZGAMMA=0.
      IF(kmi.gt.1.AND.K02.ne.1.) CALL d1sym(kmi,z,GAMMA,DZGAMMA,9)

      DO k=1,kmi
        Gsx(k)=DZGAMMA(kmi-k+1)*LX(kmi-k+1)/dsqrt(GAMMA(kmi-k+1))
        Gsy(k)=DZGAMMA(kmi-k+1)*LY(kmi-k+1)/dsqrt(GAMMA(kmi-k+1))
      ENDDO

      CALL dqtfg(zint(1:kkpint2),Gsx,Gsxint,kkpint2)
      CALL dqtfg(zint(1:kkpint2),Gsy,Gsyint,kkpint2)

      IF(KINT2(kkpint2).NE.0.) Isx=
     *  -(LX(ktap)/KINT2(kkpint2))*(1.-dsqrt(GAMMA(ktap)))
     *  -(1./KINT2(kkpint2))*(LX(kmi)*dsqrt(GAMMA(kmi)))
     *  -(0.5/KINT2(kkpint2))*Gsxint(kkpint2)

      IF(KINT2(kkpint2).NE.0.) Isy=
     *  -(LY(ktap)/KINT2(kkpint2))*(1.-dsqrt(GAMMA(ktap)))
     *  -(1./KINT2(kkpint2))*(LY(kmi)*dsqrt(GAMMA(kmi)))
     *  -(0.5/KINT2(kkpint2))*Gsyint(kkpint2)

      DEALLOCATE(Gsx,Gsy,DZGAMMA,Gsxint,Gsyint)

      SIGMAS=1.
      IF(K02.ne.1..AND.KINT2(kkpint2).ne.0.) THEN
C-- correction to ud, SIGMAS=SIGMAT/(1+SIGMAT),SIGMAT=1 -> SIGMAS=0.5
c     SIGMAS=1./(1.+1.)
      ud = UB1AVE - SIGMAS* f*(rd0**2)*(-Isy) - SIGMAS*beta*rd0**2
      vd = VB1AVE - SIGMAS* f*(rd0**2)*(+Isx)
c     ud = UB1AVE - SIGMAS* f*(rd**2)*(-Isy) - SIGMAS*beta*rd**2
c     vd = VB1AVE - SIGMAS* f*(rd**2)*(+Isx)
C--
      ENDIF
      DEALLOCATE(KINT2,GAMMAREV)

C-- Calculate factor M
      ALLOCATE(MK(kmi))
c     MK=0.
      MK=1.
c     IF(K02.ne.0.) MK=1./(1.+0.75*((U-ud)**2+(V-vd)**2)/K02)
      DO k=1,kmi
c     IF(K02*GAMMA(k).ne.0.)
      IF(K02.ne.1..AND.K02*GAMMA(k).ne.0.)
     *MK(k)=1./(1.+0.5*((U(k)-ud(k))**2+(V(k)-vd(k))**2)/(K02*GAMMA(k)))
      ENDDO

C-- Calculate KAPPAM(Z)
      ALLOCATE(KAPPAMZ(kmi))
      KAPPAMZ=rd*dsqrt(K02*GAMMA)*MK
c     KAPPAMZ=rd*dsqrt(K02*GAMMA)
      KAPPAMZa=0.
      KAPPAMZa(1:kmi)=KAPPAMZ
c     DEALLOCATE(KAPPAMZ)

C-- Rosby number: Ro=kappam/fl^2
      IF(rd*f.NE.0.) Ro(1:kmi)=sqrt(K02*GAMMA)/(rd*f)
c     IF(igrid.EQ.it.AND.jgrid.EQ.jt) THEN
c        IF(igrid.eq.it.AND.jgrid.eq.jt.OR.
c    *      igrid.eq.it+1.AND.jgrid.eq.jt.OR.
c    *      igrid.eq.it.AND.jgrid.eq.jt+1) THEN
c       WRITE(*,*)'DEBUG-Ro:',nstep,igrid,jgrid,ktap,kmi,K02,rd,f
c       DO k=1,ktap
c         WRITE(*,*)k,Ro(k)
c       ENDDO
c     ENDIF
c     IF(dabs(Ro(1)).lt.1.) WRITE(*,*)nstep,igrid,jgrid,Ro(1)


c     IF(igrid.EQ.it.AND.jgrid.EQ.jt) THEN
cc       IF(igrid.eq.it.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it+1.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it.AND.jgrid.eq.jt+1) THEN
c       WRITE(*,*)'DEBUG-KAPPAMZ-0:',nstep,igrid,jgrid,ktap,kmi,K02,rd
c       DO k=1,ktap
c         WRITE(*,114)k,KAPPAMZ(k),MK(k),rd*dsqrt(K02*GAMMA(k))
c       ENDDO
c     ENDIF
c114  FORMAT(I,3E)

C-- Calculate KAPPA vector
      ALLOCATE(KAPPAX(kmi),KAPPAY(kmi))
      KAPPAX=0.
      KAPPAY=0.
c     KAPPAX=KAPPAMZ*z*(F1X+F2X)
c     KAPPAY=KAPPAMZ*z*(F1Y+F2Y)
      KAPPAX=z*(F1X+F2X)
      KAPPAY=z*(F1Y+F2Y)

C-- Calculate the mesoscale vertical flux for buoyancy Fv(M,b)
      ALLOCATE(FvM(kmi))
      FvM = 0.
      IF(K02.GT.1.) THEN
        FvM=-KAPPAMZ*z*((F1X+F2X)*dxb+(F1Y+F2Y)*dyb)
      ENDIF
      DO k=1,kmi
        IF(FvM(k).lt.0.) THEN
          FvM(k)=0.
c         k02count=1.
c       ELSE
c         k02count=0.
        ENDIF
      ENDDO

C-- Calculate KAPPAl vector
      ALLOCATE(KAPPAlx(kmi),KAPPAly(kmi))
      KAPPAlx = 0.
      KAPPAly = 0.
      IF(K02.gt.1.) THEN
        DO k=1,kmi
          IF((dxb(k)**2+dyb(k)**2).ne.0.) THEN
            KAPPAlx(k) = -FvM(k)*dxb(k)/(dxb(k)**2+dyb(k)**2)
            KAPPAly(k) = -FvM(k)*dyb(k)/(dxb(k)**2+dyb(k)**2)
          ENDIF
        ENDDO
      ENDIF

C-- Calculate vector OMEGA
c     OMEGAx = 0.
c     OMEGAy = 0.
      OMEGAxT = 0.
      OMEGAyT = 0.
      OMEGAx1mT = 0.
      OMEGAy1mT = 0.
c     IF(K02.gt.1.) THEN
c     IF(K02.gt.1..AND.dabs(Ro(1)).LT.1.) THEN
c     IF(K02.gt.1..AND.Tap(1).GT.0..AND.n20(1).gt.0.) THEN
c     IF(K02.gt.1..AND.Tap(1).GT.0..AND.minval(n20(1:ktap)).GT.0.) THEN
      kappamsigx=0.; kappamsigy=0.
      kappamsigmsx=0.; kappamsigmsx=0.
      LX3d(1:kmi)=0.; LY3d(1:kmi)=0.
      IF(K02.gt.1..AND.minval(Tap(1:ktap)).GT.0.
     *  .AND.minval(n20(1:kmi)).GT.0.) THEN
        IF(dabs(Ro(1)).LT.1.) THEN
        IF(Tap(1).LE.1.) THEN
        IF(rd.LE.10000000.) THEN
        IF(maxval(KAPPAMZ(1:kmi)).LT.15000.*10000.) THEN
        IF(k02count.EQ.0.) THEN
c       OMEGAx(1:kmi) = Tap*KAPPAlx + (1.-Tap)*KAPPAX
c       OMEGAy(1:kmi) = Tap*KAPPAly + (1.-Tap)*KAPPAY
c       OMEGAx(1:ktap) = Tap*KAPPAlx + (1.-Tap)*KAPPAX
c       OMEGAy(1:ktap) = Tap*KAPPAly + (1.-Tap)*KAPPAY
        OMEGAxT(1:ktap) = Tap*KAPPAlx
        OMEGAyT(1:ktap) = Tap*KAPPAly
        OMEGAx1mT(1:ktap) = (1.-Tap)*KAPPAX
        OMEGAy1mT(1:ktap) = (1.-Tap)*KAPPAY
c       WRITE(*,*)'YES',nstep,igrid,jgrid,Tap(1),K02,minval(n20(1:ktap))
c       ELSE
c       WRITE(*,*)'NO',nstep,igrid,jgrid,Tap(1),K02,minval(n20(1:ktap))
C-- Calculate kappam x sigma
        ALLOCATE(DZKAPPAMZ(kmi),Iintx(kmi),Iinty(kmi),Ix(kmi),Iy(kmi),
     *           alpha(kmi),Ixd(kmi),Iyd(kmi),alphad(kmi),
     *           ung(kmi),vng(kmi),alphax(kmi),alphay(kmi))

C-- Calculate ung
        ung=0.; vng=0.
        DO k=1,kmi
          IF(SIGMAS*f*rd**2.NE.0.) THEN
            ung(k)=(KAPPAMZ(k)/(SIGMAS*f*rd**2))
     *         *(-(V(k)-(VB1AVE - SIGMAS* f*(rd**2)*(+Isx))))
            vng(k)=(KAPPAMZ(k)/(SIGMAS*f*rd**2))
     *         *(+(U(k)-(UB1AVE - SIGMAS* f*(rd**2)*(-Isy))))
          ENDIF
        ENDDO

C-- Calculate alpha
        DZKAPPAMZ=0.
        IF(kmi.gt.1) CALL d1sym(kmi,z,KAPPAMZ,DZKAPPAMZ,10)

        alphax=0.; alphay=0.
        alphax=ung-LX*DZKAPPAMZ
        alphay=vng-LY*DZKAPPAMZ

        Iintx=0.; Iinty=0.
        DO k=1,kmi
          Iintx(k)=alphax(kmi-k+1)
          Iinty(k)=alphay(kmi-k+1)
        ENDDO
        Ix=0.; Iy=0.
        CALL dqtfg(zint,Iintx,Ix,kmi)
        CALL dqtfg(zint,Iinty,Iy,kmi)
 

        Ixd=0.; Iyd=0.
        DO k=1,kmi
          Ixd(k)=Ix(kmi-k+1)
          Iyd(k)=Iy(kmi-k+1)
        ENDDO
         
        kappamsigx(1:kmi) = Ixd
        kappamsigy(1:kmi) = Iyd
        kappamsigmsx(1:kmi) = kappamsigx(1:kmi)+KAPPAMZ*LX
        kappamsigmsy(1:kmi) = kappamsigy(1:kmi)+KAPPAMZ*LY

        DEALLOCATE(DZKAPPAMZ,Iintx,Iinty,Ix,Iy,alpha,Ixd,Iyd,alphad,
     *             ung,vng)

        LX3d(1:kmi)=LX; LY3d(1:kmi)=LY
        
        ENDIF
        ENDIF
        ENDIF
        ENDIF
        ENDIF
      ENDIF


c     IF(igrid.EQ.it.AND.jgrid.EQ.jt) THEN
cc       IF(igrid.eq.it.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it+1.AND.jgrid.eq.jt.OR.
cc   *      igrid.eq.it.AND.jgrid.eq.jt+1) THEN
c       WRITE(*,*)'DEBUG-OMEGA-0:',nstep,igrid,jgrid,ktap,kmi,K02,
c    *   minval(Tap(1:ktap)),minval(n20),rd/100000.
c       DO k=1,ktap
c         WRITE(*,113)k,OMEGAxT(k),OMEGAx1mT(k),OMEGAyT(k),OMEGAy1mT(k),
c    *      Tap(k),FvM(k)
c       ENDDO
c     ENDIF
c113  FORMAT(I,6E)

      Ustar=0.
      Vstar=0.

      DEALLOCATE(KAPPAX,KAPPAY,FvM,KAPPAlx,KAPPAly)
      DEALLOCATE(KAPPAMZ)
      DEALLOCATE(Tap)
      DEALLOCATE(ud,vd,MK)
      DEALLOCATE(GAMMA)
      DEALLOCATE(n2sm)
      DEALLOCATE(n2t,n2t1,z,zm,b1,b1t,
     *  U,V,rho,n2tr,
     *  LX,LY,DZLX,DZLY,
     *  UINT,VINT,UT,VT,
     *  dxb,dyb,KP,
     *  zint,UREV,VREV,KPINT,
     *  b1rev,DZLXREV,DZLYREV,
     *  UTINT,VTINT,DZU,DZV,
     *  F1X,F1Y,F2X,F2Y,KP2,KP2INT)
 
 10   CONTINUE

 
      END  subroutine mesoscales1d


C*****
      SUBROUTINE B1AVERAGE(F,Z,B1,K,FAVE,KML)
      implicit none
      INTEGER K,KML
      REAL*8 F(K),Z(K),B1(K),FAVE,FB1AVE(K-KML+1),
     *  B1AVE(K-KML+1),a02
C
      a02=0.03
c     a02=0.5
C
      CALL dqtfg(Z(1:K-KML+1),F(1:K-KML+1)
     *           *dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           FB1AVE,K-KML+1)     
      CALL dqtfg(Z(1:K-KML+1),dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           B1AVE,K-KML+1)
C     WRITE(*,*),"B1AVE",FB1AVE/B1AVE
      FAVE=0.
      IF(B1AVE(K-KML+1).NE.0.)
     *FAVE=FB1AVE(K-KML+1)/B1AVE(K-KML+1)
C
      RETURN
      END
C*****
      SUBROUTINE B1AVERAGEG(F,Z,B1,K,FAVE,KML,GAMMA)
      implicit none
      INTEGER K,KML,kk
      REAL*8 F(K),Z(K),B1(K),FAVE,FB1AVE(K-KML+1),
     *  B1AVE(K-KML+1),a02,GAMMA(K),GAMMAREV(K)

      GAMMAREV=0.
      DO kk=1,K
        GAMMAREV(kk)=GAMMA(K-kk+1)
      ENDDO
C
      CALL dqtfg(Z(1:K-KML+1),F(1:K-KML+1)*dsqrt(GAMMAREV(1:K-KML+1)),
     *           FB1AVE,K-KML+1)
      CALL dqtfg(Z(1:K-KML+1),dsqrt(GAMMAREV(1:K-KML+1)),
     *           B1AVE,K-KML+1)
C     WRITE(*,*),"B1AVE",FB1AVE/B1AVE
      FAVE=0.
      IF(B1AVE(K-KML+1).NE.0.)
     *FAVE=FB1AVE(K-KML+1)/B1AVE(K-KML+1)
C
      RETURN
      END
C*****
        SUBROUTINE d1sym(n,x,y,dyodx,flg)
C080626 Input n,x(n),y(n)               Output dyodx(n)
C       Calculates the numerical ordinary derivative of y with respect to x, dyodx,
C       using the symmetrical approach of considering the nearest neighbor points
C       on both sides and assuming the derivative is changing linearly.
C       The ratio of differences on each side thus give the midpoint derivatives
C       on the respective sides and the averages of the midpoint derivatives
C       weighted each by the other midpoint's  distance give the derivative at the point.
C       Note for equal spacing same as ratio of neighbor differences with each other.
C       At top and bottom revert to ratio of differences with lone nearest neighbor.
        IMPLICIT NONE
        INTEGER n,flg
        REAL*8 x(n),y(n),dyodx(n)
        REAL*8 dyodxm,dyodxp,dxm,dxp,dym,dyp

        INTEGER i

        if(n.eq.1) then
          dyodx=0.
        else
        DO i=1,n
           IF(i.EQ.1) THEN
c             write(20,'(3i4,9e14.4)') flg,n,i,x(i+1),x(i),y(i+1),y(i)
c             dyodx(i) = (y(i+1)-y(i))/(x(i+1)-x(i))
              dyodx(i) = (y(i+1)-y(i))/(x(i+1)-x(i)+1d-20)
           ELSE IF(i.EQ.n) THEN
c             dyodx(i) = ((y(i)-y(i-1))/(x(i)-x(i-1)))
C-- TONY - 03/09/11 - updated the bottom calculation with a 3 points derivative
              IF(n.eq.2) THEN
                dyodx(i) = 0.
              ELSE
                dxm = x(i-1)-x(i-2)
                dxp = x(i)-x(i-1)
                dym = y(i-1)-y(i-2)
                dyp = y(i)-y(i-1)
c             write(21,'(3i4,9e14.4)') flg,n,i,x(i-1),x(i)
c    &                                        ,y(i-1),y(i)
c               dyodxm = dym/(dxm)
c               dyodxp = dyp/(dxp)
c               dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp)
                dyodxm = dym/(dxm+1d-20)
                dyodxp = dyp/(dxp+1d-20)
                dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp+1d-20)
              ENDIF
           ELSE
              dxm = x(i)-x(i-1)
              dxp = x(i+1)-x(i)
              dym = y(i)-y(i-1)
              dyp = y(i+1)-y(i)
c             write(22,'(3i4,9e14.4)') flg,n,i,x(i-1),x(i),x(i+1)
c    &                                        ,y(i-1),y(i),y(i+1)
c             dyodxm = dym/(dxm)
c             dyodxp = dyp/(dxp)
c             dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp)
              dyodxm = dym/(dxm+1d-20)
              dyodxp = dyp/(dxp+1d-20)
              dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp+1d-20)
           END IF
        END DO
        endif

        RETURN
        END
C*****
C080625-30AH Subroutine to calculate the first baroclinic mode in CD WKB approximation
C       from an input N^2 profile, artificially treating N^2 as zero where it is less than
C       zero since its square and fourth roots appear in the WKB approximation formula,
C       B_1(z)=[N(z)/N_max]^1/2 cos[(f r_d)^-1 \Integral_z^z_max N(z) dz (where z_max was
C       called "z_0" by CD and take (f r_d) = {\Integral -H^z_max N(z) dz \over pi},
C       but with B_1(z) kept at 1 above z_max as formula wasn't intended for that range,
C       corrected and adapted from Cheng's 080522 program that calls dqtfg , named B1b.f:
!@sum B1b.f based on B1a.f but only reads form one data file
!@+   and ignors the last line in the data file
!@    5-22-2008

        SUBROUTINE baroclin1st(z,n2,m,ba1,im,Nm,frd,igrid,jgrid)
C       Inputs:
C               z(m)    !Depth                                                  [m]
C               n2(m)   !Square of Brunt Vaisala frequency                      [s^-2]
C
C       Outputs:
C               ba1(m)  !First baroclinic mode
C               im      !Index of maximum N on column
C               Nm      !Maximum N on column                                    [s^-2]
C               frd     !(f r_d)                                                [m/s]
C
C
C       Internal:
C               lifout  !Logical switch set .TRUE. iff diagnostic output here

      implicit none
C080625AH
        integer igrid,jgrid
        INTEGER m       ! Number of ocean levels
        REAL*8 ba1(m)   ! Array for first baroclinic mode as calculated by this routine
        LOGICAL lifout
        PARAMETER(lifout=.TRUE.)
      real*8 n2,z,n,ni,pi,a,nm,arg,B1
      dimension n2(m)   !Array for N^2
      dimension z(m),n(m),ni(m) !Arrays for depth, SQRT(N^2), and z integral of N^2
        REAL*8 frd      !"f r_d" calculated in this routine
C******AH
C080625AH Work variables introduced.
        INTEGER im,i
C******AH
C080630AH Depth integrations are only done from maximum N level, "the pycnocline depth",
C       to the bottom so points above the depth of Maximum N are excluded unlike in B1b.f.
C       Introduce new arrays for depth and SQRT(N^2) which start at "pycnocline depth".
        REAL*8 zsubpyc(m),nsubpyc(m)
        INTEGER isubpyc(m),msubpyc
C******AH


      pi=acos(-1.d0)
      a=.1


C080625AH Canuto's WKB approximation formula contains (N^2)^1/2 and (N^2)^(1/4)
C       and therefore CANNOT HANDLE NEGATIVE N^2. Artificially set N^2=0 .
        DO i=1,m
           n2(i)=MAX(0.D0,n2(i))
        END DO
C******AH


      do i=1,m
         n(i)=sqrt(n2(i))
      end do

      ! find Nm = N_max
      Nm=N(1)
C080625AH Write out maximum N and depth when lifout=.TRUE. .
      im=1
      do i=2,m
         IF(N(i).GT.Nm) im=i
         Nm=max(N(i),Nm)
      end do
      IF(lifout) THEN
c       IF(ik.eq.177.and.jk.eq.163) THEN
c       write(*,*) "Nm=",Nm
c       WRITE(*,*) "im=",im
c       ENDIF
      END IF
C******AH

C080630AH Fill z & N arrays which go from max. N level down and integrate N on same range.
        DO i=im,m
           isubpyc(i)=i+1-im
           zsubpyc(isubpyc(i))=z(i)
           nsubpyc(isubpyc(i))=n(i)
        END DO
        msubpyc=m+1-im

      call dqtfg(zsubpyc,nsubpyc,ni,msubpyc)
C******AH

      do i=1,m ! z in meters
C80630AH Integrate N only over subpycnocline to find "f r_d". Keep B_1=1 above pycnocline.
         frd=ni(msubpyc)/pi
         IF(i.LT.im) THEN
           B1=1.D0
         ELSE
           if(frd.ne.0.d0) then           !Natassa
           arg=1./frd*ni(isubpyc(i))
           B1=(nsubpyc(isubpyc(i))/nm)**.5 * cos(arg)
           else
           B1=1.D0
           endif
         END IF
C******AH
C080625AH Store first baroclinic mode values in an array.
         ba1(i)=B1
C        Gam=1./(1+a) * (a+B1**2)
C        DMW=1./4.+3./4.*exp(-abs(z(i))/500.)
C        write(23,'(9e14.6)') z(i),Gam,sqrt(Gam)
C &        ,n2(i)/n2(1),DMW
C******AH
c        IF(ik.eq.177.and.jk.eq.163.and.i.eq.m) THEN
c          WRITE(*,*)'b1:ba1=',i,ba1(i)
c          WRITE(*,*)'b1:nsubpyc(isubpyc(i))=',nsubpyc(isubpyc(i))
c          WRITE(*,*)'b1:nm=',nm
c          WRITE(*,*)'b1:cos(arg)=',cos(arg)
c          WRITE(*,*)'b1:arg=',arg
c          WRITE(*,*)'b1:frd=',frd
c          WRITE(*,*)'b1:ni(isubpyc(i))=',ni(isubpyc(i))
c        ENDIF
         end do

c     IF(im.eq.m) WRITE(*,*)'b1:im,m,frd=',ik,jk,im,m,frd

      end
C*****
      subroutine dqtfg(x,y,z,n)
      !@sum integrate y from x(1) to x(i) and store it in z(i)
      implicit none
      real*8 sum1,sum2
      integer n,i
      real*8 x(n),y(n),z(n)
      sum2=0.
      if(n-1)4,3,1
    1 do 2 i=2,n
      sum1=sum2
      sum2=sum2+.5d0*(x(i)-x(i-1))*(y(i)+y(i-1))
    2 z(i-1)=sum1
    3 z(n)=sum2
    4 return
      end

C*********
      subroutine z121_mesosc (v0,kmtj,km)
!@sum z121 Apply 121 smoothing in k to 2-d array V(k=1,km)
!@+   top (0) value is used as a dummy
!@+   bottom (km+1) value is set to input value from above.
      IMPLICIT NONE
      REAL*8, PARAMETER :: p5=5d-1, p25=2.5d-1
      INTEGER, INTENT (IN) :: kmtj,km
      REAL*8, INTENT (INOUT) :: V0(km)  ! 2-D array to be smoothed
      REAL*8 V(0:km+1)
      INTEGER K
      REAL*8 tmp
      V(1:km)=V0
      V(0)      =  p25 * V(1)
      V(kmtj+1) =        V(kmtj)

      do k=1,kmtj
         tmp      =  V(k)
         V(k)   =  V(0)  + p5 * V(k) + p25 * V(k+1)
         V(0)   =  p25 * tmp
      end do
      V0=V(1:km)
      return
      end subroutine z121_mesosc
C************

#ifdef OCN_GISS_MESO
c     SUBROUTINE MESO_D(TRM0,TRXMO0,TRYMO0,TRM,TRXMO,TRYMO,TRZMO)
      SUBROUTINE MESO_D(TRM,TRXMO0,TRYMO0,TRXMO,TRYMO,TRZMO)
      USE MODEL_COM,  only : nstep=>itime
      USE MESOML, only: ktap,ogrid,kappam3d0,OMEGAx3d,OMEGAy3d,
     *  Uplus,Vplus,Wplus,idm,kdm,dxpo,dypo,jdm
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
     ., HALO_UPDATE, NORTH, SOUTH
      USE OCEAN, only : lmm,lmu,lmv,focean
      USE OCEAN, only : dxypo,mo,dypo,uo
      USE OCEAN,      only : flux_x_sm,flux_y_sm,flux_z_sm,kpl
      USE GM_COM, ONLY: grid,GETDomainBounds

      IMPLICIT NONE

      INTEGER i,j,i_0,i_1,j_0,j_1,ip1,im1,k,l,it,jt,itp1,jtp1,it0,jt0
      INTEGER kx,ky,kz,flux_y0(idm,ogrid%j_strt_halo:ogrid%j_stop_halo)
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     *  TRM0,TRXMO0,TRYMO0,TRM,TRXMO,TRYMO,TRZMO,Fdiffx,Fdiffy,
     *  flux_x,flux_y,flux_z,dxFdiffx,dyFdiffy
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     *  dyFdiffx,dxFdiffy
      REAL*8, DIMENSION(ogrid%j_strt_halo:ogrid%j_stop_halo) :: BYDYP
      REAL*8 MOFY
      REAL*8 STRNP(kdm)
      REAL*8 wta3m

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      CALL HALO_UPDATE(ogrid,
     *                ktap(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
     *                FROM=SOUTH+NORTH)

c     DO j=ogrid%J_STRT,ogrid%J_STOP
c     IF(j.eq.127) THEN
c       WRITE(*,*)'DEBUG-ktap-MESO_D:',ktap(42,j)
c       WRITE(*,*)'DEBUG-kappam3d-MESO_D:',kappam3d0(42,j,1)
c     ENDIF
c     ENDDO

      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP

C-- Fdiff = -(kappam*grad(tau) + OMEGA*dtaudz)
      Fdiffx = 0.
      Fdiffy = 0.
      Fdiffx = -(kappam3d0*TRXMO0 + OMEGAx3d*TRZMO)
      Fdiffy = -(kappam3d0*TRYMO0 + OMEGAy3d*TRZMO)

C--- Calculate dxFdiffx,dyFdiffy
c     dxFdiffx = 0.
c     dyFdiffy = 0.
c     CALL HALO_UPDATE(ogrid,
c    *                Fdiffy(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                dyFdiffy(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)
c     do 102 j=j_0,j_1
c     do 102 i=i_0,i_1
c       ip1=i+1
c       im1=i-1
c       IF(i.eq.1) im1=288
c       IF(i.eq.288) ip1=1
c       do k=1,ktap(i,j)
C-- dxFdiffx
c       IF (Fdiffx(i,j,k).ne.0.) THEN
c         IF(Fdiffx(ip1,j,k).ne.0.) THEN
c          dxFdiffx(i,j,k)=(Fdiffx(ip1,j,k)-Fdiffx(i,j,k))/dxpo(j)
c         ELSEIF(Fdiffx(ip1,j,k).eq.0..AND.Fdiffx(im1,j,k).ne.0.
c    *           .AND.i.ne.288) THEN
c          dxFdiffx(i,j,k)=dxFdiffx(im1,j,k)
c         ENDIF
c       ENDIF
C-- dyFdiffy
c       IF (Fdiffy(i,j,k).ne.0.) THEN
c         IF(Fdiffy(i,j+1,k).ne.0.) THEN
c          dyFdiffy(i,j,k)=(Fdiffy(i,j+1,k)-Fdiffy(i,j,k))/dypo(j)
c         ELSEIF(Fdiffy(i,j+1,k).eq.0..AND.Fdiffy(i,j-1,k).ne.0.) THEN
c         ELSEIF(Fdiffy(i,j+1,k).eq.0..AND.Fdiffy(i,j-1,k).ne.0.
c    *           .OR.j.eq.jdm) THEN
c          dyFdiffy(i,j,k)=dyFdiffy(i,j-1,k)
c         ENDIF
c       ENDIF
c       enddo
c       IF(MINVAL(OMEGAx3d(i,j,:)).EQ.0.
c    &.AND.MAXVAL(OMEGAx3d(i,j,:)).EQ.0.) THEN
c        dxFdiffx(i,j,:)=0.
c       ENDIF
c       IF(MINVAL(OMEGAy3d(i,j,:)).EQ.0.
c    &.AND.MAXVAL(OMEGAy3d(i,j,:)).EQ.0.) THEN
c        dyFdiffy(i,j,:)=0.
c       ENDIF
c       IF(ktap(i,j).eq.1) THEN
c        dxFdiffx(i,j,:)=0.
c        dyFdiffy(i,j,:)=0.
c       ENDIF
c102  CONTINUE

c     CALL get_gradients(mo,Fdiffy,0,dxFdiffyi,dyFdiffyi)
c     CALL get_gradients(mo,Fdiffy,1,dxFdiffye,dyFdiffye)

      CALL get_gradients(mo,Fdiffx,0,dxFdiffx,dyFdiffx)
      CALL get_gradients(mo,Fdiffy,0,dxFdiffy,dyFdiffy)

      do 102 j=j_0,j_1
      do 102 i=i_0,i_1
        IF(SUM(OMEGAx3d(i,j,:)).EQ.0.) THEN
          dxFdiffx(i,j,:)=0.
        ENDIF
        IF(SUM(OMEGAy3d(i,j,:)).EQ.0.) THEN
          dyFdiffy(i,j,:)=0.
        ENDIF
        IF(ktap(i,j).eq.1) THEN
          dxFdiffx(i,j,:)=0.
          dyFdiffy(i,j,:)=0.
        ENDIF
 102  CONTINUE

c     it=42; jt=127
      it0=42; jt0=127
      itp1=it0+1; jtp1=jt0+1

      it=it0; jt=jt0
c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         WRITE(*,*)'DEBUG-dyFdiffy:',nstep,it,jt,
c    *        MINVAL(OMEGAy3d(it,J,:)),MAXVAL(OMEGAy3d(it,J,:))
cc        WRITE(*,*)'DEBUG-NP:',nstep,HAVE_NORTH_POLE
c         DO L=1,ktap(it,J)
c           WRITE(*,115)L,dyFdiffy(it,J,L),Fdiffy(it,J,L),
c    *        Fdiffy(it,J+1,L),kappam3d0(it,J,L),kappam3d0(it,J+1,L),
c    *        OMEGAy3d(it,J,L),OMEGAy3d(it,J+1,L)
cc   *        ,dyFdiffyi(it,J,L),dyFdiffye(it,J,L)
c         ENDDO
c       ENDIF
c     ENDDO

 115  FORMAT(I,9E)

      

C- fluxes
      flux_x=0.; flux_y=0.; flux_z=0.

c     CALL HALO_UPDATE(ogrid,
c    *                LMV(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
c    *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                Vplus(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                TRYMO0(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                dyFdiffy(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)

      DO J=j_0,j_1
        IM1 = idm
        DO I=i_0,i_1
          DO L=1,ktap(I,J)
            IF(LMU(IM1,J).ge.L) THEN
              flux_x(I,J,L) = Uplus(I,J,L)*TRXMO0(I,J,L)+dxFdiffx(I,J,L)
            ENDIF
            IF(LMV(I,J-1).ge.L) THEN
              flux_y(I,J,L) = Vplus(I,J,L)*TRYMO0(I,J,L)+dyFdiffy(I,J,L)
c             flux_y(I,J,L) = TRZMO(I,J,L)
            ENDIF
c           IF(ktap(I,J).gt.L) THEN
            IF(ktap(I,J).ge.L.AND.L.NE.1.AND.L.NE.LMM(I,J)) THEN
              flux_z(I,J,L) = Wplus(I,J,L)*TRZMO(I,J,L)
            ENDIF
          ENDDO
          IM1=I
        ENDDO
      ENDDO


      wta3m=exp(-1./(3.*1440.))

      i=42
c     DO j=j_0,j_1
c     IF(j.eq.127) WRITE(*,*)'DEBUG-flux_x_sm before:',
c    *  nstep,i,j,flux_x(i,j,1),flux_x_sm(i,j,1)
c     ENDDO

      flux_x_sm=wta3m*flux_x_sm + (1.-wta3m)*flux_x
      flux_y_sm=wta3m*flux_y_sm + (1.-wta3m)*flux_y
      flux_z_sm=wta3m*flux_z_sm + (1.-wta3m)*flux_z

      flux_x=flux_x_sm
      flux_y=flux_y_sm
      flux_z=flux_z_sm

c     DO j=j_0,j_1
c     IF(j.eq.127) WRITE(*,*)'DEBUG-flux_x_sm after:',
c    *  nstep,i,j,flux_x(i,j,1),flux_x_sm(i,j,1)
c     ENDDO



c     CALL HALO_UPDATE(ogrid,
c    *                LMM(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
c    *                FROM=SOUTH+NORTH)

c     CALL HALO_UPDATE(ogrid,
c    *                TRM0(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH)

C??? HALO flux_x, flux_y ???


c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         IF(SUM(flux_x(it,J,:)).NE.0.) THEN
c         WRITE(*,*)'DEBUG-flux_x:OCN_mesosc',nstep,it,J,ktap(it,J),
c    *      lmm(it,J),LMU(it-1,J),focean(it,j)
c         DO L=1,ktap(it,J)
c           WRITE(*,111)L,flux_x(it,J,L),Uplus(it,J,L)*TRXMO0(it,J,L),
c    *       dxFdiffx(it,J,L),Uplus(it,J,L),TRM0(it,J,L),TRM0(it-1,J,L),
c    *       TRM(it,J,L),TRM(it-1,J,L)
c         ENDDO
c         ENDIF
c       ENDIF
c     ENDDO

c     it=itp1; jt=jt0
c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         IF(SUM(flux_x(it,J,:)).NE.0.) THEN
c         WRITE(*,*)'DEBUG-flux_x:OCN_mesosc',nstep,it,J,ktap(it,J),
c    *      lmm(it,J),LMU(it-1,J),focean(it,j)
c         DO L=1,ktap(it,J)
c           WRITE(*,111)L,flux_x(it,J,L),Uplus(it,J,L)*TRXMO0(it,J,L),
c    *       dxFdiffx(it,J,L),Uplus(it,J,L),TRM0(it,J,L),TRM0(it-1,J,L),
c    *       TRM(it,J,L),TRM(it-1,J,L)
c         ENDDO
c         ENDIF
c       ENDIF
c     ENDDO

c     it=it0; jt=jt0
c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         IF(SUM(flux_y(it,J,:)).NE.0.) THEN
c         WRITE(*,*)'DEBUG-flux_y:OCN_mesosc',nstep,it,J,ktap(it,J),
c    *      lmm(it,J),LMV(it,J-1),focean(it,j)
c         DO L=1,ktap(it,J)
c           WRITE(*,111)L,flux_y(it,J,L),Vplus(it,J,L)*TRYMO0(it,J,L),
c    *       dyFdiffy(it,J,L),Vplus(it,J,L),TRM0(it,J,L),TRM0(it,J-1,L),
c    *       TRM(it,J,L),TRM(it,J-1,L)
c         ENDDO
c         ENDIF
c       ENDIF
c     ENDDO

c     it=it0; jt=jtp1
c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         IF(SUM(flux_y(it,J,:)).NE.0.) THEN
c         WRITE(*,*)'DEBUG-flux_y:OCN_mesosc',nstep,it,J,ktap(it,J),
c    *      lmm(it,J),LMV(it,J-1),focean(it,j)
c         DO L=1,ktap(it,J)
c           WRITE(*,111)L,flux_y(it,J,L),Vplus(it,J,L)*TRYMO0(it,J,L),
c    *       dyFdiffy(it,J,L),Vplus(it,J,L),TRM0(it,J,L),TRM0(it,J-1,L),
c    *       TRM(it,J,L),TRM(it,J-1,L)
c         ENDDO
c         ENDIF
c       ENDIF
c     ENDDO
c111  FORMAT(I,8E)

c     it=it0; jt=jt0
c     DO J=J_0,J_1
c       IF(J.EQ.jt) THEN
c         IF(SUM(flux_z(it,J,:)).NE.0.) THEN
c         WRITE(*,*)'DEBUG-flux_z:OCN_mesosc',nstep,it,J,ktap(it,J),
c    *      lmm(it,J),focean(it,j)
c         DO L=1,ktap(it,J)
c           WRITE(*,121)L,flux_z(it,J,L),Wplus(it,J,L)*TRZMO(it,J,L),
c    *       Wplus(it,J,L),TRM0(it,J,L),TRM0(it,J,L+1)
c         ENDDO
c         ENDIF
c       ENDIF
c     ENDDO
c121  FORMAT(I,5E)

c     flux_y0=0
c     DO J=j_0,j_1
c       DO I=i_0,i_1
c         IF(SUM(flux_y(I,J,:)).EQ.0.) THEN
c           flux_y0(I,J)=1
c           GO TO 211
c         ENDIF
c         DO L=1,ktap(I,J)
c           IF(SUM(flux_y(I,J,:)).NE.0..AND.flux_y(I,J,L).EQ.0.) THEN
c    *         WRITE(*,*)nstep,I,J,L,flux_y(I,J,L),SUM(flux_y(I,J,:))
c             flux_y0(I,J)=1
c             GO TO 211
c           ENDIF
c         ENDDO
c211      CONTINUE
c       ENDDO
c     ENDDO

c     BYDYP=0.
c     DO J=J_0,J_1
c       BYDYP(J)=1d0/DYPO(J)
c     END DO
c     MOFY=0.
c     DO J=j_0,j_1
c       IF(J.EQ.jt.OR.J.EQ.jt+1) THEN 
c       DO I=57,63
c       DO L=1,ktap(I,J)
c       MOFY =((MO(I,J-1,L)*BYDYP(J-1)*DXYPO(J-1)) +
c    *         (MO(I,J  ,L)*BYDYP(J  )*DXYPO(J  )))*0.5
c       WRITE(*,*)'DEBUG-MO:',I,J,MOFY,MO(I,J-1,L),MO(I,J,L),
c    * BYDYP(J-1),BYDYP(J),DXYPO(J-1),DXYPO(J),focean(I,J),flux_y(I,J,L)
c       ENDDO
c       ENDDO
c       ENDIF
c     ENDDO

c     CALL HALO_UPDATE(ogrid,
c    *                TRM(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,
     *                flux_y(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                ktap(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
c    *                FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,
     *                TRM(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                FROM=SOUTH)
c     CALL HALO_UPDATE(ogrid,
c    *                LMV(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
c    *                FROM=SOUTH+NORTH)
c     CALL HALO_UPDATE(ogrid,
c    *                LMM(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
c    *                FROM=SOUTH+NORTH)

      kx=0; ky=0.; kz=0.
      DO J=j_0,j_1
        IF(J.LT.JDM) THEN
        IM1 = idm
        DO I=i_0,i_1
c         IF(I.eq.62.AND.J.EQ.16) GO TO 210
c         IF(I.GE.61.AND.I.LE.63.AND.J.GE.15.AND.J.LE.17) GO TO 210
c         IF(I.GE.62.AND.I.LE.63.AND.J.GE.16.AND.J.LE.17) GO TO 210
c         IF(I.GE.56.AND.I.LE.58.AND.J.GE.15.AND.J.LE.17) GO TO 210
          IP1 = I+1
          IF(I.EQ.288) IP1=1
c         IF(FOCEAN(I,J).gt.0.) THEN
          DO L=1,ktap(I,J)
c         DO L=1,1
          IF(LMM(I,J).gt.0.AND.TRM(I,J,L).ne.0.) THEN
C-- X-direction
            IF(LMU(IM1,J).ge.L) THEN
            IF(Uplus(I,J,L).NE.0..AND.Uplus(IP1,J,L).NE.0.) THEN
cc          IF(OMEGAx3d(I,J,L).NE.0..AND.OMEGAx3d(IM1,J,L).NE.0.) THEN
            IF(LMM(IM1,J).NE.0) THEN
C??? Still need this one ^ ???
c             TRM(I  ,J,L) = TRM0(I  ,J,L) + flux_x(I,J,L)
              IF(ktap(IM1,J).gt.L) THEN
c             IF(MAXVAL(flux_x(I,J,:)).NE.0..AND.
c    *           MINVAL(flux_x(I,J,:)).NE.0.) THEN
              IF(SUM(flux_x(I,J,:)).NE.0.) THEN
              TRM(I  ,J,L) = TRM(I  ,J,L) + flux_x(I,J,L)
              TRM(IM1,J,L) = TRM(IM1,J,L) - flux_x(I,J,L)
c             ELSEIF(MAXVAL(flux_x(I,J,:)).EQ.0..AND.
c    *           MINVAL(flux_x(I,J,:)).EQ.0..AND.L.EQ.2) THEN
c                WRITE(*,*)'X:',nstep,I,J
              kx=kx+1
              ENDIF
              ENDIF
            ENDIF
            ENDIF
            ENDIF
C-- Y-direction
C??? DOES LMV NEED TO BE HALO-UPDATED???
            IF(LMV(I,J-1).ge.L) THEN
c           IF(LMM(I,J+1).NE.0) THEN
C??? Still need this one ^ ???
c             IF(J.LE.70.OR.J.GE.110) THEN
c             IF(J.LE.65) THEN
c             IF(I.LE.100.OR.I.GE.200) THEN
c             TRM(I,J  ,L) = TRM0(I,J  ,L) + flux_y(I,J,L)
              IF(ktap(I,J-1).gt.L) THEN
c             IF(MAXVAL(flux_y(I,J,:)).NE.0..AND.
c    *           MINVAL(flux_y(I,J,:)).NE.0.) THEN
              IF(SUM(flux_y(I,J,:)).NE.0.) THEN
c             IF(flux_y0(I,J).EQ.0) THEN
c             TRM(I,J  ,L) = TRM0(I,J  ,L) + flux_y(I,J,L)
c             TRM(I,J-1,L) = TRM0(I,J-1,L) - flux_y(I,J,L)
              TRM(I,J  ,L) = TRM(I,J  ,L) + flux_y(I,J,L)
              TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
c                IF(L.EQ.2) WRITE(*,*)'Y:yep',nstep,I,J
              ky=ky+1
c             ELSEIF(MAXVAL(flux_y(I,J,:)).EQ.0..AND.
c    *           MINVAL(flux_y(I,J,:)).EQ.0..AND.L.EQ.2) THEN
c                WRITE(*,*)'Y:nope',nstep,I,J
              ENDIF
              ENDIF
c             ENDIF
c             ENDIF
c           ENDIF
            ENDIF
C-- Z-direction
            IF(ktap(I,J).gt.L) THEN
c             IF(MAXVAL(flux_z(I,J,:)).NE.0..AND.
c    *           MINVAL(flux_z(I,J,:)).NE.0.) THEN
              IF(SUM(flux_z(I,J,:)).NE.0.) THEN
              TRM(I,J,L  ) = TRM(I,J,L  ) + flux_z(I,J,L)
              TRM(I,J,L+1) = TRM(I,J,L+1) - flux_z(I,J,L)
              kz=kz+1
c             ELSEIF(MAXVAL(flux_z(I,J,:)).EQ.0..AND.
c    *           MINVAL(flux_z(I,J,:)).EQ.0..AND.L.EQ.2) THEN
c                WRITE(*,*)'Z:',nstep,I,J
              ENDIF
            ENDIF
          ENDIF
          ENDDO
c210      CONTINUE
          IM1 = I
        ENDDO
        ENDIF
      ENDDO
c     WRITE(*,*)nstep,'kx=',kx,'ky=',ky,'kz=',kz

      J = J_1S + 1
      if (J.lt.JDM) then

C**** Non-Polar boxes
c     IM1 = IDM
      DO I=1,IDM
      DO L=1,ktap(I,J)

      IF(LMM(I,J).gt.0.AND.TRM(I,J,L).ne.0.) THEN

        IF(LMV(I,J-1).ge.L) THEN
c       IF(LMM(I,J+1).NE.0) THEN
        IF(ktap(I,J-1).gt.L) THEN
        IF(SUM(flux_y(I,J,:)).NE.0.) THEN
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
        END IF
c       END IF
        END IF
        END IF

      END IF

C**** END of L and I loops
      END DO
c     IM1 = I
      END DO

      endif

c     IF(HAVE_NORTH_POLE) THEN
C****   North Polar box
c       STRNP=0.
C****   Fluxes in Y-direction
c       DO I=1,IM
c         DO L=1,ktap(I,J)
c           IF(LMV(I,JM-1).ge.L) THEN
C??? .AND.ktap(I,J-1).GE.L ??? and also to be applied above ???
C****       Add and Subtract horizontal Y fluxes
c             STRNP(L)= STRNP(L) + flux_y(I,JM,L)
c             TRM(I,JM-1,L) = TRM0(I,JM-1,L) - flux_y(I,JM,L)
c           END IF
c         ENDDO
c       END DO
C**** adjust polar box
c       DO I=1,IM
c         DO L=1,ktap(I,J)
c           TRM(1,JM,L)=TRM0(1,JM,L) + STRNP(L)/IM
c         ENDDO
c       ENDDO
c     END IF


C??? ALSO LOOP at J=J1S+1 ???

c     DO J=j_0,j_1
c       DO I=i_0,i_1
c         DO L=1,ktap(I,J)
c           TRXMO(I,J,L)=TRXMO0(I,J,L)
c           TRYMO(I,J,L)=TRYMO0(I,J,L)
c         ENDDO
c       ENDDO
c     ENDDO

C**** Non-Polar boxes
c     DO L=1,KDM
c     DO J=J_0S,J_1S
c     IM1 = IDM
c     DO I=I_0,I_1

c     IF(LMM(I,J).le.0) GO TO 613

c     DO L=1,ktap(I,J)


c       IF(LMV(I,J-1).ge.L) THEN
c         TRM(I,J  ,L) = TRM(I,J  ,L) + flux_y(I,J,L)
c         TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
c       END IF

c     ENDDO

C**** END of I and J loops
c 613 IM1 = I
c     END DO
c     END DO
c     ENDDO

c     J = J_1S + 1
c     if (J.lt.JDM) then

c     DO L=1,KDM
C**** Non-Polar boxes
c     IM1 = IDM
c     DO I=1,IDM

c     IF(LMM(I,J).le.0) GO TO 614

c       IF(LMV(I,J-1).ge.L) THEN
c         TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
c       END IF

C**** END of L and I loops
c 614 IM1 = I
c     END DO
c     END DO

c     endif

c     CALL HALO_UPDATE(ogrid,
c    *                TRM(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
c    *                FROM=SOUTH)

      RETURN
      END SUBROUTINE MESO_D
#endif

      subroutine get_gradients(mokgm2,q_in,flag,qx,qy)
      use domain_decomp_1d, only :
     &     getDomainBounds,halo_update,south,north
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : dxpo,dyvo,dxypo
      use ocean, only : lmu,lmm,
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      use ocean, only : im,jm,lmo,ivnp,sinic,cosic,sinu,cosu
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mokgm2, ! units: kg/m2
     &     q_in, ! extensive or intesive units
     &     qx,qy ! outputs have intensive units
      integer flag ! 1: q_in is extensive; 0: q_in is intensive

      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     q,qxe,qyn
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
          q(i,j,l) = q_in(i,j,l)/(mokgm2(i,j,l)*dxypo(j))
        enddo
        enddo
        enddo
        enddo
      else               ! q_in is intensive
        q=q_in
      endif

      call halo_update(grid,q,from=north)

      if(have_north_pole) then
        do l=1,lmo
          q(2:im,jm,l) = q(1,jm,l)
        enddo
      endif

      ! compute gradients at cell edges.  gradients at coastlines
      ! are zero
      qxe = 0.
      qyn = 0.
      do l=1,lmo
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            IF(q(i+1,j,l).NE.0..AND.q(i,j,l).NE.0.) THEN
            qxe(i,j,l) = (q(i+1,j,l)-q(i,j,l))/dxpo(j)
            ENDIF
          enddo
          enddo
          i=im
          if(lmu(i,j).ge.l) then
            IF(q(1,j,l).NE.0..AND.q(i,j,l).NE.0.) THEN
            qxe(i,j,l) = (q(1,j,l)-q(i,j,l))/dxpo(j)
            ENDIF
          endif
          do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            IF(q(i,j+1,l).NE.0..AND.q(i,j,l).NE.0.) THEN
            qyn(i,j,l) = (q(i,j+1,l)-q(i,j,l))/dyvo(j)
            ENDIF
          enddo
          enddo
        enddo
      enddo

      ! average cell-edge gradients to cell centers.
      ! gradients in the north polar cap temporarily left at zero.
      call halo_update(grid,qyn,from=south)

      qx = 0.
      qy = 0.
      do l=1,lmo
      do j=j_0s,j_1s
        do n=1,nbyzm(j,l)
        do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
          qx(i,j,l) = .5*(qxe(i-1,j,l)+qxe(i,j,l))
          qy(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
        enddo
        enddo
        i=1
        if(lmm(i,j).ge.l) then
          qx(i,j,l) = .5*(qxe(im,j,l)+qxe(i,j,l))
          qy(i,j,l) = .5*(qyn(i,j-1,l)+qyn(i,j,l))
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
            qx(i,jm,l) = unp*cosu(i)  + vnp*sinu(i)
            qy(i,jm,l) = vnp*cosic(i) - unp*sinic(i)
          enddo
c         qx(im,jm,l) = unp   ! as a result of the above loop
c         qx(ivnp,jm,l) = vnp ! as a result of the above loop
        enddo
      endif
      return
      end subroutine get_gradients

C**********
      SUBROUTINE WKB(z,n2,km,frd,b1)
      IMPLICIT NONE
      INTEGER km,k,j
      REAL*8 z(km),n2(km),pi,ng(km),n2g(km),frd,n2gi(km),ngi(km),cm,
     *  Bsurf,h,phi0(km),tmpa(km),n(km),B,phi(km),db1w(km),b1w(km),tmp,
     *  zg(km),b1(km)

      pi=4.*atan(1.)

      DO k=1,km
c       n2(k)=MAX(0.D0,n2(k))
        n2(k)=MAX(1.d-8,n2(k))
      END DO

      n=dsqrt(n2)

      do k=1,km
        j=km+1-k
        zg(k)=-z(j)
        n2g(k)=n2(j)
        ng(k)=n(j)
      end do

      call dqtfg(zg,n2g,n2gi,km)
      call dqtfg(zg,ng,ngi,km)

      cm=ngi(km)/pi
      Bsurf=ng(km)**(-.5)*cm**(-1.5) ! mike
      h=-zg(1)

      do k=1,km ! going upward from -h to surf
        phi0(k)=(ng(k)/cm)**(-.5)*sin(ngi(k)/cm)
        tmpa(k)=(zg(k)+h)*n2g(k)*phi0(k)
      end do
      call dqtfg(zg,tmpa,tmpa,km)

      B=h/tmpa(km)

      do k=1,km
        phi(k)=B*phi0(k)
        db1w(k)=n2g(k)*phi(k)
        tmpa(k)=phi(k)**2
      end do

      call dqtfg(zg,db1w,b1w,km)
      !output int. of phi**2 from -h to 0:
      call dqtfg(zg,tmpa,tmpa,km)
c     write(*,'(a,14e14.6)') "int. phi**2 from -h to 0=",tmpa(m)
      !end output int. of phi**2 from -h to 0

      tmp=b1w(km)
      do k=1,km
        b1w(k)=b1w(k)-tmp+1.
      end do

      DO k=1,km
        b1(k)=b1w(km-k+1)
      ENDDO

      frd=cm

      RETURN
      END SUBROUTINE WKB

C************
#ifdef OCN_GISS_MESO
      SUBROUTINE MESO_A(TRM,TRXMO,TRYMO,TRZMO)
      USE IEEE_ARITHMETIC
      USE MODEL_COM,  only : nstep=>itime
      USE MESOML, only: ktap,ogrid,
     *  Unew,Vnew,Wnew,idm,kdm,dxpo,dypo,jdm
     * ,kappamsigmsx3d,kappamsigmsy3d,LX3d,LY3d,DH,zoc
     * ,kappamsigx3d,kappamsigy3d
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
     ., HALO_UPDATE, NORTH, SOUTH
      USE OCEAN, only : lmm,lmu,lmv,focean
      USE OCEAN, only : dxypo,mo,dypo
      USE OCEAN,      only : fluxA_x_sm,fluxA_y_sm,fluxA_z_sm,kpl
      USE GM_COM, ONLY: grid,GETDomainBounds

      IMPLICIT NONE

      INTEGER i,j,i_0,i_1,j_0,j_1,ip1,im1,k,l,it,jt,itp1,jtp1,it0,jt0
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     *  TRM,TRXMO,TRYMO,TRZMO,Fdiffx,Fdiffy,
     *  flux_x,flux_y,flux_z,
     *  Fnewx,Fnewy,dzFnewx,dzFnewy,
     *  DrhoTAUx,DrhoTAUy,dzTRM
      REAL*8 wta3m
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     .   zma,zh,zt

      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      CALL HALO_UPDATE(ogrid,
     *                ktap(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
     *                FROM=SOUTH+NORTH)

      CALL HALO_UPDATE(ogrid,
     *                TRM(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                FROM=SOUTH)

      CALL HALO_UPDATE(ogrid,
     *                LMM(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
     *                FROM=SOUTH+NORTH)

      CALL HALO_UPDATE(ogrid,
     *                LMV(:,ogrid%j_strt_halo:ogrid%j_stop_halo),
     *                FROM=SOUTH+NORTH)

      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP

C-- z at tracer level
      DO j=j_0,j_1
      DO i=i_0,i_1
        zh(i,j,1)=DH(i,j,1)/2.
        zma(i,j,1)=DH(i,j,1)
        DO k=2,kdm-1
          zh(i,j,k)=zma(i,j,k-1)+DH(i,j,k)/2.
          zma(i,j,k)=zma(i,j,k-1)+DH(i,j,k)
        ENDDO
      ENDDO
      ENDDO
      zh(:,:,kdm)=zoc(kdm)

      zt=-zh*100.
C--

c     do 104 j=j_0,j_1
c     do 104 i=i_0,i_1
c     IF(i.eq.42.AND.j.eq.127) THEN
c       WRITE(*,*)'DEBUG-Unew-in MESO_A',nstep,i,j
c       DO k=1,lmm(i,j)
c         WRITE(*,*)k,Unew(i,j,k),Vnew(i,j,k),Wnew(i,j,k)
c       ENDDO
c     ENDIF
c104  CONTINUE

C-- Calculcate Fnew
C--
      dzTRM=0.
      DO j=j_0,j_1
        DO i=i_0,i_1
        IF(FOCEAN(I,J).gt.0..AND.lmm(i,j).gt.1) THEN
          CALL d1sym(lmm(i,j),zt(i,j,1:lmm(i,j))/100.,
     *     TRM(i,j,1:lmm(i,j)),dzTRM(i,j,1:lmm(i,j)),11)
        ENDIF
        ENDDO
      ENDDO
C--- Calculate gradrho(tau)
c     DrhoTAUx = TRXMO + LX3D*TRZMO
c     DrhoTAUy = TRYMO + LY3D*TRZMO
      DrhoTAUx = TRXMO + LX3D*dzTRM
      DrhoTAUy = TRYMO + LY3D*dzTRM
C--
      Fnewx=kappamsigmsx3d*DrhoTAUx
      Fnewy=kappamsigmsy3d*DrhoTAUy
C--
      dzFnewx=0.; dzFnewy=0.
      DO j=j_0,j_1
        DO i=i_0,i_1
        IF(FOCEAN(I,J).gt.0..AND.lmm(i,j).gt.3) THEN
          DO k=1,50
            CALL z121_mesosc(Fnewx,lmm(i,j),lmm(i,j))
            CALL z121_mesosc(Fnewy,lmm(i,j),lmm(i,j))
          ENDDO
          CALL d1sym(lmm(i,j),zt(i,j,1:lmm(i,j))/100.,
     *     Fnewx(i,j,1:lmm(i,j)),dzFnewx(i,j,1:lmm(i,j)),12)
          CALL d1sym(lmm(i,j),zt(i,j,1:lmm(i,j))/100.,
     *     Fnewy(i,j,1:lmm(i,j)),dzFnewy(i,j,1:lmm(i,j)),13)
        ENDIF
        ENDDO
      ENDDO

C-- NEW COND!
      do 102 j=j_0,j_1
      do 102 i=i_0,i_1
        IF(SUM(kappamsigmsx3d(i,j,:)).EQ.0.
     * .OR.SUM(kappamsigx3d(i,j,:)).EQ.0.
     * .OR.SUM(Unew(i,j,:)).EQ.0.) THEN
          dzFnewx(i,j,:)=0.
        ENDIF
        IF(SUM(kappamsigmsy3d(i,j,:)).EQ.0.
     * .OR.SUM(kappamsigy3d(i,j,:)).EQ.0.
     * .OR.SUM(Vnew(i,j,:)).EQ.0.) THEN
          dzFnewy(i,j,:)=0.
        ENDIF
cc      IF(ktap(i,j).eq.1) THEN
cc        dxFdiffx(i,j,:)=0.
cc        dyFdiffy(i,j,:)=0.
cc      ENDIF
 102  CONTINUE

c     do 103 j=j_0,j_1
c     do 103 i=i_0,i_1
c     IF(i.eq.272.AND.j.eq.45) THEN
c       WRITE(*,*)'DEBUG-Fnew-in MESO_A',nstep,i,j,lmm(i,j)
c       DO k=1,lmm(i,j)
c         WRITE(*,*)k,LX3d(i,j,k)
c       ENDDO
c     ENDIF
c103  CONTINUE

C- fluxes
      flux_x=0.; flux_y=0.; flux_z=0.

      DO J=j_0,j_1
        IM1 = idm
        DO I=i_0,i_1
          IF(FOCEAN(I,J).GT.0.) THEN
          DO L=1,lmm(I,J)
            IF(LMU(IM1,J).ge.L) THEN
C-- NC!!
c             IF(Unew(I,J,L).NE.0..AND.Unew(IP1,J,L).NE.0.) THEN
C-- NC!!
              IF(LMM(IM1,J).NE.0) THEN
C-- NC!!
              IF(lmm(IM1,J).gt.L) THEN
              flux_x(I,J,L) = Unew(I,J,L)*TRXMO(I,J,L)+dzFnewx(I,J,L)
              ENDIF
              ENDIF
c             ENDIF
c             flux_x(I,J,L) = Unew(I,J,L)*TRXMO(I,J,L)
c             IF(.NOT.IEEE_IS_FINITE(flux_x(I,J,L)).OR.
c    *        IEEE_IS_NAN(flux_x(I,J,L)))
c    *         WRITE(*,*)'DEBUG-flux_xA:',nstep,I,J,L,flux_x(I,J,L),
c    *        Unew(I,J,L),dzFnewx(I,J,L),lmm(I,J),FOCEAN(I,J)
            ENDIF
            IF(LMV(I,J-1).ge.L) THEN
C-- NC!!
c             IF(LMM(I,J+1).NE.0) THEN
C-- NC!!
              IF(lmm(I,J-1).gt.L) THEN
              flux_y(I,J,L) = Vnew(I,J,L)*TRYMO(I,J,L)+dzFnewy(I,J,L)
              ENDIF
              ENDIF
c             flux_y(I,J,L) = Vnew(I,J,L)*TRYMO(I,J,L)
c             IF(.NOT.IEEE_IS_FINITE(flux_y(I,J,L)).OR.
c    *        IEEE_IS_NAN(flux_y(I,J,L))) 
c    *         WRITE(*,*)'DEBUG-flux_yA:',nstep,I,J,L,flux_y(I,J,L),
c    *        Vnew(I,J,L),dzFnewy(I,J,L),lmm(I,J),FOCEAN(I,J)
c           ENDIF
            IF(lmm(I,J).ge.L.AND.L.NE.1.AND.L.NE.LMM(I,J)) THEN
c             flux_z(I,J,L) = Wnew(I,J,L)*TRZMO(I,J,L)
              flux_z(I,J,L) = Wnew(I,J,L)*dzTRM(I,J,L)
c             IF(.NOT.IEEE_IS_FINITE(flux_z(I,J,L)).OR.
c    *        IEEE_IS_NAN(flux_z(I,J,L))) 
c    *         WRITE(*,*)'DEBUG-flux_zA:',nstep,I,J,L,flux_z(I,J,L),
c    *        Wnew(I,J,L),lmm(I,J),FOCEAN(I,J)
            ENDIF
          ENDDO
          IM1=I
          ENDIF
        ENDDO
      ENDDO


      wta3m=exp(-1./(3.*1440.))

      fluxA_x_sm=wta3m*fluxA_x_sm + (1.-wta3m)*flux_x
      fluxA_y_sm=wta3m*fluxA_y_sm + (1.-wta3m)*flux_y
      fluxA_z_sm=wta3m*fluxA_z_sm + (1.-wta3m)*flux_z

      flux_x=fluxA_x_sm
      flux_y=fluxA_y_sm
      flux_z=fluxA_z_sm

c     it0=42;jt0=127
c     do 104 j=j_0,j_1
c     IF(j.eq.jt0) THEN
c       WRITE(*,*)'DEBUG-fluxes-in MESO_A',nstep,it0,j,ktap(it0,j)
c       DO L=1,lmm(it0,j)
c         WRITE(*,114)L,flux_x(it0,j,L),flux_y(it0,j,L),flux_z(it0,j,L)
c       ENDDO
c       WRITE(*,*)'DEBUG-fluxez-in MESO_A',nstep,it0,j,ktap(it0,j)
c       DO L=1,lmm(it0,j)
c         WRITE(*,114)L,flux_z(it0,j,L),Wnew(it0,j,L),dzTRM(it0,j,L)
c       ENDDO
c       WRITE(*,*)'DEBUG-Unew-in MESO_A',nstep,it0,j
c       DO L=1,lmm(it0,j)
c         WRITE(*,114)k,Unew(it0,j,k),Vnew(it0,j,k),Wnew(it0,j,k)
c       ENDDO
c       WRITE(*,*)'DEBUG-dzFnewx-in MESO_A',nstep,it0,j
c       DO L=1,lmm(it0,j)
c         WRITE(*,114)k,dzFnewx(it0,j,k),dzFnewy(it0,j,k),zt(it0,j,L)
c       ENDDO
c     ENDIF
c104  CONTINUE
c114  FORMAT(I,3E)

      CALL HALO_UPDATE(ogrid,
     *                flux_y(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                FROM=SOUTH+NORTH)

      DO J=j_0,j_1
        IF(J.LT.JDM) THEN
        IM1 = idm
        DO I=i_0,i_1
        IF(ktap(I,J).gt.1) THEN
          DO L=min(ktap(I,J)+1,lmm(I,J)),lmm(I,J)
c         DO L=1,lmm(I,J)
          IF(LMM(I,J).gt.0.AND.TRM(I,J,L).ne.0.) THEN
C-- X-direction
            IF(LMU(IM1,J).ge.L) THEN
C-- NC!!
c             IF(Unew(I,J,L).NE.0..AND.Unew(IP1,J,L).NE.0.) THEN
C-- NC!!
c             IF(LMM(IM1,J).NE.0) THEN
C-- NC!!
c             IF(lmm(IM1,J).gt.L) THEN
C-- NC!!
c             IF(SUM(flux_x(I,J,:)).NE.0.) THEN
C--
              TRM(I  ,J,L) = TRM(I  ,J,L) + flux_x(I,J,L)
              TRM(IM1,J,L) = TRM(IM1,J,L) - flux_x(I,J,L)
C--
c             ENDIF
c             ENDIF
c             ENDIF
c             ENDIF
C--
            ENDIF
C-- Y-direction
            IF(LMV(I,J-1).ge.L) THEN
C-- NC!!
c             IF(LMM(I,J+1).NE.0) THEN
C-- NC!!
c             IF(lmm(I,J-1).gt.L) THEN
C-- NC!!
c             IF(SUM(flux_y(I,J,:)).NE.0.) THEN
C--
              TRM(I,J  ,L) = TRM(I,J  ,L) + flux_y(I,J,L)
              TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
C--
c             ENDIF
c             ENDIF
c             ENDIF
C--
            ENDIF
C-- Z-direction
            IF(lmm(I,J).gt.L) THEN
C-- NC!!
c             IF(SUM(flux_z(I,J,:)).NE.0.) THEN
C--
              TRM(I,J,L  ) = TRM(I,J,L  ) + flux_z(I,J,L)
              TRM(I,J,L+1) = TRM(I,J,L+1) - flux_z(I,J,L)
C--
c             ENDIF
C--
            ENDIF
          ENDIF
          ENDDO
          IM1 = I
        ENDIF
        ENDDO
        ENDIF
      ENDDO

      J = J_1S + 1
      if (J.lt.JDM) then

C**** Non-Polar boxes
c     IM1 = IDM
      DO I=1,IDM
      IF(ktap(I,J).gt.1) THEN
      DO L=min(ktap(I,J)+1,lmm(I,J)),lmm(I,J)

      IF(LMM(I,J).gt.0.AND.TRM(I,J,L).ne.0.) THEN

        IF(LMV(I,J-1).ge.L) THEN
          TRM(I,J-1,L) = TRM(I,J-1,L) - flux_y(I,J,L)
        END IF

      END IF

C**** END of L and I loops
c     IM1 = I
      END DO
      ENDIF
      END DO

      endif
      

      RETURN
      END SUBROUTINE MESO_A
#endif
