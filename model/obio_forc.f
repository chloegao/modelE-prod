#include "rundeck_opts.h"

      MODULE obio_forc

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo 
#else
      USE hycom_dim_glob
#endif
      USE obio_dim
      use ocalbedo_mod, only: nlt


      implicit none


      integer, ALLOCATABLE, DIMENSION(:,:) :: ihra            !counter for daylight hours

      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe_glob      !surface iron deposition
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity in 'umol/kg'
#ifdef TRACERS_Alkalinity
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk_glob        !alkalinity in 'umol/kg'
#endif

#ifndef OBIO_RAD_coupling
      real, ALLOCATABLE, DIMENSION(:,:,:,:,:):: Eda,Esa       !direct,diffuse downwelling irradiance
#endif

      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
#ifdef OBIO_RAD_coupling 
      real eda_frac(nlt),esa_frac(nlt)
      real ovisdir_ij,ovisdif_ij,onirdir_ij,onirdif_ij
#else
      real Eda2(nlt,nhn),Esa2(nlt,nhn)
#endif
      real Ed(nlt),Es(nlt)
      real wind               !surface wind from atmos
      real tirrq(kdm)         !total mean irradiance in quanta
      real, parameter ::  tirrq_critical=10. !in quanta threshold at compensation depth
      real rmud               !downwelling irradiance average cosine
      real atmCO2
      real rhosrf             !surface air density which comes from PBL.f

      END MODULE obio_forc

!------------------------------------------------------------------------------
      subroutine alloc_obio_forc
      USE obio_forc
      USE obio_dim
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only : idm=>imo,jdm=>jmo,kdm=>lmo
#else
      USE hycom_dim_glob
      USE hycom_dim, only : ogrid,i_0h,i_1h,j_0h,j_1h
#endif


      implicit none

#ifdef OBIO_ON_GARYocean
      INTEGER :: j_0h,j_1h,i_0h,i_1h

      I_0H = ogrid%I_STRT_HALO
      I_1H = ogrid%I_STOP_HALO
      J_0H = ogrid%J_STRT_HALO
      J_1H = ogrid%J_STOP_HALO
#endif

      ALLOCATE(tirrq3d(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(   ihra(i_0h:i_1h,j_0h:j_1h))
      ALLOCATE(   avgq(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(    alk(i_0h:i_1h,j_0h:j_1h,kdm))
      ALLOCATE(atmFe(i_0h:i_1h,j_0h:j_1h,12),atmFe_glob(idm,jdm,12))
#ifdef TRACERS_Alkalinity
      ALLOCATE(alk_glob(idm,jdm,kdm))
#endif

      end subroutine alloc_obio_forc

      subroutine obio_forc_init
      use obio_forc, only : atmco2
      use dictionary_mod
      implicit none
#ifdef constCO2
      call get_param("atmCO2",atmCO2)   !need to do this here also
#ifdef OBIO_ON_GARYocean
      print*, 'OCNDYN, atmco2=',atmCO2
#else
      print*, 'OCEAN_hycom, atmco2=',atmCO2
#endif
#else
      atmCO2=0.  !progn. atmCO2, set here to zero, dummy anyway
#endif
      end subroutine obio_forc_init
