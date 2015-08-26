#include "rundeck_opts.h"

      MODULE obio_forc

#ifdef OBIO_ON_GARYocean
      USE OCEANRES, only : kdm=>lmo 
#else
      USE hycom_dim, only: kdm
#endif
      use ocalbedo_mod, only: nlt


      implicit none


      integer, ALLOCATABLE, DIMENSION(:,:) :: ihra            !counter for daylight hours

      real, ALLOCATABLE, DIMENSION(:,:,:)  :: tirrq3d
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: avgq            !mean daily irradiance in quanta
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: atmFe
      real, ALLOCATABLE, DIMENSION(:,:,:)  :: alk             !alkalinity in 'umol/kg'

      real solz               !mean cosine solar zenith angle
      real sunz               !solar zenith angle
      real Ed(nlt),Es(nlt)
      real wind               !surface wind from atmos
      real tirrq(kdm)         !total mean irradiance in quanta
      real, parameter ::  tirrq_critical=10. !in quanta threshold at compensation depth
      real rmud               !downwelling irradiance average cosine
      real rhosrf             !surface air density which comes from PBL.f

      END MODULE obio_forc

!------------------------------------------------------------------------------
      subroutine alloc_obio_forc
      USE obio_forc
#ifdef OBIO_ON_GARYocean
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES, only : kdm=>lmo
#else
      USE hycom_dim, only : ogrid,i_0,i_1,j_0,j_1
#endif


      implicit none

#ifdef OBIO_ON_GARYocean
      INTEGER :: j_0,j_1,i_0,i_1

      I_0 = ogrid%I_STRT
      I_1 = ogrid%I_STOP
      J_0 = ogrid%J_STRT
      J_1 = ogrid%J_STOP
#endif

      ALLOCATE(tirrq3d(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(   ihra(i_0:i_1,j_0:j_1))
      ALLOCATE(   avgq(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(    alk(i_0:i_1,j_0:j_1,kdm))
      ALLOCATE(atmFe(i_0:i_1,j_0:j_1,12))

      end subroutine alloc_obio_forc
