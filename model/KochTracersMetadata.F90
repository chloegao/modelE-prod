#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module KochTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  KochTracersMetadata_mod encapsulates the KOCH tracers metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: DMS_setspec, &
    SO2_setspec, H2O2_s_setspec
  use sharedTracersMetadata_mod, only: convert_HSTAR
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: nPart, nGAS
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_emisPerFireByVegType
  use OldTracer_mod, only: set_pm2p5fact
  use OldTracer_mod, only: set_pm10fact
  use OldTracer_mod, only: set_has_chemistry
  use OldTracer_mod, only: tr_RKD 
  use TRACER_COM, only:  n_MSA, n_SO4, n_DMS, &
    n_BCII,  n_BCIA,  n_BCB, n_OCII,  n_OCIA,  n_OCB, n_H2O2_s
  use TRACER_COM, only: whichEPFCs
  use TRACER_COM, only: tracers
  use Dictionary_mod, only: sync_param
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: sulf_only_aerosols
  use RunTimeControls_mod, only: tracers_special_shindell
  use RunTimeControls_mod, only: dynamic_biomass_burning
  use Tracer_mod, only: Tracer

  implicit none
  private

  public Koch_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine KOCH_InitMetadata(pTracer)
    use TRACER_COM, only: coupled_chem
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  DMS_setSpec('DMS')
    call  MSA_setSpec('MSA')
    call  SO2_setSpec('SO2')
    call  SO4_setSpec('SO4')
    if (.not. tracers_special_shindell .or. coupled_chem.eq.0) then
      call  H2O2_s_setSpec('H2O2_s')
    end if
    if (.not. sulf_only_aerosols) then
      call  BCII_setSpec('BCII')
      call  BCIA_setSpec('BCIA')
      call  BCB_setSpec('BCB')
#ifndef TRACERS_AEROSOLS_VBS
      call  OCII_setSpec('OCII')   !Insoluble industrial organic mass
      call  OCIA_setSpec('OCIA')   !Aged industrial organic mass
      call  OCB_setSpec('OCB')     !Biomass organic mass
#endif /* not TRACERS_AEROSOLS_VBS */
    end if

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine MSA_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_MSA = n
      call set_ntm_power(n, -13)
      call set_tr_mm(n, 96.d+0) !(H2O2 34;SO2 64)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 5.d-7 ) !m (SO4 3;BC 1;OC 3)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine MSA_setSpec

    subroutine SO4_setSpec(name)
      character(len=*), intent(in) :: name
      type (Tracer), pointer :: t
      n = oldAddTracer(name)
      n_SO4 = n 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, 1.7d3) !kg/m3 this is sulfate value
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)

      t => tracers%getReference(trim(name))
      call t%insert('SO4',.true.)

    end subroutine SO4_setSpec

    subroutine BCII_setSpec(name)
      character(len=*), intent(in) :: name
      type (Tracer), pointer :: t
      n = oldAddTracer(name)
      n_BCII = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)

      t => tracers%getReference(trim(name))
      call t%insert('BC',.true.)
      
    end subroutine BCII_setSpec

    subroutine BCIA_setSpec(name)
      character(len=*), intent(in) :: name
      type (Tracer), pointer :: t
      n = oldAddTracer(name)
      n_BCIA = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)

      t => tracers%getReference(trim(name))
      call t%insert('BC',.true.)

    end subroutine BCIA_setSpec

    subroutine BCB_setSpec(name)
      character(len=*), intent(in) :: name
      type (Tracer), pointer :: t
      n = oldAddTracer(name)
      n_BCB = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d0)
      call set_trpdens(n, 1.3d3) !kg/m3
      call set_trradius(n, 1.d-7 ) !m
      call set_fq_aer(n, 0.8d0 ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      t => tracers%getReference(trim(name))
      call t%insert('BC',.true.)
      
#ifdef DYNAMIC_BIOMASS_BURNING
      if (dynamic_biomass_burning) then
        ! 12 below are the 12 VDATA veg types or Ent remapped to them,
        ! from Olga Pechony's EPFC.xlsx e-mailed to Greg 1/13/2013
        call sync_param("whichEPFCs",whichEPFCs)
        select case(whichEPFCs)
        case(1) ! AR5
          call set_emisPerFireByVegType(n, [0.d0,4.80d-9,8.33d-8,4.13d-8, &
          & 1.16d-7,1.08d-7,5.84d-8,6.10d-8,0.d0,0.d0,0.d0,0.d0] )
        case(2) ! GFED3
          call set_emisPerFireByVegType(n, [0.d0,1.64d-7,6.21d-8,3.04d-8, &
          & 3.01d-8,4.20d-8,6.87d-8,7.21d-8,0.d0,0.d0,0.d0,0.d0] )
        case(3) ! GFED2
          call set_emisPerFireByVegType(n, [0.d0,6.07d-8,3.98d-8,4.53d-8, &
          & 4.80d-8,4.78d-8,5.19d-8,7.95d-8,0.d0,0.d0,0.d0,0.d0] )
        case(4) ! MOPITT
          call set_emisPerFireByVegType(n, [0.d0,2.77d-8,1.12d-7,2.58d-8, &
          & 8.94d-8,7.49d-8,5.88d-9,2.50d-8,0.d0,0.d0,0.d0,0.d0] )
        case default
          call stop_model('whichEPFCs unknown',255)
        end select
      end if
#endif
    end subroutine BCB_setSpec

    subroutine OCII_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCII = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCII_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine OCII_setSpec

    subroutine OCIA_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCIA = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCIA_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 1.d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
      call set_has_chemistry(n, .true.)
    end subroutine OCIA_setSpec

    subroutine OCB_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_OCB = n
      call set_om2oc(n, 1.4d0)
      tmp = om2oc(n)
      call sync_param("OCB_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7 ) !m
      call set_fq_aer(n, 0.8d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_pm2p5fact(n, 1.d0) ! fraction that's PM2.5
      call set_pm10fact(n, 1.d0) ! fraction that's PM10
#ifdef DYNAMIC_BIOMASS_BURNING
      if (dynamic_biomass_burning) then
        ! 12 below are the 12 VDATA veg types or Ent remapped to them,
        ! from Olga Pechony's EPFC.xlsx e-mailed to Greg 1/13/2013
        call sync_param("whichEPFCs",whichEPFCs)
        select case(whichEPFCs)
        case(1) ! AR5
          call set_emisPerFireByVegType(n, [0.d0,5.71d-7,7.11d-7,4.02d-7, &
          & 1.18d-6,1.25d-6,9.51d-7,9.80d-7,0.d0,0.d0,0.d0,0.d0] )
        case(2) ! GFED3
          call set_emisPerFireByVegType(n, [0.d0,3.37d-6,9.98d-7,1.40d-7, &
          & 2.51d-7,6.44d-7,1.34d-6,6.92d-7,0.d0,0.d0,0.d0,0.d0] )
        case(3) ! GFED2
          call set_emisPerFireByVegType(n, [0.d0,8.74d-7,3.15d-7,4.02d-7, &
          & 4.70d-7,5.40d-7,9.81d-7,9.18d-7,0.d0,0.d0,0.d0,0.d0] )
        case(4) ! MOPITT
          call set_emisPerFireByVegType(n, [0.d0,3.51d-7,1.13d-6,2.37d-7, &
          & 7.18d-7,1.05d-6,1.18d-7,2.86d-7,0.d0,0.d0,0.d0,0.d0] )
        case default
          call stop_model('whichEPFCs unknown',255)
        end select
      end if
#endif
    end subroutine OCB_setSpec

  end subroutine KOCH_InitMetadata

end module KochTracersMetadata_mod



