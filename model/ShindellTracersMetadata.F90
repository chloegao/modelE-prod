#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module ShindellTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  ShindellTracersMetadata_mod encapsulates the TRACERS_SPECIAL_Shindell
!@+    metadata.
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: CH4_setspec, &
    N2O_setspec, H2O2_setspec
  use sharedTracersMetadata_mod, only: convert_HSTAR
  use TRACER_COM, only: NTM_chem_beg, NTM_chem_end
  use TRACER_COM, only: n_CH4,  n_N2O, n_Ox,   n_NOx, & 
    n_N2O5,   n_HNO3,  n_H2O2,  n_CH3OOH,   n_HCHO,  &
    n_HO2NO2, n_CO,    n_PAN,   n_H2O17,             &
    n_Isoprene, n_AlkylNit, n_Alkenes, n_Paraffin,   &
    n_stratOx, n_Terpenes,n_codirect,                &
    n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a,         &
    n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a,         &
    n_ClOx,   n_BrOx,  n_HCl,   n_HOCl,   n_ClONO2,  &
    n_HBr,    n_HOBr,  n_BrONO2,n_CFC,    n_GLT
#ifdef TRACERS_AEROSOLS_SOA
  USE TRACERS_SOA, only: n_soa_i, n_soa_e
#endif
  use OldTracer_mod, only: nPart
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_F0
  use OldTracer_mod, only: set_tr_RKD
  use OldTracer_mod, only: set_tr_DHD
  use OldTracer_mod, only: tr_RKD 
  use OldTracer_mod, only: set_trdecay
  use OldTracer_mod, only: dodrydep
  use OldTracer_mod, only: F0
  use OldTracer_mod, only: HSTAR
  use OldTracer_mod, only: ngas, nPART
  use OldTracer_mod, only: set_emisPerFireByVegType
  use RunTimeControls_mod, only: tracers_special_shindell
  use RunTimeControls_mod, only: tracers_drydep
  use RunTimeControls_mod, only: tracers_terp
  use RunTimeControls_mod, only: tracers_aerosols_soa
  use RunTimeControls_mod, only: shindell_strat_extra
  use RunTimeControls_mod, only: accmip_like_diags
  USE CONSTANT, only: mair
#ifdef TRACERS_AEROSOLS_SOA
  USE CONSTANT, only: gasc
#endif
  use Tracer_mod, only: Tracer

  implicit none
  private

  public SHINDELL_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine SHINDELL_initMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  Ox_setSpec('Ox')
    call  NOx_setSpec('NOx')
    call  ClOx_setSpec('ClOx')
    call  BrOx_setSpec('BrOx')
    call  N2O5_setSpec('N2O5')
    call  HNO3_setSpec('HNO3')
    call  H2O2_setSpec('H2O2')
    call  CH3OOH_setSpec('CH3OOH')

    call  HCHO_setSpec('HCHO')
    call  HO2NO2_setSpec('HO2NO2')
    call  CO_setSpec('CO')
    call  CH4_setSpec('CH4')
    call  PAN_setSpec('PAN')
    call  Isoprene_setSpec('Isoprene')
    call  AlkylNit_setSpec('AlkylNit')
    call  Alkenes_setSpec('Alkenes')
    call  Paraffin_setSpec('Paraffin')

    if (tracers_terp) then
      call  Terpenes_setSpec('Terpenes')
    end if

#ifdef TRACERS_AEROSOLS_SOA
    if (tracers_aerosols_soa) then
      call  isopp1g_setSpec('isopp1g')
      call  isopp1a_setSpec('isopp1a')
      call  isopp2g_setSpec('isopp2g')
      call  isopp2a_setSpec('isopp2a')
      if (tracers_terp) then
        call  apinp1g_setSpec('apinp1g')
        call  apinp1a_setSpec('apinp1a')
        call  apinp2g_setSpec('apinp2g')
        call  apinp2a_setSpec('apinp2a')
      end if
    end if
#endif

    call  HCl_setSpec('HCl')
    call  HOCl_setSpec('HOCl')
    call  ClONO2_setSpec('ClONO2')
    call  HBr_setSpec('HBr')
    call  HOBr_setSpec('HOBr')
    call  BrONO2_setSpec('BrONO2')
    call  N2O_setSpec('N2O')
    call  CFC_setSpec('CFC')

    if (shindell_strat_extra) then
      if (accmip_like_diags) then
        call  codirect_setSpec('codirect')
        call  stratOx_setSpec('stratOx')
        call  GLT_setSpec('GLT') ! generic linear tracer
      end if
    end if

    call calculateIndexOffsets

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine calculateIndexOffsets
      use TRACER_COM, only: nn_CH4,  nn_N2O, nn_Ox,   nn_NOx, & 
           nn_N2O5,   nn_HNO3,  nn_H2O2,  nn_CH3OOH,   nn_HCHO,  &
           nn_HO2NO2, nn_CO,    nn_PAN,   nn_H2O17,             &
           nn_Isoprene, nn_AlkylNit, nn_Alkenes, nn_Paraffin,   &
           nn_stratOx, nn_Terpenes,nn_codirect,                &
           nn_isopp1g,nn_isopp1a,nn_isopp2g,nn_isopp2a,         &
           nn_apinp1g,nn_apinp1a,nn_apinp2g,nn_apinp2a,         &
           nn_ClOx,   nn_BrOx,  nn_HCl,   nn_HOCl,   nn_ClONO2,  &
           nn_HBr,    nn_HOBr,  nn_BrONO2,nn_CFC,    nn_GLT
      use TRACER_COM, only: NTM_chem_beg
      integer :: offset

     offset = ntm_chem_beg - 1
     nn_CH4 = n_CH4 - offset
     nn_N2O = n_N2O - offset
     nn_Ox = n_Ox - offset
     nn_NOx = n_NOx - offset
     nn_N2O5 = n_N2O5 - offset
     nn_HNO3 = n_HNO3 - offset
     nn_H2O2 = n_H2O2 - offset
     nn_CH3OOH = n_CH3OOH - offset
     nn_HCHO = n_HCHO - offset
     nn_HO2NO2 = n_HO2NO2 - offset
     nn_CO = n_CO - offset
     nn_PAN = n_PAN - offset
     nn_H2O17 = n_H2O17 - offset
     nn_Isoprene = n_Isoprene - offset
     nn_AlkylNit = n_AlkylNit - offset
     nn_Alkenes = n_Alkenes - offset
     nn_Paraffin = n_Paraffin - offset
     nn_stratOx = n_stratOx - offset
    if (tracers_terp) then
       nn_Terpenes = n_Terpenes - offset
     end if
     nn_codirect = n_codirect - offset
#ifdef TRACERS_AEROSOLS_SOA
     nn_isopp1g = n_isopp1g - offset
     nn_isopp1a = n_isopp1a - offset
     nn_isopp2g = n_isopp2g - offset
     nn_isopp2a = n_isopp2a - offset
     nn_apinp1g = n_apinp1g - offset
     nn_apinp1a = n_apinp1a - offset
     nn_apinp2g = n_apinp2g - offset
     nn_apinp2a = n_apinp2a - offset
#endif
     nn_ClOx = n_ClOx - offset
     nn_BrOx = n_BrOx - offset
     nn_HCl = n_HCl - offset
     nn_HOCl = n_HOCl - offset
     nn_ClONO2 = n_ClONO2 - offset
     nn_HBr = n_HBr - offset
     nn_HOBr = n_HOBr - offset
     nn_BrONO2 = n_BrONO2 - offset
     nn_CFC = n_CFC - offset
     nn_GLT = n_GLT - offset

    end subroutine calculateIndexOffsets

    subroutine Ox_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Ox = n
      NTM_chem_beg = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
      if (tracers_drydep) then
        call set_F0(n,  1.4d0)
        call set_HSTAR(n,  1.d-2)
      end if
    end subroutine Ox_setSpec

    subroutine NOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_NOx = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 14.01d0)
      if (tracers_drydep) then
        call set_F0(n,  1.d-1)
        call set_HSTAR(n,  1.d-2)
      end if
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.1378230d-07, &
      &  3.2166037d-07, 1.5559274d-07, 4.1611088d-07, 5.7316458d-07, &
      &  2.1700112d-07, 3.0054335d-07, 0.0000000d+00, 0.0000000d+00, &
      &  0.0000000d+00, 0.0000000d+00/)
#endif
      call check_aircraft_sectors(name) ! special 3D source case
    end subroutine NOx_setSpec

    subroutine ClOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_ClOx = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 51.5d0)
    end subroutine ClOx_setSpec

    subroutine BrOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BrOx = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 95.9d0)
    end subroutine BrOx_setSpec

    subroutine N2O5_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_N2O5 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 108.02d0)
    end subroutine N2O5_setSpec

    subroutine HNO3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HNO3 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 63.018d0)
      call set_tr_RKD(n, 2.073d3 ) ! in mole/J = 2.1d5 mole/(L atm)
      if (tracers_drydep) call set_HSTAR(n, 1.d14)
    end subroutine HNO3_setSpec

    subroutine CH3OOH_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CH3OOH = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 48.042d0)
      if (tracers_drydep) call set_HSTAR(n,  3.d2)
    end subroutine CH3OOH_setSpec

    subroutine HCHO_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HCHO = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 30.026d0)
      call set_tr_RKD(n, 6.218d1 ) ! mole/J = 6.3d3 mole/(L atm)
      if (tracers_drydep) call set_HSTAR(n, 6.d3)
    end subroutine HCHO_setSpec

    subroutine HO2NO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HO2NO2 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 79.018d0)
    end subroutine HO2NO2_setSpec

    subroutine CO_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CO = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 28.01d0)
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 7.0401156d-06, &
      & 1.8708386d-05, 1.0678024d-05, 2.6742857d-05, 4.0226296d-05,&
      & 2.3661527d-05, 4.4639346d-05, 0.0000000d+00, 0.0000000d+00,&
      & 0.0000000d+00, 0.0000000d+00/)
#endif
    end subroutine CO_setSpec
    subroutine PAN_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_PAN = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 121.054d0) ! assuming CH3COOONO2 = PAN)
      if (tracers_drydep) call set_HSTAR(n,  3.6d0)
    end subroutine PAN_setSpec

    subroutine Isoprene_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Isoprene = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 60.05d0) ! i.e. 5 carbons
      if (tracers_drydep) call set_HSTAR(n,  1.3d-2)
    end subroutine Isoprene_setSpec

    subroutine AlkylNit_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_AlkylNit = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)   !unknown molecular weight, so use air and make
      ! note in the diagnostics write-out...
    end subroutine AlkylNit_setSpec

    subroutine Alkenes_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Alkenes = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 1.0d0)  ! So, careful: source files now in Kmole/m2/s or
      ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 6.1516259d-09, &
      & 1.1544214d-08, 6.9501711d-09, 1.7481154d-08, 2.5840087d-08, &
      & 1.5709551d-08, 3.5913079d-08, 0.0000000d+00, 0.0000000d+00, &
      & 0.0000000d+00, 0.0000000d+00/)
#endif
    end subroutine Alkenes_setSpec

    subroutine Paraffin_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Paraffin = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 1.0d0)  ! So, careful: source files now in Kmole/m2/s or
      ! equivalently, kg/m2/s for species with tr_mm=1
#ifdef DYNAMIC_BIOMASS_BURNING
      ! 12 below are the 12 VDATA veg types or Ent remapped to them,
      ! from Olga Pechony's AR5_EPFC_factors_incl_SO2.xlsx file.
      emisPerFireByVegType(n,1:12)=(/0.0000000d+00, 1.5258348d-09, &
      & 5.6236904d-09, 3.1752858d-09, 1.0662656d-08, 1.5271524d-08, &
      & 8.0735774d-09, 2.6055675d-08, 0.0000000d+00, 0.0000000d+00, &
      & 0.0000000d+00, 0.0000000d+00/)
#endif
    end subroutine Paraffin_setSpec

    subroutine Terpenes_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_Terpenes = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
      if (tracers_drydep) call set_HSTAR(n,  1.3d-2)
    end subroutine Terpenes_setSpec

#ifdef TRACERS_AEROSOLS_SOA
    subroutine isopp1g_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_isopp1g = n
      n_soa_i = n_isopp1g       !the first from the soa species
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    end subroutine isopp1g_setSpec

    subroutine isopp1a_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_isopp1a = n
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
    end subroutine isopp1a_setSpec

    subroutine isopp2g_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_isopp2g = n
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    end subroutine isopp2g_setSpec

    subroutine isopp2a_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_isopp2a = n
      if (.not. tracers_terp) n_soa_e = n_isopp2a       !the last from the soa species
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
    end subroutine isopp2a_setSpec

    subroutine apinp1g_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_apinp1g = n
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    end subroutine apinp1g_setSpec

    subroutine apinp1a_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_apinp1a = n
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
    end subroutine apinp1a_setSpec

    subroutine apinp2g_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_apinp2g = n
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_tr_RKD(n, 1.d4 / convert_HSTAR ) !Henry; from mole/(L atm) to mole/J
      call set_tr_DHD(n, -12.d0 * gasc        ) !Henry temp dependence (J/mole), Chung and Seinfeld, 2002
      call set_tr_wd_type(n, ngas)
      if (tracers_drydep) call set_HSTAR(n, tr_RKD(n)*convert_HSTAR)
    end subroutine apinp2g_setSpec

    subroutine apinp2a_setSpec(name)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      real*8 :: tmp
      n = oldAddTracer(name)
      n_apinp2a = n
      n_soa_e = n_apinp2a       !the last from the soa species
      tmp = om2oc(n)
      call sync_param(trim(name)//"_om2oc",tmp)
      call set_om2oc(n, tmp)
      call set_ntm_power(n, -11)
      tmp = 12.d0 * om2oc(n)
      call set_tr_mm(n, tmp)
      call set_trpdens(n, 1.5d3) !kg/m3
      call set_trradius(n, 3.d-7) !m
      call set_fq_aer(n, 0.8d0) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, nPART)
    end subroutine apinp2a_setSpec
#endif  /* TRACERS_AEROSOLS_SOA */

    subroutine HCl_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HCl = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 36.5d0)
    end subroutine HCl_setSpec

    subroutine HOCl_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HOCl = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 52.5d0)
    end subroutine HOCl_setSpec

    subroutine ClONO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_ClONO2 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 97.5d0)
    end subroutine ClONO2_setSpec

    subroutine HBr_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 80.9d0)
    end subroutine HBr_setSpec

    subroutine HOBr_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_HOBr = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 96.9d0)
    end subroutine HOBr_setSpec

    subroutine BrONO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_BrONO2 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 141.9d0)
    end subroutine BrONO2_setSpec

    subroutine CFC_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CFC = n
      NTM_chem_end = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0) !CFC11
    end subroutine CFC_setSpec

    subroutine codirect_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_codirect = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 28.01d0)
      call set_trdecay(n,  2.31482d-7) ! 1/(50 days)
      ! not a radiactive decay, but functionally identical
    end subroutine codirect_setSpec

    subroutine stratOx_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_stratOx = n
      ! assumes initial Ox conditions read in for Ox tracer
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
      if (tracers_drydep) then
        call set_F0(n,  1.4d0)
        call set_HSTAR(n,  1.d-2)
      end if
    end subroutine stratOx_setSpec

    subroutine GLT_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_GLT = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, mair)
    end subroutine GLT_setSpec

  end subroutine SHINDELL_initMetadata

end Module ShindellTracersMetadata_mod

