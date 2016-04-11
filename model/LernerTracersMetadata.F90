!------------------------------------------------------------------------------
module LernerTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  LernerTracersMetadata_mod encapsulates the Lerner tracers metadata
!@auth NCCS ASTG
  use Dictionary_mod, only: sync_param
  use sharedTracersMetadata_mod
  use TRACER_COM, only: n_CH4, n_N2O, n_SF6, n_CO2, n_CFC11, n_14CO2, &
    n_O3, n_SF6_c, n_Rn222
  USE TRACERS_MPchem_COM, only: n_MPtable, tcscale
  use TRACER_COM, only: set_ntsurfsrc
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
  USE LINOZ_CHEM_COM, only: dsol
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: set_tr_mm, set_ntm_power
  use OldTracer_mod, only: set_t_qlimit
  use RunTimeControls_mod, only: tracers_special_lerner
  use Tracer_mod, only: Tracer
  implicit none

  private
  public Lerner_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Lerner_InitMetadata(pTracer, phase)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer
    integer , intent(in) :: phase

    if (phase==1) then
      call  N2O_setSpec('N2O')
      call  CH4_setSpec('CH4')
    else
      call  SF6_setSpec('SF6')
      call  Rn222_setSpec('Rn222')
      call  CO2_setSpec('CO2')
      call  CFC11_setSpec('CFC11')
      call  C_14O2_setSpec('14CO2')
      call  O3_setSpec('O3')
      call  SF6_c_setSpec('SF6_c')
    end if

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine SF6_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6 = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
    end subroutine SF6_setSpec

    subroutine CO2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CO2 = n
      call set_ntm_power(n, -6)
      call set_tr_mm(n, 44.d0)
      call set_t_qlimit(n,  .false.)
      call set_ntsurfsrc(n,  6)
    end subroutine CO2_setSpec

    subroutine CFC11_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_CFC11 = n
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 137.4d0)
      call set_ntsurfsrc(n,  1)
      if (tracers_special_lerner) then
        n_mptable(n) = 2
        tcscale(n_MPtable(n)) = 1.
      end if
    end subroutine CFC11_setSpec

    subroutine C_14O2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_14CO2 = n
      call set_ntm_power(n, -18)
      call set_tr_mm(n, 46.d0)
      call set_ntsurfsrc(n,  1)
    end subroutine C_14O2_setSpec

    subroutine O3_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_O3 = n
      call set_ntm_power(n, -8)
      call set_tr_mm(n, 48.d0)
      call set_ntsurfsrc(n,  1)
      if (tracers_special_lerner) then
      !**** Get solar variability coefficient from namelist if it exits
        dsol = 0.
        call sync_param("dsol",dsol)
      end if
    end subroutine O3_setSpec

    subroutine SF6_c_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_SF6_c = n
      call set_ntm_power(n, -14)
      call set_tr_mm(n, 146.01d0)
      call set_ntsurfsrc(n,  1)
    end subroutine SF6_c_setSpec

  end subroutine Lerner_InitMetadata

end module LernerTracersMetadata_mod



