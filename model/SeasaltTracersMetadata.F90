!------------------------------------------------------------------------------
module SeasaltTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  SeasaltTracersMetadata_mod encapsulates the sea salt tracers metadata
!@auth NCCS ASTG, Kostas Tsigaridis
  use OldTracer_mod, only: oldAddTracer
  use OldTracer_mod, only: nPart, nGAS
  use OldTracer_mod, only: set_tr_mm
  use OldTracer_mod, only: set_ntm_power
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_tr_wd_type
  use OldTracer_mod, only: set_ntisurfsrc
  use TRACER_COM, only:  n_seasalt1,  n_seasalt2
  use TRACER_COM, only: offline_dms_ss, offline_ss
  use TRACER_COM, only: set_ntsurfsrc
  use Tracer_mod, only: Tracer

  implicit none
  private

  public Seasalt_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine Seasalt_InitMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer

    call  seasalt1_setSpec('seasalt1')
    call  seasalt2_setSpec('seasalt2')
     
!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    subroutine seasalt1_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt1 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 4.4d-7 ) ! This is non-hydrated
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine seasalt1_setSpec

    subroutine seasalt2_setSpec(name)
      character(len=*), intent(in) :: name
      n = oldAddTracer(name)
      n_seasalt2 = n
      call set_ntsurfsrc(n,  0) ! ocean bubbles
      call set_ntisurfsrc(n, 1)
      call set_ntm_power(n, -9)
      call set_tr_mm(n, 75.d0)  !Na x 3.256
      call set_trpdens(n, 2.2d3) !kg/m3 This is for non-hydrated
      call set_trradius(n, 5.0d-6) ! This is non-hydrated
      if (OFFLINE_DMS_SS.ne.1 .and. OFFLINE_SS.ne.1) then
        call set_trradius(n, 1.7d-6 ) ! This is non-hydrated
      end if
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
    end subroutine seasalt2_setSpec

  end subroutine Seasalt_InitMetadata

end module SeasaltTracersMetadata_mod



