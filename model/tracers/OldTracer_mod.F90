#include "rundeck_opts.h"




























module OldTracer_mod
  use TracerBundle_mod
  use Tracer_mod
  implicit none
  private

  public :: OldTracer_type
  public :: initializeOldTracers
  public :: oldAddTracer
  public :: findtracer
!!$  public :: internalTracers
  public :: trName
  public :: MAX_LEN_NAME

  public :: nGAS, nPART, nWATER

  public :: set_tr_mm, tr_mm 
public :: set_ntm_power, ntm_power 
public :: set_t_qlimit, t_qlimit 
public :: set_needtrs, needtrs 
public :: set_trdecay, trdecay 
public :: set_itime_tr0, itime_tr0 
public :: set_trsi0, trsi0 
public :: set_trw0, trw0 
public :: set_mass2vol, mass2vol 
public :: set_vol2mass, vol2mass 
public :: set_dodrydep, dodrydep 
public :: set_F0, F0 
public :: set_HSTAR, HSTAR 
public :: set_do_fire, do_fire 
public :: set_nBBsources, nBBsources 
public :: set_emisPerFireByVegType, emisPerFireByVegType 
public :: set_trpdens, trpdens 
public :: set_trradius, trradius 
public :: set_tr_wd_TYPE, tr_wd_TYPE 
public :: set_tr_RKD, tr_RKD 
public :: set_tr_DHD, tr_DHD 
public :: set_fq_aer, fq_aer 
public :: set_rc_washt, rc_washt 
public :: set_isDust, isDust 
public :: set_tr_H2ObyCH4, tr_H2ObyCH4 
public :: set_dowetdep, dowetdep 
public :: set_ntrocn, ntrocn 
public :: set_conc_from_fw, conc_from_fw 
public :: set_trglac, trglac 
public :: set_ntisurfsrc, ntisurfsrc 
public :: set_iso_index, iso_index 
public :: set_om2oc, om2oc 
public :: set_to_volume_MixRat, to_volume_MixRat 
public :: set_to_conc, to_conc 
public :: set_TRLI0, TRLI0 
public :: set_src_dist_base, src_dist_base 
public :: set_src_dist_index, src_dist_index 

  interface tr_mm
module procedure tr_mm_s
module procedure tr_mm_all
   module procedure tr_mm_m
end interface
interface ntm_power
module procedure ntm_power_s
module procedure ntm_power_all
   module procedure ntm_power_m
end interface
interface t_qlimit
module procedure t_qlimit_s
module procedure t_qlimit_all
   module procedure t_qlimit_m
end interface
interface needtrs
module procedure needtrs_s
module procedure needtrs_all
   module procedure needtrs_m
end interface
interface trdecay
module procedure trdecay_s
module procedure trdecay_all
   module procedure trdecay_m
end interface
interface itime_tr0
module procedure itime_tr0_s
module procedure itime_tr0_all
   module procedure itime_tr0_m
end interface
interface trsi0
module procedure trsi0_s
module procedure trsi0_all
   module procedure trsi0_m
end interface
interface trw0
module procedure trw0_s
module procedure trw0_all
   module procedure trw0_m
end interface
interface mass2vol
module procedure mass2vol_s
module procedure mass2vol_all
   module procedure mass2vol_m
end interface
interface vol2mass
module procedure vol2mass_s
module procedure vol2mass_all
   module procedure vol2mass_m
end interface
interface dodrydep
module procedure dodrydep_s
module procedure dodrydep_all
   module procedure dodrydep_m
end interface
interface F0
module procedure F0_s
module procedure F0_all
   module procedure F0_m
end interface
interface HSTAR
module procedure HSTAR_s
module procedure HSTAR_all
   module procedure HSTAR_m
end interface
interface do_fire
module procedure do_fire_s
module procedure do_fire_all
   module procedure do_fire_m
end interface
interface nBBsources
module procedure nBBsources_s
module procedure nBBsources_all
   module procedure nBBsources_m
end interface
interface emisPerFireByVegType
module procedure emisPerFireByVegType_s

end interface
interface trpdens
module procedure trpdens_s
module procedure trpdens_all
   module procedure trpdens_m
end interface
interface trradius
module procedure trradius_s
module procedure trradius_all
   module procedure trradius_m
end interface
interface tr_wd_TYPE
module procedure tr_wd_TYPE_s
module procedure tr_wd_TYPE_all
   module procedure tr_wd_TYPE_m
end interface
interface tr_RKD
module procedure tr_RKD_s
module procedure tr_RKD_all
   module procedure tr_RKD_m
end interface
interface tr_DHD
module procedure tr_DHD_s
module procedure tr_DHD_all
   module procedure tr_DHD_m
end interface
interface fq_aer
module procedure fq_aer_s
module procedure fq_aer_all
   module procedure fq_aer_m
end interface
interface rc_washt
module procedure rc_washt_s
module procedure rc_washt_all
   module procedure rc_washt_m
end interface
interface isDust
module procedure isDust_s
module procedure isDust_all
   module procedure isDust_m
end interface
interface tr_H2ObyCH4
module procedure tr_H2ObyCH4_s
module procedure tr_H2ObyCH4_all
   module procedure tr_H2ObyCH4_m
end interface
interface dowetdep
module procedure dowetdep_s
module procedure dowetdep_all
   module procedure dowetdep_m
end interface
interface ntrocn
module procedure ntrocn_s
module procedure ntrocn_all
   module procedure ntrocn_m
end interface
interface conc_from_fw
module procedure conc_from_fw_s
module procedure conc_from_fw_all
   module procedure conc_from_fw_m
end interface
interface trglac
module procedure trglac_s
module procedure trglac_all
   module procedure trglac_m
end interface
interface ntisurfsrc
module procedure ntisurfsrc_s
module procedure ntisurfsrc_all
   module procedure ntisurfsrc_m
end interface
interface iso_index
module procedure iso_index_s
module procedure iso_index_all
   module procedure iso_index_m
end interface
interface om2oc
module procedure om2oc_s
module procedure om2oc_all
   module procedure om2oc_m
end interface
interface to_volume_MixRat
module procedure to_volume_MixRat_s
module procedure to_volume_MixRat_all
   module procedure to_volume_MixRat_m
end interface
interface to_conc
module procedure to_conc_s
module procedure to_conc_all
   module procedure to_conc_m
end interface
interface TRLI0
module procedure TRLI0_s
module procedure TRLI0_all
   module procedure TRLI0_m
end interface
interface src_dist_base
module procedure src_dist_base_s
module procedure src_dist_base_all
   module procedure src_dist_base_m
end interface
interface src_dist_index
module procedure src_dist_index_s
module procedure src_dist_index_all
   module procedure src_dist_index_m
end interface


  integer, parameter :: MAX_LEN_NAME = 20
!**** parameters for tr_wd_TYPE
!@param nGAS   index for wetdep tracer type = gas
!@param nPART  index for wetdep tracer type = particle/aerosol
!@param nWATER index for wetdep tracer type = water
      integer, parameter :: nGAS=1, nPART=2, nWATER=3

  type OldTracer_type
    private
    character(len=MAX_LEN_NAME) :: name

!@var tr_mm: !@var TR_MM: molecular mass of each tracer (g/mole)
    real*8 :: tr_mm  
!@var ntm_power: Power of 10 associated with each tracer (for printing)
    integer :: ntm_power  
!@var t_qlimit: if t_qlimit=.true. tracer is maintained as positive
    logical :: t_qlimit = .true. 
!@var needtrs: true if surface tracer value from PBL is required
    logical :: needtrs = .false. 
!@var trdecay: radioactive decay constant (1/s) (=0 for stable tracers)
    real*8 :: trdecay = 0.d0 
!@var itime_tr0: start time for each tracer (hours)
    integer :: itime_tr0  
!@var trsi0: conc. in sea ice (kg/m^2)
    real*8 :: trsi0 = 0.d0 
!@var trw0: concentration in water (kg/kg)
    real*8 :: trw0 = 0.d0 
!@var mass2vol: mass to volume ratio = mair/tr_mm
    real*8 :: mass2vol  
!@var vol2mass: volume to mass ratio = tr_mm/mair
    real*8 :: vol2mass  
!@var dodrydep: true if tracer should undergo dry deposition
    logical :: dodrydep = .false. 
!@var F0: reactivity factor for oxidation of biological substances
    real*8 :: F0 = 0.d0 
!@var HSTAR: Henry's Law const for tracer dry deposition. mole/(L atm)
!@+   Same as the tr_RKD wet dep variable
    real*8 :: HSTAR = 0.d0 
!@var do_fire: true if tracer should have emissions via flammability
    logical :: do_fire = .false. 
!@var nBBsources: number of sources attributed to biomass burning
    integer :: nBBsources = 0 
!@var emisPerFireByVegType: emisPerFireByVegType tracer emissions per fire count as a
!@+ function of 12 standard GISS (VDATA) vegetation types
    real*8, dimension(12) :: emisPerFireByVegType = 0.0d0 
!@var trpdens: tracer particle density (kg/m^3)
!@+               (=0 for non-particle tracers)
    real*8 :: trpdens = 0.d0 
!@var trradius: tracer effective radius (m) (=0 for non particle tracers)
    real*8 :: trradius = 0.d0 
!@var tr_wd_TYPE: tracer wet dep type (gas, particle, water)
    integer :: tr_wd_TYPE = nGas 
!@var tr_RKD: Henry's Law coefficient (in mole/Joule please !)
    real*8 :: tr_RKD = 0.d0 
!@var tr_DHD: coefficient of temperature-dependence term of Henry's
!@+   Law coefficient (in Joule/mole please !)
    real*8 :: tr_DHD = 0.d0 
!@var fq_aer: fraction of aerosol that condenses
    real*8 :: fq_aer = 0.d0 
!@var rc_washt: aerosol washout rate
    real*8 :: rc_washt = 1.d-1 
!@var isDust: index array for testing if tracer is a dust type
    integer :: isDust = 0 
!@var tr_H2ObyCH4: conc. of tracer in water from methane oxidation
    real*8 :: tr_H2ObyCH4 = 0.d0 
!@var dowetdep: true if tracer has some form of wet deposition
    logical :: dowetdep = .false. 
!@var ntrocn: scaling power factor for ocean/ice tracer concentrations
    integer :: ntrocn = 0 
!@var conc_from_fw: true if ocean conc is defined using fresh water
    logical :: conc_from_fw = .true. 
!@var trglac: tracer ratio in glacial runoff to ocean (kg/kg)
    real*8 :: trglac  
!@var ntisurfsrc: no. of interactive surface sources for each tracer
    integer :: ntisurfsrc  
!@var iso_index: indexing taking actual tracer number to isotope
    integer :: iso_index = 1 
!@var om2oc: ratio of organic matter to organic carbon
    real*8 :: om2oc = 1.4d0 
!@var to_volume_MixRat: =0: print tracer conc. by vol mix ratio; =1:mass mixing ratio
    integer :: to_volume_MixRat = 0 
!@var to_conc: =0: print 3D tracer conc. by to_volume_MixRat; =1: kg/m3
    integer :: to_conc = 0 
!@var TRLI0: default tracer conc. for land ice (kg/kg)
    real*8 :: TRLI0 = 0.d0 
!@var src_dist_base: source distribution base
    integer :: src_dist_base = 0 
!@var src_dist_index: source distribution index
    integer :: src_dist_index = 0 


  end type OldTracer_type

  integer, save :: numTracers = 0
  type (OldTracer_type), allocatable :: internalTracers(:)
  procedure(IdefaultSpec), pointer :: defaultSpec

  abstract interface
    subroutine IdefaultSpec(n_idx, trcer)
    use Tracer_mod
    use Attributes_mod, only: assignment(=)
    integer, intent(in) :: n_idx
    class (Tracer), pointer :: trcer
    end subroutine IdefaultSpec
  end interface

  interface trName
    module procedure trname_s
    module procedure trname_all
  end interface trName
  
  type (TracerBundle), pointer :: tracerReference => null()

contains

  function trName_s(idx) result(name)
    integer, intent(in) :: idx
    character(len=len_trim(internalTracers(idx)%name)) :: name
    name = trim(internalTracers(idx)%name)
  end function trName_s

  function trName_all() result(name)
    character(len=MAX_LEN_NAME) :: name(numTracers)
    name = internalTracers(:)%name
  end function trName_all

  subroutine initializeOldTracers(tracerRef, spec)
    type (TracerBundle), target :: tracerRef
    procedure(IdefaultSpec) :: spec

    tracerReference => tracerRef
    defaultSpec => spec
    allocate(internalTracers(0))
    
  end subroutine initializeOldTracers

  integer function oldAddTracer(name, initname) result(n)
     use Attributes_mod, only: assignment(=)
    character(len=*), intent(in) :: name
    character(len=*), intent(in), optional :: initname
    type (OldTracer_type), allocatable :: tmp(:)
    class (Tracer), pointer :: t, it=>null()
    integer :: i

    allocate(tmp(numTracers))
    tmp = internalTracers
    deallocate(internalTracers)

    allocate(internalTracers(numTracers+1))
    internalTracers(1:numTracers) = tmp
    deallocate(tmp)

    if (present(initname)) then
      do i=1, numtracers
        if (internaltracers(i)%name==initname) then
          internaltracers(numtracers+1)=internaltracers(i)
          exit
        endif
      end do
    endif

    numTracers = numTracers + 1
    internalTracers(numTracers)%name = trim(name)

    call tracerReference%insert(name, newTracer(name))

!    t => newTracer(name)
!    call tracerReference%insert(name, t)
    t => tracerReference%getReference(name)
    if (present(initname)) it=>tracerreference%getreference(initname)
    if (associated(it)) then
      call copyinto(t, it)
    else
      call defaultSpec(numTracers, t)
    endif
    call t%insert('index', numTracers)
    n = numTracers

  end function oldAddTracer

  function findtracer(name)
  implicit none
  character(len=*), intent(in) :: name
  integer :: findtracer, i

#ifdef NEW_TRACER_PROPERTIES
  findtracer=tracerreference%getproperty(name, 'index')
#else
  findtracer=0
  do i=1, numtracers
#ifdef MIXED_TRACER_PROPERTIES
    if (tracerreference%internaltracers(i)%getproperty('name')==name) then
#else
    if (internaltracers(i)%name==name) then
#endif
#endif
      findtracer=i
      exit
    endif
  end do
  return
  end function findtracer

    subroutine set_tr_mm(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%tr_mm = value
    call tracerReference%setAttribute(trName(oldIndex), "tr_mm", newAttribute(value))
  end subroutine set_tr_mm
  
  function tr_mm_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: tr_mm_s
#ifdef NEW_TRACER_PROPERTIES
    tr_mm_s = tracerReference%getProperty(trName(oldIndex), "tr_mm")
#else
#ifdef MIXED_TRACER_PROPERTIES
    tr_mm_s = tracerReference%internalTracers(oldIndex)%getProperty("tr_mm")
#else
    tr_mm_s = internalTracers(oldIndex)%tr_mm
#endif
#endif
  end function tr_mm_s

  function tr_mm_all()
    real*8 :: tr_mm_all(size(internalTracers))
    tr_mm_all = internalTracers(:)%tr_mm
  end function tr_mm_all

  function tr_mm_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_mm_m(size(oldIndices))
    tr_mm_m = internalTracers(oldIndices(:))%tr_mm
  end function tr_mm_m



  subroutine set_ntm_power(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%ntm_power = value
    call tracerReference%setAttribute(trName(oldIndex), "ntm_power", newAttribute(value))
  end subroutine set_ntm_power
  
  function ntm_power_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: ntm_power_s
#ifdef NEW_TRACER_PROPERTIES
    ntm_power_s = tracerReference%getProperty(trName(oldIndex), "ntm_power")
#else
#ifdef MIXED_TRACER_PROPERTIES
    ntm_power_s = tracerReference%internalTracers(oldIndex)%getProperty("ntm_power")
#else
    ntm_power_s = internalTracers(oldIndex)%ntm_power
#endif
#endif
  end function ntm_power_s

  function ntm_power_all()
    integer :: ntm_power_all(size(internalTracers))
    ntm_power_all = internalTracers(:)%ntm_power
  end function ntm_power_all

  function ntm_power_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntm_power_m(size(oldIndices))
    ntm_power_m = internalTracers(oldIndices(:))%ntm_power
  end function ntm_power_m



  subroutine set_t_qlimit(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%t_qlimit = value
    call tracerReference%setAttribute(trName(oldIndex), "t_qlimit", newAttribute(value))
  end subroutine set_t_qlimit
  
  function t_qlimit_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: t_qlimit_s
#ifdef NEW_TRACER_PROPERTIES
    t_qlimit_s = tracerReference%getProperty(trName(oldIndex), "t_qlimit")
#else
#ifdef MIXED_TRACER_PROPERTIES
    t_qlimit_s = tracerReference%internalTracers(oldIndex)%getProperty("t_qlimit")
#else
    t_qlimit_s = internalTracers(oldIndex)%t_qlimit
#endif
#endif
  end function t_qlimit_s

  function t_qlimit_all()
    logical :: t_qlimit_all(size(internalTracers))
    t_qlimit_all = internalTracers(:)%t_qlimit
  end function t_qlimit_all

  function t_qlimit_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: t_qlimit_m(size(oldIndices))
    t_qlimit_m = internalTracers(oldIndices(:))%t_qlimit
  end function t_qlimit_m



  subroutine set_needtrs(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%needtrs = value
    call tracerReference%setAttribute(trName(oldIndex), "needtrs", newAttribute(value))
  end subroutine set_needtrs
  
  function needtrs_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: needtrs_s
#ifdef NEW_TRACER_PROPERTIES
    needtrs_s = tracerReference%getProperty(trName(oldIndex), "needtrs")
#else
#ifdef MIXED_TRACER_PROPERTIES
    needtrs_s = tracerReference%internalTracers(oldIndex)%getProperty("needtrs")
#else
    needtrs_s = internalTracers(oldIndex)%needtrs
#endif
#endif
  end function needtrs_s

  function needtrs_all()
    logical :: needtrs_all(size(internalTracers))
    needtrs_all = internalTracers(:)%needtrs
  end function needtrs_all

  function needtrs_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: needtrs_m(size(oldIndices))
    needtrs_m = internalTracers(oldIndices(:))%needtrs
  end function needtrs_m



  subroutine set_trdecay(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trdecay = value
    call tracerReference%setAttribute(trName(oldIndex), "trdecay", newAttribute(value))
  end subroutine set_trdecay
  
  function trdecay_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trdecay_s
#ifdef NEW_TRACER_PROPERTIES
    trdecay_s = tracerReference%getProperty(trName(oldIndex), "trdecay")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trdecay_s = tracerReference%internalTracers(oldIndex)%getProperty("trdecay")
#else
    trdecay_s = internalTracers(oldIndex)%trdecay
#endif
#endif
  end function trdecay_s

  function trdecay_all()
    real*8 :: trdecay_all(size(internalTracers))
    trdecay_all = internalTracers(:)%trdecay
  end function trdecay_all

  function trdecay_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trdecay_m(size(oldIndices))
    trdecay_m = internalTracers(oldIndices(:))%trdecay
  end function trdecay_m



  subroutine set_itime_tr0(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%itime_tr0 = value
    call tracerReference%setAttribute(trName(oldIndex), "itime_tr0", newAttribute(value))
  end subroutine set_itime_tr0
  
  function itime_tr0_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: itime_tr0_s
#ifdef NEW_TRACER_PROPERTIES
    itime_tr0_s = tracerReference%getProperty(trName(oldIndex), "itime_tr0")
#else
#ifdef MIXED_TRACER_PROPERTIES
    itime_tr0_s = tracerReference%internalTracers(oldIndex)%getProperty("itime_tr0")
#else
    itime_tr0_s = internalTracers(oldIndex)%itime_tr0
#endif
#endif
  end function itime_tr0_s

  function itime_tr0_all()
    integer :: itime_tr0_all(size(internalTracers))
    itime_tr0_all = internalTracers(:)%itime_tr0
  end function itime_tr0_all

  function itime_tr0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: itime_tr0_m(size(oldIndices))
    itime_tr0_m = internalTracers(oldIndices(:))%itime_tr0
  end function itime_tr0_m



  subroutine set_trsi0(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trsi0 = value
    call tracerReference%setAttribute(trName(oldIndex), "trsi0", newAttribute(value))
  end subroutine set_trsi0
  
  function trsi0_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trsi0_s
#ifdef NEW_TRACER_PROPERTIES
    trsi0_s = tracerReference%getProperty(trName(oldIndex), "trsi0")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trsi0_s = tracerReference%internalTracers(oldIndex)%getProperty("trsi0")
#else
    trsi0_s = internalTracers(oldIndex)%trsi0
#endif
#endif
  end function trsi0_s

  function trsi0_all()
    real*8 :: trsi0_all(size(internalTracers))
    trsi0_all = internalTracers(:)%trsi0
  end function trsi0_all

  function trsi0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trsi0_m(size(oldIndices))
    trsi0_m = internalTracers(oldIndices(:))%trsi0
  end function trsi0_m



  subroutine set_trw0(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trw0 = value
    call tracerReference%setAttribute(trName(oldIndex), "trw0", newAttribute(value))
  end subroutine set_trw0
  
  function trw0_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trw0_s
#ifdef NEW_TRACER_PROPERTIES
    trw0_s = tracerReference%getProperty(trName(oldIndex), "trw0")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trw0_s = tracerReference%internalTracers(oldIndex)%getProperty("trw0")
#else
    trw0_s = internalTracers(oldIndex)%trw0
#endif
#endif
  end function trw0_s

  function trw0_all()
    real*8 :: trw0_all(size(internalTracers))
    trw0_all = internalTracers(:)%trw0
  end function trw0_all

  function trw0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trw0_m(size(oldIndices))
    trw0_m = internalTracers(oldIndices(:))%trw0
  end function trw0_m



  subroutine set_mass2vol(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%mass2vol = value
    call tracerReference%setAttribute(trName(oldIndex), "mass2vol", newAttribute(value))
  end subroutine set_mass2vol
  
  function mass2vol_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: mass2vol_s
#ifdef NEW_TRACER_PROPERTIES
    mass2vol_s = tracerReference%getProperty(trName(oldIndex), "mass2vol")
#else
#ifdef MIXED_TRACER_PROPERTIES
    mass2vol_s = tracerReference%internalTracers(oldIndex)%getProperty("mass2vol")
#else
    mass2vol_s = internalTracers(oldIndex)%mass2vol
#endif
#endif
  end function mass2vol_s

  function mass2vol_all()
    real*8 :: mass2vol_all(size(internalTracers))
    mass2vol_all = internalTracers(:)%mass2vol
  end function mass2vol_all

  function mass2vol_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: mass2vol_m(size(oldIndices))
    mass2vol_m = internalTracers(oldIndices(:))%mass2vol
  end function mass2vol_m



  subroutine set_vol2mass(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%vol2mass = value
    call tracerReference%setAttribute(trName(oldIndex), "vol2mass", newAttribute(value))
  end subroutine set_vol2mass
  
  function vol2mass_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: vol2mass_s
#ifdef NEW_TRACER_PROPERTIES
    vol2mass_s = tracerReference%getProperty(trName(oldIndex), "vol2mass")
#else
#ifdef MIXED_TRACER_PROPERTIES
    vol2mass_s = tracerReference%internalTracers(oldIndex)%getProperty("vol2mass")
#else
    vol2mass_s = internalTracers(oldIndex)%vol2mass
#endif
#endif
  end function vol2mass_s

  function vol2mass_all()
    real*8 :: vol2mass_all(size(internalTracers))
    vol2mass_all = internalTracers(:)%vol2mass
  end function vol2mass_all

  function vol2mass_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: vol2mass_m(size(oldIndices))
    vol2mass_m = internalTracers(oldIndices(:))%vol2mass
  end function vol2mass_m



  subroutine set_dodrydep(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%dodrydep = value
    call tracerReference%setAttribute(trName(oldIndex), "dodrydep", newAttribute(value))
  end subroutine set_dodrydep
  
  function dodrydep_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: dodrydep_s
#ifdef NEW_TRACER_PROPERTIES
    dodrydep_s = tracerReference%getProperty(trName(oldIndex), "dodrydep")
#else
#ifdef MIXED_TRACER_PROPERTIES
    dodrydep_s = tracerReference%internalTracers(oldIndex)%getProperty("dodrydep")
#else
    dodrydep_s = internalTracers(oldIndex)%dodrydep
#endif
#endif
  end function dodrydep_s

  function dodrydep_all()
    logical :: dodrydep_all(size(internalTracers))
    dodrydep_all = internalTracers(:)%dodrydep
  end function dodrydep_all

  function dodrydep_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: dodrydep_m(size(oldIndices))
    dodrydep_m = internalTracers(oldIndices(:))%dodrydep
  end function dodrydep_m



  subroutine set_F0(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%F0 = value
    call tracerReference%setAttribute(trName(oldIndex), "F0", newAttribute(value))
  end subroutine set_F0
  
  function F0_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: F0_s
#ifdef NEW_TRACER_PROPERTIES
    F0_s = tracerReference%getProperty(trName(oldIndex), "F0")
#else
#ifdef MIXED_TRACER_PROPERTIES
    F0_s = tracerReference%internalTracers(oldIndex)%getProperty("F0")
#else
    F0_s = internalTracers(oldIndex)%F0
#endif
#endif
  end function F0_s

  function F0_all()
    real*8 :: F0_all(size(internalTracers))
    F0_all = internalTracers(:)%F0
  end function F0_all

  function F0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: F0_m(size(oldIndices))
    F0_m = internalTracers(oldIndices(:))%F0
  end function F0_m



  subroutine set_HSTAR(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%HSTAR = value
    call tracerReference%setAttribute(trName(oldIndex), "HSTAR", newAttribute(value))
  end subroutine set_HSTAR
  
  function HSTAR_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: HSTAR_s
#ifdef NEW_TRACER_PROPERTIES
    HSTAR_s = tracerReference%getProperty(trName(oldIndex), "HSTAR")
#else
#ifdef MIXED_TRACER_PROPERTIES
    HSTAR_s = tracerReference%internalTracers(oldIndex)%getProperty("HSTAR")
#else
    HSTAR_s = internalTracers(oldIndex)%HSTAR
#endif
#endif
  end function HSTAR_s

  function HSTAR_all()
    real*8 :: HSTAR_all(size(internalTracers))
    HSTAR_all = internalTracers(:)%HSTAR
  end function HSTAR_all

  function HSTAR_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: HSTAR_m(size(oldIndices))
    HSTAR_m = internalTracers(oldIndices(:))%HSTAR
  end function HSTAR_m



  subroutine set_do_fire(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%do_fire = value
    call tracerReference%setAttribute(trName(oldIndex), "do_fire", newAttribute(value))
  end subroutine set_do_fire
  
  function do_fire_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: do_fire_s
#ifdef NEW_TRACER_PROPERTIES
    do_fire_s = tracerReference%getProperty(trName(oldIndex), "do_fire")
#else
#ifdef MIXED_TRACER_PROPERTIES
    do_fire_s = tracerReference%internalTracers(oldIndex)%getProperty("do_fire")
#else
    do_fire_s = internalTracers(oldIndex)%do_fire
#endif
#endif
  end function do_fire_s

  function do_fire_all()
    logical :: do_fire_all(size(internalTracers))
    do_fire_all = internalTracers(:)%do_fire
  end function do_fire_all

  function do_fire_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: do_fire_m(size(oldIndices))
    do_fire_m = internalTracers(oldIndices(:))%do_fire
  end function do_fire_m



  subroutine set_nBBsources(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%nBBsources = value
    call tracerReference%setAttribute(trName(oldIndex), "nBBsources", newAttribute(value))
  end subroutine set_nBBsources
  
  function nBBsources_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: nBBsources_s
#ifdef NEW_TRACER_PROPERTIES
    nBBsources_s = tracerReference%getProperty(trName(oldIndex), "nBBsources")
#else
#ifdef MIXED_TRACER_PROPERTIES
    nBBsources_s = tracerReference%internalTracers(oldIndex)%getProperty("nBBsources")
#else
    nBBsources_s = internalTracers(oldIndex)%nBBsources
#endif
#endif
  end function nBBsources_s

  function nBBsources_all()
    integer :: nBBsources_all(size(internalTracers))
    nBBsources_all = internalTracers(:)%nBBsources
  end function nBBsources_all

  function nBBsources_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: nBBsources_m(size(oldIndices))
    nBBsources_m = internalTracers(oldIndices(:))%nBBsources
  end function nBBsources_m



  subroutine set_emisPerFireByVegType(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, dimension(12), intent(in) :: value
    internalTracers(oldIndex)%emisPerFireByVegType = value
    call tracerReference%setAttribute(trName(oldIndex), "emisPerFireByVegType", newAttribute(value))
  end subroutine set_emisPerFireByVegType
  
  function emisPerFireByVegType_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8, dimension(12) :: emisPerFireByVegType_s
#ifdef NEW_TRACER_PROPERTIES
    emisPerFireByVegType_s = tracerReference%getProperty(trName(oldIndex), "emisPerFireByVegType")
#else
#ifdef MIXED_TRACER_PROPERTIES
    emisPerFireByVegType_s = tracerReference%internalTracers(oldIndex)%getProperty("emisPerFireByVegType")
#else
    emisPerFireByVegType_s = internalTracers(oldIndex)%emisPerFireByVegType
#endif
#endif
  end function emisPerFireByVegType_s



  subroutine set_trpdens(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trpdens = value
    call tracerReference%setAttribute(trName(oldIndex), "trpdens", newAttribute(value))
  end subroutine set_trpdens
  
  function trpdens_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trpdens_s
#ifdef NEW_TRACER_PROPERTIES
    trpdens_s = tracerReference%getProperty(trName(oldIndex), "trpdens")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trpdens_s = tracerReference%internalTracers(oldIndex)%getProperty("trpdens")
#else
    trpdens_s = internalTracers(oldIndex)%trpdens
#endif
#endif
  end function trpdens_s

  function trpdens_all()
    real*8 :: trpdens_all(size(internalTracers))
    trpdens_all = internalTracers(:)%trpdens
  end function trpdens_all

  function trpdens_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trpdens_m(size(oldIndices))
    trpdens_m = internalTracers(oldIndices(:))%trpdens
  end function trpdens_m



  subroutine set_trradius(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trradius = value
    call tracerReference%setAttribute(trName(oldIndex), "trradius", newAttribute(value))
  end subroutine set_trradius
  
  function trradius_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trradius_s
#ifdef NEW_TRACER_PROPERTIES
    trradius_s = tracerReference%getProperty(trName(oldIndex), "trradius")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trradius_s = tracerReference%internalTracers(oldIndex)%getProperty("trradius")
#else
    trradius_s = internalTracers(oldIndex)%trradius
#endif
#endif
  end function trradius_s

  function trradius_all()
    real*8 :: trradius_all(size(internalTracers))
    trradius_all = internalTracers(:)%trradius
  end function trradius_all

  function trradius_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trradius_m(size(oldIndices))
    trradius_m = internalTracers(oldIndices(:))%trradius
  end function trradius_m



  subroutine set_tr_wd_TYPE(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%tr_wd_TYPE = value
    call tracerReference%setAttribute(trName(oldIndex), "tr_wd_TYPE", newAttribute(value))
  end subroutine set_tr_wd_TYPE
  
  function tr_wd_TYPE_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: tr_wd_TYPE_s
#ifdef NEW_TRACER_PROPERTIES
    tr_wd_TYPE_s = tracerReference%getProperty(trName(oldIndex), "tr_wd_TYPE")
#else
#ifdef MIXED_TRACER_PROPERTIES
    tr_wd_TYPE_s = tracerReference%internalTracers(oldIndex)%getProperty("tr_wd_TYPE")
#else
    tr_wd_TYPE_s = internalTracers(oldIndex)%tr_wd_TYPE
#endif
#endif
  end function tr_wd_TYPE_s

  function tr_wd_TYPE_all()
    integer :: tr_wd_TYPE_all(size(internalTracers))
    tr_wd_TYPE_all = internalTracers(:)%tr_wd_TYPE
  end function tr_wd_TYPE_all

  function tr_wd_TYPE_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: tr_wd_TYPE_m(size(oldIndices))
    tr_wd_TYPE_m = internalTracers(oldIndices(:))%tr_wd_TYPE
  end function tr_wd_TYPE_m



  subroutine set_tr_RKD(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%tr_RKD = value
    call tracerReference%setAttribute(trName(oldIndex), "tr_RKD", newAttribute(value))
  end subroutine set_tr_RKD
  
  function tr_RKD_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: tr_RKD_s
#ifdef NEW_TRACER_PROPERTIES
    tr_RKD_s = tracerReference%getProperty(trName(oldIndex), "tr_RKD")
#else
#ifdef MIXED_TRACER_PROPERTIES
    tr_RKD_s = tracerReference%internalTracers(oldIndex)%getProperty("tr_RKD")
#else
    tr_RKD_s = internalTracers(oldIndex)%tr_RKD
#endif
#endif
  end function tr_RKD_s

  function tr_RKD_all()
    real*8 :: tr_RKD_all(size(internalTracers))
    tr_RKD_all = internalTracers(:)%tr_RKD
  end function tr_RKD_all

  function tr_RKD_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_RKD_m(size(oldIndices))
    tr_RKD_m = internalTracers(oldIndices(:))%tr_RKD
  end function tr_RKD_m



  subroutine set_tr_DHD(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%tr_DHD = value
    call tracerReference%setAttribute(trName(oldIndex), "tr_DHD", newAttribute(value))
  end subroutine set_tr_DHD
  
  function tr_DHD_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: tr_DHD_s
#ifdef NEW_TRACER_PROPERTIES
    tr_DHD_s = tracerReference%getProperty(trName(oldIndex), "tr_DHD")
#else
#ifdef MIXED_TRACER_PROPERTIES
    tr_DHD_s = tracerReference%internalTracers(oldIndex)%getProperty("tr_DHD")
#else
    tr_DHD_s = internalTracers(oldIndex)%tr_DHD
#endif
#endif
  end function tr_DHD_s

  function tr_DHD_all()
    real*8 :: tr_DHD_all(size(internalTracers))
    tr_DHD_all = internalTracers(:)%tr_DHD
  end function tr_DHD_all

  function tr_DHD_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_DHD_m(size(oldIndices))
    tr_DHD_m = internalTracers(oldIndices(:))%tr_DHD
  end function tr_DHD_m



  subroutine set_fq_aer(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%fq_aer = value
    call tracerReference%setAttribute(trName(oldIndex), "fq_aer", newAttribute(value))
  end subroutine set_fq_aer
  
  function fq_aer_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: fq_aer_s
#ifdef NEW_TRACER_PROPERTIES
    fq_aer_s = tracerReference%getProperty(trName(oldIndex), "fq_aer")
#else
#ifdef MIXED_TRACER_PROPERTIES
    fq_aer_s = tracerReference%internalTracers(oldIndex)%getProperty("fq_aer")
#else
    fq_aer_s = internalTracers(oldIndex)%fq_aer
#endif
#endif
  end function fq_aer_s

  function fq_aer_all()
    real*8 :: fq_aer_all(size(internalTracers))
    fq_aer_all = internalTracers(:)%fq_aer
  end function fq_aer_all

  function fq_aer_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: fq_aer_m(size(oldIndices))
    fq_aer_m = internalTracers(oldIndices(:))%fq_aer
  end function fq_aer_m



  subroutine set_rc_washt(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%rc_washt = value
    call tracerReference%setAttribute(trName(oldIndex), "rc_washt", newAttribute(value))
  end subroutine set_rc_washt
  
  function rc_washt_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: rc_washt_s
#ifdef NEW_TRACER_PROPERTIES
    rc_washt_s = tracerReference%getProperty(trName(oldIndex), "rc_washt")
#else
#ifdef MIXED_TRACER_PROPERTIES
    rc_washt_s = tracerReference%internalTracers(oldIndex)%getProperty("rc_washt")
#else
    rc_washt_s = internalTracers(oldIndex)%rc_washt
#endif
#endif
  end function rc_washt_s

  function rc_washt_all()
    real*8 :: rc_washt_all(size(internalTracers))
    rc_washt_all = internalTracers(:)%rc_washt
  end function rc_washt_all

  function rc_washt_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: rc_washt_m(size(oldIndices))
    rc_washt_m = internalTracers(oldIndices(:))%rc_washt
  end function rc_washt_m



  subroutine set_isDust(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%isDust = value
    call tracerReference%setAttribute(trName(oldIndex), "isDust", newAttribute(value))
  end subroutine set_isDust
  
  function isDust_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: isDust_s
#ifdef NEW_TRACER_PROPERTIES
    isDust_s = tracerReference%getProperty(trName(oldIndex), "isDust")
#else
#ifdef MIXED_TRACER_PROPERTIES
    isDust_s = tracerReference%internalTracers(oldIndex)%getProperty("isDust")
#else
    isDust_s = internalTracers(oldIndex)%isDust
#endif
#endif
  end function isDust_s

  function isDust_all()
    integer :: isDust_all(size(internalTracers))
    isDust_all = internalTracers(:)%isDust
  end function isDust_all

  function isDust_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: isDust_m(size(oldIndices))
    isDust_m = internalTracers(oldIndices(:))%isDust
  end function isDust_m



  subroutine set_tr_H2ObyCH4(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%tr_H2ObyCH4 = value
    call tracerReference%setAttribute(trName(oldIndex), "tr_H2ObyCH4", newAttribute(value))
  end subroutine set_tr_H2ObyCH4
  
  function tr_H2ObyCH4_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: tr_H2ObyCH4_s
#ifdef NEW_TRACER_PROPERTIES
    tr_H2ObyCH4_s = tracerReference%getProperty(trName(oldIndex), "tr_H2ObyCH4")
#else
#ifdef MIXED_TRACER_PROPERTIES
    tr_H2ObyCH4_s = tracerReference%internalTracers(oldIndex)%getProperty("tr_H2ObyCH4")
#else
    tr_H2ObyCH4_s = internalTracers(oldIndex)%tr_H2ObyCH4
#endif
#endif
  end function tr_H2ObyCH4_s

  function tr_H2ObyCH4_all()
    real*8 :: tr_H2ObyCH4_all(size(internalTracers))
    tr_H2ObyCH4_all = internalTracers(:)%tr_H2ObyCH4
  end function tr_H2ObyCH4_all

  function tr_H2ObyCH4_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: tr_H2ObyCH4_m(size(oldIndices))
    tr_H2ObyCH4_m = internalTracers(oldIndices(:))%tr_H2ObyCH4
  end function tr_H2ObyCH4_m



  subroutine set_dowetdep(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%dowetdep = value
    call tracerReference%setAttribute(trName(oldIndex), "dowetdep", newAttribute(value))
  end subroutine set_dowetdep
  
  function dowetdep_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: dowetdep_s
#ifdef NEW_TRACER_PROPERTIES
    dowetdep_s = tracerReference%getProperty(trName(oldIndex), "dowetdep")
#else
#ifdef MIXED_TRACER_PROPERTIES
    dowetdep_s = tracerReference%internalTracers(oldIndex)%getProperty("dowetdep")
#else
    dowetdep_s = internalTracers(oldIndex)%dowetdep
#endif
#endif
  end function dowetdep_s

  function dowetdep_all()
    logical :: dowetdep_all(size(internalTracers))
    dowetdep_all = internalTracers(:)%dowetdep
  end function dowetdep_all

  function dowetdep_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: dowetdep_m(size(oldIndices))
    dowetdep_m = internalTracers(oldIndices(:))%dowetdep
  end function dowetdep_m



  subroutine set_ntrocn(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%ntrocn = value
    call tracerReference%setAttribute(trName(oldIndex), "ntrocn", newAttribute(value))
  end subroutine set_ntrocn
  
  function ntrocn_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: ntrocn_s
#ifdef NEW_TRACER_PROPERTIES
    ntrocn_s = tracerReference%getProperty(trName(oldIndex), "ntrocn")
#else
#ifdef MIXED_TRACER_PROPERTIES
    ntrocn_s = tracerReference%internalTracers(oldIndex)%getProperty("ntrocn")
#else
    ntrocn_s = internalTracers(oldIndex)%ntrocn
#endif
#endif
  end function ntrocn_s

  function ntrocn_all()
    integer :: ntrocn_all(size(internalTracers))
    ntrocn_all = internalTracers(:)%ntrocn
  end function ntrocn_all

  function ntrocn_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntrocn_m(size(oldIndices))
    ntrocn_m = internalTracers(oldIndices(:))%ntrocn
  end function ntrocn_m



  subroutine set_conc_from_fw(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    logical, intent(in) :: value
    internalTracers(oldIndex)%conc_from_fw = value
    call tracerReference%setAttribute(trName(oldIndex), "conc_from_fw", newAttribute(value))
  end subroutine set_conc_from_fw
  
  function conc_from_fw_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    logical :: conc_from_fw_s
#ifdef NEW_TRACER_PROPERTIES
    conc_from_fw_s = tracerReference%getProperty(trName(oldIndex), "conc_from_fw")
#else
#ifdef MIXED_TRACER_PROPERTIES
    conc_from_fw_s = tracerReference%internalTracers(oldIndex)%getProperty("conc_from_fw")
#else
    conc_from_fw_s = internalTracers(oldIndex)%conc_from_fw
#endif
#endif
  end function conc_from_fw_s

  function conc_from_fw_all()
    logical :: conc_from_fw_all(size(internalTracers))
    conc_from_fw_all = internalTracers(:)%conc_from_fw
  end function conc_from_fw_all

  function conc_from_fw_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    logical :: conc_from_fw_m(size(oldIndices))
    conc_from_fw_m = internalTracers(oldIndices(:))%conc_from_fw
  end function conc_from_fw_m



  subroutine set_trglac(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%trglac = value
    call tracerReference%setAttribute(trName(oldIndex), "trglac", newAttribute(value))
  end subroutine set_trglac
  
  function trglac_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: trglac_s
#ifdef NEW_TRACER_PROPERTIES
    trglac_s = tracerReference%getProperty(trName(oldIndex), "trglac")
#else
#ifdef MIXED_TRACER_PROPERTIES
    trglac_s = tracerReference%internalTracers(oldIndex)%getProperty("trglac")
#else
    trglac_s = internalTracers(oldIndex)%trglac
#endif
#endif
  end function trglac_s

  function trglac_all()
    real*8 :: trglac_all(size(internalTracers))
    trglac_all = internalTracers(:)%trglac
  end function trglac_all

  function trglac_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: trglac_m(size(oldIndices))
    trglac_m = internalTracers(oldIndices(:))%trglac
  end function trglac_m



  subroutine set_ntisurfsrc(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%ntisurfsrc = value
    call tracerReference%setAttribute(trName(oldIndex), "ntisurfsrc", newAttribute(value))
  end subroutine set_ntisurfsrc
  
  function ntisurfsrc_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: ntisurfsrc_s
#ifdef NEW_TRACER_PROPERTIES
    ntisurfsrc_s = tracerReference%getProperty(trName(oldIndex), "ntisurfsrc")
#else
#ifdef MIXED_TRACER_PROPERTIES
    ntisurfsrc_s = tracerReference%internalTracers(oldIndex)%getProperty("ntisurfsrc")
#else
    ntisurfsrc_s = internalTracers(oldIndex)%ntisurfsrc
#endif
#endif
  end function ntisurfsrc_s

  function ntisurfsrc_all()
    integer :: ntisurfsrc_all(size(internalTracers))
    ntisurfsrc_all = internalTracers(:)%ntisurfsrc
  end function ntisurfsrc_all

  function ntisurfsrc_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: ntisurfsrc_m(size(oldIndices))
    ntisurfsrc_m = internalTracers(oldIndices(:))%ntisurfsrc
  end function ntisurfsrc_m



  subroutine set_iso_index(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%iso_index = value
    call tracerReference%setAttribute(trName(oldIndex), "iso_index", newAttribute(value))
  end subroutine set_iso_index
  
  function iso_index_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: iso_index_s
#ifdef NEW_TRACER_PROPERTIES
    iso_index_s = tracerReference%getProperty(trName(oldIndex), "iso_index")
#else
#ifdef MIXED_TRACER_PROPERTIES
    iso_index_s = tracerReference%internalTracers(oldIndex)%getProperty("iso_index")
#else
    iso_index_s = internalTracers(oldIndex)%iso_index
#endif
#endif
  end function iso_index_s

  function iso_index_all()
    integer :: iso_index_all(size(internalTracers))
    iso_index_all = internalTracers(:)%iso_index
  end function iso_index_all

  function iso_index_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: iso_index_m(size(oldIndices))
    iso_index_m = internalTracers(oldIndices(:))%iso_index
  end function iso_index_m



  subroutine set_om2oc(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%om2oc = value
    call tracerReference%setAttribute(trName(oldIndex), "om2oc", newAttribute(value))
  end subroutine set_om2oc
  
  function om2oc_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: om2oc_s
#ifdef NEW_TRACER_PROPERTIES
    om2oc_s = tracerReference%getProperty(trName(oldIndex), "om2oc")
#else
#ifdef MIXED_TRACER_PROPERTIES
    om2oc_s = tracerReference%internalTracers(oldIndex)%getProperty("om2oc")
#else
    om2oc_s = internalTracers(oldIndex)%om2oc
#endif
#endif
  end function om2oc_s

  function om2oc_all()
    real*8 :: om2oc_all(size(internalTracers))
    om2oc_all = internalTracers(:)%om2oc
  end function om2oc_all

  function om2oc_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: om2oc_m(size(oldIndices))
    om2oc_m = internalTracers(oldIndices(:))%om2oc
  end function om2oc_m



  subroutine set_to_volume_MixRat(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%to_volume_MixRat = value
    call tracerReference%setAttribute(trName(oldIndex), "to_volume_MixRat", newAttribute(value))
  end subroutine set_to_volume_MixRat
  
  function to_volume_MixRat_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: to_volume_MixRat_s
#ifdef NEW_TRACER_PROPERTIES
    to_volume_MixRat_s = tracerReference%getProperty(trName(oldIndex), "to_volume_MixRat")
#else
#ifdef MIXED_TRACER_PROPERTIES
    to_volume_MixRat_s = tracerReference%internalTracers(oldIndex)%getProperty("to_volume_MixRat")
#else
    to_volume_MixRat_s = internalTracers(oldIndex)%to_volume_MixRat
#endif
#endif
  end function to_volume_MixRat_s

  function to_volume_MixRat_all()
    integer :: to_volume_MixRat_all(size(internalTracers))
    to_volume_MixRat_all = internalTracers(:)%to_volume_MixRat
  end function to_volume_MixRat_all

  function to_volume_MixRat_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: to_volume_MixRat_m(size(oldIndices))
    to_volume_MixRat_m = internalTracers(oldIndices(:))%to_volume_MixRat
  end function to_volume_MixRat_m



  subroutine set_to_conc(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%to_conc = value
    call tracerReference%setAttribute(trName(oldIndex), "to_conc", newAttribute(value))
  end subroutine set_to_conc
  
  function to_conc_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: to_conc_s
#ifdef NEW_TRACER_PROPERTIES
    to_conc_s = tracerReference%getProperty(trName(oldIndex), "to_conc")
#else
#ifdef MIXED_TRACER_PROPERTIES
    to_conc_s = tracerReference%internalTracers(oldIndex)%getProperty("to_conc")
#else
    to_conc_s = internalTracers(oldIndex)%to_conc
#endif
#endif
  end function to_conc_s

  function to_conc_all()
    integer :: to_conc_all(size(internalTracers))
    to_conc_all = internalTracers(:)%to_conc
  end function to_conc_all

  function to_conc_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: to_conc_m(size(oldIndices))
    to_conc_m = internalTracers(oldIndices(:))%to_conc
  end function to_conc_m



  subroutine set_TRLI0(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    real*8, intent(in) :: value
    internalTracers(oldIndex)%TRLI0 = value
    call tracerReference%setAttribute(trName(oldIndex), "TRLI0", newAttribute(value))
  end subroutine set_TRLI0
  
  function TRLI0_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    real*8 :: TRLI0_s
#ifdef NEW_TRACER_PROPERTIES
    TRLI0_s = tracerReference%getProperty(trName(oldIndex), "TRLI0")
#else
#ifdef MIXED_TRACER_PROPERTIES
    TRLI0_s = tracerReference%internalTracers(oldIndex)%getProperty("TRLI0")
#else
    TRLI0_s = internalTracers(oldIndex)%TRLI0
#endif
#endif
  end function TRLI0_s

  function TRLI0_all()
    real*8 :: TRLI0_all(size(internalTracers))
    TRLI0_all = internalTracers(:)%TRLI0
  end function TRLI0_all

  function TRLI0_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    real*8 :: TRLI0_m(size(oldIndices))
    TRLI0_m = internalTracers(oldIndices(:))%TRLI0
  end function TRLI0_m



  subroutine set_src_dist_base(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%src_dist_base = value
    call tracerReference%setAttribute(trName(oldIndex), "src_dist_base", newAttribute(value))
  end subroutine set_src_dist_base
  
  function src_dist_base_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: src_dist_base_s
#ifdef NEW_TRACER_PROPERTIES
    src_dist_base_s = tracerReference%getProperty(trName(oldIndex), "src_dist_base")
#else
#ifdef MIXED_TRACER_PROPERTIES
    src_dist_base_s = tracerReference%internalTracers(oldIndex)%getProperty("src_dist_base")
#else
    src_dist_base_s = internalTracers(oldIndex)%src_dist_base
#endif
#endif
  end function src_dist_base_s

  function src_dist_base_all()
    integer :: src_dist_base_all(size(internalTracers))
    src_dist_base_all = internalTracers(:)%src_dist_base
  end function src_dist_base_all

  function src_dist_base_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: src_dist_base_m(size(oldIndices))
    src_dist_base_m = internalTracers(oldIndices(:))%src_dist_base
  end function src_dist_base_m



  subroutine set_src_dist_index(oldIndex, value)
   use Attributes_mod
    integer, intent(in) :: oldIndex
    integer, intent(in) :: value
    internalTracers(oldIndex)%src_dist_index = value
    call tracerReference%setAttribute(trName(oldIndex), "src_dist_index", newAttribute(value))
  end subroutine set_src_dist_index
  
  function src_dist_index_s(oldIndex)
    use GenericType_mod
    integer, intent(in) :: oldIndex
    type (Tracer), pointer :: p
    integer :: src_dist_index_s
#ifdef NEW_TRACER_PROPERTIES
    src_dist_index_s = tracerReference%getProperty(trName(oldIndex), "src_dist_index")
#else
#ifdef MIXED_TRACER_PROPERTIES
    src_dist_index_s = tracerReference%internalTracers(oldIndex)%getProperty("src_dist_index")
#else
    src_dist_index_s = internalTracers(oldIndex)%src_dist_index
#endif
#endif
  end function src_dist_index_s

  function src_dist_index_all()
    integer :: src_dist_index_all(size(internalTracers))
    src_dist_index_all = internalTracers(:)%src_dist_index
  end function src_dist_index_all

  function src_dist_index_m(oldIndices)
    integer, intent(in) :: oldIndices(:)
    integer :: src_dist_index_m(size(oldIndices))
    src_dist_index_m = internalTracers(oldIndices(:))%src_dist_index
  end function src_dist_index_m





end module OldTracer_mod
