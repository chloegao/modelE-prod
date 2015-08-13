module TracerSource_mod
  implicit none
  private

  public :: TracerSource
  public :: TracerSource3D
  public :: N_MAX_SECT


!@param n_max_sect maximum number of sectors for emissions altering
  integer, parameter :: N_MAX_SECT = 10

  type TracerSource
    integer :: num_tr_sectors = 0! number of sectors for a particular tracer and source
    integer :: tr_sect_index(N_MAX_SECT)  = 0 ! array hold the sector index for given tracer/source
    character(len=10) :: tr_sect_name(N_MAX_SECT) = ' '! array hold the sector name for given tracer/source
  end type TracerSource

  type, extends(TracerSource) :: TracerSource3D
  end type TracerSource3D

end module TracerSource_mod
