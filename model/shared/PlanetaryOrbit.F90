!-----------------------------------------------------------------------
! PlanetaryOrbit extends FixedOrbit and is intended to be used for
! planets other than the Earth.  The only functionality added here is
! a constructor for specifying orbital parameters and a trivial
! makeCalendar() method which delegates to the PlanetaryCalendar class.
!-----------------------------------------------------------------------

module PlanetaryOrbit_mod
  use AbstractOrbit_mod
  use FixedOrbit_mod
  use KindParameters_Mod, only: WP => DP, DP
  implicit none
  private

  public :: PlanetaryOrbit

  type, extends(FixedOrbit) :: PlanetaryOrbit
    private
    integer :: foo
  contains
    procedure :: makeCalendar
    procedure :: print_unit
  end type PlanetaryOrbit


  interface PlanetaryOrbit
    module procedure newPlanetaryOrbit
  end interface PlanetaryOrbit


  real(kind=WP), parameter :: EARTH_LON_AT_PERIHELION = 282.9
  real (kind=DP), parameter :: PI = 2*asin(1.d0)

contains


  ! TODO: too many parameters to constructor
  ! TODO: mean distance is missing
  function newPlanetaryOrbit(obliquity, eccentricity, longitudeOfPeriapsis, &
       & siderealPeriod, rotationPeriod, meanDistance) result(orbit)
    use Rational_mod
    use BaseTime_mod
    use TimeInterval_mod
    use OrbitUtilities_mod, only: computeMeanAnomaly
    type (PlanetaryOrbit) :: orbit
    real (kind=WP), intent(in) :: obliquity
    real (kind=WP), intent(in) :: eccentricity
    real (kind=WP), intent(in) :: longitudeOfPeriapsis
    real (kind=WP), intent(in) :: siderealPeriod
    real (kind=WP), intent(in) :: rotationPeriod
    real (kind=WP), intent(in) :: meanDistance

    real (kind=WP) :: meanDay
    type (TimeInterval) :: meanDayInterval
    integer :: daysPerYear
    type (Rational) :: q
    real (kind=WP) :: MA0

    call orbit%setLongitudeAtPeriapsis(longitudeOfPeriapsis)
    call orbit%setObliquity(obliquity)
    call orbit%setEccentricity(eccentricity)

    call orbit%setMeanDistance(meanDistance)

    !--------------------------------------------------------------------------------------
    ! Note sidereal period and rotation period are adjusted to ensure integer days per year
    ! while preserving the length of the mean day.   Other conventions are possible.
    !--------------------------------------------------------------------------------------
    meanDay = 1/(1/rotationPeriod - 1/siderealPeriod)
    daysPerYear = nint(siderealPeriod / meanDay)
    q=Rational(meanDay, tolerance=1.d-6)
    meanDayInterval = TimeInterval(q)
    call orbit%setMeanDay(meanDayInterval)
    call orbit%setSiderealPeriod(TimeInterval(daysPerYear * meanDayInterval))
    call orbit%setRotationPeriod(TimeInterval(meanDayInterval * Rational(daysPerYear, daysPerYear+1)))
    
    MA0 = computeMeanAnomaly(PI/180*(longitudeOfPeriapsis - EARTH_LON_AT_PERIHELION), &
         & eccentricity)
    call orbit%setTimeAtPeriapsis(newBaseTime(MA0/(2*PI) * (daysPerYear*meanDay)))

  end function newPlanetaryOrbit


  ! Pass instance of self to constructor for PlanetaryCalendar.
  function makeCalendar(this) result(calendar)
    use PlanetaryCalendar_mod
    use AbstractCalendar_mod
    class (AbstractCalendar), allocatable :: calendar
    class (PlanetaryOrbit), intent(in) :: this

    allocate(calendar, source=PlanetaryCalendar(this))
  end function makeCalendar


  subroutine print_unit(this, unit)
    class (PlanetaryOrbit), intent(in) :: this
    integer, intent(in) :: unit

    write(unit,*) 'Fixed orbital parameters for planet.'
    write(unit,*) '  Eccentricity:', this%getEccentricity()
    write(unit,*) '  Obliquity (degs):',this%getObliquity()
    write(unit,*) '  Longitude at periapsis (degs from ve):', &
         & this%getLongitudeAtPeriapsis()

  end subroutine print_unit

end module PlanetaryOrbit_mod
