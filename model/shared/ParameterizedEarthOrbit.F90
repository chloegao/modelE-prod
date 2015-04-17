!
!
module ParameterizedEarthOrbit_mod
  use AbstractOrbit_mod
  use AbstractCalendar_mod
  use KindParameters_mod, only: WP => DP, DP
  use TimeInterval_mod
  use BaseTime_mod
  use Rational_mod
  implicit none
  private

  public :: ParameterizedEarthOrbit
  public :: newParameterizedEarthOrbit

  type, extends(AbstractOrbit) :: ParameterizedEarthOrbit
    private
    class (AbstractCalendar), allocatable :: calendar

    integer :: yearBeforePresent = -1

    type (BaseTime) :: timeAtPeriapsis    ! seconds
    type (BaseTime) :: timeAtVernalEquinox

    real (kind=WP) :: eccentricity
    real (kind=WP) :: obliquity
    real (kind=WP) :: longitudeAtPeriapsis
    
  contains

    procedure :: setYear

    procedure :: getEccentricity
    procedure :: getObliquity
    procedure :: getLongitudeAtPeriapsis

    procedure :: setTimeAtPeriapsis
    procedure :: getTimeAtPeriapsis

    procedure :: getMeanAnomaly
    procedure :: getTrueAnomaly

    procedure :: makeCalendar
    procedure :: print_unit

  end type ParameterizedEarthOrbit

  real (kind=DP), parameter :: PI = 2*asin(1.d0)

  interface ParameterizedEarthOrbit
     module procedure newParameterizedEarthOrbit
  end interface ParameterizedEarthOrbit

contains

   function newParameterizedEarthOrbit(yearBeforePresent) result(orbit)
      use BaseTime_mod
      use TimeInterval_mod
      use Rational_mod
      type (ParameterizedEarthOrbit) :: orbit
      integer, intent(in) :: yearBeforePresent
      

      call orbit%setMeanDistance(1.0_WP) ! 1 A.U.
      orbit%timeAtVernalEquinox = newBaseTime(Rational(79*24+12)*3600)
      call orbit%setTimeAtPeriapsis(newBaseTime(Rational(2*24+5)*3600))

      call orbit%setSiderealOrbitalPeriod(TimeInterval(Rational(365*24*3600)))
      call orbit%setSiderealRotationPeriod(TimeInterval(Rational(24*3600 * 365,366)))

      orbit%yearBeforePresent = yearBeforePresent

   end function newParameterizedEarthOrbit

  
  subroutine setYear(this, year)
    class (ParameterizedEarthOrbit), intent(inout) :: this
    real(kind=WP), intent(in) :: year

    real(kind=WP) :: pYear

    pYear = year - this%yearBeforePresent
    call orbpar(pYear, this%eccentricity, this%obliquity, this%longitudeAtPeriapsis)

    if (this%getVerbose()) then
       write(6,*) 'Set orbital parameters for year ',pyear,' (CE)'
       if (this%yearBeforePresent /= 0) write(6,*) 'offset by', &
            &  this%yearBeforePresent,' years from model year'
       write(6,*) "   Eccentricity: ", this%getEccentricity()
       write(6,*) "   Obliquity (degs): ",this%getObliquity()
       write(6,*) "   Precession (degs from ve): ", &
            &         this%getLongitudeAtPeriapsis()
    end if

  end subroutine setYear
  

  function getEccentricity(this) result(eccentricity)
    real(kind=WP) :: eccentricity
    class (ParameterizedEarthOrbit), intent(in) :: this
    eccentricity = this%eccentricity
  end function getEccentricity


  function getObliquity(this) result(obliquity)
    real(kind=WP) :: obliquity
    class (ParameterizedEarthOrbit), intent(in) :: this
    obliquity = this%obliquity
  end function getObliquity


  function getLongitudeAtPeriapsis(this) result(longitudeAtPeriapsis)
    real(kind=WP) :: longitudeAtPeriapsis
    class (ParameterizedEarthOrbit), intent(in) :: this
    longitudeAtPeriapsis = this%longitudeAtPeriapsis
  end function getLongitudeAtPeriapsis


  function getMeanAnomaly(this, t) result(meanAnomaly)
    use BaseTime_mod, only: BaseTime
    use Rational_mod
    use TimeInterval_mod
    real (kind=WP) :: meanAnomaly
    class (ParameterizedEarthOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    ! TODO: "fraction" is a bad name - Fortran intrinsic
    type (Rational) :: fraction
    type (TimeInterval) :: P

    P = this%getSiderealOrbitalPeriod()
    fraction = modulo(t,P) - modulo(this%timeAtPeriapsis,P)
    fraction = fraction / P

    meanAnomaly = real(fraction) * (2*PI)

  end function getMeanAnomaly


  function getTrueAnomaly(this, t) result(trueAnomaly)
    use BaseTime_mod
    use OrbitUtilities_mod, only: computeTrueAnomaly
    real (kind=WP) :: trueAnomaly
    class (ParameterizedEarthOrbit), intent(in) :: this
    class (BaseTime), intent(in) :: t

    real(kind=WP) :: meanAnomaly

    meanAnomaly = this%getMeanAnomaly(t)
    trueAnomaly = computeTrueAnomaly(meanAnomaly, this%getEccentricity())

  end function getTrueAnomaly


  subroutine setTimeAtPeriapsis(this, timeAtPeriapsis)
    use BaseTime_mod, only: BaseTime
    class (ParameterizedEarthOrbit), intent(inout) :: this
    type (BaseTime), intent(in) :: timeAtPeriapsis
    this%timeAtPeriapsis = timeAtPeriapsis
  end subroutine setTimeAtPeriapsis


  function getTimeAtPeriapsis(this) result(timeAtPeriapsis)
    use BaseTime_mod, only: BaseTime
    class (ParameterizedEarthOrbit), intent(in) :: this
    type (BaseTime) :: timeAtPeriapsis
    timeAtPeriapsis = this%timeAtPeriapsis
  end function getTimeAtPeriapsis


  function makeCalendar(this) result(calendar)
    use AbstractCalendar_mod, only: AbstractCalendar
    use JulianCalendar_mod, only: JulianCalendar
    class (AbstractCalendar), allocatable :: calendar
    class (ParameterizedEarthOrbit), intent(in) :: this

    type (BaseTime) :: vernalEquinox
    type (BaseTime) :: autumnalEquinox
    type (BaseTime) :: winterSolstice
    type (BaseTime) :: summerSolstice

    allocate(calendar, source=JulianCalendar())

    ! Add orbital dates

    vernalEquinox  = this%timeAtVernalEquinox
    summerSolstice = this%rotate(vernalEquinox, PI/2)
    autumnalEquinox = this%rotate(vernalEquinox, PI)
    winterSolstice = this%rotate(vernalEquinox, 3*PI/2)

    call calendar%addTransitionDate('vernal equinox', &
         & calendar%getCalendarDate(vernalEquinox))
    call calendar%addTransitionDate('autumnal equinox', &
         & calendar%getCalendarDate(autumnalEquinox))
    call calendar%addTransitionDate('winter solstice', &
         & calendar%getCalendarDate(winterSolstice))
    call calendar%addTransitionDate('summer solstice', &
         & calendar%getCalendarDate(summerSolstice))

  end function makeCalendar


  subroutine print_unit(this, unit)
    class (ParameterizedEarthOrbit), intent(in) :: this
    integer, intent(in) :: unit

    write(unit,*) 'Variable orbital parameters, updated each year. '

  end subroutine print_unit


end module ParameterizedEarthOrbit_mod
