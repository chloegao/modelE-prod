module PlanetCalendar_mod
!@sum Parameterized calendar that allows for independent day length and year length.
!@+   Years are assumed to be integer number of days.  
!@+   Days are assumed to be an integer number of seconds.
!@+   All months have 30 days, except the last which is <= 30 days.
!@+   Time of day is given in h/m/s.   Note that a day may not have an integral number 
!@+ of hours/minutes.

!@auth T. Clune
  use Calendar_mod
  use Month_mod
  use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
  use TimeConstants_mod, only: INT_HOURS_PER_DAY
  use BaseTime_mod
  implicit none
  private 

  public :: PlanetCalendar
  public :: makePlanetCalendar
  
  integer, parameter :: LONG = selected_int_kind(15)
  type, extends(Calendar) :: PlanetCalendar
    private 
    type (BaseTime) :: secondsPerDay
    type (BaseTime) :: secondsPerPlanetHour
    type (BaseTime) :: secondsPerYear
    integer :: daysPerYear
    integer :: lastDayOfMonth(0:INT_MONTHS_PER_YEAR)
  contains
    procedure :: getDayOfYear
    procedure :: getYear
    procedure :: getMonth
    procedure :: getDate
    procedure :: getHour
    procedure :: getAbbreviation
    procedure :: convertToTime
    procedure :: getSecondsPerDay
    procedure :: getDaysPerYear
    procedure :: getSecondsPerHour
    procedure :: getDaysPerMonth
    procedure :: getLastDayOfMonth
    procedure :: getMidDayOfMonth
  end type PlanetCalendar

  integer, parameter :: BASE_YEAR = 1 ! there was no year "0"
  
  type (PlanetCalendar), save, target :: singletonPlanetCalendar


contains

  ! Returns the singleton instance of the Planet calendar.
  ! It does not make sense to have multiple Planet calendars.
  function makePlanetCalendar(secondsPerDay, daysPerYear) result(ptr)
    use Rational_mod
    class (Calendar), pointer :: ptr
    type (BaseTime), intent(in) :: secondsPerDay
    integer, intent(in) :: daysPerYear

    integer :: i
    type (BaseTime) :: secondsPerPlanetHour

    
    singletonPlanetCalendar%secondsPerDay = secondsPerDay
    singletonPlanetCalendar%daysPerYear = daysPerYear

    singletonPlanetCalendar%secondsPerPlanetHour = newBaseTime((secondsPerDay / INT_HOURS_PER_DAY))
    singletonPlanetCalendar%secondsPerYear = newBaseTime(secondsPerDay * daysPerYear)
    
    do i = 0, INT_MONTHS_PER_YEAR
       singletonPlanetCalendar%lastDayOfMonth(i) = &
            & i*daysPerYear/INT_MONTHS_PER_YEAR
    end do
    ptr => singletonPlanetCalendar
    
  end function makePlanetCalendar

  integer function getYear(this, t)
    use Time_mod
    use Rational_mod
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    integer(kind=LONG) :: tSeconds
    type (Rational) :: years
    
    years = t / this%secondsPerYear
    getYear = BASE_YEAR + years%getWhole()

  end function getYear

  integer function getDayOfYear(this, t)
    use Rational_mod
    use Time_mod, only: Time
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t
    
    type (Rational) :: secondsIntoYear
    type (Rational) :: daysIntoYear

    secondsIntoYear = t - this%secondsPerYear * (this%getYear(t) - BASE_YEAR)

    daysIntoYear = secondsIntoYear / this%secondsPerDay
    getDayOfYear = 1 + daysIntoYear%getWhole()

  end function getDayOfYear

  integer function getMonth(this, t) result(month)
    use Time_mod
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t
    
    integer :: day
    integer :: i

    day = this%getDayOfYear(t)
    do i = 1, INT_MONTHS_PER_YEAR
       if (day <= this%lastDayOfMonth(i)) then
          month = i
          exit
       end if
    end do

  end function getMonth

  function getAbbreviation(this, t) result(abbrev)
    use Month_mod, only: LEN_MONTH_ABBREVIATION
    use Time_mod
!!$$    character(len=LEN_MONTH_ABBREVIATION) :: abbrev
    !TODO workaround for NAG - needs literal here
    character(len=4) :: abbrev
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    select case(this%getMonth(t))
    case (1)
      abbrev = 'JAN '
    case (2)
      abbrev = 'FEB '
    case (3)
      abbrev = 'MAR '
    case (4)
      abbrev = 'APR '
    case (5)
      abbrev = 'MAY '
    case (6)
      abbrev = 'JUN '
    case (7)
      abbrev = 'JUL '
    case (8)
      abbrev = 'AUG '
    case (9)
      abbrev = 'SEP '
    case (10)
      abbrev = 'OCT '
    case (11)
      abbrev = 'NOV '
    case (12)
      abbrev = 'DEC '
    end select


  end function getAbbreviation

  integer function getDate(this, t) result(date)
    use Time_mod
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    integer :: month
    
    month = this%getMonth(t)
    date = this%getDayOfYear(t) - this%lastDayOfMonth(month-1)

  end function getDate

  integer function getHour(this, t) result(hour)
    use Time_mod
    use Rational_mod
    class (PlanetCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    type (BaseTime) :: timeOfDay
    type (Rational) :: hoursIntoDay
    
    timeOfDay = this%getTimeOfDay(t)
    hoursIntoDay = (timeOfDay / this%secondsPerPlanetHour)
    hour = hoursIntoDay%getWhole()

 end function getHour

  function convertToTime(this, year, month, date, hour) result(t)
    use Time_mod
    use Rational_mod
    class (PlanetCalendar), intent(in) :: this
    integer, intent(in) :: year, month, date, hour
    type (BaseTime) :: t

    integer :: numYears, numDays
    type (Rational) :: numSeconds
    
    numYears = year - BASE_YEAR
    numDays = numYears*this%daysPerYear + this%lastDayOfMonth(month-1) + (date-1)
    numSeconds = (this%secondsPerPlanetHour*hour) + (this%secondsPerDay*numDays)

    t = newBaseTime(numSeconds)

  end function convertToTime

  function getSecondsPerDay(this) result(secondsPerDay)
    type (BaseTime) :: secondsPerDay
    class (PlanetCalendar), intent(in) :: this
    secondsPerDay = this%secondsPerDay
  end function getSecondsPerDay

  integer function getDaysPerMonth(this, month) result(days)
    class (PlanetCalendar), intent(in) :: this
    integer, intent(in) :: month

    days = this%lastDayOfMonth(month) - this%lastDayOfMonth(month-1)
  end function getDaysPerMonth

  integer function getLastDayOfMonth(this, month) result(lastDay)
    class (PlanetCalendar), intent(in) :: this
    integer, intent(in) :: month
    lastDay = this%lastDayOfMonth(month)
  end function getLastDayOfMonth

  integer function getMidDayOfMonth(this, month) result(midDay)
    class (PlanetCalendar), intent(in) :: this
    integer, intent(in) :: month

    select case(month)
    case (0)
       midDay = (this%lastDayOfMonth(month-1) - this%lastDayOfMonth(month-1))/2 
    case (1:INT_MONTHS_PER_YEAR)
       midDay = (this%lastDayOfMonth(month) + this%lastDayOfMonth(month-1))/2
    end select
  end function getMidDayOfMonth

  integer function getDaysPerYear(this) result(days)
    class (PlanetCalendar), intent(in) :: this
    days = this%daysPerYear
  end function getDaysPerYear
 
  function getSecondsPerHour(this) result(numSeconds)
    type (BaseTime) :: numSeconds
    class (PlanetCalendar), intent(in) :: this
    numSeconds = this%secondsPerPlanetHour
  end function getSecondsPerHour

end module PlanetCalendar_mod
