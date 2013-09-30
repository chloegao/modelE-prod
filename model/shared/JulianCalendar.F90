! During the transition period, this module will continue to support the legacy
! names for constants which tend towards the terse end of the spectrum.  

module JulianCalendar_mod
!@sum Specifies the parameters for the Julian calendar used in modelE.
!@auth T. Clune
  use Calendar_mod
  use Month_mod
  use TimeConstants_mod, only: INT_DAYS_PER_YEAR
  use TimeConstants_mod, only: INT_SECONDS_PER_DAY
  use TimeConstants_mod, only: INT_SECONDS_PER_HOUR
  use Rational_mod
  use BaseTime_mod
  implicit none
  private 

  public :: JulianCalendar
  public :: makeJulianCalendar
  
  type, extends(Calendar) :: JulianCalendar
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
  end type JulianCalendar

  !         Legacy :  New
  public :: JDendOfM, LAST_JULIAN_DAY_IN_MONTH
  public :: JDmidOfM, MID_JULIAN_DAY_IN_MONTH

  ! Not entirely clear if these constants belong here, since they actually
  ! should vary for other planets and/or eras (e.g. paleo has shorter day)

  public :: Month_type, JULIAN_MONTHS
  public :: JANUARY,   FEBRUARY, MARCH,    APRIL
  public :: MAY,       JUNE,     JULY,     AUGUST
  public :: SEPTEMBER, OCTOBER,  NOVEMBER, DECEMBER

  
!@var DAYS_PER_YEAR (JDPERY)    number of days per year
  integer, parameter :: DAYS_PER_YEAR = 365, JDPERY = DAYS_PER_YEAR
!@var MONTHS_PER_YEAR (JMperY)  number of months per year
  integer, parameter :: MONTHS_PER_YEAR = 12, JMPERY = MONTHS_PER_YEAR
!@var INT_SECONDS_PER_YEAR (JMperY)  number of seconds per year (long integer to be safe for arithmetic)
  integer*8, parameter :: INT_SECONDS_PER_YEAR = DAYS_PER_YEAR * INT_SECONDS_PER_DAY


!@var LAST_JULIAN_DAY_IN_MONTH (JDendOfM, ) last Julian day in month
  integer, parameter :: LAST_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR) = (/ &
       & 0,31,59,90,120,151,181,212,243,273,304,334,365 &
       & /)
  integer, parameter :: JDendOfM(0:MONTHS_PER_YEAR) = LAST_JULIAN_DAY_IN_MONTH
!@var MID_JULIAN_DAY_IN_MONTH(0:13) (JDmidOfM(0:13)) middle Julian day in month
  integer, parameter :: MID_JULIAN_DAY_IN_MONTH(0:MONTHS_PER_YEAR+1) = (/ &
       & -15,16,45,75,106,136,167,197,228,259,289,320,350,381 &
       & /)
  integer, parameter :: JDmidOfM(0:MONTHS_PER_YEAR+1) = MID_JULIAN_DAY_IN_MONTH

  ! Months
  type (Month_type), parameter :: JANUARY   = Month_type('JAN ', 'January   ', 31,  31,  16)
  type (Month_type), parameter :: FEBRUARY  = Month_type('FEB ', 'February  ', 28,  59,  45)
  type (Month_type), parameter :: MARCH     = Month_type('MAR ', 'March     ', 31,  90,  75)
  type (Month_type), parameter :: APRIL     = Month_type('APR ', 'April     ', 30, 120,  106)
  type (Month_type), parameter :: MAY       = Month_type('MAY ', 'May       ', 31, 151,  136)
  type (Month_type), parameter :: JUNE      = Month_type('JUN ', 'June      ', 30, 181,  167)
  type (Month_type), parameter :: JULY      = Month_type('JUL ', 'July      ', 31, 212,  197)
  type (Month_type), parameter :: AUGUST    = Month_type('AUG ', 'August    ', 31, 243,  228)
  type (Month_type), parameter :: SEPTEMBER = Month_type('SEP ', 'September ', 30, 273,  259)
  type (Month_type), parameter :: OCTOBER   = Month_type('OCT ', 'October   ', 31, 304,  289)
  type (Month_type), parameter :: NOVEMBER  = Month_type('NOV ', 'November  ', 30, 334,  320)
  type (Month_type), parameter :: DECEMBER  = Month_type('DEC ', 'December  ', 31, 365,  350)

  type (Month_type), parameter :: JULIAN_MONTHS(MONTHS_PER_YEAR) = [ &
       & JANUARY,  FEBRUARY, MARCH,     &
       & APRIL,    MAY,      JUNE,      &
       & JULY,     AUGUST,   SEPTEMBER, &
       & OCTOBER,  NOVEMBER, DECEMBER   &
       & ]

!@var AMONTH(0:12)  (3-4 letter) names for all months
! AMONTH(0) = 'IC' (i.e. initial conditions) only used early in a
! model run.  Should find a way to eliminate it.
  character(len=LEN_MONTH_ABBREVIATION), parameter :: AMONTH(0:12) = [ 'IC  ', JULIAN_MONTHS%abbreviation ]

  integer, parameter :: BASE_YEAR = 1 ! there was no year "0"
  type (JulianCalendar), save, target :: singletonJulianCalendar

contains

  ! Returns the singleton instance of the Julian calendar.
  ! It does not make sense to have multiple Julian calendars.
  function makeJulianCalendar() result(ptr)
    class (Calendar), pointer :: ptr
    ptr => singletonJulianCalendar
  end function makeJulianCalendar

  integer function getYear(this, t)
    use Time_mod
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    integer*8 :: tSeconds

    tSeconds = t%getWhole()
    getYear = BASE_YEAR + (tSeconds / INT_SECONDS_PER_YEAR)

  end function getYear

  integer function getDayOfYear(this, t)
    use Time_mod, only: Time
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t
    
    integer*8 :: tSecondsInYear

    tSecondsInYear = t%getWhole() - INT_SECONDS_PER_YEAR * (this%getYear(t) - BASE_YEAR)
    getDayOfYear = 1 + (tSecondsInYear / INT_SECONDS_PER_DAY)

  end function getDayOfYear

  integer function getMonth(this, t) result(month)
    use Time_mod
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t
    
    integer :: day

    day = this%getDayOfYear(t)
    month = 1
    do while (day > LAST_JULIAN_DAY_IN_MONTH(month))
      month = month + 1
    end do

  end function getMonth

  function getAbbreviation(this, t) result(abbrev)
    use Month_mod, only: LEN_MONTH_ABBREVIATION
    use Time_mod
!!$$    character(len=LEN_MONTH_ABBREVIATION) :: abbrev
    !TODO workaround for NAG - needs literal here
    character(len=4) :: abbrev
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    integer :: month

    month = this%getMonth(t)
    abbrev = JULIAN_MONTHS(month)%abbreviation

  end function getAbbreviation

  integer function getDate(this, t) result(date)
    use Time_mod
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    date = this%getDayOfYear(t) - LAST_JULIAN_DAY_IN_MONTH(this%getMonth(t)-1)
  end function getDate

  integer function getHour(this, t) result(hour)
    use Time_mod
    class (JulianCalendar), intent(in) :: this
    class (Time), intent(in) :: t

    type (BaseTime) :: timeOfDay

    timeOfDay = this%getTimeOfDay(t)
    hour = timeOfDay%getWhole()/INT_SECONDS_PER_HOUR

 end function getHour

  function convertToTime(this, year, month, date, hour) result(t)
    use Time_mod
    use TimeConstants_mod, only: INT_SECONDS_PER_HOUR, INT_HOURS_PER_DAY, INT_DAYS_PER_YEAR
    class (JulianCalendar), intent(in) :: this
    integer, intent(in) :: year, month, date, hour
    type (BaseTime) :: t

    integer*8 :: numYears, numDays, numHours, numSeconds

    numYears = year - BASE_YEAR
    numDays = numYears * INT_DAYS_PER_YEAR + LAST_JULIAN_DAY_IN_MONTH(month-1) + (date-1)
    numHours = numDays * INT_HOURS_PER_DAY + hour
    numSeconds = numHours * INT_SECONDS_PER_HOUR
    t = newBaseTime(numSeconds)

  end function convertToTime

  type (BaseTime) function getSecondsPerDay(this) result(secondsPerDay)
    use Rational_mod
    class (JulianCalendar), intent(in) :: this 
    secondsPerDay = newBaseTime(INT_SECONDS_PER_DAY)
  end function getSecondsPerDay

  integer function getDaysPerMonth(this, month) result(days)
    class (JulianCalendar), intent(in) :: this
    integer, intent(in) :: month
    days = JULIAN_MONTHS(month)%numDays
  end function getDaysPerMonth

  integer function getLastDayOfMonth(this, month) result(lastDay)
    class (JulianCalendar), intent(in) :: this
    integer, intent(in) :: month
    lastDay = JDendOfM(month)
  end function getLastDayOfMonth

  integer function getMidDayOfMonth(this, month) result(midDay)
    class (JulianCalendar), intent(in) :: this
    integer, intent(in) :: month
    midDay = JDmidOfM(month)
  end function getMidDayOfMonth

  integer function getDaysPerYear(this) result(numDays)
    class (JulianCalendar), intent(in) :: this
    numDays = INT_DAYS_PER_YEAR
  end function getDaysPerYear

  function getSecondsPerHour(this) result(numSeconds)
    type (BaseTime) :: numSeconds
    class (JulianCalendar), intent(in) :: this
    numSeconds = newBaseTime(INT_SECONDS_PER_HOUR)
  end function getSecondsPerHour
  
end module JulianCalendar_mod
