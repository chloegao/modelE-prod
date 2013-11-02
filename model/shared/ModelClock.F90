module ModelClock_mod
  use BaseTime_mod
  use Time_mod
  implicit none
  private

  public :: ModelClock
  public :: newModelClock

  type :: ModelClock
!!$    private
    type (Time) :: currentTime
    type (Time) :: startTime
    type (BaseTime) :: dt

    ! modelE legacy representation
    integer :: tick
    integer :: stepsPerDay

  contains
    procedure :: getCurrentTime
    procedure :: getTimeInSecondsFromDate
    procedure :: getAbsoluteTimeInSeconds
    procedure :: getTimeTick
    procedure :: getDt
    procedure :: isBeginningOfDay
    procedure :: nextTick
    procedure :: getDate
    procedure :: year
    procedure :: month
    procedure :: date
    procedure :: dayOfYear
    procedure :: hour
    procedure :: abbrev ! month abbreviation
  end type ModelClock

contains

  ! constructor
  function newModelClock(startTime, startTick, stepsPerDay) result(clock)
    use AbstractCalendar_mod
    type (ModelClock) :: clock
    type (Time), intent(in) :: startTime
    integer, intent(in) :: startTick
    integer, intent(in) :: stepsPerDay

    class (AbstractCalendar), pointer :: pCalendar

    clock%tick = startTick
    clock%stepsPerDay = stepsPerDay

    clock%startTime = startTime
    clock%currentTime = startTime

    pCalendar => startTime%calendar

    clock%dt = newBaseTime(pCalendar%getSecondsPerDay() / stepsPerDay)

  end function newModelClock

  subroutine nextTick(this)
    class (ModelClock), intent(inout) :: this

    this%tick = this%tick + 1
    call this%currentTime%setBaseTime(newBaseTime(this%currentTime + this%dt))

  end subroutine nextTick

  integer function getTimeTick(this)
    class (ModelClock), intent(in) :: this
    getTimeTick = this%tick
  end function getTimeTick

  type (BaseTime) function getDt(this) result(dt)
    class (ModelClock), intent(in) :: this
    dt = this%dt
  end function getDt

  logical function isBeginningOfDay(this)
    class (ModelClock), intent(in) :: this
    
    isBeginningOfDay = mod(this%tick, this%stepsPerDay) == 0
  end function isBeginningOfDay

  function getAbsoluteTimeInSeconds(this) result (secs)
    integer*8 :: secs
    class (ModelClock), intent(inout) :: this
    secs = this%currentTime%getWhole()
  end function getAbsoluteTimeInSeconds


  function getCurrentTime(this) result(t)
    type (Time) :: t
    class (ModelClock), intent(in) :: this
    t = this%currentTime
  end function getCurrentTime


  function getTimeInSecondsFromDate(this, year, month, date, hour) result (seconds)
    use AbstractCalendar_mod, only: AbstractCalendar
    use BaseTime_mod
    type (BaseTime) :: seconds
    class (ModelClock), intent(inout) :: this
    integer, intent(in) :: year, month, date, hour
    type (Time) :: aTime
    class (AbstractCalendar), pointer :: pCalendar

    pCalendar => this%currentTime%calendar
    aTime = newTime(pCalendar)
    call aTime%setByDate(year, month, date, hour)
    seconds = newBaseTime(this%currentTime - aTime)

  end function getTimeInSecondsFromDate


  subroutine getDate(this, year, month, dayOfYear, date, hour, amn)
!@sum  getDate gets Calendar info from internal timing info
!@auth Gavin Schmidt (updated by Tom CLune)
    use TimeConstants_mod, only: INT_SECONDS_PER_HOUR
    use JulianCalendar_mod, only: JULIAN_MONTHS
    use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION, CalendarMonth

    class (ModelClock), intent(in) :: this
    integer, optional, intent(out) :: year
    integer, optional, intent(out) :: month
    integer, optional, intent(out) :: dayOfYear
    integer, optional, intent(out) :: date
    integer, optional, intent(out) :: hour
    character(len=LEN_MONTH_ABBREVIATION), optional, intent(out) :: amn
    integer :: mnth

    if (present(year)) year = this%currentTime%getYear()
    if (present(dayOfYear)) dayOfYear = this%currentTime%getDayOfYear()
    if (present(month)) month = this%currentTime%getMonth()

    if (present(amn)) then
      mnth = this%currentTime%getMonth()
      amn = this%currentTime%getAbbreviation()
    end if

    if (present(date)) date = this%currentTime%getDate()
    if (present(hour)) hour = this%currentTime%getHour()

    return
  end subroutine getDate

  integer function year(this)
    class (ModelClock), intent(in) :: this
    year = this%currentTime%getYear()
  end function year

  integer function month(this)
    class (ModelClock), intent(in) :: this
    month = this%currentTime%getMonth()
  end function month

  integer function date(this)
    class (ModelClock), intent(in) :: this

    date = this%currentTime%getDate()
  end function date

  integer function dayOfYear(this)
    class (ModelClock), intent(in) :: this
    dayOfYear = this%currentTime%getDayOfYear()
  end function dayOfYear

  integer function hour(this)
    class (ModelClock), intent(in) :: this

    hour = this%currentTime%getHour()
  end function hour

  function abbrev(this)
    use CalendarMonth_mod, only: LEN_MONTH_ABBREVIATION
    character(len=LEN_MONTH_ABBREVIATION) abbrev
    class (ModelClock), intent(in) :: this

    abbrev = this%currentTime%getAbbreviation()
  end function abbrev

end module ModelClock_mod
