! These constants are used by Clock/Calendar to manage time within the model.
! For the moment, legacy names are preserved until the model is adapted to the newer
! names.
! 
! Note that many of these "constants" are actually dependent on choice of calendar.  E.g. 
! Number of days per year.   The constants placed here correspond to the Julian calendar
! that is the primary case for the model.
!

module TimeConstants_mod
  implicit none
  public  ! just a bunch of named constants

  ! Constants that are assumed for _all_ planets
  integer, parameter :: INT_HOURS_PER_DAY = 24
  integer, parameter :: INT_MONTHS_PER_YEAR = 12

  ! Earth specific values (Julian calendar)
  integer, parameter :: INT_SECONDS_PER_MINUTE = 60
  integer, parameter :: INT_MINUTES_PER_HOUR = 60
  integer, parameter :: INT_DAYS_PER_YEAR = 365

  integer, parameter :: INT_SECONDS_PER_HOUR = INT_SECONDS_PER_MINUTE * INT_MINUTES_PER_HOUR
  integer, parameter :: INT_SECONDS_PER_DAY = INT_SECONDS_PER_HOUR * INT_HOURS_PER_DAY
  integer, parameter :: INT_MINUTES_PER_DAY = INT_MINUTES_PER_HOUR * INT_HOURS_PER_DAY
  integer, parameter :: INT_SECONDS_PER_YEAR = INT_DAYS_PER_YEAR * INT_SECONDS_PER_DAY

  real*8, parameter :: SECONDS_PER_MINUTE = INT_SECONDS_PER_MINUTE
  real*8, parameter :: SECONDS_PER_HOUR = INT_SECONDS_PER_HOUR
  real*8, parameter :: SECONDS_PER_DAY = INT_SECONDS_PER_DAY
  real*8, parameter :: SECONDS_PER_YEAR = SECONDS_PER_DAY * INT_DAYS_PER_YEAR

  real*8, parameter :: MINUTES_PER_HOUR = INT_MINUTES_PER_HOUR
  real*8, parameter :: HOURS_PER_DAY = INT_HOURS_PER_DAY
  real*8, parameter :: DAYS_PER_YEAR = INT_DAYS_PER_YEAR

  real*8, parameter :: EARTH_DAYS_PER_DAY = 1.
  real*8, parameter :: EARTH_DAYS_PER_YEAR = 365.

! Deprecated names
!-----------------

!@var JDPERY  number of days per year
  integer, parameter :: JDPERY = INT_DAYS_PER_YEAR
!@var JMperY number of months per year
  integer, parameter :: JMPERY = INT_MONTHS_PER_YEAR
!@var SDAY seconds per day
  real*8, parameter :: SDAY = SECONDS_PER_DAY
!@var SYR seconds per year
  real*8, parameter :: SYR = SECONDS_PER_YEAR
!@var HRDAY hours per day
  real*8, parameter :: HRDAY = HOURS_PER_DAY
!@var EDPERD Earth days per day
  real*8, parameter :: EDPERD = EARTH_DAYS_PER_DAY
!@var EDPERY Earth days per year
  real*8, parameter :: EDPERY = EARTH_DAYS_PER_YEAR

end module TimeConstants_mod
