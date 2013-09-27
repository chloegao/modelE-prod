! Design loosely based upon that from ESMF to support time represented as
! rational numbers with large range and high accuracy.
! Rational number ensure that subcycles align perfectly with 
! larger cycles.  E.g. dt evenly divides day.

module Rational_mod
  implicit none
  private

  public :: Rational
  public :: nint

  integer, parameter :: DEFAULT_INT = kind(1)
  integer, parameter :: LONG_INT = selected_int_kind(10)
  integer, parameter :: LONG = max(DEFAULT_INT, LONG_INT) ! in case only 32 bit available
  integer, parameter :: SP = selected_real_kind(6)
  integer, parameter :: DP = selected_real_kind(14)

  ! Numerical value is whole + numerator/denominator
  type Rational
    integer(kind=LONG) :: numerator   = 0 
    integer(kind=LONG) :: denominator = 1 ! always positive
    integer(kind=LONG) :: whole       = 0
  contains
    procedure :: getWhole
!!$    procedure, pass(this) :: toReal_sp
!!$    procedure, pass(this) :: toReal_dp

    procedure :: add_fraction
    procedure :: subtract_fraction

    procedure :: multiply_fraction
    procedure, pass(a) :: multiply_int
    procedure, pass(a) :: multiply_int2

    procedure :: divide_fraction
    procedure, pass(a) :: divide_realdp
    procedure, pass(a) :: divide_int4
    procedure, pass(a) :: divide_int8

    procedure :: equals_fraction
    procedure, pass(a) :: equals_int
    procedure, pass(b) :: equals_int2
    procedure :: lessThan_fraction
    procedure :: greaterThan_fraction
    procedure :: convertToReal

    generic :: operator(+) => add_fraction
    generic :: operator(-) => subtract_fraction
    generic :: operator(*) => multiply_fraction, multiply_int, multiply_int2
    generic :: operator(/) => divide_fraction, divide_realDP, divide_int4, divide_int8
    generic :: operator(==) => equals_fraction, equals_int, equals_int2
    generic :: operator(<) => lessThan_fraction
    generic :: operator(>) => greaterThan_fraction
!!$    generic :: assignment(=) => toReal_sp, toReal_dp
    procedure, private :: reduce

    procedure :: print

  end type Rational

  interface Rational
    module procedure newRational_defaultWhole
    module procedure newRational_default_n_over_d
    module procedure newRational_default
#ifndef LONG_UNSUPPORTED
    module procedure newRational_longWhole
    module procedure newRational_long_n_over_d
    module procedure newRational_long
#endif

    module procedure newRational_real_sp
    module procedure newRational_real_dp
  end interface Rational

  interface nint
    module procedure nintRational
  end interface nint

contains

  integer(kind=LONG) function getWhole(this) result(whole)
    class (Rational), intent(in) :: this
    whole = this%whole
  end function getWhole

  function convertToReal(this) result(x)
    real(kind=DP) :: x
    class(Rational), intent(in) :: this
    x = real(this%numerator,kind=DP) / real(this%denominator,kind=DP)
    x = x + this%whole
  end function convertToReal

  subroutine toReal_sp(x, this)
    real(kind=SP), intent(out) :: x
    class (Rational), intent(in) :: this

    real(kind=DP) :: x_dp

    x_dp = this%convertToReal()
    x = x_dp

  end subroutine toReal_sp

  subroutine toReal_dp(x, this)
    real(kind=dp), intent(out) :: x
    class (Rational), intent(in) :: this

    x = this%whole + real(this%numerator,kind=dp)/this%denominator
  end subroutine toReal_dp

  ! Return the nearest integer
  ! Rounds to even for r = n + m/2
  integer function nintRational(this) result(n)
    class (Rational), intent(in) :: this

    n = abs(this%whole)

    if (2*abs(this%numerator) > this%denominator) then
      n = n + 1
    else if (2*abs(this%numerator) == this%denominator) then
      if (mod(n,2) == 1) n = n + 1
    end if
    
    if (this%whole < 0 .or. this%numerator < 0) n = - n

  end function nintRational

  function newRational_defaultWhole(whole) result(r)
    type (Rational) :: r
    integer, intent(in) :: whole

    r%whole = whole
    r%numerator = 0
    r%denominator = 1
  end function newRational_defaultWhole

  function newRational_default_n_over_d(numerator, denominator) result(r)
    type (Rational) :: r
    integer, intent(in) :: numerator
    integer, intent(in) :: denominator

    r%whole = 0
    r%numerator = numerator
    r%denominator = denominator

    call r%reduce()

  end function newRational_default_n_over_d

  function newRational_default(whole, numerator, denominator) result(r)
    type (Rational) :: r
    integer, intent(in) :: whole
    integer, intent(in) :: numerator
    integer, intent(in) :: denominator

    ! TODO: Ifort - 14.0 does not use the correct interface for other part of constructor
    r = Rational(whole) + Rational(numerator, denominator)

 end function newRational_default

#ifndef LONG_UNSUPPORTED

  function newRational_longWhole(whole) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: whole

    r%whole = whole
    r%numerator = 0
    r%denominator = 1
  end function newRational_longWhole

  function newRational_long_n_over_d(numerator, denominator) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: numerator
    integer(kind=LONG), intent(in) :: denominator

    r%whole = 0
    r%numerator = numerator
    r%denominator = denominator
    call r%reduce()
  end function newRational_long_n_over_d

  function newRational_long(whole, numerator, denominator) result(r)
    type (Rational) :: r
    integer(kind=LONG), intent(in) :: whole
    integer(kind=LONG), intent(in) :: numerator
    integer(kind=LONG), intent(in) :: denominator

    r = Rational(whole) + Rational(numerator, denominator)
 end function newRational_long
#endif

! Use continued fractions to convert floating point to a fraction within
! a specified tolerance.
! http://en.wikipedia.org/wiki/Continued_fraction
! http://www-math.mit.edu/phase2/UJM/vol1/COLLIN~1.PDF

  function newRational_real_sp(x, tolerance) result(r)
     type (Rational) :: r
     real(kind=SP), intent(in) :: x
     real(kind=SP), intent(in) :: tolerance

     r = Rational(real(x,kind=DP), real(tolerance,kind=DP))
  end function newRational_real_sp

  function newRational_real_dp(x, tolerance) result(r)
     type (Rational) :: r
     real(kind=DP), intent(in) :: x
     real(kind=DP), intent(in) :: tolerance

     real(kind=DP) :: xx, absx
     integer(kind=LONG) :: a, w
     integer(kind=LONG) :: p_n, p_nm1, p_nm2
     integer(kind=LONG) :: q_n, q_nm1, q_nm2

     absx = abs(x)

     xx = absx
     w = floor(xx)
     xx = xx - w

     p_nm2 = 1
     q_nm2 = 0

     p_nm1 = 0
     q_nm1 = 1

     p_n = 0
     q_n = 1

     do while (xx /= 0.d0)
       xx = 1/xx
       a = floor(xx)
       xx = xx - a

        p_n = a * p_nm1 + p_nm2
        q_n = a * q_nm1 + q_nm2
        
        p_nm2 = p_nm1
        q_nm2 = q_nm1
        
        p_nm1 = p_n
        q_nm1 = q_n
        
        if (abs(q_n*(absx-w) - p_n) < q_n*tolerance) exit

     end do

     ! Restore sign of result
     if (x < 0) then
        w = -w
        p_n = -p_n
     end if

     r = newRational_long(w, p_n, q_n)

  end function newRational_real_dp

! Add two fractions and reduce to simplest form.
  function add_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    integer(kind=LONG) :: gcf

    c%whole = a%whole + b%whole 

    gcf = greatestCommonFactor(a%denominator, b%denominator)
    
    c%numerator= a%numerator*(b%denominator/gcf) + b%numerator*(a%denominator/gcf)
    c%denominator = (a%denominator / gcf) * b%denominator

    call c%reduce()

  end function add_fraction

! Subtract two fractions and reduce to simplest form.
! a - b = a + (- b)
  function subtract_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    c = a + b*(-1)

  end function subtract_fraction

! Multiply two fractions and reduce to simplest form.
  function multiply_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    c%whole = a%whole * b%whole 
    c%numerator = b%whole*a%numerator*b%denominator + &
         & a%whole*b%numerator*a%denominator + &
         & a%numerator * b%numerator
    c%denominator = a%denominator*b%denominator
    call c%reduce()

  end function multiply_fraction

! Multiply by int on right
  function multiply_int(a, i) result(c)
    class (Rational), intent(in) :: a
    integer, intent(in) :: i
    type (Rational) :: c

    c%whole = a%whole * i
    c%numerator = a%numerator * i
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_int

! Multiply by int on left
  function multiply_int2(i, a) result(c)
    integer, intent(in) :: i
    class (Rational), intent(in) :: a
    type (Rational) :: c

    c%whole = i * a%whole
    c%numerator = i * a%numerator
    c%denominator = a%denominator
    call c%reduce()

  end function multiply_int2

! Divide two fractions and reduce to simplest form.
  function divide_fraction(a, b) result(c)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b
    type (Rational) :: c

    c = Rational(a%whole*a%denominator + a%numerator, &
         & b%whole*b%denominator + b%numerator) * &
         & Rational(b%denominator, a%denominator)
    call c%reduce()

  end function divide_fraction

! Divide fractions by an double and reduce
  function divide_realDP(a, x) result(c)
    class (Rational), intent(in) :: a
    real(kind=DP), intent(in) :: x
    type (Rational) :: c

    c = a / Rational(x)
    call c%reduce()

  end function divide_realDP

! Divide fractions by an integer and reduce
  function divide_int4(a, i) result(c)
    class (Rational), intent(in) :: a
    integer(kind=DEFAULT_INT), intent(in) :: i
    type (Rational) :: c

    c = a / Rational(i)
    call c%reduce()

  end function divide_int4

! Divide fractions by an integer and reduce
  function divide_int8(a, i) result(c)
    class (Rational), intent(in) :: a
    integer(kind=LONG), intent(in) :: i
    type (Rational) :: c

    c = a / Rational(i)
    call c%reduce()

  end function divide_int8

  logical function equals_fraction(a, b) result(equals)
    class (Rational), intent(in) :: a
    class (Rational), intent(in) :: b

    equals = (a%whole == b%whole) .and. &
         & (a%numerator == b%numerator) .and. &
         & (a%denominator == b%denominator)

  end function equals_fraction

  logical function equals_int(a, b) result(equals)
    class (Rational), intent(in) :: a
    integer, intent(in) :: b

    equals = (a%whole == b) .and. (a%numerator == 0)

  end function equals_int

  logical function equals_int2(a, b) result(equals)
    integer, intent(in) :: a
    class (Rational), intent(in) :: b

    equals = (b%whole == a) .and. (b%numerator == 0)

  end function equals_int2

  logical function lessThan_fraction(r1, r2) result(lessThan)
    class (Rational), intent(in) :: r1
    class (Rational), intent(in) :: r2

    if (r1%whole < r2%whole) then
      lessThan = .true.
    else if (r1%whole > r2%whole) then
      lessThan = .false.
    else
      if (r1%numerator*r2%denominator < r2%numerator*r1%denominator) then
        lessThan = .true.
      else
        lessThan = .false.
      end if
    end if

  end function lessThan_fraction

  logical function greaterThan_fraction(r1, r2) result(greaterThan)
    class (Rational), intent(in) :: r1
    class (Rational), intent(in) :: r2

    if (r1%whole > r2%whole) then
      greaterThan = .true.
    else if (r1%whole < r2%whole) then
      greaterThan = .false.
    else
      if (r1%numerator*r2%denominator > r2%numerator*r1%denominator) then
        greaterThan = .true.
      else
        greaterThan = .false.
      end if
    end if

  end function greaterThan_fraction

  !
  ! Reduce rational number to standard form:
  !    -  |numerator| < denominator
  !    - sign(whole)*sign(numerator) >= 0
  !    - gcd(numerotor,denominator) = 1
  !
  subroutine reduce(this)
    class (Rational), intent(inout) :: this

    integer(kind=Long) :: whole, factor

    associate(w => this%whole, n => this%numerator, d => this%denominator)

      if (d < 0) then
         d = -d
         n = -n
      end if

      if (w > 0 .and. n < 0) then
         w = w - 1
         n = n + d
      else if (w < 0 .and. n > 0) then
         w =  w + 1
         n = n - d
      end if
         
      if (abs(n) >= d) then
         whole = n / d
         w = w + whole
         n = n - whole*d
      end if

      factor = greatestCommonFactor(abs(n), d)
      n = n / factor
      d = d / factor

    end associate

  end subroutine reduce

! Citation - Euclid
  function greatestCommonFactor(a, b) result(factor)
     integer(kind=LONG) :: factor
     integer(kind=LONG), intent(in) :: a
     integer(kind=LONG), intent(in) :: b

     integer(kind=LONG) :: ia, ib, f
     
     if (a==0 .and. b==0) then
        factor = 1
        return
     end if

     ia = a
     ib = b
     f = ia

     do while (f > 0)
        ia = mod(ib, f)
        ib = f
        f = ia
     end do

     factor = ib

  end function greatestCommonFactor

  subroutine print(this)
     class (Rational), intent(in) :: this
     write(*,'(a,i0," + ",i0,"/",i0)') 'Rational: ',this%whole, this%numerator, this%denominator
  end subroutine print

end module Rational_mod
