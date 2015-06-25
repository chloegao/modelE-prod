#include "rundeck_opts.h"


















module AttributeDictionary_mod
  use Attributes_mod
  use AttributeHashMap_mod
  implicit none
  private

  public :: AttributeDictionary
  public :: newAttributeDictionary
  public :: assignment(=)  !, copyIt
  public :: toPointer
  public :: clean

  type, extends(AttributeHashMap) :: AttributeDictionary
    private
  contains
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: writeFormatted
    procedure :: equals
!!$    generic :: operator(==) => equals
     procedure :: insertInteger1dAttribute
 procedure :: insertIntegerAttribute
 procedure :: insertLogical1dAttribute
 procedure :: insertLogicalAttribute
 procedure :: insertRealDP1dAttribute
 procedure :: insertRealDPAttribute
 procedure :: insertString1dAttribute
 procedure :: insertStringAttribute
     generic :: insert => insertInteger1dAttribute
 generic :: insert => insertIntegerAttribute
 generic :: insert => insertLogical1dAttribute
 generic :: insert => insertLogicalAttribute
 generic :: insert => insertRealDP1dAttribute
 generic :: insert => insertRealDPAttribute
 generic :: insert => insertString1dAttribute
 generic :: insert => insertStringAttribute
  end type AttributeDictionary

  interface assignment(=)
    module procedure copyIt
  end interface assignment(=)

  ! TODO make type-bound when ifort is fixed - compilation slows to craws for higher-level modules
  interface operator(==)
    module procedure equals
  end interface operator(==)

  interface clean
    module procedure clean_
  end interface clean

  ! Prototypes used for distinguishing types - must include all possible attribute types
   type (Integer1dAttribute), target :: Integer1dAttributePrototype
 type (IntegerAttribute), target :: IntegerAttributePrototype
 type (Logical1dAttribute), target :: Logical1dAttributePrototype
 type (LogicalAttribute), target :: LogicalAttributePrototype
 type (RealDP1dAttribute), target :: RealDP1dAttributePrototype
 type (RealDPAttribute), target :: RealDPAttributePrototype
 type (String1dAttribute), target :: String1dAttributePrototype
 type (StringAttribute), target :: StringAttributePrototype


  integer, parameter :: NUM_TYPES = 8 
  type (AttributePointer) :: prototypes(NUM_TYPES)

contains

  function newAttributeDictionary() result(this)
    type (AttributeDictionary) :: this

    this%AttributeHashMap = newAttributeHashMap(100)
    
    ! Initialize the prototypes
    
    prototypes(1)%p => Integer1dAttributePrototype
   prototypes(2)%p => IntegerAttributePrototype
   prototypes(3)%p => Logical1dAttributePrototype
   prototypes(4)%p => LogicalAttributePrototype
   prototypes(5)%p => RealDP1dAttributePrototype
   prototypes(6)%p => RealDPAttributePrototype
   prototypes(7)%p => String1dAttributePrototype
   prototypes(8)%p => StringAttributePrototype
   
   

  end function newAttributeDictionary

  logical function equals(this, b)
    class (AttributeDictionary), intent(in) :: this
    class (AttributeDictionary), intent(in) :: b

    type (AttributeHashMapIterator) :: iter
    class (AbstractAttribute), pointer :: p1
!!$    type (AbstractAttributeReference) :: p2
    class (AbstractAttribute), pointer :: p2

    equals = .true.
    if (this%size() /= b%size()) then
      equals = .false. 
      print*,'different size',this%size(), b%size()
      return
    end if

    iter = this%begin()
    do while (iter /= this%last())

      if (.not. b%has(iter%key())) then
        equals = .false.
        print*,'different key'
        return
      end if
      
      p1 => iter%value()
      p2 => b%getReference(iter%key())

!!$      if (.not. (p1%equals(p2%ptr))) then
      if (.not. (p1%equals(p2))) then
        equals = .false.
        print*,'different value for key <',trim(iter%key()),'>'
        call p1%print()
        call p2%print()
        return
      end if

      call iter%next()

    end do

  end function equals

  ! 
  subroutine readUnformatted(this, unit)
    class (AttributeDictionary),intent(inout) :: this
    integer, intent(in) :: unit

    integer :: n
    integer :: i
    class (AbstractAttribute), pointer :: p, q
    character(len=MAX_LEN_KEY) :: key
    integer :: attributeType

    read(unit) n
    do i = 1, n
      read(unit) key
      read(unit) attributeType
      p => prototypes(attributeType)%p
      q => p%readUnformatted(unit)
      call this%insert(trim(key), q)
    end do

  end subroutine readUnformatted



  subroutine insertInteger1dAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value (:)

    call this%insert(key, newAttribute(value))
  end subroutine insertInteger1dAttribute

  subroutine insertIntegerAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    integer, intent(in) :: value 

    call this%insert(key, newAttribute(value))
  end subroutine insertIntegerAttribute

  subroutine insertLogical1dAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value (:)

    call this%insert(key, newAttribute(value))
  end subroutine insertLogical1dAttribute

  subroutine insertLogicalAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    logical, intent(in) :: value 

    call this%insert(key, newAttribute(value))
  end subroutine insertLogicalAttribute

  subroutine insertRealDP1dAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(kind=DP), intent(in) :: value (:)

    call this%insert(key, newAttribute(value))
  end subroutine insertRealDP1dAttribute

  subroutine insertRealDPAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    real(kind=DP), intent(in) :: value 

    call this%insert(key, newAttribute(value))
  end subroutine insertRealDPAttribute

  subroutine insertString1dAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value (:)

    call this%insert(key, newAttribute(value))
  end subroutine insertString1dAttribute

  subroutine insertStringAttribute(this, key, value)
    class (AttributeDictionary), intent(inout) :: this
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: value 

    call this%insert(key, newAttribute(value))
  end subroutine insertStringAttribute


  subroutine copyIt(a, b)
    type (AttributeDictionary), intent(out) :: a
    type (AttributeDictionary), intent(in) :: b

    a%AttributeHashMap = b%AttributeHashMap
  end subroutine copyIt

  subroutine writeUnformatted(this, unit)
    class (AttributeDictionary), intent(in) :: this
    integer, intent(in) :: unit

    type (AttributeHashMapIterator) :: iter
    class (AbstractAttribute), pointer :: p

    write(unit) this%size()
    iter = this%begin()
    do while (iter /= this%last())

      write(unit) iter%key()
      p => iter%value()
      write(unit) getAttributeType(p)
      call p%writeUnformatted(unit)

      call iter%next()
    end do

  contains

    integer function getAttributeType(p) result(attributeType)
      class (AbstractAttribute), pointer :: p
      integer :: i

      do i = 1, NUM_TYPES
        if (same_type_as(p, prototypes(i)%p)) then
          attributeType = i
          return
        end if
      end do
      
      call stop_model('No prototype for attribute in AttributeDictionary writeUnformatted.',255)

    end function getAttributeType

  end subroutine writeUnformatted

  subroutine writeFormatted(this, unit)
    class (AttributeDictionary), intent(in) :: this
    integer, intent(in) :: unit


  end subroutine writeFormatted

  subroutine clean_(this)
    type (AttributeDictionary), intent(inout) :: this
    call clean(this%AttributeHashMap) ! parent
  end subroutine clean_

end module AttributeDictionary_mod



