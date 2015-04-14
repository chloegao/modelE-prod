#include "rundeck_opts.h"




















module integerAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  
  implicit none
  private

#define TYPE integerAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    integer :: value
  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    integer, intent(in) :: value 
    

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     integer , intent(inout) :: value 
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)
        value = entry%value



    class default
      call throwException('Illegal conversion of integerAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    integer, pointer :: ptr 
    class (AbstractAttribute), target, intent(in) :: entry
    integer :: cast 

    select type (q => entry)
    type is (integerAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of integerAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (integerAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (integerAttribute)
      if ((this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (integerAttribute), intent(in) :: this
    print*,'  Type:  ', 'integerAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (integerAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = toStringElemental(this%value)

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (integerAttribute), intent(in) :: this
    integer, intent(in) :: unit

    
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (integerAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    integer, pointer :: value 
    

    allocate(value )
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (integerAttribute), intent(inout) :: this
    
  end subroutine clean

end module integerAttribute_mod
#undef TYPE


module integer1dAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  use integerAttribute_mod
  implicit none
  private

#define TYPE integer1dAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    integer, allocatable :: value (:)

  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    integer, intent(in) :: value (:)
    allocate(entry%value (size(value,1)))

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     integer , allocatable, intent(inout) :: value (:)
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)

#ifdef __GFORTRAN__
        value = entry%value
#else
        allocate(value, source=entry%value)
#endif



    type is (VectorAttribute)
       allocate(value(entry%size()))
       do i = 1, entry%size()
          value(i) = entry%items(i)%ptr
       end do

    class default
      call throwException('Illegal conversion of integer1dAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    integer, pointer :: ptr (:)
    class (AbstractAttribute), target, intent(in) :: entry
    integer :: cast (:)

    select type (q => entry)
    type is (integer1dAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of integer1dAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (integer1dAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (integer1dAttribute)
      if (all(this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (integer1dAttribute), intent(in) :: this
    print*,'  Type:  ', 'integer1dAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (integer1dAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = join(reshape(toStringElemental(this%value),(/size(this%value)/)),', ')

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (integer1dAttribute), intent(in) :: this
    integer, intent(in) :: unit

    write(unit) shape(this%value)
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (integer1dAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    integer, pointer :: value (:)
    

    integer :: attributeShape(1)

    read(unit) attributeShape
    allocate(value(attributeShape(1)))
    
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (integer1dAttribute), intent(inout) :: this
    deallocate(this%value)
  end subroutine clean

end module integer1dAttribute_mod
#undef TYPE


module logicalAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  
  implicit none
  private

#define TYPE logicalAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    logical :: value
  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    logical, intent(in) :: value 
    

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     logical , intent(inout) :: value 
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)
        value = entry%value



    class default
      call throwException('Illegal conversion of logicalAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    logical, pointer :: ptr 
    class (AbstractAttribute), target, intent(in) :: entry
    logical :: cast 

    select type (q => entry)
    type is (logicalAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of logicalAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (logicalAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (logicalAttribute)
      if ((this%value .eqv. p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (logicalAttribute), intent(in) :: this
    print*,'  Type:  ', 'logicalAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (logicalAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = toStringElemental(this%value)

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (logicalAttribute), intent(in) :: this
    integer, intent(in) :: unit

    
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (logicalAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    logical, pointer :: value 
    

    allocate(value )
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (logicalAttribute), intent(inout) :: this
    
  end subroutine clean

end module logicalAttribute_mod
#undef TYPE


module logical1dAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  use logicalAttribute_mod
  implicit none
  private

#define TYPE logical1dAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    logical, allocatable :: value (:)

  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    logical, intent(in) :: value (:)
    allocate(entry%value (size(value,1)))

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     logical , allocatable, intent(inout) :: value (:)
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)

#ifdef __GFORTRAN__
        value = entry%value
#else
        allocate(value, source=entry%value)
#endif



    type is (VectorAttribute)
       allocate(value(entry%size()))
       do i = 1, entry%size()
          value(i) = entry%items(i)%ptr
       end do

    class default
      call throwException('Illegal conversion of logical1dAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    logical, pointer :: ptr (:)
    class (AbstractAttribute), target, intent(in) :: entry
    logical :: cast (:)

    select type (q => entry)
    type is (logical1dAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of logical1dAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (logical1dAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (logical1dAttribute)
      if (all(this%value .eqv. p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (logical1dAttribute), intent(in) :: this
    print*,'  Type:  ', 'logical1dAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (logical1dAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = join(reshape(toStringElemental(this%value),(/size(this%value)/)),', ')

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (logical1dAttribute), intent(in) :: this
    integer, intent(in) :: unit

    write(unit) shape(this%value)
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (logical1dAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    logical, pointer :: value (:)
    

    integer :: attributeShape(1)

    read(unit) attributeShape
    allocate(value(attributeShape(1)))
    
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (logical1dAttribute), intent(inout) :: this
    deallocate(this%value)
  end subroutine clean

end module logical1dAttribute_mod
#undef TYPE


module RealDPAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  
  implicit none
  private

#define TYPE RealDPAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    real(kind=DP) :: value
  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    real(kind=DP), intent(in) :: value 
    

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     real(kind=DP) , intent(inout) :: value 
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)
        value = entry%value



    class default
      call throwException('Illegal conversion of RealDPAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    real(kind=DP), pointer :: ptr 
    class (AbstractAttribute), target, intent(in) :: entry
    real(kind=DP) :: cast 

    select type (q => entry)
    type is (RealDPAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of RealDPAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (RealDPAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (RealDPAttribute)
      if ((this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (RealDPAttribute), intent(in) :: this
    print*,'  Type:  ', 'RealDPAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (RealDPAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = toStringElemental(this%value)

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (RealDPAttribute), intent(in) :: this
    integer, intent(in) :: unit

    
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (RealDPAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    real(kind=DP), pointer :: value 
    

    allocate(value )
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (RealDPAttribute), intent(inout) :: this
    
  end subroutine clean

end module RealDPAttribute_mod
#undef TYPE


module RealDP1dAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  use RealDPAttribute_mod
  implicit none
  private

#define TYPE RealDP1dAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    real(kind=DP), allocatable :: value (:)

  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    real(kind=DP), intent(in) :: value (:)
    allocate(entry%value (size(value,1)))

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     real(kind=DP) , allocatable, intent(inout) :: value (:)
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)

#ifdef __GFORTRAN__
        value = entry%value
#else
        allocate(value, source=entry%value)
#endif



    type is (VectorAttribute)
       allocate(value(entry%size()))
       do i = 1, entry%size()
          value(i) = entry%items(i)%ptr
       end do

    class default
      call throwException('Illegal conversion of RealDP1dAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    real(kind=DP), pointer :: ptr (:)
    class (AbstractAttribute), target, intent(in) :: entry
    real(kind=DP) :: cast (:)

    select type (q => entry)
    type is (RealDP1dAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of RealDP1dAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (RealDP1dAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (RealDP1dAttribute)
      if (all(this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (RealDP1dAttribute), intent(in) :: this
    print*,'  Type:  ', 'RealDP1dAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (RealDP1dAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = join(reshape(toStringElemental(this%value),(/size(this%value)/)),', ')

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (RealDP1dAttribute), intent(in) :: this
    integer, intent(in) :: unit

    write(unit) shape(this%value)
    
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    
    class (RealDP1dAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    real(kind=DP), pointer :: value (:)
    

    integer :: attributeShape(1)

    read(unit) attributeShape
    allocate(value(attributeShape(1)))
    
    

    read(unit) value

    

    allocate(new, source=newAttribute(value))
    deallocate(value)
    

  end function readUnformatted

  subroutine clean(this)
    class (RealDP1dAttribute), intent(inout) :: this
    deallocate(this%value)
  end subroutine clean

end module RealDP1dAttribute_mod
#undef TYPE


module StringAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  
  implicit none
  private

#define TYPE StringAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    character(len=MAX_LEN_ATTRIBUTE_STRING) :: value
  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    character(len=*), intent(in) :: value 
    

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     character(len=MAX_LEN_ATTRIBUTE_STRING) , intent(inout) :: value 
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)
        value = entry%value



    class default
      call throwException('Illegal conversion of StringAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: ptr 
    class (AbstractAttribute), target, intent(in) :: entry
    character(len=MAX_LEN_ATTRIBUTE_STRING) :: cast 

    select type (q => entry)
    type is (StringAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of StringAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (StringAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (StringAttribute)
      if ((this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (StringAttribute), intent(in) :: this
    print*,'  Type:  ', 'StringAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (StringAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = toStringElemental(this%value)

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (StringAttribute), intent(in) :: this
    integer, intent(in) :: unit

    
    write(unit) len_trim(this%value)
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    use StringUtilities_mod, only: forceTrim
    class (StringAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: value 
    integer, pointer :: lengths 

    allocate(value , lengths)
    read(unit) lengths

    read(unit) value

    call forceTrim(value,lengths)

    allocate(new, source=newAttribute(value))
    deallocate(value)
    deallocate(lengths)

  end function readUnformatted

  subroutine clean(this)
    class (StringAttribute), intent(inout) :: this
    
  end subroutine clean

end module StringAttribute_mod
#undef TYPE


module String1dAttribute_mod
  use AbstractAttribute_mod
  use AttributeReference_mod
  use StringAttribute_mod
  implicit none
  private

#define TYPE String1dAttribute

  public :: TYPE
  public :: newAttribute
  public :: toType
  public :: toPointer
  public :: assignment(=)

  type, extends(AbstractAttribute) :: TYPE
    character(len=MAX_LEN_ATTRIBUTE_STRING), allocatable :: value (:)

  contains   
    procedure :: equals
    procedure :: clean
    procedure :: print => printIt ! gfortran workaround
    procedure :: writeUnformatted
    procedure :: readUnformatted
    procedure :: toString
  end type TYPE

  interface newAttribute
    module procedure constructor
  end interface

  interface assignment(=)
     module procedure toType_
  end interface

  interface toType
     module procedure toType_
  end interface

  interface toPointer
     module procedure toPointerType
  end interface toPointer

contains

  function constructor(value) result(entry)
    type (TYPE) :: entry
    character(len=*), intent(in) :: value (:)
    allocate(entry%value (size(value,1)))

    entry%value = value

  end function constructor

  subroutine toType_(value, entry)
     character(len=MAX_LEN_ATTRIBUTE_STRING) , allocatable, intent(inout) :: value (:)
    class (AbstractAttribute), intent(in) :: entry

    integer :: i

    select type (entry)
    type is (TYPE)

#ifdef __GFORTRAN__
        value = entry%value
#else
        allocate(value, source=entry%value)
#endif



    type is (VectorAttribute)
       allocate(value(entry%size()))
       do i = 1, entry%size()
          value(i) = entry%items(i)%ptr
       end do

    class default
      call throwException('Illegal conversion of String1dAttribute.',255)
    end select
  end subroutine toType_

  function toPointerType(entry, cast) result(ptr)
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: ptr (:)
    class (AbstractAttribute), target, intent(in) :: entry
    character(len=MAX_LEN_ATTRIBUTE_STRING) :: cast (:)

    select type (q => entry)
    type is (String1dAttribute)
       ptr => q%value
    class default
      call throwException('Illegal association of String1dAttribute.',255)
    end select
  end function toPointerType


  logical function equals(this, b)
    class (String1dAttribute), intent(in) :: this
    class (AbstractAttribute), intent(in) :: b

    select type (p => b)
    class is (String1dAttribute)
      if (all(this%value == p%value)) then
        equals = .true.
      else
        equals = .false.
      end if
    class default
      equals = .false.
    end select

  end function equals

  subroutine printIt(this)
    class (String1dAttribute), intent(in) :: this
    print*,'  Type:  ', 'String1dAttribute'
    print*,'  Value: <', this%value,'>'
    print*,'--------------'
  end subroutine printIt

  function toString(this) result(string)
    use StringUtilities_mod, only: toStringElemental => toString
    class (String1dAttribute), intent(in) :: this
    character(len=MAX_LEN_LINE) :: string

    string = join(reshape(toStringElemental(this%value),(/size(this%value)/)),', ')

  contains

    function join(strArray, separator) result(string)
      character(len=*), intent(in) :: strArray(:)
      character(len=*), intent(in) :: separator
      character(len=MAX_LEN_LINE) :: string

      integer :: i
      string = trim(strArray(1))
      do i = 2, size(strArray)
        string = trim(string) // trim(separator) // trim(strArray(i))
      end do
    end function join

  end function toString

  subroutine writeUnformatted(this, unit)
    class (String1dAttribute), intent(in) :: this
    integer, intent(in) :: unit

    write(unit) shape(this%value)
    write(unit) len_trim(this%value)
    write(unit) this%value

  end subroutine writeUnformatted

  function readUnformatted(this, unit) result(new)
    use StringUtilities_mod, only: forceTrim
    class (String1dAttribute), intent(in) :: this
    integer, intent(in) :: unit
    class (AbstractAttribute), pointer :: new

    integer :: rank
    character(len=MAX_LEN_ATTRIBUTE_STRING), pointer :: value (:)
    integer, pointer :: lengths (:)

    integer :: attributeShape(1)

    read(unit) attributeShape
    allocate(value(attributeShape(1)))
    allocate(lengths(attributeShape(1)))
    
    read(unit) lengths

    read(unit) value

    call forceTrim(value,lengths)

    allocate(new, source=newAttribute(value))
    deallocate(value)
    deallocate(lengths)

  end function readUnformatted

  subroutine clean(this)
    class (String1dAttribute), intent(inout) :: this
    deallocate(this%value)
  end subroutine clean

end module String1dAttribute_mod
#undef TYPE



module Attributes_mod
  use AbstractAttribute_mod
  use integerAttribute_mod
  use integer1dAttribute_mod
  use logicalAttribute_mod
  use logical1dAttribute_mod
  use RealDPAttribute_mod
  use RealDP1dAttribute_mod
  use StringAttribute_mod
  use String1dAttribute_mod
  implicit none

  public :: assignment(=)
  public :: toPointer
  public :: toType

end module Attributes_mod  
