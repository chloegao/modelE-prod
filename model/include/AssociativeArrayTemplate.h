#define IDENTITY(A) A
#define CONCAT(A,B) IDENTITY(A)IDENTITY(B)

#define DEFAULT_CONTAINER_TYPE(type) CONCAT(type,AssociativeArray)
#define DEFAULT_ITERATOR_TYPE(type) CONCAT(type,Iterator)
#define DEFAULT_MODULE_NAME(type) CONCAT(type,_mod)
#define DEFAULT_CONSTRUCTOR(type) CONCAT(new,type)

#ifndef USE_MODULE
#define USE_MODULE DEFAULT_MODULE_NAME(TYPE)
#endif

#ifndef CONTAINER_TYPE
#define CONTAINER_TYPE DEFAULT_CONTAINER_TYPE(TYPE)
#endif

#ifndef CONSTRUCTOR
#define CONSTRUCTOR DEFAULT_CONSTRUCTOR(CONTAINER_TYPE)
#endif

#ifndef MODULE_NAME
#define MODULE_NAME DEFAULT_MODULE_NAME(CONTAINER_TYPE)
#endif

#ifndef TYPE_NAME
#define TYPE_NAME TYPE
#endif

#ifndef ITERATOR_TYPE
#define ITERATOR_TYPE DEFAULT_ITERATOR_TYPE(CONTAINER_TYPE)
#endif

module MODULE_NAME
   use USE_MODULE, only: TYPE_NAME
   implicit none
   private

   public :: CONTAINER_TYPE ! The type for the Associative Array that is exported
   public :: CONSTRUCTOR
   public :: ITERATOR_TYPE
!!$   public :: assignment(=)
   public :: operator(/=)
   public :: operator(==)
   public :: clean
   public :: MAX_LEN_KEY

  integer, parameter :: MAX_LEN_KEY = 32
  integer, parameter :: INDEX_NOT_FOUND = -1

  type KeyValue
    character(len=MAX_LEN_KEY) :: key
    class (TYPE_NAME), allocatable :: value ! => null()
  end type KeyValue

  type :: CONTAINER_TYPE
    private
    integer :: numEntries
    type (KeyValue), allocatable :: entries(:)
  contains
    procedure :: size => getSize
    procedure :: getReference
    procedure :: insertEntry
    generic :: insert => insertEntry
    procedure :: insertReference
    procedure :: has => hasIt
    procedure :: merge
    procedure :: print
    ! iterator operations
    procedure :: begin
    procedure :: last
  end type CONTAINER_TYPE

  type :: ITERATOR_TYPE
!!$    private
    integer :: iter
    integer :: iterStop = 0
    class (CONTAINER_TYPE), pointer :: reference => null()
  contains
    procedure :: hasNext
    procedure :: next
    procedure :: key
    procedure :: getKey
    procedure :: value
  end type ITERATOR_TYPE

  interface clean
    module procedure clean_container
    module procedure clean_iterator
  end interface clean

!!$  interface assignment(=)
!!$    module procedure copy
!!$    ! TODO - thas copyIter procedure should not be necessary. Bug in intel?
!!$    module procedure copyIter
!!$  end interface assignment(=)

  interface operator(/=)
    module procedure notEqual
  end interface operator(/=)

  interface operator(==)
    module procedure equal
  end interface operator(==)

contains

  function CONSTRUCTOR () result(dictionary)
    type (CONTAINER_TYPE) :: dictionary

    allocate(dictionary%entries(0))
    dictionary%numEntries = 0
  end function CONSTRUCTOR

  function ITERATOR_CONSTRUCTOR (dictionary) result(iterator)
    type (CONTAINER_TYPE), target :: dictionary
    type (ITERATOR_TYPE), pointer :: iterator

    allocate(iterator)
    allocate(iterator%reference, SOURCE=dictionary)
  end function ITERATOR_CONSTRUCTOR

  integer function getSize(this) 
    class (CONTAINER_TYPE), intent(in) :: this
    getSize = this%numEntries
  end function getSize

  subroutine insertEntry(this, key, value)
    use StringUtilities_mod, only: toLowerCase
    class (CONTAINER_TYPE), target, intent(inout) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME) :: value

    type (KeyValue), allocatable :: tmpList(:)
    integer :: newCount
    integer :: i

    if (this%has(key)) then
      i = getIndex(this, key)
      deallocate(this%entries(i)%value)
      allocate(this%entries(i)%value, source=value)
      return
    end if
    
    if (this%numEntries > 0) then
       call move_alloc(from=this%entries, to=tmpList)
    else
       deallocate(this%entries)
    end if

    newCount = this%numEntries + 1

    allocate(this%entries(newCount))

    if (this%numEntries > 0) then
      do i = 1, this%numEntries
        this%entries(i)%key = tmpList(i)%key
#ifdef COMPILER_Intel8
        call kludge_move_alloc(tmplist(i)%value, this%entries(i)%value)
#else
        call move_alloc(from=tmplist(i)%value, to=this%entries(i)%value)
#endif
      end do
    end if

    if (this%numEntries > 0) then
       deallocate(tmpList)
    end if
    
    this%entries(newCount)%key = trim(toLowerCase(key))
    call kludge_allocate(this%entries(newCount)%value, source=value)
    this%numEntries = newCount

  contains

    subroutine kludge_move_alloc(from,to)
      class (TYPE_NAME), allocatable, intent(inout) :: from
      class (TYPE_NAME), allocatable, intent(out) :: to

      call move_alloc(from,to)
    end subroutine kludge_move_alloc

    subroutine kludge_allocate(element, source)
      class (TYPE_NAME), allocatable, intent(inout) :: element
      class (TYPE_NAME), intent(in) :: source

      allocate(element, source=source)
   end subroutine kludge_allocate

  end subroutine insertEntry

  subroutine insertReference(this, key, value)
    use StringUtilities_mod, only: toLowerCase
    class (CONTAINER_TYPE), intent(inout) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME), target :: value

!!$    type (KeyValue), allocatable :: tmpList(:)
!!$    integer :: newCount
!!$    integer :: i
!!$
!!$    if (this%has(key)) then
!!$      i = getIndex(this, key)
!!$      deallocate(this%entries(i)%value)
!!$      allocate(this%entries(i)%value, source=value)
!!$      return
!!$    end if
!!$    
!!$    call move_alloc(this%entries, tmpList)
!!$
!!$    newCount = this%numEntries + 1
!!$
!!$    allocate(this%entries(newCount))
!!$
!!$    if (this%numEntries > 0) then
!!$      do i = 1, this%numEntries
!!$        this%entries(i)%key = tmpList(i)%key
!!$        call move_alloc(tmplist(i)%value, this%entries(i)%value)
!!$      end do
!!$    end if
!!$
!!$    select case (trim(key))
!!$    case ('n2o5','N2O5','dms','DMS','SO4_d3','so4_d3')
!!$       write(77,*)'insert ', this%size()
!!$       do i = 1, this%size()
!!$          write(77,*)'      i=',i,trim(this%entries(i)%key)
!!$       end do
!!$    end select
!!$
!!$    deallocate(tmpList)
!!$    
!!$    this%entries(newCount)%key = trim(toLowerCase(key))
!!$    allocate(this%entries(newCount)%value, source=value)
!!$    this%numEntries = newCount
!!$
!!$    select case (trim(key))
!!$    case ('n2o5','N2O5')
!!$       write(77,*)'ha ha new: ','N2O5', this%size()
!!$       call this%print()
!!$    end select

  end subroutine insertReference

  function getReference(this, key) result(ptr)
    use StringUtilities_mod, only: toLowerCase
    class (CONTAINER_TYPE), target, intent(in) :: this
    character(len=*), intent(in) :: key
    class (TYPE_NAME), pointer :: ptr

    type (ITERATOR_TYPE) :: iter
    character(len=MAX_LEN_KEY) :: lowerCaseKey

    lowerCaseKey = trim(toLowerCase(key))

    iter = this%begin()

    do while (iter /= this%last())
      if (trim(iter%key()) == trim(lowerCaseKey)) then
        ptr => iter%value()
        return
      end if
      call iter%next()
    end do

    ptr => null()

  end function getReference

  logical function hasIt(this, key)
    use StringUtilities_mod, only: toLowerCase
    class (CONTAINER_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key
    
    integer :: i

    do i = 1, this%numEntries
      if (trim(toLowerCase(key)) == trim(this%entries(i)%key)) then
        hasIt = .true.
        return
      end if
    end do

    hasIt = .false.

  end function hasIt

  integer function getIndex(this, key) result(index)
    use StringUtilities_mod, only: toLowerCase
    class (CONTAINER_TYPE), intent(in) :: this
    character(len=*), intent(in) :: key
    
    integer :: i

    do i = 1, this%numEntries
      if (trim(toLowerCase(key)) == trim(this%entries(i)%key)) then
        index = i
        return
      end if
    end do

    index = INDEX_NOT_FOUND

  end function getIndex

  subroutine copy(a, b)
    type (CONTAINER_TYPE), intent(out) :: a
    type (CONTAINER_TYPE), intent(in)  :: b

    type (ITERATOR_TYPE) :: iter

    a%numEntries = 0
    allocate(a%entries(0))
    
    iter = b%begin()
    do while (iter /= b%last())
      call a%insert(iter%key(), iter%value())
      call iter%next()
    end do

  end subroutine copy

  subroutine copyIter(a, b)
     type (ITERATOR_TYPE), intent(inout) :: a
     type (ITERATOR_TYPE), intenT(in) :: b
     a%reference => b%reference
     a%iter = b%iter
     a%iterStop = b%iterStop
  end subroutine copyIter

  subroutine merge(this, b)
    class (CONTAINER_TYPE), intent(inout) :: this
    class (CONTAINER_TYPE), intent(in) :: b

    type (ITERATOR_TYPE) :: iter

    iter = b%begin()
    do while (iter /= b%last())
       if (.not. this%has(iter%key())) then
          call this%insert(iter%key(), iter%value())
       else
          call throwException('AssociativeArray::merge() failed due to duplicate keys: <' &
               & // trim(iter%key()) // '>.', 255)
       end if
      call iter%next()
    end do

  end subroutine merge

  subroutine print(this)
    class (CONTAINER_TYPE), intent(in) :: this
    type (ITERATOR_TYPE) :: iter
    class (TYPE_NAME), pointer :: t

    print*,'--------------------------'
    print*,' AssociativeArray: '
    print*,'--------------------------'

    iter = this%begin()
    do while (iter /= this%last())
      print*,'   key: <',trim(iter%key()),'>'
#ifdef HAS_PRINT
      t => iter%value()
      call t%print()
#endif
      call iter%next()
    end do
    print*,'--------------------------'
    print*,'--------------------------'
    print*,' '

  end subroutine print

  type (ITERATOR_TYPE) function begin(this) result(iterator)
    class (CONTAINER_TYPE), target, intent(in) :: this
    iterator%reference => this
    iterator%iter = 1
    iterator%iterStop = this%numEntries
  end function begin

  type (ITERATOR_TYPE) function last(this) result(iterator)
    class (CONTAINER_TYPE), target, intent(in) :: this
    iterator%reference => this
    iterator%iter = this%numEntries + 1
    iterator%iterStop = this%numEntries
  end function last

  logical function notEqual(a, b)
    class (ITERATOR_TYPE), intent(in) :: a
    class (ITERATOR_TYPE), intent(in) :: b
    ! TODO: throw exception if a and b do not have the same reference
    notEqual = .not. (a == b)
  end function notEqual

  logical function equal(a, b)
    class (ITERATOR_TYPE), intent(in) :: a
    class (ITERATOR_TYPE), intent(in) :: b
    ! TODO: throw exception if a and b do not have the same reference
    equal = (a%iter == b%iter)
  end function equal

  logical function hasNext(this)
    class (ITERATOR_TYPE), intent(in) :: this
    hasNext = (this%iter < this%iterStop)
  end function hasNext

  subroutine next(this)
    class (ITERATOR_TYPE), intent(inout) :: this
    this%iter = this%iter + 1
  end subroutine next

  function key(this)
    class (ITERATOR_TYPE), target, intent(in) :: this
    character(len=MAX_LEN_KEY), pointer :: key
    type (KeyValue), pointer :: p
    p => this%reference%entries(this%iter)
    key => p%key
!!$    key => this%reference%entries(this%iter)%key
  end function key

  subroutine getKey(this, key)
     class (ITERATOR_TYPE), target, intent(in) :: this
     character(len=MAX_LEN_KEY), pointer :: key
     type (KeyValue), pointer :: p
     p => this%reference%entries(this%iter)
     key => p%key
!!$    key => this%reference%entries(this%iter)%key
     
  end subroutine getKey

  function value(this)
    class (ITERATOR_TYPE), target, intent(in) :: this
    class (TYPE_NAME), pointer :: value
    type (KeyValue), pointer :: p
    p => this%reference%entries(this%iter)
    value => p%value
!!$    value => this%reference%entries(this%iter)%value
  end function value

  subroutine clean_container(this)
    type (CONTAINER_TYPE), intent(inout) :: this
    deallocate(this%entries)
  end subroutine clean_container

  subroutine clean_iterator(this)
    type (ITERATOR_TYPE), intent(inout) :: this
  end subroutine clean_iterator

end module MODULE_NAME

#undef USE_MODULE
#undef CONTAINER_TYPE
#undef CONSTRUCTOR
#undef MODULE_NAME
#undef TYPE_NAME
