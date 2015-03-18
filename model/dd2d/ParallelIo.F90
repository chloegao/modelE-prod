module ParallelIo_mod
   use dd2d_utils, only : dist_grid
   use pario
   implicit none
   private

   public :: ParallelIo ! type and constructor
   public :: doVar      ! procedure and method
   
   type ParallelIo
      private
      type (dist_grid), pointer :: grid
      integer :: fileId
      character(len=:), allocatable :: compareString
   contains
      generic :: doVar => doVar_2d_real64
      generic :: doVar => doVar_3d_real64
      generic :: doVar => doVar_4d_real64
      generic :: doVar => doVar_5d_real64
      generic :: doVar => doVar_2d_integer
      generic :: doVar => doVar_3d_integer
      generic :: doVar => doVar_4d_integer
      generic :: doVar => doVar_2d_logical
      procedure :: doVar_2d_real64
      procedure :: doVar_3d_real64
      procedure :: doVar_4d_real64
      procedure :: doVar_5d_real64
      procedure :: doVar_2d_integer
      procedure :: doVar_3d_integer
      procedure :: doVar_4d_integer
      procedure :: doVar_2d_logical
   end type ParallelIo

   ! This type is used to force keyword association
   ! for optional arguments in argument lists.
   ! it is intentionally private
   type UnusableArgument
   end type UnusableArgument


   interface doVar
      module procedure doVar_2d_real64
      module procedure doVar_3d_real64
      module procedure doVar_4d_real64
      module procedure doVar_5d_real64
      module procedure doVar_2d_integer
      module procedure doVar_3d_integer
      module procedure doVar_4d_integer
      module procedure doVar_2d_logical
   end interface doVar

   interface ParallelIo
      module procedure newParallelIo
   end interface ParallelIo


contains


   function newParallelIo(grid, fid, compareString) result(handle)
      type (ParallelIo) :: handle
      type (dist_grid), target :: grid
      integer, intent(in) :: fid
      character(len=*), optional, intent(in) :: compareString
      
      handle%grid => grid
      handle%fileId = fid
      if (present(compareString)) then
         handle%compareString = compareString
      end if

   end function newParallelIo

#define _DECLARE_(arr) real(kind=real64), dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:) :: arr
#define _PROC_NAME_ doVar_3d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:,:) :: arr
#define _PROC_NAME_ doVar_4d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:,:,:) :: arr
#define _PROC_NAME_ doVar_5d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:,:) :: arr
#define _PROC_NAME_ doVar_3d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:,:,:) :: arr
#define _PROC_NAME_ doVar_4d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) logical, dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_logical
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_


   ! Extract variable name from varInfo string
   ! Such strings are all of the form    <varName>[(dims)]
   function getVariableNameFrom(varInfo) result(variableName)
      character(len=:), allocatable :: variableName
      character(len=*), intent(in) :: varInfo

      integer :: idx
      
      idx = index(varInfo, '(')
      if (idx == 0) idx = len_trim(varInfo)
      variableName = trim(adjustl(varInfo(1:idx-1)))
      
   end function getVariableNameFrom


end module ParallelIo_mod

