# Get ARCH information
execute_process( 
   COMMAND uname -m 
   COMMAND tr -d '\n' 
   OUTPUT_VARIABLE ARCHITECTURE 
)

message(ARCH ${ARCHITECTURE})
message("FORTRAN " ${CMAKE_Fortran_COMPILER_ID})

#if (MPI MATCHES YES)
#   # Add CPP flags for MPI=YES
#   add_definitions(-DUSE_MPI)
#   find_package(MPI REQUIRED)
#   include_directories(${MPI_INCLUDE_PATH})
#   set (COMPILER_WRAPPER mpif90)
#endif()

## Hashmap templates are here:
#include_directories(${CMAKE_SOURCE_DIR}/model/include)
#include_directories(${CMAKE_SOURCE_DIR}/model/shared)
## rundeck_opts are here:
#include_directories(${PROJECT_BINARY_DIR}/model/include)

## Add flags for PFUNIT
#if(EXISTS $ENV{PFUNIT})
#   include_directories($ENV{PFUNIT}/mod)
#   include_directories($ENV{PFUNIT}/include)
#   link_directories($ENV{PFUNIT}/lib)
#   set(CPPFLAGS -DUSE_PFUNIT)
#   set(WITH_PFUNIT YES)
#else()
#   set(WITH_PFUNIT NO)
#endif()
#
#set(CPP gcc)
set(CFLAGS "-O2 -m64")

# Intel compiler flags
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

   if (CMAKE_SYSTEM_NAME MATCHES Linux)
      set (CPPFLAGS 
         "${CPPFLAGS} -DCOMPILER_Intel8 -DCONVERT_BIGENDIAN -DMACHINE_Linux"
      )
   else()
      set (CPPFLAGS 
         "${CPPFLAGS} -DCOMPILER_Intel8 -DCONVERT_BIGENDIAN -DMACHINE_MAC"
      )
   endif()

   set (FFLAGS_RELEASE 
      "-assume protect_parens -fp-model strict -warn nousage -assume realloc_lhs"
   )
   # Base flags
   if (WITH_PFUNIT MATCHES NO)
      set(FFLAGS 
         "${CPPFLAGS} ${FFLAGS_RELEASE} -fpp -g -O2 -ftz -convert big_endian"
      )
      set(F90FLAGS 
         "${CPPFLAGS} ${FFLAGS_RELEASE} -fpp -g -O2 -ftz -convert big_endian -free"
      )
   else()
      set(FFLAGS 
         "${CPPFLAGS} -fpp -g -O2 -ftz -convert big_endian"
      )
      set(F90FLAGS 
         "${CPPFLAGS} -fpp -g -O2 -ftz -convert big_endian -free"
      )
   endif()

   if (COMPILE_WITH_DEBUG MATCHES YES)
      set(FFLAGS 
         "${FFLAGS} -g -O0 -traceback"
      )
      set(F90FLAGS 
         "${F90FLAGS} -g -O0 -traceback"
      )
      set (LFLAGS  "-O2 -ftz")
   endif()

   if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(FFLAGS 
         "${FFLAGS} -CB -fpe0 -check uninit -ftrapuv -traceback"
      )
      set(F90FLAGS 
         "${F90FLAGS} -CB -fpe0 -check uninit -ftrapuv -traceback"
      )
      set(LFLAGS "${LFLAGS} -CB -fpe0 -check uninit -ftrapuv")
   endif()

   set(R8 "-r8")
   set(EXTENDED_SOURCE "-extend_source")


# ===================================== GNU compiler flags
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

   set (CPP ${CMAKE_C_COMPILER} -E)

   # Base compiler
   set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -fbacktrace -DCOMPILER_G95")
   if (CMAKE_SYSTEM_NAME MATCHES Linux)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMACHINE_Linux")
   else()
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DMACHINE_MAC")
   endif()

   # ----------------- Base Flags
   # -fno-range-check:
   #    Disable range checking on results of simplification
   #    of constant expressions during compilation. For example, GNU
   #    Fortran will give an error at compile time when simplifying a =
   #    1. / 0. With this option, no error will be given and a will be
   #    assigned the value +Infinity. If an expression evaluates to a value
   #    outside of the relevant range of [-HUGE():HUGE()], then the
   #    expression will be replaced by -Inf or +Inf as
   #    appropriate. Similarly, DATA i/Z'FFFFFFFF'/ will result in an
   #    integer overflow on most systems, but with -fno-range-check the
   #    value will “wrap around” and i will be initialized to -1 instead.
   #  -fconvert=conversion
   #     Specify the representation of data for unformatted files. Valid
   #     values for conversion are: ‘native’, the default; ‘swap’, swap
   #     between big- and little-endian; ‘big-endian’, use big-endian
   #     representation for unformatted files; ‘little-endian’, use
   #     little-endian representation for unformatted files.
   #     This option has an effect only when used in the main program. The
   #     CONVERT specifier and the GFORTRAN_CONVERT_UNIT environment
   #     variable override the default specified by -fconvert.

   #  -cpp
   #     Enable preprocessing. The preprocessor is automatically invoked if
   #     the file extension is .fpp, .FPP, .F, .FOR, .FTN, .F90, .F95, .F03
   #     or .F08. Use this option to manually enable preprocessing of any
   #     kind of Fortran file.

   #     The preprocessor is run in traditional mode. Any restrictions of
   #     the file-format, especially the limits on line length, apply for
   #     preprocessed output as well, so it might be advisable to use the
   #     -ffree-line-length-none or -ffixed-line-length-none options.
   set(CMAKE_Fortran_FLAGS 
      "${CMAKE_Fortran_FLAGS} -cpp -fconvert=big-endian -fno-range-check -ffree-line-length-none -DUSE_MPI -DMPITYPE_LOOKUP_HACK"
   )




   # https://gcc.gnu.org/onlinedocs/gfortran/Code-Gen-Options.html#Code-Gen-Options
   if ("${COMPILE_WITH_DEBUG}" STREQUAL "YES")
      # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O1 -fcheck=all")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O1 -fcheck=bounds -fcheck=array-temps -fcheck=do -fcheck=mem -fcheck=recursion")
      # ModelE crashes under: -fcheck=pointer
   else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
   endif()

   if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow")
      # ModelE crashes under: -finit-real=snan 
   endif()
  
   set(R8 "-fdefault-real-8 -fdefault-double-8")
   set(EXTENDED_SOURCE "-ffixed-line-length-132")

endif()

#if (MPI MATCHES YES)
#if (NOT "${PNETCDFHOME}" STREQUAL "")
#  link_directories(${PNETCDFHOME}/lib)
#  include_directories(${PNETCDFHOME}/include)
#endif()
#endif()
#
#if (NOT "${NETCDFHOME}" STREQUAL "")
#  link_directories(${NETCDFHOME}/lib)
#  include_directories(${NETCDFHOME}/include)
#endif()
#
## Will this be an ENV variable?
#if (defined $ENV{GLINT2INCLUDEDIR})
#  include_directories(${GLINT2INCLUDEDIR}/include)
#endif()
#if (defined $ENV{GLINT2LIBDIR})
#  link_directories(${GLINT2LIBDIR}/lib)
## What about 
## -L/opt/local/lib -lproj -lblitz -lCGAL -lmpfr -lgmp
#endif()
#
## TODO:
#
## FVCUBED OPTIONS
#
#set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")
#set(CMAKE_SKIP_RPATH ON)

