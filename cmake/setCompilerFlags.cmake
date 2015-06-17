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


# GNU compiler flags
elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

   set (CPP ${CMAKE_C_COMPILER} -E)

   if (CMAKE_SYSTEM_NAME MATCHES Linux)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -DCOMPILER_G95 -DMACHINE_Linux")
   else()
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -DCOMPILER_G95 -DMACHINE_MAC")
   endif()

   # Base flags
   # https://gcc.gnu.org/onlinedocs/gfortran/Code-Gen-Options.html
   set(CMAKE_Fortran_FLAGS 
      "${CMAKE_Fortran_FLAGS} -g -cpp -fconvert=big-endian -O2 -fno-range-check"
   )

   if ("${COMPILE_WITH_DEBUG}" STREQUAL "YES")
      set(CMAKE_Fortran_FLAGS 
      "${CMAKE_Fortran_FLAGS} -g -cpp -fconvert=big-endian -O1 -fcheck=mem -fcheck=pointer -fcheck=bounds -fno-range-check"
      )
   endif()

   if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(CMAKE_Fortran_FLAGS 
         "${CMAKE_Fortran_FLAGS} -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace"
      )
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

