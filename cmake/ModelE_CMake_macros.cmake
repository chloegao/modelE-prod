# This file contains CMake macros used in the root CMakeLists.txt


macro(modele_find_prerequisites)
  if (${USE_FEXCEPTION})
    find_package(FException REQUIRED)
  endif()

  find_package (MPI REQUIRED)

  # Other required libraries
message(${CMAKE_MODULE_PATH})
  find_package (NetCDF REQUIRED)

  # Use option values to set compiler and linker flags
  set (ModelE_EXTERNAL_LIBS "")

endmacro()

macro(modele_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories (
    ${NETCDF_INCLUDES}
    ${MPI_Fortran_INCLUDE_PATH})

  list (APPEND ModelE_EXTERNAL_LIBS
    ${NETCDF_LIBRARIES}
    ${MPI_Fortran_LIBRARIES})

  if (${USE_FEXCEPTION})
    include_directories(${FEXCEPTION_INCLUDE_DIR})
    list(APPEND ModelE_EXTERNAL_LIBS ${FEXCEPTION_LIBRARY})
  endif()

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE TAO_DIR TAO_INCLUDE_DIRS NETCDF_PAR_H)

endmacro()
# ------------------------------------------
macro(modele_set_flags)

   # This need to be set from the rundeck and passed through CMake config...
   add_definitions(-DNEW_IO)

   if (CMAKE_SYSTEM_NAME MATCHES Linux)
      set (CPPFLAGS "${CPPFLAGS} -DCOMPILER_G95 -DMACHINE_Linux")
   else()
      set (CPPFLAGS "${CPPFLAGS} -DCOMPILER_G95 -DMACHINE_MAC")
   endif()

   # Base flags
   set(FFLAGS 
      "${CPPFLAGS} -cpp -fconvert=big-endian -O2 -fno-range-check"
   )
   set(F90FLAGS 
      "${CPPFLAGS} -cpp -fconvert=big-endian -O2 -fno-range-check -ffree-line-length-none"
   )
   if (COMPILE_WITH_DEBUG MATCHES YES)
      set(FFLAGS 
         "${FFLAGS} -g -O0 -fbacktrace"
      )
      set(F90FLAGS 
         "${F90FLAGS} -g -O0 -fbacktrace"
      )
   endif()

   if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(FFLAGS 
         "${FFLAGS} -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace"
      )
      set(F90FLAGS 
         "${F90FLAGS} -fbounds-check -fcheck-array-temporaries  -ffpe-trap=invalid,zero,overflow -finit-real=snan -fbacktrace"
      )
   endif()
  
   set(LFLAGS  "")
   set(R8 "-fdefault-real-8 -fdefault-double-8")
   set(EXTENDED_SOURCE "-ffixed-line-length-132")
endmacro()