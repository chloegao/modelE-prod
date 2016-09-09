# This file contains CMake macros used in the root CMakeLists.txt

macro(modele_find_prerequisites)
  if (${USE_FEXCEPTION})
    find_package(FException REQUIRED)
  endif()

  if (COMPILE_MODEL)
    find_package (MPI REQUIRED)

    if (${USE_PNETCDF})
        find_package(PNetCDF REQUIRED)
#        find_package(NetCDF_CXX REQUIRED)
    endif()

  endif()
  
  # Other required libraries
  find_package (NetCDF4_Fortran REQUIRED)

  # Use option values to set compiler and linker flags
  set (ModelE_EXTERNAL_LIBS "")

endmacro()

macro(modele_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories(${NETCDF4_FORTRAN_INCLUDE_DIR})

  if (COMPILE_MODEL)
    include_directories(${MPI_Fortran_INCLUDE_PATH})
    list (APPEND ModelE_EXTERNAL_LIBS ${MPI_Fortran_LIBRARIES})
  endif()

  list (APPEND ModelE_EXTERNAL_LIBS ${NETCDF4_FORTRAN_LIBRARY})

  if (${USE_FEXCEPTION})
    include_directories(${FEXCEPTION_INCLUDE_DIR})
    list(APPEND ModelE_EXTERNAL_LIBS ${FEXCEPTION_LIBRARY})
  endif()

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE TAO_DIR TAO_INCLUDE_DIRS NETCDF_PAR_H)

  if (USE_PNETCDF MATCHES YES)
    include_directories(${PNETCDF_INCLUDE_DIR})
#    include_directories(${NETCDF_CXX_INCLUDE_DIR})
    list (APPEND ModelE_EXTERNAL_LIBS ${PNETCDF_LIBRARY})# ${NETCDF_CXX_LIBRARY})
  endif()

endmacro()


macro(aux_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories(${NETCDF4_FORTRAN_INCLUDE_DIR})

  list (APPEND AUX_EXTERNAL_LIBS ${NETCDF4_FORTRAN_LIBRARY})

endmacro()

macro(mkdiags_set_dependencies)

  # Set include and library directories for *required* libraries.
  include_directories(${NETCDF4_FORTRAN_INCLUDE_DIR})

  list (APPEND MKDIAGS_EXTERNAL_LIBS ${NETCDF4_FORTRAN_LIBRARY})

  # Hide distracting CMake variables
  mark_as_advanced(file_cmd MPI_LIBRARY MPI_EXTRA_LIBRARY
    CMAKE_OSX_ARCHITECTURES CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT
    MAKE_EXECUTABLE TAO_DIR TAO_INCLUDE_DIRS NETCDF_PAR_H)

endmacro()

# ------------------------------------------
macro(modele_set_flags)
   # Get ARCH information
  execute_process( 
    COMMAND uname -m 
    COMMAND tr -d '\n' 
    OUTPUT_VARIABLE ARCHITECTURE)

  set(CFLAGS "-O2 -m64")
  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    add_definitions(-DMACHINE_Linux)
  else()
    add_definitions(-DMACHINE_MAC)
  endif()
  if (MPI MATCHES "YES")
    add_definitions(-DMPI_LOOKUP_HACK)
  endif()

# ===================================== Intel compiler flags
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")

    set(CPP /usr/local/other/SLES11/gcc/4.9.1/bin/gcc -E)

    add_definitions(-DCOMPILER_Intel8 -DCONVERT_BIGENDIAN)

    set (CMAKE_Fortran_FLAGS_RELEASE "${CPPFLAGS} -fpp -O2 -g -convert big_endian \
        -ftz -assume protect_parens -fp-model strict -warn nousage -assume realloc_lhs")
    set (CMAKE_Fortran_FLAGS_DEBUG   "${CPPFLAGS} -fpp -O0 -g -traceback \
        -ftz -convert big_endian \
        -assume protect_parens -fp-model strict -warn nousage -assume realloc_lhs")

    if (CMAKE_BUILD_TYPE MATCHES Release)
      set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
    else()
      set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
      set (LFLAGS  "-O2 -ftz")
    endif()

    if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CB -fpe0 \
           -check uninit -ftrapuv -traceback")
      set(LFLAGS "${LFLAGS} -CB -fpe0 -check uninit -ftrapuv")
    endif()

    set(R8 "-r8")
    set(EXTENDED_SOURCE "-extend_source")

# ===================================== GNU compiler flags
  elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL GNU)

    set (CPP ${CMAKE_C_COMPILER} -E)

    set(CPPFLAGS "${CPPFLAGS} -cpp")
    add_definitions(-DCOMPILER_G95)
    # This breaks if you try to wrap the long line in the obvious way.
    set (CMAKE_Fortran_FLAGS_RELEASE "${CPPFLAGS} -O2 -g -fconvert=big-endian -fno-range-check -ffree-line-length-none")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CPPFLAGS} -O -g -fbacktrace -fconvert=big-endian -fno-range-check -ffree-line-length-none -fcheck=bounds -fcheck=do -fcheck=mem -fcheck=recursion")


    if (CMAKE_BUILD_TYPE MATCHES Release)
      set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_RELEASE})
    else()
      set(CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS_DEBUG})
    endif()

    if ("${COMPILE_WITH_TRAPS}" STREQUAL "YES")
      set(CMAKE_Fortran_FLAGS
           "${CMAKE_Fortran_FLAGS} -fbounds-check -fcheck-array-temporaries -ffpe-trap=invalid,zero,overflow")
    endif()

    set(LFLAGS  "")
    set(R8 "-fdefault-real-8 -fdefault-double-8")
    set(EXTENDED_SOURCE "-ffixed-line-length-132")

  endif()

endmacro()
