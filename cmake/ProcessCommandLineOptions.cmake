# USE_FEXCEPTION
if (NOT DEFINED USE_FEXCEPTION)
   # Use the fexception package in stop_model, rather than simply terminating
   # the process. Allows for better debugging from Python.
   set(USE_FEXCEPTION NO)
endif()

if(${USE_FEXCEPTION})
    # This will require the fexception library
    add_definitions(-DUSE_FEXCEPTION)
endif()

# BUILD_TYPE
if (NOT DEFINED CMAKE_BUILD_TYPE)
   # Use the fexception package in stop_model, rather than simply terminating
   # the process. Allows for better debugging from Python.
   message(STATUS "Setting build type to 'Debug' as none was specified.")
   set(CMAKE_BUILD_TYPE Debug)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "Traps")
endif()


# USE_MPI
if (NOT DEFINED MPI)
   set(MPI YES)
endif()
if(${MPI})
    add_definitions(-DUSE_MPI)
endif()

if (NOT DEFINED USE_PNETCDF)
    set(USE_PNETCDF YES)
endif()

if (NOT DEFINED COMPILE_MODEL)
    set(COMPILE_MODEL YES)
endif()

if (NOT DEFINED COMPILE_AUX)
    set(COMPILE_AUX NO)
endif()

if (NOT DEFINED COMPILE_DIAGS)
    set(COMPILE_DIAGS NO)
endif()

if (NOT DEFINED COMPILE_IC)
    set(COMPILE_IC NO)
endif()

