# ------------------------------------------------------------------------ #
# modelE cmake setup
# ------------------------------------------------------------------------ #
#
# Usage: cmake -DRUN=<modelE run> /path/to/my/repo [optional arguments]  
#
# Required argument:
#
#	-DRUN=<modelE run>
#       /path/to/my/repo
#
# Optional arguments (when possible a default is assigned):
#
# 	-DCMAKE_CXX_COMPILER=<c++ compiler> (default is system's compiler)
#	-DCMAKE_C_COMPILER=<c compiler> (default is system's)
#	-DCMAKE_Fortran_COMPILER=<fortran compiler> (default is system's)
#	-DCMAKE_BUILD_TYPE=Debug (or Release)
#	-DCOMPILE_WITH_TRAPS=NO (or YES)
#	-DCOMPILE_IC=NO (or YES to compile init_cond directory)
#	-DCOMPILE_DIAGS=NO (or YES to compile mk_diags directory)
#	-DCOMPILE_AUX=NO (or YES to compile aux directory)
#	-DMPI=YES (or NO)
#	-DRUNSRC=<modelE template> (default is none)
#       -DMODEL_DECKS_DIR=/path/to/my/decks (default is none)
#
# NOTE: RUN must be specified with full path.
#
# Examples:
#
#  Run cmake in a temporary directory: 
#
#    1) cmake -DRUN=/path/my/decks/ef40cc.R /path/to/my/repo
#
#    Uses specified e4f40cc.R rundeck using source code in /path/to/my/repo
#
#  It is common to create a "build" subdirectory in the source directory ,
#  cd to it and run cmake as follows:
#
#    2) cmake -DRUN=e4f0.R -DRUNSRC=E4F40.R ..
#
#    Creates e4f0.R rundeck in current directory using the specified template.
#
#  If you want to store or share your rundeck, specify where it should be 
#  created using the MODEL_DECKS_DIR option:
#
#    3) cmake -DMODEL_DECKS_DIR=/path/to/my/decks -DRUN=e4f0.R -DRUNSRC=E4F40.R ..
#
#
# ------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.2)
project (ModelE)
enable_language(Fortran)

# Where to look first for project cmake modules,
# before ${CMAKE_ROOT}/Modules/ is checked
set(CMAKE_MODULE_PATH
   ${CMAKE_SOURCE_DIR}/cmake
)

include(PreventInSourceBuild)

# Process command line options
include(ProcessCommandLineOptions)

#
# Set some CMake defaults
include(DefineCMakeDefaults)

# Set platform settings and check for required compiler versions
include(DefinePlatformDefaults)

# Set compiler flags
#include ("cmake/PISM_CMake_macros.cmake")
#include ("cmake/GLINT2_CMake_macros.cmake")
include (ModelE_CMake_macros)
modele_set_flags()

# Find required packages
#pism_find_prerequisites()
#glint2_find_prerequisites()
modele_find_prerequisites()

# Specify all the subdirectories to build
if (COMPILE_MODEL MATCHES YES)
    add_subdirectory(model)
endif()

# Not yet on the master (dev) branch...
# add_subdirectory (pyext)

# Do we want to build aux?
if(COMPILE_AUX MATCHES YES)
   add_subdirectory(aux)
endif()

# Do we want to build mk_diags?
if(COMPILE_DIAGS MATCHES YES)
   add_subdirectory(model/mk_diags)
endif()

# Do we want to build init_cond?
if(COMPILE_IC MATCHES YES)
   add_subdirectory(init_cond)
endif()

# If PFUNIT environment variable is set then WITH_PFUNIT is set to YES
if(WITH_PFUNIT MATCHES YES)
   add_subdirectory(tests)
endif()

message(STATUS "********************************************")
message(STATUS "********** PROJECT: ${PROJECT_NAME} **********")
message(STATUS "Architecture: ${ARCHITECTURE}")
message(STATUS "System:       ${CMAKE_SYSTEM_NAME}")
message(STATUS "MODELERC:     $ENV{MODELERC}")
message(STATUS "COMPILER:     ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "RUNSRC:       ${RUNSRC}")
message(STATUS "RUN:          ${RUN}")
message(STATUS "MPI:          ${MPI}")
message(STATUS "WITH_PFUNIT:  ${WITH_PFUNIT}")
if (COMPILE_WITH_DEBUG MATCHES YES)
message(STATUS "COMPILE_WITH_DEBUG: ${COMPILE_WITH_DEBUG}")
endif()
if (COMPILE_WITH_TRAPS MATCHES YES)
message(STATUS "COMPILE_WITH_TRAPS: ${COMPILE_WITH_TRAPS}")
endif()
message(STATUS "********************************************")

