# The name of our project is "HELLO". CMakeLists files in this project can
# refer to the root source directory of the project as ${HELLO_SOURCE_DIR} and
# to the root binary directory of the project as ${HELLO_BINARY_DIR}.
cmake_minimum_required (VERSION 3.2)
project (ModelE)
enable_language(Fortran)

list (APPEND CMAKE_MODULE_PATH "${ModelE_SOURCE_DIR}/cmake")
#include(setup_rpath)


# ------- Process Configure Options
if (NOT DEFINED USE_FEXCEPTION)
	# Use the fexception package in stop_model, rather than simply terminating the process.
	# Allows for better debugging from Python
	set(USE_FEXCEPTION NO)
endif()

if(${USE_FEXCEPTION})
	# This will require the fexception library
	add_definitions(-DUSE_FEXCEPTION)
endif()

# ---------------------------------


#include ("cmake/PISM_CMake_macros.cmake")
#include ("cmake/GLINT2_CMake_macros.cmake")
include ("cmake/ModelE_CMake_macros.cmake")
include ("cmake/setCompilerFlags.cmake")

# ==========================================================
set(NETCDF_F77 "YES")
set(NETCDF_CXX "YES")
#pism_find_prerequisites()
#glint2_find_prerequisites()
modele_find_prerequisites()

modele_set_flags()



# Recurse into the "Hello" and "Demo" subdirectories. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (model)

# Not yet on the master (dev) branch...
# add_subdirectory (pyext)
