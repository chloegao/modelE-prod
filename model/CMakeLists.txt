# To turn of GCC errors in one file...
# http://stackoverflow.com/questions/59670/how-to-get-rid-of-deprecated-conversion-from-string-constant-to-char-warnin/541079#541079
#pragma GCC diagnostic ignored "-Wwrite-strings"


# Set Pism_EXTERNAL_LIBS and include directories.
#pism_set_dependencies()

message("XX PROJECT_SOURCE_DIR " ${PROJECT_SOURCE_DIR})
message("XX NETCDF_INCLUDES " ${NETCDF_INCLUDES})
message("XX PROJECT_BUILD_DIR " ${PROJECT_BUILD_DIR})
message("XX RUNDECK " ${RUNDECK})
modele_set_dependencies()

# -----------------------------------------------------------
# Get lists of source files for all ModelE components.
# They may or may not be used in the final ModelE library.
set(BUILD_CPP_OPTIONS )		# Accumulate options here

include(${CMAKE_CURRENT_SOURCE_DIR}/shared/shared_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/MPI_Support_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/solvers/solvers_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/giss_LSM/giss_LSM_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/dd2d/dd2d_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/Ent/Ent_SOURCES.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/profiler/profiler_SOURCES.cmake)
#include(${CMAKE_CURRENT_SOURCE_DIR}/landice/landice_SOURCES.cmake)

# -----------------------------------------------------------

execute_process(
	COMMAND
	python ${PROJECT_SOURCE_DIR}/cmake/scripts/rundeck_to_cmake.py
	${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${RUNDECK} ${BUILD_CPP_OPTIONS})

#list (APPEND CMAKE_MODULE_PATH " ${CMAKE_CURRENT_BINARY_DIR}")
message("XX CMAKE_MODULE_PATH " ${CMAKE_MODULE_PATH})

include(${CMAKE_CURRENT_BINARY_DIR}/modele_SOURCES.cmake)

# Test with just one source
#set(modele_SOURCES CONST.f)

add_library (modele SHARED ${modele_SOURCES})
target_link_libraries(modele ${ModelE_EXTERNAL_LIBS})

#add_library (modele ${modele_SOURCES})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/model/shared)


# Set RPATH in the installed library
# (although we are unlikely to ever install here)
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510
set_target_properties(modele PROPERTIES
	INSTALL_RPATH_USE_LINK_PATH TRUE)


install(TARGETS modele DESTINATION lib)



# ---------------------------------------------

# ---------------------------------------------
add_executable(modelexe main.F90)
target_link_libraries(modelexe modele ${ModelE_EXTERNAL_LIBS})

