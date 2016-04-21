# main model CMakeLists.txt

get_filename_component(DIRNAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
set(LIB ${DIRNAME})

# In modelE, the final model executable is built upon specifications in
# a so-called 'rundeck' file. The rundeck contains information about all
# the components, model parameters and data files needed to compile and
# run the model.
# As a result one has to somehow parse the template file to extract the names
# of the components and the corresponding files, the compile and run options
# and the runtime parameters and data files. 

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
# This is a file included by model files.
set(RUNDECK_OPTS ${CMAKE_CURRENT_BINARY_DIR}/include/rundeck_opts.h)

# Export rundeck file settings into CMake
execute_process(
   COMMAND
   python ${PROJECT_SOURCE_DIR}/cmake/scripts/rundeck_to_cmake.py
   ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} ${RUN} ${BUILD_CPP_OPTIONS}
)

# include containers
include_directories(${PROJECT_SOURCE_DIR}/model/shared)

# include generated rundeck_opts
#include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include(${CMAKE_CURRENT_BINARY_DIR}/modele_SOURCES.cmake)

# Specify TARGETS. 
add_library(modele SHARED ${modele_SOURCES})
target_link_libraries(modele   ${ModelE_EXTERNAL_LIBS}  )

# Set RPATH in the installed library
# (although we are unlikely to ever install here)
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510

set_target_properties(modele PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS modele DESTINATION lib)

# modelE executable
add_executable(modelexe main.F90)

# Set additional link library dependencies
target_link_libraries(modelexe modele ${ModelE_EXTERNAL_LIBS})

