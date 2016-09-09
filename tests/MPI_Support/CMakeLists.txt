# modelE unit tests for MPI_Support component

get_filename_component(DIRNAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)
set(THIS ${DIRNAME}_tests)

# Pre-process pf files - output is in build directory
file(GLOB pfiles "*.pf")
foreach(file ${pfiles})
   get_filename_component(basename ${file} NAME_WE)
   add_custom_command (
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${basename}.F90
      COMMAND python
      ARGS $ENV{PFUNIT}/bin/pFUnitParser.py ${file} 
      ${CMAKE_CURRENT_BINARY_DIR}/${basename}.F90
      DEPENDS ${file}
   )
endforeach()

set (SOURCES
   Test_DistGrid.F90
   Test_GatherScatter.F90
   Test_GlobalSum.F90
   Test_Halo.F90
)

set_source_files_properties(
   ${SOURCES}
   PROPERTIES COMPILE_FLAGS
   "${F90FLAGS}"
   GENERATED TRUE
)

include_directories(
   ${CMAKE_BINARY_DIR}/model/shared
   ${CMAKE_BINARY_DIR}/model/MPI_Support
   $ENV{PFUNIT}/include
   $ENV{PFUNIT}/mod
)

# Targets

add_library(${THIS} ${SOURCES})
target_link_libraries(${THIS} shared MPI_Support)


