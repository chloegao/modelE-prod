# modelE unit tests for profiler component

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
   Test_ProfileReport.F90
   Test_ReportColumn.F90
   Test_TimeFormatUtilities.F90
   Test_Timer.F90
   Test_TimerList.F90
)

if (MPI MATCHES YES)
   set (SOURCES
      ${SOURCES}
      Test_TimerParallel.F90
   )
endif()

set_source_files_properties(
   ${SOURCES}
   PROPERTIES COMPILE_FLAGS
   "${F90FLAGS}"
   GENERATED TRUE
)

include_directories(
   ${CMAKE_BINARY_DIR}/model/shared
   ${CMAKE_BINARY_DIR}/model/profiler
   $ENV{PFUNIT}/include
   $ENV{PFUNIT}/mod
)

# Targets
add_library(${THIS} ${SOURCES})
target_link_libraries(${THIS} shared profiler)


