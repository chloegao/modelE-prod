# modelE unit tests for shared component

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

# Preprocess Foo modules
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/FooHashMap.F90
   COMMAND ${CPP}
   ARGS -E -P -I${CMAKE_SOURCE_DIR}/model/include -I${CMAKE_SOURCE_DIR}/model/shared
   ${CMAKE_CURRENT_SOURCE_DIR}/FooHashMap.F90 >
   ${CMAKE_CURRENT_BINARY_DIR}/FooHashMap.F90
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/FooHashMap.F90
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/FooAssociativeArray.F90
   COMMAND ${CPP}
   ARGS -E -P -I${CMAKE_SOURCE_DIR}/model/include -I${CMAKE_SOURCE_DIR}/model/shared
   ${CMAKE_CURRENT_SOURCE_DIR}/FooAssociativeArray.F90 >
   ${CMAKE_CURRENT_BINARY_DIR}/FooAssociativeArray.F90
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/FooAssociativeArray.F90 
)

set (SOURCES
   ${CMAKE_CURRENT_BINARY_DIR}/FooAssociativeArray.F90 
   ${CMAKE_CURRENT_BINARY_DIR}/FooHashMap.F90
)

set (PFSOURCES
   Test_AssociativeArray.F90
   Test_AssociativeArrayIterator.F90
   Test_AttributeAssociativeArray.F90
   Test_AttributeDictionary.F90
   Test_BaseTime.F90
   Test_CalendarDateHashMap.F90
   Test_Dictionary.F90
   Test_GenericType.F90
   Test_Geometry.F90
   Test_HashMap.F90
   Test_HashMapIterator.F90
   Test_JulianCalendar.F90
   Test_KeyValuePair.F90
   Test_OrbitUtilities.F90
   Test_Parser.F90
   Test_Random.F90
   Test_Rational.F90
   Test_StringUtilities.F90
   Test_Time.F90
)

set(ALLSOURCES
   ${PFSOURCES}
   ${SOURCES}
)
set_source_files_properties(
   ${SOURCES}
   PROPERTIES COMPILE_FLAGS
   "${F90FLAGS}"
)

set_source_files_properties(
   ${PFSOURCES}
   PROPERTIES COMPILE_FLAGS
   "${F90FLAGS}"
   GENERATED TRUE
)

include_directories(
   ${CMAKE_BINARY_DIR}/model
   ${CMAKE_BINARY_DIR}/model/shared
   ${CMAKE_SOURCE_DIR}/model/include
   $ENV{PFUNIT}/include
   $ENV{PFUNIT}/mod
)

# Targets
add_library(${THIS} ${ALLSOURCES})
target_link_libraries(${THIS} shared)



