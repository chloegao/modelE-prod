# modelE unit tests

set (SOURCE_DIR 
   ${CMAKE_SOURCE_DIR}/model
)

# Unfortunately COMPILER is defined differently in modelE (gfortran)
# than in PFUNIT (GNU). So, we rename here to inherit the correct PFUNIT
# configuration.
if (${COMPILER} MATCHES gfortran)
   set(COMPILER GNU)
endif()
# NAG note: PFUNIT must be compiled with -kind=byte
if (${COMPILER} MATCHES gfortran)
   set(F90FLAGS "${F90FLAGS} -kind=byte")
endif()

list(APPEND TESTDIRS shared)
list(APPEND TESTDIRS profiler)
if (${USE_TRACERS} MATCHES YES)
   list(APPEND TESTDIRS tracers)
endif()
if (MPI MATCHES YES)
   list(APPEND TESTDIRS MPI_Support)
   
endif()

foreach(dir ${TESTDIRS})
   list(APPEND TESTLIBS ${dir}_tests)
   list(APPEND MODELLIBS ${dir})
   list(APPEND INCDIRS ${CMAKE_SOURCE_DIR}/tests/${dir})
   add_subdirectory(${dir})
endforeach()

include_directories(${INCDIRS})

# TARGETS

set(TEST_EXE tests.x)
add_executable(${TEST_EXE} $ENV{PFUNIT}/include/driver.F90 testSuites.inc)
target_link_libraries(${TEST_EXE}
   ${TESTLIBS}
   ${MODELLIBS}
   netcdf
   pfunit
)

if (MPI MATCHES YES)
   target_link_libraries(${TEST_EXE} pnetcdf ${MPI_Fortran_LIBRARIES})
   add_custom_target(tests
      COMMAND mpirun -np 12 ${TEST_EXE}
      DEPENDS ${TEST_EXE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "MPI run ${TEST_EXE}"
   )
else()
   add_custom_target(tests
      COMMAND ${TEST_EXE}
      DEPENDS ${TEST_EXE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "run ${TEST_EXE}"
   )
endif()

