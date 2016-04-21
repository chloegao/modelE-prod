set(dd2d_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/cdl_mod.f
   ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/timestream_mod.f
   ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/ParallelIo.F90
)

if (USE_PNETCDF)
   if(MPI)
      list(APPEND dd2d_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/pario_pnc.f)
   else()
      list(APPEND dd2d_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/pario_nc.f)
   endif()
else()
   list(APPEND dd2d_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/dd2d/pario_nc.f)
endif()
