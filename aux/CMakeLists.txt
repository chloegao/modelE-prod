# aux

aux_set_dependencies()

include_directories(
   ${PROJECT_BINARY_DIR}/model
   ${PROJECT_BINARY_DIR}/model/include
   ${PROJECT_BINARY_DIR}/model/profiler
   ${PROJECT_BINARY_DIR}/model/shared
   ${PROJECT_BINARY_DIR}/model/giss_LSM
   ${PROJECT_BINARY_DIR}/model/Ent
   ${PROJECT_BINARY_DIR}/model/dd2d
)

# Targets

add_executable(RMS RMS.f)
add_executable(mkdeep.exe mkdeep.f)
add_executable(RDijk2llEM RDijk2llEM.f)
add_executable(inputll2cs inputLL2CS.f regrid.f ncio.f)
target_link_libraries(inputll2cs  ${AUX_EXTERNAL_LIBS})
add_executable(ll2cs ll2cs.f regrid.f ncio.f)
target_link_libraries(ll2cs  ${AUX_EXTERNAL_LIBS})
add_executable(ncll2cs ncll2cs.f regrid.f ncio.f)
target_link_libraries(ncll2cs  ${AUX_EXTERNAL_LIBS})

list(APPEND apps RMS)
list(APPEND apps mkdeep.exe) 

foreach(app ${apps})
   #message(STATUS " ---> Add executable ${app}")
   target_link_libraries(${app} modele shared profiler giss_LSM Ent dd2d
      ${AUX_EXTERNAL_LIBS})
endforeach()


