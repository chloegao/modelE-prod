# init_cond

include_directories(${PROJECT_BINARY_DIR}/model)

# Targets

add_executable(mkAIC.bin AIC.D771201.f HNTRPS.f)
target_link_libraries(mkAIC.bin modele)
add_executable(OIC.WOA98.exe OIC.WOA98.f)
add_executable(OPF.exe OPF.f)

install(TARGETS mkAIC.bin OIC.WOA98.exe OPF.exe DESTINATION bin)


