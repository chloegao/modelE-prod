# init_cond

include_directories(${PROJECT_BINARY_DIR}/model)

# Targets

# This requires a Rundeck to build...
#add_executable(mkAIC.bin AIC.D771201.f HNTRPS.f)
#target_link_libraries(mkAIC.bin modele)

# Does not build; something wrong with FOCEAN/ZOCEAN common block
#add_executable(OIC.WOA98.exe OIC.WOA98.f)
add_executable(OPF.exe OPF.f)

install(TARGETS OPF.exe DESTINATION bin)
# mkAIC.bin

