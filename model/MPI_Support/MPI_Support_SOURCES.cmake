# Set sources
set(MPI_Support_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/dd2d_utils.f
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/DomainDecompLatLon.f
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/pario_fbsa.f
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/assert.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/dist_grid_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/DomainDecomposition_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/Domain_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/GatherScatter_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/GlobalSum_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/Halo_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/Hidden_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/MpiSupport_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/ProcessTopology_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/MPI_Support/SpecialIO_mod.F90
)
