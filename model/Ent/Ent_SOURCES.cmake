
# Pre-process m4F90 files
file(GLOB m4files "Ent/*.m4f")
foreach(file ${m4files})
   get_filename_component (name_without_extension ${file} NAME_WE)

   add_custom_command (
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/Ent

      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/Ent/${name_without_extension}.f
      COMMAND m4
      ARGS -s -I${CMAKE_CURRENT_SOURCE_DIR}/Ent ${file} > ${CMAKE_CURRENT_BINARY_DIR}/Ent/${name_without_extension}.f
      DEPENDS ${file}
   )
endforeach()


# Set Sources
set(Ent_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_prescribed_drv.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_prescribed_drv_geo.f
   ${CMAKE_CURRENT_BINARY_DIR}/Ent/ent_mod.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/cohorts.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/patches.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/entcells.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/physutil.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/allometryfn.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/reproduction.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/phenology.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/respauto_physio.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/disturbance.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/soilbgc.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_const.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_types.f 
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_prescr_veg.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_prescribed_updates.f
   ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_debug.f)


# Set sources based on rundeck options
if (PFT_MODEL MATCHES "ENT")
   if (FLUXNET MATCHES "YES")
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_pfts_ENT_FLUXNET.f)
   else()
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_pfts_ENT.f)
   endif()
else()
   if (FLUXNET MATCHES "YES")
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_pfts_FLUXNET.f)
   else()
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_pfts.f)
   endif()
endif()

# Set Ent_SOURCES based on rundeck options
if (PS_MODEL MATCHES "FBB")
   list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/FBBphotosynthesis.f)
   if (RAD_MODEL MATCHES "GORT")
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/canopyradiation.f ${CMAKE_CURRENT_SOURCE_DIR}/Ent/canopygort.f)
   else()
      list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/canopyspitters.f)
   endif()
   if (PFT_MODEL MATCHES "ENT")
      if (FLUXNET MATCHES "YES")
         list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/FBBpfts_ENT_FLUXNET.f)
      else()
         list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/FBBpfts_ENT.f)
      endif()
   else()
      if (FLUXNET MATCHES "YES")
         list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/FBBpfts_FLUXNET.f)
      else()
         list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/FBBpfts.f)
      endif()
   endif()
else()
   list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/biophysics.f)
endif()

if (MIXED_CANOPY_OPT MATCHES "ENT")
   list(APPEND Ent_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/Ent/ent_make_struct.f)
endif()

# TODO: This needs to get fixed up...
# Would like to write all preprocessor flags into a .h file
# for easy inspection.  Also, preprocessor flag should match
# the Ent option name.

# Set CPPFLAGS based on rundeck options
if (MIXED_CANOPY_OPT MATCHES "YES")
   set(CPPFLAGS "${CPPFLAGS} -DMIXED_CANOPY")
endif()

if (ENT_STANDALONE_DIAG MATCHES "YES")
   set(CPPFLAGS "${CPPFLAGS} -DENT_STANDALONE_DIAG")
endif()

if (SITE MATCHES "YES")
   set(CPPFLAGS "${CPPFLAGS} -DSOILCARB_SITE")
endif()

if (PFT_MODEL MATCHES "ENT")
    add_definitions(-DPFT_MODEL_ENT)
endif()

if (FLUXNET MATCHES "YES")
    add_definitions(-DSOILCARB_SITE)
endif()

if (PFT_MODEL MATCHES "YES")
    add_definitions(-DENT_STANDALONE_DIAG)
endif()
