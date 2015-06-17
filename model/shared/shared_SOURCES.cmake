# Remember that some .F90 files were preprocessed from .m4F90 files.
# TODO: We need to add in this code generation at some point.


# Create a post-processed AttributeHashMap.F90
add_custom_command(
   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/shared


   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/shared/AttributeHashMap.F90
   COMMAND ${CPP}
   ARGS -P -I${CMAKE_SOURCE_DIR}/model/include -I${CMAKE_SOURCE_DIR}/model/shared
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AttributeHashMap.F90 -o
   ${CMAKE_CURRENT_BINARY_DIR}/shared/AttributeHashMap.F90
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/shared/AttributeHashMap.F90


   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/shared/AbstractTimeStamp.F90
   COMMAND ${CPP}
   ARGS -P -I${CMAKE_SOURCE_DIR}/model/include -I${CMAKE_SOURCE_DIR}/model/shared
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AbstractTimeStamp.F90 -o
   ${CMAKE_CURRENT_BINARY_DIR}/shared/AbstractTimeStamp.F90
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/shared/AbstractTimeStamp.F90


   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/shared/CalendarDate.F90
   COMMAND ${CPP}
   ARGS -P -I${CMAKE_SOURCE_DIR}/model/include -I${CMAKE_SOURCE_DIR}/model/shared
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/CalendarDate.F90 -o
   ${CMAKE_CURRENT_BINARY_DIR}/shared/CalendarDate.F90
   DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/shared/CalendarDate.F90
)



set(shared_SOURCES
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/system_tools.c
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/libmodele_refaddr.c
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/orbpar.f
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AbstractAttribute.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/ArrayBundle_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AttributeDictionary.F90
   ${CMAKE_CURRENT_BINARY_DIR}/shared/AttributeHashMap.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AttributeReference.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Attributes.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Constants_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/CubicEquation_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Dictionary_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/FileManager_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/GaussianQuadrature.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/GenericType_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Geometry_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/GetTime_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/KeyValuePair_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Parser_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/MathematicalConstants.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/PlanetaryParams.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/PlanetParams_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/PolynomialInterpolator.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Precision_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Random_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/RootFinding_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/RunTimeControls_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/SpecialFunctions.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/stop_model.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/modele_python.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/StringUtilities_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/System.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/dast.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/SystemTimers_mod.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/TimeConstants.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Time.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Utilities.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/SystemTools.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/KindParameters.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Rational.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/BaseTime.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/TimeInterval.F90
   ${CMAKE_CURRENT_BINARY_DIR}/shared/AbstractTimeStamp.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AnniversaryDate.F90
   ${CMAKE_CURRENT_BINARY_DIR}/shared/CalendarDate.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/CalendarMonth.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AbstractCalendar.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/FixedCalendar.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/JulianCalendar.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/PlanetaryCalendar.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Time.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/OrbitUtilities.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/AbstractOrbit.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/FixedOrbit.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/Earth365DayOrbit.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/ParameterizedEarthOrbit.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/PlanetaryOrbit.F90
   ${CMAKE_CURRENT_SOURCE_DIR}/shared/ModelClock.F90
)
