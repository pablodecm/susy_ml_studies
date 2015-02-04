
# Create ROOT dictionary
root_generate_dictionary("GeneralSkimmer.h" "" "GeneralSkimmerDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )
root_generate_dictionary("EventData.h" "" "EventDataDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )
root_generate_dictionary("EventTopology.h" "" "EventTopologyDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )
root_generate_dictionary("EventRecObjects.h" "" "EventRecObjectsDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )
root_generate_dictionary("EventGenObjects.h" "" "EventGenObjectsDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )
root_generate_dictionary("EventHighLevel.h" "" "EventHighLevelDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )

# Specify .cpp files
set(SRC GeneralSkimmer.cpp mt2_bisect.cpp process_tchain.cpp GeneralSkimmerDict.cxx
  EventDataDict.cxx EventTopologyDict.cxx EventRecObjectsDict.cxx EventGenObjectsDict.cxx
  EventHighLevelDict.cxx)

# Set headers
set(HEADERS GeneralSkimmer.h mt2bisect.h)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Specify executable
add_executable(susy_ml_studies-bin process_tchain.cpp ${SRC} ${HEADERS})
target_link_libraries(susy_ml_studies-bin ${LIBRARIES})
set_target_properties(susy_ml_studies-bin PROPERTIES OUTPUT_NAME process_tchain)

# Build a library
add_library(susy_ml_studies ${SRC})
target_link_libraries(susy_ml_studies ${LIBRARIES})

# vim:ft=cmake

