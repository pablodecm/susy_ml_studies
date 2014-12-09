
# Create ROOT dictionary
root_generate_dictionary("GeneralSkimmer.h" "" "GeneralSkimmerDict.cxx" ${CMAKE_CURRENT_SOURCE_DIR} )

# Specify .cpp files
set(SRC GeneralSkimmer.cpp one_ntuple.cpp GeneralSkimmerDict.cxx)

# Set headers
set(HEADERS GeneralSkimmer.h)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Specify executable
add_executable(susy_ml_studies-bin one_ntuple.cpp ${SRC} ${HEADERS})
target_link_libraries(susy_ml_studies-bin ${LIBRARIES})
set_target_properties(susy_ml_studies-bin PROPERTIES OUTPUT_NAME one_ntuple)

# Build a library
add_library(susy_ml_studies ${SRC})
target_link_libraries(susy_ml_studies ${LIBRARIES})

# vim:ft=cmake

