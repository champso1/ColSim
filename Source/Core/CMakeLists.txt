file(GLOB_RECURSE sources inc/* src/*)

add_library(ColSim_Core STATIC ${sources})
add_library(ColSim::Core ALIAS ColSim_Core)

target_include_directories(ColSim_Core PUBLIC inc)

# lhapdf
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/share/cmake")
find_package(LHAPDF REQUIRED)
target_include_directories(ColSim_Core PUBLIC ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(ColSim_Core PUBLIC ${LHAPDF_LIBRARIES})
