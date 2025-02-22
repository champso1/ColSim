file(GLOB_RECURSE sources src/* inc/*)

add_library(ColSim_Main STATIC ${sources})
add_library(ColSim ALIAS ColSim_Main)

target_include_directories(ColSim_Main PUBLIC inc)
target_link_libraries(ColSim_Main PUBLIC ColSim::Core ColSim::Physics)
