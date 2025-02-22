file(GLOB_RECURSE sources src/* inc/*)

add_library(ColSim_Physics STATIC ${sources})
add_library(ColSim::Physics ALIAS ColSim_Physics)

target_include_directories(ColSim_Physics PUBLIC inc)
target_link_libraries(ColSim_Physics PUBLIC ColSim::Core ColSim::Plotting)
