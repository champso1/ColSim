file(GLOB_RECURSE sources src/* inc/*)

add_library(ColSim_Plotting STATIC ${sources})
add_library(ColSim::Plotting ALIAS ColSim_Plotting)

target_include_directories(ColSim_Plotting PUBLIC inc)
target_link_libraries(ColSim_Plotting PUBLIC ColSim::Core) # we need core
