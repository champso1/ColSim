cmake_minimum_required(VERSION 3.20)
project(ColSimTest_Milestone3 LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(milestone3 milestone3.cpp)
target_compile_options(milestone3 PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(milestone3 PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(milestone3 PUBLIC ${ColSim_LIBRARIES})
