cmake_minimum_required(VERSION 3.20)
project(ColSimTest_four-vector LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(four-vector four-vector.cpp)
target_compile_options(four-vector PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(four-vector PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(four-vector PUBLIC ${ColSim_LIBRARIES})
