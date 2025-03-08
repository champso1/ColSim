cmake_minimum_required(VERSION 3.20)
project(ColSimTest_Basic LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(basic basic.cpp)
target_compile_options(basic PUBLIC -O3)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(basic PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(basic PUBLIC ${ColSim_LIBRARIES})
