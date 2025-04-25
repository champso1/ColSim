cmake_minimum_required(VERSION 3.20)
project(ColSimTest_gnuplot LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(gnuplot gnuplot.cpp)
target_compile_options(gnuplot PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(gnuplot PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(gnuplot PUBLIC ${ColSim_LIBRARIES})
