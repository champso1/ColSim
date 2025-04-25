cmake_minimum_required(VERSION 3.20)
project(ColSimTest_monte-carlo LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(monte-carlo monte-carlo.cpp)
target_compile_options(monte-carlo PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(monte-carlo PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(monte-carlo PUBLIC ${ColSim_LIBRARIES})
