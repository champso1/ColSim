cmake_minimum_required(VERSION 3.20)
project(ColSimTest_settings LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(settings settings.cpp)
target_compile_options(settings PUBLIC -O2)

list(APPEND CMAKE_PREFIX_PATH "../../install/lib/cmake/")
find_package(ColSim REQUIRED)
target_include_directories(settings PUBLIC ${ColSim_INCLUDE_DIRS})
target_link_libraries(settings PUBLIC ${ColSim_LIBRARIES})
