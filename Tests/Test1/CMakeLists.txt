add_executable(main main.cpp)
target_link_libraries(main PRIVATE ColSim)
target_compile_options(main PRIVATE -O3)
