add_executable(basic basic.cpp)
target_link_libraries(basic PRIVATE colsim)
target_compile_options(basic PRIVATE -O3)
