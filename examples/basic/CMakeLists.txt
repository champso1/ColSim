add_executable(basic basic.cpp)
target_compile_options(basic PRIVATE -Wall -Wextra)

target_include_directories(basic PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(basic PRIVATE colsim)
