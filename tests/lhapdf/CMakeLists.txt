add_executable(lhapdf lhapdf.cpp)
target_include_directories(lhapdf PUBLIC ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(lhapdf PUBLIC ${LHAPDF_LIBRARIES})
