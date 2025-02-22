cmake_minimum_required(VERSION 3.20)
project(ColSim LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(Source)
add_subdirectory(Tests)



# copy the `res` directory to build folder
# so we can handle relative paths like "res/xx.xx" in the code
# once the executable is produced
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/res DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


# --------------------------
#   Doxygen/Documentation
# --------------------------

# set default option to no docs
option(BUILD_DOCS "Build documentation" OFF)

if(BUILD_DOCS)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
	set(DOXYGEN_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/docs)
	set(DOXYGEN_GENERATE_TREEVIEW YES)
	set(DOXYGEN_HAVE_DOT YES)
	file(GLOB_RECURSE manualfiles ${CMAKE_BINARY_SOURCE_DIR}/res/manual/*)
	doxygen_add_docs(doxygen
	  ${sources} ${manualfiles}
	  ALL)
  endif(DOXYGEN_FOUND)
endif(BUILD_DOCS)
