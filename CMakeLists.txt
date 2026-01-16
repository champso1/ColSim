cmake_minimum_required(VERSION 3.20)
project(ColSim LANGUAGES CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/install)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(sources
  src/alphas.cpp
  src/colsim.cpp
  src/fourvector.cpp
  src/hard_process.cpp
  src/math.cpp
  src/parton_shower.cpp
  src/phase_space.cpp
  src/settings.cpp)
set(headers
  include/colsim/alphas.hpp
  include/colsim/colsim.hpp
  include/colsim/common.hpp
  include/colsim/event.hpp
  include/colsim/fourvector.hpp
  include/colsim/gnuplot.hpp
  include/colsim/hard_process.hpp
  include/colsim/math.hpp
  include/colsim/particle.hpp
  include/colsim/parton_shower.hpp
  include/colsim/phase_space.hpp
  include/colsim/settings.hpp
  include/colsim/utils.hpp)

add_library(colsim STATIC)
target_sources(colsim PRIVATE ${sources})
target_sources(colsim PUBLIC
  FILE_SET colsim_headers
  TYPE HEADERS
  BASE_DIRS ${CMAKE_SOURCE_DIR}/include
  FILES ${headers})

target_compile_options(colsim PUBLIC -Wall -Wextra)

target_include_directories(colsim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)


list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/share/cmake")
find_package(LHAPDF REQUIRED)
target_include_directories(colsim PUBLIC ${LHAPDF_INCLUDE_DIRS})
target_link_libraries(colsim PUBLIC ${LHAPDF_LIBRARIES})

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

# --- installing --- #


include(GNUInstallDirs)
install(TARGETS colsim
  EXPORT colsimTargets
  FILE_SET colsim_headers
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT colsimTargets
  FILE colsimTargets.cmake
  NAMESPACE colsim::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/colsim)

include(CMakePackageConfigHelpers)
set(INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/share/cmake/Config.cmake.in
  ${CMAKE_BINARY_DIR}/colsimConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/colsim
  PATH_VARS INSTALL_INCLUDE_DIR)
install(FILES
  ${CMAKE_BINARY_DIR}/colsimConfig.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/colsim)


add_subdirectory(examples)
