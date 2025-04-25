#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "ColSim::ColSim" for configuration ""
set_property(TARGET ColSim::ColSim APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(ColSim::ColSim PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libColSim.a"
  )

list(APPEND _cmake_import_check_targets ColSim::ColSim )
list(APPEND _cmake_import_check_files_for_ColSim::ColSim "${_IMPORT_PREFIX}/lib/libColSim.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
