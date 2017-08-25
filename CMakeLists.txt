cmake_minimum_required (VERSION 3.0)
project (hirshfeld)

include(GNUInstallDirs)

set(VERSION "0.2.5")
set(PACKAGE_NAME "hirshfeld")

option (DEBUG 
        "turn on debug to get more information" OFF) 

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )



add_subdirectory(src)
add_subdirectory(data)


include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set (CPACK_PACKAGE_VERSION_MAJOR "0")
set (CPACK_PACKAGE_VERSION_MINOR "2")
set (CPACK_PACKAGE_VERSION_PATCH "5")
include (CPack)
