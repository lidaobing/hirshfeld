cmake_minimum_required (VERSION 3.0)
project (hirshfeld
  VERSION "0.2.5")

include(GNUInstallDirs)

option (DEBUG 
        "turn on debug to get more information" OFF) 

configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
  )

SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_CXX_FLAGS_DEBUG "-O0")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(CMAKE_BUILD_TYPE MATCHES "")
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "CMAKE_BUILD_TYPE is " ${CMAKE_BUILD_TYPE})

add_subdirectory(src)
add_subdirectory(data)

include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
include (CPack)
