#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(test)
set(EXE_NAME "test") 

#----------------------------------------------------------------------------
# Setup GNUInstallDirs for installation
include(GNUInstallDirs)

#----------------------------------------------------------------------------
# The following lines are required so that the executable can dynamically locate 
# the library dependencies independent of use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_FULL_LIBDIR}" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_FULL_LIBDIR}")
ENDIF("${isSystemDir}" STREQUAL "-1")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

#----------------------------------------------------------------------------
# Force use of c++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------------
# Directories relevant for this project 
include_directories(
${PROJECT_SOURCE_DIR}/include 
${PROJECT_SOURCE_DIR}/src 
${CMAKE_CURRENT_BINARY_DIR}/include
)

#----------------------------------------------------------------------------
# Sources and headers for the project 
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# link against UMass-ZC
set(UMZC_LIB "UMass-ZC")  
include_directories($ENV{UMZCSYS}/include) 
link_directories($ENV{UMZCSYS}/lib)

#----------------------------------------------------------------------------
# Generate the executable 
add_executable(${EXE_NAME} ${SOURCES} ${PROJECT_SOURCE_DIR}/${EXE_NAME}.cc)
target_link_libraries(${EXE_NAME} ${UMZC_LIB})

#----------------------------------------------------------------------------
# Install items 
install (TARGETS ${EXE_NAME} DESTINATION ${INSTALL_DIR})
install (FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
# install(DIRECTORY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${CMAKE_INSTALL_PREFIX})
