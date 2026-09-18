set(CMAKE_CXX_STANDARD 11)

# define build type
IF( DEFINED CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ELSE()
  SET( CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ENDIF()

# In CMake 4, 3.10 is deprecated and 3.5 has been removed.
# use 3.1...3.10 so all of these versions are acceptable as min. version.
# https://cmake.org/cmake/help/latest/command/cmake_minimum_required.html
cmake_minimum_required(VERSION 3.1...3.10)

PROJECT(openEMS CXX C)

# Comparing the engine output of two builds bit by bit is a way to check that a
# change which should not alter the result really does not. Two things get in
# the way of that, both enabled by default because both are wanted in a normal
# build. Turn them off in localConfig.cmake or on the cmake command line.

# Flushing denormals to zero drops the lowest bits of a decaying field.
OPTION(ENABLE_FLUSH_TO_ZERO "Flush denormal floating point values to zero in the engines" ON)

# Floating-point contraction fuses a multiply and an add into an FMA, which
# rounds once instead of twice, so the result depends on where the compiler
# chose to contract. Debug builds turn it off, they are not built for speed.
IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ffp-contract=off")
ENDIF()

# default
set(LIB_VERSION_MAJOR 0)
set(LIB_VERSION_MINOR 37)
set(LIB_VERSION_PATCH 0)
set(LIB_VERSION_STRING ${LIB_VERSION_MAJOR}.${LIB_VERSION_MINOR}.${LIB_VERSION_PATCH})

set(VERSION "v${LIB_VERSION_STRING}")

IF(EXISTS ${PROJECT_SOURCE_DIR}/localConfig.cmake)
   include(${PROJECT_SOURCE_DIR}/localConfig.cmake)
ENDIF()

#ADD_DEFINITIONS( -D__SSE2__ )

# Metal GPU backend for the GPU engine (--engine=gpu), macOS only
IF(APPLE)
  OPTION(ENABLE_METAL "Build the Metal GPU backend" ON)
ENDIF()
IF(APPLE AND ENABLE_METAL)
  message(STATUS "Metal GPU backend enabled")
  enable_language(OBJCXX)
  set(CMAKE_OBJCXX_STANDARD ${CMAKE_CXX_STANDARD})
  set(CMAKE_OBJCXX_FLAGS "${CMAKE_OBJCXX_FLAGS} -fobjc-arc")
  ADD_DEFINITIONS( -DOPENEMS_WITH_METAL )
ENDIF()

IF(NOT ENABLE_FLUSH_TO_ZERO)
  message(STATUS "Flush-to-zero disabled, the engines keep denormal values")
  ADD_DEFINITIONS( -DOPENEMS_NO_FLUSH_TO_ZERO )
ENDIF()

# add git revision
IF(EXISTS ${PROJECT_SOURCE_DIR}/.git )
  FIND_PACKAGE(Git)
  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --always
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GITREV
  )
  set(VERSION ${GITREV})
  string(STRIP ${VERSION} VERSION)
  message(STATUS "Found Git repository, ${PROJECT_NAME} version tag: ${VERSION}")
ENDIF()

ADD_DEFINITIONS(-DGIT_VERSION=\"${VERSION}\")

#
# Set up RPATH for the project
#
option(ENABLE_RPATH "Enable rpath support on Linux and Mac" ON)
if(NOT CMAKE_INSTALL_RPATH)
  # the RPATH to be used when installing, but only if it's not a system directory
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
  ENDIF("${isSystemDir}" STREQUAL "-1")
endif()
if(APPLE AND NOT CMAKE_INSTALL_NAME_DIR)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}")
endif()
if(UNIX AND ENABLE_RPATH)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# Set locations of extra CMake modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${openEMS_SOURCE_DIR}/cmake/Modules/")

# fparser
# $ cmake -D FPARSER_ROOT_DIR=~/opt/openEMS .
# SET(FPARSER_ROOT_DIR ~/opt/openEMS)
find_library(fparser_LIBRARIES 
  NAMES fparser
  HINTS ${FPARSER_ROOT_DIR}/lib${LIB_SUFFIX}
  NO_CMAKE_FIND_ROOT_PATH
)
message(STATUS "fparser: ${fparser_LIBRARIES}" )
#TODO test if fparser was found
INCLUDE_DIRECTORIES( ${FPARSER_ROOT_DIR}/include ) 

# CSXCAD
# $ cmake -D CSXCAD_ROOT_DIR=~/opt/openEMS .
# SET(CSXCAD_ROOT_DIR ~/opt/openEMS)
find_library(CSXCAD_LIBRARIES 
  NAMES CSXCAD
  HINTS ${CSXCAD_ROOT_DIR}/lib${LIB_SUFFIX}
  NO_CMAKE_FIND_ROOT_PATH
)
message(STATUS "CSXCAD_LIBRARIES: ${CSXCAD_LIBRARIES}" )
#TODO test if CSXCADs was found
find_path(CSXCAD_INCLUDE_DIR 
    NAMES ContinuousStructure.h 
    HINTS ${CSXCAD_ROOT_DIR}/include
    PATH_SUFFIXES "CSXCAD" ${CSXCAD_INCLUDE_DIR}
    NO_CMAKE_FIND_ROOT_PATH
)
message(STATUS "CSXCAD_INCLUDE_DIR: ${CSXCAD_INCLUDE_DIR}" )
INCLUDE_DIRECTORIES( ${CSXCAD_INCLUDE_DIR} ) 

# TinyXML module from https://github.com/ros/cmake_modules
find_package(TinyXML REQUIRED)
ADD_DEFINITIONS( -DTIXML_USE_STL )
message(STATUS "TinyXML_INCLUDE_DIR: ${TinyXML_INCLUDE_DIR}")
message(STATUS "TinyXML_LIBRARY: ${TinyXML_LIBRARY}")

INCLUDE_DIRECTORIES( ${TinyXML_INCLUDE_DIR} )

# hdf5
find_package(HDF5 COMPONENTS C HL REQUIRED)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIRS})
# hdf5 compat
#ADD_DEFINITIONS( -DH5_USE_16_API )
#ADD_DEFINITIONS( -DH5_BUILT_AS_DYNAMIC_LIB )

# boost
find_package(Boost 1.46 COMPONENTS
  thread
  date_time
  serialization
  chrono
  program_options
  REQUIRED
)
if (${Boost_MINOR_VERSION} LESS 69)
    # Since Boost 1.69, "system" is header only, so a check
    # is needed only in Boost 1.68 and earlier.
    find_package(Boost 1.46 COMPONENTS system thread date_time serialization chrono program_options REQUIRED)
endif()
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})

# vtk
find_package(VTK COMPONENTS vtkCommonCore NO_MODULE QUIET)
if ("${VTK_VERSION}" VERSION_GREATER "9")
    find_package(VTK REQUIRED COMPONENTS IOXML IOGeometry IOLegacy IOPLY NO_MODULE REQUIRED)
else()
    find_package(VTK REQUIRED COMPONENTS vtkIOXML vtkIOGeometry vtkIOLegacy vtkIOPLY NO_MODULE REQUIRED)
    include(${VTK_USE_FILE})
endif()

message(STATUS "Found package VTK. Using version " ${VTK_VERSION})
set( vtk_LIBS ${VTK_LIBRARIES} )
message(STATUS "vtk libraries " ${vtk_LIBS})

INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIR})

# independent tool
ADD_SUBDIRECTORY( nf2ff )

set(SOURCES
  openems.cpp
)

set(PUB_HEADERS openems.h openems_global.h)

# libs
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( FDTD )
ADD_SUBDIRECTORY( FDTD/extensions )
ADD_SUBDIRECTORY( Common )

INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )    #find tools

message(STATUS "Sources:"  ${SOURCES})

add_library( openEMS SHARED ${SOURCES})
#ADD_EXECUTABLE( openEMS main.cpp ${SOURCES})
set_target_properties(openEMS PROPERTIES VERSION ${LIB_VERSION_STRING} SOVERSION ${LIB_VERSION_MAJOR} )
TARGET_LINK_LIBRARIES( openEMS
  ${CSXCAD_LIBRARIES}
  ${fparser_LIBRARIES}
  ${TinyXML_LIBRARY}
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${Boost_LIBRARIES}
  ${vtk_LIBS}
)
IF(APPLE AND ENABLE_METAL)
  TARGET_LINK_LIBRARIES( openEMS "-framework Metal" "-framework Foundation" )
ENDIF()

if (WIN32)
    # make sure "BUILD_OPENEMS_LIB" is only set for openEMS (dll) not openEMS_bin (exe)
    target_compile_definitions(openEMS PRIVATE -DBUILD_OPENEMS_LIB )
endif (WIN32)

# main program
ADD_EXECUTABLE( openEMS_bin main.cpp )
SET_TARGET_PROPERTIES(openEMS_bin PROPERTIES OUTPUT_NAME openEMS)
TARGET_LINK_LIBRARIES(openEMS_bin openEMS)
 
if (WIN32)
    INSTALL(TARGETS openEMS RUNTIME DESTINATION bin ARCHIVE DESTINATION lib${LIB_SUFFIX})
else()
    INSTALL(TARGETS openEMS DESTINATION lib${LIB_SUFFIX})
endif()
INSTALL(TARGETS openEMS_bin DESTINATION bin)

if (UNIX)
    INSTALL( FILES openEMS.sh 
         DESTINATION bin 
         PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
endif ()
INSTALL(FILES ${PUB_HEADERS} DESTINATION include/openEMS)
INSTALL( DIRECTORY matlab DESTINATION share/openEMS )
INSTALL( DIRECTORY python/Tutorials DESTINATION share/openEMS/python )
INSTALL( DIRECTORY resources DESTINATION share/openEMS )
# TODO tarball, debug, release
