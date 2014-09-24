
PROJECT(openEMS CXX)
cmake_minimum_required(VERSION 2.8)

IF(EXISTS ${PROJECT_SOURCE_DIR}/localConfig.cmake)
   include(${PROJECT_SOURCE_DIR}/localConfig.cmake)
ENDIF()

# add git revision
IF(EXISTS ${PROJECT_SOURCE_DIR}/.git )
  FIND_PACKAGE(Git)
  # Get the latest abbreviated commit hash of the working branch
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
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
  LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
  IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
  ENDIF("${isSystemDir}" STREQUAL "-1")
endif()
if(APPLE AND NOT CMAKE_INSTALL_NAME_DIR)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
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
find_library(fparser_LIBRARIES NAMES fparser
  HINTS
  ${FPARSER_ROOT_DIR}/lib
)
message(STATUS "fparser: ${fparser_LIBRARIES}" )
#TODO test if fparser was found
INCLUDE_DIRECTORIES( ${FPARSER_ROOT_DIR}/include ) 

# CSXCAD
# $ cmake -D CSXCAD_ROOT_DIR=~/opt/openEMS .
# SET(CSXCAD_ROOT_DIR ~/opt/openEMS)
find_library(CSXCAD_LIBRARIES NAMES CSXCAD
  HINTS
  ${CSXCAD_ROOT_DIR}/lib
)
message(STATUS "CSXCAD: ${CSXCAD_LIBRARIES}" )
#TODO test if CSXCADs was found
INCLUDE_DIRECTORIES( ${CSXCAD_ROOT_DIR}/include/CSXCAD ) 

# TinyXML module from https://github.com/ros/cmake_modules
find_package(TinyXML REQUIRED)
ADD_DEFINITIONS( -DTIXML_USE_STL )

find_package(HDF5 REQUIRED)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})

# hdf5 compat
ADD_DEFINITIONS( -DH5_USE_16_API )

find_package(Boost 1.46 COMPONENTS
  thread
  system
  date_time
  serialization
  chrono
)

find_package(VTK REQUIRED)

if("${VTK_MAJOR_VERSION}" EQUAL 5)
  set( vtk_LIBS
    vtkCommon
    vtkFiltering
    vtkGraphics
    vtkHybrid
    vtkIO
    vtkRendering
    vtkWidgets
    QVTK
  )
endif()


if("${VTK_MAJOR_VERSION}" EQUAL 6)
  set( vtk_LIBS
    vtkCommonCore
    vtkCommonDataModel
    vtkIOLegacy
    vtkIOXML
    vtkIOGeometry
    vtkIOPLY
    vtksys
    vtkIOCore
)
endif()

if(VTK_FOUND)
	message(STATUS "Found package VTK. Using " ${VTK_USE_FILE})
    include(${VTK_USE_FILE})
else(VTK_FOUND)
    message(FATAL_ERROR "VTK not found. Please set VTK_DIR.")
endif(VTK_FOUND)

INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIR})

#message(status "vtk libs: ${vtk_LIBS}")

# independent tool
ADD_SUBDIRECTORY( nf2ff )

# libs
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( FDTD )
ADD_SUBDIRECTORY( FDTD/extensions )

# depend on CSXCAD
ADD_SUBDIRECTORY( Common )

# main program
set(SOURCES
  main.cpp
  openems.cpp
)

INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )    #find tools

ADD_EXECUTABLE( openEMS ${SOURCES} )

TARGET_LINK_LIBRARIES( openEMS
  FDTD
  Common
  extensions
  tools
  ${CSXCAD_LIBRARIES}
  ${fparser_LIBRARIES}
  tinyxml
  hdf5
  hdf5_hl
  #${HDF_LIBRARIES}
  ${Boost_LIBRARIES}
  #${VTK_LIBRARIES}
  ${vtk_LIBS}
)

INSTALL( TARGETS openEMS DESTINATION bin )

INSTALL( DIRECTORY matlab DESTINATION share/openEMS )

# TODO mpi, tarball, debug, release
