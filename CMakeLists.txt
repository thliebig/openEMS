
# define build type
IF( DEFINED CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ELSE()
  SET( CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ENDIF()

PROJECT(openEMS CXX C)
cmake_minimum_required(VERSION 2.8)

# default
set(LIB_VERSION_MAJOR 0)
set(LIB_VERSION_MINOR 0)
set(LIB_VERSION_PATCH 35)
set(LIB_VERSION_STRING ${LIB_VERSION_MAJOR}.${LIB_VERSION_MINOR}.${LIB_VERSION_PATCH})

set(VERSION "v${LIB_VERSION_STRING}")

IF(EXISTS ${PROJECT_SOURCE_DIR}/localConfig.cmake)
   include(${PROJECT_SOURCE_DIR}/localConfig.cmake)
ENDIF()

#ADD_DEFINITIONS( -D__SSE2__ )

set(VERSION "v0.0.35")

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


if (WITH_MPI)
	ADD_DEFINITIONS(-DMPI_SUPPORT)
	# Require MPI for this project:
	find_package(MPI REQUIRED)
	INCLUDE_DIRECTORIES( ${MPI_INCLUDE_PATH} )

	message(STATUS "Found MPI. INCLUDE " ${MPI_INCLUDE_PATH})
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
find_package(HDF5 1.8 COMPONENTS C HL REQUIRED)
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARIES})
# hdf5 compat
#ADD_DEFINITIONS( -DH5_USE_16_API )
#ADD_DEFINITIONS( -DH5_BUILT_AS_DYNAMIC_LIB )

# boost
find_package(Boost 1.46 COMPONENTS
  thread
  system
  date_time
  serialization
  chrono
  REQUIRED
)
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES: ${Boost_LIBRARIES}")
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})

# vtk
find_package(VTK COMPONENTS vtkIOXML vtkIOGeometry vtkIOLegacy vtkIOPLY NO_MODULE REQUIRED)

message(STATUS "Found package VTK. Using version " ${VTK_VERSION})
if("${VTK_MAJOR_VERSION}" GREATER 5)
    set( vtk_LIBS ${VTK_LIBRARIES} )
else()
    set( vtk_LIBS
        vtkCommon
        vtkIO
    )    
endif()
message(STATUS "vtk libraries " ${vtk_LIBS})

include(${VTK_USE_FILE})
INCLUDE_DIRECTORIES (${VTK_INCLUDE_DIR})

if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
  set(ARCH "x86_64")
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "amd64")
  set(ARCH "x86_64")
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "AMD64")
  set(ARCH "x86_64")
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64")
  set(ARCH "ppc64")
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le")
  set(ARCH "ppc64")
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "unknown")
  set(ARCH "unknown")
  message(FATAL_ERROR "Unable to determine target architecture!  Try setting CMAKE_SYSTEM_PROCESSOR to a valid value.")
else()
  set(ARCH "unsupported")
endif()

if(${ARCH} STREQUAL "x86_64")
  message(STATUS "Detected 64-bit x86 target")
  #set(CMAKE_CXX_FLAGS "-msse -march=native")
elseif(${ARCH} STREQUAL "ppc64")
  message(STATUS "Detected 64-bit POWER target")
  set(CMAKE_CXX_FLAGS "-DNO_WARN_X86_INTRINSICS -DSSE_CORRECT_DENORMALS")
elseif(${ARCH} STREQUAL "unsupported")
  message(FATAL_ERROR "Unsupported target architecture!  Try porting openEMS to your architecture...")
else()
  message(FATAL_ERROR "Unable to determine target architecture!  Aborting.")
endif()

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

if (${MPI_CXX_FOUND})
	set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
	set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
	SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
endif()

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
  ${MPI_LIBRARIES}
)

if (WIN32)
    # make sure "BUILD_OPENEMS_LIB" is only set for openEMS (dll) not openEMS_bin (exe)
    target_compile_definitions(openEMS PRIVATE -DBUILD_OPENEMS_LIB )
endif (WIN32)

# main program
ADD_EXECUTABLE( openEMS_bin main.cpp )
SET_TARGET_PROPERTIES(openEMS_bin PROPERTIES OUTPUT_NAME openEMS)
TARGET_LINK_LIBRARIES(openEMS_bin openEMS)
 
if (WIN32)
    INSTALL(TARGETS openEMS DESTINATION bin)
else()
    INSTALL(TARGETS openEMS DESTINATION lib${LIB_SUFFIX})
endif()
INSTALL(TARGETS openEMS_bin DESTINATION bin)

if (UNIX)
    INSTALL( FILES openEMS.sh 
         DESTINATION bin 
         PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
	if (WITH_MPI)
	INSTALL( FILES openEMS_MPI.sh
		DESTINATION bin
		PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ WORLD_READ OWNER_EXECUTE GROUP_EXECUTE WORLD_EXECUTE)
	endif()
endif ()
INSTALL(FILES ${PUB_HEADERS} DESTINATION include/openEMS)
INSTALL( DIRECTORY matlab DESTINATION share/openEMS )
# TODO mpi, tarball, debug, release
