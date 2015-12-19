
# define build type
IF( DEFINED CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ELSE()
  SET( CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ENDIF()

PROJECT(nf2ff CXX)
cmake_minimum_required(VERSION 2.8)

set(SOURCES
  main.cpp
  nf2ff.cpp
  nf2ff_calc.cpp
  ../tools/array_ops.cpp
  ../tools/useful.cpp
  ../tools/hdf5_file_reader.cpp
  ../tools/hdf5_file_writer.cpp
)

#ADD_SUBDIRECTORY( ../tools )
set(HEADERS
  nf2ff.h
  nf2ff_calc.h
)

ADD_EXECUTABLE( nf2ff ${SOURCES} )

TARGET_LINK_LIBRARIES( nf2ff
  tinyxml
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
)

INSTALL(TARGETS nf2ff DESTINATION bin)


#TODO tarball, debug, release
