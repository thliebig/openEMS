
# define build type
IF( DEFINED CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ELSE()
  SET( CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ENDIF()

PROJECT(nf2ff CXX)
cmake_minimum_required(VERSION 2.8)

set(LIB_VERSION_MAJOR 0)
set(LIB_VERSION_MINOR 1)
set(LIB_VERSION_PATCH 0)
set(LIB_VERSION_STRING ${LIB_VERSION_MAJOR}.${LIB_VERSION_MINOR}.${LIB_VERSION_PATCH})

set(VERSION "v${LIB_VERSION_STRING}")

set(SOURCES
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

add_library( nf2ff SHARED ${SOURCES})
set_target_properties(nf2ff PROPERTIES VERSION ${LIB_VERSION_STRING} SOVERSION ${LIB_VERSION_MAJOR})
if (WIN32)
    target_compile_definitions(nf2ff PRIVATE -DBUILD_NF2FF_LIB )
endif (WIN32)

TARGET_LINK_LIBRARIES( nf2ff
  ${TinyXML_LIBRARY}
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
  ${MPI_LIBRARIES}
)

ADD_EXECUTABLE( nf2ff_bin main.cpp )
SET_TARGET_PROPERTIES(nf2ff_bin PROPERTIES OUTPUT_NAME nf2ff)
TARGET_LINK_LIBRARIES(nf2ff_bin nf2ff)

INSTALL(TARGETS nf2ff_bin DESTINATION bin)
if (WIN32)
    INSTALL(TARGETS nf2ff     DESTINATION bin)
else()
    INSTALL(TARGETS nf2ff     DESTINATION lib${LIB_SUFFIX})
endif()

INSTALL(FILES nf2ff.h DESTINATION include/openEMS)


#TODO tarball, debug, release
