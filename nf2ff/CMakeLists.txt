
# define build type
IF( DEFINED CMAKE_BUILD_TYPE )
  SET( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ELSE()
  SET( CMAKE_BUILD_TYPE Release CACHE STRING "Set to either \"Release\" or \"Debug\"" )
ENDIF()

PROJECT(nf2ff CXX)

# In CMake 4, 3.10 is deprecated and 3.5 has been removed.
# use 3.0...3.10 so all of these versions are acceptable as min. version.
# https://cmake.org/cmake/help/latest/command/cmake_minimum_required.html
cmake_minimum_required(VERSION 3.0...3.10)

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

# Optional GPU backend of the surface integral, written in HIP so that it runs
# on AMD and on NVIDIA GPUs. Off by default, so a plain build has no GPU
# dependency; the CPU code is always built and stays the fallback (see
# nf2ff::SetBackend).
option(NF2FF_HIP "Build the HIP (AMD or NVIDIA GPU) backend of nf2ff" OFF)

if (NF2FF_HIP)
  if (CMAKE_VERSION VERSION_LESS 3.21)
    message(FATAL_ERROR "NF2FF_HIP needs CMake 3.21 or newer (3.28 for NVIDIA GPUs)")
  endif()
  # select the GPUs with CMAKE_HIP_ARCHITECTURES, the platform with CMAKE_HIP_PLATFORM (amd or nvidia)
  enable_language(HIP)
  find_package(hip REQUIRED)
  message(STATUS "nf2ff: HIP backend enabled (${CMAKE_HIP_COMPILER_ID} ${CMAKE_HIP_COMPILER_VERSION})")
  set(SOURCES ${SOURCES} nf2ff_gpu.hip)
  set(HEADERS ${HEADERS} nf2ff_gpu.h)
endif()

add_library( nf2ff SHARED ${SOURCES})
set_target_properties(nf2ff PROPERTIES VERSION ${LIB_VERSION_STRING} SOVERSION ${LIB_VERSION_MAJOR})
set_target_properties(nf2ff PROPERTIES CXX_STANDARD 11)
if (WIN32)
    target_compile_definitions(nf2ff PRIVATE -DBUILD_NF2FF_LIB )
endif (WIN32)

TARGET_LINK_LIBRARIES( nf2ff
  ${TinyXML_LIBRARY}
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
)

if (NF2FF_HIP)
  target_compile_definitions(nf2ff PRIVATE NF2FF_WITH_GPU)
  if (MSVC)
    # keep std::mutex usable with the older msvcp140.dll another install may have loaded first
    target_compile_definitions(nf2ff PRIVATE _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
  endif()
  target_link_libraries(nf2ff hip::host)
endif()

ADD_EXECUTABLE( nf2ff_bin main.cpp )
SET_TARGET_PROPERTIES(nf2ff_bin PROPERTIES OUTPUT_NAME nf2ff)
TARGET_LINK_LIBRARIES(nf2ff_bin nf2ff)

INSTALL(TARGETS nf2ff_bin DESTINATION bin)
if (WIN32)
    INSTALL(TARGETS nf2ff RUNTIME DESTINATION bin ARCHIVE DESTINATION lib${LIB_SUFFIX})
else()
    INSTALL(TARGETS nf2ff     DESTINATION lib${LIB_SUFFIX})
endif()

INSTALL(FILES nf2ff.h DESTINATION include/openEMS)


#TODO tarball, debug, release
