set(CMAKE_CXX_STANDARD 17)

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

# zlib, for the compressed VTK files
find_package(ZLIB REQUIRED)
INCLUDE_DIRECTORIES (${ZLIB_INCLUDE_DIRS})

# independent tool
ADD_SUBDIRECTORY( nf2ff )

set(SOURCES
  openems.cpp
)

set(PUB_HEADERS openems.h openems_global.h)

# HIP GPU backend for the GPU engine (--engine=gpu), built if a HIP compiler is found.
# HIP runs on AMD GPUs through ROCm and on NVIDIA GPUs through CUDA, from the same
# sources: on the NVIDIA platform the compiler is nvcc and the binary needs the CUDA
# runtime, not a ROCm one.
#
# With OPENEMS_HIP_AS_CUDA the same sources are compiled as CUDA instead, for NVIDIA GPUs
# without a HIP compiler: it is the same compiler (nvcc) over the same headers, since the
# HIP headers of the NVIDIA platform map the calls to CUDA inline. This needs a CUDA
# compiler and the HIP headers (HIP_INCLUDE_DIR, looked for in ROCM_PATH and /opt/rocm;
# the hip-dev package of AMD holds the headers of every system).
# It is the default on Windows, where the HIP language of CMake is not usable for NVIDIA:
# it refuses every generator but Make and Ninja, and its compiler test fails to link with
# nvcc. For AMD it works there, with ROCm's clang and Ninja.
#
# Whatever is missing for the backend switches it off, with a warning if it was asked
# for: openEMS builds without it.
IF(POLICY CMP0077)
  # the options below keep a value set in localConfig.cmake
  cmake_policy(SET CMP0077 NEW)
ENDIF()
include(CheckLanguage)
set(HIP_MISSING "")
IF(WIN32)
  OPTION(OPENEMS_HIP_AS_CUDA "Compile the HIP GPU backend as CUDA sources (NVIDIA GPUs)" ON)
ELSE()
  OPTION(OPENEMS_HIP_AS_CUDA "Compile the HIP GPU backend as CUDA sources (NVIDIA GPUs)" OFF)
ENDIF()
IF(OPENEMS_HIP_AS_CUDA)
  IF(CMAKE_VERSION VERSION_LESS 3.24)
    set(HIP_MISSING "CMake 3.24 or newer")
  ELSE()
    check_language(CUDA)
    find_path(HIP_INCLUDE_DIR hip/hip_runtime.h HINTS ${ROCM_PATH}/include $ENV{ROCM_PATH}/include /opt/rocm/include)
    IF(NOT CMAKE_CUDA_COMPILER)
      set(HIP_MISSING "a CUDA compiler")
    ELSEIF(NOT HIP_INCLUDE_DIR)
      set(HIP_MISSING "the HIP headers (HIP_INCLUDE_DIR)")
    ENDIF()
  ENDIF()
ELSE()
  # CMAKE_HIP_PLATFORM, which tells AMD from NVIDIA below, came with CMake 3.28
  IF(CMAKE_VERSION VERSION_LESS 3.28)
    set(HIP_MISSING "CMake 3.28 or newer")
  ELSE()
    check_language(HIP)
    IF(NOT CMAKE_HIP_COMPILER)
      set(HIP_MISSING "a HIP compiler")
    ENDIF()
  ENDIF()
ENDIF()
IF(HIP_MISSING)
  OPTION(ENABLE_HIP "Build the HIP GPU backend" OFF)
ELSE()
  OPTION(ENABLE_HIP "Build the HIP GPU backend" ON)
ENDIF()
IF(ENABLE_HIP AND HIP_MISSING)
  message(WARNING "The HIP GPU backend needs ${HIP_MISSING}, building without it")
  set(ENABLE_HIP OFF)
ENDIF()
# The GPU architectures to build for, unless CMAKE_HIP_ARCHITECTURES /
# CMAKE_CUDA_ARCHITECTURES are given: the GPU of the build machine if it has one,
# else a list of all GPUs the backend is meant for, so that a build machine without
# a GPU (CI, packages) builds for the machines of the users. A GPU missing from the
# build is reported when the GPU engine starts; the CPU engines are unaffected.
set(OPENEMS_HIP_AMD_ARCHITECTURES
  "gfx908;gfx90a;gfx942;gfx950;gfx1010;gfx1011;gfx1012;gfx1030;gfx1031;gfx1032;gfx1033;gfx1034;gfx1035;gfx1036;gfx1100;gfx1101;gfx1102;gfx1103;gfx1150;gfx1151;gfx1152;gfx1153;gfx1200;gfx1201"
  CACHE STRING "AMD GPU architectures of the HIP backend if the build machine has no AMD GPU")
IF(ENABLE_HIP AND OPENEMS_HIP_AS_CUDA)
  IF(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    find_program(NVIDIA_SMI nvidia-smi)
    set(_nvidia_gpus "")
    IF(NVIDIA_SMI)
      execute_process(COMMAND ${NVIDIA_SMI} -L OUTPUT_VARIABLE _nvidia_gpus ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    ENDIF()
    IF(_nvidia_gpus MATCHES "GPU")
      set(CMAKE_CUDA_ARCHITECTURES native)
    ELSE()
      # every major architecture the CUDA compiler knows
      set(CMAKE_CUDA_ARCHITECTURES all-major)
    ENDIF()
  ENDIF()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  # the same rounding as the HIP flags below give on the other systems
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --fmad=false")
  IF(ENABLE_FLUSH_TO_ZERO)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ftz=true")
  ENDIF()
  message(STATUS "HIP GPU backend enabled as CUDA sources, architectures: ${CMAKE_CUDA_ARCHITECTURES}")
ELSEIF(ENABLE_HIP)
  set(_hip_nvidia OFF)
  IF((CMAKE_HIP_PLATFORM STREQUAL "nvidia") OR ("$ENV{HIP_PLATFORM}" STREQUAL "nvidia"))
    set(_hip_nvidia ON)
  ENDIF()
  # the AMD GPUs of this machine, for the default and for a given "native"
  set(_amd_gpus "")
  IF(NOT _hip_nvidia AND ((NOT DEFINED CMAKE_HIP_ARCHITECTURES) OR ("native" IN_LIST CMAKE_HIP_ARCHITECTURES)))
    # amdgpu-arch prints them one per line
    find_program(AMDGPU_ARCH amdgpu-arch HINTS ${ROCM_PATH}/llvm/bin $ENV{ROCM_PATH}/llvm/bin
                 ${ROCM_PATH}/lib/llvm/bin $ENV{ROCM_PATH}/lib/llvm/bin /opt/rocm/llvm/bin)
    IF(AMDGPU_ARCH)
      execute_process(COMMAND ${AMDGPU_ARCH} OUTPUT_VARIABLE _amd_gpus ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    ENDIF()
    IF(_amd_gpus MATCHES "gfx")
      string(REPLACE "\n" ";" _amd_gpus "${_amd_gpus}")
      list(REMOVE_DUPLICATES _amd_gpus)
    ELSE()
      set(_amd_gpus "")
    ENDIF()
  ENDIF()
  IF(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    IF(_hip_nvidia)
      set(CMAKE_HIP_ARCHITECTURES native)
    ELSEIF(_amd_gpus)
      # by name rather than native: the backend checks a device against the list
      set(CMAKE_HIP_ARCHITECTURES "${_amd_gpus}")
    ELSE()
      set(CMAKE_HIP_ARCHITECTURES "${OPENEMS_HIP_AMD_ARCHITECTURES}")
    ENDIF()
  ENDIF()
  enable_language(HIP)
  set(CMAKE_HIP_STANDARD 17)
  # No fused multiply-add and the same denormal handling, so the kernels round like the
  # CPU engines (python/Tests/GPU_Engine.py compares them value by value). The flags are
  # those of the compiler behind hipcc: hipcc on the NVIDIA platform, clang on AMD, and
  # the clang spelling would only reach the host compiler under hipcc.
  IF(CMAKE_HIP_PLATFORM STREQUAL "nvidia")
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} --fmad=false")
    IF(ENABLE_FLUSH_TO_ZERO)
      set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} --ftz=true")
    ENDIF()
  ELSE()
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -ffp-contract=off")
    IF(ENABLE_FLUSH_TO_ZERO)
      set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -fgpu-flush-denormals-to-zero")
    ENDIF()
  ENDIF()
  IF(CMAKE_HIP_PLATFORM STREQUAL "amd")
    # On the AMD platform HIP is a runtime library the kernels need at link time. Not
    # find_package(hip): its config file requires the ROCm clang as the C++ compiler,
    # while openEMS builds its C++ with the system one. On the NVIDIA platform hipcc
    # compiles to CUDA and the CUDA runtime comes with the HIP sources.
    find_library(HIP_RUNTIME_LIBRARY amdhip64 HINTS ${ROCM_PATH}/lib $ENV{ROCM_PATH}/lib /opt/rocm/lib)
    IF(NOT HIP_RUNTIME_LIBRARY)
      message(WARNING "The HIP GPU backend needs libamdhip64 on AMD (set ROCM_PATH), building without it")
      set(ENABLE_HIP OFF)
    ENDIF()
  ENDIF()
ENDIF()
IF(ENABLE_HIP AND NOT OPENEMS_HIP_AS_CUDA)
  message(STATUS "HIP GPU backend enabled, architectures: ${CMAKE_HIP_ARCHITECTURES}")
ENDIF()

# The HIP backend is a plugin loaded at runtime (FDTD/gpu_plugin.h): libopenEMS then
# needs no GPU runtime, a missing driver or device only makes the backend unavailable.
# OFF links it into libopenEMS, e.g. for GPU debuggers and profilers.
OPTION(OPENEMS_GPU_PLUGINS "Build the HIP GPU backend as a plugin loaded at runtime" ON)
IF(ENABLE_HIP AND NOT OPENEMS_GPU_PLUGINS)
  ADD_DEFINITIONS( -DOPENEMS_WITH_HIP )
ENDIF()

# SYCL GPU backend for the GPU engine (--engine=gpu), for Intel GPUs through Level Zero.
# It needs the Intel oneAPI DPC++ compiler (icpx) and is only built as a plugin:
# libopenEMS then needs no SYCL runtime. With icpx as the C++ compiler of the whole
# build the plugin is a target of it; with another one (e.g. GCC) icpx builds the plugin
# as a project of its own (FDTD/sycl) against the libopenEMS just built, as the HIP
# plugin is built by hipcc: OPENEMS_ICPX is that icpx, a path that does not exist turns
# it off. Whatever is missing switches the backend off, with a warning if it was asked for.
include(CheckCXXSourceCompiles)
set(SYCL_MISSING "")
set(_sycl_subbuild OFF)
IF(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_REQUIRED_FLAGS "-fsycl")
  set(CMAKE_REQUIRED_LINK_OPTIONS "-fsycl")
  check_cxx_source_compiles("#include <sycl/sycl.hpp>
int main() { sycl::range<1> r(1); return (int)r.size() - 1; }" OPENEMS_SYCL_COMPILES)
  unset(CMAKE_REQUIRED_FLAGS)
  unset(CMAKE_REQUIRED_LINK_OPTIONS)
  IF(NOT OPENEMS_SYCL_COMPILES)
    set(SYCL_MISSING "a working -fsycl compiler")
  ENDIF()
ELSE()
  # searched once: a value given by the user, even OPENEMS_ICPX-NOTFOUND, is kept
  IF(NOT DEFINED OPENEMS_ICPX)
    find_program(OPENEMS_ICPX icpx HINTS ${ONEAPI_ROOT}/compiler/latest/bin
                 $ENV{ONEAPI_ROOT}/compiler/latest/bin /opt/intel/oneapi/compiler/latest/bin)
  ENDIF()
  IF(OPENEMS_ICPX AND EXISTS "${OPENEMS_ICPX}")
    set(_sycl_subbuild ON)
  ELSE()
    set(SYCL_MISSING "the Intel oneAPI DPC++ compiler (icpx)")
  ENDIF()
ENDIF()
IF(SYCL_MISSING)
  OPTION(ENABLE_SYCL "Build the SYCL GPU backend (Intel GPUs)" OFF)
ELSE()
  OPTION(ENABLE_SYCL "Build the SYCL GPU backend (Intel GPUs)" ON)
ENDIF()
IF(ENABLE_SYCL AND SYCL_MISSING)
  message(WARNING "The SYCL GPU backend needs ${SYCL_MISSING}, building without it")
  set(ENABLE_SYCL OFF)
ENDIF()
IF(ENABLE_SYCL AND NOT OPENEMS_GPU_PLUGINS)
  message(WARNING "The SYCL GPU backend is only built as a plugin (OPENEMS_GPU_PLUGINS), building without it")
  set(ENABLE_SYCL OFF)
ENDIF()
# The Intel GPUs compiled ahead of time; any other Intel GPU compiles the SPIR-V image
# when the GPU engine starts.
set(OPENEMS_SYCL_AOT_DEVICES "bmg-g21" CACHE STRING
  "Intel GPU targets compiled ahead of time by the SYCL backend (ocloc -device syntax); empty: SPIR-V only, compiled on the device at run time")
IF(ENABLE_SYCL AND _sycl_subbuild)
  message(STATUS "SYCL GPU backend enabled, built by ${OPENEMS_ICPX} in the sub-build sycl-plugin, AOT devices: ${OPENEMS_SYCL_AOT_DEVICES}")
ELSEIF(ENABLE_SYCL)
  message(STATUS "SYCL GPU backend enabled, AOT devices: ${OPENEMS_SYCL_AOT_DEVICES}")
ENDIF()

# libopenEMS refuses a plugin of another build (see FDTD/gpu_plugin.h). In a git
# checkout the version does not tell uncommitted changes apart: a random id per build
# directory, kept across reconfigures. A release (no .git) is told by its version
# alone, so that its build is reproducible.
IF(NOT OPENEMS_GPU_BUILD_ID)
  IF(EXISTS ${PROJECT_SOURCE_DIR}/.git)
    string(RANDOM LENGTH 12 _gpu_build_id)
  ELSE()
    set(_gpu_build_id release)
  ENDIF()
  set(OPENEMS_GPU_BUILD_ID ${_gpu_build_id} CACHE INTERNAL "id of this build, shared by libopenEMS and its GPU plugins")
ENDIF()
configure_file(FDTD/gpu_build_id.h.in ${CMAKE_CURRENT_BINARY_DIR}/gpu_build_id.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Metal GPU backend for the GPU engine (--engine=gpu), macOS only
IF(APPLE)
  # Objective-C++ is a language of CMake since 3.16
  IF(CMAKE_VERSION VERSION_LESS 3.16)
    OPTION(ENABLE_METAL "Build the Metal GPU backend" OFF)
    IF(ENABLE_METAL)
      message(WARNING "The Metal GPU backend needs CMake 3.16 or newer, building without it")
      set(ENABLE_METAL OFF)
    ENDIF()
  ELSE()
    OPTION(ENABLE_METAL "Build the Metal GPU backend" ON)
  ENDIF()
ENDIF()
IF(APPLE AND ENABLE_METAL)
  message(STATUS "Metal GPU backend enabled")
  enable_language(OBJCXX)
  set(CMAKE_OBJCXX_STANDARD ${CMAKE_CXX_STANDARD})
  set(CMAKE_OBJCXX_FLAGS "${CMAKE_OBJCXX_FLAGS} -fobjc-arc")
  ADD_DEFINITIONS( -DOPENEMS_WITH_METAL )
ENDIF()

# libs
ADD_SUBDIRECTORY( tools )
ADD_SUBDIRECTORY( FDTD )
ADD_SUBDIRECTORY( FDTD/extensions )
ADD_SUBDIRECTORY( Common )

INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )    #find tools

message(STATUS "Sources:"  ${SOURCES})

IF(ENABLE_HIP AND OPENEMS_HIP_AS_CUDA)
  # the backend sources are HIP, compiled by nvcc as CUDA (see the HIP block above);
  # nvcc does not know the .hip suffix, hence -x cu
  set_source_files_properties(${HIP_SOURCES} PROPERTIES LANGUAGE CUDA COMPILE_OPTIONS "-x;cu")
ENDIF()

add_library( openEMS SHARED ${SOURCES})

# the target the HIP sources are built into: libopenEMS or the plugin
IF(ENABLE_HIP)
  IF(OPENEMS_HIP_AS_CUDA OR (CMAKE_HIP_PLATFORM STREQUAL "nvidia"))
    set(GPU_PLUGIN_NAME hip_nvidia)
  ELSE()
    set(GPU_PLUGIN_NAME hip_amd)
  ENDIF()
ENDIF()
IF(ENABLE_HIP AND OPENEMS_GPU_PLUGINS)
  set(HIP_TARGET openEMS_gpu_${GPU_PLUGIN_NAME})
  message(STATUS "HIP GPU backend built as the plugin ${HIP_TARGET}")
  add_library( ${HIP_TARGET} MODULE ${HIP_PLUGIN_SOURCES} )
  # installed to openEMS/plugins next to libopenEMS, where it is looked for; in the
  # build tree, where the openEMS binary takes that name, OPENEMS_PLUGIN_PATH=plugins
  set_target_properties( ${HIP_TARGET} PROPERTIES PREFIX ""
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugins
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/plugins )
  TARGET_LINK_LIBRARIES( ${HIP_TARGET} openEMS )
ELSEIF(ENABLE_HIP)
  set(HIP_TARGET openEMS)
ENDIF()
IF(ENABLE_HIP)
  # the name of the backend in the device list, e.g. "hip-amd"
  string(REPLACE "_" "-" _gpu_plugin_label ${GPU_PLUGIN_NAME})
  target_compile_definitions( ${HIP_TARGET} PRIVATE OPENEMS_GPU_PLUGIN_NAME="${_gpu_plugin_label}" )
  IF(NOT OPENEMS_HIP_AS_CUDA AND CMAKE_HIP_PLATFORM STREQUAL "amd")
    # the AMD runtime may abort on a device the build has no code for, see GPU_Backend_HIP::New()
    # "native" are the GPUs of this machine; no check if amdgpu-arch did not name them
    set(_hip_archs ${CMAKE_HIP_ARCHITECTURES})
    IF("native" IN_LIST _hip_archs)
      list(REMOVE_ITEM _hip_archs native)
      IF(_amd_gpus)
        list(APPEND _hip_archs ${_amd_gpus})
      ELSE()
        set(_hip_archs "")
      ENDIF()
    ENDIF()
    IF(_hip_archs)
      # comma separated: a semicolon would split the definition
      string(REPLACE ";" "," _hip_archs "${_hip_archs}")
      target_compile_definitions( ${HIP_TARGET} PRIVATE OPENEMS_HIP_ARCHITECTURES="${_hip_archs}" )
    ENDIF()
    IF(WIN32)
      # the device compile ignores the dllimport of every exported class, once per architecture
      target_compile_options( ${HIP_TARGET} PRIVATE $<$<COMPILE_LANGUAGE:HIP>:-Wno-ignored-attributes> )
    ENDIF()
  ENDIF()
ENDIF()
IF(ENABLE_HIP AND OPENEMS_HIP_AS_CUDA)
  target_include_directories( ${HIP_TARGET} PRIVATE ${HIP_INCLUDE_DIR} )
  target_compile_definitions( ${HIP_TARGET} PRIVATE __HIP_PLATFORM_NVIDIA__=1 )
ENDIF()
IF(ENABLE_SYCL AND NOT _sycl_subbuild)
  ADD_SUBDIRECTORY( FDTD/sycl )
ELSEIF(ENABLE_SYCL)
  # FDTD/sycl as a project of its own, compiled by icpx with the include directories and
  # definitions of libopenEMS; its plugin lands in plugins/ as the in-tree one does
  include(ExternalProject)
  get_directory_property(_sycl_includes INCLUDE_DIRECTORIES)
  get_directory_property(_sycl_definitions COMPILE_DEFINITIONS)
  # GIT_VERSION: its quotes do not survive the command line, and the plugin has no use for it
  set(_sycl_defs "")
  foreach(_def ${_sycl_definitions})
    IF(NOT _def MATCHES "^GIT_VERSION=")
      list(APPEND _sycl_defs ${_def})
    ENDIF()
  endforeach()
  string(REPLACE ";" "|" _sycl_includes "${_sycl_includes}")
  string(REPLACE ";" "|" _sycl_defs "${_sycl_defs}")
  set(_sycl_args
    -DCMAKE_CXX_COMPILER=${OPENEMS_ICPX}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DOPENEMS_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
    -DOPENEMS_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}
    -DOPENEMS_LIBRARY=$<TARGET_FILE:openEMS>
    -DOPENEMS_INCLUDE_DIRS=${_sycl_includes}
    -DOPENEMS_COMPILE_DEFINITIONS=${_sycl_defs}
    -DOPENEMS_SYCL_AOT_DEVICES=${OPENEMS_SYCL_AOT_DEVICES}
    -DENABLE_FLUSH_TO_ZERO=${ENABLE_FLUSH_TO_ZERO})
  IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # the GCC of libopenEMS, whose libstdc++ the plugin is compiled against (icpx would
    # take the newest GCC it finds)
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -print-libgcc-file-name
      OUTPUT_VARIABLE _libgcc OUTPUT_STRIP_TRAILING_WHITESPACE)
    get_filename_component(_gcc_install_dir "${_libgcc}" DIRECTORY)
    list(APPEND _sycl_args -DOPENEMS_GCC_INSTALL_DIR=${_gcc_install_dir})
  ENDIF()
  # BUILD_ALWAYS: the sub-build sees changes of the SYCL sources, the outer one does not
  ExternalProject_Add( openEMS_gpu_sycl
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/FDTD/sycl
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/sycl-plugin
    DEPENDS openEMS
    LIST_SEPARATOR |
    CMAKE_ARGS ${_sycl_args}
    BUILD_ALWAYS ON
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/plugins/openEMS_gpu_sycl.so
    INSTALL_COMMAND "" )
ENDIF()
#ADD_EXECUTABLE( openEMS main.cpp ${SOURCES})
IF(ENABLE_HIP AND NOT OPENEMS_HIP_AS_CUDA AND CMAKE_HIP_PLATFORM STREQUAL "amd")
  TARGET_LINK_LIBRARIES( ${HIP_TARGET} ${HIP_RUNTIME_LIBRARY} )
ENDIF()
# the plugin loader (FDTD/gpu_plugins.cpp): dlopen, and std::filesystem, which GCC 8
# keeps in a library of its own
TARGET_LINK_LIBRARIES( openEMS ${CMAKE_DL_LIBS} )
IF(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
  TARGET_LINK_LIBRARIES( openEMS stdc++fs )
ENDIF()
IF(APPLE AND ENABLE_METAL)
  TARGET_LINK_LIBRARIES( openEMS "-framework Metal" "-framework Foundation" )
ENDIF()

set_target_properties(openEMS PROPERTIES VERSION ${LIB_VERSION_STRING} SOVERSION ${LIB_VERSION_MAJOR} )
TARGET_LINK_LIBRARIES( openEMS
  ${CSXCAD_LIBRARIES}
  ${fparser_LIBRARIES}
  ${TinyXML_LIBRARY}
  ${HDF5_LIBRARIES}
  ${HDF5_HL_LIBRARIES}
  ${Boost_LIBRARIES}
  ${ZLIB_LIBRARIES}
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
    INSTALL(TARGETS openEMS RUNTIME DESTINATION bin ARCHIVE DESTINATION lib${LIB_SUFFIX})
else()
    INSTALL(TARGETS openEMS DESTINATION lib${LIB_SUFFIX})
endif()
INSTALL(TARGETS openEMS_bin DESTINATION bin)
IF(ENABLE_HIP AND OPENEMS_GPU_PLUGINS)
  # next to libopenEMS, where it is looked for: bin on Windows (a module is a
  # LIBRARY to install() there too)
  IF(WIN32)
    INSTALL(TARGETS ${HIP_TARGET} LIBRARY DESTINATION bin/openEMS/plugins)
  ELSE()
    INSTALL(TARGETS ${HIP_TARGET} LIBRARY DESTINATION lib${LIB_SUFFIX}/openEMS/plugins)
  ENDIF()
ENDIF()

IF(ENABLE_SYCL AND _sycl_subbuild)
  # the file of the sub-build; the in-tree target installs itself (FDTD/sycl)
  INSTALL(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/plugins/openEMS_gpu_sycl.so
    DESTINATION lib${LIB_SUFFIX}/openEMS/plugins)
ENDIF()

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
