# SYCL GPU backend (Intel GPUs), always a plugin (see the SYCL block of the top level).
# Added by the top level when icpx is the C++ compiler of the whole build; otherwise a
# project of its own, built by icpx against the libopenEMS just built (the top level
# passes the OPENEMS_* variables below and installs the plugin).
IF(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.1...3.10)
  # the GCC of libopenEMS (see below) for project() too, whose compiler identification
  # and check take CXXFLAGS: the GCC icpx picks by itself, the newest, need not have
  # C++ headers, and icpx would not be recognized (no -fPIC) or fail the check
  IF(OPENEMS_GCC_INSTALL_DIR)
    set(ENV{CXXFLAGS} "$ENV{CXXFLAGS} --gcc-install-dir=${OPENEMS_GCC_INSTALL_DIR}")
  ENDIF()
  project(openEMS_gpu_sycl CXX)
  set(CMAKE_CXX_STANDARD 17)
  set(_sycl_standalone ON)
  # what the openEMS sources are compiled with, from the top level
  include_directories(${OPENEMS_SOURCE_DIR} ${OPENEMS_BINARY_DIR} ${OPENEMS_INCLUDE_DIRS})
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${OPENEMS_COMPILE_DEFINITIONS})
  set(_openems_lib ${OPENEMS_LIBRARY})
ELSE()
  set(_sycl_standalone OFF)
  set(OPENEMS_BINARY_DIR ${openEMS_BINARY_DIR})
  set(_openems_lib openEMS)
ENDIF()

# the linker's dependency file names the temporary objects of the icpx driver, gone
# after the link: every build would link the plugin again
set(CMAKE_LINK_DEPENDS_USE_LINKER FALSE)

add_library( openEMS_gpu_sycl MODULE
  gpu_backend_sycl.cpp
  sycl_ext_excitation.cpp
  sycl_ext_upml.cpp
  sycl_ext_mur_abc.cpp
  sycl_ext_lorentzmaterial.cpp
  sycl_ext_debyematerial.cpp
  sycl_ext_lumpedRLC.cpp
  sycl_ext_tfsf.cpp
  sycl_ext_absorbing_bc.cpp
  sycl_ext_steadystate.cpp
  sycl_ext_cylinder.cpp
  sycl_multigrid.cpp
  sycl_plugin.cpp
)
# installed to openEMS/plugins next to libopenEMS, where it is looked for; in the
# build tree, where the openEMS binary takes that name, OPENEMS_PLUGIN_PATH=plugins
set_target_properties( openEMS_gpu_sycl PROPERTIES PREFIX ""
  LIBRARY_OUTPUT_DIRECTORY ${OPENEMS_BINARY_DIR}/plugins )
TARGET_LINK_LIBRARIES( openEMS_gpu_sycl ${_openems_lib} )
target_compile_definitions( openEMS_gpu_sycl PRIVATE OPENEMS_GPU_PLUGIN_NAME="sycl" )

# The kernels round like the CPU engines: no fused multiply-add, correctly rounded
# division and square root, flush-to-zero as ENABLE_FLUSH_TO_ZERO. On the compile and
# the link line: the device code is generated at link time, where the
# -foffload-fp32-prec-* flags act; -fno-fast-math because the icpx host default is fast math.
IF(OPENEMS_SYCL_AOT_DEVICES)
  set(_sycl_targets spir64,spir64_gen)
ELSE()
  set(_sycl_targets spir64)
ENDIF()
set(_sycl_flags -fsycl -fsycl-targets=${_sycl_targets} -fno-fast-math -ffp-contract=off
  -foffload-fp32-prec-div -foffload-fp32-prec-sqrt)
# the libstdc++ of the GCC that built libopenEMS, not whichever GCC icpx would pick:
# the plugin shares its classes and standard library objects
IF(OPENEMS_GCC_INSTALL_DIR)
  list(APPEND _sycl_flags --gcc-install-dir=${OPENEMS_GCC_INSTALL_DIR})
ENDIF()
target_compile_options( openEMS_gpu_sycl PRIVATE ${_sycl_flags} )
target_link_options( openEMS_gpu_sycl PRIVATE ${_sycl_flags} )
# options of the device back ends, link line only (SHELL: keeps a quoted string one
# argument and the repeated -Xsycl-target-backend from being de-duplicated)
IF(OPENEMS_SYCL_AOT_DEVICES)
  IF(ENABLE_FLUSH_TO_ZERO)
    target_link_options( openEMS_gpu_sycl PRIVATE
      "SHELL:-Xsycl-target-backend=spir64_gen \"-device ${OPENEMS_SYCL_AOT_DEVICES} -options -cl-denorms-are-zero\"" )
  ELSE()
    target_link_options( openEMS_gpu_sycl PRIVATE
      "SHELL:-Xsycl-target-backend=spir64_gen \"-device ${OPENEMS_SYCL_AOT_DEVICES}\"" )
  ENDIF()
ENDIF()
IF(ENABLE_FLUSH_TO_ZERO)
  target_link_options( openEMS_gpu_sycl PRIVATE "SHELL:-Xsycl-target-backend=spir64 -cl-denorms-are-zero" )
ENDIF()

# The SYCL runtime (and the Intel math libraries icpx links) sit next to the compiler,
# which the driver links implicitly: no RPATH entry of CMake's own. With it in the
# RPATH the plugin loads without the oneAPI environment (setvars.sh).
get_filename_component(_icpx_dir ${CMAKE_CXX_COMPILER} REALPATH)
get_filename_component(_icpx_dir ${_icpx_dir} DIRECTORY)
find_library(OPENEMS_SYCL_LIBRARY sycl HINTS ${_icpx_dir}/../lib)
IF(OPENEMS_SYCL_LIBRARY)
  get_filename_component(_sycl_lib_dir ${OPENEMS_SYCL_LIBRARY} REALPATH)
  get_filename_component(_sycl_lib_dir ${_sycl_lib_dir} DIRECTORY)
  set_property(TARGET openEMS_gpu_sycl PROPERTY BUILD_RPATH ${_sycl_lib_dir})
  set_property(TARGET openEMS_gpu_sycl PROPERTY INSTALL_RPATH ${CMAKE_INSTALL_RPATH} ${_sycl_lib_dir})
ENDIF()

IF(_sycl_standalone)
  # the top level installs this build-tree file as it is: libopenEMS is loaded already
  # when the plugin is, so the build RPATH needs no more than the SYCL runtime
  set_property(TARGET openEMS_gpu_sycl PROPERTY BUILD_WITH_INSTALL_RPATH ON)
ELSE()
  INSTALL(TARGETS openEMS_gpu_sycl LIBRARY DESTINATION lib${LIB_SUFFIX}/openEMS/plugins)
ENDIF()
