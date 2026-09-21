# Metal GPU backend (Objective-C++)
if(APPLE AND ENABLE_METAL)
  set(METAL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/gpu_backend_metal.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_excitation.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_upml.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_mur_abc.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_lorentzmaterial.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_lumpedRLC.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_tfsf.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_absorbing_bc.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_steadystate.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_ext_cylinder.mm
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_multigrid.mm
  )
endif()

# HIP GPU backend
if(ENABLE_HIP)
  set(HIP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/gpu_backend_hip.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_excitation.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_upml.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_mur_abc.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_lorentzmaterial.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_lumpedRLC.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_tfsf.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_absorbing_bc.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_steadystate.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_cylinder.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_multigrid.hip
  )
  # the top level needs the list by itself: source file properties belong to the
  # directory of the target, and the openEMS target is defined there
  set(HIP_SOURCES ${HIP_SOURCES} PARENT_SCOPE)
endif()

set(SOURCES
  ${SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_multithread.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_cylinder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_cylinder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_sse_compressed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_sse_compressed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_multithread.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/excitation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_cylindermultigrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_cylindermultigrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_sse_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_cylindrical_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_gpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_gpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_gpu_cylindermultigrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_coeff_sets.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend_reference.cpp
  ${METAL_SOURCES}
  ${HIP_SOURCES}
  PARENT_SCOPE
)

