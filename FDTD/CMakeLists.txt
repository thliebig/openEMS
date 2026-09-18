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
  )
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
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend_reference.cpp
  ${METAL_SOURCES}
  PARENT_SCOPE
)

