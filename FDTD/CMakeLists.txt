# HIP GPU backend
if(ENABLE_HIP)
  set(HIP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/gpu_backend_hip.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_cell_model.hip
    ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_ext_ops.hip
  )
  # the top level needs the list by itself: source file properties belong to the
  # directory of the target, and the openEMS and plugin targets are defined there
  set(HIP_SOURCES ${HIP_SOURCES} PARENT_SCOPE)
  if(OPENEMS_GPU_PLUGINS)
    set(HIP_PLUGIN_SOURCES ${HIP_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/hip/hip_plugin.cpp PARENT_SCOPE)
    # in a plugin, not in libopenEMS
    set(HIP_SOURCES "")
  endif()
endif()

# Metal GPU backend (Objective-C++)
if(APPLE AND ENABLE_METAL)
  set(METAL_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/gpu_backend_metal.mm
  )

  # All kernels of the Metal backend in Metal Shading Language: the shared code and the entry
  # points, compiled by the backend at run time from metal_cm_source.h. metal_cm_kernels.metal
  # is the same source for an offline compile: `cmake --build . --target metal_check` (CI).
  set(METAL_CM_PARTS
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels/cell_update.h
    ${CMAKE_CURRENT_SOURCE_DIR}/kernels/ext_ops.h
    ${CMAKE_CURRENT_SOURCE_DIR}/metal/metal_kernels.metal)
  set(METAL_CM_SOURCE "#include <metal_stdlib>\nusing namespace metal;\n#pragma STDC FP_CONTRACT OFF\n")
  foreach(part ${METAL_CM_PARTS})
    file(READ ${part} part_text)
    string(APPEND METAL_CM_SOURCE "${part_text}\n")
  endforeach()
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${METAL_CM_PARTS})
  file(WRITE ${CMAKE_BINARY_DIR}/metal_cm_kernels.metal "${METAL_CM_SOURCE}")
  file(WRITE ${CMAKE_BINARY_DIR}/metal_cm_source.h
       "// generated from ${METAL_CM_PARTS}\nstatic const char* METAL_CM_SOURCE = R\"CMSRC(${METAL_CM_SOURCE})CMSRC\";\n")
  add_custom_target(metal_check
    COMMAND xcrun -sdk macosx metal -std=metal3.0 -c ${CMAKE_BINARY_DIR}/metal_cm_kernels.metal -o ${CMAKE_BINARY_DIR}/metal_cm_kernels.air
    COMMENT "Compiling the Metal kernels offline"
    VERBATIM)
endif()

set(SOURCES
  ${SOURCES}
  ${HIP_SOURCES}
  ${METAL_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_threadpool.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_plan.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/cell_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_cell_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_ext.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_cpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_cpu_kernels.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_coeffs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_cylinder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_sse.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/excitation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_cylindermultigrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_subgrid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_sse_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_gpu_fdtd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_field_gather.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_gpu.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_plugins.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/gpu_backend_reference.cpp
  PARENT_SCOPE
)
