
INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )

set(SOURCES
  engine.cpp
  operator.cpp
  engine_multithread.cpp
  operator_cylinder.cpp
  engine_cylinder.cpp
  engine_sse.cpp
  operator_sse.cpp
  operator_sse_compressed.cpp
  engine_sse_compressed.cpp
  operator_multithread.cpp
  excitation.cpp
  operator_cylindermultigrid.cpp
  engine_cylindermultigrid.cpp
  engine_interface_fdtd.cpp
  engine_interface_sse_fdtd.cpp
  engine_interface_cylindrical_fdtd.cpp
)

# FDTD lib
add_library( FDTD ${SOURCES} )

