
INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )
INCLUDE_DIRECTORIES( ${CSXCAD_SOURCE_DIR}/src )

set(SOURCES
  engine_extension.cpp
  operator_ext_dispersive.cpp
  operator_ext_lorentzmaterial.cpp
  operator_ext_conductingsheet.cpp
  engine_ext_dispersive.cpp
  engine_ext_lorentzmaterial.cpp
  engine_ext_cylindermultigrid.cpp
  operator_ext_upml.cpp
  engine_ext_upml.cpp
  operator_extension.cpp
  engine_ext_mur_abc.cpp
  operator_ext_mur_abc.cpp
  operator_ext_cylinder.cpp
  engine_ext_cylinder.cpp
  operator_ext_excitation.cpp
  engine_ext_excitation.cpp
  operator_ext_tfsf.cpp
  engine_ext_tfsf.cpp
  operator_ext_steadystate.cpp
  engine_ext_steadystate.cpp
)

# FDTD/extensions lib
add_library( extensions ${SOURCES} )

