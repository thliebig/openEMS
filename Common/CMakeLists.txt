
INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )

set( SOURCES
  engine_interface_base.cpp
  operator_base.cpp
  processcurrent.cpp
  processfieldprobe.cpp
  processfields.cpp
  processfields_fd.cpp
  processfields_sar.cpp
  processfields_td.cpp
  processing.cpp
  processintegral.cpp
  processmodematch.cpp
  processvoltage.cpp
)

# Common lib
add_library( Common STATIC ${SOURCES} )

