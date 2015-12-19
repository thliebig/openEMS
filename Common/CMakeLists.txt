
#INCLUDE_DIRECTORIES( ${openEMS_SOURCE_DIR} )

set(SOURCES
  ${SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/engine_interface_base.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/operator_base.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processcurrent.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processfieldprobe.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processfields.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processfields_fd.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processfields_sar.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processfields_td.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processing.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processintegral.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processmodematch.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/processvoltage.cpp
  PARENT_SCOPE
)

# Common lib
#add_library( Common STATIC ${SOURCES} )

