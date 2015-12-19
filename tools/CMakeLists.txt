
set(SOURCES
  ${SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/AdrOp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ErrorMsg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/array_ops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/global.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_file_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_file_writer.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sar_calculation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/useful.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/vtk_file_writer.cpp
  PARENT_SCOPE
)

#set(HEADERS
#  constants.h
#  array_ops.h
#  global.h
#  useful.h
#  aligned_allocator.h
#  hdf5_file_reader.h
#  hdf5_file_writer.h
#)

# tools lib
#add_library(tools STATIC ${SOURCES} )

