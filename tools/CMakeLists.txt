
set(SOURCES
  ${SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/AdrOp.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ErrorMsg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/array_ops.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/signal.cpp
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



PROJECT(sar CXX)
cmake_minimum_required(VERSION 2.8)

ADD_EXECUTABLE( sar_bin sar_calc.cpp)
SET_TARGET_PROPERTIES(sar_bin PROPERTIES OUTPUT_NAME sar_calc)
TARGET_LINK_LIBRARIES(sar_bin openEMS)

INSTALL(TARGETS sar_bin DESTINATION bin)

