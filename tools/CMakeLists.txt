
set(SOURCES
  AdrOp.cpp
  ErrorMsg.cpp
  array_ops.cpp
  global.cpp
  hdf5_file_reader.cpp
  hdf5_file_writer.cpp
  sar_calculation.cpp
  useful.cpp
  vtk_file_writer.cpp
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
add_library(tools ${SOURCES} )

