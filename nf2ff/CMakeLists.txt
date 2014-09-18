
PROJECT(nf2ff CXX)
cmake_minimum_required(VERSION 2.8)

set(SOURCES
  main.cpp
	nf2ff.cpp
	nf2ff_calc.cpp
)

set(HEADERS
  nf2ff.h
  nf2ff_calc.h
)

ADD_EXECUTABLE( nf2ff ${SOURCES} )

TARGET_LINK_LIBRARIES( nf2ff
  tools
  tinyxml
  ${HDF5_LIBRARIES}
  ${Boost_LIBRARIES}
)

INSTALL(TARGETS nf2ff DESTINATION bin)


#TODO tarball, debug, release
