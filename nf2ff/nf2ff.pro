TARGET = nf2ff
CONFIG += console
CONFIG -= app_bundle qt
TEMPLATE = app
OBJECTS_DIR = obj
INCLUDEPATH += .
INCLUDEPATH += ../../tinyxml
CONFIG += debug_and_release

win32 {
    QMAKE_CXXFLAGS += -DH5_USE_16_API
    INCLUDEPATH += ../../hdf5/include ../../hdf5/include/cpp ../../boost/include/boost-1_42
	LIBS +=  ../../hdf5/lib/hdf5.lib
    LIBS += ../../boost/lib/libboost_thread-mgw44-mt.lib
	LIBS += ../../tinyxml/release/libtinyxml2.a
}
!win32 {
    LIBS += -lboost_thread
	LIBS += -lhdf5
	LIBS += ../../tinyxml/libtinyxml.so
}
QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../../tinyxml\'

TOOLSPATH = ../tools

#### SOURCES ################################################################
SOURCES += main.cpp \
	nf2ff.cpp \
	nf2ff_calc.cpp

# tools
 SOURCES += $$TOOLSPATH/global.cpp \
	$$TOOLSPATH/useful.cpp \
	$$TOOLSPATH/array_ops.cpp \
	$$TOOLSPATH/hdf5_file_reader.cpp \
	$$TOOLSPATH/hdf5_file_writer.cpp

#### HEADERS ################################################################
HEADERS += nf2ff.h \
	nf2ff_calc.h

# tools
HEADERS += $$TOOLSPATH/constants.h \
	$$TOOLSPATH/array_ops.h \
	$$TOOLSPATH/global.h \
	$$TOOLSPATH/useful.h \
	$$TOOLSPATH/aligned_allocator.h \
	$$TOOLSPATH/hdf5_file_reader.h \
	$$TOOLSPATH/hdf5_file_writer.h

QMAKE_CXXFLAGS_RELEASE = -O3 \
	-g \
	-march=native
QMAKE_CXXFLAGS_DEBUG = -O0 \
	-g \
	-march=native

# add git revision
# QMAKE_CXXFLAGS += -DGIT_VERSION=\\\"`git describe --tags`\\\"


