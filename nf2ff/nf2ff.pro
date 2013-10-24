CONFIG -= app_bundle qt
TEMPLATE = app
OBJECTS_DIR = obj
CONFIG += debug_and_release

VERSION = 0.1.0

exists(localPaths.pri) {
    include(localPaths.pri)
}

win32 {
    CONFIG += console

    isEmpty(WIN32_LIB_ROOT) {
	WIN32_LIB_ROOT = ../..
    }

    # #3rd party libraries#
    # tinyxml
    DEFINES += TIXML_USE_STL
    INCLUDEPATH += $$WIN32_LIB_ROOT/tinyxml/include
    LIBS += -L$$WIN32_LIB_ROOT/tinyxml/bin -ltinyxml2

    # hdf5
    INCLUDEPATH += $$WIN32_LIB_ROOT/hdf5/include
    LIBS += -L$$WIN32_LIB_ROOT/hdf5/lib -lhdf5
    # zlib
    LIBS += -L$$WIN32_LIB_ROOT/zlib/lib -lz

    # boost
    DEFINES += BOOST_THREAD_USE_LIB
    INCLUDEPATH += $$WIN32_LIB_ROOT/boost/include
    LIBS += -L$$WIN32_LIB_ROOT/boost/lib -lboost_thread -lboost_chrono -lboost_system
}
!win32 {
    LIBS += -lboost_thread -lboost_system
    LIBS += -ltinyxml
    #vtk
    isEmpty(VTK_LIBRARYPATH){
    } else {
    LIBS +=-L$$VTK_LIBRARYPATH
    }
    LIBS += -lhdf5_hl -lhdf5
}

# hdf5 compat
DEFINES += H5_USE_16_API

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

!packaging {
# if packaging is not set in CONFIG, set some default flags
# if packaging is enabled, give the flags on the qmake comandline
QMAKE_CXXFLAGS_RELEASE = -O3 \
	-g \
	-march=native
QMAKE_CXXFLAGS_DEBUG = -O0 \
	-g \
	-march=native
}

# add git revision
# QMAKE_CXXFLAGS += -DGIT_VERSION=\\\"`git describe --tags`\\\"




#
# INSTALL (only the nf2ff executable)
#
isEmpty(PREFIX) {
 PREFIX = /usr/local
}
install.target = install
install.commands = mkdir -p \"$$PREFIX/bin\"
unix:install.commands += && cp -at \"$$PREFIX/bin/\" nf2ff
win32:install.commands += && cp -at \"$$PREFIX/bin/\" release/nf2ff.exe
QMAKE_EXTRA_TARGETS += install


#
# create .PHONY target
#
phony.target = .PHONY
phony.depends = $$QMAKE_EXTRA_TARGETS
QMAKE_EXTRA_TARGETS += phony

