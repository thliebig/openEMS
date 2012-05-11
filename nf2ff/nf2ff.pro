CONFIG -= app_bundle qt
TEMPLATE = app
OBJECTS_DIR = obj
CONFIG += debug_and_release

win32 {
    CONFIG += console
    INCLUDEPATH += ../../hdf5/include ../../hdf5/include/cpp ../../boost/include/boost-1_42
    INCLUDEPATH += ../../tinyxml
    LIBS +=  ../../hdf5/lib/hdf5.lib
    LIBS += ../../boost/lib/libboost_thread-mgw44-mt.lib
    LIBS += ../../tinyxml/release/libtinyxml2.a
}
!win32 {
    LIBS += -lboost_thread-mt
    LIBS += -lhdf5
    LIBS += -ltinyxml
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
install.target = install
install.commands = mkdir -p \"$(INSTALL_ROOT)/usr/bin\"
install.commands += && cp -at \"$(INSTALL_ROOT)/usr/bin/\" nf2ff
QMAKE_EXTRA_TARGETS += install


#
# create .PHONY target
#
phony.target = .PHONY
phony.depends = $$QMAKE_EXTRA_TARGETS
QMAKE_EXTRA_TARGETS += phony

