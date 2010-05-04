# -------------------------------------------------
# Project created by QtCreator 2010-02-26T22:34:51
# -------------------------------------------------
TARGET = openEMS
CONFIG += console
CONFIG -= app_bundle qt
TEMPLATE = app
OBJECTS_DIR = obj
INCLUDEPATH += .
INCLUDEPATH += ../CSXCAD \
    ../fparser \
    ../tinyxml
LIBS += -L../CSXCAD -lCSXCAD

win32 {
    INCLUDEPATH += ../hdf5/include ../boost/include/boost-1_42
    LIBS +=  ../hdf5/lib/libhdf5_cpp.a ../hdf5/lib/libhdf5.a
    LIBS += ../boost/lib/libboost_thread-mgw44-mt.lib
    LIBS += -L../CSXCAD/release
    LIBS += ../fparser/release/libfparser4.a
    LIBS += ../tinyxml/release/libtinyxml2.a
}
!win32 {
    LIBS += ../fparser/libfparser.so
    LIBS += ../tinyxml/libtinyxml.so
    LIBS += -lboost_thread
    LIBS += -lhdf5 -lhdf5_cpp
}

QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../CSXCAD\'
QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../fparser\'
QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../tinyxml\'
SOURCES += main.cpp \
    tools/ErrorMsg.cpp \
    tools/AdrOp.cpp \
    FDTD/engine.cpp \
    FDTD/operator.cpp \
    tools/array_ops.cpp \
    FDTD/processvoltage.cpp \
    FDTD/processing.cpp \
    FDTD/processfields.cpp \
    FDTD/processfields_td.cpp \
    FDTD/processcurrent.cpp \
    examples/FDTD_examples.cpp \
    openems.cpp \
    FDTD/engine_multithread.cpp \
    FDTD/operator_cylinder.cpp \
    FDTD/engine_cylinder.cpp \
    FDTD/engine_sse.cpp \
    FDTD/operator_sse.cpp \
    FDTD/operator_extension.cpp \
    FDTD/engine_extension.cpp \
    FDTD/engine_ext_mur_abc.cpp \
    FDTD/operator_ext_mur_abc.cpp \
    FDTD/excitation.cpp
HEADERS += tools/ErrorMsg.h \
    tools/AdrOp.h \
    tools/constants.h \
    FDTD/engine.h \
    FDTD/operator.h \
    tools/array_ops.h \
    FDTD/processvoltage.h \
    FDTD/processing.h \
    FDTD/processfields.h \
    FDTD/processfields_td.h \
    FDTD/processcurrent.h \
    examples/FDTD_examples.h \
    openems.h \
    FDTD/engine_multithread.h \
    FDTD/operator_cylinder.h \
    FDTD/engine_sse.h \
    FDTD/operator_sse.h \
    FDTD/engine_cylinder.h \
    FDTD/operator_extension.h \
    FDTD/engine_extension.h \
    FDTD/engine_ext_mur_abc.h \
    FDTD/operator_ext_mur_abc.h \
    FDTD/excitation.h
QMAKE_CXXFLAGS_RELEASE = -O3 \
    -g \
	-march=native
QMAKE_CXXFLAGS_DEBUG = -O0 \
    -g \
	-march=native

# add git revision
QMAKE_CXXFLAGS += -DGIT_VERSION=\\\"`git describe --tags`\\\"

# to use ABI2 target:
# qmake CONFIG+="ABI2 bits64" -o Makefile.ABI2-64 openEMS.pro
# make -fMakefile.ABI2-64
ABI2 { 
    CONFIG -= debug \
        debug_and_release
    CONFIG += release
    QMAKE_CFLAGS_RELEASE = -O2 \
        -fabi-version=2
    QMAKE_CXXFLAGS_RELEASE = -O2 \
        -fabi-version=2
    QMAKE_CC = apgcc
    QMAKE_CXX = apg++
    QMAKE_LINK = apg++
    QMAKE_LINK_SHLIB = apg++
    QMAKE_LFLAGS_RPATH = 
    QMAKE_LFLAGS = \'-Wl,-rpath,\$$ORIGIN/lib\'
}
bits64 { 
    QMAKE_CXXFLAGS_RELEASE += -m64 \
        -march=athlon64
    QMAKE_LFLAGS_RELEASE += -m64 \
        -march=athlon64
    OBJECTS_DIR = ABI2-64
    LIBS = ../CSXCAD/ABI2-64/libCSXCAD.so
    LIBS += ../fparser/ABI2-64/libfparser.so
    LIBS += ../tinyxml/ABI2-64/libtinyxml.so
    LIBS += ../boost-64/lib/libboost_thread.so
    LIBS += ../hdf5-64/lib/libhdf5.so
    LIBS += ../hdf5-64/lib/libhdf5_cpp.so \
        -lpthread
    INCLUDEPATH += ../hdf5-64/include
    INCLUDEPATH += ../boost-64/include
}
bits32 { 
    QMAKE_CXXFLAGS_RELEASE += -m32 \
        -march=i686
    QMAKE_LFLAGS_RELEASE += -m32 \
        -march=i686
    OBJECTS_DIR = ABI2-32
    LIBS = ../CSXCAD/ABI2-32/libCSXCAD.so
    LIBS += ../fparser/ABI2-32/libfparser.so
    LIBS += ../tinyxml/ABI2-32/libtinyxml.so
    LIBS += ../boost-32/lib/libboost_thread.so
    LIBS += ../hdf5-32/lib/libhdf5.so
    LIBS += ../hdf5-32/lib/libhdf5_cpp.so
    INCLUDEPATH += ../hdf5-32/include
    INCLUDEPATH += ../boost-32/include
}
ABI2 { 
    DESTDIR = $$OBJECTS_DIR
    MOC_DIR = $$OBJECTS_DIR
    UI_DIR = $$OBJECTS_DIR
    RCC_DIR = $$OBJECTS_DIR
}
