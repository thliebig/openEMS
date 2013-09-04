# -------------------------------------------------
# Project created by QtCreator 2010-02-26T22:34:51
# -------------------------------------------------
TARGET = openEMS
CONFIG -= app_bundle qt
CONFIG += rtti exceptions
TEMPLATE = app
OBJECTS_DIR = obj
INCLUDEPATH += .
CONFIG += debug_and_release

#
# VERSION
#
VERSION=0.0.30

# add git revision
GITREV = $$system(git describe --tags)
isEmpty(GITREV):GITREV=$$VERSION
DEFINES += GIT_VERSION=\\\"$$GITREV\\\"

# remove unnecessary webkit define
DEFINES -= QT_WEBKIT

exists(localPaths.pri) {
    include(localPaths.pri)
}

###############################################################################
# CONFIG SECTION

# the SSE engine defaults to flush-to-zero mode, because of speed advantages
# to restore the correct handling of denormals and to comply to IEEE 754 uncomment:
# DEFINES += SSE_CORRECT_DENORMALS

# openEMS defaults to output length in unit meters; to recover the old behaviour
# to output length in terms of the drawing unit, uncomment:
# DEFINES += OUTPUT_IN_DRAWINGUNITS

# CONFIG SECTION
###############################################################################

win32 {
    CONFIG += console
    isEmpty(WIN32_LIB_ROOT) {
        WIN32_LIB_ROOT = ..
    }
    isEmpty(CSXCAD_ROOT) {
        CSXCAD_ROOT = $$WIN32_LIB_ROOT/CSXCAD
    }
    # CSXCAD
    INCLUDEPATH += $$CSXCAD_ROOT/include/CSXCAD
    LIBS += -L$$CSXCAD_ROOT/lib -lCSXCAD0

    # #3rd party libraries#
    # tinyxml
    INCLUDEPATH += $$WIN32_LIB_ROOT/tinyxml/include
    LIBS += -L$$WIN32_LIB_ROOT/tinyxml/bin -ltinyxml2
    DEFINES += TIXML_USE_STL
    # fparser
    INCLUDEPATH += $$WIN32_LIB_ROOT/fparser/include
    LIBS += -L$$WIN32_LIB_ROOT/fparser/bin -lfparser4
    # hdf5
    INCLUDEPATH += $$WIN32_LIB_ROOT/hdf5/include $$WIN32_LIB_ROOT/hdf5/include/cpp
    LIBS += -L$$WIN32_LIB_ROOT/hdf5/lib -lhdf5
    # zlib
    LIBS += -L$$WIN32_LIB_ROOT/zlib/lib -lz
    # boost
    DEFINES += BOOST_THREAD_USE_LIB
    INCLUDEPATH += $$WIN32_LIB_ROOT/boost/include
    LIBS += -L$$WIN32_LIB_ROOT/boost/lib -lboost_thread -lboost_chrono -lboost_system
    # vtk
    INCLUDEPATH +=   $$WIN32_LIB_ROOT/vtk/include/vtk-5.10
    LIBS += -L$$WIN32_LIB_ROOT/vtk/bin -lvtkCommon -lvtkIO -lvtkFiltering
}
!win32 {
    # CSXCAD
    isEmpty(CSXCAD_ROOT) {
        CSXCAD_ROOT = /usr
    } else {
        QMAKE_LFLAGS += \'-Wl,-rpath,$$CSXCAD_ROOT/lib\'
    }
    INCLUDEPATH += $$CSXCAD_ROOT/include/CSXCAD
    LIBS += -L$$CSXCAD_ROOT/lib -lCSXCAD

    # #3rd party libraries#
    #fparser
    isEmpty(FPARSER_ROOT) {
        FPARSER_ROOT = /usr
    } else {
        INCLUDEPATH += $$FPARSER_ROOT/include
        LIBS += -L$$FPARSER_ROOT/lib
        QMAKE_LFLAGS += \'-Wl,-rpath,$$FPARSER_ROOT/lib\'
    }
    LIBS += -lfparser

    LIBS += -ltinyxml
    DEFINES += TIXML_USE_STL
    LIBS += -lboost_thread-mt
    # hdf5 (and mpi for parallel hdf5)
    LIBS += -lhdf5_hl -lhdf5
    LIBS += -lmpi -lmpi_cxx
    INCLUDEPATH += /usr/include/mpi
    ### vtk ###
    INCLUDEPATH += /usr/include/vtk-5.2 \
        /usr/include/vtk-5.4 \
        /usr/include/vtk-5.6 \
        /usr/include/vtk-5.8 \
        /usr/include/vtk-5.10 \
        /usr/include/vtk
    LIBS += -lvtkCommon \
        -lvtkIO \
        -lvtksys \
        -lvtkFiltering
}

# vtk includes deprecated header files; silence the corresponding warning
QMAKE_CXXFLAGS += -Wno-deprecated

# hdf5 compat
DEFINES += H5_USE_16_API


#### SOURCES ################################################################
SOURCES += main.cpp \
    openems.cpp

# FDTD
SOURCES += FDTD/engine.cpp \
    FDTD/operator.cpp \
    FDTD/engine_multithread.cpp \
    FDTD/operator_cylinder.cpp \
    FDTD/engine_cylinder.cpp \
    FDTD/engine_sse.cpp \
    FDTD/operator_sse.cpp \
    FDTD/operator_sse_compressed.cpp \
    FDTD/engine_sse_compressed.cpp \
    FDTD/operator_multithread.cpp \
    FDTD/excitation.cpp \
    FDTD/operator_cylindermultigrid.cpp \
    FDTD/engine_cylindermultigrid.cpp \
    FDTD/engine_interface_fdtd.cpp \
    FDTD/engine_interface_sse_fdtd.cpp \
    FDTD/engine_interface_cylindrical_fdtd.cpp

# FDTD/extensions source files
SOURCES += FDTD/extensions/engine_extension.cpp \
    FDTD/extensions/operator_ext_dispersive.cpp \
    FDTD/extensions/operator_ext_lorentzmaterial.cpp \
    FDTD/extensions/operator_ext_conductingsheet.cpp \
    FDTD/extensions/engine_ext_dispersive.cpp \
    FDTD/extensions/engine_ext_lorentzmaterial.cpp \
    FDTD/extensions/engine_ext_cylindermultigrid.cpp \
    FDTD/extensions/operator_ext_upml.cpp \
    FDTD/extensions/engine_ext_upml.cpp \
    FDTD/extensions/operator_extension.cpp \
    FDTD/extensions/engine_ext_mur_abc.cpp \
    FDTD/extensions/operator_ext_mur_abc.cpp \
    FDTD/extensions/operator_ext_cylinder.cpp \
    FDTD/extensions/engine_ext_cylinder.cpp \
    FDTD/extensions/operator_ext_excitation.cpp \
    FDTD/extensions/engine_ext_excitation.cpp \
    FDTD/extensions/operator_ext_tfsf.cpp \
    FDTD/extensions/engine_ext_tfsf.cpp

# Common source files
SOURCES += Common/operator_base.cpp \
    Common/engine_interface_base.cpp \
    Common/processmodematch.cpp \
    Common/processvoltage.cpp \
    Common/processing.cpp \
    Common/processintegral.cpp \
    Common/processfields.cpp \
    Common/processfields_td.cpp \
    Common/processcurrent.cpp \
    Common/processfields_fd.cpp \
    Common/processfieldprobe.cpp \
    Common/processfields_sar.cpp

# tools
SOURCES += tools/global.cpp \
    tools/useful.cpp \
    tools/array_ops.cpp \
    tools/ErrorMsg.cpp \
    tools/AdrOp.cpp \
    tools/sar_calculation.cpp \
    tools/vtk_file_writer.cpp \
    tools/hdf5_file_writer.cpp

#### HEADERS ################################################################
HEADERS += openems.h

# FDTD
HEADERS += FDTD/engine.h \
    FDTD/operator.h \
    FDTD/engine_multithread.h \
    FDTD/operator_cylinder.h \
    FDTD/engine_cylinder.h \
    FDTD/engine_sse.h \
    FDTD/operator_sse.h \
    FDTD/excitation.h \
    FDTD/operator_sse_compressed.h \
    FDTD/engine_sse_compressed.h \
    FDTD/operator_multithread.h \
    FDTD/operator_cylindermultigrid.h \
    FDTD/engine_cylindermultigrid.h \
    FDTD/engine_interface_fdtd.h \
    FDTD/engine_interface_sse_fdtd.h \
    FDTD/engine_interface_cylindrical_fdtd.h

# FDTD/extensions header files
HEADERS += FDTD/extensions/operator_extension.h \
    FDTD/extensions/engine_extension.h \
    FDTD/extensions/engine_ext_mur_abc.h \
    FDTD/extensions/operator_ext_mur_abc.h \
    FDTD/extensions/operator_ext_cylinder.h \
    FDTD/extensions/engine_ext_cylinder.h \
    FDTD/extensions/operator_ext_dispersive.h \
    FDTD/extensions/operator_ext_lorentzmaterial.h \
    FDTD/extensions/operator_ext_conductingsheet.h \
    FDTD/extensions/cond_sheet_parameter.h \
    FDTD/extensions/engine_ext_dispersive.h \
    FDTD/extensions/engine_ext_lorentzmaterial.h \
    FDTD/extensions/engine_ext_cylindermultigrid.h \
    FDTD/extensions/operator_ext_upml.h \
    FDTD/extensions/engine_ext_upml.h \
    FDTD/extensions/operator_ext_excitation.h \
    FDTD/extensions/engine_ext_excitation.h \
    FDTD/extensions/operator_ext_tfsf.h \
    FDTD/extensions/engine_ext_tfsf.h

# Common header files
HEADERS += Common/operator_base.h \
    Common/engine_interface_base.h \
    Common/processvoltage.h \
    Common/processing.h \
    Common/processintegral.h \
    Common/processfields.h \
    Common/processfields_td.h \
    Common/processcurrent.h \
    Common/processmodematch.h \
    Common/processfields_fd.h \
    Common/processfieldprobe.h \
    Common/processfields_sar.h

# tools
HEADERS += tools/ErrorMsg.h \
    tools/AdrOp.h \
    tools/constants.h \
    tools/array_ops.h \
    tools/global.h \
    tools/useful.h \
    tools/aligned_allocator.h \
    tools/sar_calculation.h \
    tools/vtk_file_writer.h \
    tools/hdf5_file_writer.h

!packaging {
    # if packaging is not set in CONFIG, set some default flags
    # if packaging is enabled, give the flags on the qmake comandline
    QMAKE_CXXFLAGS_RELEASE = -O3 -g -march=native
    QMAKE_CXXFLAGS_DEBUG = -O0 -g -march=native
}

MPI_SUPPORT {
    DEFINES += MPI_SUPPORT
    QMAKE_CC     = mpicc
    QMAKE_CXX    = mpicxx
    QMAKE_LINK   = mpicxx
    QMAKE_LINK_C = mpicc
    HEADERS += FDTD/operator_mpi.h \
               FDTD/engine_mpi.h \
               FDTD/openems_fdtd_mpi.h
    SOURCES += FDTD/operator_mpi.cpp \
               FDTD/engine_mpi.cpp \
               FDTD/openems_fdtd_mpi.cpp

    QMAKE_CXXFLAGS_RELEASE += -Wno-unused-parameter #needed because mpich2 produces tons of unused parameter
}



#
# create tar file (for the whole openEMS project)
#
tarball.target = tarball
tarball.commands = git archive --format=tar --prefix=openEMS-$$VERSION/ HEAD | bzip2 > openEMS-$${VERSION}.tar.bz2
QMAKE_EXTRA_TARGETS += tarball


#
# INSTALL (only the openEMS executable and matlab scripts)
#
isEmpty(PREFIX) {
    PREFIX = /usr/local
}
install.target = install
install.commands = mkdir -p \"$$PREFIX/bin\"
install.commands += && mkdir -p \"$$PREFIX/share/openEMS/matlab\"
unix:install.commands += && cp -at \"$$PREFIX/bin/\" openEMS.sh openEMS_MPI.sh openEMS
win32:install.commands += && cp -at \"$$PREFIX/bin/\" release/openEMS.exe
install.commands += && cp -at \"$$PREFIX/share/openEMS/\" matlab/
QMAKE_EXTRA_TARGETS += install


#
# create .PHONY target
#
phony.target = .PHONY
phony.depends = $$QMAKE_EXTRA_TARGETS
QMAKE_EXTRA_TARGETS += phony
