# -------------------------------------------------
# Project created by QtCreator 2010-02-26T22:34:51
# -------------------------------------------------
TARGET = openEMS
CONFIG -= app_bundle qt
TEMPLATE = app
OBJECTS_DIR = obj
INCLUDEPATH += .
INCLUDEPATH += ../CSXCAD ../fparser
CONFIG += debug_and_release


#
# VERSION
#
VERSION=0.0.28


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
    WIN32_LIB_ROOT = ..
    # tinyxml
    INCLUDEPATH += $$WIN32_LIB_ROOT/tinyxml
    LIBS += -L../tinyxml/release -ltinyxml2
    # fparser
    LIBS += -L../fparser/release -lfparser4
    # CSXCAD
    LIBS += -L../CSXCAD/release  -lCSXCAD0
    # hdf5
    INCLUDEPATH += $$WIN32_LIB_ROOT/hdf5/include $$WIN32_LIB_ROOT/hdf5/include/cpp
    LIBS += -L$$WIN32_LIB_ROOT/hdf5/lib -lhdf5
    # boost
    DEFINES += BOOST_THREAD_USE_LIB
    INCLUDEPATH += $$WIN32_LIB_ROOT/boost/include
    LIBS += $$WIN32_LIB_ROOT/boost/lib/libboost_thread-mgw44-mt.lib
    # vtk
    INCLUDEPATH +=   $$WIN32_LIB_ROOT/vtk \
        $$WIN32_LIB_ROOT/vtk/Common \
        $$WIN32_LIB_ROOT/vtk/Filtering \
        $$WIN32_LIB_ROOT/vtk/IO
     LIBS += -L$$WIN32_LIB_ROOT/vtk/bin -lvtkCommon -lvtkIO -lvtkFiltering
}
!win32 {
    LIBS += -L../fparser -lfparser
    LIBS += -ltinyxml
    LIBS += -L../CSXCAD -lCSXCAD
    LIBS += -lboost_thread-mt
    LIBS += -lhdf5 -lhdf5_cpp
    ### vtk ###
    INCLUDEPATH += /usr/include/vtk-5.2 \
        /usr/include/vtk-5.4 \
        /usr/include/vtk-5.6 \
        /usr/include/vtk-5.8 \
        /usr/include/vtk-5.10 \
        /usr/include/vtk
    INCLUDEPATH += /usr/include/CSXCAD
    LIBS += -lvtkCommon \
        -lvtkIO \
        -lvtksys \
        -lvtkFiltering
    QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../CSXCAD\'
    QMAKE_LFLAGS += \'-Wl,-rpath,\$$ORIGIN/../fparser\'
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
    FDTD/extensions/operator_ext_pml_sf.cpp \
    FDTD/extensions/engine_ext_pml_sf.cpp \
    FDTD/extensions/engine_ext_cylindermultigrid.cpp \
    FDTD/extensions/operator_ext_upml.cpp \
    FDTD/extensions/engine_ext_upml.cpp \
    FDTD/extensions/operator_extension.cpp \
    FDTD/extensions/engine_ext_mur_abc.cpp \
    FDTD/extensions/operator_ext_mur_abc.cpp \
    FDTD/extensions/operator_ext_cylinder.cpp \
    FDTD/extensions/engine_ext_cylinder.cpp \
    FDTD/extensions/operator_ext_excitation.cpp \
    FDTD/extensions/engine_ext_excitation.cpp

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
    FDTD/extensions/operator_ext_pml_sf.h \
    FDTD/extensions/engine_ext_pml_sf.h \
    FDTD/extensions/engine_ext_cylindermultigrid.h \
    FDTD/extensions/operator_ext_upml.h \
    FDTD/extensions/engine_ext_upml.h \
    FDTD/extensions/operator_ext_excitation.h \
    FDTD/extensions/engine_ext_excitation.h

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

# add git revision
GITREV = $$system(git describe --tags)
DEFINES += GIT_VERSION=\\\"$$GITREV\\\"





#
# create tar file (for the whole openEMS project)
#
tarball.target = tarball
tarball.commands = git archive --format=tar --prefix=openEMS-$$VERSION/ HEAD | bzip2 > openEMS-$${VERSION}.tar.bz2
QMAKE_EXTRA_TARGETS += tarball


#
# INSTALL (only the openEMS executable and matlab scripts)
#
install.target = install
install.commands = mkdir -p \"$(INSTALL_ROOT)/usr/bin\"
install.commands += && mkdir -p \"$(INSTALL_ROOT)/usr/share/openEMS/matlab\"
install.commands += && cp -at \"$(INSTALL_ROOT)/usr/bin/\" openEMS.sh openEMS_MPI.sh openEMS
install.commands += && cp -at \"$(INSTALL_ROOT)/usr/share/openEMS/\" matlab/
QMAKE_EXTRA_TARGETS += install


#
# create .PHONY target
#
phony.target = .PHONY
phony.depends = $$QMAKE_EXTRA_TARGETS
QMAKE_EXTRA_TARGETS += phony
