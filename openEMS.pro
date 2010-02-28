# -------------------------------------------------
# Project created by QtCreator 2010-02-26T22:34:51
# -------------------------------------------------
QT -= gui
TARGET = openEMS
CONFIG += console
CONFIG -= app_bundle
TEMPLATE = app
OBJECTS_DIR = obj
INCLUDEPATH += ../CSXCAD
LIBS += -L../CSXCAD \
    -lCSXCAD
SOURCES += main.cpp \
    FDTD/cartoperator.cpp \
    tools/ErrorMsg.cpp \
    tools/AdrOp.cpp
HEADERS += FDTD/cartoperator.h \
    tools/ErrorMsg.h \
    tools/AdrOp.h \
    tools/constants.h
