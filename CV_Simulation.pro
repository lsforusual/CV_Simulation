#-------------------------------------------------
#
# Project created by QtCreator 2016-09-17T21:09:08
#
#-------------------------------------------------

QT       += core
QT       -= gui
CONFIG   += console
CONFIG   += app_buddle

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = CV_Simulation
TEMPLATE = app


SOURCES += main.cpp\
    datarw.cpp \
    fitfunc.cpp \
    iodata.cpp \
    datascan.cpp

HEADERS  += \
    inputdata.h \
    outputdata.h

FORMS    +=
