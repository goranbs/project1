TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp

LIBS += -larmadillo -lblas -llapack

HEADERS += \
    lib.h

