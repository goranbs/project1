TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    lib.cpp \
    hfsolve.cpp

LIBS += -larmadillo -lblas -llapack

HEADERS += \
    lib.h \
    hfsolve.h \
    testingHFSolve.h

CXX_FLAGS += -O3 -std=c++0x
