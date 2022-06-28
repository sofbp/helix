TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -O2 -march=native -std=c++11 -Wno-unused-parameter -Wno-sign-compare

SOURCES += \
        main.cpp

HEADERS += \
    residue.h \
    peptide.h \
    energy_calculator.h \
    atom.h
