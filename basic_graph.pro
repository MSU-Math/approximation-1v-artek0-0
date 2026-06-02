CONFIG        += c++11
QMAKE_CXXFLAGS += -Werror
HEADERS       = window.h \
                newt.h \
                bessel.h
SOURCES       = main.cpp \
                window.cpp \
                newt.cpp \
                bessel.cpp
QT += widgets
