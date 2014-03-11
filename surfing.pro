QT -= core gui
QMAKE_CXXFLAGS += -std=gnu++0x
TARGET = surfing
TEMPLATE = app

HEADERS += \
    include/*.h \
    ../core/optionManager.h \
    include/cloudVolume.h \
    include/timer.h

SOURCES += \
    src/*.cpp

INCLUDEPATH += include/ ../core

CONFIG(debug, debug|release) {
    DESTDIR = .
    OBJECTS_DIR = build
} else {
    DESTDIR = .
    OBJECTS_DIR = build
}

