QT -= core gui
QMAKE_CXXFLAGS += -std=gnu++0x
TARGET = surfing
TEMPLATE = app

HEADERS += \
    include/*.h \
    ../core/optionManager.h

SOURCES += \
    src/*.cpp \
#    src/*.cu

#SOURCES -= \
#    src/*.cu

INCLUDEPATH += include/ ../core

#CUDA_SOURCES += src/computeIsoValues_CUDA.cu
#CUDA_DIR = /usr/local/cuda-5.5
#INCLUDEPATH += $$CUDA_DIR/include
#LIBS += -L$$CUDA_DIR/lib -lcudart -lcuda
#CUDA_ARCH = sm_20
#NVCCFLAGS = --compiler-options -fno-strict-aliasing -use_fast_math --ptxas-options=-v -std=gnu++0x

#CUDA_INC = $$join(INCLUDEPATH,' -I','-I',' ')

#cuda.commands = $$CUDA_DIR/bin/nvcc -m64 -O3 -arch=$$CUDA_ARCH -c $$NVCCFLAGS \
#                $$CUDA_INC $$LIBS  ${QMAKE_FILE_NAME} -o ${QMAKE_FILE_OUT} \
#                2>&1 | sed -r \"s/\\(([0-9]+)\\)/:\\1/g\" 1>&2
#cuda.dependency_type = TYPE_C
#cuda.depend_command = $$CUDA_DIR/bin/nvcc -O3 -M $$CUDA_INC $$NVCCFLAGS   ${QMAKE_FILE_NAME}
#cuda.input = CUDA_SOURCES
#cuda.output = ${OBJECTS_DIR}${QMAKE_FILE_BASE}.o
#QMAKE_EXTRA_COMPILERS += cuda

CONFIG(debug, debug|release) {
    DESTDIR = .
    OBJECTS_DIR = build
} else {
    DESTDIR = .
    OBJECTS_DIR = build
}

