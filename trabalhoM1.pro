TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        ../FDA_Original_C_unsigned/diff2d.c \
        diff2d.c \
        main.c \
        pgmfiles.c \
        pgmtolist.c

DISTFILES += \
    lena_noise.pgm \
    lena_noise_out.pgm

HEADERS += \
    diff2d.h \
    pgmfiles.h
