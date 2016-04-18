CC = gcc
CFLAGS = -std=gnu99 -Wall -pedantic

OS = $(shell "uname")

ifeq ($(OS), Darwin)
	LIBMACTIME = libmactime.a
	LFLAGS= -lpthread -lm -lmactime
else
	LIBMACTIME = 
	LFLAGS= -lpthread -lm -lrt
endif

.PHONY: all
all: blurc thresc

.PHONY: clean
clean:
	-$(RM) *.o blurc thresc *.a

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o $(LIBMACTIME)
	$(CC) -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o -L. $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o $(LIBMACTIME)
	$(CC) -fast -o $@ thresmain.o ppmio.o thresfilter.o -L. $(LFLAGS)

libmactime.a: mactime.o
	ar cr libmactime.a mactime.o

arc:
	tar cf - *.c *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
