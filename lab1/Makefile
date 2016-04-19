CC = gcc
CFLAGS = -std=gnu99 -Wall -pedantic

OS = $(shell "uname")

ifeq ($(OS), Darwin)
	MACTIME = mactime
	LFLAGS= -lpthread -lm -lmactime
else
	MACTIME = 
	LFLAGS= -lpthread -lm -lrt
endif

.PHONY: all
all: blur

.PHONY: clean
clean:
	cd blur && $(MAKE) clean
	cd mactime && $(MAKE) clean

.PHONY: blur
blur: $(MACTIME)
	cd blur	&& $(MAKE)

.PHONY: mactime
mactime: 
	cd mactime && $(MAKE)
