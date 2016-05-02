CC = gcc
CFLAGS = -std=gnu99 -Wall -pedantic

OS = $(shell "uname")

.PHONY: all
all: blur

.PHONY: clean
clean:
	cd blur && $(MAKE) clean

.PHONY: blur
blur:
	cd blur	&& $(MAKE)

