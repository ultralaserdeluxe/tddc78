CC = gcc
CFLAGS = -std=gnu99 -Wall -pedantic

SRC_DIR = src
OBJ_DIR = obj

BIN_DIR = bin
BIN_NAME = threshold
BIN_FULL_NAME = $(BIN_DIR)/$(BIN_NAME)

.PHONY: all
all: $(BIN_FULL_NAME)

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

$(BIN_FULL_NAME): $(OBJ_DIR)/$(BIN_NAME).o $(OBJ_DIR)/ppm.o
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ -l mpi

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<
