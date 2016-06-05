CC = gcc-4.9
CFLAGS = -fopenmp -std=c99 -Wall -pedantic

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
BIN = laplsolv

$(BIN_DIR)/$(BIN): $(OBJS)
	mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
