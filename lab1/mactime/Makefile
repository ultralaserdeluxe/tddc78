CC = gcc
CFLAGS = -std=gnu99 -Wall -pedantic

SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj

LIB_DIR = lib
LIB_NAME = mactime
LIB_FULL_NAME = $(LIB_DIR)/lib$(LIB_NAME).a

.PHONY: all
all: $(LIB_FULL_NAME)

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) $(LIB_DIR)

$(LIB_FULL_NAME): $(OBJ_DIR)/$(LIB_NAME).o
	mkdir -p $(LIB_DIR)
	ar cr $(LIB_FULL_NAME) $(OBJ_DIR)/$(LIB_NAME).o

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(INC_DIR)/%.h
	mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -I $(INC_DIR) -c -o $@ $<
