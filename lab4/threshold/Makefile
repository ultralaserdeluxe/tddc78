CC = mpicc
CFLAGS = -std=gnu99 -Wall -pedantic -I$(VT_ROOT)/include
LDFLAGS = -L$(VT_LIB_DIR)
LDLIBS =  -lVT $(VT_ADD_LIBS)

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
BIN = $(BIN_DIR)/threshold

$(BIN): $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS)  -o $@ $^ $(LDFLAGS) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJS): | $(OBJ_DIR)
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
