CXX = mpicxx

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))
BIN = $(BIN_DIR)/gas

$(BIN): $(OBJS) | $(BIN_DIR)
	$(CXX) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) -c -o $@ $<

$(OBJS): | $(OBJ_DIR)
$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

.PHONY: clean
clean:
	rm -rf $(BIN_DIR)
	rm -rf $(OBJ_DIR)
