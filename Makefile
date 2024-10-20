# Makefile for ENGDRAW_COPY_2 project

# Compiler and linker
CXX = g++
CC = gcc

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build
DEPENDENCIES_DIR = dependencies

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
OBJS += $(BUILD_DIR)/glad.o

# Executable
TARGET = $(BUILD_DIR)/app

# Include paths
INCLUDE_PATHS = -I$(INCLUDE_DIR) \
                -I$(DEPENDENCIES_DIR)/include \
                -I/opt/homebrew/opt/glfw/include \
                -I/opt/homebrew/opt/glew/include \
                -I/opt/homebrew/opt/glm/include

# Library paths
LIBRARY_PATHS = -L$(DEPENDENCIES_DIR)/library \
                -L/opt/homebrew/opt/glfw/lib \
                -L/opt/homebrew/opt/glew/lib \
                -L/opt/homebrew/opt/glm/lib

# Libraries to link
LIBS = -lglfw -lGLEW -framework OpenGL -framework Cocoa -framework IOKit -framework CoreVideo

# Compiler flags for C++
CXXFLAGS = -std=c++17 -Wall -g -DGLM_ENABLE_EXPERIMENTAL -Wno-deprecated-declarations $(INCLUDE_PATHS)

# Compiler flags for C
CFLAGS = -Wall -g -Wno-deprecated-declarations $(INCLUDE_PATHS)

# Build rules
all: $(TARGET)

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Compile .cpp files to .o files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile glad.c separately using C compiler and CFLAGS
$(BUILD_DIR)/glad.o: $(SRC_DIR)/glad.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Link the object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBRARY_PATHS) $^ -o $@ $(LIBS)

# Clean up build artifacts
clean:
	rm -rf $(BUILD_DIR)/*.o $(TARGET)

.PHONY: all clean