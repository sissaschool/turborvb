# Load settings

include ../../make.inc

# Variables

SRC_DIR := .
QLAPACK_BUILD_DIR := $(BUILD_DIR)/qlapack
MODULE_DIR := $(BUILD_DIR)/modules
TARGET := $(BUILD_DIR)/qlapack.a

# Set rest files
SRCS := $(wildcard $(SRC_DIR)/*.f)
OBJS := $(SRCS:%=$(QLAPACK_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(QLAPACK_BUILD_DIR)
	echo $(OBJS)

# Generic rule for compiling Fortran files
$(QLAPACK_BUILD_DIR)/%.f.o: %.f 
	$(FC) $(FCFLAGS) -J$(MODULE_DIR) -c $< -o $@

# Rule for making the target
$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

