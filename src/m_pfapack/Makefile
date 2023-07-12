# Load settings

include ../../make.inc

# Variables

PFAPACK_SRC_DIR := .
PFAPACK_BUILD_DIR := $(BUILD_DIR)/pfapack
TARGET := $(BUILD_DIR)/pfapack.a

# File Lists
SRCS := $(wildcard $(PFAPACK_SRC_DIR)/*.f)
OBJS := $(SRCS:%=$(PFAPACK_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(PFAPACK_BUILD_DIR)

# Rule for making the target
$(TARGET): $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)

# Generic rule for compiling CPP files
$(PFAPACK_BUILD_DIR)/%.f.o: %.f
	$(FC) $(FCFLAGS) -c $< -o $@

# Cleaning everything
clean:
	rm -rf $(PFAPACK_BUILD_DIR)

.PHONY: all dirs clean

