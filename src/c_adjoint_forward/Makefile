# Load settings

include ../../make.inc

# Variables

SRC_DIR := .
ADJOINT_FORWARD_BUILD_DIR := $(BUILD_DIR)/adjoint_forward
MODULE_DIR := $(BUILD_DIR)/modules
TARGET := $(BUILD_DIR)/adjoint_forward.a
TARGET_MODULES := $(BUILD_DIR)/adjoint_forward_module.a

# File Lists
#
# First set modules:

MODULE_SRCS := $(SRC_DIR)/dspev_drv.f90

MODULE_OBJS := $(MODULE_SRCS:%=$(ADJOINT_FORWARD_BUILD_DIR)/%.o)
MODULE_MODS := $(MODULE_SRCS:%=$(BUILD_DIR)/%.mod)

# Set rest files
SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRCS += $(wildcard $(SRC_DIR)/*.c)
SRCS := $(filter-out $(MODULE_SRCS), $(SRCS))
OBJS := $(SRCS:%=$(ADJOINT_FORWARD_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(ADJOINT_FORWARD_BUILD_DIR)
	mkdir -p $(MODULE_DIR)
	echo $(MODULE_OBJ)

# Compile modules
$(ADJOINT_FORWARD_BUILD_DIR)/%.f90.o: %.f90
	$(FC) $(FCFLAGS) -J $(MODULE_DIR) -I $(MODULE_DIR) -c $< -o $@

# Common module
$(TARGET_MODULES): $(MODULE_OBJS)
	$(AR) $(ARFLAGS) $(TARGET_MODULES) $(MODULE_OBJS)

# Rule for making the target
$(TARGET): $(TARGET_MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) $(TARGET_MODULES)

################################################################################

# Generic rule for compiling Fortran files
$(ADJOINT_FORWARD_BUILD_DIR)/%.f90.o: %.f90 
	$(FC) $(FCFLAGS) -J$(MODULE_DIR) -c $< -o $@

# Generic rule for compiling C files
$(ADJOINT_FORWARD_BUILD_DIR)/%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

