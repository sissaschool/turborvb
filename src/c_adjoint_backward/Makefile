# Load settings

include ../../make.inc

# Variables

SRC_DIR := .
ADJOINT_BACKWARD_BUILD_DIR := $(BUILD_DIR)/adjoint_backward
MODULE_DIR := $(BUILD_DIR)/modules
TARGET := $(BUILD_DIR)/adjoint_backward.a

# File Lists
#
# First set modules:

# Set rest files
SRCS := $(wildcard $(SRC_DIR)/*.f90)
OBJS := $(SRCS:%=$(ADJOINT_BACKWARD_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(ADJOINT_BACKWARD_BUILD_DIR)
	mkdir -p $(MODULE_DIR)
	echo $(MODULE_OBJ)

# Rule for making the target
$(TARGET): $(TARGET_MODULES) $(OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) $(TARGET_MODULES)

# Generic rule for compiling Fortran files
$(ADJOINT_BACKWARD_BUILD_DIR)/%.f90.o: %.f90 

	$(if $(filter $<, $(SOURCES_AGGRESSIVE)), \
	@echo "\e[31m Compiling $< with aggressive optimization\e[0m", \
	@echo "\e[32m Compiling $<\e[0m" \
	)

	$(eval FCFLAGS_SPECIAL=$(if $(filter $<, $(SOURCES_AGGRESSIVE)), \
	$(FCFLAGS_AGGRESSIVE), \
	$(FCFLAGS_PASSIVE)) \
	)

	$(FC) $(FCFLAGS) $(FCFLAGS_SPECIAL) \
	        $(MODULE_STORE) $(MODULE_DIR) \
	        $(MODULE_INCLUDE) $(MODULE_DIR) \
	        -c $< -o $@

# Generic rule for compiling C files
$(ADJOINT_BACKWARD_BUILD_DIR)/%.c.o: %.c
	@echo "\e[30m Compiling $<\e[0m"
	$(CC) $(CFLAGS) -c $< -o $@

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

