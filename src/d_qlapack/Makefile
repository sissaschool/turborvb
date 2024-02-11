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

# Rule for making the target
$(TARGET): $(OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS)

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

