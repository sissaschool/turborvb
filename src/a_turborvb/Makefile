# Load settings

include ../../make.inc

# Variables

SRC_DIR := .
MAIN_BUILD_DIR := $(BUILD_DIR)/main
MODULE_DIR := $(BUILD_DIR)/modules
TARGET := $(BUILD_DIR)/turborvb.x
TARGET_NAME := $(BUILD_DIR)/turborvb$(SUFFIX).x
TARGET_MODULES := $(BUILD_DIR)/turborvb$(SUFFIX)_module.a

# File Lists
#
# First set modules:

MODULE_SRCS := $(SRC_DIR)/scal_lins.f90

MODULE_OBJS := $(MODULE_SRCS:%=$(MAIN_BUILD_DIR)/%.o)
MODULE_MODS := $(MODULE_SRCS:%=$(BUILD_DIR)/%.mod)

# Set rest files
SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRCS += $(wildcard $(SRC_DIR)/*.c)
SRCS := $(filter-out $(MODULE_SRCS), $(SRCS))
OBJS := $(SRCS:%=$(MAIN_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(MAIN_BUILD_DIR)
	mkdir -p $(MODULE_DIR)
	echo $(MODULE_OBJS)

# Compile modules
$(MAIN_BUILD_DIR)/%.f90.o: %.f90
	$(FC) $(FCFLAGS) \
		$(MODULE_STORE) $(MODULE_DIR) \
 	        $(MODULE_INCLUDE) $(MODULE_DIR) \
		-c $< -o $@


# Common module
$(TARGET_MODULES): $(MODULE_OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(AR) $(ARFLAGS) $(TARGET_MODULES) $(MODULE_OBJS)

# Rule for making the target
$(TARGET): $(TARGET_MODULES) $(OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(FC) $(FLINK) \
	$(OBJS) \
	$(MODULE_OBJS) \
	$(BUILD_DIR)/complex.a \
	$(BUILD_DIR)/complex_module.a \
	$(BUILD_DIR)/adjoint_forward.a \
	$(BUILD_DIR)/adjoint_forward_module.a \
	$(BUILD_DIR)/adjoint_backward.a \
	$(BUILD_DIR)/qlapack.a \
	$(BUILD_DIR)/common.a \
	$(BUILD_DIR)/common_module.a \
	$(BUILD_DIR)/pfapack.a \
	$(BUILD_DIR)/complex.a \
	$(BUILD_DIR)/complex_module.a \
	$(BUILD_DIR)/adjoint_forward.a \
	$(BUILD_DIR)/adjoint_forward_module.a \
	$(BUILD_DIR)/adjoint_backward.a \
	$(BUILD_DIR)/qlapack.a \
	$(BUILD_DIR)/common.a \
	$(BUILD_DIR)/common_module.a \
	$(BUILD_DIR)/pfapack.a \
	-J $(MODULE_DIR) \
	$(LINK_LIBS) \
       	-o $(TARGET)

################################################################################

# Generic rule for compiling Fortran files
$(MAIN_BUILD_DIR)/%.f90.o: %.f90 

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
$(MAIN_BUILD_DIR)/%.c.o: %.c
	@echo "\e[32m Compiling $<\e[0m"
	$(CC) $(CFLAGS) -c $< -o $(SUFFIX)

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

