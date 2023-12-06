# Load settings

include ../../make.inc
include ../../make.txt

# Variables

SRC_DIR := .
COMMON_BUILD_DIR := $(BUILD_DIR)/common
COMMON_MODULES_BUILD_DIR := $(BUILD_DIR)/common
MODULE_DIR := $(BUILD_DIR)/modules
TARGET := $(BUILD_DIR)/common.a
TARGET_MODULES := $(BUILD_DIR)/common_module.a

# File Lists
#
# First set modules:

MODULE_SRCS := $(SRC_DIR)/constants.f90    \
	       $(SRC_DIR)/symmetries.f90   \
	       $(SRC_DIR)/cell.f90         \
	       $(SRC_DIR)/dielectric.f90   \
	       $(SRC_DIR)/ewald.f90        \
	       $(SRC_DIR)/types.f90        \
	       $(SRC_DIR)/kpoints.f90      \
	       $(SRC_DIR)/mod_IO.f90       \
	       $(SRC_DIR)/mod_extpot.f90   \
	       $(SRC_DIR)/sub_comm.f90     \
	       $(SRC_DIR)/mpiio.f90        \
	       $(SRC_DIR)/kind.f90         \
	       $(SRC_DIR)/allio.f90        \

MODULE_OBJS := $(MODULE_SRCS:%=$(COMMON_MODULES_BUILD_DIR)/%.o)
MODULE_MODS := $(MODULE_SRCS:%=$(BUILD_DIR)/%.mod)

# Set rest files
SRCS := $(wildcard $(SRC_DIR)/*.f90)
SRCS += $(wildcard $(SRC_DIR)/*.c)
SRCS := $(filter-out $(MODULE_SRCS), $(SRCS))
OBJS := $(SRCS:%=$(COMMON_BUILD_DIR)/%.o)

# Default rule
all: dirs $(TARGET)

# Rule for making the build directory
dirs:
	mkdir -p $(COMMON_BUILD_DIR)
	mkdir -p $(COMMON_MODULES_BUILD_DIR)
	mkdir -p $(MODULE_DIR)

# Common module
$(TARGET_MODULES): $(MODULE_OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(AR) $(ARFLAGS) $(TARGET_MODULES) $(MODULE_OBJS)

# Rule for making the target
$(TARGET): $(TARGET_MODULES) $(OBJS)
	@echo "\e[34m Linking $(notdir $(TARGET))\e[0m"
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) $(TARGET_MODULES)

################################################################################

# Generic rule for compiling Fortran files
$(COMMON_BUILD_DIR)/%.f90.o: %.f90

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
$(COMMON_BUILD_DIR)/%.c.o: %.c
	@echo "\e[32m Compiling $<\e[0m"
	$(CC) $(CFLAGS) \
	-c $< -o $@

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

