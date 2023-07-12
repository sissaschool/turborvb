# Load settings

include ../../make.inc

# Variables

SRC_DIR := .
COMMON_BUILD_DIR := $(BUILD_DIR)/common
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

MODULE_OBJS := $(MODULE_SRCS:%=$(COMMON_BUILD_DIR)/%.o)
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
	mkdir -p $(MODULE_DIR)
	echo $(MODULE_OBJ)

# Compile modules
$(COMMON_BUILD_DIR)/%.f90.o: %.f90
	$(FC) $(FCFLAGS) \
       	-I $(MODULE_DIR) \
       	-J $(MODULE_DIR) \
       	-c $< -o $@

# Common module
$(TARGET_MODULES): $(MODULE_OBJS)
	$(AR) $(ARFLAGS) $(TARGET_MODULES) $(MODULE_OBJS)

# Rule for making the target
$(TARGET): $(TARGET_MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJS) $(TARGET_MODULES)

################################################################################

# Generic rule for compiling Fortran files
$(COMMON_BUILD_DIR)/%.f90.o: %.f90 
	$(FC) $(FCFLAGS) -I $(MODULE_DIR) -c $< -o $@

# Generic rule for compiling C files
$(COMMON_BUILD_DIR)/%.c.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Cleaning everything
clean:
	rm -rf $(BUILD_DIR)

.PHONY: all dirs clean

