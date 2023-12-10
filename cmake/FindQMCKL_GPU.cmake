#===========================================

# Try to find the QMCkl GPU library;
# If found, it will define the following variables (note the plural form):
#  QMCKL_GPU_FOUND - System has libqmckl;
#  QMCKL_GPU_INCLUDE_DIRS - The QMCKL GPU include directories;
#  QMCKL_GPU_LIBRARIES - The libraries needed to use QMCKL GPU;

# If QMCKL GPU has been installed in a non-standard location, one can set an
# environment variable $QMCKL_GPU_DIR in the current shell:
# $ export QMCKL_GPU_DIR=<custom_path>
# to indicate the prefix used during the QMCKL installation
# (typically `./configure prefix=<custom_path> ..` or `cmake -DCMAKE_INSTALL_DIR=<custom_path> ..`)

# This file should be located WITHIN your project source tree.
# (e.g. in cmake/FindQMCKL_GPU.cmake)
# How to use it in your project CMakeLists.txt:

# This is needed to locate FindQMCKL_GPU.cmake file, modify it according to your source tree.
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# find_package(QMCKL_GPU)
# if (QMCKL_GPU_FOUND)
#   include_directories(${QMCKL_GPU_INCLUDE_DIRS})
#   target_link_libraries(your_target ${QMCKL_GPU_LIBRARIES})
# endif()

#===========================================

# This file is distirbuted under the BSD 3-Clause License.
# Copyright (c) 2021, TREX Center of Excellence

#===========================================

message(" ")
message("Looking for the QMCKL_GPU library:")

set(QMCKL_GPU_SEARCH_PATHS
        ~/Library/Frameworks
        /Library/Frameworks
        /usr/local
        /usr
        /sw # Fink
        /opt/local # DarwinPorts
        /opt/csw # Blastwave
        /opt
)

find_path(QMCKL_GPU_INCLUDE_DIR
          NAMES qmckl_gpu_f.F90
          HINTS $ENV{QMCKL_GPU_DIR}
          PATH_SUFFIXES include
          PATHS ${QMCKL_GPU_SEARCH_PATHS}
          )


# No need to specify platform-specific prefix (e.g. libqmckl_gpu on Unix) or
# suffix (e.g. .so on Unix or .dylib on MacOS) in NAMES. CMake takes care of that.
find_library(QMCKL_GPU_LIBRARY
             NAMES qmckl_gpu
             HINTS $ENV{QMCKL_GPU_DIR}
             PATH_SUFFIXES lib64 lib
             PATHS ${QMCKL_GPU_SEARCH_PATHS}
             )

get_filename_component(QMCKL_GPU_LIBRARY_DIR ${QMCKL_GPU_LIBRARY} DIRECTORY)

set(QMCKL_GPU_LIBRARIES "-L${QMCKL_GPU_LIBRARY_DIR} -lqmckl_gpu")

# Handle the QUIETLY and REQUIRED arguments and set QMCKL_GPU_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QMCKL_GPU DEFAULT_GPU_MSG QMCKL_GPU_LIBRARY QMCKL_GPU_INCLUDE_DIR )
MARK_AS_ADVANCED(QMCKL_GPU_INCLUDE_DIR QMCKL_GPU_LIBRARY)

# Mot setting _INCLUDE_DIR and _LIBRARIES is considered a bug,
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Find-Libraries
set(QMCKL_GPU_INCLUDE_DIRS ${QMCKL_GPU_INCLUDE_DIR})

