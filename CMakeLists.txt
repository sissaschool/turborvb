
######################################################################
#
# Copyright (C) 2022 TurboRVB group
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
######################################################################

######################################################################
#
# Here are some definitions of variables used in this CMakeLists.txt:
#
# EXECUTABLES_S_L   - list of all serial executables
# LIBRARIES_S_L     - list of all serial libraries
# EXECUTABLES_P_L   - list of all parallel executables
# LIBRARIES_P_L     - list of all parallel libraries
# AGRESSIVE_L       - list of aggressively optimizable serial sources
# LIBRARIES_L       - list of all libraries
# PPDIRECTIVES      - stores list of preprocessor flags common
#                     to all targets
# PARALLEL_F        - list of parallel precompiler flags
# SERIAL_F          - list of serial precompiler flags
# AGRESIVE_FINAL_L  - final list of agresive optimized object files
# PASIVE_FINAL_L    - final list of non-agresive optimized object files
# AGRESIVE_FLAG_S   - flag for agressive optimizations
# PASIVE_FLAG_S     - flag for non-agressive optimizations
# NATIVE_OMP_B      - flag if true CMake native OMP linking is done
#
######################################################################

# Geting the git revision:

find_package(Git QUIET)

if(GIT_FOUND)
  execute_process(
      COMMAND
          ${GIT_EXECUTABLE} rev-parse --short HEAD
      RESULT_VARIABLE
          SHORT_HASH_RESULT
      OUTPUT_VARIABLE
          SHORT_HASH
      )
  if ( SHORT_HASH STREQUAL "" )
    set( SHORT_HASH "N/A" )
  endif()
  string(STRIP ${SHORT_HASH} SHORT_HASH)
else()
  set( SHORT_HASH "N/A" )
endif()

message("")
message("-------------------------------------------------------------------------")
message("   TurboRVB CMake build (git rev. ${SHORT_HASH})")
message("-------------------------------------------------------------------------")
message("")

cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

project(TurboRVB LANGUAGES Fortran C) 

cmake_policy(SET CMP0057 NEW) # Adding IN_LIST
cmake_policy(SET CMP0054 NEW) # Adding IN_LIST

######################################################################
#
# Define special target properties
#

define_property(TARGET 
                PROPERTY NOFLAG 
                BRIEF_DOCS "List of precompiler flags \
                            that has to be removed from target"
                FULL_DOCS  "List of precompiler flags \
                            that has to be removed from target.\
                            It works for serial and parallel flags." )
                                                   
#
######################################################################

######################################################################
#
# Set PARALLEL_F
#

list( APPEND PARALLEL_F PARALLEL )

#
######################################################################

######################################################################
#
# Here we load modules, there is one important for ScaLAPACK.
# It was obtained from here:
#
# https://github.com/ecrc/ecrc_cmake/blob/master/modules/find/FindSCALAPACK.cmake
#
# Is published under French CeCILL-C license.
#

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#
######################################################################

######################################################################
#
# I have a suggestion to nomenclature where all externally set
# variables will have EXT_ prefix. This might not be a good strategy
# but for the time being ...
#

option(EXT_SERIAL     "Compile serial version"                ON)
option(EXT_PARALLEL   "Compile parallel version"              ON)
option(EXT_OPT        "Turn on optimization"                  ON)
option(EXT_DEBUG      "Turn on debug (even in release)"      OFF)
option(EXT_TIME       "Turn on internal time counter"        OFF)
option(EXT_DFT        "Compile DFT code"                      ON)
option(EXT_QMC        "Compile QMC code"                      ON)
option(EXT_MODTEST    "Compile module test tools"             ON)
option(EXT_TOOLS      "Compile auxiliary tools"               ON)
option(EXT_GPU        "Compile GPU version"                   ON)
option(EXT_LTO        "Enable LTO (if optimizations are on)" OFF)
option(EXT_SPEEDTEST  "Enable speed tests targets"           OFF)
option(EXT_STATICPACK "Produce static packaging if possible" OFF)
option(EXT_DETAIL     "Produce more detailed config output"  OFF)
option(EXT_NVTX       "Turn on Nvidia NVTX ranges"           OFF)

set(EXT_BLAS_LIB   ""    CACHE STRING    "Blas libraries")
set(EXT_LAPACK_LIB ""    CACHE STRING    "Lapack libraries")
set(EXT_OTHER_LIB  ""    CACHE STRING    "Other libraries to link")

set(EXT_FLAGS      ""    CACHE STRING    "Common preprocessor flags")
set(EXT_BLKL_FLAGS ""    CACHE STRING    "Blacklist common preprocessor flags")

set(EXT_AGRS       ""    CACHE STRING    "Add to agressive")
set(EXT_BLKL_AGRS  ""    CACHE STRING    "Remove from agressive")

set(EXT_GPUTYPE    ""    CACHE STRING    "GPU type (like comp. cap.)")

#
######################################################################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)

    set(CMAKE_C_FLAGS_RELEASE "")
    set(CMAKE_Fortran_FLAGS_RELEASE "")
endif()

######################################################################
#
# CPack part
#

if(EXT_STATIC_PACK)
  set(CMAKE_INSTALL_MFC_LIBRARIES ON)
  set(CMAKE_INSTALL_OPENMP_LIBRARIES ON)
  include(InstallRequiredSystemLibraries)
endif()

set(CPACK_PACKAGE_VENDOR "prof. Sorella's group SISSA Trieste")
set(CPACK_PACKAGE_DESCRIPTION "TurboRVB is Quantum Monte Carlo code")
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/License)
include(CPack)

#
######################################################################

# Here goes all of the includes

include(CheckFortranCompilerFlag)
include(CheckCCompilerFlag)

# Add subdirectory with helper macros and functions:

add_subdirectory(cmake)

list( APPEND PPDIRECTIVES _SIMD )

######################################################################
#
# This variable will make BLAS/LAPACK statically linked if set ON
# However, there is a problem with statically linked libraries and
# MKL. Basically CMake has problem correctly solve cyclic dependecies
# of MKL. There are some tricks with repeated linking, but they dis
# not work in my case.
#

set( BLA_STATIC OFF )

#
######################################################################

if(EXT_PARALLEL)
  find_package(MPI COMPONENTS Fortran C)

  if(NOT MPI_FOUND)
    message( WARNING "MPI not found, turning off compilation of parallel code." )
    set(EXT_PARALLEL OFF)
  endif()
endif()

list( APPEND MKL_ALLOWED_COMPILERS "GNU" )
list( APPEND MKL_ALLOWED_COMPILERS "Intel" )
list( APPEND MKL_ALLOWED_COMPILERS "IntelLLVM" )
list( APPEND MKL_ALLOWED_COMPILERS "NVHPC" )

if( CMAKE_Fortran_COMPILER_ID IN_LIST MKL_ALLOWED_COMPILERS)

  set(MKL_LINK "static")
  set(MKL_INTERFACE lp64)

  if(EXT_PARALLEL)
    set(MKL_MPI openmpi)
    if( CMAKE_Fortran_COMPILER_ID MATCHES "Intel*.")
      set(MKL_MPI intelmpi)
    endif()
    set( ENABLE_SCALAPACK ON )
  endif()

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(MKL_THREADING pgi_thread)
  endif()

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MKL_THREADING gnu_thread)
  endif()
  
  if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel.*")
    set(MKL_THREADING intel_thread)
  endif()
  
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
      set( CMAKE_Fortran_COMPILER_ID "PGI" )
  endif()
  
  find_package(MKL QUIET)

  if (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
      set( CMAKE_Fortran_COMPILER_ID "NVHPC" )
  endif()

endif()

######################################################################
#
# The main idea behind this to first check the if math libs are set,
# if not then try to find them and if it is not possibile to find
# them throw error
#

if(EXT_BLAS_LIB)
  link_libraries(${EXT_BLAS_LIB})
else()
  if(MKL_FOUND)
    link_libraries(MKL::MKL)
  else()
    find_package(BLAS REQUIRED)
    link_libraries(${BLAS_LIBRARIES})
  endif()
endif()

if(EXT_LAPACK_LIB)
  link_libraries(${EXT_LAPACK_LIB})
else()
  if(MKL_FOUND)
    link_libraries(MKL::MKL)
  else()
    find_package(LAPACK)
    if(LAPACK_FOUND)
        link_libraries(${LAPACK_LIBRARIES})
    else()
      if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
        message( FATAL_ERROR "LAPACK not found!" )
      endif()
    endif()
  endif()
endif()

#
######################################################################

foreach(OTHER_LIB ${EXT_OTHER_LIB})
  link_libraries(${OTHER_LIB})
endforeach()

######################################################################
#
# Setting native OMP to off
#

set( NATIVE_OMP_B OFF )

#
######################################################################

find_package(OpenMP)
find_package(Threads)

if(EXT_GPU)
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message( WARNING "Due to a bug in gfortran compiler, it is not possible to compile GPU version.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    message( WARNING "Turbo uses CUDA libraries, which are not supported by Intel compilers.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message( WARNING "Turbo uses CUDA libraries, which are not supported by Intel compilers.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    message( WARNING "Due to a bug in gfortran compiler, it is not possible to compile GPU version.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
  if(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    message( WARNING "Turbo uses CUDA libraries, which are not supported by Intel compilers.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
  if(CMAKE_C_COMPILER_ID STREQUAL "Intel")
    message( WARNING "Turbo uses CUDA libraries, which are not supported by Intel compilers.\n"
                     "Turning off GPU support")
    set( EXT_GPU OFF )
  endif()
endif()

if(EXT_GPU)

  find_package(CUDAToolkit)
  
  if(CUDAToolkit_FOUND)
    list(APPEND PPDIRECTIVES _OFFLOAD)
    list(APPEND PPDIRECTIVES _CUBLAS)
    list(APPEND PPDIRECTIVES _CUSOLVER)
    #link_libraries(${CUDA_cusolver_LIBRARY})
    #link_libraries(${CUDA_cublas_LIBRARY})
    #link_libraries(${CUDA_LIBRARIES})
  else()
      message( WARNING "No GPU math libraries found code will be very inefficient." )
    list(APPEND PPDIRECTIVES _OFFLOAD)
    #set(EXT_GPU OFF)
  endif()

endif()

if(OpenMP_FOUND AND NATIVE_OMP_B)
    link_libraries(${OpenMP_Fortran_LIBRARIES})
    link_libraries(${OpenMP_C_LIBRARIES})
endif()

if(Threads_FOUND)
    link_libraries(${Threads_Fortran_LIBRARIES})
    link_libraries(${Threads_C_LIBRARIES})
endif()

# Setup libraries and executables:

if (EXT_QMC)
  if(EXT_SERIAL)
    add_target_lib_serial_wrapper(pfapack)
    add_target_lib_serial_wrapper(qlapack-serial)
    add_target_lib_serial_wrapper(common-serial)
    add_target_lib_serial_wrapper(adjoint-forward-serial)
    add_target_lib_serial_wrapper(adjoint-backward-serial)
    add_target_lib_serial_wrapper(complex-serial)
    add_target_exe_serial_wrapper(turborvb-serial)
  endif()
  
  if(EXT_PARALLEL)
    add_target_lib_serial_wrapper(pfapack)
    add_target_lib_parallel_wrapper(qlapack-mpi)
    add_target_lib_parallel_wrapper(common-mpi)
    add_target_lib_parallel_wrapper(adjoint-forward-mpi)
    add_target_lib_parallel_wrapper(adjoint-backward-mpi)
    add_target_lib_parallel_wrapper(complex-mpi)
    add_target_exe_parallel_wrapper(turborvb-mpi)
  endif()
endif()

if(EXT_TIME)
  list( APPEND PPDIRECTIVES _TIME )
endif()

if(EXT_NVTX)
  list( APPEND PPDIRECTIVES _NVTX )
endif()

if(EXT_DEBUG)
  list( APPEND PPDIRECTIVES _DEBUG )
endif()

if(EXT_DFT)
  if(EXT_SERIAL)
    add_target_lib_serial_wrapper(pfapack)
    add_target_lib_serial_wrapper(qlapack-serial)
    add_target_lib_serial_wrapper(common-serial)
    add_target_lib_serial_wrapper(adjoint-forward-serial)
    add_target_lib_serial_wrapper(adjoint-backward-serial)
    add_target_lib_serial_wrapper(complex-serial)
    add_target_exe_serial_wrapper(prep-serial)
  endif()

  if(EXT_PARALLEL)
    add_target_lib_serial_wrapper(pfapack)
    add_target_lib_parallel_wrapper(qlapack-mpi)
    add_target_lib_parallel_wrapper(common-mpi)
    add_target_lib_parallel_wrapper(adjoint-forward-mpi)
    add_target_lib_parallel_wrapper(adjoint-backward-mpi)
    add_target_lib_parallel_wrapper(complex-mpi)
    add_target_exe_parallel_wrapper(prep-mpi)
  endif()

endif()

######################################################################
#
# One can easilly setup prefix/name/suffix with these options:
#
# https://stackoverflow.com/questions/31038963/how-do-you-rename-a-library-filename-in-cmake
#
# Cheking for module can be done this way:
#
# https://stackoverflow.com/questions/59507199/how-to-check-if-a-module-included-with-include-by-cmake-exists
#
# This solution was abandoned, but is interesting so I will keep it here comented
#

if(EXT_PARALLEL)
  if(MPI_FOUND)
    foreach( TARGET IN LISTS LIBRARIES_P_L )
      target_include_directories( ${TARGET} PUBLIC ${MPI_INCLUDE_PATH} )
      target_link_libraries( ${TARGET} PUBLIC MPI::MPI_Fortran )
    endforeach()
  else()
    set(EXT_PARALLEL OFF)
    message( WARNING "MPI not found, turning off compilation of parallel code." )
  endif()

endif()

#
######################################################################

######################################################################
#
# https://stackoverflow.com/questions/26144327/how-can-i-suppress-the-warning-that-a-package-configuration-file-was-not-found
#
# This solution was abandoned
#

if(EXT_PARALLEL)
  if(MKL_FOUND)
    list(APPEND PARALLEL_F __SCALAPACK)
    message( STATUS "Linking MKL ScaLAPACK" )
    foreach( EXECUTABLE IN LISTS EXECUTABLES_P_L )
      target_link_libraries( ${EXECUTABLE} PRIVATE MKL::MKL )
    endforeach()
  else()
    find_package(SCALAPACK)
    if(SCALAPACK_FOUND)
      list(APPEND PARALLEL_F __SCALAPACK)
      message( STATUS "Linking ScaLAPACK: ${SCALAPACK_LIBRARIES}" )
      foreach( EXECUTABLE IN LISTS EXECUTABLES_P_L )
        target_link_libraries( ${EXECUTABLE} PRIVATE ${SCALAPACK_LIBRARIES} )
      endforeach()
    endif()
  endif()
endif()

#
######################################################################

if(EXT_TOOLS)

  # These tools have to be threaded in special way:
  
  # serial tools
  add_target_exe_serial_wrapper( makefort10 )
  add_target_exe_serial_wrapper( readforward-serial )
  add_target_exe_serial_wrapper( convertfort10-serial )
  add_target_exe_serial_wrapper( convertfort10mol-serial )
  
  # parallel tools
  if(EXT_PARALLEL)
      add_target_exe_parallel_wrapper( readforward-mpi )
      add_target_exe_parallel_wrapper( convertfort10-mpi )
      add_target_exe_parallel_wrapper( convertfort10mol-mpi )
  endif()

  # simple serial tools
  list(APPEND TOOLS_LIST_SIMPLE corrvar )
  list(APPEND TOOLS_LIST_SIMPLE corrforza )
  list(APPEND TOOLS_LIST_SIMPLE corrforzap )
  list(APPEND TOOLS_LIST_SIMPLE corrforzap_complex )
  list(APPEND TOOLS_LIST_SIMPLE readf )
  list(APPEND TOOLS_LIST_SIMPLE readffn )
  list(APPEND TOOLS_LIST_SIMPLE readalles )
  list(APPEND TOOLS_LIST_SIMPLE funvsa )
  list(APPEND TOOLS_LIST_SIMPLE copyjas )
  list(APPEND TOOLS_LIST_SIMPLE copydet )
  list(APPEND TOOLS_LIST_SIMPLE pseudo )
  list(APPEND TOOLS_LIST_SIMPLE plot_orbitals )
  list(APPEND TOOLS_LIST_SIMPLE assembling_pseudo)
  list(APPEND TOOLS_LIST_SIMPLE orthomol )

  add_target_lib_serial_wrapper(pfapack)

  if(EXT_GPU)
      #serial
      add_target_lib_serial_wrapper(common-nogpu-serial)
      add_target_lib_serial_wrapper(complex-nogpu-serial)
      add_target_lib_serial_wrapper(qlapack-nogpu-serial)
      add_target_lib_serial_wrapper(adjoint-forward-nogpu-serial)
      add_target_lib_serial_wrapper(adjoint-backward-nogpu-serial)
      #parallel
      if(EXT_PARALLEL)
          add_target_lib_parallel_wrapper(common-mpi)
          add_target_lib_parallel_wrapper(complex-mpi)
          add_target_lib_parallel_wrapper(qlapack-mpi)
          add_target_lib_parallel_wrapper(adjoint-forward-mpi)
          add_target_lib_parallel_wrapper(adjoint-backward-mpi)
      endif()
  else()
      #serial
      add_target_lib_serial_wrapper(common-serial)
      add_target_lib_serial_wrapper(complex-serial)
      add_target_lib_serial_wrapper(qlapack-serial)
      add_target_lib_serial_wrapper(adjoint-forward-serial)
      add_target_lib_serial_wrapper(adjoint-backward-serial)
      #parallel
      if(EXT_PARALLEL)
          add_target_lib_parallel_wrapper(common-mpi)
          add_target_lib_parallel_wrapper(complex-mpi)
          add_target_lib_parallel_wrapper(qlapack-mpi)
          add_target_lib_parallel_wrapper(adjoint-forward-mpi)
          add_target_lib_parallel_wrapper(adjoint-backward-mpi)
      endif()
  endif()

  foreach(TOOL_NAME ${TOOLS_LIST_SIMPLE})
      add_target_exe_serial_wrapper(${TOOL_NAME})
  endforeach()

endif()

if(EXT_MODTEST)

  add_target_exe_serial_wrapper(testad test_tools)
  add_target_exe_serial_wrapper(test_zdotc test_tools)
  add_target_exe_serial_wrapper(test_dger test_tools)
  add_target_exe_serial_wrapper(test_dger2 test_tools)
  add_target_exe_serial_wrapper(test_zgeru test_tools)
  add_target_exe_serial_wrapper(test_zger2 test_tools)
  add_target_exe_serial_wrapper(test_dgemm test_tools)
  add_target_exe_serial_wrapper(test_zgemm test_tools)
  add_target_exe_serial_wrapper(test_dgemm_b test_tools)
  add_target_exe_serial_wrapper(test_zgemm_b test_tools)
  add_target_exe_serial_wrapper(test_dtrsm test_tools)
  add_target_exe_serial_wrapper(test_ztrsm test_tools)
  add_target_exe_serial_wrapper(test_dgemv test_tools)
  add_target_exe_serial_wrapper(test_dgemv_ test_tools)
  add_target_exe_serial_wrapper(test_dgemv_offload test_tools)
  add_target_exe_serial_wrapper(test_zgemv test_tools)
  add_target_exe_serial_wrapper(test_zgemv_offload test_tools)
  add_target_exe_serial_wrapper(test_dskmv test_tools)
  add_target_exe_serial_wrapper(test_zskmv test_tools)
  add_target_exe_serial_wrapper(test_dgetrf test_tools)
  add_target_exe_serial_wrapper(test_zgetrf test_tools)
  add_target_exe_serial_wrapper(test_dgetrfi test_tools)
  add_target_exe_serial_wrapper(test_zgetrfi test_tools)
  add_target_exe_serial_wrapper(test_openmp_reduction_real test_tools)
  add_target_exe_serial_wrapper(test_openmp_reduction_complex test_tools)

  if(EXT_GPU)

    add_target_exe_serial_wrapper(test_offload_pointer_transcription test_tools)
    add_target_exe_serial_wrapper(test_offload_data_real test_tools)
    add_target_exe_serial_wrapper(test_offload_data_complex test_tools)
    add_target_exe_serial_wrapper(test_offload_code_assessment_real test_tools)
    add_target_exe_serial_wrapper(test_offload_code_assessment_complex test_tools)
    add_target_exe_serial_wrapper(test_offload_if test_tools)
    add_target_exe_serial_wrapper(test_offload_code_forcycle_real test_tools)
    add_target_exe_serial_wrapper(test_offload_code_forcycle_complex test_tools)

  endif()

  # These tools should be threated in the same way:
  # Also compile these tools only if serial is enabled, other it does not make sense

  if(EXT_SPEEDTEST)

    if( TARGET test_dgemm )
      add_custom_target( run_speed_dgemm
          COMMAND ${BASH_EXECUTABLE} cat ss | $<TARGET_FILE:test_dgemm>
          COMMAND ${BASH_EXECUTABLE} echo 512 N N 1000 | $<TARGET_FILE:test_dgemm>
          COMMAND ${BASH_EXECUTABLE} echo 1024 N N 100 | $<TARGET_FILE:test_dgemm>
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_dgemm_NN128/
        )
      list(APPEND SPEED_TESTS run_speed_dgemm)
    endif()

    if( TARGET test_zgemm )
      add_custom_target( run_speed_zgemm
          COMMAND ${BASH_EXECUTABLE} cat ss | $<TARGET_FILE:test_zgemm>
          COMMAND ${BASH_EXECUTABLE} echo 512 N N 1000 | $<TARGET_FILE:test_zgemm>
          COMMAND ${BASH_EXECUTABLE} echo 1024 N N 100 | $<TARGET_FILE:test_zgemm>
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_zgemm_NN128/
        )
      list(APPEND SPEED_TESTS run_speed_zgemm)
    endif()

    if( TARGET test_dgemm_b )
      add_custom_target( run_speed_dgemm_b
          COMMAND ${BASH_EXECUTABLE} cat ss | $<TARGET_FILE:test_dgemm_b>
          COMMAND ${BASH_EXECUTABLE} echo 512 N N 500 | $<TARGET_FILE:test_dgemm_b>
          COMMAND ${BASH_EXECUTABLE} echo 1024 N N 50 | $<TARGET_FILE:test_dgemm_b>
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_dgemm_b_NN128/
        )
      list(APPEND SPEED_TESTS run_speed_dgemm_b)
    endif()

    if( TARGET test_zgemm_b )
      add_custom_target( run_speed_zgemm_b
          COMMAND ${BASH_EXECUTABLE} cat ss | $<TARGET_FILE:test_zgemm_b>
          COMMAND ${BASH_EXECUTABLE} echo 512 N N 500 | $<TARGET_FILE:test_zgemm_b>
          COMMAND ${BASH_EXECUTABLE} echo 1024 N N 50 | $<TARGET_FILE:test_zgemm_b>
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_zgemm_b_NN128/
        )
      list(APPEND SPEED_TESTS run_speed_zgemm_b)
    endif()

    if( TARGET turborvb-serial )
      add_custom_target( run_speed_dft
          COMMAND ${BASH_EXECUTABLE} ./cm.stest.sh $<TARGET_FILE:prep-serial> prep.speed.o
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_dft_pbc_twist_complex/
        )
      list(APPEND SPEED_TESTS run_speed_dft)

      add_custom_target( run_speed_vmc
          COMMAND ${BASH_EXECUTABLE} ./cm.stest.sh $<TARGET_FILE:turborvb-serial> vmc.speed.o
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_vmc/
        )
      list(APPEND SPEED_TESTS run_speed_vmc)

      add_custom_target( run_speed_lrdmc
          COMMAND ${BASH_EXECUTABLE} ./cm.stest.sh $<TARGET_FILE:turborvb-serial> lrdmc.speed.o
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/test_lrdmc_new/
        )
      list(APPEND SPEED_TESTS run_speed_lrdmc)
    endif()
    
    add_custom_target( speed_tests
        DEPENDS "${SPEED_TESTS}"
        )

  endif()

endif()


######################################################################
#
# Add subdirectories
#
######################################################################

add_subdirectory(src/m_pfapack)
add_subdirectory(src/m_common)

add_subdirectory(src/d_qlapack)

add_subdirectory(src/c_adjoint_forward)
add_subdirectory(src/c_adjoint_backward)

add_subdirectory(src/b_complex)

add_subdirectory(src/a_prep)
add_subdirectory(src/a_turborvb)

add_subdirectory(src/a_tools)
add_subdirectory(src/a_readforward)
add_subdirectory(src/a_module_tests)
add_subdirectory(src/a_makefort10)
add_subdirectory(src/a_convertfort10mol)
add_subdirectory(src/a_convertfort10)

#prepare_agressive()

foreach( FLG ${EXT_AGRS} )
    list(APPEND AGRESSIVE_L ${FLG})
endforeach()
list(REMOVE_DUPLICATES AGRESSIVE_L)

foreach( FLG ${EXT_BLKL_AGRS} )
    list(REMOVE_ITEM AGRESSIVE_L ${FLG})
endforeach()


if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")

    prepare_agressive_gnu()
    set_gnu_fortran_compiler()
  
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")

    prepare_agressive_intel()
    set_intel_fortran_compiler()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")

    prepare_agressive_oneapi()
    set_oneapi_fortran_compiler()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")

    prepare_agressive_xl()
    set_xl_fortran_compiler()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")

    prepare_agressive_nvhpc()
    set_nvhpc_fortran_compiler()

else()

    message( WARNING "Unsupported Fortran compiler!" )

endif()
if(CMAKE_C_COMPILER_ID STREQUAL "GNU")

  ######################################################################
  #
  # Here GCC part starts
  #

  set_nvhpc_c_compiler()

  #
  # Here GCC part ends
  #
  ######################################################################
  
elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel")

  ######################################################################
  #
  # Here Intel (legacy) part starts
  #

  set_intel_c_compiler()
  
  #
  # Here Intel (legacy) part ends
  #
  ######################################################################

elseif(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")

  ######################################################################
  #
  # Here Intel oneAPI part starts
  #

  set_oneapi_c_compiler()
  
  #
  # Here Intel oneAPI part ends
  #
  ######################################################################

elseif(CMAKE_C_COMPILER_ID STREQUAL "XLClang")

  ######################################################################
  #
  # Here IBM XL part starts
  #

  set_xl_c_compiler()
  
  #
  # Here IBM part ends
  #
  ######################################################################

elseif(CMAKE_C_COMPILER_ID STREQUAL "NVHPC")

  ######################################################################
  #
  # Here Nvidia XL part starts
  #

  set_nvhpc_c_compiler()

  #
  # Here Nvidia part ends
  #
  ######################################################################

elseif(CMAKE_C_COMPILER_ID STREQUAL "AppleClang")

  ######################################################################
  #
  # Here Nvidia XL part starts
  #

  set_appleclang_c_compiler()

  #
  # Here Nvidia part ends
  #
  ######################################################################

else()

    message( FATAL_ERROR "This C compiler is not supportted")

endif()

# Adding preprocessor flags from input

foreach( FLG ${EXT_FLAGS} )
    list(APPEND PPDIRECTIVES ${FLG})
endforeach()
list(REMOVE_DUPLICATES PPDIRECTIVES)

# Removing blacklisted flags

foreach( FLG ${EXT_BLKL_FLAGS} )
    list(REMOVE_ITEM PPDIRECTIVES ${FLG})
endforeach()

foreach( FLG ${EXT_BLKL_FLAGS} )
	list(REMOVE_ITEM PARALLEL_F ${FLG})
endforeach()

list(SORT PPDIRECTIVES)

######################################################################
#
# There is a problem C preprocessor directives in Fortran,
# this solution is based on:
#
# https://stackoverflow.com/questions/16553227/cmake-missing-fortran-module-dependency-in-preprocessor-directive
#

foreach( NAME_OF_TARGET IN LISTS LIBRARIES_S_L)
    set(PREPROCESSOR_FLAGS ${PPDIRECTIVES})
    get_target_property(NOTTHIS ${NAME_OF_TARGET} NOFLAG)
    foreach( FL ${NOTTHIS} )
      list(REMOVE_ITEM PREPROCESSOR_FLAGS ${FL})
    endforeach()
    set_target_properties ("${NAME_OF_TARGET}" PROPERTIES COMPILE_DEFINITIONS "${PREPROCESSOR_FLAGS}")
endforeach()

foreach( NAME_OF_TARGET IN LISTS LIBRARIES_P_L)
    list(APPEND PREPROCESSOR_FLAGS ${PPDIRECTIVES} ${PARALLEL_F} )
    get_target_property(NOTTHIS ${NAME_OF_TARGET} NOFLAG)
    foreach( FL ${NOTTHIS} )
      list(REMOVE_ITEM PREPROCESSOR_FLAGS ${FL})
    endforeach()
    set_target_properties ("${NAME_OF_TARGET}" PROPERTIES COMPILE_DEFINITIONS "${PREPROCESSOR_FLAGS}")
endforeach()

foreach( NAME_OF_TARGET IN LISTS EXECUTABLES_S_L )
    set(PREPROCESSOR_FLAGS ${PPDIRECTIVES})
    get_target_property(NOTTHIS ${NAME_OF_TARGET} NOFLAG)
    foreach( FL NOTTHIS )
      list(REMOVE_ITEM PREPROCESSOR_FLAGS ${FL})
    endforeach()
    set_target_properties ("${NAME_OF_TARGET}" PROPERTIES COMPILE_DEFINITIONS "${PREPROCESSOR_FLAGS}")
endforeach()

foreach( NAME_OF_TARGET IN LISTS EXECUTABLES_P_L )
    list(APPEND PREPROCESSOR_FLAGS ${PPDIRECTIVES} ${PARALLEL_F} )
    get_target_property(NOTTHIS ${NAME_OF_TARGET} NOFLAG)
    foreach( FL NOTTHIS )
      list(REMOVE_ITEM PREPROCESSOR_FLAGS ${FL})
    endforeach()
    set_target_properties ("${NAME_OF_TARGET}" PROPERTIES COMPILE_DEFINITIONS "${PREPROCESSOR_FLAGS}")
endforeach()

#
# Important thing is you dont have to specify these flags in flags,
# they will be added automaticaly, very nice!
#
######################################################################

# Prepare testing suite

find_package(Python3)

if(Python3_FOUND)
  # Here should be python version test
  message( STATUS "Since python is available loading tests" )
  enable_testing()
  add_subdirectory(test)
endif()

summary()

