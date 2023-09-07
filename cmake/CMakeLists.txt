######################################################################
#

add_subdirectory(compilers)

#
######################################################################

######################################################################
#
# These wrappers add targets to lists for futher processing
# but only if they are missing there
#

macro(ADD_TARGET_LIB_SERIAL_WRAPPER)
  if(NOT ${ARGV0} IN_LIST LIBRARIES_S_L)
    list( APPEND LIBRARIES_S_L ${ARGV0} )
    add_library(${ARGV0})
    install( TARGETS ${ARGV0} )
  endif()
endmacro()

macro(ADD_TARGET_EXE_SERIAL_WRAPPER)
  if(NOT ${ARGV0} IN_LIST EXECUTABLES_S_L)
    list( APPEND EXECUTABLES_S_L ${ARGV0} )
    add_executable(${ARGV0})
    if("${ARGV1}" STREQUAL "")
      install( TARGETS ${ARGV0} )
    else()
      install( TARGETS ${ARGV0}
               RUNTIME DESTINATION ${ARGV1})
    endif()
  endif()
endmacro()

macro(ADD_TARGET_LIB_PARALLEL_WRAPPER)
  if(NOT ${ARGV0} IN_LIST LIBRARIES_P_L)
    list( APPEND LIBRARIES_P_L ${ARGV0} )
    add_library(${ARGV0})
    install( TARGETS ${ARGV0} )
  endif()
endmacro()

macro(ADD_TARGET_EXE_PARALLEL_WRAPPER)
  if(NOT ${ARGV0} IN_LIST EXECUTABLES_P_L)
    list( APPEND EXECUTABLES_P_L ${ARGV0} )
    add_executable(${ARGV0})
    if("${ARGV1}" STREQUAL "")
      install( TARGETS ${ARGV0} )
    else()
      install( TARGETS ${ARGV0}
               RUNTIME DESTINATION ${ARGV1})
    endif()
  endif()
endmacro()

#
######################################################################

macro(ADD_TARGET_LINK_FLAG)
    foreach( TARGET IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L
                             LIBRARIES_S_L
                             LIBRARIES_P_L)
        target_link_options( ${TARGET} PRIVATE ${ARGV0} )
    endforeach()
endmacro()

macro(SET_MODULE_PATH_DIRECTORY)
    if( ${ARGV1} MATCHES "serial")
        set_target_properties(${ARGV0} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/serial_modules)
        target_include_directories(${ARGV0} PUBLIC ${CMAKE_BINARY_DIR}/serial_modules)
    else()
        if( ${ARGV0} MATCHES ".*-serial" )
            set_target_properties(${ARGV0} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/serial_modules)
            target_include_directories(${ARGV0} PUBLIC ${CMAKE_BINARY_DIR}/serial_modules)
        endif()

        if( ${ARGV0} MATCHES ".*-mpi" )
          set_target_properties(${ARGV0} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mpi_modules)
          target_include_directories(${ARGV0} PUBLIC ${CMAKE_BINARY_DIR}/mpi_modules)
        endif()
    endif()
endmacro()

# We have to find a way how to make summaries

macro(SUMMARY)
  
  message( "" )
  message( " ######################### SUMMARY ######################### " )
  message( "" )

  get_filename_component( FC "${CMAKE_Fortran_COMPILER}" NAME )
  get_filename_component( CC "${CMAKE_C_COMPILER}" NAME )
  message( "   C compiler:         ${CMAKE_C_COMPILER_ID} (${CC} ver. ${CMAKE_C_COMPILER_VERSION})" )
  message( "   Fortran compiler:   ${CMAKE_Fortran_COMPILER_ID} (${FC} ver. ${CMAKE_Fortran_COMPILER_VERSION})" )

  message( "" )
  message( "   BLAS used:" )
  message( "   ${BLAS_LIBRARIES}" )

  message( "" )
  message( "   Compiling QMC code:          ${EXT_QMC}" )
  message( "   Compiling DFT code:          ${EXT_DFT}" )
  message( "   Compiling tools:             ${EXT_TOOLS}" )
  message( "" )
  message( "   Compiling serial version:    ${EXT_SERIAL}" )
  message( "   Compiling parallel version:  ${EXT_PARALLEL}" )
  message( "   Compiling GPU:               ${EXT_GPU}" )
  message( "" )
  message( "   Optimizations:               ${EXT_OPT}" )
  if(EXT_LTO)
    set( LTO_WARNING " (Beware, turning LTO on can make compilation very long (hours))" )
  endif()
  message( "   LTO optimizations:           ${EXT_LTO}${LTO_WARNING}" )

  # Showing preprocessor flags

  message( "" )
  message( "   Base preprocessor flags are set to:")
  foreach(DIRECTIVE IN LISTS PPDIRECTIVES)
    message("     ::  ${DIRECTIVE}")
  endforeach()

  message( "" )
  message( "   Parallel specific preprocessor flags are set to:")
  foreach(DIRECTIVE IN LISTS PARALLEL_F)
    message("     ::  ${DIRECTIVE}")
  endforeach()

  message( "" )

  # Showing compiler flags

  message( "   Base Fortran flags:" )
  foreach(FLAG ${CMAKE_Fortran_FLAGS})
  	message("     :: ${FLAG}")
  endforeach()

  if(NOT EXT_DETAIL)
    list( LENGTH AGRESSIVE_FINAL_L NUM_AGR )
    list( LENGTH PASIVE_FINAL_L NUM_NONAGR )
    message( "" )
    message( "   Agressive flags:         ${AGRESSIVE_FLAG_S}  (${NUM_AGR} files)" )
    message( "   Non-agressive flags:     ${PASIVE_FLAG_S}  (${NUM_NONAGR} files)" )
  else()
    message( "" )
    message( "<Here should be list of agresive and non-agressive object files>" )
    foreach(X ${AGRESSIVE_FINAL_L})
        message(${X})
    endforeach()
  endif()

  message( "" )
  message( " ######################### xxxxxxx ######################### " )
  message( "" )

endmacro()

# This macro sould apply SAFE flags
macro(MAKE_SAVE)
  set(SAFE_OPTIONS ${ARGV0})
  set(SOURCE_FILES ${ARGV1})
  foreach(FL ${SAFE_FILES})
    set_source_files_properties(${FL} PROPERTIES COMPILE_FLAGS "${SAFE_OPTIONS}")
  endforeach()
endmacro()

######################################################################
#
# This macro will add Fortran compiler flags to the targets defined in
# lists LIBRARIES_(S|P)_L and EXECUTABLES_(S|P)_L. If a second 
# argument is defined macro will consider this flag as mandatory. 
# Then the second argument will be considered as the error message 
# if flag does not exists.
#

macro(ADD_FORTRAN_FLAG_IF_AVAIL)
  set(FLAG ${ARGV0})
  set(MESSAGE ${ARGV1})
  if(NOT "${FLAG}" STREQUAL "")
    string(REGEX REPLACE "^-" "fortran_" FLAG_NAME "${FLAG}")
    string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
    string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
    string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")
    
    # runs check to see flag is supported by compiler
    check_fortran_compiler_flag("${FLAG}" ${FLAG_NAME})
    if(${FLAG_NAME})
      foreach(NAME_OF_TARGET IN LISTS LIBRARIES_S_L
                                      LIBRARIES_P_L)
        target_compile_options(${NAME_OF_TARGET} PRIVATE
                               $<$<COMPILE_LANGUAGE:Fortran>:${FLAG}>)
      endforeach()
      foreach(NAME_OF_TARGET IN LISTS EXECUTABLES_S_L
                                      EXECUTABLES_P_L)
        target_compile_options(${NAME_OF_TARGET} PRIVATE
                               $<$<COMPILE_LANGUAGE:Fortran>:${FLAG}>)
      endforeach()
    else()
      if(NOT "${MESSAGE}" STREQUAL "")
        message( FATAL_ERROR ${MESSAGE} )
      endif()
    endif()
  endif()
endmacro()

#
######################################################################

######################################################################
#
# This macro will add C compiler flags to the targets defined in
# lists LIBRARIES_(S|P)_L and EXECUTABLES_(S|P)_L. If a second 
# argument is defined macro will consider this flag as mandatory. 
# Then the second argument will be considered as the error message 
# if flag does not exists.
#

macro(ADD_C_FLAG_IF_AVAIL)
  set(FLAG ${ARGV0})
  set(MESSAGE ${ARGV1})
  if(NOT "${FLAG}" STREQUAL "")
    string(REGEX REPLACE "^-" "c_" FLAG_NAME "${FLAG}")
    string(REPLACE "-" "_" FLAG_NAME "${FLAG_NAME}")
    string(REPLACE " " "_" FLAG_NAME "${FLAG_NAME}")
    string(REPLACE "=" "_" FLAG_NAME "${FLAG_NAME}")
    
    # runs check to see flag is supported by compiler
    check_c_compiler_flag("${FLAG}" ${FLAG_NAME})
    if(${FLAG_NAME})
      foreach(NAME_OF_TARGET IN LISTS LIBRARIES_S_L
                                      LIBRARIES_P_L)
        target_compile_options(${NAME_OF_TARGET} PRIVATE
                               $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
      endforeach()
      foreach(NAME_OF_TARGET IN LISTS EXECUTABLES_S_L
                                      EXECUTABLES_P_L)
        target_compile_options(${NAME_OF_TARGET} PRIVATE
                               $<$<COMPILE_LANGUAGE:C>:${FLAG}>)
      endforeach()
    else()
      if(NOT "${MESSAGE}" STREQUAL "")
        message( FATAL_ERROR ${MESSAGE} )
      endif()
    endif()
  endif()
endmacro()

#
######################################################################

######################################################################
#
# Check if library is added by hand if not try to find it:
#
# This has to be changed
#

macro( ADD_OR_FIND_PACKAGE LIBR)
  set(VAR_EXT "EXT_${LIBR}_LIB")
  set(VAR_EXT "EXT_${LIBR}_LIB")
  cmake_language( EVAL CODE "set(LOC_LIB_PATH ${VAR_EXT})" )
  message("${LOC_LIB_PATH} ${VAR_EXT}")
  #set(LIBR_ [[${VAR_EXT}]] )
  #;message(cmake_language(EVAL CODE "\"${LIBR_}\""))
  #if( re_evaluated[[${VAR_EXT}]] )
  #  link_libraries( ${VAR_EXT} )
  #else()
  #  find_package(LAPACK)
  #endif()
endmacro()

#
######################################################################

macro( PREPARE_AGRESSIVE )
    list( APPEND AGRESSIVE_L 
        cell.f90
        dielectric.f90
        ewald.f90 
        makefun.f90 
        makefun0.f90 
        makefun_pbc.f90 
        makefun0_pbc.f90 
        makefun_bump.f90 
        makefun0_bump.f90  
        upnewwf_new.f90 
        upinvhop.f90
        upinvhop_fnf_new.f90 
        uptabpip_new.f90 
        ratiovar.f90 
        scratchdet.f90 
        scalevect.f90
        jastrow_exp.f90 
        jastrowgrad_exp.f90 
        jastrowgrad_exp_pbc.f90 
        up2bodypsi_pbc_exp.f90 
        upwinv.f90 
        upwinvp.f90
        jastrow_ei.f90 
        jastrow_ee_exp.f90 
        conjginv_prep.f90 
        hopping_pbc.f90
        mapping_sz.f90
        updiag.f90
        upsim.f90
        upsimp.f90
        uptable.f90
        upvpot_pbc_exp.f90
        t_lrdmc.f90
        ngivej_pbc.f90
        dgemm_my.f90
        subener.f90
        ratiofn_psi.f90
        ratio_psi.f90
        dsktri.f90
        dsktrs.f90
        zsktrs.f90
        zsktri.f90
        )
endmacro()
