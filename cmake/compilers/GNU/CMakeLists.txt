
macro(SET_GNU_FORTRAN_COMPILER)

  ######################################################################
  #
  # Here GCC part starts
  #

  message( STATUS "Using GCC Compiler" )

  set( AGRESSIVE_FLAG_S "-Ofast" )
  set( AGRESSIVE_FLAG_S "${AGRESSIVE_FLAG_S} ${EXT_AGRS}" )
  if( ${EXT_OPT} )
    set( PASIVE_FLAG_S "-O2" )
  else()
    set( PASIVE_FLAG_S "-O0" )
  endif()

  # Maybe this can be done by check_Fortran_compiler_flag

  if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "10.0.0")
    message( STATUS "GNU Fortran with \"no-argument-mismatch\" flag")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-argument-mismatch ")
  else()
    message( STATUS "GNU Fortran with \"allow-argument-mismatch\" flag")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch ")
  endif()

  # Checking if compiler knows these flags:

  check_Fortran_compiler_flag("-ffree-line-length-none" _free_line_length_none_works)

  if (${EXT_PROFILED})
    add_fortran_flag_if_avail( "-pg" )
  endif()

  if(EXT_GPU)
	  add_fortran_flag_if_avail( "-fno-whole-file" "GNU Fortran compiler does not supports no-whole-flag, turn off GPU support" )
	  add_fortran_flag_if_avail( "-fopenmp" )
	  add_fortran_flag_if_avail( "-foffload=nvptx-none" )
	  add_fortran_flag_if_avail( "-foffload-options=-lgfortran" )
	  add_fortran_flag_if_avail( "-foffload-options=-lm" )
	  add_fortran_flag_if_avail( "-foffload-options=-latomic" )
	  add_c_flag_if_avail( "-fopenmp" )
	  add_c_flag_if_avail( "-foffload=nvptx-none" )
	  add_c_flag_if_avail( "-foffload-options=-lgfortran" )
	  add_c_flag_if_avail( "-foffload-options=-lm" )
	  add_c_flag_if_avail( "-foffload-options=-latomic" )
      foreach( TARGET IN LISTS EXECUTABLES_S_L
                               EXECUTABLES_P_L
                               LIBRARIES_S_L
                               LIBRARIES_P_L)
        target_link_options( ${TARGET} PRIVATE "-fopenmp" )
        target_link_options( ${TARGET} PRIVATE "-foffload=nvptx-none" )
	      target_link_options( ${TARGET} PRIVATE "-foffload-options=-lgfortran" )
	      target_link_options( ${TARGET} PRIVATE "-foffload-options=-lm" )
	      target_link_options( ${TARGET} PRIVATE "-foffload-options=-latomic" )
        if(CUDAToolkit_FOUND)
          target_link_libraries( ${TARGET} PUBLIC CUDA::cudart)
          target_link_libraries( ${TARGET} PUBLIC CUDA::cublas)
          target_link_libraries( ${TARGET} PUBLIC CUDA::cusolver)
        endif()
      endforeach()
  endif()
  
  foreach( TARGET IN LISTS EXECUTABLES_S_L
                           EXECUTABLES_P_L
                           LIBRARIES_S_L
                           LIBRARIES_P_L)
    if (${EXT_PROFILED})
      target_link_options( ${TARGET} PRIVATE "-pg" )
    endif()
  endforeach()

  add_fortran_flag_if_avail( "-march=native" )
  add_fortran_flag_if_avail( "-mtune=native" )
  add_fortran_flag_if_avail( "-w" )
  add_fortran_flag_if_avail( "-g" )
  add_fortran_flag_if_avail( "-ggdb3" )
  add_fortran_flag_if_avail( "-cpp"                    "GNU Fortran compiler does not supports c preprocessor" )
  add_fortran_flag_if_avail( "-ffree-line-length-none" "GNU Fortran compiler does not supports arbitrary line length" )
  add_link_options( "-ggdb3" )
  add_link_options( "-static" )

  # Set optimization flags:
  foreach( LIBRARY IN LISTS LIBRARIES_S_L
                            LIBRARIES_P_L
                            EXECUTABLES_S_L
                            EXECUTABLES_P_L
                            )

      get_target_property(TURBO_SOURCES ${LIBRARY} SOURCES)
      foreach(SRC ${TURBO_SOURCES})
        set( AGR "OFF" )
        if(${LIBRARY} MATCHES pfapack)
            set( AGR "ON" )
        endif()
        foreach(ARGSRC ${AGRESSIVE_L})
          if(${SRC} MATCHES ".*${ARGSRC}")
              set( AGR "ON" )
          endif()
        endforeach()
        if(${AGR} AND ${EXT_OPT})
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS ${AGRESSIVE_FLAG_S})
            if(NOT ${SRC} IN_LIST AGRESSIVE_FINAL_L)
              list(APPEND AGRESSIVE_FINAL_L ${SRC})
            endif()
        else()
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS ${PASIVE_FLAG_S})
            if(NOT ${SRC} IN_LIST PASIVE_FINAL_L)
              list(APPEND PASIVE_FINAL_L ${SRC})
            endif()
        endif()
      endforeach()
  endforeach()

  if(${EXT_LTO} AND ${EXT_OPT} )
      add_fortran_flag_if_avail( "-flto" )
      foreach( LIBRARY IN LISTS LIBRARIES_S_L
                                LIBRARIES_P_L )
                            target_link_options( ${LIBRARY} PRIVATE "-flto -O3" )
      endforeach()
  endif()

  #
  # Here GCC part ends
  #
  ######################################################################

endmacro()

macro(SET_GNU_C_COMPILER)

  list( APPEND PPDIRECTIVES "__FFTW" )
  list( APPEND PPDIRECTIVES "__USE_INTERNAL_FFTW" )

  add_c_flag_if_avail( "-g" )

  if (${EXT_PROFILED})
    add_c_flag_if_avail( "-pg" )
  endif()

  if( ${EXT_OPT} )
    add_c_flag_if_avail( "-O3" )
  else()
    add_c_flag_if_avail( "-O0" )
  endif()

endmacro()

macro( PREPARE_AGRESSIVE_GNU )
  set(AGRESSIVE_L "")
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

  # The cell.o object has to be removed from agressive flags for gcc 11
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 11
      AND
      CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 12.0.0)
    list(REMOVE_ITEM AGRESSIVE_L cell.f90)
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 11.3)
      list(REMOVE_ITEM AGRESSIVE_L jastrowgrad_exp_pbc.f90)
    endif()
  endif()

endmacro()
