
macro(SET_ONEAPI_FORTRAN_COMPILER)

  ######################################################################
  #
  # Here Intel oneAPI part starts
  #

  message( STATUS "Using Intel oneAPI Compiler" )

  add_fortran_flag_if_avail( "-fpp"   "Intel LLVM Fortran compiler does not supports C preprocessor" )
  add_fortran_flag_if_avail( "-xHost" )
  add_fortran_flag_if_avail( "-finline-functions" )
  add_fortran_flag_if_avail( "-qopt-report" )

  set(AGRESSIVE_F "-Ofast")
  if( ${EXT_OPT} )
    set(PASIVE_F "-O1")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 2023.1)
      set(PASIVE_F "-O0")
    endif()
  else()
    set(PASIVE_F "-O0")
  endif()

  if( ${EXT_PROFILED} )
    add_fortran_flag_if_avail( "-pg" )
  endif()

  foreach( LIBRARY IN LISTS LIBRARIES_S_L
                            LIBRARIES_P_L
                            EXECUTABLES_S_L
                            EXECUTABLES_P_L
                            )
    if( ${EXT_PROFILED} )
      target_link_options( ${LIBRARY} PRIVATE "-pg" )
    endif()
    if( ${EXT_OMP} )
      target_link_options( ${LIBRARY} PRIVATE "-fopenmp" )
    endif()
  endforeach()

  # Set optimization flags:
  foreach( LIBRARY IN LISTS LIBRARIES_S_L
                            LIBRARIES_P_L
                            EXECUTABLES_S_L
                            EXECUTABLES_P_L
                            )
      if(${EXT_LTO} AND ${EXT_OPT})
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
      endif()
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
            if(NOT ${SRC} IN_LIST AGRESSIVE_FINAL_L)
              list(APPEND AGRESSIVE_FINAL_L ${SRC})
            endif()
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS ${AGRESSIVE_F})
        else()
            if(NOT ${SRC} IN_LIST PASIVE_FINAL_L)
              list(APPEND PASIVE_FINAL_L ${SRC})
            endif()
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS ${PASIVE_F})
        endif()
      endforeach()
  endforeach()
  
  #
  # Here Intel oneAPI part ends
  #
  ######################################################################

endmacro()

macro(SET_ONEAPI_C_COMPILER)

  list( APPEND PPDIRECTIVES "__FFTW" )
  list( APPEND PPDIRECTIVES "__USE_INTERNAL_FFTW" )

  add_c_flag_if_avail( "-g" )

  if( ${EXT_OPT} )
    add_c_flag_if_avail( "-O3" )
  else()
    add_c_flag_if_avail( "-O0" )
  endif()

endmacro()

macro( PREPARE_AGRESSIVE_ONEAPI )
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

  if (CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL 2023.1)
    list(REMOVE_ITEM AGRESSIVE_L ewald.f90)
  endif()
endmacro()
