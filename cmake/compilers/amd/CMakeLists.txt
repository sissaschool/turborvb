
macro(SET_AMD_FORTRAN_COMPILER)

  ######################################################################
  #
  # Here AMD part starts
  #

  set( AGRESSIVE_FLAG_S "-Ofast" )
  set( AGRESSIVE_FLAG_S "${AGRESSIVE_FLAG_S} ${EXT_AGRS}" )
  if( ${EXT_OPT} )
    set( PASIVE_FLAG_S "-O2" )
  else()
    set( PASIVE_FLAG_S "-O0" )
  endif()

  list( APPEND PPDIRECTIVES "__AMD" )

  add_fortran_flag_if_avail( "-cpp"                    "GNU Fortran compiler does not supports c preprocessor" )
  add_fortran_flag_if_avail( "-g" )
  add_fortran_flag_if_avail( "-march=native" )
  add_fortran_flag_if_avail( "-funroll-loops" )

  if( ${EXT_OPT} )
    add_fortran_flag_if_avail( "-fvectorize" )
    add_fortran_flag_if_avail( "-function-specialize" )
    add_fortran_flag_if_avail( "-finline-aggressive" )
    add_fortran_flag_if_avail( "-ffast-math" )
    add_fortran_flag_if_avail( "-freciprocal-math" )
    add_fortran_flag_if_avail( "-reduce-array-computation=3" )
    add_fortran_flag_if_avail( "-fprofile-instr-generate" )
    add_fortran_flag_if_avail( "-fprofile-instr-use" )
    add_fortran_flag_if_avail( "-fnt-store" )
  endif()

  add_fortran_flag_if_avail( "-openmp" )

  add_link_options("-qopenmp")

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

  #
  # Here GCC part ends
  #
  ######################################################################

endmacro()

macro(SET_AMD_C_COMPILER)

  list( APPEND PPDIRECTIVES "__FFTW" )
  list( APPEND PPDIRECTIVES "__USE_INTERNAL_FFTW" )

  add_c_flag_if_avail( "-g" )
  add_c_flag_if_avail( "-march=native" )
  add_c_flag_if_avail( "-openmp" )

  if( ${EXT_OPT} )
    add_c_flag_if_avail( "-O3" )
  else()
    add_c_flag_if_avail( "-O0" )
  endif()

endmacro()

macro( PREPARE_AGRESSIVE_AMD )
  set(AGRESSIVE_L "")
  list( APPEND AGRESSIVE_L 
	  #cell.f90
	  #dielectric.f90
	  #ewald.f90
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
      #jastrow_ei.f90
      #jastrow_ee_exp.f90
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
