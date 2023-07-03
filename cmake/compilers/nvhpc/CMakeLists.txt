
macro(SET_NVHPC_FORTRAN_COMPILER)

  ######################################################################
  #
  # Here Nvidia compiler part starts
  #

  message( STATUS "Using Nvidia Compiler" )

  set( AGRESSIVE_FLAG_S "-O3" )
  set( PASIVE_FLAG_S "-O1/-O0" )

  add_fortran_flag_if_avail( "-g")
  add_fortran_flag_if_avail( "-Mpreprocess" "Fortran compiler does not supports -Mpreprocess")

  if(EXT_DEBUG)
      # Option -g is now always activate
      #add_fortran_flag_if_avail( "-g" )
  endif()

  list( APPEND PPDIRECTIVES "__PORT" )

  foreach( EXE IN LISTS EXECUTABLES_S_L
                        EXECUTABLES_P_L
                        LIBRARIES_S_L
                        LIBRARIES_P_L)
    if(CUDAToolkit_FOUND)
      target_link_libraries( ${EXE} PUBLIC CUDA::cudart)
      target_link_libraries( ${EXE} PUBLIC CUDA::cublas)
      target_link_libraries( ${EXE} PUBLIC CUDA::cusolver)
    endif()

    if(EXT_GPU)
      if( NOT NATIVE_OMP_B )
        add_fortran_flag_if_avail( "-mp=gpu" )
        if(NOT EXT_GPUTYPE STREQUAL "")
          add_fortran_flag_if_avail( "-gpu=${EXT_GPUTYPE}" )
        endif()
        target_link_options( ${EXE} PUBLIC "-mp=gpu" )
        if(NOT EXT_GPUTYPE STREQUAL "")
          target_link_options( ${EXE} PUBLIC "-gpu=${EXT_GPUTYPE}" )
        endif()
      endif()
    else()
      add_fortran_flag_if_avail( "-mp" )
      target_link_options( ${EXE} PUBLIC "-mp" )
    endif()


    ######################################################################
    #
    # There is a problem with nvidia compiler. When readf is link with
    # BLAS libraries it will segfault. We dont know the reason for this
    # this is a quick workaround
    #
    # After restructuring this is not longer a problem, but I am keeping
    # this lines commented for a future memo anout this problem.
    #

    #if( ${EXE} MATCHES readf)
    #    set_property(TARGET ${EXE} PROPERTY LINK_LIBRARIES "-L./" "libturbo-serial.a")
    #endif()

    #
    ######################################################################

  endforeach()

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
            else()
            endif()
        endforeach()

        ######################################################################
        #
        #

        if(${EXT_OPT})
          if(${AGR})
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS -O3)
            if(NOT ${SRC} IN_LIST AGRESSIVE_FINAL_L)
              list(APPEND AGRESSIVE_FINAL_L ${SRC})
            endif()
          else()
            set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS -O1)
            if(NOT ${SRC} IN_LIST PASIVE_FINAL_L)
              list(APPEND PASIVE_FINAL_L ${SRC})
            endif()
          endif()
        else()
          set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS -O0)
          if(NOT ${SRC} IN_LIST PASIVE_FINAL_L)
            list(APPEND PASIVE_FINAL_L ${SRC})
          endif()
        endif()

        #
        #
        ######################################################################

    endforeach()

  endforeach()

  #
  # Here Nvidia part ends
  #
  ######################################################################

endmacro()

macro(SET_NVHPC_C_COMPILER)

  list( APPEND PPDIRECTIVES "__FFTW" )
  list( APPEND PPDIRECTIVES "__USE_INTERNAL_FFTW" )

  add_c_flag_if_avail( "-g")

  if( ${EXT_OPT} )
    add_c_flag_if_avail( "-O3")
  else()
    add_c_flag_if_avail( "-O0")
  endif()

  if( ${EXT_GPU} )
    add_c_flag_if_avail( "-mp=gpu" )
    if(NOT EXT_GPUTYPE STREQUAL "")
      add_c_flag_if_avail( "-gpu=${EXT_GPUTYPE}" )
    endif()
  else()
    add_c_flag_if_avail( "-mp" )
  endif()

endmacro()

macro( PREPARE_AGRESSIVE_NVHPC )
    set(AGRESSIVE_L "")
    list( APPEND AGRESSIVE_L
      cell.f90
      conjginv_prep.f90
      dgemm_my.f90
      dielectric.f90
      dsktri.f90
      dsktrs.f90
      ewald.f90
      hopping_pbc.f90
      jastrow_ee_exp.f90
      jastrow_ei.f90
      jastrow_exp.f90
      jastrowgrad_exp_pbc.f90
      makefun.f90
      makefun0.f90
      makefun0_bump.f90
      makefun0_pbc.f90
      makefun_bump.f90
      makefun_pbc.f90
      mapping_sz.f90
      ngivej_pbc.f90
      ratio_psi.f90
      ratiofn_psi.f90
      ratiovar.f90
      scalevect.f90
      scratchdet.f90
      subener.f90
      t_lrdmc.f90
      up2bodypsi_pbc_exp.f90
      updiag.f90
      upinvhop.f90
      upinvhop_fnf_new.f90
      upnewwf_new.f90
      upsim.f90
      upsimp.f90
      uptable.f90
      uptabpip_new.f90
      uptabtot_new.f90
      upvpot_pbc_exp.f90
      upvpotaa.f90
      upvpotdiag_pbc.f90
      upwinv.f90
      upwinvp.f90
      zsktri.f90
      zsktrs.f90
        )

  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 23)
      list(REMOVE_ITEM AGRESSIVE_L upinvhop_fnf_new.f90)
  endif()
endmacro()
