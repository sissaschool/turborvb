
macro(SET_XL_FORTRAN_COMPILER)

  ######################################################################
  #
  # Here IBM XL part starts
  #

  message( STATUS "Using IBM Compiler" )
  message( WARNING "In order to compile TurboRVB properly with IBM compilers one has to precompile and link netlib LAPACK and link it with rest of the application." )

  # Use threadsafe versions of xl compiler set

  set( CMAKE_C_COMPILER "xlc_r" )
  set( CMAKE_Fortran_COMPILER "xlf_r" )

  # This globing is bad practice, we should find a way how to set these in the respective directories
  #
  # Maybe here interesting anwsers:
  # https://cmake.org/cmake/help/latest/command/set_property.html#command:set_property

  file( GLOB_RECURSE FS "AD/**/*.f" )
  file( GLOB_RECURSE F90AD "AD/**/*.f90" )
  file( GLOB_RECURSE F90SRC "src/*.f90" )
  file( GLOB_RECURSE F90DFT "DFT/*.f90" )
  file( GLOB_RECURSE F90TOOL "tool/*.f90" )
  file( GLOB_RECURSE F90TOOLPLUS "tool/**/*.f90" )
  list( APPEND F90S ${F90AD} ${F90SRC} ${F90DFT} ${F90TOOL} ${F90TOOLPLUS} )

  list( APPEND LOC_F_FLAGS "-qfixed" )
  list( APPEND LOC_F_FLAGS "-qsuffix=cpp=f" )
  list( JOIN LOC_F_FLAGS " " LOC_F )

  list( APPEND LOC_F90_FLAGS " -qfree" )
  list( APPEND LOC_F90_FLAGS " -qsuffix=cpp=f90" )
  list( JOIN LOC_F90_FLAGS " " LOC_F90 )

  if(OpenMP_FOUND)
      link_libraries(${OpenMP_Fortran_LIBRARIES})
      link_libraries(${OpenMP_C_LIBRARIES})
  endif()

  foreach ( EXE IN LISTS EXECUTABLES_S_L
                         EXECUTABLES_P_L)
    if(EXT_GPU)
      add_fortran_flag_if_avail( "-qoffload" )
      target_link_options( ${EXE} PUBLIC "-qoffload" )
    endif()
    set( _essl_path $ENV{ESSL_LIB} )
    target_link_libraries( ${EXE} PUBLIC "${_essl_path}/libesslsmp.so" )
    target_link_options( ${EXE} PUBLIC "-qsmp=omp:noauto" )
    #target_link_options( ${EXE} PUBLIC "-L${_essl_path}" "-lesslsmp" )
  endforeach()

  add_fortran_flag_if_avail( "-q64" )
  #set_source_files_properties( ${FS} PROPERTIES COMPILE_FLAGS "${LOC_F}" )
  #set_source_files_properties( ${F90S} PROPERTIES COMPILE_FLAGS "${LOC_F90}" )
  add_fortran_flag_if_avail( "-qfloat=norngchk" )
  #add_fortran_flag_if_avail( "-qsmp=omp:noauto" )
  add_fortran_flag_if_avail( "-qsmp=omp" )
  add_fortran_flag_if_avail( "-qsimd=auto" )
  add_fortran_flag_if_avail( "-Q" )
  add_fortran_flag_if_avail( "-qstrict" )
  add_fortran_flag_if_avail( "-qcheck=all" )
  add_fortran_flag_if_avail( "-qmaxmem=-1" )

  list( APPEND PPDIRECTIVES "RISC" )

  # Set optimization flags:
  foreach( TARGET IN LISTS LIBRARIES_S_L
                           LIBRARIES_P_L
		           EXECUTABLES_S_L
		           EXECUTABLES_P_L )
      if(CUDAToolkit_FOUND)
	target_link_libraries( ${TARGET} PUBLIC CUDA::cudart)
        target_link_libraries( ${TARGET} PUBLIC CUDA::cublas)
        target_link_libraries( ${TARGET} PUBLIC CUDA::cusolver)
      endif()
      get_target_property(TURBO_SOURCES ${TARGET} SOURCES)
      foreach(SRC ${TURBO_SOURCES})

	  set( SP "" )
	  if ( ${SRC} MATCHES ".*\.f$" )
	    set( SP ${LOC_F} )
	  endif()
	  if ( ${SRC} MATCHES ".*\.f90$" )
	    set( SP ${LOC_F90} )
	  endif()

          set( AGR "OFF" )
	  if(${TARGET} MATCHES pfapack)
              set( AGR "ON" )
          endif()
          foreach(ARGSRC ${AGRESSIVE_L})
              if(${SRC} MATCHES ".*${ARGSRC}")
                  set( AGR "ON" )
              endif()
          endforeach()
          if(${AGR} AND ${EXT_OPT})
	      set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "${SP} -Ofast")
              if(NOT ${SRC} IN_LIST AGRESSIVE_FINAL_L)
                list(APPEND AGRESSIVE_FINAL_L ${SRC})
              endif()
          else()
	      set_source_files_properties(${SRC} PROPERTIES COMPILE_FLAGS "${SP} -O2")
              if(NOT ${SRC} IN_LIST PASIVE_FINAL_L)
                list(APPEND PASIVE_FINAL_L ${SRC})
              endif()
          endif()
      endforeach()
  endforeach()

  if(${EXT_LTO} AND ${EXT_OPT})
      add_fortran_flag_if_avail( "-flto" )
      add_link_options( "-flto" )
  endif()

  #
  # Here IBM part ends
  #
  ######################################################################

endmacro()

macro(SET_XL_C_COMPILER)

  list( APPEND PPDIRECTIVES "__FFTW" )
  list( APPEND PPDIRECTIVES "__USE_INTERNAL_FFTW" )

  foreach ( EXE IN LISTS EXECUTABLES_S_L
                         EXECUTABLES_P_L)
    if(EXT_GPU)
      add_c_flag_if_avail( "-qoffload" )
    endif()

    add_c_flag_if_avail( "-qsmp=omp:noauto" )
    add_c_flag_if_avail( "-q64" )
    add_c_flag_if_avail( "-O2" )

  endforeach()

endmacro()

macro( PREPARE_AGRESSIVE_XL )
    set(AGRESSIVE_L "")
    list( APPEND AGRESSIVE_L 
      cell.o 
      dielectric.o 
      ewald.o 
      makefun.o 
      makefun_pbc.o 
      makefun_bump.f90
      makefun0.f90
      makefun0_pbc.f90
      makefun0_bump.f90
      upnewwf_new.f90
      upinvhop.f90
      upinvhop_fnf_new.f90
      uptabpip_new.f90
      uptabtot_new.f90
      ratiovar.f90
      scratchdet.f90
      scalevect.f90
      jastrow_exp.f90
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
      upvpotdiag_pbc.f90
      t_lrdmc.f90
      ngivej_pbc.f90
      upvpotaa.f90
      dgemm_my.f90
      subener.f90
      ratiofn_psi.f90
      ratio_psi.f90
      dsktri.f90
      dsktrs.f90
      zsktrs.f90
      zsktri.o 
        )
endmacro()
