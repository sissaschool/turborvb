message( STATUS "Processing subdirectory tool")

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )

  if(${EXECUTABLE} STREQUAL testad)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE} 
      PRIVATE
        ${EXECUTABLE}.f90
      )

    if(EXT_GPU)
      target_link_libraries( ${EXECUTABLE}
        PRIVATE
            complex-nogpu-serial
        )
    else()
      target_link_libraries( ${EXECUTABLE}
        PRIVATE
            complex-serial
        )
    endif()
  endif()

  ######################################################################
  #
  # Add sources and link libraries for zdotc wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_zdotc)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    set_module_path_directory( ${EXECUTABLE} serial )

    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )

    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        qlapack-serial
        common-serial
      )
  endif()

  #
  ######################################################################

  ######################################################################
  #
  # Add sources and link libraries for [d,z]ger2? wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_[dz]geru?)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    set_module_path_directory( ${EXECUTABLE} serial )

    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )

    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        common-serial
        qlapack-serial
      )
  endif()

  #
  ######################################################################

  ######################################################################
  #
  # Add sources and link libraries for [dz]gem[mv](_b)?_? wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_[dz]gem[mv][\(_b\)\(_\)\(_offload\)]?)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE} 
      PRIVATE
      ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        qlapack-serial
        common-serial
        adjoint-backward-serial
      )
  endif()

  #
  ######################################################################

  ######################################################################
  #
  # Add sources and link libraries for [dz]skmv wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_[dz]skmv)
	  set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
  
	  target_sources( ${EXECUTABLE}
      PRIVATE
        ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        pfapack
        common-serial
      )
  endif()

  #
  ######################################################################

  ######################################################################
  #
  # Add sources and link libraries for [dz]getrfi? wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_[dz]getrfi?)
      set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        qlapack-serial
        common-serial
      )
  endif()

  #
  ######################################################################

  ######################################################################
  #
  # Add sources and link libraries for dgetrsm wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_[dz]trsm)
      set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        qlapack-serial
        common-serial
      )
  endif()

  ######################################################################
  #
  # Add sources and link libraries for zsktr[ri] wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_zsktr[si])
      set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        common-serial
      )
  endif()

  ######################################################################
  #
  # Add sources and link libraries for upwinvp wrapper tests
  #

  if(${EXECUTABLE} MATCHES test_upwinvp.*)
      set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
      ${EXECUTABLE}.f90
      )
    
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        common-serial
      )
  endif()

  ######################################################################
  #
  #
  
  if(${EXECUTABLE} STREQUAL test_offload_pointer_transcription)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE} 
      PRIVATE
        offload_pointer_transcription.f90
        transcribe.c
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_data_real)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
  
    target_sources( ${EXECUTABLE} 
      PRIVATE
        offload_data_real.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_data_complex)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
  
    target_sources( ${EXECUTABLE} 
      PRIVATE
        offload_data_complex.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_code_assessment_real)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        offload_code_assessment_real.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_code_assessment_complex)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        offload_code_assessment_complex.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_if)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE} 
      PRIVATE
        offload_if.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_code_forcycle_real)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        offload_code_forcycle_real.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_offload_code_forcycle_complex)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        offload_code_forcycle_complex.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_openmp_reduction_real)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        openmp_reduction_real.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_openmp_reduction_complex)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        openmp_reduction_complex.f90
      )
  endif()

  if(${EXECUTABLE} STREQUAL test_makefun)
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")
    
    target_sources( ${EXECUTABLE}
      PRIVATE
        ${EXECUTABLE}.f90
        c_lister.h
        c_lister.c
      )

    if(EXT_GPU)
      target_link_libraries( ${EXECUTABLE}
        PRIVATE
            complex-nogpu-serial
        )
    else()
      target_link_libraries( ${EXECUTABLE}
        PRIVATE
            complex-serial
        )
    endif()
  endif()

endforeach()
