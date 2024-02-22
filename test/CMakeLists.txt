# Testing tests:

  ######################################################################
  #
  # This is "add_test" wrapper. Basically it works as normal add_test
  # function, however, I accepts an extra keyword variable DEPENDENCY_TREE
  # This variable act as a queue for test. Every thes with the same
  # DEPENDENCY_TREE will be executed sequentualy as the were added
  #

  macro(ADD_TEST_DEPENDENCY_TREE)

    set(__KEY_VALUES "NAME;COMMAND;WORKING_DIRECTORY;DEPENDENCY_TREE")

    foreach( __KEY IN ITEMS ${__KEY_VALUES} )
      set("__${__KEY}" "")
    endforeach()

    set(__POINTER "")
    set(__LINE "")

    foreach(__ARG IN ITEMS ${ARGN})
      if( "${__ARG}" IN_LIST __KEY_VALUES )
        set( __POINTER ${__ARG} )
      else()
	list(APPEND "__${__POINTER}" "${__ARG}")
      endif()
    endforeach()

    list(REMOVE_ITEM __KEY_VALUES DEPENDENCY_TREE)

    foreach( __KEY IN ITEMS ${__KEY_VALUES} )
      if( NOT "${__${__KEY}}" STREQUAL "" )
        string(REPLACE ";" " " __VALUES "${__${__KEY}}")
        string(CONCAT __LINE "${__LINE} " "${__KEY} " "\"${__VALUES}\"")
      endif()
    endforeach()
    
    foreach( __KEY NAME;COMMAND;DEPENDENCY_TREE ) 
      if( "${__${__KEY}}" STREQUAL "" )
	message( FATAL_ERROR "${__KEY} is missing" )
      endif()
    endforeach()

    add_test( NAME ${__NAME}
              COMMAND ${__COMMAND})

    list(REMOVE_ITEM __KEY_VALUES NAME)
    list(REMOVE_ITEM __KEY_VALUES COMMAND)

    foreach( __KEY IN ITEMS ${__KEY_VALUES} )
      if( NOT "${__${__KEY}}" STREQUAL "" )
        set_tests_properties(${__NAME} PROPERTIES ${__KEY} "${__${__KEY}}")
      endif()
    endforeach()

    if( NOT "${__HIDDEN_${__DEPENDENCY_TREE}}" STREQUAL "" )
      set_tests_properties(${__NAME} PROPERTIES DEPENDS "${__HIDDEN_${__DEPENDENCY_TREE}}")
    endif()

    list(APPEND "__HIDDEN_${__DEPENDENCY_TREE}" "${__NAME}" )

  endmacro()

  #
  # end of add_test_dependency_tree macro
  #
  ######################################################################

  ######################################################################
  #
  # Testing functionalities of wrappers and routines
  #

if(EXT_MODTEST)

  if(EXT_GPU)

    add_subdirectory(test_offload_pointer_trancription)
    add_subdirectory(test_offload_data)
    add_subdirectory(test_offload_code_assessment)
    add_subdirectory(test_offload_if)
    add_subdirectory(test_offload_code_forcycle)

  endif()

  add_subdirectory(test_openmp_reduction)

  add_subdirectory(test_zdotc_128)

  add_subdirectory(test_dger_256)
  add_subdirectory(test_dger2_256)
  add_subdirectory(test_zgeru_256)
  add_subdirectory(test_zger2_256)

  add_subdirectory(test_dgemm_NN128)
  add_subdirectory(test_dgemm_NT128)
  add_subdirectory(test_dgemm_TT128)

  add_subdirectory(test_zgemm_NN128)
  add_subdirectory(test_zgemm_NT128)
  add_subdirectory(test_zgemm_CT128)

  add_subdirectory(test_dgemm_b_NN128)
  add_subdirectory(test_dgemm_b_NT128)
  add_subdirectory(test_dgemm_b_TT128)

  add_subdirectory(test_zgemm_b_NN128)
  add_subdirectory(test_zgemm_b_NT128)
  add_subdirectory(test_zgemm_b_TT128)

  add_subdirectory(test_dtrsm_ULN128)
  add_subdirectory(test_dtrsm_LRT128)

  add_subdirectory(test_ztrsm_URN128)
  add_subdirectory(test_ztrsm_LLC128)

  add_subdirectory(test_dgemv_N256)
  add_subdirectory(test_dgemv_T256)

  add_subdirectory(test_zgemv_N256)
  add_subdirectory(test_zgemv_T256)
  add_subdirectory(test_zgemv_C256)

  add_subdirectory(test_dskmv_X256)
  add_subdirectory(test_zskmv_X256)

  add_subdirectory(test_dgetrf_128)
  add_subdirectory(test_zgetrf_128)
  add_subdirectory(test_dgetrfi_128)
  add_subdirectory(test_zgetrfi_128)
  
  add_subdirectory(test_zsktri_128_L)
  add_subdirectory(test_zsktri_128_U)

  add_subdirectory(test_zsktrs_128_L)
  add_subdirectory(test_zsktrs_128_U)

  add_subdirectory(test_upwinvp)
  add_subdirectory(test_upwinvp_complex)
  #add_subdirectory(test_upwinvp_pfaff)
  #add_subdirectory(test_upwinvp_pfaff_complex)
  
  #
  ######################################################################

endif()

if(EXT_DFT)
  
  ######################################################################
  #
  # These test are the same for serial and parallel version
  #
  
  add_subdirectory(test_dft_open)
  add_subdirectory(test_dft_pbc_gamma)
  add_subdirectory(test_dft_pbc_twist_real)
  add_subdirectory(test_dft_pbc_twist_complex)

  #
  ######################################################################

endif()

if(EXT_TOOLS AND EXT_QMC)
  
  ######################################################################
  #
  # These tests are only if tools are available
  #
  
  # VMC
  add_subdirectory(test_vmc)
  add_subdirectory(test_vmc_open_100_molecules)
  add_subdirectory(test_vmcpbc_tilted)
  add_subdirectory(test_vmc_pfaff)
  add_subdirectory(test_vmc_pbc_gamma_8_crystals_insulator)
  add_subdirectory(test_vmc_pbc_complex_8_crystals_insulator)  
  add_subdirectory(test_vmc_pbc_trim_8_crystals_insulator)

  # VMCopt
  add_subdirectory(test_vmcopt_jastrow_1-3b_linear_method)
  add_subdirectory(test_vmcopt_jastrow_1-3b_stochastic_reconfiguration)
  add_subdirectory(test_vmcopt_agp_matrix_linear_method)
  add_subdirectory(test_vmcopt_agp_matrix_stochastic_reconfiguration)
  #add_subdirectory(test_vmcopt_sd_MOs_linear_method)

  # LRDMC
  add_subdirectory(test_lrdmc)
  add_subdirectory(test_lrdmc_new)
  add_subdirectory(test_lrdmc_pfaff)
  add_subdirectory(test_lrdmc_tilted)
  add_subdirectory(test_lrdmc_dla)
  add_subdirectory(test_lrdmc_dltm)
  add_subdirectory(test_lrdmc_la)
  add_subdirectory(test_lrdmc_tmove)
 
  #DMC-Tmoves
  add_subdirectory(test_dmc_tmoves)
 
  # MD
  add_subdirectory(test_MD_classical_vmc)
  add_subdirectory(test_MD_quantum_vmc)
  
  ######################################################################

endif()
