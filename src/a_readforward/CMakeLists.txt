foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES readforward-.* )

    set_target_properties( ${EXECUTABLE} PROPERTIES OUTPUT_NAME ${EXECUTABLE} )
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${EXECUTABLE} )

    target_sources( ${EXECUTABLE}
      PRIVATE
        CompPW.f90
        CompRhok.f90
        compute_corr_fun.f90
        cossin.f90
        dsort.f90
        el_ion_distance_2d.f90
        el_ion_distance_for_dipole.f90
        el_ion_distance.f90
        read_corr_fun.f90
        readforward.f90
        shells.f90
        write_corr_fun.f90
        berry_phase.f90
        compute_dipole.f90
        compute_grid.f90
        compute_qpwf.f90
        compute_rho_assar.f90
        compute_rho.f90
        compute_spin2.f90
    )

    string(REGEX REPLACE "^readforward" "common" COMMON_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^readforward" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^readforward" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^readforward" "complex" COMPLEX_LIBRARY_TO_LINK ${EXECUTABLE})
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${COMPLEX_LIBRARY_TO_LINK}
    )
    
    target_include_directories( ${EXECUTABLE}
      PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/b_complex
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/c_adjoint_backward
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/c_adjoint_forward
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/d_qlapack
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/m_common
        ${CMAKE_CURRENT_LIST_DIR}/CMakeFiles/${EXECUTABLE}.dir/src/m_pfapack
      )

  endif()
endforeach()
