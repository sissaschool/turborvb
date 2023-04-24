# Explenation for these command can be found in the root CMakeLists.txt file

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES prep-.* )

    set_target_properties( ${EXECUTABLE} PROPERTIES OUTPUT_NAME ${EXECUTABLE} )
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${EXECUTABLE} )

    target_sources( ${EXECUTABLE}
      PRIVATE
        hamiltonian.f90
        initialize_matrices.f90
        writeout_finalize.f90
        upocc_kpoints.f90
        graham_scalapack.f90
        eval_hamilt.f90
        fillbuffer.f90
        eval_hamilt_complex.f90
        invsymn.f90
        prep.f90
        lda.f90
        lsda.f90
        elec_density.f90
        printden.f90
        update_fort10.f90
        fill_matrix_scalapack.f90
        improve_diag.f90
        real_complex_routines.f90
        compute_bands.f90
        initialize_environment.f90
        buffers.f90
        parallel_mod.f90
        fourier_mod.f90
        freeelmod_complex_kpoints.f90
        eval_efermi.f90
        setup.f90
        )

    string(REGEX REPLACE "^prep" "common" COMMON_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^prep" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^prep" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^prep" "complex" COMPLEX_LIBRARY_TO_LINK ${EXECUTABLE})
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${COMPLEX_LIBRARY_TO_LINK}
    )
    
  endif()
endforeach()

