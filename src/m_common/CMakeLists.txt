foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )

if( ${LIBRARY} MATCHES common-.* )

    set_module_path_directory( ${LIBRARY} )

    # replace path info. in _help_online.f90
    set(HELP_FILE ${CMAKE_CURRENT_BINARY_DIR}/help_online.f90)
    set(ORG_HELP_FILE ${CMAKE_CURRENT_SOURCE_DIR}/_help_online.f90)
    add_custom_command(
        OUTPUT ${HELP_FILE}
        COMMAND sed "s%NAME_DIR%${CMAKE_SOURCE_DIR}%" ${ORG_HELP_FILE} > ${HELP_FILE}
        DEPENDS ${ORG_HELP_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    # replace version info. in _version.f90
    set(VER_FILE ${CMAKE_CURRENT_BINARY_DIR}/version.f90)
    set(ORG_VER_FILE ${CMAKE_CURRENT_SOURCE_DIR}/_version.f90)
    add_custom_command(
        OUTPUT ${VER_FILE}
        COMMAND sed "s%unknown%${SHORT_HASH}%" ${ORG_VER_FILE} > ${VER_FILE}
        DEPENDS ${ORG_VER_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    
    if( ${LIBRARY} MATCHES ".*-serial" )
        if( NOT TARGET help_file_target_serial )
            add_custom_target(help_file_target_serial DEPENDS ${HELP_FILE})
        endif()
        if( NOT TARGET ver_file_target_serial )
            add_custom_target(ver_file_target_serial DEPENDS ${VER_FILE})
        endif()
        add_dependencies(${LIBRARY} help_file_target_serial)
        add_dependencies(${LIBRARY} ver_file_target_serial)
    endif()
    if( ${LIBRARY} MATCHES ".*-mpi" )
        if( NOT TARGET help_file_target_mpi )
            add_custom_target(help_file_target_mpi DEPENDS ${HELP_FILE})
        endif()
        if( NOT TARGET ver_file_target_mpi )
            add_custom_target(ver_file_target_mpi DEPENDS ${VER_FILE})
        endif()
        add_dependencies(${LIBRARY} help_file_target_mpi)
        add_dependencies(${LIBRARY} ver_file_target_mpi)
    endif()
    
    target_sources( ${LIBRARY}
        PRIVATE
        # modules
        allio.f90
        buffer.f90
        cell.f90
        constants.f90
        damping.f90
        dielectric.f90
        evaluate_invariant.f90
        ewald.f90
        kpoints.f90
        mod_extpot.f90
        mod_IO.f90
        mpiio.f90
        nvtx.f90
        qmckl.f90
        sub_comm.f90
        symmetries.f90
        trexio.f90
        types.f90
        # subroutines
        atom_weight.f90
        bconstraint.f90
        bconstrbra.f90
        bconstrbr_complex.f90
        bconstrbr.f90
        bspline90_22.f90
        constrbra_complex.f90
        constrbra.f90
        constrbr_complex.f90
        constrbr.f90
        convertdec.f90
        copy_eagp.f90
        definition.f90
        descriptors.f90
        dgemm_my.f90
        dscalzero.f90
        dsktri.f90
        dsktrs.f90
        dsortx.f90
        error.f90
        extv.f90
        fft_scalar.f90
        fft_stick.c
        fftw.c
        fileOP.f90
        fillmatrix.f90
        findmaxmat.f90
        findrionref.f90
        forces_ext.f90
        fortran.c
        fort11_io.f90
        graham.f90
        grahamo.f90
        checkmatrix.f90
        iesdr1iesd.f90
        iesdr2iesd.f90
        initconf_pbc.f90
        invsymeps.f90
        io.c
        ioptorbcontr.f90
        kind.f90
        max_ovlp.f90
        memOP.f90
        molec_pfaff.f90
        pareff.f90
        plot_3d_data.f90
        plot_3d_data_tilted.f90
        print_eigenvalues.f90
        random.f90
        randomnumber.c
        ran.f90
        read_pseudo.f90
        rotate_tools.f90
        save_jall_fn.f90
        scalevect.f90
        slaterorb.f90
        symmetrize_agp.f90
        Thomas_Fermi_model.f90
        updatedwarp.f90
        update_jastrowall.f90
        upsim.f90
        upsimp.f90
        upvpotaa.f90
        upvpot_ei.f90
        upwinv.f90
        upwinvp.f90
        write_type_orb.f90
        zgemm_my.f90
        zsktri.f90
        zsktrs.f90
        zsktrs_qp.f90
        ${HELP_FILE}
        ${VER_FILE}
    )

    #target_include_directories( ${LIBRARY}
    #  PRIVATE
    #    ${CMAKE_CURRENT_LIST_DIR}
    #)

    if( ${LIBRARY} MATCHES .*-nogpu-.* )
      set_target_properties( ${LIBRARY} PROPERTIES NOFLAG "_OFFLOAD;_CUSOLVER;_CUBLAS" )
    endif()

endif()
endforeach()
