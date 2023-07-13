# Explenation for these command can be found in the root CMakeLists.txt file

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES turborvb-.* )

    set_target_properties( ${EXECUTABLE} PROPERTIES OUTPUT_NAME ${EXECUTABLE} )
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${EXECUTABLE} )

    target_sources( ${EXECUTABLE}
      PRIVATE
        adam_opt.f90
        add_sincos.f90
        adjust_tpar.f90
        atom_weight.f90
        bootforcecov.f90
        bootparameter.f90
        bootpress.f90
        branchingo.f90
        comp_econf.f90
        conjginv.f90
        conjginv_prep.f90
        conjginvs.f90
        constraint_complex.f90
        cutminzmaxz.f90
        dgelscut_new.f90
        diffus.f90
        dsyev_sing_sz.f90
        enercont.f90
        enercut.f90
        enercutoff.f90
        energy.f90
        enforce_detail.f90
        eval_pulay.f90
        evalgren.f90
        gather_avgs.f90
        gemm.f90
        hopping_pbc.f90
        jastrow_exp.f90
        main.f90
        makeu.f90
        mapping_sz.f90
        periodic.f90
        ratio_psi.f90
        ratiofn_psi.f90
        reshuffhub_sz_new.f90
        reweight0_k6.f90
        root2mat.f90
        rotation.f90
        scal_lins.f90
        scratchdet.f90
        set_turboq.f90
        transpsip.f90
        transpsip_complex.f90
        upcellkel.f90
        updiag.f90
        upgradcont.f90
        upjbra.f90
        upkinc.f90
        uptable.f90
        uptable_complex.f90
        uptabpip_new.f90
        uptabtot_new.f90
        upvpot_pbc_exp.f90
        wf_sign.f90
        writeoutput.f90
        )

    string(REGEX REPLACE "^turborvb" "common" COMMON_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^turborvb" "qlapack" QLAPACK_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^turborvb" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^turborvb" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^turborvb" "complex" COMPLEX_LIBRARY_TO_LINK ${EXECUTABLE})

    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${QLAPACK_LIBRARY_TO_LINK}
        pfapack
        ${ADBACKWARD_LIBRARY_TO_LINK}
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${COMPLEX_LIBRARY_TO_LINK}
    )

  endif()
endforeach()

