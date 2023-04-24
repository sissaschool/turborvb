foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )
  if( ${LIBRARY} MATCHES adjoint-forward-.* )

    set_module_path_directory( ${LIBRARY} )

    target_sources( ${LIBRARY}
      PRIVATE
        #adjoint forward
        compute_fast.f90
        eval_iond.f90
        jastrow_ee.f90
        jastrow_ei.f90
        jastrowgrad.f90
        jastrowgrad_pbc.f90
        legfun.f90
        makefun.f90
        makefun_bump.f90
        makefun_pbc.f90
        makefun0.f90
        makefun0_bump.f90
        makefun0_pbc.f90
        pseudofun.f90
        pseudoset.f90
        subener.f90
        t_lrdmc.f90
        up2bodypsi.f90
        upvpotdiag.f90
        upwf.f90
        ngivej.f90
        #else
        dgesvd.f90
        dsyev.f90
        evaldet_complex.f90
        evaldet.f90
        ratiovar.f90
        upinvhop_fnf_new.f90
        upinvhop.f90
        dspev_drv.f90
        zgesvd.f90
        zsyev.f90
        )

    string(REGEX REPLACE "^adjoint-forward" "common" COMMON_LIBRARY_TO_LINK ${LIBRARY})
    string(REGEX REPLACE "^adjoint-forward" "qlapack" QLAPACK_LIBRARY_TO_LINK ${LIBRARY})

    target_link_libraries( ${LIBRARY}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${QLAPACK_LIBRARY_TO_LINK}
        pfapack
        )
    
  endif()
endforeach()
