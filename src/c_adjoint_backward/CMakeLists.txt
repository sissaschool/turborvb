foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )
  if( ${LIBRARY} MATCHES adjoint-backward-.* )

    set_module_path_directory( ${LIBRARY} )

    target_sources( ${LIBRARY}
      PRIVATE
        adBuffer.f90              # only in backward
        compute_fast_b.f90        #
        ddot_b.f90                # only in backward
        derfc_b.f90               # only in backward
        dgemm_b.f90               # only in backward
        dgemv_b.f90               # only in backward
        eval_iond_b.f90           #
        jastrow_ee_b.f90          #
        jastrow_ei_b.f90          #
        jastrowgrad_b.f90         #
        jastrowgrad_pbc_b.f90     #
        legfun_b.f90              #
        makefun_b.f90             #
        makefun_bump_b.f90        #
        makefun_pbc_b.f90         #
        makefun0_b.f90            #
        makefun0_bump_b.f90       #
        makefun0_pbc_b.f90        #
        ngivej_b.f90              # only in backward
        pseudofun_b.f90           #
        pseudoset_b.f90           #
        subener_b.f90             #
        t_lrdmc_b.f90             #
        up2bodypsi_b.f90          #
        upvpotdiag_b.f90          #
        upwf_b.f90                #
        zgemm_b.f90               # only in backward
        )
    
    string(REGEX REPLACE "^adjoint-backward" "common" COMMON_LIBRARY_TO_LINK ${LIBRARY})
    string(REGEX REPLACE "^adjoint-backward" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${LIBRARY})

    target_link_libraries( ${LIBRARY}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        qlapack
        pfapack
        ${ADFORWARD_LIBRARY_TO_LINK}
        )
    
  endif()
endforeach()
