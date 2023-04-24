foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )

  if( ${LIBRARY} MATCHES complex-.* )
    target_sources( ${LIBRARY}
      PRIVATE
        fort10_io.f90
        convertmod.f90
        eval_molec_unpaired.f90
        eval_molec_epsdgel.f90
        read_datas.f90
        )

    set_module_path_directory( ${LIBRARY} )

    string(REGEX REPLACE "^complex" "common" COMMON_LIBRARY_TO_LINK ${LIBRARY})
    string(REGEX REPLACE "^complex" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${LIBRARY})
    string(REGEX REPLACE "^complex" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${LIBRARY})

    target_link_libraries( ${LIBRARY}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        qlapack
        pfapack
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${ADBACKWARD_LIBRARY_TO_LINK}
        )
    
    if( ${LIBRARY} MATCHES .*-nogpu-.* )
      set_target_properties( ${LIBRARY} PROPERTIES NOFLAG "_OFFLOAD;_CUSOLVER;_CUBLAS" )
    endif()

  endif()
endforeach()
