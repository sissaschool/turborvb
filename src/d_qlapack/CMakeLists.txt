foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )
  if( ${LIBRARY} MATCHES qlapack-.* )

    set_module_path_directory( ${LIBRARY} )

    target_sources( ${LIBRARY}
        PRIVATE
           dgemm_tn.f
           dgemv.f
           dgemvn.f
           dgemvt.f
           dger.f
           dger2.f
           getrfi.f
           zgemm_tn.f
           zgemv.f
           zgemvn.f
           zgemvt.f
           zgeru.f
           zger2.f
           zdotc.f
    )
    
  string(REGEX REPLACE "^qlapack" "common" COMMON_LIBRARY_TO_LINK ${LIBRARY})

    target_link_libraries( ${LIBRARY}
        PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        )

    if( ${LIBRARY} MATCHES .*-nogpu-.* )
      set_target_properties( ${LIBRARY} PROPERTIES NOFLAG "_OFFLOAD;_CUSOLVER;_CUBLAS" )
    endif()

  endif()
endforeach()
