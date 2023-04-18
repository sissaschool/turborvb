foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES convertfort10mol-.* )

    set_target_properties( ${EXECUTABLE} PROPERTIES OUTPUT_NAME ${EXECUTABLE} )
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${EXECUTABLE} )

    target_sources( ${EXECUTABLE}
      PRIVATE
        convertfort10mol.f90
    )

    string(REGEX REPLACE "^convertfort10mol" "common" COMMON_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10mol" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10mol" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10mol" "complex" COMPLEX_LIBRARY_TO_LINK ${EXECUTABLE})
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${COMPLEX_LIBRARY_TO_LINK}
    )

  endif()
endforeach()
