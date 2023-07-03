foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES convertfort10-.* )

    set_target_properties( ${EXECUTABLE} PROPERTIES OUTPUT_NAME ${EXECUTABLE} )
    set_target_properties( ${EXECUTABLE} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${EXECUTABLE} )

    target_sources( ${EXECUTABLE}
      PRIVATE
        convertfort10.f90
    )

    string(REGEX REPLACE "^convertfort10" "common" COMMON_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10" "adjoint-forward" ADFORWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10" "adjoint-backward" ADBACKWARD_LIBRARY_TO_LINK ${EXECUTABLE})
    string(REGEX REPLACE "^convertfort10" "complex" COMPLEX_LIBRARY_TO_LINK ${EXECUTABLE})
    target_link_libraries( ${EXECUTABLE}
      PRIVATE
        ${COMMON_LIBRARY_TO_LINK}
        ${ADFORWARD_LIBRARY_TO_LINK}
        ${COMPLEX_LIBRARY_TO_LINK}
    )

  endif()
endforeach()
