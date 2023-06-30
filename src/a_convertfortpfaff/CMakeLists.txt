if(convertpfaff IN_LIST EXECUTABLES_S_L)

  set_module_path_directory( convertpfaff serial )
  
  target_sources( convertpfaff
    PRIVATE
      convertpfaff.f90
    )
  
  set_target_properties( convertpfaff PROPERTIES SUFFIX ".x")

  target_link_libraries( convertpfaff
    PRIVATE
      common-serial
      adjoint-forward-serial
      complex-serial
    )

endif()