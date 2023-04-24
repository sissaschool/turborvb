if(makefort10 IN_LIST EXECUTABLES_S_L)

  set_module_path_directory( makefort10 serial )
  
  target_sources( makefort10
    PRIVATE
      kind.f90
      constants2.f90
      symm_data.f90
      orbitals.f90
      findsection.f90
      makelambda.f90
      makeforces.f90
      cubicsym.f90
      invmat.f90
      intc.f90
      makepfaff.f90
      makefort10.f90
    )
  
  set_target_properties( makefort10 PROPERTIES SUFFIX ".x")

  target_link_libraries( makefort10
    PRIVATE
      common-serial
      adjoint-forward-serial
      complex-serial
    )

endif()
