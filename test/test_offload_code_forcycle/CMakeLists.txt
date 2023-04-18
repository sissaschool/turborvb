add_test(
    NAME "Test offload code forcycle (real)"
  COMMAND $<TARGET_FILE:test_offload_code_forcycle_real>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_test(
    NAME "Test offload code forcycle (complex)"
  COMMAND $<TARGET_FILE:test_offload_code_forcycle_complex>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
