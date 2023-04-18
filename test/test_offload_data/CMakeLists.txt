add_test(
    NAME "Test offload data (real)"
  COMMAND $<TARGET_FILE:test_offload_data_real>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_test(
    NAME "Test offload data (complex)"
  COMMAND $<TARGET_FILE:test_offload_data_complex>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
