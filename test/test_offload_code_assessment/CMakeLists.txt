add_test(
    NAME "Test offload code assessment (real)"
  COMMAND $<TARGET_FILE:test_offload_code_assessment_real>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_test(
    NAME "Test offload code assessment (complex)"
  COMMAND $<TARGET_FILE:test_offload_code_assessment_complex>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
