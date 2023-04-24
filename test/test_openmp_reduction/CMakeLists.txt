add_test(
  NAME "Test OpenMP reduction (GPU if needed, real)"
  COMMAND $<TARGET_FILE:test_openmp_reduction_real>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_test(
  NAME "Test OpenMP reduction (GPU if needed, complex)"
  COMMAND $<TARGET_FILE:test_openmp_reduction_complex>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
