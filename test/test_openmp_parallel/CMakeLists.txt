add_test(
  NAME "Test OpenMP"
  COMMAND $<TARGET_FILE:test_openmp_parallel>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
