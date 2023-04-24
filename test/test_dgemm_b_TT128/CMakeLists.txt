add_test(
  NAME "Test dgemm_b TT (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgemm_b> test_dgemm_b.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
