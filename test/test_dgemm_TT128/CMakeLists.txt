add_test(
  NAME "Test dgemm_ wrapper TT (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgemm> test_dgemm.out 6
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
