add_test(
  NAME "Test dtrsm_ wrapper ULN (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dtrsm> test_dtrsm.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
