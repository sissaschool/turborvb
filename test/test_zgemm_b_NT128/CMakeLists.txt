add_test(
  NAME "Test zgemm_b NN (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgemm_b> test_zgemm_b.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
