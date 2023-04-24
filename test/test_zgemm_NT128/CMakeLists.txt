add_test(
  NAME "Test zgemm_ wrapper NT (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgemm> test_zgemm.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
