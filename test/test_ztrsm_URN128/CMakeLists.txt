add_test(
  NAME "Test ztrsm_ wrapper URN (128 x 128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_ztrsm> test_ztrsm.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
