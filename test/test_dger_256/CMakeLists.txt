add_test(
  NAME "Test dger_ wrapper (256 x 256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dger> test_dger.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
