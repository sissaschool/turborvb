add_test(
  NAME "Test dger2 wrapper (256 x 256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dger2> test_dger2.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
