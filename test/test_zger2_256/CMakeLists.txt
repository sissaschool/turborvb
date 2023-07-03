add_test(
  NAME "Test zger2 wrapper (256 x 256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zger2> test_zger2.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
