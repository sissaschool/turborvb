add_test(
  NAME "Test zgeru_ wrapper (256 x 256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgeru> test_zgeru.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
