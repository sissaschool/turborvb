add_test(
  NAME "Test zdotc_ wrapper (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zdotc> test_zdotc.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
