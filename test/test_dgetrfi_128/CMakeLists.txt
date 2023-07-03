add_test(
  NAME "Test dgetrfi_ wrapper (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgetrfi> test_dgetrfi.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
