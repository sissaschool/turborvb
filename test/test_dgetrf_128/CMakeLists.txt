add_test(
  NAME "Test dgetrf_ wrapper (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgetrf> test_dgetrf.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
