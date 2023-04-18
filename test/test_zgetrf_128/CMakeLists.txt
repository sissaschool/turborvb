add_test(
  NAME "Test zgetrf_ wrapper (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgetrf> test_zgetrf.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
