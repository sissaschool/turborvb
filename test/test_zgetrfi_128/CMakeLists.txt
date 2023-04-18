add_test(
  NAME "Test zgetrfi_ wrapper (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgetrfi> test_zgetrfi.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
