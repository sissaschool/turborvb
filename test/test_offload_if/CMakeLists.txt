add_test(
  NAME "Test offload if"
  COMMAND $<TARGET_FILE:test_offload_if>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
