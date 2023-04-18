add_test(
  NAME "Test offload pointer transcription"
  COMMAND $<TARGET_FILE:test_offload_pointer_transcription>
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
