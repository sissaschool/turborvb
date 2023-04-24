add_test(
  NAME "Test zgemv_ wrapper T (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgemv> test_zgemv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
if(EXT_GPU)
  add_test(
    NAME "Test zgemv_offload subroutine T (256)"
    COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zgemv_offload> test_zgemv_offload.out
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
