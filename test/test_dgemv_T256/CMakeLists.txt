add_test(
  NAME "Test dgemv_ wrapper T (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgemv> test_dgemv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
add_test(
  NAME "Test dgemv__ wrapper T (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgemv_> test_dgemv_.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  )
if(EXT_GPU)
  add_test(
    NAME "Test dgemv_offload subroutine T (256)"
    COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dgemv_offload> test_dgemv_offload.out
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
