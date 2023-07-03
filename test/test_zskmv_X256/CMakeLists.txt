get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test_dependency_tree(
  NAME "Test zskmv_ wrapper U (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zskmv> U test_zskmv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDENCY_TREE ${PARENT_DIR}
  )

add_test_dependency_tree(
  NAME "Test zskmv_ wrapper L (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zskmv> L test_zskmv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDENCY_TREE ${PARENT_DIR}
  )
