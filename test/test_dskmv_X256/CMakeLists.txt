get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test_dependency_tree(
  NAME "Test dskmv_ wrapper U (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dskmv> U test_dskmv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDENCY_TREE ${PARENT_DIR}
  )

add_test_dependency_tree(
  NAME "Test dskmv_ wrapper L (256)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_dskmv> L test_dskmv.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDENCY_TREE ${PARENT_DIR}
  )
