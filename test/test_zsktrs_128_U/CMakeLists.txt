get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_test_dependency_tree(
  NAME "Test zsktrs wrapper U (128)"
  COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:test_zsktrs> U test_zsktrs.out
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDENCY_TREE ${PARENT_DIR}
  )
