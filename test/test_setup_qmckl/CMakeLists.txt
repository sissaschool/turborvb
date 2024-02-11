get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES turborvb-serial)
  add_test(
    NAME "test qmckl setup (${EXECUTABLE})"
    COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:${EXECUTABLE}> $<TARGET_FILE:readf> datasvmc.d out_vmc fort.21_true 6
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
  endif()

endforeach()
