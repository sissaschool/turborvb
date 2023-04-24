get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES prep-.* )
    set( PREFIX "" )
    if( ${EXECUTABLE} MATCHES ".*-mpi" )
      foreach(N 1;2)
        set( PREFIX "mpirun -np ${N}" )
        set( _TEST_NAME "Test DFT 38 open systems (num mpi = ${N} ${EXECUTABLE})" )
        add_test_dependency_tree(
          NAME ${_TEST_NAME}
          COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:${EXECUTABLE}> prep.out out_true.o REFERENCE_fort.10_new 5 ${PREFIX}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          DEPENDENCY_TREE ${PARENT_DIR}
          )
      endforeach()
    else()
      set( _TEST_NAME "Test DFT 38 open systems (num mpi = ${N} ${EXECUTABLE})" )
      add_test_dependency_tree(
        NAME ${_TEST_NAME}
        COMMAND ${BASH_EXECUTABLE} cm.test.sh $<TARGET_FILE:${EXECUTABLE}> prep.out out_true.o REFERENCE_fort.10_new 5 ${PREFIX}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDENCY_TREE ${PARENT_DIR}
        )
    endif()
  endif()
endforeach()
