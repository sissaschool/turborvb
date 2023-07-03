get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES turborvb-.*)
  if( ${EXECUTABLE} MATCHES ".*-serial" )
     set( PREFIX "" )
     add_test_dependency_tree(
      NAME "Test MD classical-VMC run (${EXECUTABLE})"
	  COMMAND ${BASH_EXECUTABLE} cm.test1.sh $<TARGET_FILE:${EXECUTABLE}> $<TARGET_FILE:readf> NA serial.out out_true.o REFERENCE_fortXXI 6 ${PREFIX}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDENCY_TREE ${PARENT_DIR}
      )
  endif()
  if( ${EXECUTABLE} MATCHES ".*-mpi" )
     set( PREFIX "mpirun -np 2" )
     add_test_dependency_tree(
      NAME "Test MD classical-VMC run (${EXECUTABLE})"
	  COMMAND ${BASH_EXECUTABLE} cm.test1.sh $<TARGET_FILE:${EXECUTABLE}> $<TARGET_FILE:readf> NA mpi_np2.out out_true_mpi_np2.o REFERENCE_fortXXI_mpi_np2 6 "${PREFIX}"
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDENCY_TREE ${PARENT_DIR}
      )
  endif()
  endif()
endforeach()
