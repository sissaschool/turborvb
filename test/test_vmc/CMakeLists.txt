get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} NAME)

foreach( EXECUTABLE IN LISTS EXECUTABLES_S_L
                             EXECUTABLES_P_L )
  if( ${EXECUTABLE} MATCHES turborvb-.*)
  if( ${EXECUTABLE} MATCHES ".*-serial" )
     set( PREFIX "" )
     add_test_dependency_tree(
      NAME "Test VMC run (${EXECUTABLE})"
      COMMAND ${BASH_EXECUTABLE} cm.test1.sh TURBORVB=$<TARGET_FILE:${EXECUTABLE}> 
                                             READF=$<TARGET_FILE:readf> 
                                             INPUT_TURBORVB=datasvmc.d
                                             OUTPUT_TURBORVB=serial.out
                                             OUTPUT_REF=out_true.o
                                             FORT21_REF=REFERENCE_fortXXI
                                             ROUND_OFF=6
                                             EXECUTE_PREFIX=${PREFIX}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDENCY_TREE ${PARENT_DIR}
      )
  endif()
  if( ${EXECUTABLE} MATCHES ".*-mpi" )
     set( PREFIX "mpirun -np 2" )
     add_test_dependency_tree(
      NAME "Test VMC run (${EXECUTABLE})"
      COMMAND ${BASH_EXECUTABLE} cm.test1.sh TURBORVB=$<TARGET_FILE:${EXECUTABLE}>
                                             READF=$<TARGET_FILE:readf>
                                             INPUT_TURBORVB=datasvmc.d
                                             OUTPUT_TURBORVB=mpi_np2.out
                                             OUTPUT_REF=out_true_mpi_np2.o
                                             FORT21_REF=REFERENCE_fortXXI_mpi_np2
                                             ROUND_OFF=6
                                             EXECUTE_PREFIX=${PREFIX}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      DEPENDENCY_TREE ${PARENT_DIR}
      )
  endif()
  if( ${EXT_GPU} )
      if( ${EXECUTABLE} MATCHES ".*-serial" )
         set( PREFIX "" )
         add_test_dependency_tree(
          NAME "Test VMC run, forcing yes_ontarget (${EXECUTABLE})"
          COMMAND ${BASH_EXECUTABLE} cm.test1.sh TURBORVB=$<TARGET_FILE:${EXECUTABLE}> 
                                                 READF=$<TARGET_FILE:readf> 
                                                 INPUT_TURBORVB=datasvmc_yot.d
                                                 OUTPUT_TURBORVB=serial_yot.out
                                                 OUTPUT_REF=out_true.o
                                                 FORT21_REF=REFERENCE_fortXXI
                                                 ROUND_OFF=6
                                                 EXECUTE_PREFIX=${PREFIX}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          DEPENDENCY_TREE ${PARENT_DIR}
          )
      endif()
      if( ${EXECUTABLE} MATCHES ".*-mpi" )
         set( PREFIX "mpirun -np 2" )
         add_test_dependency_tree(
          NAME "Test VMC run, forcing yes_ontarget (${EXECUTABLE})"
          COMMAND ${BASH_EXECUTABLE} cm.test1.sh TURBORVB=$<TARGET_FILE:${EXECUTABLE}>
                                                 READF=$<TARGET_FILE:readf>
                                                 INPUT_TURBORVB=datasvmc_yot.d
                                                 OUTPUT_TURBORVB=mpi_np2_yot.out
                                                 OUTPUT_REF=out_true_mpi_np2.o
                                                 FORT21_REF=REFERENCE_fortXXI_mpi_np2
                                                 ROUND_OFF=6
                                                 EXECUTE_PREFIX=${PREFIX}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          DEPENDENCY_TREE ${PARENT_DIR}
          )
      endif()
  endif()
  endif()

  if( ${EXECUTABLE} MATCHES testad*)
    add_test_dependency_tree(
    	NAME "Test VMC (AD) run"
    	COMMAND ${BASH_EXECUTABLE} cm.test_ad.sh $<TARGET_FILE:testad>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDENCY_TREE ${PARENT_DIR}
      )
  endif()

endforeach()

