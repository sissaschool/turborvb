
  foreach(TOOL_NAME ${TOOLS_LIST_SIMPLE})

    set_target_properties( ${TOOL_NAME} PROPERTIES SUFFIX ".x")

    set_module_path_directory( ${TOOL_NAME} )
    
    if( ${TOOL_NAME} MATCHES "assembling_pseudo" )
    # replace path info. in _assembling_pseudo.f90
    set(PP_FILE ${CMAKE_CURRENT_BINARY_DIR}/assembling_pseudo.f90)
    set(ORG_PP_FILE ${CMAKE_CURRENT_SOURCE_DIR}/_assembling_pseudo.f90)
    add_custom_command(
        OUTPUT ${PP_FILE}
        COMMAND sed "s%NAME_DIR%${CMAKE_SOURCE_DIR}%" ${ORG_PP_FILE} > ${PP_FILE}
        DEPENDS ${ORG_PP_FILE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_custom_target(pp_file_target DEPENDS ${PP_FILE})
    add_dependencies(${TOOL_NAME} pp_file_target)
    target_sources( ${TOOL_NAME}
      PRIVATE
        ${PP_FILE}
      )
    else()
    target_sources( ${TOOL_NAME}
      PRIVATE
        "${TOOL_NAME}.f90"
      )
    endif()
      
    if(EXT_GPU)
        target_link_libraries( ${TOOL_NAME}
         PRIVATE
         complex-nogpu-serial
        )
    else()
        target_link_libraries( ${TOOL_NAME}
         PRIVATE
         complex-serial
        )
    endif()
    install( TARGETS ${TOOL_NAME} )

  endforeach()
