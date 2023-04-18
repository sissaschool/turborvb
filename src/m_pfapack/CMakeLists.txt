foreach( LIBRARY IN LISTS LIBRARIES_S_L
                          LIBRARIES_P_L )

if( ${LIBRARY} MATCHES pfapack )

    set_module_path_directory( ${LIBRARY} )
    
    target_sources( ${LIBRARY}
        PUBLIC
        dsktrf.f
        dlasktrf.f
        dsktf2.f
        dskr2.f
        dskr2k.f
        zsktrf.f
        zlasktrf.f
        zsktf2.f
        zskr2.f
        zskr2k.f
        zsktrd.f
        dsktrd.f
        dlasktrd.f
        dskmv.f
        dsktd2.f
        zlasktrd.f
        zskmv.f
        zsktd2.f
    )
    
    target_include_directories( ${LIBRARY}
    PUBLIC
        ${CMAKE_CURRENT_LIST_DIR}
    )

endif()
endforeach()
