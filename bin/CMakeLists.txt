#file( READ turbo-genius.sh turbogenius )
#file( READ forcevm.sh turbogenius )

install( FILES forcevmc.sh
	 PERMISSIONS OWNER_READ OWNER_EXECUTE GROUP_READ WORLD_READ
         DESTINATION bin )

install( FILES check_fort12_columns_length.py
	 PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
         DESTINATION bin )

install( FILES read_columns_fort21.py
	 PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
         DESTINATION bin )
