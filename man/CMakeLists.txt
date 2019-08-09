#########################################################################################
## man files

# Find all the input files
FILE(GLOB ManFiles1 RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/*.1")

INSTALL(FILES  ${ManFiles1}
	DESTINATION ${MAN1_PATH}
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT Documentation
)


FILE(GLOB ManFiles5 RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
	"${CMAKE_CURRENT_SOURCE_DIR}/*.5")

INSTALL(FILES  ${ManFiles5}
	DESTINATION ${MAN5_PATH}
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT Documentation
)

