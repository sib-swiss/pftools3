#########################################################################################
## Documentation

SET(Documents
    profile.pdf
    profile.txt
    poster/posterPfsearchV3SIBdays2013_v0.3.pdf
    Article/btt129.pdf
)

INSTALL(FILES  ${Documents}
	DESTINATION ${DOC_PATH}
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT Documentation
)
