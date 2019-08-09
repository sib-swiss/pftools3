FILE(GLOB PerlScripts *.pl)
SET(PERL_SYNTAX_OK TRUE)
MESSAGE(STATUS "Check Perl script syntax")
FOREACH(target ${PerlScripts})
	EXECUTE_PROCESS(COMMAND echo -n "       ")
	EXECUTE_PROCESS(COMMAND basename ${target})
	EXECUTE_PROCESS(COMMAND perl -c INPUT_FILE ${target} RESULT_VARIABLE
			retcode ERROR_FILE stderr)
	IF(NOT "${retcode}" STREQUAL "0")
		EXECUTE_PROCESS(COMMAND cat stderr)
		#MESSAGE(FATAL_ERROR "")
		SET(PERL_SYNTAX_OK FALSE)
	ELSE()
		EXECUTE_PROCESS(COMMAND echo -n "           ")
		EXECUTE_PROCESS(COMMAND cat stderr)
	ENDIF()
ENDFOREACH()


#########################################################################################
# Installation
IF(PERL_SYNTAX_OK)
	INSTALL(FILES ${PerlScripts}
					DESTINATION bin
					COMPONENT Applications
					PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
	)
ELSE (PERL_SYNTAX_OK)
	MESSAGE(STATUS "Perl scripts contains error, not installing them")
ENDIF(PERL_SYNTAX_OK)
