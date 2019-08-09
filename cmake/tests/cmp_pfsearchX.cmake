###############################################################################################
# Common output format to fortan and C
# 1 replicates pfsearch -lxz output
# 2 replicates pfsearch -zkx output
# 3 replicates pfsearch   -x output
# 4 xPSA output
# 5 tsv output (single line tab delimited)
SET(COMMON_FORMAT "")
LIST(APPEND COMMON_FORMAT "-lxz")
LIST(APPEND COMMON_FORMAT "-zkx")
LIST(APPEND COMMON_FORMAT "-x")

###############################################################################################
MACRO(EXEC_RUN CMD opt format outfile)
    IF(EXISTS ${TESTS_DIRECTORY}/${outfile})
      MESSAGE("Deleting existing ${TESTS_DIRECTORY}/${outfile}")
      FILE(REMOVE ${TESTS_DIRECTORY}/${outfile})
    ENDIF(EXISTS ${TESTS_DIRECTORY}/${outfile})
    
    STRING(REPLACE ";" " " OPTIONS "${${opt}}" )
    MESSAGE("Executing ${CMD} ${format} ${OPTIONS}")
    EXECUTE_PROCESS(COMMAND "${CMD}" ${format} ${${opt}}
		                RESULT_VARIABLE CMD_RESULT
		                OUTPUT_FILE "${TESTS_DIRECTORY}/${outfile}")
    MESSAGE("Output in  ${TESTS_DIRECTORY}/${outfile}") 
    IF(CMD_RESULT)
        MESSAGE(FATAL_ERROR "Error running ${CMD} ${OPTIONS}\nCMAKE RESULTS: ${CMD_RESULT}")
    ENDIF(CMD_RESULT)
ENDMACRO(EXEC_RUN)


SET(ARGUMENTS3 "-t 1")
SET(ARGUMENTS "-f")

IF(DEFINED OPTIONS3)
	IF(NOT ${OPTIONS3} STREQUAL "")
		FOREACH(val IN ITEMS ${OPTIONS3})
			LIST(APPEND ARGUMENTS3 ${val})
		ENDFOREACH()
	ENDIF(NOT ${OPTIONS3} STREQUAL "")
ENDIF(DEFINED OPTIONS3)

IF(DEFINED OPTIONS)
	IF(NOT ${OPTIONS} STREQUAL "")
		FOREACH(val IN ITEMS ${OPTIONS})
			LIST(APPEND ARGUMENTS3 ${val})
		ENDFOREACH()
	ENDIF(NOT ${OPTIONS} STREQUAL "")
ENDIF(DEFINED OPTIONS)

IF(NOT ${PRF} STREQUAL "")
  LIST(APPEND ARGUMENTS3 "${PRF}")
  LIST(APPEND ARGUMENTS "${PRF}")
ENDIF(NOT ${PRF} STREQUAL "")

LIST(APPEND ARGUMENTS3 "${SEQDB}")
LIST(APPEND ARGUMENTS "${SEQDB}")

FOREACH(frmt_id RANGE 1 3)
	MESSAGE(" Args: -o\ ${frmt_id};${ARGUMENTS3}")
	exec_run(${EXECUTABLE3} ARGUMENTS3 "-o\ ${frmt_id}" "${OUT_FILE}_C_${frmt_id}.test")
	IF(NOT EXISTS  "${TESTS_DIRECTORY}/${OUT_FILE}_C_${frmt_id}.test")
		 MESSAGE(FATAL_ERROR "No output file ${TESTS_DIRECTORY}/${OUT_FILE}_C_${frmt_id}.test")
	ENDIF(NOT EXISTS "${TESTS_DIRECTORY}/${OUT_FILE}_C_${frmt_id}.test")
	
	MATH(EXPR index "${frmt_id}-1")
	LIST(GET COMMON_FORMAT ${index} frmt)
	MESSAGE(" Args: ${frmt};${ARGUMENTS}")
	exec_run(${EXECUTABLE} ARGUMENTS "${frmt}" "${OUT_FILE}_Fortran_${frmt_id}.test")
	IF(NOT EXISTS  "${TESTS_DIRECTORY}/${OUT_FILE}_Fortran_${frmt_id}.test")
		 MESSAGE(FATAL_ERROR "No output file ${TESTS_DIRECTORY}/${OUT_FILE}_Fortran_${frmt_id}.test")
	ENDIF(NOT EXISTS "${TESTS_DIRECTORY}/${OUT_FILE}_Fortran_${frmt_id}.test")
	
	
	EXECUTE_PROCESS(COMMAND diff "${TESTS_DIRECTORY}/${OUT_FILE}_C_${frmt_id}.test" "${TESTS_DIRECTORY}/${OUT_FILE}_Fortran_${frmt_id}.test"
									RESULT_VARIABLE CMD_RESULT
									OUTPUT_FILE "${TESTS_DIRECTORY}/${OUT_FILE}.diff"
			)
	IF(CMD_RESULT)
		MESSAGE(FATAL_ERROR "Difference found between ${TESTS_DIRECTORY}/${OUT_FILE}_C_${frmt_id}.test and ${TESTS_DIRECTORY}/${OUT_FILE}_Fortran_${frmt_id}.test")
	ELSE(CMD_RESULT)
		FILE(REMOVE "${TESTS_DIRECTORY}/${OUT_FILE}.diff")
	ENDIF(CMD_RESULT)
	
ENDFOREACH(frmt_id RANGE 1 3)

