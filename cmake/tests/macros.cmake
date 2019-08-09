MACRO(TEST_PFSEARCHV3 subdir name prf seqdb options output_format output_filename)
	IF (NOT EXISTS ${TESTS_DIRECTORY}/${subdir}/)
		FILE(MAKE_DIRECTORY ${TESTS_DIRECTORY}/${subdir}/)
	ENDIF(NOT EXISTS ${TESTS_DIRECTORY}/${subdir}/)

	ADD_TEST(NAME pfsearchV3_${subdir}_${name}
		WORKING_DIRECTORY "${TESTS_DIRECTORY}"
		COMMAND ${CMAKE_COMMAND}
					-DEXECUTABLE3=$<TARGET_FILE:pfsearchV3>
					-DTESTS_DIRECTORY=${TESTS_DIRECTORY}
					-DRESULTS_DIR=${TESTS_DATA_DIR}/cmake/Truth
					-DPRF=${TESTS_DATA_DIR}/cmake/${subdir}/${prf}
					-DSEQDB=${TESTS_DATA_DIR}/cmake/${subdir}/${seqdb}
					-DOUTPUT_FORMAT=${output_format}
					-DOPTIONS=${options}
					-DOUT_FILE=${subdir}/${output_filename}.test
					-P ${PROJECT_SOURCE_DIR}/cmake/tests/run_pfsearch_test.cmake
	)
ENDMACRO(TEST_PFSEARCHV3)

MACRO(COMPARE_PFSEARCH subdir name prf seqdb opt3 opt)
	IF (NOT EXISTS ${TESTS_DIRECTORY}/${subdir}/)
		FILE(MAKE_DIRECTORY ${TESTS_DIRECTORY}/${subdir}/)
	ENDIF(NOT EXISTS ${TESTS_DIRECTORY}/${subdir}/)

	ADD_TEST(NAME compare_${subdir}_${name}
		WORKING_DIRECTORY "${TESTS_DIRECTORY}"
		COMMAND ${CMAKE_COMMAND}
					-DEXECUTABLE3=$<TARGET_FILE:pfsearchV3>
					-DEXECUTABLE=$<TARGET_FILE:pfsearch>
					-DTESTS_DIRECTORY=${TESTS_DIRECTORY}
					-DRESULTS_DIR=${TESTS_DATA_DIR}/cmake/Truth
					-DPRF=${TESTS_DATA_DIR}/cmake/${subdir}/${prf}
					-DSEQDB=${TESTS_DATA_DIR}/cmake/${subdir}/${seqdb}
					-DOPTIONS3=${opt3}
					-DOUT_FILE=${subdir}/${name}
					-P ${PROJECT_SOURCE_DIR}/cmake/tests/cmp_pfsearchX.cmake
	)
ENDMACRO(COMPARE_PFSEARCH)
