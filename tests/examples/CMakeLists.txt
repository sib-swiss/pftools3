ENABLE_TESTING()
ADD_TEST(
		NAME "execute_test_V2.sh"
		WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
		COMMAND
			sh -ve test_V2.sh
)
ADD_TEST(
	NAME "execute_test_V3.sh"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND
		"sh" "-ve" "test_V3.sh"
)
ADD_TEST(
	NAME "check_output_of_test_V3.sh"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND
		"sh" "-e" "test_V3_output.sh"
)
ADD_TEST(
	NAME "execute_test_scan_search.pl"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND "perl" "test_scan_search.pl" "all.prf" "all.seq" "/tmp/test_scan_search"
)
ADD_TEST(
	NAME "execute_test_pfsearchV3_iupac.pl"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND "perl" "test_pfsearchV3_iupac.pl" "/tmp/test_pfsearchV3_iupac"
)

