ENABLE_TESTING()
ADD_TEST(
		NAME "execute_test_V2.sh"
		WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
		COMMAND
			sh -ve test_V2.sh
#		[ 2>&1 | perl -ne '$ok=1 if /#!/;print $_ if $ok' ] # see comment below
#        	> "${TESTS_DIRECTORY}/examples/test_V3.tmp"
)

#ADD_TEST(
#    NAME "output_of_test_V2.sh"
#	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
#	${CMAKE_COMMAND} -E compare_files "test_V2.out" "execute_test_V2.sh"
#)
# ADD_TEST(
#   NAME "output_of_test_V2.sh"
#   WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
#   COMMAND
#       "diff" "test_V2.out" "/tmp/test_V2.tmp"
#)
ADD_TEST(
	NAME "execute_test_V3.sh"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND
		"sh" "-ve" "test_V3.sh"
#		[ 2>&1 | perl -ne '$ok=1 if /#!/;print $_ if $ok' ] # see comment below
#        "> /tmp/test_V3.tmp "
)
# ADD_TEST(
#	NAME "output_of_test_V3.sh"
#	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
#	COMMAND
#		"diff" "test_V3.out" "/tmp/test_V3.tmp"
#)

ADD_TEST(
	NAME "execution_of_test_scan_search.pl"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND perl test_scan_search.pl all.prf all.seq /tmp
)

ADD_TEST(
	NAME "execution_of_test_pfsearchV3_iupac.pl"
	WORKING_DIRECTORY "${TESTS_DIRECTORY}/examples/"
	COMMAND "perl" "test_pfsearchV3_iupac.pl" "/tmp"
)
# Comment: In our development environment, the command
# "sh -ve ... 2>&1" was printing an unexpected header,
# that was suppressed by the supplied Perl one-liner.
# This should not impact overwise the test results

