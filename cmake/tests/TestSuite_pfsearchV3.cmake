MESSAGE(STATUS "Building test suite for pfsearch in ${TESTS_DIRECTORY}")
INCLUDE(${PROJECT_SOURCE_DIR}/cmake/tests/macros.cmake)

##############################################################################################
## SH2 tests
IF (NOT EXISTS ${TESTS_DIRECTORY}/SH2/)
  FILE(MAKE_DIRECTORY ${TESTS_DIRECTORY}/SH2/)
ENDIF(NOT EXISTS ${TESTS_DIRECTORY}/SH2/)

# MESSAGE(STATUS "PFSEARCH STAGE 0: output format verification")

FOREACH(frmt RANGE 0 6)
	IF ( NOT ${frmt} STREQUAL "7" )
		TEST_PFSEARCHV3(SH2 "format_${frmt}" sh2.prf sequence.seq "-n" "-o\ ${frmt}" output_format${frmt})
	ENDIF( NOT ${frmt} STREQUAL "7" )
ENDFOREACH(frmt RANGE 0 8)

# MESSAGE(STATUS "PFSEARCH STAGE 1: reverse sequence verification")
TEST_PFSEARCHV3(SH2 Reverse_seq sh2.prf sequence_reverse.seq "-n --reverse" "-o\ 1" output_format_rev)

COMPARE_PFSEARCH(SH2 compare sh2.prf sequence.seq "-n" "")

##############################################################################################
## 
