IF (NOT EXISTS "${TESTS_DIRECTORY}/examples")
	FILE(MAKE_DIRECTORY "${TESTS_DIRECTORY}/examples")
ENDIF (NOT EXISTS "${TESTS_DIRECTORY}/examples")
SET(EXAMPLES_TMP_DIR "${TESTS_DIRECTORY}/examples")
SET(EXAMPLES_DATA_DIR "${TESTS_DATA_DIR}/examples")
SET(PERL_SCRIPT_DIR "${PROJECT_SOURCE_DIR}/src/Perl")


FILE(GLOB EXAMPLES_FILES examples/*)

FILE(COPY ${EXAMPLES_FILES} DESTINATION ${TESTS_DIRECTORY}/examples/)

CONFIGURE_FILE(${TESTS_DIRECTORY}/examples/test_V2.sh.cmake
               ${TESTS_DIRECTORY}/examples/test_V2.sh.tmp @ONLY)
FILE(GENERATE OUTPUT ${TESTS_DIRECTORY}/examples/test_V2.sh
              INPUT ${TESTS_DIRECTORY}/examples/test_V2.sh.tmp)
FILE(REMOVE ${TESTS_DIRECTORY}/examples/test_V2.sh.cmake)


CONFIGURE_FILE(${TESTS_DIRECTORY}/examples/test_V3.sh.cmake
               ${TESTS_DIRECTORY}/examples/test_V3.sh.tmp @ONLY)
FILE(GENERATE OUTPUT ${TESTS_DIRECTORY}/examples/test_V3.sh
              INPUT ${TESTS_DIRECTORY}/examples/test_V3.sh.tmp)
FILE(REMOVE ${TESTS_DIRECTORY}/examples/test_V3.sh.cmake)


CONFIGURE_FILE(${TESTS_DIRECTORY}/examples/test_pfsearchV3_iupac.pl.cmake
               ${TESTS_DIRECTORY}/examples/test_pfsearchV3_iupac.pl.tmp @ONLY)
FILE(GENERATE OUTPUT ${TESTS_DIRECTORY}/examples/test_pfsearchV3_iupac.pl
              INPUT ${TESTS_DIRECTORY}/examples/test_pfsearchV3_iupac.pl.tmp)
FILE(REMOVE ${TESTS_DIRECTORY}/examples/test_pfsearchV3_iupac.pl.cmake)

CONFIGURE_FILE(${TESTS_DIRECTORY}/examples/test_scan_search.pl.cmake
               ${TESTS_DIRECTORY}/examples/test_scan_search.pl.tmp @ONLY)
FILE(GENERATE OUTPUT ${TESTS_DIRECTORY}/examples/test_scan_search.pl
              INPUT ${TESTS_DIRECTORY}/examples/test_scan_search.pl.tmp)
FILE(REMOVE ${TESTS_DIRECTORY}/examples/test_scan_search.pl.cmake)


INSTALL(DIRECTORY examples/
        DESTINATION share/examples
        PATTERN "CMakeLists.txt" EXCLUDE
)
INSTALL(CODE "EXECUTE_PROCESS(COMMAND find ${CMAKE_INSTALL_PREFIX}/share/examples/ -type f -name \*.cmake -exec chmod a+rx {} \;)")
INSTALL(CODE "EXECUTE_PROCESS(COMMAND find ${CMAKE_INSTALL_PREFIX}/share/examples/ -type f -name \*.cmake -exec sed -i \"s|@DATA_DIR@/Matrices|${CMAKE_INSTALL_PREFIX}/share/data/cmp|; s|@PERL_SCRIPT_DIR@|${CMAKE_INSTALL_PREFIX}/bin|\" {} \;)")
INSTALL(CODE "EXECUTE_PROCESS(COMMAND find ${CMAKE_INSTALL_PREFIX}/share/examples/ -type f -name \*.cmake -exec perl -i -pe \"s|.<TARGET_FILE:(.+?)>|${CMAKE_INSTALL_PREFIX}/bin/$1|\" {} \;)")
INSTALL(CODE "FILE(RENAME ${CMAKE_INSTALL_PREFIX}/share/examples/test_V2.sh.cmake               ${CMAKE_INSTALL_PREFIX}/share/examples/test_V2.sh)")
INSTALL(CODE "FILE(RENAME ${CMAKE_INSTALL_PREFIX}/share/examples/test_V3.sh.cmake               ${CMAKE_INSTALL_PREFIX}/share/examples/test_V3.sh)")
INSTALL(CODE "FILE(RENAME ${CMAKE_INSTALL_PREFIX}/share/examples/test_scan_search.pl.cmake      ${CMAKE_INSTALL_PREFIX}/share/examples/test_scan_search.pl)")
INSTALL(CODE "FILE(RENAME ${CMAKE_INSTALL_PREFIX}/share/examples/test_pfsearchV3_iupac.pl.cmake ${CMAKE_INSTALL_PREFIX}/share/examples/test_pfsearchV3_iupac.pl)")

