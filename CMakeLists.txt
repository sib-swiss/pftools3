#########################################################################################
# PfTools versioning
#
SET(VERSION "3.2.12")
SET(MAJOR_VERSION "3")
SET(MINOR_VERSION "2")
SET(PATCH_VERSION "12")
SET(DATA_PATH "share/data"     CACHE STRING "Installation data folder" )
SET(DOC_PATH  "share/doc"      CACHE STRING "Installation doc folder" )
SET(MAN1_PATH "share/man/man1" CACHE STRING "Installation man1 folder" )
SET(MAN5_PATH "share/man/man5" CACHE STRING "Installation man5 folder" )

#########################################################################################
# Header Message
#
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "|                          PfTools   v${VERSION}                          |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")
MESSAGE(STATUS "|     (C) Copyright SIB Swiss Institute of Bioinformatics            |")
MESSAGE(STATUS "|         Thierry Schuepbach (projects@vital-it.ch)                  |")
MESSAGE(STATUS "|                                                                    |")
MESSAGE(STATUS "|     PfTools is available from                                      |")
MESSAGE(STATUS "|         https://github.com/sib-swiss/pftools3                      |")
MESSAGE(STATUS "|     under the GPL v2. See LICENSE.                                 |")
MESSAGE(STATUS "|                                                                    |")
MESSAGE(STATUS "+--------------------------------------------------------------------+")

#########################################################################################
# CMAKE
#
CMAKE_MINIMUM_REQUIRED(VERSION 3.7)
PROJECT(PFTOOLS Fortran C)
IF(COMMAND cmake_policy)
	CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)
SET(CMAKE_MODULE_PATH ${PFTOOLS_SOURCE_DIR}/cmake)

#########################################################################################
# Architectures
#
INCLUDE(CMakeDetermineSystem)
#INCLUDE(CMakePrintSystemInformation)
IF(CMAKE_SYSTEM_PROCESSOR)
	MESSAGE(STATUS "Compilation on architecture ${CMAKE_SYSTEM_PROCESSOR}.")
ENDIF(CMAKE_SYSTEM_PROCESSOR)

#########################################################################################
# What type of build
#
IF(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#########################################################################################
# Determine C Compiler
#
IF(NOT CMAKE_COMPILER_IS_GNUCC)
	EXEC_PROGRAM(${CMAKE_C_COMPILER} ARGS ${CMAKE_BOOT_C_FLAGS} -E "\"${PFTOOLS_SOURCE_DIR}/cmake/CMakeTestIntel.c\"" OUTPUT_VARIABLE CMAKE_COMPILER_OUTPUT RETURN_VALUE CMAKE_COMPILER_RETURN)
	SET(CMAKE_COMPILER_IS_GNUCC_RUN 1)
	IF(NOT CMAKE_COMPILER_RETURN)
		IF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
			INCLUDE(${PFTOOLS_SOURCE_DIR}/cmake/Linux-icc.cmake)
		ENDIF("${CMAKE_COMPILER_OUTPUT}" MATCHES ".*THIS_IS_INTEL.*" )
	ENDIF(NOT CMAKE_COMPILER_RETURN)
ENDIF(NOT CMAKE_COMPILER_IS_GNUCC)

INCLUDE(CheckCCompilerFlag)
CHECK_C_COMPILER_FLAG(-msse HANDLE_SSE)
IF(HANDLE_SSE)
	SET(CMAKE_C_SSE_FLAGS "-msse")
ELSE(HANDLE_SSE)
	MESSAGE(FATAL_ERROR "C compiler does not support SSE")
ENDIF(HANDLE_SSE)
CHECK_C_COMPILER_FLAG(-msse2 HANDLE_SSE2)
IF(HANDLE_SSE2)
	SET(CMAKE_C_SSE2_FLAGS "-msse2")
ELSE(HANDLE_SSE2)
	MESSAGE(FATAL_ERROR "C compiler does not support SSE 2")
ENDIF(HANDLE_SSE2)
SET(SSE41_FLAGS "-msse4.1" "-msse41")
UNSET(CMAKE_C_SSE41_FLAGS)
FOREACH(flag IN LISTS SSE41_FLAGS)
	MESSAGE(STATUS "testing flag ${flag}...")
	CHECK_C_COMPILER_FLAG("${flag}" HANDLE_SSE41)
	IF(HANDLE_SSE41)
		SET(CMAKE_C_SSE41_FLAGS "${flag}")
		BREAK()
	ENDIF(HANDLE_SSE41)
ENDFOREACH(flag)
IF(NOT DEFINED CMAKE_C_SSE41_FLAGS)
	MESSAGE(FATAL_ERROR "C compiler does not support SSE 4.1")
ENDIF(NOT DEFINED CMAKE_C_SSE41_FLAGS)

CHECK_C_COMPILER_FLAG("-std=c99" HANDLE_C99)
IF(NOT HANDLE_C99)
	MESSAGE(FATAL_ERROR "PFTOOLS require a C compiler with C99 support")
ELSE(NOT HANDLE_C99)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
	MESSAGE(STATUS "Add -std=c99 to C compiler options")
ENDIF(NOT HANDLE_C99)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_SSE2_FLAGS}")
MESSAGE(STATUS "Add SSE2 to C compiler options")

IF(CMAKE_COMPILER_IS_GNUCC)
	CHECK_C_COMPILER_FLAG("-fms-extensions" MS_EXTENSION)
	IF(MS_EXTENSION)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fms-extensions")
	ENDIF(MS_EXTENSION)
ENDIF(CMAKE_COMPILER_IS_GNUCC)

#########################################################################################
# Determine Fortran Compiler
#

#########################################################################################
# External Dependencies
#
SET(FIND_LIBRARY_USE_LIB64_PATHS 1)

# Check for pthread
FIND_PACKAGE(Threads REQUIRED)
IF( ${CMAKE_USE_WIN32_THREADS_INIT})
	MESSAGE(STATUS "Using Windows API for threads")
	SET(USE_WINAPI 1)
ENDIF(${CMAKE_USE_WIN32_THREADS_INIT})

#Check for presence of Perl for profile generation
FIND_PACKAGE(Perl)
IF(PERL_FOUND)
	MESSAGE(STATUS "Perl found - full test suite usable")
ELSE(PERL_FOUND)
	MESSAGE(STATUS "Perl not found - test suite only partially usable")
ENDIF(PERL_FOUND)
### insert here condition for which we do full test in case Perl is found
### or a slim version if Perl is absent

#########################################################################################
# Potentially Internal Dependencies
#
OPTION(STANDALONE "This will compile executable in static mode." OFF)
OPTION(USE_PDF "Enable output in pdf file format." OFF)
OPTION(USE_GRAPHICS "Enable graphics generation." OFF)
OPTION(USE_PCRE "Enable Perl Regex." OFF)
OPTION(USE_PCRE2 "Enable Perl Regex." ON)

# Check for PCRE2 library
IF(USE_PCRE2)
	SET(USE_PCRE OFF)
	FIND_PACKAGE(PCRE2)
	IF (NOT PCRE2_FOUND)
		MESSAGE(STATUS "libpcre2 not found")
	ENDIF(NOT PCRE2_FOUND)
ENDIF(USE_PCRE2)

# Check for PCRE library
IF(USE_PCRE)
	IF(STANDALONE)
		SET(BUILD_PCRE ON)
	ELSE(STANDALONE)
		OPTION(BUILD_PCRE "Build libpcre from internal source file" OFF)
		IF (NOT BUILD_PCRE)
			FIND_PACKAGE(PCRE)
			IF (NOT PCRE_FOUND)
				MESSAGE(STATUS "libpcre not found, enforcing internal version to be built...")
				SET(BUILD_PCRE ON CACHE BOOL "" FORCE)
			ENDIF(NOT PCRE_FOUND)
		ENDIF(NOT BUILD_PCRE)
	ENDIF(STANDALONE)
ENDIF(USE_PCRE)

#Check for HARU PDF library
IF(USE_PDF)
	IF(STANDALONE)
		SET(BUILD_hpdf ON)
		SET(hpdf_FOUND ON)
	ELSE(STANDALONE)
		OPTION(BUILD_hpdf "Build libHaru from internal source file" OFF)
		IF(NOT BUILD_hpdf)
			FIND_PACKAGE(hpdf)
			IF(NOT hpdf_FOUND)
				MESSAGE(STATUS "libHaru not found, enforcing internal version to be built...")
				SET(BUILD_hpdf ON CACHE BOOL "" FORCE)
				SET(hpdf_FOUND ON)
			ENDIF(NOT hpdf_FOUND)
		ELSE (NOT BUILD_hpdf)
			MESSAGE(STATUS "libHaru taken from internal source files")
		ENDIF(NOT BUILD_hpdf)
	ENDIF(STANDALONE)
ENDIF(USE_PDF)

IF(USE_GRAPHICS)
	IF(STANDALONE)
		SET(BUILD_GD ON)
		SET(BUILD_PLplot ON)
		SET(PLplot_FOUND ON)
		SET(USE_PLPLOT ON)
	ELSE(STANDALONE)
		OPTION(BUILD_GD "Build libGD from internal source file" OFF)
		#Check for GD library
		IF( NOT BUILD_GD)
			FIND_PACKAGE(GD)
			IF(GD_FOUND)
				IF(GD_SUPPORTS_PNG AND GD_SUPPORTS_JPEG)
					MESSAGE(STATUS "System GD supports PNG and JPEG")
					SET(USE_GD ON)
				ELSE(GD_SUPPORTS_PNG AND GD_SUPPORTS_JPEG)
					IF(GD_SUPPORTS_JPEG)
						MESSAGE(STATUS "System libGD does not support PNG file format")
					ELSE(GD_SUPPORTS_JPEG)
						MESSAGE(STATUS "System libGD does not supprt JPEG file format")
					ENDIF(GD_SUPPORTS_JPEG)
					MESSAGE(STATUS "libGD does not provides all features, enforcing internal version to be built...")
					SET(BUILD_GD ON)
				ENDIF(GD_SUPPORTS_PNG AND GD_SUPPORTS_JPEG)
			ELSE(GD_FOUND)
				MESSAGE(STATUS "libGD not found,  enforcing internal version to be built...")
				SET(BUILD_GD ON CACHE BOOL "" FORCE)
			ENDIF(GD_FOUND)
		ELSE( NOT BUILD_GD)
			MESSAGE(STATUS "libGD taken from internal source file")
		ENDIF( NOT BUILD_GD)

		OPTION(BUILD_PLplot "Build PLplot from internal source file" OFF)
		#Check for PLplot library and its dependency with GD
		IF(NOT BUILD_PLplot)
			FIND_PACKAGE(PLplot)
			IF(PLplot_FOUND)
				check_library_exists(
			"${PLplot_LIBRARIES}"
			"plD_dispatch_init_jpeg"
			"${PLplot_LIBRARY_DIR}"
			PLplot_SUPPORTS_JPEG
			)
				IF(PLplot_SUPPORTS_JPEG)
					SET(USE_PLPLOT ON)
					MESSAGE(STATUS "PLplot supports JPEG file format")
				ELSE(PLplot_SUPPORTS_JPEG)
					MESSAGE(STATUS "PLplot does not support JPEG file format")
					UNSET(PLplot_FOUND)
				ENDIF(PLplot_SUPPORTS_JPEG)
			ENDIF(PLplot_FOUND)
			IF(NOT PLplot_FOUND)
				MESSAGE(STATUS "PLplot not functional, enforcing internal version to be built...")
				SET(BUILD_PLplot ON CACHE BOOL "" FORCE)
				SET(PLplot_FOUND ON)
				SET(USE_PLPLOT ON)
			ENDIF(NOT PLplot_FOUND)
		ELSE(NOT BUILD_PLplot)
			MESSAGE(STATUS "PLplot not found, enforcing internal version to be built...")
		ENDIF(NOT BUILD_PLplot)
	ENDIF(STANDALONE)
ENDIF(USE_GRAPHICS)

#########################################################################################
# EXTRA HEADER FILE
#
INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (emmintrin.h HAVE_SSE2_H)
CHECK_INCLUDE_FILES (smmintrin.h HAVE_SSE_41_H)
IF (NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)
#	MESSAGE(FATAL_ERROR "PFTOOLS require a compiler with both SSE 2 and 4.1 support")
ENDIF(NOT HAVE_SSE2_H OR NOT HAVE_SSE_41_H)

CHECK_INCLUDE_FILES (mm_malloc.h HAVE_MM_MALLOC_H)
CHECK_INCLUDE_FILES (alloca.h HAVE_ALLOCA_H)

SET(MMAP ON)
IF(MINGW)
	CHECK_INCLUDE_FILES (sys/mman.h HAVE_MMAN_H)
	IF(NOT HAVE_MMAN_H)
		MESSAGE(WARNING "MinGW package mman is missing, memory mapping cannot be used")
		SET(MMAP OFF)
	ENDIF(NOT HAVE_MMAN_H)
ENDIF(MINGW)

#########################################################################################
# Compiling options
#
# Use RPATH?
IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	# No rpath on Darwin. Setting it will only cause trouble.
ELSE(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	OPTION(USE_RPATH "Use -rpath when linking libraries, executables" ON)
	IF(USE_RPATH)
		SET(CMAKE_EXE_LINKER_FLAGS "-Wl,--enable-new-dtags")
		# use, i.e. don't skip the full RPATH for the build tree
		SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

		# when building, don't use the install RPATH already
		# (but later on when installing)
		SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

		# add the automatically determined parts of the RPATH
		# which point to directories outside the build tree to the install RPATH
		SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

		# the RPATH to be used when installing, but only if it's not a system directory
		LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
		IF("${isSystemDir}" STREQUAL "-1")
			SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
		ENDIF("${isSystemDir}" STREQUAL "-1")
	ENDIF(USE_RPATH)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

OPTION(USE_AFFINITY "CPU core affinity will be used for threads." ON)
IF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
	SET(USE_AFFINITY OFF)
ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Darwin")

IF(MMAP)
	OPTION(USE_MMAP "DB files wil be mapped to memory for speed." ON)
ENDIF(MMAP)

OPTION(BUILD_SHARED_LIBS "Builds shared development library." OFF)
OPTION(BUILD_STATIC_LIBS "Builds static development library." OFF)

OPTION(USE_32BIT_INTEGER "Builds using 32bit integer format rather than 16 bits." OFF)

#########################################################################################
# CONFIGURATION FILE
#
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/C/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/C/include/pfConfig.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/pfConfig.h)
ADD_DEFINITIONS("-I ${CMAKE_CURRENT_BINARY_DIR}/")

#########################################################################################
# STATIC BUILD
#
IF(STANDALONE)
	FIND_LIBRARY(ZLIB_STATIC libz.a)
	IF(NOT ZLIB_STATIC)
		MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library zlib\n")
	ENDIF(NOT ZLIB_STATIC)

	IF(USE_PCRE2)
		#NOTE For PCRE2 with PCRE2_CODE_UNIT_WIDTH = 8
		FIND_LIBRARY(PCRE2_STATIC libpcre2-8.a)
		IF(NOT PCRE2_STATIC)
			MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library PCRE2\n")
		ENDIF(NOT PCRE2_STATIC)
	ENDIF(USE_PCRE2)
	IF(USE_PCRE)
		FIND_LIBRARY(PCRE_STATIC libpcre.a)
		IF(NOT PCRE_STATIC)
			MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library PCRE\n")
		ENDIF(NOT PCRE_STATIC)
	ENDIF(USE_PCRE)

	IF(USE_GRAPHICS)
		FIND_LIBRARY(PNG_STATIC libpng.a)
		IF(NOT PNG_STATIC)
			MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library PNG\n")
		ENDIF(NOT PNG_STATIC)

		FIND_LIBRARY(JPEG_STATIC libjpeg.a)
		IF(NOT JPEG_STATIC)
			MESSAGE(FATAL_ERROR "Standalone static building is not possible without static library JPEG\n")
		ENDIF(NOT JPEG_STATIC)
	ENDIF(USE_GRAPHICS)

	ADD_DEFINITIONS("-DNO_USERNAME")
	SET(USE_RPATH OFF)
	SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
#	SET (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
	SET (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "-static")
	SET (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "-static")

	SET(CMAKE_SKIP_INSTALL_RPATH ON)
	SET(CMAKE_SKIP_RPATH ON)

ENDIF(STANDALONE)

#########################################################################################
# Compiling external tools
#
IF(USE_PDF)
	IF(BUILD_hpdf)
		INCLUDE(cmake/build_internal_hpdf.cmake)
	ENDIF(BUILD_hpdf)
ENDIF(USE_PDF)
IF(USE_GRAPHICS)
	IF(BUILD_GD)
		INCLUDE(cmake/build_internal_GD.cmake)
		SET(USE_GD ON)
	ENDIF(BUILD_GD)
	IF(BUILD_PLplot)
		INCLUDE(cmake/build_internal_PLplot.cmake)
		SET(PLplot_FOUND ON)
		SET(USE_PLPLOT ON)
	ENDIF(BUILD_PLplot)
ENDIF(USE_GRAPHICS)

IF(BUILD_PCRE)
	INCLUDE(cmake/build_internal_pcre.cmake)
ENDIF(BUILD_PCRE)

#########################################################################################
# CONFIGURATION FILE
#
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/C/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/C/include/pfConfig.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/pfConfig.h)
ADD_DEFINITIONS("-I ${CMAKE_CURRENT_BINARY_DIR}/")

#########################################################################################
# Now treats subdirectories
#
ADD_SUBDIRECTORY(${PFTOOLS_SOURCE_DIR}/src)
LINK_DIRECTORIES(${PFTOOLS_SOURCE_DIR}/Build/${CONF}/)

#########################################################################################
## Data files
SET(DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)
ADD_SUBDIRECTORY(data)

#########################################################################################
# Enabling TESTING
include(CTest)
SET(TESTS_DIRECTORY "${PROJECT_SOURCE_DIR}/testing" CACHE PATH "Directory for tests")
SET(TESTS_DATA_DIR "${PROJECT_SOURCE_DIR}/tests")
IF (NOT EXISTS "${TESTS_DIRECTORY}" )
	FILE(MAKE_DIRECTORY "${TESTS_DIRECTORY}" )
ENDIF (NOT EXISTS "${TESTS_DIRECTORY}" )

IF(BUILD_TESTING)
	ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTING)

ENABLE_TESTING()
ADD_SUBDIRECTORY( "${TESTS_DIRECTORY}/examples" )

#########################################################################################
## Documentation
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(man)

#########################################################################################
## Packaging
SET (CPACK_PACKAGE_NAME PFtools)
SET (CPACK_PACKAGE_CONTACT "Thierry Schuepbach")
SET (CPACK_PACKAGE_VENDOR "SIB Swiss Institute of Bioinformatics")
SET (CPACK_PACKAGE_VERSION "${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH-VERSION}")
SET (CPACK_PACKAGE_VERSION_MAJOR   ${MAJOR_VERSION})
SET (CPACK_PACKAGE_VERSION_MINOR  ${MINOR_VERSION})
SET (CPACK_PACKAGE_VERSION_PATCH  ${PATCH-VERSION})
SET (CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}-${CMAKE_SYSTEM_PROCESSOR}")
SET (CPACK_PACKAGE_INSTALL_DIRECTORY "${CPACK_PACKAGE_NAME}/${CPACK_PACKAGE_VERSION}")
SET (CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)

SET(CPACK_COMPONENT_APPLICATIONS_DISPLAY_NAME "PFtools suite")
SET(CPACK_COMPONENT_APPLICATIONS_DESCRIPTION  "Binary executables")

SET(CPACK_COMPONENT_DATAPRF_DISPLAY_NAME "Profiles files")
SET(CPACK_COMPONENT_DATAPRF_DESCRIPTION  "Profiles data files")
SET(CPACK_COMPONENT_DATAPRF_DISABLED TRUE)
SET(CPACK_COMPONENT_DATAPRF_DEPENDS Applications)
SET(CPACK_COMPONENT_DATAPRF_GROUP "data")

SET(CPACK_COMPONENT_DATAPRFSEQ_DISPLAY_NAME "Profiles sequence files")
SET(CPACK_COMPONENT_DATAPRFSEQ_DESCRIPTION  "Profiles sequence data files")
SET(CPACK_COMPONENT_DATAPRFSEQ_DISABLED TRUE)
SET(CPACK_COMPONENT_DATAPRFSEQ_DEPENDS Applications)
SET(CPACK_COMPONENT_DATAPRFSEQ_GROUP "data")

SET(CPACK_COMPONENT_DATACALIBRATION_DISPLAY_NAME "Calibration database files")
SET(CPACK_COMPONENT_DATACALIBRATION_DESCRIPTION  "Calibration database files for pfcalibration")
SET(CPACK_COMPONENT_DATACALIBRATION_DISABLED TRUE)
SET(CPACK_COMPONENT_DATACALIBRATION_DEPENDS Applications)
SET(CPACK_COMPONENT_DATACALIBRATION_GROUP "data")

SET(CPACK_COMPONENT_DATABLOSSUM_DISPLAY_NAME "Substitution matrix files")
SET(CPACK_COMPONENT_DATABLOSSUM_DESCRIPTION  "Substitution matrix files")
SET(CPACK_COMPONENT_DATABLOSSUM_DISABLED TRUE)
SET(CPACK_COMPONENT_DATABLOSSUM_DEPENDS Applications)
SET(CPACK_COMPONENT_DATABLOSSUM_GROUP "data")

SET(CPACK_COMPONENT_DATAUNIPROT_DISPLAY_NAME "Uniprot sequence database")
SET(CPACK_COMPONENT_DATAUNIPROT_DESCRIPTION  "Uniprot sequence database")
SET(CPACK_COMPONENT_DATAUNIPROT_DISABLED TRUE)
SET(CPACK_COMPONENT_DATAUNIPROT_DEPENDS Applications)

SET(CPACK_COMPONENT_DOCUMENTATION_DISPLAY_NAME "Documentation")
SET(CPACK_COMPONENT_DOCUMENTATION_DESCRIPTION  "Documentation")
SET(CPACK_COMPONENT_DOCUMENTATION_GROUP "doc")

SET(CPACK_COMPONENT_DEV_DISPLAY_NAME "Development library")
SET(CPACK_COMPONENT_DEV_DESCRIPTION  "Development library")
SET(CPACK_COMPONENT_DEV_DISABLED TRUE)
SET(CPACK_COMPONENT_DEV_GROUP "dev")

SET (CPACK_COMPONENTS_ALL Applications DataPrf DataPrfSeq DataCalibration
		DataBlossum DataReverse DataFamPrf Documentation)


SET (CPACK_RESOURCE_FILE_LICENSE "${PFTOOLS_SOURCE_DIR}/LICENSE")

UNSET(CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE)
IF(UNIX AND NOT WIN32)
	SET (CPACK_GENERATOR "TBZ2")
	SET (CPACK_DEB_COMPONENT_INSTALL ON)
	INSTALL(FILES "${PFTOOLS_SOURCE_DIR}/LICENSE"
	              "${PFTOOLS_SOURCE_DIR}/README.md"
	              "${PFTOOLS_SOURCE_DIR}/README.pftools.v2.3.md"
	              "${PFTOOLS_SOURCE_DIR}/ChangeLog"
	        DESTINATION .
	        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	        COMPONENT Applications)

	IF(USE_PDF)
		SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libhpdf")
	ENDIF(USE_PDF)
	IF(USE_PCRE2)
		SET(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS},pcre2")
	ENDIF(USE_PCRE2)
	IF(USE_PCRE)
		SET(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS},pcre")
	ENDIF(USE_PCRE)
ELSE(UNIX AND NOT WIN32)
	#SET (CPACK_COMPONENTS_ALL_IN_ONE_PACKAGE TRUE)
ENDIF(UNIX AND NOT WIN32)

SET(CPACK_SOURCE_GENERATOR "TBZ2")
#SET(CPACK_SOURCE_OUTPUT_CONFIG_FILE "${PFTOOLS_SOURCE_DIR}/CPackSourceConfig.cmake")
SET(CPACK_SOURCE_IGNORE_FILES
   "/CVS/"
   "/Build/*"
   "/build/*"
   "/\\\\.svn/*"
   "/doc/Article/pfsearch3_v0.8.odt"
   "/doc/poster/"
   "/data/Benchmarks/"
   ".*~"
   "/doc/Presentation/"
   "/src/C/experimental/"
   "/src/C/external/"
   "/src/C/generic/"
   "/data/.*\\\\.index"
   "/data/.*\\\\.fas"
   ".*\\\\.am"
   ".*\\\\.kdev4"
   "/\\\\.git/*"
   "/generic-cmake.config"
   "/output.c.diff"
)
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "PFTools-${VERSION}")
SET(CPACK_SOURCE_STRIP_FILES "")


INCLUDE(CPack)

#########################################################################################
## Post install script
IF(BUILD_PLplot)
	INSTALL (SCRIPT ${PFTOOLS_SOURCE_DIR}/cmake/PostInstall.cmake)
ENDIF(BUILD_PLplot)

#########################################################################################
## Output global build information
include(summary)
summary()
