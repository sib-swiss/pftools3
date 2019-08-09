#########################################################################################
## Data files

SET(CalibrationFiles
    Calibration/window20.seq
    Calibration/reversed.seq
)


IF(BUILD_TESTING)
# 	IF (NOT EXISTS ${PROJECT_BINARY_DIR}/prosite.dat)
# 		MESSAGE(STATUS "Downloading prosite.dat")
# 		FILE(DOWNLOAD ftp://ftp.expasy.org/databases/prosite/prosite.dat ${PROJECT_BINARY_DIR}/prosite.dat INACTIVITY_TIMEOUT 60 STATUS PROSITE_OK)
# 		IF(NOT PROSITE_OK)
# 			MESSAGE(WARNING "Unable to download prosite.dat from the Internet.")
# 		ENDIF(NOT PROSITE_OK)
# 	ENDIF (NOT EXISTS ${PROJECT_BINARY_DIR}/prosite.dat)
#
# 	IF (NOT EXISTS ${PROJECT_BINARY_DIR}/uniprot_sprot.fasta.gz)
# 		MESSAGE(STATUS "Downloading uniprot_sprot.fasta.gz")
# 		FILE(DOWNLOAD ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz ${PROJECT_BINARY_DIR}/uniprot_sprot.fasta.gz INACTIVITY_TIMEOUT 60 STATUS SWISSPROT_OK)
# 		IF(NOT SWISSPROT_OK)
# 			MESSAGE(WARNING "Unable to download uniprot_sprot.fasta.gz from the Internet.")
# 		ENDIF(NOT SWISSPROT_OK)
# 		execute_process(
# 			COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_BINARY_DIR}/uniprot_sprot.fasta.gz
# 			WORKING_DIRECTORY${PROJECT_BINARY_DIR}
# 		)
# 	ENDIF (NOT EXISTS ${PROJECT_BINARY_DIR}/uniprot_sprot.fasta.gz)
ENDIF(BUILD_TESTING)

# Find all the input files
# FILE(GLOB PrfFiles RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/prf/*.prf")
#
# FILE(GLOB PrfSeqFiles RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/prfseq/*.msa")

FILE(GLOB blossum   Matrices/*.cmp)
FILE(GLOB reverse   reverse/*)
FILE(GLOB familyPrf FamilyProfile/*)

INSTALL(FILES  ${CalibrationFiles}
	DESTINATION ${DATA_PATH}/calibration
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT DataCalibration
)

INSTALL(FILES  ${blossum}
	DESTINATION ${DATA_PATH}/cmp
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT DataBlossum
)

INSTALL(FILES  ${reverse}
	DESTINATION ${DATA_PATH}/reverse
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT DataReverse
)

INSTALL(FILES  ${familyPrf}
	DESTINATION ${DATA_PATH}/familyProfile
	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
	COMPONENT DataFamPrf
)

#INSTALL(FILES  ${PrfFiles}
# 	DESTINATION ${DATA_PATH}/prf
# 	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
# 	COMPONENT DataPrf
#)
#
#INSTALL(FILES  ${PrfSeqFiles}
# 	DESTINATION ${DATA_PATH}/prfseq
# 	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
# 	COMPONENT DataPrfSeq
#)
#
# INSTALL(FILES uniprot.fas
# 	DESTINATION ${DATA_PATH}/uniprot/
# 	PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
# 	COMPONENT DataUniprot
#)

