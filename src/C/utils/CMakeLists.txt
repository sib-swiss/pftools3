ADD_DEFINITIONS("-I ${PFTOOLS_SOURCE_DIR}/src/C/include/")

############################### PROFILES ####################################
SET(PROFILE_SOURCES
    io.c
    profile.c
    Normalization.c
)
ADD_LIBRARY(PROFILES OBJECT ${PROFILE_SOURCES})

################################ SYSTEM ######################################
ADD_LIBRARY(SYSTEM OBJECT system.c)

############################### SEQUENCES ####################################
SET(SEQUENCE_SOURCES
    ReadSequence.c
)
ADD_LIBRARY(SEQUENCES OBJECT ${SEQUENCE_SOURCES})

SET(SEQUENCE_EXTRA_SOURCES
    pam.c
)
ADD_LIBRARY(SEQUENCES_EXTRA OBJECT ${SEQUENCE_EXTRA_SOURCES})

############################# OUTPUT FORMAT ##################################
SET(OUTPUT_FORMAT_SOURCES
    output.c
)
ADD_LIBRARY( OUTPUT_FORMAT OBJECT ${OUTPUT_FORMAT_SOURCES})

############################ CALIBRATION TOOLS ###############################
SET(CALIBRATION_SOURCES
    random.c
    histogram.c
    quantileregression.c
    generatesequences.c
)
ADD_LIBRARY(CALIBRATION OBJECT ${CALIBRATION_SOURCES})
ADD_DEPENDENCIES(CALIBRATION quantile)


SET(LIBRARY_PROFILE_SOURCES "${PROFILE_SOURCES}" PARENT_SCOPE)
SET(LIBRARY_SEQUENCE_SOURCES "${SEQUENCE_SOURCES}" PARENT_SCOPE)
SET(LIBRARY_OUTPUT_FORMAT_SOURCES "${OUTPUT_FORMAT_SOURCES}" PARENT_SCOPE)
SET(LIBRARY_CALIBRATION_SOURCES "${CALIBRATION_SOURCES}" "../../Fortran/rq0.f" PARENT_SCOPE)
