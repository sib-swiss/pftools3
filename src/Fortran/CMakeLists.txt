SET_PROPERTY(SOURCE gtop.f htop.f ptoh.f ptof.f 2ft.f 6ft.f psa2msa.f pfscan.f pfsearch.f
                    pfw.f pfscale.f pfmake.f
             APPEND PROPERTY COMPILE_FLAGS "-fno-automatic" )

######################################## DATA #######################################
SET(DAT_OBJECTS psdat.f gsdat.f djdat.f nodat.f codat.f pfdat.f dfdat.f)

######################################## GTOP #######################################
ADD_EXECUTABLE(gtop gtop.f io.c)

######################################## HTOP #######################################
ADD_EXECUTABLE(htop htop.f)

######################################## PTOH #######################################
ADD_EXECUTABLE(ptoh ptoh.f)

######################################## PTOF #######################################
ADD_EXECUTABLE(ptof ptof.f)

######################################## PFW ########################################
ADD_EXECUTABLE(pfw io.c pfw.f)

######################################## 2FT ########################################
ADD_EXECUTABLE(2ft 2ft.f)

######################################## 6FT ########################################
ADD_EXECUTABLE(6ft 6ft.f)

###################################### PSA2MSA ######################################
ADD_EXECUTABLE(psa2msa io.c psa2msa.f)

###################################### PFSCAN #######################################
ADD_EXECUTABLE(pfscan pfscan.f)

###################################### PFMAKE #######################################
ADD_EXECUTABLE(pfmake pfmake.f)

##################################### PFSCALE #######################################
ADD_EXECUTABLE(pfscale pfscale.f)

################################### OLD PFSEARCH ####################################
ADD_EXECUTABLE(pfsearch pfsearch.f)

############################### QUANTILE REGRESSION #################################
ADD_LIBRARY(quantile OBJECT rq0.f)
set_target_properties( quantile PROPERTIES LINKER_LANGUAGE C)

#########################################################################################
# Installation
INSTALL(TARGETS gtop htop ptoh ptof pfw 2ft 6ft psa2msa pfscale pfsearch pfscan pfmake
        DESTINATION bin
        COMPONENT Applications
)
