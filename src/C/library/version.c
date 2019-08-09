/*******************************************************
                        PFTOOLS
 *******************************************************
  Aug 23, 2013 version.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/

#include "config.h"

char Version[] = PF_VERSION;

const char * GetVersion() 
{
  return Version;
}