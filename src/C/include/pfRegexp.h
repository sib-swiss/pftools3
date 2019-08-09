/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Jun 26, 2013 pfregexp.h
 ***************************************************************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#ifndef _PFREGEX_H
#define _PFREGEX_H
#include "pfConfig.h"
#include "pfProfile.h"

#if defined(USE_PCRE)
# include "pcre.h"
#else
# if defined(__GNUC__)
#   define _POSIX_C_SOURCE
# endif
# include <sys/types.h>
# include <regex.h>
#endif

struct RegEx {
   char * * regexString;
#if defined(USE_PCRE)
   pcre * * regexCompiled;
#else
   regex_t * regexCompiled;
#endif
   size_t count;
   size_t maxMatchCount;
   char * memory;
   size_t memorySize;
};

int InitRegExFromString(struct RegEx * const regex, const char * const regexSource, const _Bool IsAPattern, const size_t MaxMatchCount);
int InitRegExFromProfile(struct RegEx * const regex, const struct Profile * * PatternPrfs, const int PatternProfileCount, const size_t MaxMatchCount);
void FreeRegEx(struct RegEx * const regex);

char * PatternToRegex(const char * const Pattern);

#if defined(USE_PCRE)
void PrintRegex(const char * const restrict regexString, const char * const restrict Sequence,
		const int Matches[], char * const restrict Header, const size_t SeqLength);
#else
void PrintRegex(const char * const restrict regexString, const char * const restrict Sequence,
		const regmatch_t Matches[], char * const restrict Header, const size_t SeqLength);
#endif
#endif /* _PFREGEX_H */