/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Jun 26, 2013 pfregexp.c
 ***************************************************************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "../include/pfRegexp.h"

int InitRegExFromString(struct RegEx * const regex, const char * const regexSource, const _Bool IsAPattern, const size_t MaxMatchCount)
{
  size_t localCount = 0, localmemorySize;
  char * localMemory;
  char * * localregexString;

  localregexString = (char**) malloc(sizeof(char*));
  if (localregexString == NULL) return -4;
  if (!IsAPattern) {
    localmemorySize = strlen(regexSource)+1;
    localMemory = (char*) malloc(localmemorySize*sizeof(char));
    if (localMemory == NULL) {
      free(localregexString);
      return -4;
    }
    memcpy(localMemory, regexSource, localmemorySize*sizeof(char));
    localMemory[localmemorySize] = '\0';
  }
  else {
    localMemory = PatternToRegex(regexSource);
    if (localMemory == NULL) {
      free(localregexString);
      return -4;
    }
  }
  localregexString[0] = localMemory;
  localCount = 1;


  /* Allocate memory for the compiled regex */
#if defined(USE_PCRE2)
  pcre2_code * * const localregexCompiled = (pcre2_code**) malloc(localCount*sizeof(pcre2_code*));
#elif defined(USE_PCRE)
  pcre * * const localregexCompiled = (pcre**) malloc(localCount*sizeof(pcre*));
#else
  regex_t * const localregexCompiled = (regex_t*) malloc(localCount*sizeof(regex_t));
#endif
  if ( localregexCompiled == NULL) {
      free(localregexString);
      free(localMemory);
      return -4;
  }
  /* Assign pointers to structure*/
  regex->maxMatchCount = MaxMatchCount;
  regex->count         = localCount;
  regex->memory        = localMemory;
  regex->memorySize    = localmemorySize;
  regex->regexCompiled = localregexCompiled;
  regex->regexString   = localregexString;

  /* Compile the regex sources */
  char errorBuffer[1024];
  for (size_t i=0; i<localCount; ++i) {
#if defined(USE_PCRE2)
    const char *error;
    int erroffset;
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    localregexCompiled[i] = pcre2_compile ((PCRE2_SPTR)localregexString[i],  /* the pattern */
                      PCRE2_ZERO_TERMINATED,
                      /*PCRE2_MULTILINE*/ 0,
                      &error,           /* for error message */
                      &erroffset,       /* for error offset */
                      0);               /* use default character tables */
    if (!localregexCompiled[i])
    {
        fprintf(stderr, "pcre2_compile failed (offset: %d), %s\n", erroffset, error);
        return -6;
    }
#elif defined(USE_PCRE)
    const char *error;
    int erroffset;
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    localregexCompiled[i] = pcre_compile (localregexString[i],  /* the pattern */
					  /*PCRE_MULTILINE*/ 0,
					  &error,         	/* for error message */
					  &erroffset,     	/* for error offset */
					  0);             	/* use default character tables */
    if (!localregexCompiled[i])
    {
        fprintf(stderr, "pcre_compile failed (offset: %d), %s\n", erroffset, error);
        return -6;
    }
#else
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    const int error = regcomp(&localregexCompiled[i], localregexString[i], REG_EXTENDED /*| REG_NEWLINE*/);
    if (error != 0) {
      regerror(error, &localregexCompiled[i], errorBuffer, 512);
      fprintf(stderr,"REGEX error: %s\n", errorBuffer);
      regex->count = i;
      FreeRegEx(regex);
      return -6;
    }
#endif
  }
  return 0;
}

int InitRegExFromProfile(struct RegEx * const regex, const struct Profile * * PatternPrfs, const int PatternProfileCount, const size_t MaxMatchCount)
{
  char * * localregexString;

  localregexString = (char**) malloc(PatternProfileCount*sizeof(char*));
  if (localregexString == NULL) return -4;
  for (int iPatternPrf=0; iPatternPrf<PatternProfileCount; ++iPatternPrf) {
    localregexString[iPatternPrf] = PatternToRegex(PatternPrfs[iPatternPrf]->Pattern);
    if (localregexString[iPatternPrf] == NULL) {
      while (--iPatternPrf > 0) free(localregexString[iPatternPrf]);
      free(localregexString);
      return -4;
    }
  }

  /* Allocate memory for the compiled regex */
#if defined(USE_PCRE2)
  pcre2_code * * const localregexCompiled = (pcre2_code**) malloc(PatternProfileCount*sizeof(pcre2_code*));
#elif defined(USE_PCRE)
  pcre * * const localregexCompiled = (pcre**) malloc(PatternProfileCount*sizeof(pcre*));
#else
  regex_t * const localregexCompiled = (regex_t*) malloc(PatternProfileCount*sizeof(regex_t));
#endif
  if ( localregexCompiled == NULL) {
      for (int iPatternPrf=0; iPatternPrf<PatternProfileCount; ++iPatternPrf) {
	free(localregexString[iPatternPrf]);
      }
      free(localregexString);
      return -4;
  }
  /* Assign pointers to structure*/
  regex->maxMatchCount = MaxMatchCount;
  regex->count         = PatternProfileCount;
  regex->memory        = NULL;
  regex->memorySize    = 0;
  regex->regexCompiled = localregexCompiled;
  regex->regexString   = localregexString;

  /* Compile the regex sources */
  char errorBuffer[1024];
  for (int i=0; i<PatternProfileCount; ++i) {
#if defined(USE_PCRE2)
    const char *error;
    int erroffset;
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    localregexCompiled[i] = pcre2_compile ((PCRE2_SPTR)localregexString[i],  /* the pattern */
                      PCRE2_ZERO_TERMINATED,
                      /*PCRE2_MULTILINE*/ 0,
                      &error,           /* for error message */
                      &erroffset,       /* for error offset */
                      0);               /* use default character tables */
    if (!localregexCompiled[i])
    {
        fprintf(stderr, "pcre2_compile failed (offset: %d), %s\nProfile ID: %s\npattern: %s\n", erroffset, error,
        PatternPrfs[i]->Identification, PatternPrfs[i]->Pattern);

        return -6;
    }
#elif defined(USE_PCRE)
    const char *error;
    int erroffset;
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    localregexCompiled[i] = pcre_compile (localregexString[i],  /* the pattern */
					  /*PCRE_MULTILINE*/ 0,
					  &error,         	/* for error message */
					  &erroffset,     	/* for error offset */
					  0);             	/* use default character tables */
    if (!localregexCompiled[i])
    {
        fprintf(stderr, "pcre_compile failed (offset: %d), %s\nProfile ID: %s\npattern: %s\n", erroffset, error,
		PatternPrfs[i]->Identification, PatternPrfs[i]->Pattern);

        return -6;
    }
#else
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
    const int error = regcomp(&localregexCompiled[i], localregexString[i], REG_EXTENDED /*| REG_NEWLINE*/);
    if (error != 0) {
      regerror(error, &localregexCompiled[i], errorBuffer, 512);
      fprintf(stderr,"REGEX error: %s\n", errorBuffer);
      regex->count = i;
      FreeRegEx(regex);
      return -6;
    }
#endif
  }
  return 0;
}

void FreeRegEx(struct RegEx * const regex)
{
    if (regex->memory) free(regex->memory);
    if (regex->regexCompiled) {
#if defined(USE_PCRE2)
      for (size_t i=0; i<regex->count; ++i) pcre2_code_free(regex->regexCompiled[i]);
#elif defined(USE_PCRE)
      for (size_t i=0; i<regex->count; ++i) pcre_free(regex->regexCompiled[i]);
#else
      for (size_t i=0; i<regex->count; ++i) regfree(&regex->regexCompiled[i]);
#endif
      free(regex->regexCompiled);
    }
    if (regex->regexString) free(regex->regexString);
    regex->regexString = NULL;
    regex->regexCompiled = NULL;
    regex->memory = NULL;
    regex->memorySize = 0;
    regex->count = 0;
    regex->maxMatchCount = 0;
}

char * PatternToRegex(const char * const Pattern)
{
  const size_t length = strlen(Pattern);
  size_t RegexLength = 0;
  _Bool OpenMisc = false; // To treat misc pattern [..>] or [<...]
  for (size_t i=0; i<length; ++i) {
    switch (Pattern[i]) {
      case '(':
      case ')':
      case '}':
      case 'x':
	++RegexLength;
	break;
      case '{':
	RegexLength += 2;
	break;
      case 'B':
	RegexLength += 3;
	break;
      case 'z':
	RegexLength += 3;
	break;
      case '[':
	if (Pattern[i+1] == '<') {
	 OpenMisc = true;
	 RegexLength += 6;
	}
	else {
	  ++RegexLength;
	}
	break;
      case '>':
	if (Pattern[i+1] == ']') {
	  RegexLength += 8;
	}
	else {
	 ++RegexLength;
	}
	break;
      case ']':
	if (OpenMisc) {
	  OpenMisc = false;
	  RegexLength += 2;
	}
	else {
	 ++RegexLength;
	}
	break;
      case '-':
	break;
      default:
	++RegexLength;
    }
  }
  ++RegexLength;

  char * const regex = (char*) malloc(RegexLength*sizeof(char));
  char * Buffer = (char*) alloca(RegexLength*sizeof(char));
  char * restrict cptr = regex;
  if (regex != NULL) {
    OpenMisc = false;
    size_t i = 0;
    for (size_t i=0; i<length; ++i) {
      switch (Pattern[i]) {
	case '(': *cptr++ = '{'; break;
	case ')': *cptr++ = '}'; break;
	case '}': *cptr++ = ']'; break;
	case 'x': *cptr++ = '.'; break;
	case '{': cptr[0] = '['; cptr[1] = '^'; cptr += 2; break;
	case 'B': cptr[0] = 'N'; cptr[1] = 'D'; cptr[2] = 'B'; cptr += 3; break;
	case 'z': cptr[0] = 'Q'; cptr[1] = 'E'; cptr[2] = 'Z'; cptr += 3; break;
	case '[':
	  if (Pattern[i+1] == '<') {
	    OpenMisc = true;
	    cptr[0] = '('; cptr[1] = '?'; cptr[2] = ':'; cptr[3] = '^'; cptr[4] = '|'; cptr[5] = '['; cptr += 6;
	    i += 1;
	  }
	  else {
	    *cptr++ = '[';
	  }
	  break;
	case '>':
	  if (Pattern[i+1] == ']') {
	    /* Go back to opening [, copying string to buffer */
	    char * restrict Dest = &Buffer[RegexLength-1];
	    *Dest = '\0';
	    while (--cptr >= regex) {
	      if (*cptr != '[') {
		--Dest;
		*Dest = *cptr;
	      }
	      else
		break;
	    }
	    cptr[0] = '('; cptr[1] = '?'; cptr[2] = ':'; cptr[3] = '['; cptr += 4;
	    while (*Dest != '\0') {
	      *cptr++ = *Dest++;
	    }
	    cptr[0] = ']'; cptr[1] = '|'; cptr[2] = '$'; cptr[3] = ')'; cptr += 4;
	    i += 1;
	  }
	  else {
	    *cptr++ = '$';
	  }
	  break;
	case ']':
	  if (OpenMisc) {
	    OpenMisc = false;
	    cptr[0] = ']'; cptr[1] = ')';
	    cptr += 2;
	  }
	  else {
	    *cptr++ = ']';
	  }
	  break;
	case '-':
	  break;
	default:
	  *cptr++ = Pattern[i];
      }
    }

  }
  *cptr = '\0';
  return regex;
}

#if defined(USE_PCRE2)
void PrintRegex(const char * const restrict regexString, const char * const restrict Sequence,
        const int Matches[], char * const restrict Header, const size_t SeqLength)
{
  char * cptr = Header;
  while ( *cptr != ' ' && *cptr != '\0') ++cptr;
  *cptr = '\0';

  size_t count = 0;
  while (Matches[2*count] != -1) {
    fprintf(stdout, "%s_%lu %i-%i %s\n%.*s\n",
        Header, 1+count, 1+Matches[2*count], Matches[2*count+1], regexString,  Matches[2*count+1]-Matches[2*count], Sequence + Matches[2*count] );
    ++count;
  }
}
#elif defined(USE_PCRE)
void PrintRegex(const char * const restrict regexString, const char * const restrict Sequence,
		const int Matches[], char * const restrict Header, const size_t SeqLength)
{
  char * cptr = Header;
  while ( *cptr != ' ' && *cptr != '\0') ++cptr;
  *cptr = '\0';

  size_t count = 0;
  while (Matches[2*count] != -1) {
    fprintf(stdout, "%s_%lu %i-%i %s\n%.*s\n",
	    Header, 1+count, 1+Matches[2*count], Matches[2*count+1], regexString,  Matches[2*count+1]-Matches[2*count], Sequence + Matches[2*count] );
    ++count;
  }
}
#else
void PrintRegex(const char * const restrict regexString, const char * const restrict Sequence,
		const regmatch_t Matches[], char * const restrict Header, const size_t SeqLength)
{
  char * cptr = Header;
  while ( *cptr != ' ' && *cptr != '\0') ++cptr;
  *cptr = '\0';

  size_t count = 0;
  while (Matches[count].rm_so != -1) {
    const size_t s = (size_t) Matches[count].rm_so;
    const size_t e = (size_t) Matches[count].rm_eo;
    fprintf(stdout, "%s_%lu %lu-%lu\n", Header, 1+count, s, e-1);
    const char * restrict start = &Sequence[s];
    const char * const restrict end = &Sequence[e];
    while(start < end) fputc(*start++, stdout);
    fputc('\n', stdout);
    ++count;
  }
}

// int InitRegExFromString(struct RegEx * const regex, const char * const regexSource, const _Bool IsAPattern, const size_t MaxMatchCount)
// {
//   size_t localCount = 0, localmemorySize;
//   char * localMemory;
//   char * * localregexString;
//   struct stat st;
//
//   /* Create a copy of the source on the stack */
//   const size_t regexSourceSize = 1+strlen(regexSource);
//   char * ctmp = (char *) alloca(regexSourceSize*sizeof(char));
//   for (size_t i=0; i<regexSourceSize; ++i) ctmp[i] = regexSource[i];
//   ctmp[regexSourceSize] = '\0';
//
//   /* Check format extract '{' '}' border */
//   char * Start;
//   char * ptr = ctmp;
//   char * const End = &ctmp[regexSourceSize-1];
//
//   while ( *ptr != '{' && ptr <= End) ++ptr;
//   if (ptr == End)
//     return -1;
//   else
//     Start = ++ptr;
//
//   ptr = End;
//   while ( *ptr != '}' && ptr <= End) --ptr;
//   if (ptr == Start)
//     return -1;
//   else
//     *ptr = '\0';
//
//   const size_t id = IsAPattern ? 7 : 5;
//
//   /* What kind of source? file (regexfile()) or string (regex())*/
//   if (regexSource[id] == 'f') {
//     int fd = open(Start, O_RDONLY);
//     if (fd == -1) return -2;
//     if (fstat(fd, &st) != 0) return -3;
//
//     localmemorySize = st.st_size+1;
//     localMemory = (char*) malloc(localmemorySize*sizeof(char));
//     if (localMemory == NULL) return -4;
//
//     ssize_t size = read(fd, localMemory, localmemorySize);
//     if (size != localmemorySize) {
//       close(fd);
//       free(localMemory);
//       return -5;
//     }
//     close(fd);
//     localMemory[localmemorySize] = '\0';
//
//     /* Now count and replace \n with \0 */
//     for (size_t i=0; i<localmemorySize; ++i) {
//       if (localMemory[i] == '\n') ++localCount;
//     }
//     localregexString = (char**) malloc(localCount*sizeof(char*));
//     if (localregexString == NULL) {
//       free(localMemory);
//       return -4;
//     }
//     {
//       size_t id = 0;
//       localregexString[id++] = localMemory;
//       for (size_t i=0; i<localmemorySize; ++i) {
// 	if (localMemory[i] == '\n') {
// 	  localMemory[i] = '\0';
// 	  if (strlen(localregexString[id-1]) > 0) localregexString[id++] = &localMemory[i+1];
// 	}
//       }
//       localCount = id;
//     }
//   }
//   else {
//     localregexString = (char**) malloc(sizeof(char*));
//     if (localregexString == NULL) return -4;
//     if (!IsAPattern) {
//       localmemorySize = strlen(Start)+1;
//       localMemory = (char*) malloc(localmemorySize*sizeof(char));
//       if (localMemory == NULL) {
// 	free(localregexString);
// 	return -4;
//       }
//       memcpy(localMemory, Start, localmemorySize*sizeof(char));
//       localMemory[localmemorySize] = '\0';
//     }
//     else {
//       localMemory = PatternToRegex(Start);
//       if (localMemory == NULL) {
// 	free(localregexString);
// 	return -4;
//       }
//     }
//     localregexString[0] = localMemory;
//     localCount = 1;
//   }
//
//   /* Allocate memory for the compiled regex */
// #if defined(USE_PCRE)
//   pcre * * const localregexCompiled = (pcre**) malloc(localCount*sizeof(pcre*));
// #else
//   regex_t * const localregexCompiled = (regex_t*) malloc(localCount*sizeof(regex_t));
// #endif
//   if ( localregexCompiled == NULL) {
//       free(localregexString);
//       free(localMemory);
//       return -4;
//   }
//   /* Assign pointers to structure*/
//   regex->maxMatchCount = MaxMatchCount;
//   regex->count         = localCount;
//   regex->memory        = localMemory;
//   regex->memorySize    = localmemorySize;
//   regex->regexCompiled = localregexCompiled;
//   regex->regexString   = localregexString;
//
//   /* Compile the regex sources */
//   char errorBuffer[1024];
//   for (size_t i=0; i<localCount; ++i) {
// #if defined(USE_PCRE)
//     const char *error;
//     int erroffset;
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
//     localregexCompiled[i] = pcre_compile (localregexString[i],  /* the pattern */
// 					  /*PCRE_MULTILINE*/ 0,
// 					  &error,         	/* for error message */
// 					  &erroffset,     	/* for error offset */
// 					  0);             	/* use default character tables */
//     if (!localregexCompiled[i])
//     {
//         fprintf(stderr, "pcre_compile failed (offset: %d), %s\n", erroffset, error);
//         return -6;
//     }
// #else
//     fprintf(stderr, "Regex: \'%s\'\n", localregexString[i]);
//     const int error = regcomp(&localregexCompiled[i], localregexString[i], REG_EXTENDED /*| REG_NEWLINE*/);
//     if (error != 0) {
//       regerror(error, &localregexCompiled[i], errorBuffer, 512);
//       fprintf(stderr,"REGEX error: %s\n", errorBuffer);
//       regex->count = i;
//       FreeRegEx(regex);
//       return -6;
//     }
// #endif
//   }
//   return 0;
// }

#endif
