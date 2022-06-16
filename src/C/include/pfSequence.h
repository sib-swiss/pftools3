/*******************************************************
                        PFTOOLS
 *******************************************************
  May 29, 2013 pfSequence.h
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#ifndef _SEQUENCE_H
#define _SEQUENCE_H
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include "pfConfig.h"

#ifndef PFSEQ
#define PFSEQ
/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
typedef struct PFSequence {
  unsigned char * ProfileIndex;
  unsigned char * OriginalSequence;
  size_t Length;
} PFSequence;
#endif

typedef struct ss_Data {
	off_t Offset;
  unsigned int HeaderLength;
  unsigned int SequenceLength;
} ss_Data;

typedef struct s_Data {
	ss_Data Sequence;
	ss_Data Quality;
} s_Data;

typedef struct DBSequence {
  char FileName[256];
  s_Data  *DataPtr;
  off_t FileSize;
	struct timespec LastModification;
  size_t SequenceCount;
  size_t MaxSequenceSize; // this includes carriage returns and header
} DBSequence_t;

typedef struct Sequence {
  union {
    char * Header;
    void * Memory;
  } Data;
  size_t Size;
  PFSequence ProfileData;
} Sequence;

/*
 ************************************************************************************************
 *                             SEQUENCE FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
int AnalyzeFASTAStructure(const char * const FileName, DBSequence_t * const Info);
int AnalyzeFASTQStructure(const char * const FileName, DBSequence_t * const Info);
int AnalyzeEMBLStructure(const char * const FileName, DBSequence_t * const Info);
int ExportDBStructure(FILE* const stream, const DBSequence_t * const Info);
int ImportDBStructure(FILE* const stream, const char * const DBFileName, DBSequence_t * const Info);

#ifndef __USE_INLINE_FUNCTIONS__
char * CleanSequence(PFSequence * const Sequence);
PFSequence * ReadSequenceIndex(Sequence * const Seq, FILE * const stream, const ss_Data * const DataPtr);
PFSequence * MMAP_ReadSequenceIndex(Sequence * const Seq, const char * const restrict Array,
				    const ss_Data * const DataPtr, const off_t InitialArrayOffset
#ifdef MMAP_DEBUG
				    ,const size_t ThreadId, const size_t NodeId, const size_t length
#endif
);
void ReadSequenceNameIndex(char * const Name, FILE * const stream, const ss_Data * const DataPtr);
unsigned char TranslateCharToIndex(const char letter, const unsigned char * restrict const Alphabet);
PFSequence * TranslateSequenceToIndex(PFSequence * const Sequence, const unsigned char * restrict const Alphabet, const int complement );
void ReverseTranslatedSequence(PFSequence * const Sequence);
#else
#include "pfSequenceInline.h"
#endif

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */
static inline void FreeDBStructure(DBSequence_t * Info)
{
    free(Info->DataPtr);
}


#endif /* _SEQUENCE_H */
