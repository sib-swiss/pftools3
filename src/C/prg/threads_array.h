/************************************************************************************************
                        PFTOOLS
 ************************************************************************************************
  Jun 26, 2013 threads_array.h
 ************************************************************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ************************************************************************************************
  THIS FILE IS MEANT TO BE INCLUDED IN THE MAIN PROGRAM DECLARATION FILE
 ************************************************************************************************/

/*
 ************************************************************************************************
 *                                       STRUCTURES                                             *
 ************************************************************************************************
 */
struct ID {
  int PrfId;
  unsigned int SeqId;
};

struct ThreadData {
  const struct Profile * * prf;
  const DBSequence_t * FASTA;
  const char * restrict SequenceFile;
#ifdef _NEEDS_REGEX_
  const struct RegEx * regex;
#endif
#if defined(__USE_MMAP__) && defined(MMAP_DEBUG)
  size_t * maplength;
#endif
  struct ID * Array;
  TransposeMatrix * TransposeMatch;
  size_t MaxProfileSize;
  size_t profileCount;
  size_t start;
  size_t stop;
  size_t threadId;
	float NormalizedCutoff;
  unsigned int counter;
  enum Version version;
};
/*
 ************************************************************************************************
 *                                         GLOBALS                                              *
 ************************************************************************************************
 */
#ifdef _NEEDS_FILTER_
static int (*xali1_ptr)(const struct Profile * const restrict, const unsigned char * const restrict,
			int * const, const size_t, const size_t, const int, const _Bool);
#endif
#ifdef _NEEDS_ALIGNMENT_
static int (*xalip_ptr)( const struct Profile * const restrict, const unsigned char * const restrict,
           union lScores * const restrict, union Positions * const restrict,
           union Positions * const restrict, const int, const int,
           struct Alignment * const restrict, _Bool * const restrict, const size_t, const size_t, const _Bool,
           const int, const size_t);
static int (*xalit_ptr)(const struct Profile * const restrict, const size_t, const size_t, const size_t, const size_t,
          const unsigned char * const restrict, char * const restrict, union lScores * const restrict,
          struct Alignment * const restrict, const _Bool * const restrict);
static PrintFunctionPtr PrintFunction = &PrintPfscan;
static pthread_mutex_t PrintLock;
#endif
#ifdef _NEEDS_REGEX_
static size_t MaxRegexNumber = 16;
#endif
/*
 ************************************************************************************************
 *                                        FUNCTIONS                                             *
 ************************************************************************************************
 */
#ifdef _NEEDS_HEURISTIC_
static void *thread_heuristic(void * _Data) 
{
  Sequence SeqData;
  const struct Profile * const * prfs         = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict FASTA   = ((struct ThreadData*) _Data)->FASTA;
  const TransposeMatrix * const restrict TransposeMatch = ((struct ThreadData*) _Data)->TransposeMatch;
  struct ID * const restrict Array            = ((struct ThreadData*) _Data)->Array;
  PFSequence * PFSeq;
  
  const HeuristicFunctionPtr heuristic = GetHeuristicVersion(((struct ThreadData*) _Data)->version);

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(2*(1+FASTA->MaxSequenceSize)*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fputs("Thread Cannot allocate memory for sequence.\n", stderr);
    return (void*) 1;
  }
  char * const Buffer = &SeqData.Data.Header[1+FASTA->MaxSequenceSize];

  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start              = ((struct ThreadData*) _Data)->start;
  size_t Stop               = ((struct ThreadData*) _Data)->stop;

  /* LOOPS ON SEQUENCES */
  register const size_t Nprf = ((struct ThreadData*) _Data)->profileCount;
  for (size_t i=Start; i<Stop; ++i) {
    PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(FASTA->DataPtr[i].Sequence));
       
    for (size_t k=0; k<Nprf; ++k) {
      memcpy(Buffer, PFSeq->ProfileIndex, PFSeq->Length);
      Buffer[PFSeq->Length] = '\0';
      PFSequence lPFSeq = {Buffer, PFSeq->Length}; 
      const struct Profile * const restrict prf = prfs[k];
			
			
      unsigned int lHeuristicCutOff;
			if (((struct ThreadData*) _Data)->NormalizedCutoff > 0.0f) {
				if (prf->HeuristicModeIndex < 0) {
					fprintf(stderr, "Cannot compute heuristic cutoff from normalized score, profile %s (%s) is not calibrated!\n", prf->AC_Number, prf->Description);
					exit(1);
				}
				const int CutOff = prf->NormalizedToRaw(((struct ThreadData*) _Data)->NormalizedCutoff,
																								prf->NormalizationCoefs, 0.0f, 0);
				const int itmp = N2R_1((float) CutOff, (const float*) &(prf->NormalizationData.Values[prf->HeuristicModeIndex].RNOP),
																 0.0f, 0);
				lHeuristicCutOff = (itmp >0) ? itmp : 0U;
			}
			else {
				lHeuristicCutOff = prf->HeuristicCutOff;
			}
			
      if (lHeuristicCutOff > 0U) {
				/* Translate first sequence */
				TranslateSequenceToIndex(&lPFSeq, prf->Alphabet_Mapping);
				const unsigned int score = heuristic(TransposeMatch[k], prf->Alphabet_Length, prf->Length, &lPFSeq);
				Array[i*Nprf+k].PrfId = (score >= (unsigned int) lHeuristicCutOff) ? (int) k : -1 ;	
      } 
      else {
				/* No heuristic cutoff, accept all */
				Array[i*Nprf+k].PrfId = (int) k;
      }
      Array[i*Nprf+k].SeqId = (unsigned int) i;
    }
  }

  /* close sequence file */
  UNSET_DATABASE_ACCESS();
  
  /* Free Memory */
  free(SeqData.Data.Memory);
  pthread_exit(0);
};
#endif /* _NEEDS_HEURISTIC_ */

#ifdef _NEEDS_FILTER_
static void *thread_xali1( void * _Data)
{
  Sequence SeqData;
  const struct Profile * const * prfs         = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict FASTA   = ((struct ThreadData*) _Data)->FASTA;
  struct ID * const restrict Array            = ((struct ThreadData*) _Data)->Array;
  const size_t MaxProfileSize                 = ((struct ThreadData*) _Data)->MaxProfileSize;
  PFSequence * PFSeq;

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(FASTA->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fputs("Thread Cannot allocate memory for sequence.\n", stderr);
    return (void*) 1;
  }
  /* Allocate work aligned memory for xali1 */
  int * Work = _mm_malloc((1+MaxProfileSize)*4*sizeof(int)+63,64);
  if (Work == NULL) pthread_exit((void*) 1);
  
  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);
  
  size_t Start  = ((struct ThreadData*) _Data)->start;
  size_t Stop   = ((struct ThreadData*) _Data)->stop;
//   fprintf(stderr,"Thread %lu - %lu\n", Start, Stop); 
  /* LOOPS ON SEQUENCES AND PROFILES */
  for (size_t i=Start; i<Stop; ++i) {
    PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(FASTA->DataPtr[Array[i].SeqId].Sequence));
    
    const struct Profile * const restrict prf = prfs[Array[i].PrfId];
    
    /* Translate first sequence */
    PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
		
    const float lNormalizedCutoff = (((struct ThreadData*) _Data)->NormalizedCutoff > 0.0f) ?
                                    ((struct ThreadData*) _Data)->NormalizedCutoff : prf->NormalizedCutOff;   
    int CutOff;
    if (prf->NormalizationType == GLE_ZSCAVE) {
      NormalizedToRawFunctionPtr NormalizedToRawFunction = prf->NormalizedToRaw;
      const float * const restrict NormCoefs = prf->NormalizationCoefs;
      const float RAVE = ComputeAverageFrequencies(PFSeq, &prf->Average);
      CutOff = NormalizedToRawFunction(lNormalizedCutoff, NormCoefs, RAVE, PFSeq->Length);
    }
    else if ( ((struct ThreadData*) _Data)->NormalizedCutoff > 0.0f ) {
			CutOff = prf->NormalizedToRaw(lNormalizedCutoff, prf->NormalizationCoefs, 0.0f, 0);
		}
    else {
      CutOff = prf->CutOff;
    }
    
    const int score = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, CutOff, false);
    if (score < CutOff) Array[i].PrfId = -1;
  }
  
  /* close sequence file */
  UNSET_DATABASE_ACCESS();
  
  /* Free Memory */
  free(SeqData.Data.Memory);
  _mm_free(Work);
    
  pthread_exit(0);
}
#endif /* _NEEDS_FILTER_ */

#ifdef _NEEDS_ALIGNMENT_
static void *thread_xaliPT( void * _Data)
{ 
  Sequence SeqData;
  const struct Profile * const * prfs         = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict FASTA   = ((struct ThreadData*) _Data)->FASTA; 
  struct ID * const restrict Array            = ((struct ThreadData*) _Data)->Array;
  const size_t MaxProfileSize                 = ((struct ThreadData*) _Data)->MaxProfileSize;
  PFSequence * PFSeq;

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(FASTA->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fputs("Thread cannot allocate memory for sequence.\n", stderr);
    return (void*) 1;
  }
  /* Allocate work aligned memory for xali1 */
  union lScores * const restrict iop          = _mm_malloc((1+MaxProfileSize)*sizeof(union lScores), 16);
  union Positions * const restrict iom        = _mm_malloc((1+MaxProfileSize)*sizeof(union Positions), 16);
  union Positions * const restrict ioi        = _mm_malloc((1+MaxProfileSize)*sizeof(union Positions), 16);
  struct Alignment * const restrict alignment = _mm_malloc(NALI*sizeof(struct Alignment),16);
  _Bool * const restrict Lock                 = _mm_malloc(FASTA->MaxSequenceSize*sizeof(_Bool), 16);
  const size_t MaxAlignmentSize               = FASTA->MaxSequenceSize+(1+MaxProfileSize)+1;
  char * const restrict Sequences             = _mm_malloc(NALI*MaxAlignmentSize*sizeof(char),16);
  if ( iop == NULL || iom == NULL || ioi == NULL || alignment == NULL || Lock == NULL || Sequences == NULL) {
    if (iop) _mm_free(iop);
    if (iom) _mm_free(iom);
    if (ioi) _mm_free(ioi);
    if (alignment) _mm_free(alignment);
    if (Lock) _mm_free(Lock);
    if (Sequences) _mm_free(Sequences);
    pthread_exit((void*)1);
  }
  
  /* Open sequence file */
  PrintInput_t FASTQ;
	FASTQ.DB = FASTA;
	FASTQ.SequenceFile = ((struct ThreadData*) _Data)->SequenceFile;
	SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start = ((struct ThreadData*) _Data)->start;
  size_t Stop  = ((struct ThreadData*) _Data)->stop;

  unsigned int AlignedSeqCounter = 0;
  
  // Allocate on the stack for maximum NALI alignment
  char ** const restrict AlignedSequences = (char **) alloca((NALI)*sizeof(char *));
  AlignedSequences[0] = &Sequences[0];
  for (size_t i=1; i<NALI; ++i) AlignedSequences[i] = &Sequences[i*MaxAlignmentSize];
   
  /* LOOPS ON SEQUENCES AND PROFILES */
  for (size_t i=Start; i<Stop; ++i) {
    PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(FASTA->DataPtr[Array[i].SeqId].Sequence));
    FASTQ.SeqId =  Array[i].SeqId;
		
    const struct Profile * const restrict prf = prfs[Array[i].PrfId];

    /* Translate first sequence */
    PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);

    /* Clear Lock */
    memset(Lock, 0, FASTA->MaxSequenceSize*sizeof(_Bool));
    
		const float lNormalizedCutoff = (((struct ThreadData*) _Data)->NormalizedCutoff > 0.0f) ?
                                    ((struct ThreadData*) _Data)->NormalizedCutoff : prf->NormalizedCutOff;
    int CutOff;
    if (prf->NormalizationType == GLE_ZSCAVE) {
      NormalizedToRawFunctionPtr NormalizedToRawFunction = prf->NormalizedToRaw;
      const float * const restrict NormCoefs = prf->NormalizationCoefs;
      const float RAVE = ComputeAverageFrequencies(PFSeq, &prf->Average);
      CutOff = NormalizedToRawFunction(lNormalizedCutoff, NormCoefs, RAVE, PFSeq->Length);
    }
    else if ( ((struct ThreadData*) _Data)->NormalizedCutoff > 0.0f ) {
			CutOff = prf->NormalizedToRaw(lNormalizedCutoff, prf->NormalizationCoefs, 0.0f, 0);
		}
    else {
      CutOff = prf->CutOff;
    }
    
    // It seems we must have sequence starting from 1 here
    const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
				                       Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false, 
				                       CutOff, NALI); 
			  
    if (nali <= 0) {
      fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu(%u) (nali=%i)!\n%s\n",
	      ((struct ThreadData*) _Data)->threadId, i, Array[i].SeqId, nali, SeqData.Data.Header);
//       pthread_exit((void*)1);          
    }
    
    // Alignement is not filled from start !!!
    for ( int j=1; j<=nali; j++) {
    
      /* Remove lock for aligned sequence generation */
      memset(Lock, 0, FASTA->MaxSequenceSize*sizeof(_Bool));
      
      if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
		    AlignedSequences[j-1], iop, &alignment[j], Lock) < 0 ) {
				fputs("Internal error within xalit!\n", stderr);
				pthread_exit((void*)1);
      }
    }
    pthread_mutex_lock(&PrintLock);
    PrintFunction(prf, (const char ** const restrict) AlignedSequences, &alignment[1], SeqData.Data.Header, PFSeq->Length, 0.0f, nali, &FASTQ);
    pthread_mutex_unlock(&PrintLock);
    
    AlignedSeqCounter += nali;
  }
  

  /* Set the number of aligned sequences */
  ((struct ThreadData*) _Data)->counter = AlignedSeqCounter; 

  /* close sequence file */
  UNSET_DATABASE_ACCESS();
  
  /* Free Memory */
  free(SeqData.Data.Memory);
  _mm_free(iop);
  _mm_free(iom);
  _mm_free(ioi);
  _mm_free(alignment);
  _mm_free(Lock);
  _mm_free(Sequences);

  pthread_exit(0);
}
#endif /* _NEEDS_ALIGNMENT_ */

#ifdef _NEEDS_REGEX_
static void *thread_regex( void * _Data)
{
  Sequence SeqData;
  const struct RegEx * const restrict regex   = ((struct ThreadData*) _Data)->regex;
  const struct Profile * const * prfs         = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict FASTA   = ((struct ThreadData*) _Data)->FASTA;
  PFSequence * PFSeq;
  
  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(FASTA->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fputs("Thread Cannot allocate memory for sequence.\n", stderr);
    return (void*) 1;
  }
  
   /* Allocate the memory to hold the matches */
  const size_t nmatch = regex->maxMatchCount;
#if defined(USE_PCRE)
  int * Matches = (int*) malloc(2*(1+nmatch)*sizeof(int));
#else
  regmatch_t * const restrict Matches = (regmatch_t*) malloc(nmatch*sizeof(regmatch_t));
#endif
  if (Matches == NULL) {
    fputs("Thread Cannot allocate memory for regex match structure.\n", stderr);
    pthread_exit((void*) 1);
  }

  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start              = ((struct ThreadData*) _Data)->start;
  size_t Stop               = ((struct ThreadData*) _Data)->stop;

//   fprintf(stderr,"Thread %lu - %lu\n", Start, Stop);
  /* LOOPS ON SEQUENCES */
  register const size_t Nprf = regex->count;
  
  for (size_t i=Start; i<Stop; ++i) {
    PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(FASTA->DataPtr[i].Sequence));
       
    /* Translate first sequence */
    const char * const CleanSeq = CleanSequence(PFSeq);
      
    for (size_t iPatternPrf=0; iPatternPrf<Nprf; ++iPatternPrf) {
      /* Run the regex engine */
#if defined(USE_PCRE)
      const pcre * const restrict rg = regex->regexCompiled[iPatternPrf];
      int ovector[2*8];
      unsigned int offset = 0;
      const unsigned int len = PFSeq->Length;
      int rc;
      size_t count = 0;
      while (offset < len && (rc = pcre_exec(rg, 0, CleanSeq, len, offset, 0, ovector, 8)) >= 0)
      {
	  for(int k = 0; k < rc; ++k)
	  {
	      if (count < nmatch) {
		  Matches[2*count] = ovector[2*k];
		  Matches[2*count+1] = ovector[2*k+1];
		  ++count;
	      }
	      else {
		  fprintf(stderr, "Warning: maximum number of matches reached for %s with %s\n", 
			prfs[iPatternPrf]->Identification, prfs[iPatternPrf]->Pattern);
	      }
	  }
	  offset = ovector[1];
      }

      Matches[2*count] = -1;
      if (count)
#else
      const regex_t * const restrict rg = &(regex->regexCompiled[iPatternPrf]);
      memset(Matches, 0, nmatch*sizeof(regmatch_t));
      const int res = regexec(rg, CleanSeq, nmatch, Matches, 0);
      if ( res == REG_ESPACE) {
	fputs("Regex ran out of memory\n", stderr);
	exit(1);
      }
      else if (res == 0)
#endif /* defined(USE_PCRE) */
      {
#if !defined(__USE_WINAPI__)
	pthread_mutex_lock(&PrintLock);
#else
	EnterCriticalSection(&PrintLock);
#endif
	PrintRegex(/*prfs[iPatternPrf]->Pattern*/ regex->regexString[iPatternPrf], CleanSeq, Matches, SeqData.Data.Header, PFSeq->Length);
#if !defined(__USE_WINAPI__)
	pthread_mutex_unlock(&PrintLock);
#else
	LeaveCriticalSection(&PrintLock);
#endif
      }
    }
  }
  /* close sequence file */
  UNSET_DATABASE_ACCESS();
  
  /* Free Memory */
  free(SeqData.Data.Memory);
  free(Matches);
    
  pthread_exit(0);
}
#endif /* _NEEDS_REGEX_ */
