/************************************************************************************************
                        PFTOOLS
 ************************************************************************************************
  Oct 3, 2011 threads.h
 ************************************************************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ************************************************************************************************
  THIS FILE IS MEANT TO BE INCLUDED IN THE MAIN PROGRAM DECLARATION FILE
 ************************************************************************************************/
#if defined(__USE_WINAPI__)
#define pthread_exit(x) return (LPDWORD) x
#define THREAD_FUNCTION(F) static LPDWORD WINAPI F(_In_ LPVOID _Data)
#else
#define THREAD_FUNCTION(F) static void* F(void * _Data)
#endif
/*
 ************************************************************************************************
 *                                       STRUCTURES                                             *
 ************************************************************************************************
 */
struct ThreadData {
  const struct Profile * prf;
#ifdef _NEEDS_REGEX_
  const struct RegEx * regex;
#endif
  const DBSequence_t * DB;
  const char * restrict SequenceFile;
#if defined(__USE_MMAP__) && defined(MMAP_DEBUG)
  size_t * maplength;
#endif
  unsigned int * SequenceID;
  int * FilterScores;
  unsigned int * HeuristicScores;
  TransposeMatrix TransposeMatch;
#ifdef _NEEDS_HEURISTIC_CALIBRATION_
  unsigned char * PamTypes;
#endif
  size_t start;
  size_t stop;
  unsigned long counter; /* WARNING: on input in filter phase counter == 0 means stop when cutoff reached, 1 -> compute real filter value */
  size_t threadId;
	enum Strand strand;
  enum Version version;
};
/*
 ************************************************************************************************
 *                                         GLOBALS                                              *
 ************************************************************************************************
 */
static pthread_t *Threads = NULL;
static size_t nThreads = 0;
#if defined(_NEEDS_FILTER_) || defined(_NEEDS_SCORES_ONLY_)
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
#endif
#ifdef _NEEDS_HEURISTIC_CALIBRATION_
static long int Seed = 1234L;
static size_t PamDistanceStart = 1;
static size_t PamDistanceStop = 100;
static size_t PamDistanceStep = 20;
static size_t PamDistanceCounter = 5;	/* Number of pam distances to execute */
static size_t PamSampling = 5;
#endif

#if defined(_NEEDS_ALIGNMENT_)
static PrintFunctionPtr PrintFunction = &PrintDefault;
#if !defined(__USE_WINAPI__)
static pthread_mutex_t PrintLock;
#else
static CRITICAL_SECTION PrintLock;
#endif
#endif

#ifdef _NEEDS_OPTIMAL_ALIGNMENT_
int SatisfyCutoff = false;
#endif

/*
 ************************************************************************************************
 *                                        FUNCTIONS                                             *
 ************************************************************************************************
 */
#ifdef _NEEDS_HEURISTIC_
THREAD_FUNCTION(thread_heuristic)
{
  Sequence SeqData;
  const struct Profile * const restrict prf     = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict DB     = ((struct ThreadData*) _Data)->DB;
  const TransposeMatrix TransposeMatch          = ((struct ThreadData*) _Data)->TransposeMatch;
  unsigned int * const restrict HeuristicScores = ((struct ThreadData*) _Data)->HeuristicScores;
  unsigned int * const restrict SeqID           = ((struct ThreadData*) _Data)->SequenceID;
  PFSequence * PFSeq;

  const HeuristicFunctionPtr heuristic = GetHeuristicVersion(((struct ThreadData*) _Data)->version);

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    for (size_t i=0; i< nThreads; ++i) pthread_cancel(Threads[i]);
    pthread_exit( (void*) 1);
  }

  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start              = ((struct ThreadData*) _Data)->start;
  size_t Stop               = ((struct ThreadData*) _Data)->stop;

  if (SeqID == NULL) {
		switch(((struct ThreadData*) _Data)->strand) {
			case FORWARD:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
					HeuristicScores[i] = heuristic(TransposeMatch, prf->Alphabet_Length, prf->Length, PFSeq);
				}
				break;
			case REVERSE_COMPLEMENT:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);
					HeuristicScores[i] = heuristic(TransposeMatch, prf->Alphabet_Length, prf->Length, PFSeq);
				}
				break;
			case BOTH:
				for (size_t i=Start; i<Stop; i++) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
					HeuristicScores[2*i] = heuristic(TransposeMatch, prf->Alphabet_Length, prf->Length, PFSeq);
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);
					HeuristicScores[2*i+1] = heuristic(TransposeMatch, prf->Alphabet_Length, prf->Length, PFSeq);
				}
		}
  }
  else {
		/* LOOPS ON SEQUENCES */
		for (size_t i=Start; i<Stop; ++i) {
			const unsigned int index = SeqID[i] & 0x3FFFFFFF;
			PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));
			const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;

			PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
			if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

			HeuristicScores[i] = heuristic(TransposeMatch, prf->Alphabet_Length, prf->Length, PFSeq);
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
THREAD_FUNCTION(thread_xali1)
{
  Sequence SeqData;
  const struct Profile * const restrict prf = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict DB    = ((struct ThreadData*) _Data)->DB;
  const unsigned int * const restrict SeqID = ((struct ThreadData*) _Data)->SequenceID;
  int * const restrict Scores               = ((struct ThreadData*) _Data)->FilterScores;
  PFSequence * PFSeq;
  unsigned long TotalSequenceLength = 0L;

  /* Allocate memory to hold sequence and aligned data for xali1 */
  SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  int * Work = _mm_malloc((1+prf->Length)*4*sizeof(int)+63,64);
  if (SeqData.Data.Memory == NULL || Work == NULL) {
     fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    for (size_t i=0; i< nThreads; ++i) pthread_cancel(Threads[i]);
    pthread_exit((void*) 1);
  }
  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start  = ((struct ThreadData*) _Data)->start;
  size_t Stop   = ((struct ThreadData*) _Data)->stop;

  const _Bool RealFilterScore = (((struct ThreadData*) _Data)->counter > 0L ) ? true : false;

  /* Compute all sequences, Sequence index is not given */
  if (SeqID == NULL) {
		switch(((struct ThreadData*) _Data)->strand) {
			case FORWARD:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
					TotalSequenceLength += PFSeq->Length;
					Scores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, RealFilterScore);
				}
				break;
			case REVERSE_COMPLEMENT:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);
					TotalSequenceLength += PFSeq->Length;
					Scores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, RealFilterScore);
				}
				break;
			case BOTH:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
					TotalSequenceLength += PFSeq->Length;
					Scores[2*i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, RealFilterScore);

					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);
					Scores[2*i+1] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, RealFilterScore);
				}
		}
  }
  else { /* Compute all sequences, Sequence index is given */
    /* Do we need to compute the cutoff */
    if (prf->NormalizationType == GLE_ZSCAVE) {
      const float lNormalizedCutoff = prf->NormalizedCutOff;
      NormalizedToRawFunctionPtr NormalizedToRawFunction = prf->NormalizedToRaw;
      const float * const restrict NormCoefs = prf->NormalizationCoefs;

			/* LOOPS ON SEQUENCES */
			for (size_t i=Start; i<Stop; ++i) {
				const unsigned int index = SeqID[i] & 0x3FFFFFFF;
				PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));
				const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;

				PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
				if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

				TotalSequenceLength += PFSeq->Length;
				const float RAVE = ComputeAverageFrequencies(PFSeq, &prf->Average);
				const int CutOff = NormalizedToRawFunction(lNormalizedCutoff, NormCoefs, RAVE, PFSeq->Length);
				Scores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, CutOff, RealFilterScore);
			}
    }
    else {
      const int CutOff = prf->CutOff;
      /* LOOPS ON SEQUENCES */
			for (size_t i=Start; i<Stop; ++i) {
				const unsigned int index = SeqID[i] & 0x3FFFFFFF;
				PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));

				const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;

				PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
				if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

				TotalSequenceLength += PFSeq->Length;
				Scores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, CutOff, RealFilterScore);
			}
    }
  }

  /* close sequence file */
  UNSET_DATABASE_ACCESS();

  /* Free Memory */
  free(SeqData.Data.Memory);
  _mm_free(Work);

  /* Set total sequence length */
  ((struct ThreadData*) _Data)->counter = TotalSequenceLength;

  pthread_exit(0);
}
#endif /* _NEEDS_FILTER_ */

#ifdef _NEEDS_SCORES_ONLY_
THREAD_FUNCTION(thread_scores_only)
{
  Sequence SeqData;
  const struct Profile * const restrict prf     = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict DB     = ((struct ThreadData*) _Data)->DB;
  const TransposeMatrix TransposeMatch          = ((struct ThreadData*) _Data)->TransposeMatch;
  unsigned int * const restrict HeuristicScores = ((struct ThreadData*) _Data)->HeuristicScores;
  int * const restrict FilterScores             = ((struct ThreadData*) _Data)->FilterScores;
  unsigned int * const restrict SeqID           = ((struct ThreadData*) _Data)->SequenceID;
  PFSequence * PFSeq;

  const HeuristicFunctionPtr heuristic = GetHeuristicVersion(((struct ThreadData*) _Data)->version);

  /* Allocate memory to hold sequence and aligned data for xali1*/
  SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  int * Work = _mm_malloc((1+prf->Length)*4*sizeof(int)+63,64);
  if (SeqData.Data.Memory == NULL || Work == NULL) {
     fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    for (size_t i=0; i< nThreads; ++i) pthread_cancel(Threads[i]);
    pthread_exit((void*) 1);
  }
  SeqData.Size = DB->MaxSequenceSize;

  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start              = ((struct ThreadData*) _Data)->start;
  size_t Stop               = ((struct ThreadData*) _Data)->stop;

  if (SeqID) {
    /* LOOPS ON SEQUENCES */
    const size_t AlphabetLength = prf->Alphabet_Length;
    const size_t prfLength = prf->Length;
    for (size_t i=Start; i<Stop; ++i) {
			const unsigned int index = SeqID[i] & 0x3FFFFFFF;
      PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));

			const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;

			PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
			if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

      /* Compute filter */
      FilterScores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, true);

      /* Compute heuristic */
      HeuristicScores[i] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeq);
    }
  }
  else {
    /* LOOPS ON SEQUENCES */
    const size_t AlphabetLength = prf->Alphabet_Length;
    const size_t prfLength = prf->Length;
		switch(((struct ThreadData*) _Data)->strand) {
			case FORWARD:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);

					/* Compute filter */
					FilterScores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, true);

					/* Compute heuristic */
					HeuristicScores[i] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeq);
				}
				break;
			case REVERSE_COMPLEMENT:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);

					/* Compute filter */
					FilterScores[i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, true);

					/* Compute heuristic */
					HeuristicScores[i] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeq);
				}
				break;
			case BOTH:
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);

					/* Compute filter */
					FilterScores[2*i] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, true);

					/* Compute heuristic */
					HeuristicScores[2*i] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeq);
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);

					/* Compute filter */
					FilterScores[2*i+1] = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0.0, true);

					/* Compute heuristic */
					HeuristicScores[2*i+1] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeq);
				}
		}
  }


  /* close sequence file */
  UNSET_DATABASE_ACCESS();

  /* Free Memory */
  free(SeqData.Data.Memory);
  _mm_free(Work);

  pthread_exit(0);
}
#endif /* _NEEDS_SCORES_ONLY_ */

#ifdef _NEEDS_ALIGNMENT_
THREAD_FUNCTION(thread_xaliPT)
{
  Sequence SeqData;
  const struct Profile * const restrict prf = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict DB    = ((struct ThreadData*) _Data)->DB;
  const unsigned int * const restrict SeqID = ((struct ThreadData*) _Data)->SequenceID;
  PFSequence * PFSeq;
	PrintInput_t FASTQ;
	FASTQ.DB = DB;

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
     fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    return (void*) 1;
  }
  /* Allocate work aligned memory for xali1 */
  union lScores * const restrict iop          = _mm_malloc((1+prf->Length)*sizeof(union lScores), 16);
  union Positions * const restrict iom        = _mm_malloc((1+prf->Length)*sizeof(union Positions), 16);
  union Positions * const restrict ioi        = _mm_malloc((1+prf->Length)*sizeof(union Positions), 16);
  struct Alignment * const restrict alignment = _mm_malloc(NALI*sizeof(struct Alignment),16);
  _Bool * const restrict Lock                 = _mm_malloc(DB->MaxSequenceSize*sizeof(_Bool), 16);
  const size_t MaxAlignmentSize               = DB->MaxSequenceSize+(1+prf->Length)+1;
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
	FASTQ.SequenceFile = ((struct ThreadData*) _Data)->SequenceFile;
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start = ((struct ThreadData*) _Data)->start;
  size_t Stop  = ((struct ThreadData*) _Data)->stop;

  unsigned int AlignedSeqCounter = 0;

  // Allocate on the stack for maximum NALI alignment
  char ** const restrict AlignedSequences = (char **) alloca((NALI)*sizeof(char *));
  AlignedSequences[0] = &Sequences[0];
  for (size_t i=1; i<NALI; ++i) AlignedSequences[i] = &Sequences[i*MaxAlignmentSize];

  if (prf->NormalizationType != GLE_ZSCAVE) {
    register const int CutOff = prf->CutOff;
		/* LOOPS ON SEQUENCES */
		for (size_t i=Start; i<Stop; ++i) {
			const unsigned int  index = SeqID[i] & 0x3FFFFFFF;
			PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));
			FASTQ.SeqId = SeqID[i];

			/* Translate first sequence */
			const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;
			PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
			if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

			/* Clear Lock */
			memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

			// It seems we must have sequence starting from 1 here
			const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
							Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
							CutOff, NALI);

			if (nali <= 0) {
				fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu(%u) (nali=%i)!\n%s\n",
					((struct ThreadData*) _Data)->threadId, i, SeqID[i], nali, SeqData.Data.Header);
				((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
				pthread_exit((void*)1);
			}

			// Alignement is not filled from start !!!
			for ( int j=1; j<=nali; j++) {
				/* Remove lock for aligned sequence generation */
				memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));
				if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
					AlignedSequences[j-1], iop, &alignment[j], Lock) < 0 ) {
					fputs("Internal error within xalit!\n", stderr);
					((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
					pthread_exit((void*)1);
				}
				++AlignedSeqCounter;
			}
#if !defined(__USE_WINAPI__)
			pthread_mutex_lock(&PrintLock);
#else
			EnterCriticalSection(&PrintLock);
#endif
			PrintFunction(prf, (const char**) AlignedSequences, &alignment[1], SeqData.Data.Header, PFSeq->Length, 0.0f, nali, &FASTQ);
#if !defined(__USE_WINAPI__)
			pthread_mutex_unlock(&PrintLock);
#else
			LeaveCriticalSection(&PrintLock);
#endif
		}
  }
  else {
    const float lNormalizedCutoff = prf->NormalizedCutOff;
    NormalizedToRawFunctionPtr NormalizedToRawFunction = prf->NormalizedToRaw;
    const float * const restrict NormCoefs = prf->NormalizationCoefs;

		/* LOOPS ON SEQUENCES */
		for (size_t i=Start; i<Stop; ++i) {
			const unsigned int  index = SeqID[i] & 0x3FFFFFFF;
			PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[index].Sequence));
			FASTQ.SeqId = SeqID[i];

			/* Translate first sequence */
			const unsigned char * AlphabetPtr = (SeqID[i] & 0x80000000) ? prf->Complement_Alphabet_Mapping : prf->Alphabet_Mapping;
			PFSeq = TranslateSequenceToIndex(PFSeq, AlphabetPtr);
			if (SeqID[i] & 0x40000000) ReverseTranslatedSequence(PFSeq);

			const float RAVE = ComputeAverageFrequencies(PFSeq, &prf->Average);
			const int CutOff = NormalizedToRawFunction(lNormalizedCutoff, NormCoefs, RAVE, PFSeq->Length);

			/* Clear Lock */
			memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

			// It seems we must have sequence starting from 1 here
			const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
							Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
							CutOff, NALI);

			if (nali <= 0) {
				fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu(%u) (nali=%i)!\n%s\n",
					((struct ThreadData*) _Data)->threadId, i, SeqID[i], nali, SeqData.Data.Header);
				((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
				pthread_exit((void*)1);
			}

			// Alignement is not filled from start !!!
			for ( int j=1; j<=nali; j++) {
				/* Remove lock for aligned sequence generation */
				memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));
				if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
					AlignedSequences[j-1], iop, &alignment[j], Lock) < 0 ) {
						fputs("Internal error within xalit!\n", stderr);
						((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
						pthread_exit((void*)1);
				}
				++AlignedSeqCounter;
			}
#if !defined(__USE_WINAPI__)
			pthread_mutex_lock(&PrintLock);
#else
			EnterCriticalSection(&PrintLock);
#endif
			PrintFunction(prf, (const char**) AlignedSequences, &alignment[1], SeqData.Data.Header, PFSeq->Length, RAVE, nali, &FASTQ);
#if !defined(__USE_WINAPI__)
			pthread_mutex_unlock(&PrintLock);
#else
			LeaveCriticalSection(&PrintLock);
#endif
		}
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

#ifdef _NEEDS_OPTIMAL_ALIGNMENT_
THREAD_FUNCTION(thread_optimal)
{
	Sequence SeqData;
	const struct Profile * const restrict prf   = ((struct ThreadData*) _Data)->prf;
	const DBSequence_t * const restrict DB = ((struct ThreadData*) _Data)->DB;
	PFSequence * PFSeq;
	PrintInput_t FASTQ;
	FASTQ.DB = DB;

	/* Allocate memory to hold sequence */
	SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
	if (SeqData.Data.Memory == NULL) {
		 fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
		return (void*) 1;
	}

	/* Allocate work aligned memory for xali1 */
	int * const restrict Work                   = _mm_malloc((1+prf->Length)*4*sizeof(int)+63,64);
	union lScores * const restrict iop          = _mm_malloc((1+prf->Length)*sizeof(union lScores), 16);
	union Positions * const restrict iom        = _mm_malloc((1+prf->Length)*sizeof(union Positions), 16);
	union Positions * const restrict ioi        = _mm_malloc((1+prf->Length)*sizeof(union Positions), 16);
	struct Alignment * const restrict alignment = _mm_malloc(NALI*sizeof(struct Alignment),16);
	_Bool * const restrict Lock                 = _mm_malloc(DB->MaxSequenceSize*sizeof(_Bool), 16);
	const size_t MaxAlignmentSize               = DB->MaxSequenceSize+(1+prf->Length)+1;
	char * const restrict Sequences             = _mm_malloc(MaxAlignmentSize*sizeof(char),16);
	if ( iop == NULL || iom == NULL || ioi == NULL || alignment == NULL || Lock == NULL || Sequences == NULL || Work == NULL) {
		if (iop) _mm_free(iop);
		if (iom) _mm_free(iom);
		if (ioi) _mm_free(ioi);
		if (alignment) _mm_free(alignment);
		if (Lock) _mm_free(Lock);
		if (Sequences) _mm_free(Sequences);
		if (Work) _mm_free(Work);
		pthread_exit((void*)1);
	}

	/* Open sequence file */
	SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

	size_t Start = ((struct ThreadData*) _Data)->start;
	size_t Stop  = ((struct ThreadData*) _Data)->stop;

	unsigned int AlignedSeqCounter = 0;

	// Allocate on the stack for maximum NALI alignment
	{
		char ** const restrict AlignedSequences = (char **) alloca(1*sizeof(char *));
		AlignedSequences[0] = &Sequences[0];

		switch(((struct ThreadData*) _Data)->strand) {
			case FORWARD:
				FASTQ.SeqId = 0U;
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);

					/* Get the optimal alignment score */
					const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

					/* Clear Lock */
					memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

					// It seems we must have sequence starting from 1 here
					const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																			Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																		CutOff, NALI);

					if (nali <= 0) {
						fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
										((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
						((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
						pthread_exit((void*)1);
					}

					// Get best alignment
					int best = 1;
					int BestScore = alignment[1].Score;
					for (int i=2; i<=nali; i++) {
						if (BestScore < alignment[i].Score ) {
							BestScore = alignment[i].Score;
							best = i;
						}
					}

					if (BestScore < prf->CutOff && SatisfyCutoff) continue;

					// Alignement is not filled from start !!!
					{
						/* Remove lock for aligned sequence generation */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
							AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
							fputs("Internal error within xalit!\n", stderr);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}
						++AlignedSeqCounter;
					}

					#if !defined(__USE_WINAPI__)
					pthread_mutex_lock(&PrintLock);
					#else
					EnterCriticalSection(&PrintLock);
					#endif
					PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
					#if !defined(__USE_WINAPI__)
					pthread_mutex_unlock(&PrintLock);
					#else
					LeaveCriticalSection(&PrintLock);
					#endif
				}
				break;
			case REVERSE:
				FASTQ.SeqId = 0x40000000;
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);

					/* Get the optimal alignment score */
					const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

					/* Clear Lock */
					memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

					// It seems we must have sequence starting from 1 here
					const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																			Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																		CutOff, NALI);

					if (nali <= 0) {
						fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
										((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
						((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
						pthread_exit((void*)1);
					}

					// Get best alignment
					int best = 1;
					int BestScore = alignment[1].Score;
					for (int i=2; i<=nali; i++) {
						if (BestScore < alignment[i].Score ) {
							BestScore = alignment[i].Score;
							best = i;
						}
					}

					if (BestScore < prf->CutOff && SatisfyCutoff) continue;

					// Alignement is not filled from start !!!
					{
						/* Remove lock for aligned sequence generation */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
							AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
							fputs("Internal error within xalit!\n", stderr);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}
						++AlignedSeqCounter;
					}

					#if !defined(__USE_WINAPI__)
					pthread_mutex_lock(&PrintLock);
					#else
					EnterCriticalSection(&PrintLock);
					#endif
					PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
					#if !defined(__USE_WINAPI__)
					pthread_mutex_unlock(&PrintLock);
					#else
					LeaveCriticalSection(&PrintLock);
					#endif
				}
				break;
			case COMPLEMENT:
				FASTQ.SeqId = 0x80000000;
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);

					/* Get the optimal alignment score */
					const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

					/* Clear Lock */
					memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

					// It seems we must have sequence starting from 1 here
					const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																			Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																		CutOff, NALI);

					if (nali <= 0) {
						fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
										((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
						((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
						pthread_exit((void*)1);
					}

					// Get best alignment
					int best = 1;
					int BestScore = alignment[1].Score;
					for (int i=2; i<=nali; i++) {
						if (BestScore < alignment[i].Score ) {
							BestScore = alignment[i].Score;
							best = i;
						}
					}

					if (BestScore < prf->CutOff && SatisfyCutoff) continue;

					// Alignement is not filled from start !!!
					{
						/* Remove lock for aligned sequence generation */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
							AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
							fputs("Internal error within xalit!\n", stderr);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}
						++AlignedSeqCounter;
					}

					#if !defined(__USE_WINAPI__)
					pthread_mutex_lock(&PrintLock);
					#else
					EnterCriticalSection(&PrintLock);
					#endif
					PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
					#if !defined(__USE_WINAPI__)
					pthread_mutex_unlock(&PrintLock);
					#else
					LeaveCriticalSection(&PrintLock);
					#endif
				}
				break;
			case REVERSE_COMPLEMENT:
				FASTQ.SeqId = 0xC0000000;
				for (size_t i=Start; i<Stop; ++i) {
					PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
					ReverseTranslatedSequence(PFSeq);

					/* Get the optimal alignment score */
					const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

					/* Clear Lock */
					memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

					// It seems we must have sequence starting from 1 here
					const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																			Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																		CutOff, NALI);

					if (nali <= 0) {
						fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
										((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
						((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
						pthread_exit((void*)1);
					}

					// Get best alignment
					int best = 1;
					int BestScore = alignment[1].Score;
					for (int i=2; i<=nali; i++) {
						if (BestScore < alignment[i].Score ) {
							BestScore = alignment[i].Score;
							best = i;
						}
					}

					if (BestScore < prf->CutOff && SatisfyCutoff) continue;

					// Alignement is not filled from start !!!
					{
						/* Remove lock for aligned sequence generation */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
							AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
							fputs("Internal error within xalit!\n", stderr);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}
						++AlignedSeqCounter;
					}

					#if !defined(__USE_WINAPI__)
					pthread_mutex_lock(&PrintLock);
					#else
					EnterCriticalSection(&PrintLock);
					#endif
					PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
					#if !defined(__USE_WINAPI__)
					pthread_mutex_unlock(&PrintLock);
					#else
					LeaveCriticalSection(&PrintLock);
					#endif
				}
				break;
			case BOTH:
				FASTQ.SeqId = 0U;
				for (size_t i=Start; i<Stop; ++i) {
					{
						PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

						/* Translate first sequence */
						PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping);

						/* Get the optimal alignment score */
						const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

						/* Clear Lock */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						// It seems we must have sequence starting from 1 here
						const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																				Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																			CutOff, NALI);

						if (nali <= 0) {
							fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
											((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}

						// Get best alignment
						int best = 1;
						int BestScore = alignment[1].Score;
						for (int i=2; i<=nali; i++) {
							if (BestScore < alignment[i].Score ) {
								BestScore = alignment[i].Score;
								best = i;
							}
						}

						if (BestScore < prf->CutOff && SatisfyCutoff) continue;

						// Alignement is not filled from start !!!
						{
							/* Remove lock for aligned sequence generation */
							memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

							if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
								AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
								fputs("Internal error within xalit!\n", stderr);
								((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
								pthread_exit((void*)1);
							}
							++AlignedSeqCounter;
						}

						#if !defined(__USE_WINAPI__)
						pthread_mutex_lock(&PrintLock);
						#else
						EnterCriticalSection(&PrintLock);
						#endif
						PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
						#if !defined(__USE_WINAPI__)
						pthread_mutex_unlock(&PrintLock);
						#else
						LeaveCriticalSection(&PrintLock);
						#endif
					}
					FASTQ.SeqId = 0xC0000000;
					{
						PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

						/* Translate first sequence */
						PFSeq = TranslateSequenceToIndex(PFSeq, prf->Complement_Alphabet_Mapping);
						ReverseTranslatedSequence(PFSeq);

						/* Get the optimal alignment score */
						const int CutOff = xali1_ptr(prf, PFSeq->ProfileIndex, Work, 0, PFSeq->Length, 0, true);

						/* Clear Lock */
						memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

						// It seems we must have sequence starting from 1 here
						const int nali = xalip_ptr(prf, PFSeq->ProfileIndex, iop, iom, ioi, 1, PFSeq->Length, alignment,
																				Lock, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], false,
																			CutOff, NALI);

						if (nali <= 0) {
							fprintf(stderr,"Thread %lu : Internal error xalip reported no possible alignment for sequence %lu (nali=%i)!\n%s\n",
											((struct ThreadData*) _Data)->threadId, i, nali, SeqData.Data.Header);
							((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
							pthread_exit((void*)1);
						}

						// Get best alignment
						int best = 1;
						int BestScore = alignment[1].Score;
						for (int i=2; i<=nali; i++) {
							if (BestScore < alignment[i].Score ) {
								BestScore = alignment[i].Score;
								best = i;
							}
						}

						if (BestScore < prf->CutOff && SatisfyCutoff) continue;

						// Alignement is not filled from start !!!
						{
							/* Remove lock for aligned sequence generation */
							memset(Lock, 0, DB->MaxSequenceSize*sizeof(_Bool));

							if (xalit_ptr(prf, prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1], 1, PFSeq->Length, &(PFSeq->ProfileIndex[0]),
								AlignedSequences[0], iop, &alignment[best], Lock) < 0 ) {
								fputs("Internal error within xalit!\n", stderr);
								((struct ThreadData*) _Data)->counter = AlignedSeqCounter;
								pthread_exit((void*)1);
							}
							++AlignedSeqCounter;
						}

						#if !defined(__USE_WINAPI__)
						pthread_mutex_lock(&PrintLock);
						#else
						EnterCriticalSection(&PrintLock);
						#endif
						PrintFunction(prf, (const char**) AlignedSequences, &alignment[best], SeqData.Data.Header, PFSeq->Length, 0.0f, 1, &FASTQ);
						#if !defined(__USE_WINAPI__)
						pthread_mutex_unlock(&PrintLock);
						#else
						LeaveCriticalSection(&PrintLock);
						#endif
					}
				}
		}
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
	_mm_free(Work);

	pthread_exit(0);
}
#endif /* _NEEDS_OPTIMAL_ALIGNMENT_ */

#ifdef _NEEDS_HEURISTIC_CALIBRATION_
#include "random.h"
extern void mutate(struct Random * const restrict Generator, const struct Profile * const restrict prf,
		   const PFSequence * const restrict inSequence, Sequence * const restrict outSequence,
	           const size_t MutationNumber);

THREAD_FUNCTION(thread_heuristic_calibration)
{
  Sequence SeqDataIn, SeqDataOut;
  struct Random Generator;
  const struct Profile * const restrict prf   = ((struct ThreadData*) _Data)->prf;
  const DBSequence_t * const restrict DB   = ((struct ThreadData*) _Data)->DB;
  const TransposeMatrix TransposeMatch        = ((struct ThreadData*) _Data)->TransposeMatch;
  unsigned int * const HeuristicScores        = ((struct ThreadData*) _Data)->HeuristicScores;
  int * const FilterScores                    = ((struct ThreadData*) _Data)->FilterScores;
  PFSequence * PFSeqIn, * restrict PFSeqOut = &(SeqDataOut.ProfileData);

  const HeuristicFunctionPtr heuristic = GetHeuristicVersion(((struct ThreadData*) _Data)->version);

  /* Initialize random generator */
  const long int MySeed = Seed + ((long int) ((struct ThreadData*) _Data)->threadId);
  InitializeGenerator(MySeed, &Generator);

  /* Allocate memory to hold sequence */
  SeqDataIn.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  SeqDataOut.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  if (SeqDataOut.Data.Memory == NULL || SeqDataIn.Data.Memory == NULL) {
     fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    return (void*) 1;
  }
  SeqDataOut.Size = DB->MaxSequenceSize;
  SeqDataIn.Size = DB->MaxSequenceSize;

  /* Allocate work aligned memory for xali1 */
  int * Work = _mm_malloc((1+prf->Length)*4*sizeof(int)+63,64);
  if (Work == NULL) {
    fprintf(stderr, "Thread %lu could not allocate memory\n", ((struct ThreadData*) _Data)->threadId);
    fflush(stderr);
    pthread_exit((void*) 1);
  }

  /* Open sequence file*/
  SETUP_DATABASE_ACCESS(((struct ThreadData*) _Data)->SequenceFile);

  size_t Start              = ((struct ThreadData*) _Data)->start;
  size_t Stop               = ((struct ThreadData*) _Data)->stop;
  //const unsigned int CutOff = prf->HeuristicCutOff;

//   fprintf(stderr,"Thread %lu - %lu\n", Start, Stop);
	const size_t PamMemoryStride = 1 + PamSampling*PamDistanceCounter;
	{
		const size_t PAMOffset = Start*PamMemoryStride;
		/* LOOPS ON SEQUENCES */
		unsigned int * restrict HeuristicScoresPtr = &HeuristicScores[PAMOffset];
		int * restrict FilterScoresPtr = &FilterScores[PAMOffset];
		unsigned char * restrict PamTypesPtr = &(((struct ThreadData*) _Data)->PamTypes)[PAMOffset];
		const size_t AlphabetLength = prf->Alphabet_Length;
		const size_t prfLength = prf->Length;
		for (size_t i=Start; i<Stop; ++i) {
      PFSeqIn = GET_DATABASE_SEQUENCE(&SeqDataIn, &(DB->DataPtr[i].Sequence));

			/* Translate first sequence */
			PFSeqIn = TranslateSequenceToIndex(PFSeqIn, prf->Alphabet_Mapping);

			/* Treat original sequence */
			FilterScoresPtr[0]    = xali1_ptr(prf, PFSeqIn->ProfileIndex, Work, 0, PFSeqIn->Length, 0.0, true);
			HeuristicScoresPtr[0] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeqIn);
			PamTypesPtr[0]        = 0;

			HeuristicScoresPtr += 1;
			FilterScoresPtr    += 1;
			PamTypesPtr        += 1;

			/* Loop on Pam distances */
			unsigned char ipamtype = PamDistanceStart;
			for (size_t Pam=PamDistanceStart; Pam<=PamDistanceStop; Pam+=PamDistanceStep) {
				/* Loop on the sequence mutation transformation */
				for (size_t iTransform=0; iTransform<PamSampling; ++iTransform) {
					/* Transform the sequence */
					mutate(&Generator, prf, PFSeqIn, &SeqDataOut, Pam);
					PFSeqOut = &(SeqDataOut.ProfileData);

					/* Compute filter */
					FilterScoresPtr[iTransform] = xali1_ptr(prf, PFSeqOut->ProfileIndex, Work, 0, PFSeqOut->Length, 0.0, true);

					/* Compute heuristic */
					HeuristicScoresPtr[iTransform] = heuristic(TransposeMatch, AlphabetLength, prfLength, PFSeqOut);

					/* Set pam type */
					PamTypesPtr[iTransform] = ipamtype;

				}
				ipamtype += PamDistanceStep;
				HeuristicScoresPtr += PamSampling;
				FilterScoresPtr    += PamSampling;
				PamTypesPtr        += PamSampling;
			}
		}
	}

  /* close sequence file */
  UNSET_DATABASE_ACCESS();

  /* Free Memory */
  free(SeqDataIn.Data.Memory);
  free(SeqDataOut.Data.Memory);
  _mm_free(Work);

  pthread_exit(0);
}
#endif /* _NEEDS_HEURISTIC_CALIBRATION_ */

#ifdef _NEEDS_REGEX_
THREAD_FUNCTION(thread_regex)
{
  Sequence SeqData;
  const struct RegEx * const restrict regex = ((struct ThreadData*) _Data)->regex;
  const DBSequence_t * const restrict DB = ((struct ThreadData*) _Data)->DB;
  PFSequence * PFSeq;

  /* Allocate memory to hold sequence */
  SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
  if (SeqData.Data.Memory == NULL) {
    fprintf(stderr, "Thread %s cannot allocate memory for sequence.\n", __FUNCTION__);
    pthread_exit((void*) 1);
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

  size_t Start  = ((struct ThreadData*) _Data)->start;
  size_t Stop   = ((struct ThreadData*) _Data)->stop;

  /*
   *  SEARCH IN SEQUENCE
   */
  size_t TotalCount=0;
#if defined(USE_PCRE)
  const pcre * const restrict rg = regex->regexCompiled[0];
#else
  const regex_t * const restrict rg = &regex->regexCompiled[0];
#endif

  if (((struct ThreadData*) _Data)->counter == 0L) {
    /* LOOPS ON SEQUENCES */
    for (size_t i=Start; i<Stop; ++i) {
      PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

      /* Translate first sequence */
      const char * const CleanSeq = CleanSequence(PFSeq);

      /* Run the regex engine */
#if defined(USE_PCRE)
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
					fputs("Warning: maximum number of matches reached,"
					"you may miss some unless you extend using --max-regex-match.\n",
						stderr);
					//goto OUT;
						}
				}
				offset = ovector[1];
      }

//OUT: ;

      TotalCount += count;
      Matches[2*count] = -1;
      if (count)
#else
      memset(Matches, 0, nmatch*sizeof(regmatch_t));
      const int res = regexec(rg, CleanSeq, nmatch, Matches, 0);
      if ( res == REG_ESPACE) {
				fputs("Regex ran out of memory\n", stderr);
				exit(1);
      }
      else if (res == 0)
#endif
      {
#if !defined(__USE_WINAPI__)
				pthread_mutex_lock(&PrintLock);
#else
				EnterCriticalSection(&PrintLock);
#endif
				PrintRegex(regex->regexString[0], CleanSeq, Matches, SeqData.Data.Header, PFSeq->Length);
#if !defined(__USE_WINAPI__)
				pthread_mutex_unlock(&PrintLock);
#else
				LeaveCriticalSection(&PrintLock);
#endif
      }
    }
  }
  else {
     /* LOOPS ON SEQUENCES */
    for (size_t i=Start; i<Stop; ++i) {
      PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));

      /* Put an end to the header string */
      const char * const EndedHeader = SeqData.Data.Header;

      /* Run the regex engine */
#if defined(USE_PCRE)
      int ovector[2*8];
      unsigned int offset = 0;
      const unsigned int len = (unsigned int) ((uintptr_t) &PFSeq->ProfileIndex[0] - (uintptr_t) &EndedHeader[0]);
      int rc;
      size_t count = 0;
      while (offset < len && (rc = pcre_exec(rg, 0, EndedHeader, len, offset, 0, ovector, 8)) >= 0)
      {
				for(int k = 0; k < rc; ++k)
				{
						if (count < nmatch) {
					Matches[2*count] = ovector[2*k];
					Matches[2*count+1] = ovector[2*k+1];
					++count;
						}
						else {
					fputs("Warning: maximum number of matches reached,"
					"you may miss some unless you extend using --max-regex-match.\n",
						stderr);
					//goto OUT;
						}
				}
				offset = ovector[1];
      }

//OUT: ;

      TotalCount += count;
      Matches[2*count] = -1;
      if (count)
#else
      memset(Matches, 0, nmatch*sizeof(regmatch_t));
      const int res = regexec(rg, EndedHeader, nmatch, Matches, 0);
      if ( res == REG_ESPACE) {
				fputs("Regex ran out of memory\n", stderr);
				exit(1);
      }
      else if (res == 0)
#endif
      {
#if !defined(__USE_WINAPI__)
				pthread_mutex_lock(&PrintLock);
#else
				EnterCriticalSection(&PrintLock);
#endif
				fprintf(stdout, "%s\n%s\n", EndedHeader, &PFSeq->ProfileIndex[0]);
#if !defined(__USE_WINAPI__)
				pthread_mutex_unlock(&PrintLock);
#else
				LeaveCriticalSection(&PrintLock);
#endif
      }
    }
  }
  ((struct ThreadData*) _Data)->counter = (unsigned long) TotalCount;

  /* close sequence file */
  UNSET_DATABASE_ACCESS();

  /* Free Memory */
  free(SeqData.Data.Memory);
  free(Matches);
  pthread_exit(0);
};

#endif /* _NEEDS_REGEX_ */
