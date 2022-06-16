/*******************************************************
                        PFTOOLS
 *******************************************************
  Oct 12, 2012 xali1_print_sse41.c
 *******************************************************
 (C) 2012 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <smmintrin.h>
#include "../include/pfProfile.h"
#include "sse41_inline_fcts.h"
#include "../include/pfplot_functions.h"

#define MAX(a,b) ((a>b) ? a : b)

union sIOP {
    struct {
      int Match;
      int Insertion;
    } Element;
    __m64 mm;
};

void xali1_print_sse41(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
		       int * const WORK, union lScores * const matrix, const size_t BSEQ, const size_t LSEQ,
		       const FuncPtr operation)
/*
 * WARNING: for SSE version, WORK should be aligned on cache size (64b) and 4 times the (profile size + 1)*sizeof(int)
 *	    + 63 to align to cache line
 *           matrix should be of size at least (profile size + 1)*(sequence size) and aligned on 16b.
 */
{
  int KOPD;
  const union sIOP * restrict IOP_R;
  union sIOP * restrict IOP_W = (union sIOP*) WORK;
  const __m128i __MatchMask     = _mm_set1_epi32(PRIORITY_MATCH);
  const __m128i __InsertionMask = _mm_set1_epi32(PRIORITY_INSERTION);
  const __m128i __DeletionMask  = _mm_set1_epi32(PRIORITY_DELETION);
  const __m128i __ExtraMask     = _mm_set1_epi32(PRIORITY_EXTRA);
  const __m128i __ClearMask     = _mm_set1_epi32(0xC0000000); /* two major bits on */

  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                                       = prf->Scores.Match.AlignStep;
  const size_t prfLength = prf->Length;

  /* Set matrix ptr according to BSEQ */
  __m128i * const restrict MatrixPtr = (BSEQ == 0) ? &matrix[1+prfLength].xmm : &matrix[BSEQ].xmm;

  /* NOTE: The following part could be replaced and performed only once for a profile as it
   *       is profile dependent. Nevertheless it does a good job loading Match and Transition
   *       matrices into the cache hierarchy.
   */

  /*
   * Initialize Insertion and Match Entrance Line using FirstSequenceProtein
   */
  {
    register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    /*
     * PROFILE COLUMN 0 entrance
     */
    IOP_W[0].Element.Match     = (int) FirstSequenceProtein[0].To[MATCH];
    IOP_W[0].Element.Insertion = (int) FirstSequenceProtein[0].To[INSERTION];
    KOPD                       = (int) FirstSequenceProtein[0].To[DELETION];

#ifdef TAG
    {
	__m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));
	__FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, 2);
	__FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
	// Store all scores to matrix
	operation(__FirstSequenceProtein, MatrixPtr - (1+prfLength), 0, 0, 1+prfLength);
    }
#endif

    FirstSequenceProtein++;
    register const TransitionScores (* restrict pTransitions) = &Transitions[1];
    register union sIOP * restrict pIOP = &IOP_W[1];

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      register const int KD = KOPD + (int) *lMatch;
      lMatch += AlignStep;

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and  Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(pTransitions->From[DELETION]));

      // Add KD to Transitions
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      // Move to next profile transitions
      pTransitions++;

      // Load FirstSequenceProtein
      __m128i __FirstSequenceProtein = LoadStoredIntegerVector(&(FirstSequenceProtein[0]));

      // Move to next profile First Sequence
      FirstSequenceProtein++;

#ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __FirstSequenceProtein = _mm_slli_epi32(__FirstSequenceProtein, 2);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
      __FirstSequenceProtein = _mm_or_si128(__FirstSequenceProtein, __ExtraMask);
#endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__TransitionsD, __FirstSequenceProtein);

      // Store all scores to matrix
      operation(__max, MatrixPtr - (1+prfLength), 0, iprf, 1+prfLength);

#ifdef TAG
      // Clean extra bits
      __m128i __sign = _mm_cmplt_epi32(__max, _mm_setzero_si128());
      __max = _mm_srai_epi32(__max, 2);
      __max = _mm_or_si128(__max,_mm_and_si128(__sign, __ClearMask));
#endif

      // Store IOPI and IOPM
      StoreMatchInsertion( &(pIOP->mm), (__m128) __max);
      pIOP++;

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max, DELETION);
    }
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = (union sIOP*) (((uintptr_t) &WORK[2*(prf->Length+1)] + 63) & ~63);

  /*
   * LOOP THROUGH THE SEQUENCE STRING
   */
  for ( size_t iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
    register const size_t j1 = (size_t) Sequence[iseq];
    int KOPM = IOP_R[0].Element.Match;
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    /*
     * PROFILE COLUMN 0 entrance
     */
    {
      register const int KI = IOP_R[0].Element.Insertion + (int) lInsertion[j1];

      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[0].From[INSERTION]));
     // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

       // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[0].From[EXTRA]));

      // Insert lScore into __TransitionsX not necessary as profile loading store NLOW there automatically
      //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);

#ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
#endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__TransitionsI, __TransitionsX);

      // Store all scores to matrix
      operation(__max, MatrixPtr, iseq, 0, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

#ifdef TAG
      // Clean extra bits
      __m128i __sign = _mm_cmplt_epi32(__max, _mm_setzero_si128());
      __max = _mm_srai_epi32(__max, 2);
      __max = _mm_or_si128(__max,_mm_and_si128(__sign, __ClearMask));
#endif
      // Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __max);

      // Store KOPD
      KOPD = _mm_extract_epi32(__max, DELETION);

      // Backup new score to xmm register
      //lScore = _mm_extract_epi32(__max, DUMMY);
    }

    lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;

    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      const int KM = KOPM                          + (int) lMatch[j1];
      const int KI = IOP_R[iprf].Element.Insertion + (int) lInsertion[j1];
      const int KD = KOPD                          + (int) lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;

      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

#ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsM = _mm_slli_epi32(__TransitionsM, 2);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
#endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
       // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

#ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
#endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, iseq, iprf, 1+prfLength); //_mm_store_si128(&pmatrix[iprf-1], __max1);

#ifdef TAG
      // Clean extra bits
      __m128i __sign = _mm_cmplt_epi32(__max1, _mm_setzero_si128());
      __max1 = _mm_srai_epi32(__max1, 2);
      __max1 = _mm_or_si128(__max1,_mm_and_si128(__sign, __ClearMask));
#endif

      // Store IOPI and IOPM
      StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max1);


      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);

//       lScore = _mm_extract_epi32(__max1, DUMMY);

    }

    // Swap Read and Write pointers
    const register union sIOP * const ptr = IOP_W;
    IOP_W = (union sIOP *) IOP_R;
    IOP_R = ptr;
  }

  /*
   * Last position on the Sequence using LastSequenceProtein
   */
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    /*
     * PROFILE COLUMN 0 entrance
     */
    int KOPM     = IOP_R[0].Element.Match;
    int KI       = IOP_R[0].Element.Insertion + (int) lInsertion[j1];

    KOPD = MAX( KI + (int) Transitions[0].Element[_ID], (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;

    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;
#ifdef TAG
    __m128i __Scores = _mm_set1_epi32(NLOW<<2);
#else
    __m128i __Scores = _mm_set1_epi32(NLOW);
#endif
    /*
     * LOOP THROUGH THE REST OF THE PROFILE
     */
    for (size_t iprf=1; iprf<=prfLength; ++iprf) {
      const int KM = KOPM                          + lMatch[j1];
      KI           = IOP_R[iprf].Element.Insertion + lInsertion[j1];
      const int KD = KOPD                          + lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;
#ifdef ONLY_NECESSARY
#  ifdef TAG
      int tM = ((KM + (int) Transitions[iprf].Element[_MD]) << 2 ) | PRIORITY_MATCH;
      int tX = ((     (int) Transitions[iprf].Element[_XD]) << 2 ) | PRIORITY_EXTRA;
      int tI = ((KI + (int) Transitions[iprf].Element[_ID]) << 2 ) | PRIORITY_INSERTION;
      int tD = ((KD + (int) Transitions[iprf].Element[_DD]) << 2 ) | PRIORITY_DELETION;
      int tIOPD1 = MAX(tM, tX);
      int tIOPD2 = MAX(tI, tD);
#  else
      const int tIOPD1 = MAX( KM + (int) Transitions[iprf].Element[_MD],      (int) Transitions[iprf].Element[_XD] );
      const int tIOPD2 = MAX( KI + (int) Transitions[iprf].Element[_ID], KD + (int) Transitions[iprf].Element[_DD] );
#  endif
      KOPD     = MAX( tIOPD1, tIOPD2);
      __Scores = _mm_insert_epi32(__Scores, KOPD, DELETION);

#  ifdef TAG
      tM = ((KM + (int) LastSequenceProtein[iprf].From[MATCH])     << 2 ) | PRIORITY_MATCH;
      tI = ((KI + (int) LastSequenceProtein[iprf].From[INSERTION]) << 2 ) | PRIORITY_INSERTION;
      tD = ((KD + (int) LastSequenceProtein[iprf].From[DELETION])  << 2 ) | PRIORITY_DELETION;
      const int tIOPT1 = MAX(tM, tI);
      const int lScore = MAX(tIOPT1, tD);
#  else
      const int tIOPT1 = MAX(KM + (int) LastSequenceProtein[iprf].From[MATCH] , KI + (int) LastSequenceProtein[iprf].From[INSERTION] );
      const int lScore = MAX( tIOPT1                                          , KD + (int) LastSequenceProtein[iprf].From[DELETION]  );
#  endif
      __Scores = _mm_insert_epi32(__Scores, lScore, DUMMY);

      operation(__Scores, MatrixPtr, LSEQ-1, iprf, prfLength+1); //_mm_store_si128(&pmatrix[iprf-1], __Scores);
#  ifdef TAG
      // Clean KOPD extra 2 bits
      KOPD = (KOPD < 0) ? (KOPD >> 2) | 0xC0000000 : (KOPD >> 2);
#  endif
#else
      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Insert LastProteinSequence
      __TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastSequenceProtein[iprf].From[MATCH], DUMMY);
      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Insert LastProteinSequence
      __TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastSequenceProtein[iprf].From[INSERTION], DUMMY);
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

#  ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsM = _mm_slli_epi32(__TransitionsM, 2);
      __TransitionsI = _mm_slli_epi32(__TransitionsI, 2);
      __TransitionsM = _mm_or_si128(__TransitionsM, __MatchMask);
      __TransitionsI = _mm_or_si128(__TransitionsI, __InsertionMask);
#  endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX not necessary as profile loading should have already done it
      // __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastSequenceProtein[iprf].From[DELETION], DUMMY);
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

#  ifdef TAG
      // Paste index in lowest 2 bits
      __TransitionsX = _mm_slli_epi32(__TransitionsX, 2);
      __TransitionsD = _mm_slli_epi32(__TransitionsD, 2);
      __TransitionsX = _mm_or_si128(__TransitionsX, __ExtraMask);
      __TransitionsD = _mm_or_si128(__TransitionsD, __DeletionMask);
#  endif

      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max1, __max2);

      // Store all scores to matrix
      operation(__max1, MatrixPtr, LSEQ-1, iprf, prfLength+1); //_mm_store_si128(&pmatrix[iprf-1], __max1);

#  ifdef TAG
      // Clean extra bits
      __m128i __sign = _mm_cmplt_epi32(__max1, _mm_setzero_si128());
      __max1 = _mm_srai_epi32(__max1, 2);
      __max1 = _mm_or_si128(__max1,_mm_and_si128(__sign, __ClearMask));
#  endif

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);
#endif
    }
  }
}

#undef MAX
