/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 30, 2011 xali1_sse41.c
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <smmintrin.h>
#include "../include/pfProfile.h"
#include "sse41_inline_fcts.h"

#define MAX(a,b) (a>b) ? a : b

union sIOP {
    struct { int Match; int Insertion; } Element;
    __m64 mm;
};
 
int xali1_sse41(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                int * const WORK, const size_t BSEQ, const size_t LSEQ, const int CutOff, const _Bool LOPT)
/*
 * WARNING: for SSE version, WORK should be aligned on 64b and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 */
{
  int KOPD, lScore = (int) NLOW;
  
  const union sIOP * restrict IOP_R;
  union sIOP * restrict IOP_W = (union sIOP*) WORK;

  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * const restrict Match = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * const restrict Insertion = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep = prf->Scores.Match.AlignStep;

  /* NOTE: The following part could be replaced and performed only once for a profile as it
   *       is profile dependent. Nevertheless it does a good job loading Match and Transition
   *       matrices into the cache hierarchy.
   */
  {
    /* 
     * Initialize Profile Entrance Line 
     */
    register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    IOP_W[0].Element.Match     = (int) FirstSequenceProtein[0].To[MATCH];
    IOP_W[0].Element.Insertion = (int) FirstSequenceProtein[0].To[INSERTION];
    KOPD                       = (int) FirstSequenceProtein[0].To[DELETION];
    FirstSequenceProtein++;
    register const TransitionScores (* restrict pTransitions) = &Transitions[1];
    register union sIOP * restrict pIOP = &IOP_W[1];
    register int Length = - (int) prf->Length;

    do {
      register const int KD = KOPD + (int) *lMatch;
      lMatch += AlignStep;
      
      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __Transitions = LoadStoredIntegerVector(&(pTransitions->From[DELETION]));
            
      // Add KD to Transitions
      __Transitions = _mm_add_epi32(__Transitions, __KD);
      
      // Move to next profile transitions
      pTransitions++;

      // Load FirstSequenceProtein
      __m128i __FirstSequenceProtein = LoadStoredIntegerVector(&FirstSequenceProtein[0]);


      // Move to next profile First Sequence
      FirstSequenceProtein++;
      
      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__Transitions, __FirstSequenceProtein);

      // Store IOPI and IOPM
      StoreMatchInsertion( &(pIOP->mm), (__m128) __max);
      pIOP++;
      
      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max, DELETION);

      Length++;
    } while (Length < 0);
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = (union sIOP*) (((uintptr_t) &WORK[2*(prf->Length+1)] + 63) & ~63);

  const size_t prfLength = prf->Length;
  for ( int iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
    register const size_t j1 = (size_t) Sequence[iseq];
    int KOPM = IOP_R[0].Element.Match;
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
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
     
      // Insert lScore into __TransitionsX
      __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max = _mm_max_epi32(__TransitionsI, __TransitionsX);

      // Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __max);
      
      // Store KOPD
      KOPD = _mm_extract_epi32(__max, DELETION);

      // Backup new score to xmm register
      lScore = _mm_extract_epi32(__max, DUMMY);
    }
    
    lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;
    
    for (size_t iprf=1; iprf<=prfLength; ++iprf ) {
      const int KM = KOPM                          + (int) lMatch[j1];
      const int KI = IOP_R[iprf].Element.Insertion + (int) lInsertion[j1];
      const int KD = KOPD                          + (int) lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;

      // Transform KM into a vector
      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));

      // Add KM to Transition
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);

    
      // Transform KI into a vector
      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Add KI to Transition
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      // Get maximum ( this is SSE 4.1 )
      __m128i __max1 = _mm_max_epi32(__TransitionsM, __TransitionsI);

      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));
      // Insert lscore into TransitionX
      __TransitionsX = _mm_insert_epi32(__TransitionsX, lScore, DUMMY);
      
      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Add KD to Transition
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);
      
      // Get maximum ( this is SSE 4.1 )
      __m128i __max2 = _mm_max_epi32(__TransitionsD, __TransitionsX);
      __max1 = _mm_max_epi32(__max1, __max2);

      // Store IOPI and IOPM
      StoreMatchInsertion( &IOP_W[iprf].mm, (__m128) __max1);

      // Set KOPD ( this is SSE 4.1 )
      KOPD = _mm_extract_epi32(__max1, DELETION);

      lScore = _mm_extract_epi32(__max1, DUMMY);
#ifdef XALI1_DEBUG
       printf("%i %i\t\t%i\t%i\t\t%i\t%i\t%i\t\t%i\t%i\t\t%i\n",
              iseq, iprf,
              IOP_W[iprf].Element.Match, IOP_W[iprf].Element.Insertion,
              KM, KI, KD, KOPM, KOPD,
              lScore);
#endif
    } //while (++iprf <= prf->Length);

    // Swap Read and Write pointers
    const union sIOP * const ptr = IOP_W;
    IOP_W = (union sIOP*) IOP_R;
    IOP_R = ptr;

    if ( ! LOPT && lScore >= CutOff) return lScore;
  } 
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    int KOPM     = IOP_R[0].Element.Match;
    int KI       = IOP_R[0].Element.Insertion + (int) lInsertion[j1];
    
    KOPD   = MAX( KI + (int) Transitions[0].Element[_ID],      (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;
    lScore = MAX( lScore, KI + (int) LastSequenceProtein[0].From[INSERTION] );
  
    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;
    
    for (size_t iprf=1; iprf<=prfLength; ++iprf) {
      const int KM = KOPM                          + lMatch[j1];
      KI           = IOP_R[iprf].Element.Insertion + lInsertion[j1];
      const int KD = KOPD                          + lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      KOPM = IOP_R[iprf].Element.Match;

      const int tIOPD1 = MAX( KM + (int) Transitions[iprf].Element[_MD],      (int) Transitions[iprf].Element[_XD] );
      const int tIOPD2 = MAX( KI + (int) Transitions[iprf].Element[_ID], KD + (int) Transitions[iprf].Element[_DD] );
      KOPD             = MAX( tIOPD1, tIOPD2);

      const int tIOPT1 = MAX( KM + (int) LastSequenceProtein[iprf].From[MATCH], KI + (int) LastSequenceProtein[iprf].From[INSERTION] );
      const int tIOPT2 = MAX( lScore                                          , KD + (int) LastSequenceProtein[iprf].From[DELETION] );
      lScore           = MAX( tIOPT1, tIOPT2);
    }
  }
  return lScore;
}

#undef MAX
