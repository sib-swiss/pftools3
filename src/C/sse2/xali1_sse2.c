/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 30, 2011 xali1_sse2.c
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "../include/pfProfile.h"
#include "sse2_inline_fcts.h"
#define MAX(a,b) (a>b) ? a : b

union sIOP {
    struct { float Match; float Insertion; } Element;
    __m64 mm;
};

int xali1_sse2(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                int * const WORK, const size_t BSEQ, const size_t LSEQ, const int CutOff, const _Bool LOPT)
/*
 * WARNING: for SSE version, WORK should be aligned on 64b and 4 times the (profile size + 1)*sizeof(int)
 *          + 63 to align to cache line
 */
{
  float lScore = (float) NLOW;
  union { int i; float f;} Transfer;
  union lScores KOPD;
  __m128 __lScore;

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
    register const StoredIntegerFormat * restrict lMatch = (const StoredIntegerFormat *) &Match[_D];
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    __m128 __temp       = LoadAndConvertToFloatStoredIntegerVector(&(FirstSequenceProtein[0]));
    KOPD.xmmf           = _mm_move_ss(_mm_setzero_ps(), __temp);
    StoreMatchInsertion(&(IOP_W[0].mm), __temp);
    FirstSequenceProtein++;
    register const TransitionScores (* restrict pTransitions) = &Transitions[1];
    register union sIOP * restrict pIOP = &IOP_W[1];
    register int Length = - (int) prf->Length;

    do {
      register __m128 __KD = _mm_add_ss(KOPD.xmmf, _my_cvtdw_ss(*lMatch));
      lMatch += AlignStep;

      // dispatch value to all tuple
      __KD = _mm_shuffle_ps(__KD, __KD, _MM_SHUFFLE(0,0,0,0));

      // Load Transitions / Convert signed WORD into float / Add KD to Transitions
      __m128 __Transitions = LoadAndConvertToFloatAndAddStoredIntegerVector(&(pTransitions->From[DELETION]), __KD);

      // Move to next profile transitions
      pTransitions++;

      // Load FirstSequenceProtein / Convert signed WORD into float
      __m128 __FirstSequenceProtein = LoadAndConvertToFloatStoredIntegerVector(&(FirstSequenceProtein[0]));

      // Move to next profile First Sequence
      FirstSequenceProtein++;

      // Get maximum
      __m128 __max = _mm_max_ps(__Transitions, __FirstSequenceProtein);

      // Store IOPI and IOPM
      StoreMatchInsertion( &(pIOP->mm), (__m128) __max);
      pIOP++;

      // Set KOPD
      KOPD.xmmf = _mm_move_ss(_mm_setzero_ps(), __max);

      Length++;
    } while (Length < 0);
  }

  // Swap and assign Read and write pointers
  IOP_R = IOP_W;
  IOP_W = (union sIOP*) (((uintptr_t) &WORK[2*(prf->Length+1)] + 63) & ~63);

#ifdef XALI1_DEBUG
  fprintf(stdout,"XALI1 SCORE %12i\n", (int) lScore);
#endif
  __asm__ (" movss %1,%0" : "=x"(__lScore) : "m"(lScore));

  for ( int iseq=BSEQ; iseq < LSEQ-1; ++iseq) {
//       printf("%i %i\t", iseq+1, (int) lScore);
    register const size_t j1 = (size_t) Sequence[iseq];
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    {
      register float KI = IOP_R[0].Element.Insertion + (float) lInsertion[j1];

      // Transform KI into a vector
      __m128 __KI = _mm_set1_ps(KI);
      // Load Transitions / Convert signed WORD into float / Add KI to Transition
      __m128 __TransitionsI = LoadAndConvertToFloatAndAddStoredIntegerVector(&(Transitions[0].From[INSERTION]), __KI);

      // Load Transitions and Convert signed WORD into float
      union lScores __TransitionsX = { xmmf: LoadAndConvertToFloatStoredIntegerVector(&(Transitions[0].From[EXTRA]))};
      // Insert lScore into __TransitionsX
      __TransitionsX.xmmf = _my_insert_ss_ps_POS1(__TransitionsX.xmmf, __lScore);

      // Get maximum
      __TransitionsX.xmmf = _mm_max_ps(__TransitionsI, __TransitionsX.xmmf);

      // Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[0].mm), (__m128) __TransitionsX.xmmf);

      // Store KOPD
      KOPD.xmmf = _mm_move_ss(KOPD.xmmf, __TransitionsX.xmmf);

      __lScore = _mm_shuffle_ps(__TransitionsX.xmmf, __TransitionsX.xmmf, _MM_SHUFFLE(1,1,1,1));
#ifdef XALI1_DEBUG
      fprintf(stdout,"XALI1 SCORE SEQ %i %12i\n", iseq, _mm_cvttss_si32(__lScore));
#endif
    }

    lInsertion += AlignStep;
    register const StoredIntegerFormat * restrict lMatch = Match;


    for (int iprf=1; iprf<=prf->Length; ++iprf ) {
      const float KM = IOP_R[iprf-1].Element.Match + (float) lMatch[j1];
      const float KI = IOP_R[iprf].Element.Insertion   + (float) lInsertion[j1];
      union lScores KD    = { xmmf: _mm_add_ss(KOPD.xmmf, _my_cvtdw_ss(lMatch[_D]))} ;

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      // Transform KM into a vector
      __m128 __KM = _mm_set1_ps(KM);
      // Load Transitions / Convert signed WORD into float / Add KM to Transition
      __m128 __TransitionsM = LoadAndConvertToFloatAndAddStoredIntegerVector(&(Transitions[iprf].From[MATCH]), __KM);

      // Transform KI into a vector
      __m128 __KI = _mm_set1_ps(KI);
      // Load Transitions / Convert signed WORD into float / Add KI to Transition
      __m128 __TransitionsI = LoadAndConvertToFloatAndAddStoredIntegerVector(&(Transitions[iprf].From[INSERTION]), __KI);

#ifdef XALI1_DEBUG
      fprintf(stdout,"XALI1 M    SEQ %4i %12i %12i %12i %12i\n", iseq, Transitions[iprf].From[MATCH].To[MATCH],
	Transitions[iprf].From[MATCH].To[INSERTION],Transitions[iprf].From[MATCH].To[DELETION],Transitions[iprf].From[MATCH].To[EXTRA]);
      fprintf(stdout,"XALI1 I    SEQ %4i %12i %12i %12i %12i\n", iseq, Transitions[iprf].From[INSERTION].To[MATCH],
	Transitions[iprf].From[INSERTION].To[INSERTION],Transitions[iprf].From[INSERTION].To[DELETION],Transitions[iprf].From[INSERTION].To[EXTRA]);
      fprintf(stdout,"XALI1 D    SEQ %4i %12i %12i %12i %12i\n", iseq, Transitions[iprf].From[DELETION].To[MATCH],
	Transitions[iprf].From[DELETION].To[INSERTION],Transitions[iprf].From[DELETION].To[DELETION],Transitions[iprf].From[DELETION].To[EXTRA]);
      fprintf(stdout,"XALI1 X    SEQ %4i %12i %12i %12i %12i\n", iseq, Transitions[iprf].From[EXTRA].To[MATCH],
	Transitions[iprf].From[EXTRA].To[INSERTION],Transitions[iprf].From[EXTRA].To[DELETION],Transitions[iprf].From[EXTRA].To[EXTRA]);
      fprintf(stdout,"XALI1 KMID SEQ %4i %12i %12i %12i %12i\n", iseq, (int) KM, (int) KI, (int) KD.Elementf[0], 0);
#endif
      // Get maximum
      const __m128 __max1 = _mm_max_ps(__TransitionsM, __TransitionsI);

       // Load Transitions and Convert signed WORD into float
      union lScores __TransitionsX = { xmmf: LoadAndConvertToFloatStoredIntegerVector(&(Transitions[iprf].From[EXTRA]))};
      // Insert lScore into __TransitionsX
      __TransitionsX.xmmf = _my_insert_ss_ps_POS1(__TransitionsX.xmmf, __lScore);

      // Transform KD into a vector
      KD.xmmf = _mm_shuffle_ps(KD.xmmf, KD.xmmf, _MM_SHUFFLE(0,0,0,0));
      // Load Transitions / Convert signed WORD into float / Add KD to Transition
      __m128 __TransitionsD = LoadAndConvertToFloatAndAddStoredIntegerVector(&(Transitions[iprf].From[DELETION]), KD.xmmf);


      // Get maximum / Set KOPD
      __TransitionsX.xmmf = _mm_max_ps(__TransitionsD, __TransitionsX.xmmf);
      KOPD.xmmf           = _mm_max_ps(__max1, __TransitionsX.xmmf);

      // Store IOPI and IOPM
      StoreMatchInsertion( &(IOP_W[iprf].mm), (__m128) KOPD.xmmf);

      __lScore = _mm_shuffle_ps(KOPD.xmmf, KOPD.xmmf, _MM_SHUFFLE(1,1,1,1));
#ifdef XALI1_DEBUG
      fprintf(stdout,"XALI1 SCORE SEQ %4i %12i %8.1f %8.1f %8.1f %8.1f\n", iseq, _mm_cvttss_si32(__lScore),
	      KOPD.Elementf[2], KOPD.Elementf[3], KOPD.Elementf[0], KOPD.Elementf[1]);
#endif
//       printf("%i %i\t\t%i\t%i\t\t%i\t%i\t%i\t\t%i\t%i\t\t%i\n",
//              iseq, iprf,
//              (int) IOP_W[iprf].Element.Match, (int) IOP_W[iprf].Element.Insertion,
//              (int)KM, (int)KI, (int)KD.Elementf[0], (int) 0, (int) KOPD.Elementf[DELETION],
//              (int)lScore);

    } //while (++iprf <= prf->Length);

    // Swap Read and Write pointers
    const union sIOP * const ptr = IOP_W;
    IOP_W = (union sIOP*) IOP_R;
    IOP_R = ptr;

    if ( ! LOPT ) {
      __m128 __dummy;
      __asm__ ("cvtsi2ss  %1, %0" : "=x"(__dummy) : "m"(CutOff));
      if (_mm_comige_ss(__lScore, __dummy)) {
	return _mm_cvttss_si32(__lScore);
      }
    }
  }
  int iScore = _mm_cvttss_si32(__lScore);
  {
    register const StoredIntegerFormat * restrict lInsertion = Insertion;
    const int j1 = (int) Sequence[LSEQ-1];
    int iKOPM    = (int) IOP_R[0].Element.Match;
    int KI       = (int) IOP_R[0].Element.Insertion + (int) lInsertion[j1];

    int iKOPD    = MAX( KI + (int) Transitions[0].Element[_ID],      (int) Transitions[0].Element[_XD] );
    register const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;
    iScore = MAX( iScore, KI + (int) LastSequenceProtein[0].From[INSERTION] );

    register const StoredIntegerFormat * restrict lMatch = Match;
    lInsertion += AlignStep;
    register const int itmp = prf->Length;
    for (int iprf=1; iprf<=itmp; ++iprf) {
      const int KM = iKOPM                               + (int) lMatch[j1];
      KI           = (int) IOP_R[iprf].Element.Insertion + (int) lInsertion[j1];
      const int KD = iKOPD                               + (int) lMatch[_D];

      lMatch     += AlignStep;
      lInsertion += AlignStep;

      iKOPM = (int) IOP_R[iprf].Element.Match;

      const float tIOPD1 = MAX( KM + (int) Transitions[iprf].Element[_MD],      (int) Transitions[iprf].Element[_XD] );
      const float tIOPD2 = MAX( KI + (int) Transitions[iprf].Element[_ID], KD + (int) Transitions[iprf].Element[_DD] );
      iKOPD              = MAX( tIOPD1, tIOPD2);

      const int tIOPT1 = MAX( KM + (int) LastSequenceProtein[iprf].From[MATCH], KI + (int) LastSequenceProtein[iprf].From[INSERTION] );
      const int tIOPT2 = MAX( iScore                                          , KD + (int) LastSequenceProtein[iprf].From[DELETION] );
      iScore           = MAX( tIOPT1, tIOPT2);
#ifdef XALI1_DEBUG
      fprintf(stdout,"XALI1 SCORE %12i\n", iScore);
#endif
    }
  }
  return iScore;
}


#undef MAX
