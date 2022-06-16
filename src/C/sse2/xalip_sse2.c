/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 30, 2011 xalip_sse2.c
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
#define MIN(a,b) (a<b) ? a : b

struct Minima {
   int _a;
   int _b;
   int _c;
};


static void InitR(const int iseq, const size_t N1, const size_t N2, const size_t bseq, const size_t lseq,
                  union lScores * const restrict iop, union Positions * const restrict iom,
                  union Positions * const restrict ioi,  const struct Profile * const restrict prf)
{
  int KOPD;
#ifdef XALIP_DEBUG
  fprintf(stdout,"XALIP InitR called with index %i bseq %lu\n",iseq, bseq);
#endif
  // Are we treating sequence index below given start?
  register const ScoreTuple * restrict FirstSequenceProtein = ( iseq < bseq )
      ? &(prf->Scores.Insertion.FirstSequenceProtein[0])
      : &(prf->Scores.Insertion.Transitions->From[EXTRA]);
  const size_t FirstSequenceProteinAlignStep = (iseq < bseq) ? 1 : 4;

  register int (* restrict pIOP)[4] = &(iop[0].Element);

  for (int i=0; i<4; ++i) {
    pIOP[0][i] = (int) FirstSequenceProtein->To[i];
  }

  // Set KOPD
  KOPD = pIOP[0][DELETION];

  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  register const StoredIntegerFormat * restrict pMatch       = &prf->Scores.Match.Alphabet[_D];
  const size_t AlignStep = prf->Scores.Match.AlignStep;

  // Move to next profile First Sequence
  FirstSequenceProtein += FirstSequenceProteinAlignStep;

  for (size_t iprf=1; iprf<=prf->Length; ++iprf) {
    register const int KD = KOPD + (int) *pMatch;
    pMatch += AlignStep;

    // Load Transitions
    int __Transitions[4];
    __Transitions[MATCH]     = KD + (int) Transitions[iprf].From[DELETION].To[MATCH];
    __Transitions[INSERTION] = KD + (int) Transitions[iprf].From[DELETION].To[INSERTION];
    __Transitions[DELETION]  = KD + (int) Transitions[iprf].From[DELETION].To[DELETION];
    __Transitions[EXTRA]     = KD + (int) Transitions[iprf].From[DELETION].To[EXTRA];

    // Move to next profile First Sequence
    FirstSequenceProtein += FirstSequenceProteinAlignStep;

    // Get maximum
    int __FirstSequenceProtein[4];
    for (int i=0; i<4; ++i) {
      __FirstSequenceProtein[i] = (int) FirstSequenceProtein->To[i];
      pIOP[iprf][i] = (__Transitions[i] > __FirstSequenceProtein[i]) ? __Transitions[i] : __FirstSequenceProtein[i];
    }

    // Set KOPD ( this is SSE 4.1 )
    KOPD = pIOP[iprf][DELETION];
  }

  if (N1 > 0) {
    union Positions TPOS __attribute__((aligned(16)));
    TPOS.Region.Protected.Begin = lseq + 1;
    TPOS.Region.Protected.End   = 0;
    TPOS.Region.Sequence.Begin  = iseq + 1;
    TPOS.Region.Sequence.End    = 0; // this is dummy
    for (size_t iprf=0; iprf<(N1-1); ++iprf) {
      _mm_store_si128(&(iom[iprf].xmm), TPOS.xmm);
      _mm_store_si128(&(ioi[iprf].xmm), TPOS.xmm);
    }
    ioi[N1-1].xmm = TPOS.xmm;

    TPOS.Region.Protected.Begin = TPOS.Region.Sequence.Begin;
    TPOS.Region.Protected.End = TPOS.Region.Sequence.Begin;
    iom[N1-1].xmm = TPOS.xmm;
    for (size_t iprf=N1; iprf<N2; ++iprf) {
      _mm_store_si128(&(iom[iprf].xmm), TPOS.xmm);
      _mm_store_si128(&(ioi[iprf].xmm), TPOS.xmm);
    }

    TPOS.Region.Protected.Begin = lseq + 1;
    TPOS.Region.Protected.End = 0;
    for (size_t iprf=N2; iprf<=prf->Length; ++iprf) {
      _mm_store_si128(&(iom[iprf].xmm), TPOS.xmm);
      _mm_store_si128(&(ioi[iprf].xmm), TPOS.xmm);
    }
  } else {
    fputs("BUG HERE N1 is NULL\n", stderr);
    exit(1);
  }
}

static void nextR(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                  const int iseq, union lScores * const restrict iop, union Positions * const restrict iom,
                  union Positions * const restrict ioi,const int lseq, struct Alignment * const restrict alignment,
                  struct Minima * const restrict ifer, const _Bool Lock, const size_t N1, const size_t N2)
{
#ifdef XALIP_DEBUG
   fprintf(stdout,"XALIP NextR called with iseq %i\n",iseq);
#endif
   // Initialization
   const unsigned int In = iseq + 1;
   // WARNING: Fortran uses a 1 based index for sequence
   const unsigned int SequenceIndex = (unsigned int) Sequence[iseq-1];
#ifdef XALIP_DEBUG
   if ( iseq >= lseq) {
      fputs("nextR_last should have been called\n", stderr);
      exit(1);
   }
#endif
   // Disable match and insert vertices of protected region
   if (Lock) {
      iop[N1-1].Element[MATCH] = NLOW;
      for (size_t iprf=N1; iprf<N2; ++iprf) {
         iop[iprf].Element[MATCH]     = NLOW;
         iop[iprf].Element[INSERTION] = NLOW;
      }
   }
   ////////////////////////////////////////////////////////////////////////////////////////////
   // Profile position 0
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Save previous match position
   int Kopm = iop[0].Element[MATCH];
   union Positions Kpos __attribute__((aligned(16)));
   Kpos.xmm = iom[0].xmm;

   const union Positions TEMPpos __attribute__((aligned(16))) = { lseq+1, 0, In, 0 };

   const union Positions * restrict PTRpos[4];
   union Positions Kiod;
   PTRpos[0] = &TEMPpos;
   PTRpos[1] = &Kpos;
   PTRpos[3] = &Kiod;

   // Get pointers to score data
   const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
   const StoredIntegerFormat * restrict Insertion = prf->Scores.Insertion.Alphabet;
   const size_t AlignStep = prf->Scores.Insertion.AlignStep;

   int Ki = iop[0].Element[INSERTION] + (int) Insertion[SequenceIndex];

   // Match position
   int Ji   = Ki + (int) Transitions[0].From[INSERTION].To[MATCH];
   int itmp = (int) Transitions[0].From[EXTRA].To[MATCH];
   if ( Ji > itmp) {
      iop[0].Element[MATCH] = Ji;
      iom[0].xmm = ioi[0].xmm;
   } else {
      iop[0].Element[MATCH] = itmp;
      iom[0].xmm = TEMPpos.xmm;
   }

   // Deletion position
   int Kopd;
   Ji   = Ki + (int) Transitions[0].From[INSERTION].To[DELETION];
   itmp = (int) Transitions[0].From[EXTRA].To[DELETION];
   if ( Ji > itmp ) {
      Kopd     = Ji;
      Kiod.xmm = ioi[0].xmm;
   } else {
      Kopd     = itmp;
      Kiod.xmm = TEMPpos.xmm;
   }

   // Insertion position
   Ji   = Ki + (int) Transitions[0].From[INSERTION].To[INSERTION];
   itmp = (int) Transitions[0].From[EXTRA].To[INSERTION];
   if ( Ji > itmp ) {
      iop[0].Element[INSERTION] = Ji;
   } else {
      iop[0].Element[INSERTION] = itmp;
      ioi[0].xmm = TEMPpos.xmm;
   }

   // Initialize minima
   ifer->_a = iseq;
   ifer->_b = iseq;
   itmp     = MIN(ioi[0].Region.Sequence.Begin, iom[0].Region.Sequence.Begin);
   ifer->_c = MIN(itmp, Kiod.Region.Sequence.Begin);

   // Initialize alignment
   union Positions Fpos __attribute__((aligned(16)));;
   Fpos.Region.Protected.Begin = alignment->Region.Protected.Begin;
   Fpos.Region.Protected.End   = alignment->Region.Protected.End;
   Fpos.Region.Sequence.Begin  = alignment->Region.Sequence.Begin;
   Fpos.Region.Sequence.End    = 0; // this is dummy

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Loop through rest of profile
   ////////////////////////////////////////////////////////////////////////////////////////////
   const StoredIntegerFormat * restrict Match = prf->Scores.Match.Alphabet;
   Insertion += AlignStep;

   for (size_t iprf=1; iprf<=prf->Length; ++iprf) {
      /////////////////////////////////////////////////////////////////////////////////////////
      // Match
      const register int KM = Kopm + (int) Match[SequenceIndex];
      Kopm = iop[iprf].Element[MATCH];

      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Add KM to Transitions
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      /////////////////////////////////////////////////////////////////////////////////////////
      // Insertion
      const register int KI = iop[iprf].Element[INSERTION] + (int) Insertion[SequenceIndex];
      // one could move on the seq index instead

      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Add KM to Transitions
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Deletion
      const register int KD = Kopd + (int) Match[_D];

      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
     // Add KM to Transitions
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Extensions
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));

      // Insert NLOW into Extension vector -> done in io.c
      //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);

      // Move to next profile index
      Match     += AlignStep;
      Insertion += AlignStep;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Compute Maxima (Fortran Nstep function)
      union lScores Index = { xmm: (__m128i) _mm_setzero_si128()};
      __m128i __three = _mm_set_epi32(3,3,3,3);

      __m128i  __Mask = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
      __TransitionsX  = _my_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
      Index.xmm       = _my_blendv_epi8(Index.xmm, __three, __Mask);

      __m128i __One  = _mm_set_epi32(1,1,1,1);
      __Mask         = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
      __TransitionsX = _my_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
      Index.xmm      = _my_blendv_epi8(Index.xmm, __One, __Mask);

      __m128i __Two  = _mm_add_epi32(__One, __One);
      __Mask         = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
      __TransitionsX = _my_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
      Index.xmm      = _my_blendv_epi8(Index.xmm, __Two, __Mask);

      // Set new data
      iop[iprf].xmm = __TransitionsX;
      //StoreMatchInsertion((__m64*) &iop[iprf],(__m128) __TransitionsX);
      Kopd = _mm_cvtsi128_si32 ( __TransitionsX);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Check for new maxima
      __TransitionsX = _mm_srli_si128(__TransitionsX, 4);
      const int KE = _mm_cvtsi128_si32 ( __TransitionsX);

      if (KE > alignment->Score) {
         alignment->Score = KE;
         alignment->Region.Sequence.End = iseq;
         Fpos.xmm = ioi[iprf].xmm;

         if (Index.Element[DUMMY] == 1)  { // KM is max
               Fpos.xmm = Kpos.xmm;
         } else if (Index.Element[DUMMY] == 3) { // KD is max
               Fpos.xmm = Kiod.xmm;
         }
      }
#ifdef XALIP_DEBUG
      printf("XALIP FPOS %8i %8i %8i\n", Fpos.Region.Protected.Begin, Fpos.Region.Protected.End, Fpos.Region.Sequence.Begin);
#endif
      /////////////////////////////////////////////////////////////////////////////////////////
      // Update alignment positions
      union Positions Jpos __attribute__((aligned(16)));
      Jpos.xmm  = iom[iprf].xmm;
      PTRpos[2] = &ioi[iprf];

      iom[iprf].xmm  = PTRpos[Index.Element[MATCH]]->xmm;
      Kiod.xmm       = PTRpos[Index.Element[DELETION]]->xmm;
      ioi[iprf].xmm  = PTRpos[Index.Element[INSERTION]]->xmm;

      Kpos.xmm = Jpos.xmm;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Update minima

      const int t1 = MIN(ioi[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.Begin);
      const int t2 = MIN(t1, Kiod.Region.Protected.Begin);
      ifer->_a     = MIN(ifer->_a, t2);

      if (iprf > N1) {
         ifer->_b = MIN(ifer->_b, t2);
      }
      const int t3 = MIN(ioi[iprf].Region.Sequence.Begin, iom[iprf].Region.Sequence.Begin);
      const int t4 = MIN(t3, Kiod.Region.Sequence.Begin);
      ifer->_c     = MIN(ifer->_c, t4);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Epilogue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Finish updating alignment positions
   alignment->Region.Protected.Begin = Fpos.Region.Protected.Begin;
   alignment->Region.Protected.End = Fpos.Region.Protected.End;
   alignment->Region.Sequence.Begin = Fpos.Region.Sequence.Begin;

   // Entry and exit from protected regions
   iom[N1-1].Region.Protected.Begin = MIN(iom[N1-1].Region.Protected.Begin, In);
   iom[N1-1].Region.Protected.End = In;

   for (size_t iprf=N1; iprf<N2; ++iprf) {
      iom[iprf].Region.Protected.Begin = MIN(iom[iprf].Region.Protected.Begin, In);
      iom[iprf].Region.Protected.End = In;

      ioi[iprf].Region.Protected.Begin = MIN(ioi[iprf].Region.Protected.Begin, In);
      ioi[iprf].Region.Protected.End = In;
   }

#ifdef XALIP_DEBUG
   for (size_t iprf=0; iprf<=prf->Length; ++iprf) {
      printf("XALIP NEXTR IOP %4.4lu %15i %15i %15i\n", iprf, iop[iprf].Element[MATCH], iop[iprf].Element[INSERTION],
             iop[iprf].Element[DELETION]);
      printf("XALIP NEXTR IOM %4.4lu %15i %15i %15i\n", iprf, iom[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.End,
             iom[iprf].Region.Sequence.Begin);
      printf("XALIP NEXTR IOI %4.4lu %15i %15i %15i\n", iprf, ioi[iprf].Region.Protected.Begin, ioi[iprf].Region.Protected.End,
             ioi[iprf].Region.Sequence.Begin);

   }

   printf("XALIP NEXTR ALIGN %4i %4i %4i %4i %4i\n",
                   alignment->Region.Protected.Begin, alignment->Region.Protected.End, alignment->Score, alignment->Region.Sequence.Begin,
                   alignment->Region.Sequence.End);
#endif
}

static void nextR_last(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                       const int iseq, union lScores * const restrict iop, union Positions * const restrict iom,
                       union Positions * const restrict ioi,const int lseq, struct Alignment * const restrict alignment,
                       struct Minima * const restrict ifer, const _Bool Lock, const size_t N1, const size_t N2)
{
#ifdef XALIP_DEBUG
   fprintf(stdout,"XALIP NextR_last called with iseq %i\n",iseq);
#endif
   // Initialization
   const unsigned int In = iseq + 1;
   // WARNING: Fortran uses a 1 based index for sequence
   const unsigned int SequenceIndex = (unsigned int) Sequence[iseq-1];
#ifdef XALIP_DEBUG
   if ( iseq < lseq) {
      fputs("nextR should have been called\n", stderr);
      exit(1);
   }
#endif
   // Disable match and insert vertices of protected region
   if (Lock) {
      iop[N1-1].Element[MATCH] = NLOW;
      for (size_t iprf=N1; iprf<N2; ++iprf) {
         iop[iprf].Element[MATCH]     = NLOW;
         iop[iprf].Element[INSERTION] = NLOW;
      }
   }
   ////////////////////////////////////////////////////////////////////////////////////////////
   // Profile position 0
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Save previous match position
   int Kopm = iop[0].Element[MATCH];
   union Positions Kpos __attribute__((aligned(16)));
   Kpos.xmm = iom[0].xmm;

   const union Positions TEMPpos __attribute__((aligned(16))) = { lseq+1, 0, In, 0 };

   const union Positions * restrict PTRpos[4];
   union Positions Kiod;
   PTRpos[0] = &TEMPpos;
   PTRpos[1] = &Kpos;
   PTRpos[3] = &Kiod;

   // Get pointers to score data
   const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
   const StoredIntegerFormat * restrict Insertion = prf->Scores.Insertion.Alphabet;
   const size_t AlignStep = prf->Scores.Insertion.AlignStep;

   int Ki = iop[0].Element[INSERTION] + (int) Insertion[SequenceIndex];

   // Match position
   int Ji   = Ki + (int) Transitions[0].From[INSERTION].To[MATCH];
   int itmp = (int) Transitions[0].From[EXTRA].To[MATCH];
   if ( Ji > itmp) {
      iop[0].Element[MATCH] = Ji;
      iom[0].xmm = ioi[0].xmm;
   } else {
      iop[0].Element[MATCH] = itmp;
      iom[0].xmm = TEMPpos.xmm;
   }

   // Deletion position
   int Kopd;
   Ji   = Ki + (int) Transitions[0].From[INSERTION].To[DELETION];
   itmp = (int) Transitions[0].From[EXTRA].To[DELETION];
   if ( Ji > itmp ) {
      Kopd     = Ji;
      Kiod.xmm = ioi[0].xmm;
   } else {
      Kopd     = itmp;
      Kiod.xmm = TEMPpos.xmm;
   }

   // Insertion position
   Ji   = Ki + (int) Transitions[0].From[INSERTION].To[INSERTION];
   itmp = (int) Transitions[0].From[EXTRA].To[INSERTION];
   if ( Ji > itmp ) {
      iop[0].Element[INSERTION] = Ji;
   } else {
      iop[0].Element[INSERTION] = itmp;
      ioi[0].xmm = TEMPpos.xmm;
   }

   // Initialize minima
   ifer->_a = iseq;
   ifer->_b = iseq;
   itmp     = MIN(ioi[0].Region.Sequence.Begin, iom[0].Region.Sequence.Begin);
   ifer->_c = MIN(itmp, Kiod.Region.Sequence.Begin);

   // Initialize alignment
   union Positions Fpos __attribute__((aligned(16)));
   Fpos.Region.Protected.Begin = alignment->Region.Protected.Begin;
   Fpos.Region.Protected.End   = alignment->Region.Protected.End;
   Fpos.Region.Sequence.Begin  = alignment->Region.Sequence.Begin;
   Fpos.Region.Sequence.End    = 0; // this is dummy

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Loop through rest of profile
   ////////////////////////////////////////////////////////////////////////////////////////////
   const StoredIntegerFormat * restrict Match = prf->Scores.Match.Alphabet;
   const ScoreTuple * const restrict LastProteinSequence = prf->Scores.Insertion.LastSequenceProtein;
   Insertion += AlignStep;

   for (size_t iprf=1; iprf<=prf->Length; ++iprf) {
      /////////////////////////////////////////////////////////////////////////////////////////
      // Match
      const register int KM = Kopm + (int) Match[SequenceIndex];
      Kopm = iop[iprf].Element[MATCH];

      __m128i __KM = _mm_set1_epi32(KM);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsM = LoadStoredIntegerVector(&(Transitions[iprf].From[MATCH]));
      // Insert LastProteinSequence
      __TransitionsM = _my_insert_epi32_POS1(__TransitionsM, (int) LastProteinSequence[iprf].From[MATCH]);
      // Add KM to Transitions
      __TransitionsM = _mm_add_epi32(__TransitionsM, __KM);


      /////////////////////////////////////////////////////////////////////////////////////////
      // Insertion
      const register int KI = iop[iprf].Element[INSERTION] + (int) Insertion[SequenceIndex];
      // one could move on the seq index instead

      __m128i __KI = _mm_set1_epi32(KI);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsI = LoadStoredIntegerVector(&(Transitions[iprf].From[INSERTION]));
      // Insert LastProteinSequence
      __TransitionsI = _my_insert_epi32_POS1(__TransitionsI, (int) LastProteinSequence[iprf].From[INSERTION]);
      // Add KM to Transitions
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Deletion
      const register int KD = Kopd + (int) Match[_D];

      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _my_insert_epi32_POS1(__TransitionsD, (int) LastProteinSequence[iprf].From[DELETION]);
      // Add KM to Transitions
      __TransitionsD = _mm_add_epi32(__TransitionsD, __KD);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Extensions
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsX = LoadStoredIntegerVector(&(Transitions[iprf].From[EXTRA]));

      // Insert NLOW into Extension vector -> done in io.c
      //__TransitionsX = _mm_insert_epi32(__TransitionsX, NLOW, DUMMY);

      // Move to next profile index
      Match     += AlignStep;
      Insertion += AlignStep;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Compute Maxima (Fortran Nstep function)
      union lScores Index = { xmm: (__m128i) _mm_setzero_si128()} ;
      __m128i __three = _mm_set_epi32(3,3,3,3);

      __m128i  __Mask = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
      __TransitionsX  = _my_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
      Index.xmm       = _my_blendv_epi8(Index.xmm, __three, __Mask);

      __m128i __One  = _mm_set_epi32(1,1,1,1);
      __Mask         = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
      __TransitionsX = _my_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
      Index.xmm      = _my_blendv_epi8(Index.xmm, __One, __Mask);

      __m128i __Two  = _mm_add_epi32(__One, __One);
      __Mask         = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
      __TransitionsX = _my_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
      Index.xmm      = _my_blendv_epi8(Index.xmm, __Two, __Mask);

      // Set new data
      iop[iprf].xmm = __TransitionsX;
      //_mm_storel_pi((__m64*) &iop[iprf],(__m128) __TransitionsX);
      Kopd = _mm_cvtsi128_si32(__TransitionsX);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Check for new maxima
       __TransitionsX = _mm_srli_si128(__TransitionsX, 4);
      const int KE = _mm_cvtsi128_si32 ( __TransitionsX);
      if (KE > alignment->Score) {
         alignment->Score = KE;
         alignment->Region.Sequence.End = iseq;
         Fpos.xmm = ioi[iprf].xmm;
         if (Index.Element[DUMMY] == 1)  { // KM is max
               Fpos.xmm = Kpos.xmm;
         } else if (Index.Element[DUMMY] == 3) { // KD is max
               Fpos.xmm = Kiod.xmm;
         }
      }
#ifdef XALIP_DEBUG
      printf("XALIP LAST FPOS %8i %8i %8i\n", Fpos.Region.Protected.Begin, Fpos.Region.Protected.End, Fpos.Region.Sequence.Begin);
#endif
      /////////////////////////////////////////////////////////////////////////////////////////
      // Update alignment positions
      union Positions Jpos __attribute__((aligned(16)));
      Jpos.xmm  = iom[iprf].xmm;
      PTRpos[2] = &ioi[iprf];

      iom[iprf].xmm  = PTRpos[Index.Element[MATCH]]->xmm;
      Kiod.xmm       = PTRpos[Index.Element[DELETION]]->xmm;
      ioi[iprf].xmm  = PTRpos[Index.Element[INSERTION]]->xmm;

      Kpos.xmm = Jpos.xmm;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Update minima

      const int t1 = MIN(ioi[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.Begin);
      const int t2 = MIN(t1, Kiod.Region.Protected.Begin);
      ifer->_a     = MIN(ifer->_a, t2);

      if (iprf > N1) {
         ifer->_b = MIN(ifer->_b, t2);
      }
      const int t3 = MIN(ioi[iprf].Region.Sequence.Begin, iom[iprf].Region.Sequence.Begin);
      const int t4 = MIN(t3, Kiod.Region.Sequence.Begin);
      ifer->_c     = MIN(ifer->_c, t4);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Epilogue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Finish updating alignment positions
   alignment->Region.Protected.Begin = Fpos.Region.Protected.Begin;
   alignment->Region.Protected.End = Fpos.Region.Protected.End;
   alignment->Region.Sequence.Begin = Fpos.Region.Sequence.Begin;

   // Entry and exit from protected regions
   iom[N1-1].Region.Protected.Begin = MIN(iom[N1-1].Region.Protected.Begin, In);
   iom[N1-1].Region.Protected.End = In;

   for (size_t iprf=N1; iprf<N2; ++iprf) {
      iom[iprf].Region.Protected.Begin = MIN(iom[iprf].Region.Protected.Begin, In);
      iom[iprf].Region.Protected.End = In;

      ioi[iprf].Region.Protected.Begin = MIN(ioi[iprf].Region.Protected.Begin, In);
      ioi[iprf].Region.Protected.End = In;
   }

#ifdef XALIP_DEBUG
   for (size_t iprf=0; iprf<=prf->Length; ++iprf) {
      printf("XALIP NEXTR LAST IOP %4.4lu %15i %15i %15i\n", iprf, iop[iprf].Element[MATCH], iop[iprf].Element[INSERTION],
             iop[iprf].Element[DELETION]);
      printf("XALIP NEXTR LAST IOM %4.4lu %15i %15i %15i\n", iprf, iom[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.End,
             iom[iprf].Region.Sequence.Begin);
      printf("XALIP NEXTR LAST IOI %4.4lu %15i %15i %15i\n", iprf, ioi[iprf].Region.Protected.Begin, ioi[iprf].Region.Protected.End,
             ioi[iprf].Region.Sequence.Begin);

   }

   printf("XALIP NEXTR LAST ALIGN %4i %4i %4i %4i %4i\n",
                   alignment->Region.Protected.Begin, alignment->Region.Protected.End, alignment->Score, alignment->Region.Sequence.Begin,
                   alignment->Region.Sequence.End);
#endif
}

int xalip_sse2( const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
           union lScores * const restrict iop, union Positions * const restrict iom,
           union Positions * const restrict ioi, const int bseq, const int lseq,
           struct Alignment * const restrict alignment,
           _Bool * const restrict Lock, const size_t N1, const size_t N2, const _Bool Lopt,
           const int kcut, const size_t MaxNumberOfAlignment)
{
   int iseq;
   ////////////////////////////////////////////////////////////////////////////////////////////
   // Prelogue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Alignment list
   size_t AlignmentNumber = 0;

   // Search control fields
   int ibeg    = bseq-1;
   size_t jlcp = prf->Length / 2;

   // Stack Memory
   struct Minima ifer __attribute__((aligned(16)));

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Two step forward one step backward loop
   ////////////////////////////////////////////////////////////////////////////////////////////
   MajorLoop:

   iseq = ibeg;
   struct Alignment lAlignment __attribute__((aligned(16))); ;
//    lAlignment.Region.Protected.Begin = 0;
//    lAlignment.Region.Protected.End = 0;
//    lAlignment.Region.Sequence.Begin = 0;
//    lAlignment.SequenceEnd   = 0;
   _mm_store_si128( (__m128i*) &lAlignment.Region.Protected.Begin, _mm_setzero_si128());
   lAlignment.Score = NLOW;

   // Initiate work array
   InitR(iseq, N1, N2, bseq, lseq, iop, iom, ioi, prf);
#ifdef XALIP_DEBUG
   for (size_t i=0; i<=prf->Length; ++i) {
      printf("XALIP IOP %8i %8i %8i %8i\n",
             iop[i].Element[MATCH], iop[i].Element[INSERTION], iop[i].Element[DELETION],
             iop[i].Element[DUMMY]);
      printf("XALIP IOM %8i %8i %8i %8i\n",
             iom[i].Region.Protected.Begin, iom[i].Region.Protected.End, iom[i].Region.Sequence.Begin,
             iom[i].Region.Sequence.End);
      printf("XALIP IOI %8i %8i %8i %8i\n",
             ioi[i].Region.Protected.Begin, ioi[i].Region.Protected.End, ioi[i].Region.Sequence.Begin,
             ioi[i].Region.Sequence.End);
   }
#endif

   // Initiate search control values
   int LastCheckPoint = iseq;
   int FirstCheckPoint = iseq+1;
   int NextCheckPoint = LastCheckPoint + jlcp;

   // Move one sequence position forward
   ++iseq;

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Loop over sequence positions
   ////////////////////////////////////////////////////////////////////////////////////////////
   SeqPosLoop:
   {
      /////////////////////////////////////////////////////////////////////////////////////////
      // Update work array
      if (iseq < lseq ) {
      #pragma noinline
         nextR(prf, Sequence, iseq, iop, iom, ioi, lseq, &lAlignment, &ifer, Lock[iseq], N1, N2);
      } else {
      #pragma noinline
         nextR_last(prf, Sequence, iseq, iop, iom, ioi, lseq, &lAlignment, &ifer, Lock[iseq], N1, N2);
      }
      /////////////////////////////////////////////////////////////////////////////////////////
      // If Match found
      if (lAlignment.Score >= kcut) {
         // Determine first entry of current row
         if ( (ifer._a > lAlignment.Region.Protected.End) || (iseq == lseq) ) {
            // Fill in missing alignment data
            lAlignment.Region.Protected.Begin = lAlignment.Region.Protected.Begin == 0 ? lAlignment.Region.Sequence.Begin : lAlignment.Region.Protected.Begin;
            lAlignment.Region.Protected.End = lAlignment.Region.Protected.End == 0     ? iseq                             : lAlignment.Region.Protected.End;


            // Check for errors
            if (lAlignment.Region.Protected.End < lAlignment.Region.Protected.Begin) {
	       fprintf(stderr, "Error: Illegal alignment found in alignment %lu - no list produced.\n"
	                      "       Alignement should be from %i to %i!\n",1+AlignmentNumber, lAlignment.Region.Protected.Begin, lAlignment.Region.Protected.End );
               return -1;
            }

            if (++AlignmentNumber >= MaxNumberOfAlignment) {
               fputs("Warning: Too many alignments found - list may be incomplete.\n", stderr);
               return -2;
            }

            // Accept alignment
            struct Alignment * const restrict pAlignment = &alignment[AlignmentNumber];
            pAlignment->Region.Sequence.Begin  = lAlignment.Region.Sequence.Begin;
            pAlignment->Region.Protected.Begin = lAlignment.Region.Protected.Begin;
            pAlignment->Region.Protected.End   = lAlignment.Region.Protected.End;
            pAlignment->Region.Sequence.End    = lAlignment.Region.Sequence.End;
            pAlignment->Score                  = lAlignment.Score;
#ifdef XALIP_DEBUG
            printf("XALIP ALIGN %lu %4.4i %4.4i %4.4i %4.4i %4.4i\n",AlignmentNumber,
                   lAlignment.Region.Protected.Begin, lAlignment.Region.Protected.End, lAlignment.Score, lAlignment.Region.Sequence.Begin,
                   lAlignment.Region.Sequence.End);
#endif
            // Protect sequence region
            for (int jseq=lAlignment.Region.Protected.Begin; jseq<=lAlignment.Region.Protected.End; ++jseq) {
               Lock[jseq] = true;
            }

            // Exit if only searching for optimal alignment
            if (AlignmentNumber>0 && Lopt)
               return AlignmentNumber;
            else
               goto MajorLoop;
         } else {
            if ( ++iseq <= lseq ) goto SeqPosLoop;
         }
      } else {
         // Have we reached next check point ?
         if (iseq >= NextCheckPoint) {
            // Determine firdst entry of current row
            if (ifer._b >= LastCheckPoint) {
               ibeg = FirstCheckPoint - 1;
               FirstCheckPoint = ifer._c;
               LastCheckPoint = iseq;
            }

            // Calculate next check point
            NextCheckPoint += jlcp;
         }

         // Move one sequence position forward
         if ( ++iseq <= lseq ) goto SeqPosLoop;
      }
   }
   return AlignmentNumber;
}

