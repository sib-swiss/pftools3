/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 26, 2011 xalip_sse41.c
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
#include <assert.h>

#include "../include/pfProfile.h"
#include "sse41_inline_fcts.h"

#define MAX(a,b) (a>b) ? a : b
#define MIN(a,b) (a<b) ? a : b

struct Minima {
   int ProtectedRegionBegin;
   int ProtectedRegionBeginAfterN1;
   int SequenceBegin;
};

static void InitR(const int iseq, const size_t N1, const size_t N2, const size_t bseq, const size_t lseq,
                  union lScores * const restrict iop, union Positions * const restrict iom,
                  union Positions * const restrict ioi,  const struct Profile * const restrict prf)
{
  int KOPD;
#ifdef XALIP_DEBUG
  fprintf(stdout,"InitR called with index %i bseq %lu\n",iseq, bseq);
#endif
  // Are we treating sequence index below given start?
  register const ScoreTuple * restrict FirstSequenceProtein = ( iseq < bseq )
      ? &(prf->Scores.Insertion.FirstSequenceProtein[0])
      : &(prf->Scores.Insertion.Transitions->From[EXTRA]);
  const size_t FirstSequenceProteinAlignStep = (iseq < bseq) ? 1 : 4;

//   register __m128i * restrict pIOP = &(iop[0].xmm);

  // Load FirstSequenceProtein
  __m128i __FirstSequenceProtein = LoadStoredIntegerVector(FirstSequenceProtein);

  // Store into iop
  iop[0].xmm = __FirstSequenceProtein;

  // Set KOPD
  KOPD = _mm_extract_epi32(__FirstSequenceProtein, DELETION);

  register const TransitionScores * const restrict Transitions = prf->Scores.Insertion.Transitions;
  register const StoredIntegerFormat * restrict pMatch       = &prf->Scores.Match.Alphabet[_D];
  const size_t AlignStep = prf->Scores.Match.AlignStep;

  // Move to next profile First Sequence
  FirstSequenceProtein += FirstSequenceProteinAlignStep;

  for (size_t iprf=1; iprf<=prf->Length; ++iprf) {
      register const int KD = KOPD + (int) *pMatch;
      pMatch += AlignStep;

      // Transform KD into a vector
      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __Transitions = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Add KD to Transitions
      __Transitions = _mm_add_epi32(__Transitions, __KD);

      // Load FirstSequenceProtein and Convert signed WORD into signed DWORD
      __m128i __FirstSequenceProtein = LoadStoredIntegerVector(FirstSequenceProtein);

      // Move to next profile First Sequence
      FirstSequenceProtein += FirstSequenceProteinAlignStep;

      // Get maximum
      __m128i __max = _mm_max_epi32(__Transitions, __FirstSequenceProtein);

      // Store IOPI and IOPM
      iop[iprf].xmm = __max;


      // Set KOPD
      KOPD = _mm_extract_epi32(__max, DELETION);
  }
#ifdef XALIP_DEBUG
  for (size_t iprf=0; iprf<=prf->Length; ++iprf) {
      printf("INITR IOP %5lu %15i %15i %15i %15i\n", iprf,
	     iop[iprf].Element[MATCH],
	     iop[iprf].Element[INSERTION],
	     iop[iprf].Element[DELETION],
	     iop[iprf].Element[EXTRA]
      );
  }
#endif

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
    TPOS.Region.Protected.End   = TPOS.Region.Sequence.Begin;
    iom[N1-1].xmm               = TPOS.xmm;
    for (size_t iprf=N1; iprf<N2; ++iprf) {
      _mm_store_si128(&(iom[iprf].xmm), TPOS.xmm);
      _mm_store_si128(&(ioi[iprf].xmm), TPOS.xmm);
    }
    TPOS.Region.Protected.Begin = lseq + 1;
    TPOS.Region.Protected.End   = 0;

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
                  struct Minima * const restrict minima, const _Bool isProtected, const size_t N1, const size_t N2)
{
#ifdef XALIP_DEBUG
   fprintf(stdout,"NextR called with iseq %i\n",iseq);
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
   if (isProtected) {
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
	 union Positions Kiod;
   const union Positions * PTRpos[4] = {&TEMPpos,  &Kpos, NULL, &Kiod };
//    union Positions Kiod;
//    PTRpos[0] = &TEMPpos;
//    PTRpos[1] = &Kpos;
//    PTRpos[3] = &Kiod;

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
      iom[0].xmm            = ioi[0].xmm;
   } else {
      iop[0].Element[MATCH] = itmp;
      iom[0].xmm            = TEMPpos.xmm;
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
      ioi[0].xmm                = TEMPpos.xmm;
   }

   // Initialize minima
   minima->ProtectedRegionBegin = iseq;
   minima->ProtectedRegionBeginAfterN1 = iseq;
   itmp                  = MIN(ioi[0].Region.Sequence.Begin, iom[0].Region.Sequence.Begin);
   minima->SequenceBegin = MIN(itmp, Kiod.Region.Sequence.Begin);

   // Initialize alignment
   union Positions Fpos __attribute__((aligned(16)));;
//    Fpos.Region.Protected.Begin = alignment->Region.Protected.Begin;
//    Fpos.Region.Protected.End   = alignment->Region.Protected.End;
   Fpos.Region.Protected.mm    = alignment->Region.Protected.mm;
   Fpos.Region.Sequence.Begin  = alignment->Region.Sequence.Begin;
   Fpos.Region.Sequence.End    = 0; // this is dummy

#ifdef XALIP_DEBUG
    printf("NEXTR IOP      0 %15i %15i %15i %15i\n",
	   iop[0].Element[MATCH],
	   iop[0].Element[INSERTION],
	   iop[0].Element[DELETION],
	   iop[0].Element[EXTRA]
    );
#endif

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

#ifdef XALIP_DEBUG
      _mm_empty();
      printf("NEXTR IOP M %4lu %15i %15i %15i %15i %15i\n", iprf,
	     _mm_extract_epi32(__TransitionsM, MATCH),
	     _mm_extract_epi32(__TransitionsI, MATCH),
	     _mm_extract_epi32(__TransitionsD, MATCH),
	     _mm_extract_epi32(__TransitionsX, MATCH),
	     (int) Match[SequenceIndex]
	    );
      printf("NEXTR IOP I %4lu %15i %15i %15i %15i %15i\n", iprf,
	     _mm_extract_epi32(__TransitionsM, INSERTION),
	     _mm_extract_epi32(__TransitionsI, INSERTION),
	     _mm_extract_epi32(__TransitionsD, INSERTION),
	     _mm_extract_epi32(__TransitionsX, INSERTION),
	     (int) Insertion[SequenceIndex]
	    );
      printf("NEXTR IOP D %4lu %15i %15i %15i %15i %15i\n", iprf,
	     _mm_extract_epi32(__TransitionsM, DELETION),
	     _mm_extract_epi32(__TransitionsI, DELETION),
	     _mm_extract_epi32(__TransitionsD, DELETION),
	     _mm_extract_epi32(__TransitionsX, DELETION),
	     (int) Match[_D]
	    );
      printf("NEXTR IOP X %4lu %15i %15i %15i %15i %15i\n", iprf,
	     _mm_extract_epi32(__TransitionsM, EXTRA),
	     _mm_extract_epi32(__TransitionsI, EXTRA),
	     _mm_extract_epi32(__TransitionsD, EXTRA),
	     _mm_extract_epi32(__TransitionsX, EXTRA),
	     0
	    );
#endif

      // Move to next profile index
      Match     += AlignStep;
      Insertion += AlignStep;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Compute Maxima (Fortran Nstep function)
      __m128i __Index = _mm_setzero_si128();
      __m128i __three = _mm_set_epi32(3,3,3,3);

      __m128i  __Mask = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
      __TransitionsX  = _mm_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
      __Index         = _mm_blendv_epi8(__Index, __three, __Mask);

      __m128i __One  = _mm_set_epi32(1,1,1,1);
      __Mask         = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
      __TransitionsX = _mm_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
      __Index        = _mm_blendv_epi8(__Index, __One, __Mask);

      __m128i __Two  = _mm_add_epi32(__One, __One);
      __Mask         = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
      __TransitionsX = _mm_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
      __Index        = _mm_blendv_epi8(__Index, __Two, __Mask);

      // Set new data
      iop[iprf].xmm = __TransitionsX;
      Kopd = _mm_extract_epi32(__TransitionsX, DELETION);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Check for new maxima
      const int KE = _mm_extract_epi32(__TransitionsX, DUMMY);
      if (KE > alignment->Score) {
         alignment->Score = KE;
         alignment->Region.Sequence.End = iseq;
//          Fpos.xmm = ioi[iprf].xmm;
         const unsigned int Id = (unsigned int) _mm_extract_epi32(__Index, DUMMY);
				 assert(Id <4);
         if (Id == 1)  { // KM is max
						Fpos.xmm = Kpos.xmm;
         } else if (Id == 3) { // KD is max
						Fpos.xmm = Kiod.xmm;
         } else if (Id == 2) {
						Fpos.xmm = ioi[iprf].xmm;
				 } else {
						Fpos.xmm = TEMPpos.xmm;
				 }
      }
#ifdef XALIP_DEBUG
      _mm_empty();
      printf("FPOS %6i %4lu   %8i %8i %8i\n", iseq, iprf, Fpos.Region.Protected.Begin, Fpos.Region.Protected.End, Fpos.Region.Sequence.Begin);
#endif
      /////////////////////////////////////////////////////////////////////////////////////////
      // Update alignment positions
      union Positions Jpos __attribute__((aligned(16)));
      Jpos.xmm  = iom[iprf].xmm;
      PTRpos[2] = &ioi[iprf];

      const int NewM = _mm_extract_epi32(__Index, MATCH);
      iom[iprf].xmm  = PTRpos[NewM]->xmm;

      const int NewD = _mm_extract_epi32(__Index, DELETION);
      /* There is a need to protect Kiod in case of use later for IOI */
      const __m128i __tmp = PTRpos[NewD]->xmm;

      const int NewI = _mm_extract_epi32(__Index, INSERTION);
      ioi[iprf].xmm  = PTRpos[NewI]->xmm;

      Kpos.xmm = Jpos.xmm;
      Kiod.xmm = __tmp;

#ifdef XALIP_DEBUG
      _mm_empty();
			const char *Prov[4] = { "JUMPIN", "MATCH", "INSERTION", "DELETION" };
			printf("PROV %6i %4lu %s %s %s\n", iseq, iprf, Prov[NewM], Prov[NewI], Prov[NewD]);
      printf("FSEQBEG %6i %4lu   %8i %8i %8i\n", iseq, iprf,
	     iom[iprf].Region.Sequence.Begin, ioi[iprf].Region.Sequence.Begin, Kiod.Region.Sequence.Begin);
#endif

      /////////////////////////////////////////////////////////////////////////////////////////
      // Update minima

      const int t1 = MIN(ioi[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.Begin);
      const int t2 = MIN(t1, Kiod.Region.Protected.Begin);
      minima->ProtectedRegionBegin = MIN(minima->ProtectedRegionBegin, t2);

      if (iprf > N1) {
         minima->ProtectedRegionBeginAfterN1 = MIN(minima->ProtectedRegionBeginAfterN1, t2);
      }
      const int t3 = MIN(ioi[iprf].Region.Sequence.Begin, iom[iprf].Region.Sequence.Begin);
      const int t4 = MIN(t3, Kiod.Region.Sequence.Begin);
      minima->SequenceBegin = MIN(minima->SequenceBegin, t4);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Epilogue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Finish updating alignment positions
   alignment->Region.Protected.mm    = Fpos.Region.Protected.mm;
   alignment->Region.Sequence.Begin  = Fpos.Region.Sequence.Begin;

   // Entry and exit from protected regions
   iom[N1-1].Region.Protected.Begin = MIN(iom[N1-1].Region.Protected.Begin, In);
   iom[N1-1].Region.Protected.End   = In;

   for (size_t iprf=N1; iprf<N2; ++iprf) {
      iom[iprf].Region.Protected.Begin = MIN(iom[iprf].Region.Protected.Begin, In);
      iom[iprf].Region.Protected.End   = In;

      ioi[iprf].Region.Protected.Begin = MIN(ioi[iprf].Region.Protected.Begin, In);
      ioi[iprf].Region.Protected.End   = In;
   }
	_mm_empty();

#ifdef XALIP_DEBUG
   for (size_t iprf=0; iprf<=prf->Length; ++iprf) {
      printf("NEXTR IOP %4.4lu %15i %15i %15i\n", iprf, iop[iprf].Element[MATCH], iop[iprf].Element[INSERTION],
             iop[iprf].Element[DELETION]);
      printf("NEXTR IOM %4.4lu %15i %15i %15i\n", iprf, iom[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.End,
             iom[iprf].Region.Sequence.Begin);
      printf("NEXTR IOI %4.4lu %15i %15i %15i\n", iprf, ioi[iprf].Region.Protected.Begin, ioi[iprf].Region.Protected.End,
             ioi[iprf].Region.Sequence.Begin);

   }

   printf("NEXTR ALIGN %5i %4i %4i %7i %4i %4i\n", iseq,
          alignment->Region.Protected.Begin, alignment->Region.Protected.End, alignment->Score,
		      alignment->Region.Sequence.Begin, alignment->Region.Sequence.End);
#endif

}

static void nextR_last(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                       const int iseq, union lScores * const restrict iop, union Positions * const restrict iom,
                       union Positions * const restrict ioi,const int lseq, struct Alignment * const restrict alignment,
                       struct Minima * const restrict minima, const _Bool isProtected, const size_t N1, const size_t N2)
{
#ifdef XALIP_DEBUG
   fprintf(stdout,"NextR_last called with iseq %i\n",iseq);
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
   if (isProtected) {
      iop[N1-1].Element[MATCH] = NLOW;
      for (int iprf=N1; iprf<N2; ++iprf) {
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
   union Positions Kiod __attribute__((aligned(16)));
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
   minima->ProtectedRegionBegin = iseq;
   minima->ProtectedRegionBeginAfterN1 = iseq;
   itmp     = MIN(ioi[0].Region.Sequence.Begin, iom[0].Region.Sequence.Begin);
   minima->SequenceBegin = MIN(itmp, Kiod.Region.Sequence.Begin);

   // Initialize alignment
   union Positions Fpos __attribute__((aligned(16)));
   Fpos.Region.Protected.mm    = alignment->Region.Protected.mm;
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
      __TransitionsM = _mm_insert_epi32(__TransitionsM, (int) LastProteinSequence[iprf].From[MATCH], DUMMY);
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
      __TransitionsI = _mm_insert_epi32(__TransitionsI, (int) LastProteinSequence[iprf].From[INSERTION], DUMMY);
      // Add KM to Transitions
      __TransitionsI = _mm_add_epi32(__TransitionsI, __KI);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Deletion
      const register int KD = Kopd + (int) Match[_D];

      __m128i __KD = _mm_set1_epi32(KD);
      // Load Transitions and Convert signed WORD into signed DWORD
      __m128i __TransitionsD = LoadStoredIntegerVector(&(Transitions[iprf].From[DELETION]));
      // Insert LastProteinSequence
      __TransitionsD = _mm_insert_epi32(__TransitionsD, (int) LastProteinSequence[iprf].From[DELETION], DUMMY);
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
      __m128i __Index = (__m128i) _mm_setzero_ps();
      __m128i __three = _mm_set_epi32(3,3,3,3);

      __m128i  __Mask = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
      __TransitionsX  = _mm_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
      __Index         = _mm_blendv_epi8(__Index, __three, __Mask);

      __m128i __One  = _mm_set_epi32(1,1,1,1);
      __Mask         = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
      __TransitionsX = _mm_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
      __Index        = _mm_blendv_epi8(__Index, __One, __Mask);

      __m128i __Two  = _mm_add_epi32(__One, __One);
      __Mask         = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
      __TransitionsX = _mm_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
      __Index        = _mm_blendv_epi8(__Index, __Two, __Mask);

      // Set new data
      iop[iprf].xmm = __TransitionsX;
      Kopd = _mm_extract_epi32(__TransitionsX, DELETION);

      /////////////////////////////////////////////////////////////////////////////////////////
      // Check for new maxima
      const int KE = _mm_extract_epi32(__TransitionsX, DUMMY);
      if (KE > alignment->Score) {
         alignment->Score = KE;
         alignment->Region.Sequence.End = iseq;
         Fpos.xmm = ioi[iprf].xmm;
				 const unsigned int Id = (unsigned int) _mm_extract_epi32(__Index, DUMMY);
				 assert(Id <4);
         if (Id == 1)  { // KM is max
						Fpos.xmm = Kpos.xmm;
         } else if (Id == 3) { // KD is max
						Fpos.xmm = Kiod.xmm;
         } else if (Id == 2) {
						Fpos.xmm = ioi[iprf].xmm;
				 } else {
						Fpos.xmm = TEMPpos.xmm;
				 }
      }
#ifdef XALIP_DEBUG
      _mm_empty();
      printf("FPOS %8i %8i %8i\n", Fpos.Region.Protected.Begin, Fpos.Region.Protected.End, Fpos.Region.Sequence.Begin);
#endif
      /////////////////////////////////////////////////////////////////////////////////////////
      // Update alignment positions
      union Positions Jpos __attribute__((aligned(16)));
      Jpos.xmm  = iom[iprf].xmm;
      PTRpos[2] = &ioi[iprf];

      const int NewM = _mm_extract_epi32(__Index, MATCH);
      iom[iprf].xmm  = PTRpos[NewM]->xmm;

      const int NewD = _mm_extract_epi32(__Index, DELETION);
      Kiod.xmm       = PTRpos[NewD]->xmm;

      const int NewI = _mm_extract_epi32(__Index, INSERTION);
      ioi[iprf].xmm  = PTRpos[NewI]->xmm;

      Kpos.xmm = Jpos.xmm;

      /////////////////////////////////////////////////////////////////////////////////////////
      // Update minima

      const int t1 = MIN(ioi[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.Begin);
      const int t2 = MIN(t1, Kiod.Region.Protected.Begin);
      minima->ProtectedRegionBegin     = MIN(minima->ProtectedRegionBegin, t2);

      if (iprf > N1) {
         minima->ProtectedRegionBeginAfterN1 = MIN(minima->ProtectedRegionBeginAfterN1, t2);
      }
      const int t3 = MIN(ioi[iprf].Region.Sequence.Begin, iom[iprf].Region.Sequence.Begin);
      const int t4 = MIN(t3, Kiod.Region.Sequence.Begin);
      minima->SequenceBegin = MIN(minima->SequenceBegin, t4);
   }

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Epilogue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Finish updating alignment positions
//    alignment->Region.Protected.Begin = Fpos.Region.Protected.Begin;
//    alignment->Region.Protected.End   = Fpos.Region.Protected.End;
   alignment->Region.Protected.mm    = Fpos.Region.Protected.mm;
   alignment->Region.Sequence.Begin  = Fpos.Region.Sequence.Begin;

   // Entry and exit from protected regions
   iom[N1-1].Region.Protected.Begin = MIN(iom[N1-1].Region.Protected.Begin, In);
   iom[N1-1].Region.Protected.End = In;

   for (int iprf=N1; iprf<N2; ++iprf) {
      iom[iprf].Region.Protected.Begin = MIN(iom[iprf].Region.Protected.Begin, In);
      iom[iprf].Region.Protected.End   = In;

      ioi[iprf].Region.Protected.Begin = MIN(ioi[iprf].Region.Protected.Begin, In);
      ioi[iprf].Region.Protected.End   = In;
   }
   _mm_empty();

#ifdef XALIP_DEBUG
   for (size_t iprf=0; iprf<=prf->Length; ++iprf) {
      printf("NEXTR LAST IOP %4.4lu %15i %15i %15i\n", iprf, iop[iprf].Element[MATCH], iop[iprf].Element[INSERTION],
             iop[iprf].Element[DELETION]);
      printf("NEXTR LAST IOM %4.4lu %15i %15i %15i\n", iprf, iom[iprf].Region.Protected.Begin, iom[iprf].Region.Protected.End,
             iom[iprf].Region.Sequence.Begin);
      printf("NEXTR LAST IOI %4.4lu %15i %15i %15i\n", iprf, ioi[iprf].Region.Protected.Begin, ioi[iprf].Region.Protected.End,
             ioi[iprf].Region.Sequence.Begin);

   }

   printf("NEXTR LAST ALIGN %4i %4i %7i %4i %4i\n",
                   alignment->Region.Protected.Begin, alignment->Region.Protected.End, alignment->Score,
		   alignment->Region.Sequence.Begin, alignment->Region.Sequence.End);
#endif

}

int xalip_sse41( const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
           union lScores * const restrict iop, union Positions * const restrict iom,
           union Positions * const restrict ioi, const int bseq, const int lseq,
           struct Alignment * const restrict alignment,
           _Bool * const restrict ProtectedRegion, const size_t N1, const size_t N2, const _Bool Lopt,
           const int Cutoff, const size_t MaxNumberOfAlignment)
{
   int iseq;
   ////////////////////////////////////////////////////////////////////////////////////////////
   // Prologue
   ////////////////////////////////////////////////////////////////////////////////////////////

   // Alignment list
   size_t AlignmentNumber = 0;

   // Search control fields
   int SearchSequenceBegin = bseq-1;
   const int SearchSequenceJump = ((int)prf->Length) / 2;

   // Stack Memory
   struct Minima minima __attribute__((aligned(16)));

   ////////////////////////////////////////////////////////////////////////////////////////////
   // Two step forward one step backward loop
   ////////////////////////////////////////////////////////////////////////////////////////////
   MajorLoop:

   iseq = SearchSequenceBegin;
   struct Alignment lAlignment __attribute__((aligned(16)));
   lAlignment.Region.xmm = _mm_setzero_si128();
   lAlignment.Score = NLOW;


   // Initiate work array
   InitR(iseq, N1, N2, bseq, lseq, iop, iom, ioi, prf);
#ifdef XALIP_DEBUG
   for (size_t i=0; i<=prf->Length; ++i) {
      printf("IOP %8i %8i %8i %8i\n",
             iop[i].Element[MATCH], iop[i].Element[INSERTION], iop[i].Element[DELETION],
             iop[i].Element[DUMMY]);
      printf("IOM %8i %8i %8i %8i\n",
             iom[i].Region.Protected.Begin, iom[i].Region.Protected.End, iom[i].Region.Sequence.Begin,
             iom[i].Region.Sequence.End);
      printf("IOI %8i %8i %8i %8i\n",
             ioi[i].Region.Protected.Begin, ioi[i].Region.Protected.End, ioi[i].Region.Sequence.Begin,
             ioi[i].Region.Sequence.End);
   }
#endif

   // Initiate search control values
   int LastCheckPoint  = SearchSequenceBegin;
   int FirstCheckPoint = SearchSequenceBegin + 1;
   int NextCheckPoint  = SearchSequenceBegin + SearchSequenceJump;

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
				nextR(prf, Sequence, iseq, iop, iom, ioi, lseq, &lAlignment, &minima, ProtectedRegion[iseq], N1, N2);
      } else {
				#pragma noinline
				nextR_last(prf, Sequence, iseq, iop, iom, ioi, lseq, &lAlignment, &minima, ProtectedRegion[iseq], N1, N2);
      }
#ifdef XALIP_DEBUG
			printf("MINIMA %6i %i\n", iseq, minima.SequenceBegin);
#endif
      /////////////////////////////////////////////////////////////////////////////////////////
      // If Match found
      if (lAlignment.Score >= Cutoff) {
         // Determine first entry of current row
         if ( (minima.ProtectedRegionBegin > lAlignment.Region.Protected.End) || (iseq == lseq) ) {
            // Fill in missing alignment data
            lAlignment.Region.Protected.Begin = lAlignment.Region.Protected.Begin == 0 ? lAlignment.Region.Sequence.Begin : lAlignment.Region.Protected.Begin;
            lAlignment.Region.Protected.End   = lAlignment.Region.Protected.End == 0   ? iseq                             : lAlignment.Region.Protected.End;

            // Check for errors
            if (lAlignment.Region.Protected.End < lAlignment.Region.Protected.Begin) {
							fprintf(stderr, "Error: Inconsistent alignment found in alignment %lu - no list produced.\n"
	                            "       Alignement should be from %i to %i!\n",
							        1+AlignmentNumber, lAlignment.Region.Protected.Begin, lAlignment.Region.Protected.End );
#ifdef XALIP_DEBUG
							printf("XALIP ERROR N1=%lu N2=%lu\n", N1, N2);
							printf("XALIP ERROR TOTAL SEQUENCE LENGTH %lu\n", lseq);
							printf("XALIP ERROR CUTOFF VALUE %i\n", Cutoff);
							printf("XALIP ERROR PROTECTED REGION %s\n", ProtectedRegion[iseq] ? "TRUE" : "FALSE");
							printf("XALIP ERROR SEQUENCE INDEX %lu\n", iseq);
							printf("XALIP ERROR MINIMA PROTECTED REGION BEGIN %i\n"
										 "XALIP ERROR MINIMA PROTECTED REGION END %i\n"
										 "XALIP ERROR MINIMA SEQUENCE BEGIN %i\n",
								     minima.ProtectedRegionBegin, minima.ProtectedRegionBeginAfterN1, minima.SequenceBegin);
							printf("XALIP ERROR ALIGNMENT PROTECTED REGION BEGIN %i\n"
								     "XALIP ERROR ALIGNMENT PROTECTED REGION END %i\n"
								     "XALIP ERROR ALIGNMENT SEQUENCE BEGIN %i\n"
								     "XALIP ERROR ALIGNMENT SEQUENCE END %i\n"
								     "XALIP ERROR ALIGNMENT SCORE %i\n"
										 "XALIP ERROR ALIGNMENT IPMB %i\n"
								     "XALIP ERROR ALIGNMENT IPME %i\n",
								lAlignment.Region.Protected.Begin, lAlignment.Region.Protected.End,
								lAlignment.Region.Sequence.Begin, lAlignment.Region.Sequence.End,
								lAlignment.Score, lAlignment.IPMB, lAlignment.IPME);
#endif
               return -1;
            }

            if (++AlignmentNumber >= MaxNumberOfAlignment) {
               fputs("Warning: Too many alignments found - list may be incomplete.\n", stderr);
                return (int) (MaxNumberOfAlignment-1);
            }

            // Accept alignment
            struct Alignment * const restrict pAlignment = &alignment[AlignmentNumber];
            pAlignment->Region.Sequence.Begin  = lAlignment.Region.Sequence.Begin;
            pAlignment->Region.Protected.Begin = lAlignment.Region.Protected.Begin;
            pAlignment->Region.Protected.End   = lAlignment.Region.Protected.End;
            pAlignment->Region.Sequence.End    = lAlignment.Region.Sequence.End;
            pAlignment->Score                  = lAlignment.Score;
#ifdef XALIP_DEBUG
            printf("XALIP ALIGN %lu %4.4i %4.4i %4.4i %4.4i %4.4i\n", AlignmentNumber,
                   lAlignment.Region.Protected.Begin, lAlignment.Region.Protected.End, lAlignment.Score, lAlignment.Region.Sequence.Begin, lAlignment.Region.Sequence.End);
#endif
            // Protect sequence region
            for (int jseq=lAlignment.Region.Protected.Begin; jseq<=lAlignment.Region.Protected.End; ++jseq) {
               ProtectedRegion[jseq] = true;
            }

            // Exit if only searching for optimal alignment
            if (AlignmentNumber>0 && Lopt)
               return (int) AlignmentNumber;
            else
               goto MajorLoop;
         }
         else {
            if ( ++iseq <= lseq ) goto SeqPosLoop;
         }
      }
      else {
         // Have we reached next check point ?
         if (iseq >= NextCheckPoint) {
            // Determine first entry of current row
            if (minima.ProtectedRegionBeginAfterN1 >= LastCheckPoint) {
               SearchSequenceBegin = FirstCheckPoint - 1;
               FirstCheckPoint     = minima.SequenceBegin;
               LastCheckPoint      = iseq;
            }

            // Calculate next check point
            NextCheckPoint += SearchSequenceJump;
         }

         // Move one sequence position forward
         if ( ++iseq <= lseq ) goto SeqPosLoop;
      }
   }
   return (int) AlignmentNumber;
}

