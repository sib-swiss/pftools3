/*******************************************************
                        PFTOOLS
 *******************************************************
  Dec 1, 2011 heuristic.c
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <emmintrin.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include "../include/pfProfile.h"

unsigned int TransposeHeuristic_sse2(const TransposeMatrix TransposeMatch, const size_t Alphabet_Length,
				     const size_t Profile_Length, const PFSequence * const restrict Sequence)
{
  size_t iprf;
  // Allocate vectors on the stack ( one for read, one for write )
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;

  float * restrict v0 = (float *) (( (uintptr_t) alloca(3*Aligned_Profile_Length*sizeof(float) + 63)) & ~63);
  float * restrict v1 = v0 + Aligned_Profile_Length;
  float * restrict Sc = v0 + 2*Aligned_Profile_Length;

//   float * restrict v0 = (float *) ( ( (uintptr_t) alloca(Aligned_Profile_Length*sizeof(float) + 63)) & ~63);
//   float * restrict v1 = (float *) ( ( (uintptr_t) alloca(Aligned_Profile_Length*sizeof(float) + 63)) & ~63);
//   float * restrict Sc = (float *) ( ( (uintptr_t) alloca(Aligned_Profile_Length*sizeof(float) + 63)) & ~63);

  register float * restrict v_w = v0;

  // Initialize v with first sequence
  register size_t Index = (size_t) Sequence->ProfileIndex[0];
  register const float * restrict lMatch = &TransposeMatch.f[Aligned_Profile_Length*Index];

  iprf = 0;
  register __m128 __Zero1, __Zero2, __Zero3, __Zero4;
//    __asm__ __volatile__ (
// 	  "xorps %0, %0;"
// 	  "xorps %1, %1;"
// 	  "xorps %2, %2;"
// 	  "xorps %3, %3;"
// 	  : "=x"(__Zero1) , "=x"(__Zero2), "=x"(__Zero3), "=x"(__Zero4)
//     );
  do {
    __m128 __sc1, __sc2, __sc3, __sc4;
    __m128 __m1, __m2, __m3, __m4;

    __m1  = _mm_load_ps(&lMatch[iprf]);
    __m2  = _mm_load_ps(&lMatch[iprf+4]);
    __m3  = _mm_load_ps(&lMatch[iprf+8]);
    __m4  = _mm_load_ps(&lMatch[iprf+12]);

    __asm__ __volatile__ ("xorps %0, %0;" : "=x"(__Zero1) );
    __asm__ __volatile__ ("movaps %1, %0; maxps %2, %0;" : "=x"(__sc1) : "x"(__m1), "x"(__Zero1));
    __asm__ __volatile__ ("movaps %1, %0; maxps %2, %0;" : "=x"(__sc2) : "x"(__m2), "x"(__Zero1));
    __asm__ __volatile__ ("movaps %1, %0; maxps %2, %0;" : "=x"(__sc3) : "x"(__m3), "x"(__Zero1));
    __asm__ __volatile__ ("movaps %1, %0; maxps %2, %0;" : "=x"(__sc4) : "x"(__m4), "x"(__Zero1));

//     __sc1 = _mm_max_ps(__Zero1, __m1);
//     __sc2 = _mm_max_ps(__Zero1, __m2);
//     __sc3 = _mm_max_ps(__Zero1, __m3);
//     __sc4 = _mm_max_ps(__Zero1, __m4);

//     __Zero1 = _mm_max_ps(__Zero1, __m1);
//     __Zero2 = _mm_max_ps(__Zero2, __m2);
//     __Zero3 = _mm_max_ps(__Zero3, __m3);
//     __Zero4 = _mm_max_ps(__Zero4, __m4);

    _mm_store_ps(&v_w[iprf   ], __m1);
    _mm_store_ps(&v_w[iprf+4 ], __m2);
    _mm_store_ps(&v_w[iprf+8 ], __m3);
    _mm_store_ps(&v_w[iprf+12], __m4);

    _mm_store_ps(&Sc[iprf   ], __sc1);
    _mm_store_ps(&Sc[iprf+4 ], __sc2);
    _mm_store_ps(&Sc[iprf+8 ], __sc3);
    _mm_store_ps(&Sc[iprf+12], __sc4);

//     _mm_store_ps(&Sc[iprf   ], __Zero1);
//     _mm_store_ps(&Sc[iprf+4 ], __Zero2);
//     _mm_store_ps(&Sc[iprf+8 ], __Zero3);
//     _mm_store_ps(&Sc[iprf+12], __Zero4);
//
//      __asm__ __volatile__ (
// 	  "xorps %0, %0;"
// 	  "xorps %1, %1;"
// 	  "xorps %2, %2;"
// 	  "xorps %3, %3;"
// 	  : "=x"(__Zero1) , "=x"(__Zero2), "=x"(__Zero3), "=x"(__Zero4)
//     );

    iprf += 16;
  } while (iprf < Profile_Length);

  // Set read v pointer
  register const float * v_r = v0;
  v_w = v1;

  // Run through the rest of the profile
  for (unsigned int iseq=1; iseq<(unsigned int) Sequence->Length; ++iseq) {
    Index = (size_t) Sequence->ProfileIndex[iseq];
    lMatch = &TransposeMatch.f[Aligned_Profile_Length*Index];
//     for (iprf=Profile_Length; iprf<Aligned_Profile_Length; ++iprf) fprintf(stderr,"%2u %lu %lf %lf\n", iseq, iprf, lMatch[iprf], Sc[iprf]);
#if 0
    v_w[0] = lMatch[0] > 0.0f ? lMatch[0] : 0.0f;
    if (lMatch[0] > Sc[0] ) Sc[0] = lMatch[0];

    for (iprf=1; iprf<Profile_Length; ++iprf) {
      float tmp   = v_r[iprf-1] + lMatch[iprf];
      tmp = tmp > 0.0f ? tmp : 0.0f;
      v_w[iprf] = tmp;

      if (tmp > Sc[iprf]) Sc[iprf] = tmp;
    }
#else
    __m128 __V_R_0 = _mm_load_ps(&v_r[0]);
    __V_R_0        = (__m128) _mm_slli_si128((__m128i) __V_R_0, 4);

    iprf=0;
    goto Insert;

    Loop:
      __V_R_0        = _mm_loadu_ps(&v_r[iprf-1     ]);

    Insert:
    ;
      __m128 __V_R_1 = _mm_loadu_ps(&v_r[iprf-1 +  4]);
      __m128 __V_R_2 = _mm_loadu_ps(&v_r[iprf-1 +  8]);
      __m128 __V_R_3 = _mm_loadu_ps(&v_r[iprf-1 + 12]);

//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero1) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero2) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero3) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero4) );

      __V_R_0        = _mm_add_ps(__V_R_0, *((__m128*)&lMatch[iprf   ]));
      __V_R_1        = _mm_add_ps(__V_R_1, *((__m128*)&lMatch[iprf+ 4]));
      __V_R_2        = _mm_add_ps(__V_R_2, *((__m128*)&lMatch[iprf+ 8]));
      __V_R_3        = _mm_add_ps(__V_R_3, *((__m128*)&lMatch[iprf+12]));

//       __V_R_0        = _mm_max_ps( __Zero1, __V_R_0);
//       __V_R_1        = _mm_max_ps( __Zero1, __V_R_1);
//       __V_R_2        = _mm_max_ps( __Zero1, __V_R_2);
//       __V_R_3        = _mm_max_ps( __Zero1, __V_R_3);
      __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero1) );
      __asm__ __volatile__ ("maxps %2, %0;" : "=x"(__V_R_0) : "0"(__V_R_0), "x"(__Zero1) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero2) );
      __asm__ __volatile__ ("maxps %2, %0;" : "=x"(__V_R_1) : "0"(__V_R_1), "x"(__Zero1) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero3) );
      __asm__ __volatile__ ("maxps %2, %0;" : "=x"(__V_R_2) : "0"(__V_R_2), "x"(__Zero1) );
//       __asm__ __volatile__ ( "xorps %0, %0;" : "=x"(__Zero4) );
      __asm__ __volatile__ ("maxps %2, %0;" : "=x"(__V_R_3) : "0"(__V_R_3), "x"(__Zero1) );

//       __V_R_0        = _mm_max_ps( __Zero1, __V_R_0);
//       __V_R_1        = _mm_max_ps( __Zero2, __V_R_1);
//       __V_R_2        = _mm_max_ps( __Zero3, __V_R_2);
//       __V_R_3        = _mm_max_ps( __Zero4, __V_R_3);

      _mm_store_ps(&v_w[iprf   ], __V_R_0);
      _mm_store_ps(&v_w[iprf+4 ], __V_R_1);
      _mm_store_ps(&v_w[iprf+8 ], __V_R_2);
      _mm_store_ps(&v_w[iprf+12], __V_R_3);

//       __m128 __SC_0  = _mm_load_ps(&Sc[iprf   ]);
//       __m128 __SC_1  = _mm_load_ps(&Sc[iprf+ 4]);
//       __m128 __SC_2  = _mm_load_ps(&Sc[iprf+ 8]);
//       __m128 __SC_3  = _mm_load_ps(&Sc[iprf+12]);
//
//       __SC_0 = _mm_max_ps(__SC_0, __V_R_0);
//       __SC_1 = _mm_max_ps(__SC_1, __V_R_1);
//       __SC_2 = _mm_max_ps(__SC_2, __V_R_2);
//       __SC_3 = _mm_max_ps(__SC_3, __V_R_3);
//
//       _mm_store_ps(&Sc[iprf   ], __SC_0);
//       _mm_store_ps(&Sc[iprf+4 ], __SC_1);
//       _mm_store_ps(&Sc[iprf+8 ], __SC_2);
//       _mm_store_ps(&Sc[iprf+12], __SC_3);

      __asm__ __volatile__ (
	    "maxps   (%8,%9,4), %0;"
	    "maxps 16(%8,%9,4), %1;"
	    "maxps 32(%8,%9,4), %2;"
	    "maxps 48(%8,%9,4), %3;"
	    : "=x"(__V_R_0), "=x"(__V_R_1), "=x"(__V_R_2), "=x"(__V_R_3)
	    : "0"(__V_R_0), "1"(__V_R_1), "2"(__V_R_2), "3"(__V_R_3), "r"(Sc), "r"(iprf)
      );

      _mm_store_ps(&Sc[iprf   ], __V_R_0);
      _mm_store_ps(&Sc[iprf+4 ], __V_R_1);
      _mm_store_ps(&Sc[iprf+8 ], __V_R_2);
      _mm_store_ps(&Sc[iprf+12], __V_R_3);


      iprf += 16;

      if ( iprf<Profile_Length) goto Loop;

#endif

    // Swap pointers
    const float * ptr = v_w;
    v_w = (float*) v_r;
    v_r = ptr;
  }
#if 1
  unsigned int Score;
  {
    register __m128 __Sum, __Sum1, __Sum2, __Sum3;
    __asm__ __volatile__ ("xorps %0, %0;"
			  "xorps %1, %1;"
			  "xorps %2, %2;"
			  "xorps %3, %3;"
			  : "=x"(__Sum) , "=x"(__Sum1), "=x"(__Sum2), "=x"(__Sum3) );
    iprf = 0;
    do  {
      __Sum  = _mm_add_ps(__Sum,  *((__m128*) &Sc[iprf]));
      __Sum1 = _mm_add_ps(__Sum1, *((__m128*) &Sc[iprf+4]));
      __Sum2 = _mm_add_ps(__Sum2, *((__m128*) &Sc[iprf+8]));
      __Sum3 = _mm_add_ps(__Sum3, *((__m128*) &Sc[iprf+12]));
      iprf+=16;
    } while (iprf < (Profile_Length & ~0xF));

//     while ( iprf + 4 < Profile_Length & ~0x7) {
//       __Sum  = _mm_add_ps(__Sum,  *((__m128*) &Sc[iprf]));
//       iprf  += 4;
//     }

    __Sum  = _mm_add_ps(__Sum,  __Sum1);
    __Sum2 = _mm_add_ps(__Sum2, __Sum3);
    __Sum  = _mm_add_ps(__Sum,  __Sum2);

    __asm__ __volatile__ (
	      "movaps    %1, %2   ;"
	      "movhlps   %1, %2   ;"
	      "addps     %2, %1   ;"
	      "movaps    %1, %3   ;"
	      "shufps    $245, %1, %3 ;"
	      "addss     %3, %0 ;"
	      : "=x"(__Sum)
	      : "0"(__Sum), "x"(__Sum2), "x"(__Sum3)
	    );

    while ( iprf < Profile_Length ) {
      __asm__ __volatile__ ("addss (%1,%2,4), %0;" : "=x"(__Sum) : "r"(Sc), "r"(iprf) );
      ++iprf;
    }

    Score = (unsigned int) _mm_cvttss_si32(__Sum);
  }
  return Score;
#else
  float fScore = 0.0f;
  for (iprf=0; iprf<Profile_Length; ++iprf) fScore += Sc[iprf];
  return (unsigned int) fScore;
#endif

}

unsigned int TransposeHeuristicGivenMemory_sse2(const float * const restrict TransposeMatch, float * const Memory,
					        const size_t Alphabet_Length, const size_t Profile_Length,
						const PFSequence * const restrict Sequence)
{
  size_t iprf;
  float Score = 0.0f;

  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  float * restrict v0 = Memory;
  float * restrict v1 = Memory + Aligned_Profile_Length;
  float * restrict Sc = Memory + 2*Aligned_Profile_Length;

  register float * restrict v_w = v0;

  // Initialize v with first sequence
  register size_t Index = (size_t) Sequence->ProfileIndex[0];
  register const float * restrict lMatch = &TransposeMatch[Aligned_Profile_Length*Index];

  const register __m128 __Zero = _mm_setzero_ps();
  iprf = 0;
  do {
    __m128 __m1  = _mm_load_ps(&lMatch[iprf]);
    __m128 __m2  = _mm_load_ps(&lMatch[iprf+4]);
    __m128 __m3  = _mm_load_ps(&lMatch[iprf+8]);
    __m128 __m4  = _mm_load_ps(&lMatch[iprf+12]);
    __m128 __sc1 = _mm_max_ps(__Zero, __m1);
    _mm_store_ps(&v_w[iprf], __m1);
    __m128 __sc2 = _mm_max_ps(__Zero, __m2);
    _mm_store_ps(&v_w[iprf+4], __m1);
    __m128 __sc3 = _mm_max_ps(__Zero, __m3);
    _mm_store_ps(&v_w[iprf+8], __m1);
    __m128 __sc4 = _mm_max_ps(__Zero, __m4);
    _mm_store_ps(&v_w[iprf+12], __m1);
    _mm_store_ps(&Sc[iprf   ], __sc1);
    _mm_store_ps(&Sc[iprf+4 ], __sc2);
    _mm_store_ps(&Sc[iprf+8 ], __sc3);
    _mm_store_ps(&Sc[iprf+12], __sc4);
    iprf += 16;
  } while (iprf < Profile_Length);


  // Set read v pointer
  register const float * v_r = v0;
  v_w = v1;

  // Run through the rest of the profile
  for (unsigned int iseq=1; iseq<(unsigned int) Sequence->Length; ++iseq) {
    Index = (size_t) Sequence->ProfileIndex[iseq];
    lMatch = &TransposeMatch[Aligned_Profile_Length*Index];
#if 0
    v_w[0] = lMatch[0] > 0.0f ? lMatch[0] : 0.0f;
    if (lMatch[0] > Sc[0] ) Sc[0] = lMatch[0];

    __assume_aligned(lMatch, 16);
    __assume_aligned(v_r, 16);
    __assume_aligned(v_w, 16);
    __assume_aligned(Sc, 16);

    for (iprf=1; iprf<Profile_Length; ++iprf) {
      float tmp   = v_r[iprf-1] + lMatch[iprf];
      tmp = tmp > 0.0f ? tmp : 0.0f;
      v_w[iprf] = tmp;

      if (tmp > Sc[iprf]) Sc[iprf] = tmp;
    }
#else

    __m128 __V_R_0 = _mm_load_ps(&v_r[0]);
    __V_R_0        = (__m128) _mm_slli_si128((__m128i) __V_R_0, 4);

    iprf=0;
    goto Insert;

    Loop:
      __V_R_0                  = _mm_loadu_ps(&v_r[iprf-1]);

    Insert:
    ;
      __V_R_0        = _mm_add_ps(__V_R_0, *((__m128*)&lMatch[iprf]));

      __m128 __V_R_1 = _mm_loadu_ps(&v_r[iprf-1 + 4]);
      __V_R_1        = _mm_add_ps(__V_R_1, *((__m128*)&lMatch[iprf+4]));
      _mm_store_ps(&v_w[iprf   ], __V_R_0);

      __m128 __V_R_2 = _mm_loadu_ps(&v_r[iprf-1 + 8]);
      __V_R_2        = _mm_add_ps(__V_R_2, *((__m128*)&lMatch[iprf+8]));
       _mm_store_ps(&v_w[iprf+4 ], __V_R_1);

      __m128 __V_R_3 = _mm_loadu_ps(&v_r[iprf-1 + 12]);
      __V_R_3        = _mm_add_ps(__V_R_3, *((__m128*)&lMatch[iprf+12]));
      _mm_store_ps(&v_w[iprf+8 ], __V_R_2);

      __m128 __SC_0  = _mm_load_ps(&Sc[iprf]);
      __V_R_0        = _mm_max_ps( __Zero, __V_R_0);
      _mm_store_ps(&v_w[iprf+12], __V_R_3);


      __m128 __SC_1  = _mm_load_ps(&Sc[iprf + 4]);
      __V_R_1        = _mm_max_ps( __Zero, __V_R_1);
      __SC_0         = _mm_max_ps(__SC_0, __V_R_0);
      _mm_store_ps(&Sc[iprf   ], __SC_0);

      __m128 __SC_2  = _mm_load_ps(&Sc[iprf + 8]);
      __V_R_2        = _mm_max_ps( __Zero, __V_R_2);
       __SC_1        = _mm_max_ps(__SC_1, __V_R_1);
       _mm_store_ps(&Sc[iprf+4 ], __SC_1);

      __m128 __SC_3  = _mm_load_ps(&Sc[iprf + 12]);
      __V_R_3        = _mm_max_ps( __Zero, __V_R_3);
      __SC_2         = _mm_max_ps(__SC_2, __V_R_2);
       _mm_store_ps(&Sc[iprf+8 ], __SC_2);


      __SC_3 = _mm_max_ps(__SC_3, __V_R_3);
      _mm_store_ps(&Sc[iprf+12], __SC_3);

      iprf+= 16;

      if ( iprf<Profile_Length) goto Loop;

#endif

    // Swap pointers
    float * ptr = v_w;
    v_w = (float*) v_r;
    v_r = (const float*) ptr;

    // Update Score
    //Score += max;
  }

  for (iprf=0; iprf<Profile_Length; ++iprf)
   Score += Sc[iprf];

  return (unsigned int) Score;
}
