/*******************************************************
                        PFTOOLS
 *******************************************************
  Dec 1, 2011 heuristic_sse41.c
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
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include "../include/pfProfile.h"

unsigned int TransposeHeuristic_sse41(const TransposeMatrix TransposeMatch, const size_t Alphabet_Length,
				      const size_t Profile_Length, const PFSequence * const restrict Sequence)
// WARNING: Creation of the transpose matrix took care of zeroing extra data on cache line.
//          Remember that if you happen to change the code.
{
  register size_t iprf;
  // Allocate vectors on the stack ( one for read, one for write )
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;

  int * restrict v0 = (int *) ( ( (uintptr_t) alloca(3*Aligned_Profile_Length*sizeof(int) + 63)) & ~63);
  int * restrict v1 = v0 + Aligned_Profile_Length;
  int * restrict Sc = v0 + 2*Aligned_Profile_Length;

  register int * restrict v_w = v0;

  // Initialize v with first sequence
  register size_t Index = (size_t) Sequence->ProfileIndex[0];
  register const int * restrict lMatch = &TransposeMatch.i[Aligned_Profile_Length*Index];

#if 0
  for ( iprf=0; iprf<Profile_Length; ++iprf) {
    v_w[iprf] = lMatch[iprf];
    Sc[iprf]  = lMatch[iprf] > 0 ? lMatch[iprf] : 0;
  }
#else
  iprf = 0;
  do {
    register __m128i __Zero1; //, __Zero2, __Zero3, __Zero4;
    register __m128i __m1, __m2, __m3, __m4;
    register __m128i __sc1, __sc2, __sc3, __sc4;

    __m1  = _mm_load_si128((__m128i*)&lMatch[iprf   ]);
    __m2  = _mm_load_si128((__m128i*)&lMatch[iprf+ 4]);
    __m3  = _mm_load_si128((__m128i*)&lMatch[iprf+ 8]);
    __m4  = _mm_load_si128((__m128i*)&lMatch[iprf+12]);

    __asm__ __volatile__ ("pxor %0, %0;" : "=x"(__Zero1) );
    __asm__ __volatile__ ("movdqa %1, %0; pmaxsd %2, %0;" : "=x"(__sc1) : "x"(__m1), "x"(__Zero1));
    __asm__ __volatile__ ("movdqa %1, %0; pmaxsd %2, %0;" : "=x"(__sc2) : "x"(__m2), "x"(__Zero1));
    __asm__ __volatile__ ("movdqa %1, %0; pmaxsd %2, %0;" : "=x"(__sc3) : "x"(__m3), "x"(__Zero1));
    __asm__ __volatile__ ("movdqa %1, %0; pmaxsd %2, %0;" : "=x"(__sc4) : "x"(__m4), "x"(__Zero1));

    _mm_store_si128((__m128i*)&v_w[iprf   ], __m1);
    _mm_store_si128((__m128i*)&v_w[iprf+4 ], __m2);
    _mm_store_si128((__m128i*)&v_w[iprf+8 ], __m3);
    _mm_store_si128((__m128i*)&v_w[iprf+12], __m4);

    _mm_store_si128((__m128i*)&Sc[iprf   ], __sc1);
    _mm_store_si128((__m128i*)&Sc[iprf+4 ], __sc2);
    _mm_store_si128((__m128i*)&Sc[iprf+8 ], __sc3);
    _mm_store_si128((__m128i*)&Sc[iprf+12], __sc4);
    iprf += 16;
  } while (iprf < Profile_Length);

#endif
  // Set read v pointer
  register const int * v_r = v0;
  v_w = v1;

  // Run through the rest of the profile
  for (unsigned int iseq=1; iseq<(unsigned int) Sequence->Length; ++iseq) {
    Index = (size_t) Sequence->ProfileIndex[iseq];
    lMatch = &TransposeMatch.i[Aligned_Profile_Length*Index];

#if 0
    v_w[0] = lMatch[0] > 0 ? lMatch[0] : 0;
    if (lMatch[0] > Sc[0] ) Sc[0] = lMatch[0];

    for (size_t iprf=1; iprf<Profile_Length; ++iprf) {
      register int tmp   = v_r[iprf-1] + lMatch[iprf];
      tmp = tmp > 0 ? tmp : 0;
      v_w[iprf] = tmp;

      if (tmp > Sc[iprf]) Sc[iprf] = tmp;
    }
#else
    register __m128i __Zero1; //, __Zero2, __Zero3, __Zero4;
    register __m128i __V_R_1, __V_R_2, __V_R_3, __V_R_4;
    __V_R_1 = _mm_load_si128((__m128i*) &v_r[0]);
    __V_R_1 = _mm_slli_si128(__V_R_1, 4);

    iprf=0;
    goto Insert;

    Loop:
    __asm__ __volatile__ (".align 16; ");
      __V_R_1         = _mm_loadu_si128((__m128i*) &v_r[iprf-1]);

    Insert:
    ;
      __V_R_2 = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 4]);
      __V_R_3 = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 8]);
      __V_R_4 = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 12]);

      __asm__ __volatile__ ( "pxor %0, %0;" : "=x"(__Zero1) );

      __V_R_1 = _mm_add_epi32(__V_R_1, *((__m128i*)&lMatch[iprf]));
      __V_R_2 = _mm_add_epi32(__V_R_2, *((__m128i*)&lMatch[iprf+4]));
      __V_R_3 = _mm_add_epi32(__V_R_3, *((__m128i*)&lMatch[iprf+8]));
      __V_R_4 = _mm_add_epi32(__V_R_4, *((__m128i*)&lMatch[iprf+12]));

      __V_R_1 = _mm_max_epi32( __Zero1, __V_R_1);
      __V_R_2 = _mm_max_epi32( __Zero1, __V_R_2);
      __V_R_3 = _mm_max_epi32( __Zero1, __V_R_3);
      __V_R_4 = _mm_max_epi32( __Zero1, __V_R_4);

      _mm_store_si128((__m128i*)&v_w[iprf   ], __V_R_1);
      _mm_store_si128((__m128i*)&v_w[iprf+4 ], __V_R_2);
      _mm_store_si128((__m128i*)&v_w[iprf+8 ], __V_R_3);
      _mm_store_si128((__m128i*)&v_w[iprf+12], __V_R_4);

      __asm__ __volatile__ (
	    "pmaxsd   (%8,%9,4), %0;"
	    "pmaxsd 16(%8,%9,4), %1;"
	    "pmaxsd 32(%8,%9,4), %2;"
	    "pmaxsd 48(%8,%9,4), %3;"
	    : "=x"(__V_R_1), "=x"(__V_R_2), "=x"(__V_R_3), "=x"(__V_R_4)
	    : "0"(__V_R_1), "1"(__V_R_2), "2"(__V_R_3), "3"(__V_R_4), "r"(Sc), "r"(iprf)
      );

      _mm_store_si128((__m128i*)&Sc[iprf   ], __V_R_1);
      _mm_store_si128((__m128i*)&Sc[iprf+ 4], __V_R_2);
      _mm_store_si128((__m128i*)&Sc[iprf+ 8], __V_R_3);
      _mm_store_si128((__m128i*)&Sc[iprf+12], __V_R_4);

      iprf+= 16;

      if ( iprf<Profile_Length) goto Loop;

#endif

    // Swap pointers
    int * ptr = v_w;
    v_w = (int*) v_r;
    v_r = (const int*) ptr;
  }

  unsigned int Score = 0;
#if 0
  for (iprf=0; iprf<Profile_Length; ++iprf) Score += Sc[iprf];
  return Score;
#else
  {
    register __m128i __Sum, __Sum1, __Sum2, __Sum3;
    __asm__ __volatile__ ("pxor %0, %0;"
			  "pxor %1, %1;"
			  "pxor %2, %2;"
			  "pxor %3, %3;"
			  : "=x"(__Sum) , "=x"(__Sum1), "=x"(__Sum2), "=x"(__Sum3) );
    iprf = 0;
    do  {
      __Sum  = _mm_add_epi32(__Sum,  *((__m128i*) &Sc[iprf]));
      __Sum1 = _mm_add_epi32(__Sum1, *((__m128i*) &Sc[iprf+4]));
      __Sum2 = _mm_add_epi32(__Sum2, *((__m128i*) &Sc[iprf+8]));
      __Sum3 = _mm_add_epi32(__Sum3, *((__m128i*) &Sc[iprf+12]));
      iprf+=16;
    } while (iprf < (Profile_Length & ~0xF));

    __Sum  = _mm_add_epi32(__Sum,  __Sum1);
    __Sum2 = _mm_add_epi32(__Sum2, __Sum3);
    __Sum  = _mm_add_epi32(__Sum,  __Sum2);

    __asm__ __volatile__ (
	      "phaddd    %1, %1;"
	      "pshufd    $225, %1, %2;"
	      "paddd     %2, %1 ;"
	      "movd      %1, %0;"
	      : "=r"(Score)
	      : "x"(__Sum), "x"(__Sum2)
	    );

    while ( iprf < Profile_Length ) {
      Score += Sc[iprf];
      ++iprf;
    }
    return Score;
  }
#endif
}

unsigned int TransposeHeuristicGivenMemory_sse41(const int * const restrict TransposeMatch, int * const Memory,
						 const size_t Alphabet_Length, const size_t Profile_Length,
						 const PFSequence * const restrict Sequence)
// WARNING: Creation of the transpose matrix took care of zeroing extra data on cache line.
//          Remember that if you happen to change the code.
{
  unsigned int Score = 0;
  size_t iprf;
  // Allocate vectors on the stack ( one for read, one for write )
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  int * restrict v0 = Memory;
  int * restrict v1 = Memory + Aligned_Profile_Length;
  int * restrict Sc = Memory + 2*Aligned_Profile_Length;

  register int * restrict v_w = v0;

  // Initialize v with first sequence
  register size_t Index = (size_t) Sequence->ProfileIndex[0];
  register const int * restrict lMatch = &TransposeMatch[Aligned_Profile_Length*Index];

#if 0
  for ( iprf=0; iprf<Profile_Length; ++iprf) {
    v_w[iprf] = lMatch[iprf];
    Sc[iprf]  = lMatch[iprf] > 0 ? lMatch[iprf] : 0;
  }
  for ( iprf=Profile_Length; iprf<Aligned_Profile_Length; ++iprf) {
    Sc[iprf] = 0;
  }
#else
  const register __m128i __Zero  = _mm_setzero_si128();
  iprf = 0;
  do {
    __m128i __m1  = _mm_load_si128((__m128i*)&lMatch[iprf]);
    __m128i __m2  = _mm_load_si128((__m128i*)&lMatch[iprf+4]);
    __m128i __m3  = _mm_load_si128((__m128i*)&lMatch[iprf+8]);
    __m128i __m4  = _mm_load_si128((__m128i*)&lMatch[iprf+12]);
    __m128i __sc1 = _mm_max_epi32(__Zero, __m1);
    _mm_store_si128((__m128i*)&v_w[iprf], __m1);
    __m128i __sc2 = _mm_max_epi32(__Zero, __m2);
    _mm_store_si128((__m128i*)&v_w[iprf+4], __m1);
    __m128i __sc3 = _mm_max_epi32(__Zero, __m3);
    _mm_store_si128((__m128i*)&v_w[iprf+8], __m1);
    __m128i __sc4 = _mm_max_epi32(__Zero, __m4);
    _mm_store_si128((__m128i*)&v_w[iprf+12], __m1);
    _mm_store_si128((__m128i*)&Sc[iprf   ], __sc1);
    _mm_store_si128((__m128i*)&Sc[iprf+4 ], __sc2);
    _mm_store_si128((__m128i*)&Sc[iprf+8 ], __sc3);
    _mm_store_si128((__m128i*)&Sc[iprf+12], __sc4);
    iprf += 16;
  } while (iprf < Profile_Length);
  
  while (iprf < Aligned_Profile_Length) {
    _mm_store_si128((__m128i*)&Sc[iprf   ], __Zero);
    _mm_store_si128((__m128i*)&Sc[iprf+4 ], __Zero);
    _mm_store_si128((__m128i*)&Sc[iprf+8 ], __Zero);
    _mm_store_si128((__m128i*)&Sc[iprf+12], __Zero);
    iprf += 16;
  }
#endif
  // Set read v pointer
  register const int * v_r = v0;
  v_w = v1;

  // Run through the rest of the profile
  for (unsigned int iseq=1; iseq<(unsigned int) Sequence->Length; ++iseq) {
    Index = (size_t) Sequence->ProfileIndex[iseq];
    lMatch = &TransposeMatch[Aligned_Profile_Length*Index];

#if 0
    v_w[0] = lMatch[0] > 0 ? lMatch[0] : 0;
    if (lMatch[0] > Sc[0] ) Sc[0] = lMatch[0];

    for (size_t iprf=1; iprf<Profile_Length; ++iprf) {
      register int tmp   = v_r[iprf-1] + lMatch[iprf];
      tmp = tmp > 0 ? tmp : 0;
      v_w[iprf] = tmp;

      if (tmp > Sc[iprf]) Sc[iprf] = tmp;
    }
#else

    __m128i __V_R_0 = _mm_load_si128((__m128i*) &v_r[0]);
    __V_R_0         = _mm_slli_si128(__V_R_0, 4);

    iprf=0;
    goto Insert;

    Loop:
      __V_R_0                  = _mm_loadu_si128((__m128i*) &v_r[iprf-1]);

    Insert:
    ;

    __V_R_0                  = _mm_add_epi32(__V_R_0, *((__m128i*)&lMatch[iprf]));

    __m128i __V_R_1          = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 4]);
    __V_R_1                  = _mm_add_epi32(__V_R_1, *((__m128i*)&lMatch[iprf+4]));

    __m128i __V_R_2          = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 8]);
    __V_R_2                  = _mm_add_epi32(__V_R_2, *((__m128i*)&lMatch[iprf+8]));

    __m128i __V_R_3          = _mm_loadu_si128((__m128i*)&v_r[iprf-1 + 12]);
    __V_R_3                  = _mm_add_epi32(__V_R_3, *((__m128i*)&lMatch[iprf+12]));

    __m128i __SC_0           = _mm_load_si128((__m128i*)&Sc[iprf]);
    __V_R_0                  = _mm_max_epi32( __Zero, __V_R_0);

    __m128i __SC_1           = _mm_load_si128((__m128i*)&Sc[iprf + 4]);
    __V_R_1                  = _mm_max_epi32( __Zero, __V_R_1);

    __m128i __SC_2           = _mm_load_si128((__m128i*)&Sc[iprf + 8]);
    __V_R_2                  = _mm_max_epi32( __Zero, __V_R_2);

    __m128i __SC_3           = _mm_load_si128((__m128i*)&Sc[iprf + 12]);
    __V_R_3                  = _mm_max_epi32( __Zero, __V_R_3);

    _mm_store_si128((__m128i*)&v_w[iprf   ], __V_R_0);
    _mm_store_si128((__m128i*)&v_w[iprf+4 ], __V_R_1);
    _mm_store_si128((__m128i*)&v_w[iprf+8 ], __V_R_2);
    _mm_store_si128((__m128i*)&v_w[iprf+12], __V_R_3);

    __SC_0 = _mm_max_epi32(__SC_0, __V_R_0);
    __SC_1 = _mm_max_epi32(__SC_1, __V_R_1);
    __SC_2 = _mm_max_epi32(__SC_2, __V_R_2);
    __SC_3 = _mm_max_epi32(__SC_3, __V_R_3);

    _mm_store_si128((__m128i*)&Sc[iprf   ], __SC_0);
    _mm_store_si128((__m128i*)&Sc[iprf+4 ], __SC_1);
    _mm_store_si128((__m128i*)&Sc[iprf+8 ], __SC_2);
    _mm_store_si128((__m128i*)&Sc[iprf+12], __SC_3);
    iprf+= 16;

    if ( iprf<Profile_Length) goto Loop;

#endif

    // Swap pointers
    int * ptr = v_w;
    v_w = (int*) v_r;
    v_r = (const int*) ptr;

    // Update Score
    //Score += max;
  }

#if 1
  for (iprf=0; iprf<Profile_Length; ++iprf)
   Score += Sc[iprf];
#else
   // WARNING : THERE IS AN ERROR SOMEWHERE !!!
  iprf = 0;
  __m128i __s1 = _mm_setzero_si128();
  __m128i __s2 = _mm_setzero_si128();
  __m128i __s3 = _mm_setzero_si128();
  __m128i __s4 = _mm_setzero_si128();
  
  do {
    __m128i __sc1 = _mm_load_si128((__m128i*)&Sc[iprf   ]);
    __m128i __sc2 = _mm_load_si128((__m128i*)&Sc[iprf+ 4]);
    __m128i __sc3 = _mm_load_si128((__m128i*)&Sc[iprf+ 8]);
    __m128i __sc4 = _mm_load_si128((__m128i*)&Sc[iprf+12]);
    __s1 = _mm_add_epi32(__s1, __sc1);
    __s2 = _mm_add_epi32(__s2, __sc2);
    __s3 = _mm_add_epi32(__s3, __sc3);
    __s4 = _mm_add_epi32(__s4, __sc4);
    iprf +=16;
  } while (iprf < Profile_Length);

  __s1 = _mm_add_epi32(__s1, __s2);
  __s3 = _mm_add_epi32(__s3, __s4);
  __s1 = _mm_add_epi32(__s1, __s3);

  __asm__ __volatile__ ( "pshufd    $14, %1, %2 \n\t"
			 "paddd     %2, %1      \n\t"
			 "pshufd    $57, %1, %3 \n\t"
			 "paddd     %3, %1      \n\t"
			 "pextrd    $0, %1, %0  \n\t"
			 : "=r"(Score)
			 : "x"(__s1), "x"(__s2), "x"(__s3)
	  );
#endif

  return Score;

}

