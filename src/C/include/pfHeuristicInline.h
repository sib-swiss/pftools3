/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfHeuristicInline.h
 ***************************************************************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#include <stdlib.h>
#ifdef BUILD_LIBRARY
#define PF_STATIC_INLINE(type) PFEXPORT type
#define PF_EXTERN_INLINE(type) PFEXPORT type
#else
#define PF_STATIC_INLINE(type) static inline type __ALWAYS_INLINE
#define PF_EXTERN_INLINE(type) extern inline type __ALWAYS_INLINE
#endif


#ifndef _SPECIAL_HEURISTIC_
PF_STATIC_INLINE(const int *) TransposeAndConvertMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
									 const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache line
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  int * restrict const TIMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(int),64);

  if (TIMatch == NULL) return TIMatch;
  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));
  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha <Alphabet_Length; ++alpha) {
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (int) lMatch[alpha];
    }
    lMatch += step;
  }
  return TIMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertMatchMatrixGivenMemory(int * const restrict TIMatch, const union Scores * const Matrices,
							     const size_t Alphabet_Length, const size_t Profile_Length,
							     const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));
  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (int) lMatch[alpha];
    }
    lMatch += step;
  }
}

PF_STATIC_INLINE(float * ) TransposeAndConvertToFloatMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
                                                            const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache line boundary
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  float * restrict const TFMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(float),64);

  if (TFMatch == NULL) return TFMatch;
  memset(TFMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(float));

  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TFMatch[alpha*Aligned_Profile_Length+iprf] = (float) lMatch[alpha];
    }
    lMatch += step;
  }
  return TFMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertToFloatMatchMatrixGivenMemory(float * const restrict TIMatch, const union Scores * const Matrices,
							            const size_t Alphabet_Length, const size_t Profile_Length,
							            const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));

  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) {
    for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (float) lMatch[alpha];
    }
    lMatch += step;
  }
}
#else
/*
 * These are special heuristic functions in the sense the limit negative effect of wrong match
 * replacing them where benefical by either an insertion or deletion score.
 *
 */

#define CHANGE(type) \
  const register StoredIntegerFormat * restrict lMatch = Matrices->Match.Alphabet;\
  const TransitionScores * restrict InsertionLine = Matrices->Insertion.Transitions;\
  \
  /* First line stay identical */\
  for (size_t alpha=0; alpha <Alphabet_Length; ++alpha) {\
    TIMatch[alpha*Aligned_Profile_Length] = (type) lMatch[alpha];\
  }\
  lMatch += step;\
  /* from 1 to n-1 we use special treatment */\
  for (size_t iprf = 1; iprf<Profile_Length-1; ++iprf) {\
    const StoredIntegerFormat MDDM    = InsertionLine[iprf].From[MATCH].To[DELETION]+InsertionLine[iprf+1].From[DELETION].To[MATCH];\
    const StoredIntegerFormat MIIM    = InsertionLine[iprf].From[MATCH].To[INSERTION]+InsertionLine[iprf+1].From[INSERTION].To[MATCH];\
    const StoredIntegerFormat Minimum = (MDDM>MIIM) ? MDDM : MIIM;\
    \
    for (size_t alpha=0; alpha <Alphabet_Length; ++alpha) {\
      const StoredIntegerFormat value = lMatch[alpha] < Minimum ? Minimum : lMatch[alpha];\
      TIMatch[alpha*Aligned_Profile_Length+iprf] = (type) value;\
    }\
    lMatch += step;\
  }\
  /* Last line */\
  for (size_t alpha=0; alpha <Alphabet_Length; ++alpha) {\
    TIMatch[alpha*Aligned_Profile_Length+Profile_Length-1] = (type) lMatch[alpha];\
  }

PF_EXTERN_INLINE(const int *) TransposeAndConvertMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
                                                         const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache line
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  int * restrict const TIMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(int),64);

  if (TIMatch == NULL) return TIMatch;
  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));

  CHANGE(int);

  return TIMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertMatchMatrixGivenMemory(int * const restrict TIMatch, const union Scores * const Matrices,
							     const size_t Alphabet_Length, const size_t Profile_Length,
							     const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(int));
  CHANGE(int);
}

PF_EXTERN_INLINE(float *) TransposeAndConvertToFloatMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
							    const size_t Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;
  // Profile size rounded to cache linerices
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;
  float * restrict const TIMatch = _mm_malloc(Aligned_Profile_Length*Alphabet_Length*sizeof(float),64);

  if (TIMatch == NULL) return TIMatch;
  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(float));

  CHANGE(float);

  return TIMatch;
}

PF_EXTERN_INLINE(void) TransposeAndConvertToFloatMatchMatrixGivenMemory(float * const restrict TIMatch, const union Scores * const Matrices,
							            const size_t Alphabet_Length, const size_t Profile_Length,
							            const size_t Aligned_Profile_Length)
{
  const size_t step = Matrices->Match.AlignStep;

  memset(TIMatch, 0, Aligned_Profile_Length*Alphabet_Length*sizeof(float));

  CHANGE(float);
}
#endif

#if 0
extern inline void __ALWAYS_INLINE WeightMatchMatrix(float * const restrict TIMatch, const struct Profile * prf,
				     const size_t Alphabet_Length, const size_t Profile_Length)
{
  float DatabaseWeightedScores[] = {
    0.0783f, 0.0000f, 0.0236f, 0.0496f, 0.0581f, 0.0463f, 0.0749f, 0.0248f,
    0.0621f, 0.0000f, 0.0602f, 0.0939f, 0.0251f, 0.0421f, 0.0000f, 0.0426f,
    0.0360f, 0.0524f, 0.0612f, 0.0523f, 0.0000f, 0.0700f, 0.0132f, 0.0000f,
    0.0339f, 0.0000f
  };

  /*
   * Allocate on the stack a Database scaling factor for each letter in profile alphabet
   * alphabet + 1 (+1 is for the zero position being the unknown)
   */
  float * const restrict DatabaseFactor = (float*) alloca((Alphabet_Length)*sizeof(float));
  for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
    const unsigned char letter = (unsigned char) prf->CABC[alpha];
    register size_t index = (size_t) ( (letter >= (unsigned char) 'a' ) ? letter - ((unsigned char) 'a') : letter - ((unsigned char) 'A') );

//     fprintf(stderr,"Alphabet %lu : %1c\t%lu\n",alpha,prf->CABC[alpha],index);
    if ( index < 26 ) {
      DatabaseFactor[alpha] = DatabaseWeightedScores[index];
    } else {
      fprintf(stderr,"Database weight scores does not contain all alphabet letters, %1c is missing\n", prf->CABC[alpha]);
      exit(1);
    }

  }
  /* Alignment of profile pn cache line (64 bytes) */
  const size_t Aligned_Profile_Length = (Profile_Length+1 + 15) & ~15;

  /* We allocate here on the stack, so beware */
  //float * const restrict AlphabetSum = (float*) alloca(Alphabet_Length*sizeof(float));
  float * const restrict ProfileSum = (float*) alloca(Aligned_Profile_Length*sizeof(float));
  memset(ProfileSum, 0, Aligned_Profile_Length*sizeof(float));



  /* Compute overall sum and sum per letter */
  register float Sum = 0.0f;
  for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
    register float tSum = 0.0f;
    //#pragma unroll(4)
    //for (size_t iprf = 0; iprf<Profile_Length; ++iprf) tSum += TIMatch[alpha*Aligned_Profile_Length+iprf];
    //AlphabetSum[alpha] = tSum;
#pragma unroll(4)
    for (size_t iprf=0; iprf<Profile_Length; ++iprf) {
	ProfileSum[iprf] += TIMatch[alpha*Aligned_Profile_Length+iprf];
	tSum             += TIMatch[alpha*Aligned_Profile_Length+iprf];
    }
    Sum += tSum;
  }

  /* Normalize factor and include database frequency as well */
  const register float Norm = 1.0f/Sum;
//   for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
//     AlphabetSum[alpha] *= Norm*DatabaseFactor[alpha];
//   }
  for (size_t iprf = 0; iprf<Profile_Length; ++iprf) ProfileSum[iprf] *= Norm;

  /* Apply the scaling to the Match matrix */
  for (size_t alpha=0; alpha<Alphabet_Length; ++alpha) {
    //const register float scale = AlphabetSum[alpha];
    //#pragma unroll(4)
    //for (size_t iprf = 0; iprf<Profile_Length; ++iprf) TIMatch[alpha*Aligned_Profile_Length+iprf] *= scale;
    const register float scale = DatabaseFactor[alpha];
    for (size_t iprf = 0; iprf<Profile_Length; ++iprf) TIMatch[alpha*Aligned_Profile_Length+iprf] *= (1.0f - ProfileSum[iprf]*scale);
  }
}
#endif