/*******************************************************
                        PFTOOLS
 *******************************************************
  Apr 29, 2014 sse41_inline_fcts.h
 *******************************************************
 (C) 2011-14 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#ifndef SSE41_INLINE_FCTS_H_
#define SSE41_INLINE_FCTS_H_
#ifdef __SSE4_1__
#ifndef __USE_32BIT_INTEGER__
extern __inline __m128i __ALWAYS_INLINE
LoadStoredIntegerVector(const ScoreTuple * const address)
{
      return _mm_cvtepi16_epi32(_mm_loadl_epi64((__m128i*) &(address->vector)));
}
#else  /*__USE_32BIT_INTEGER__*/
extern __inline __m128i __ALWAYS_INLINE
LoadStoredIntegerVector(const ScoreTuple * const address)
{
      return _mm_load_si128(&(address->vector));
}
#endif /*__USE_32BIT_INTEGER__*/
#else  /*__SSE4_1__*/
#error "sse41_inline_fcts.h should not be included in non sse 4.1 files"
#endif /*__SSE4_1__*/

#endif