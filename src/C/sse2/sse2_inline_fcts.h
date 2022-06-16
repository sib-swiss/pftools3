/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 30, 2011 sse2_inline_fcts.h
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/

#ifndef SSE2_INLINE_FCTS_H_
#define SSE2_INLINE_FCTS_H_
#include <emmintrin.h>
// DOES NOT ZERO EXTEND
//MOVSS __m128 _mm_load_ss(float * p)
//MOVSS void_mm_store_ss(float * p, __m128 a)
//MOVSS __m128 _mm_move_ss(__m128 a, __m128 b)
// DOES ZERO EXTEND, ALLOW from reg to xmm
//MOVD __m64 _mm_cvtsi32_si64 (int i )
//MOVD int _mm_cvtsi64_si32 ( __m64m )
//MOVD __m128i _mm_cvtsi32_si128 (int a)
//MOVD int _mm_cvtsi128_si32 ( __m128i a)

#ifndef __USE_32BIT_INTEGER__
extern __inline __m128i __ALWAYS_INLINE
LoadStoredIntegerVector(const ScoreTuple * const address)
{
    // Load short integer
    __m128i __A = _mm_loadl_epi64((__m128i*) &(address->vector));
    // Convert signed WORD into signed DWORD
    const __m128i __sign = _mm_cmpgt_epi16((__m128i) _mm_setzero_si128(), __A);
    // Interleave sign with data to produce a 128 bit (4 x DWORD)
    return _mm_unpacklo_epi16 (__A, __sign);
}
extern __inline __m128 __ALWAYS_INLINE
LoadAndConvertToFloatStoredIntegerVector(const ScoreTuple * const address)
{
  // Load short integer
  __m128i __A = _mm_loadl_epi64((__m128i*) &(address->vector));
  // Convert signed WORD into signed DWORD
  const __m128i __sign = _mm_cmpgt_epi16((__m128i) _mm_setzero_si128(), __A);
  // Interleave sign with data to produce a 128 bit (4 x DWORD)
  __A = _mm_unpacklo_epi16 (__A, __sign);
  // Convert the doublewords to floating point two at a time.
  return _mm_cvtepi32_ps (__A);
}

extern  __inline __m128 __ALWAYS_INLINE
LoadAndConvertToFloatAndAddStoredIntegerVector(const ScoreTuple * const address, const __m128 __B)
{
  // Load short integer
  __m128i __A = _mm_loadl_epi64((__m128i*) &(address->vector));
  // Convert signed WORD into signed DWORD
  const __m128i __sign = _mm_cmpgt_epi16((__m128i) _mm_setzero_si128(), __A);
  // Interleave sign with data to produce a 128 bit (4 x DWORD)
  __A = _mm_unpacklo_epi16 (__A, __sign);
  // Convert the doublewords to floating point two at a time.
  __m128 __Af = _mm_cvtepi32_ps(__A);
  // Convert the doublewords to floating point two at a time.
  return _mm_add_ps(__Af, __B);
}

#else /*************************************** __USE_32BIT_INTEGER__*******************************************/
extern __inline __m128i __ALWAYS_INLINE
LoadStoredIntegerVector(const ScoreTuple * const address)
{
  return _mm_load_si128(&(address->vector));
}

extern __inline __m128 __ALWAYS_INLINE
LoadAndConvertToFloatStoredIntegerVector(const ScoreTuple * const address)
{
  // Convert the doublewords to floating point two at a time.
  return _mm_cvtepi32_ps (address->vector);
}

extern  __inline __m128 __ALWAYS_INLINE
LoadAndConvertToFloatAndAddStoredIntegerVector(const ScoreTuple * const address, const __m128 __B)
{
  // Convert the doublewords to floating point two at a time.
  const __m128 __Af = _mm_cvtepi32_ps(address->vector);
  // Convert the doublewords to floating point two at a time.
  return _mm_add_ps(__Af, __B);
}

#endif /*************************************** __USE_32BIT_INTEGER__*******************************************/

extern __inline __m128 __ALWAYS_INLINE
_my_cvtpi16_epi32_add_cvt_ps (const short int * const address, const __m128i __B)
{
  // Load short integer
  __m128i __A = _mm_loadl_epi64((__m128i*) address);
  // Convert signed WORD into signed DWORD
  const __m128i __sign = _mm_cmpgt_epi16((__m128i) _mm_setzero_si128(), __A);
  // Interleave sign with data to produce a 128 bit (4 x DWORD)
  __A = _mm_unpacklo_epi16 (__A, __sign);
  // Add __B to __A
  __A = _mm_add_epi32(__A, __B);
  // Convert the doublewords to floating point two at a time.
  return _mm_cvtepi32_ps(__A);
}

extern __inline __m128 __ALWAYS_INLINE
_my_cvtdw_ss(const short int value)
{
  return (__m128) _mm_cvtsi32_ss(_mm_setzero_ps(),(int) value);
}

extern __inline __m128i __ALWAYS_INLINE
_my_cvtepi16_epi32(const __m128i __Value)
{
  // Convert signed WORD into signed DWORD
  const __m128i __sign = _mm_cmpgt_epi16((__m128i) _mm_setzero_si128(), __Value);
  // Interleave sign with data to produce a 128 bit (4 x DWORD)
  return  _mm_unpacklo_epi16 (__Value, __sign);
}

extern __inline __m128i __ALWAYS_INLINE
_my_blendv_epi8(const __m128i __A, const __m128i __B, const __m128i __Mask)
{
    __m128i __tmpA = _mm_andnot_si128(__Mask, __A);
    __m128i __tmpB = _mm_and_si128(__B, __Mask);
    return _mm_or_si128(__tmpA, __tmpB);
}

extern __inline __m128i __ALWAYS_INLINE
_my_insert_epi32_POS1(__m128i __A, const int B)
{
  __m128i __x = _mm_cvtsi32_si128(B);
  __m128i __y = _mm_unpacklo_epi32(__A, __x);
  __asm__ __volatile__ (" movsd %1, %0 " : "=x"(__A) : "x"(__y), "0"(__A) );
  return __A;
}

extern __inline __m128 __ALWAYS_INLINE
_my_insert_ss_ps_POS1(__m128 __A, const __m128 __B)
{
  __m128 __y = _mm_unpacklo_ps(__A, __B);
  __A        = (__m128) _mm_move_sd ((__m128d)__A, (__m128d) __y);
  return __A;
}

extern __inline __m128i __ALWAYS_INLINE
_my_insert_epi32_POS0(__m128i __A, const int B)
{
  __m128i __tmp = _mm_cvtsi32_si128(B);
  return _mm_or_si128(_mm_and_si128(__A, _mm_set_epi32(-1,-1,-1,0)), __tmp);
}
#endif
