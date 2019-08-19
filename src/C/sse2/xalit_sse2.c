/*******************************************************
                        PFTOOLS
 *******************************************************
  Sep 26, 2011 xalit_sse2.c
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

static inline unsigned char Ncode(__m128i index)
{
#ifdef XALIT_DEBUG
    static int count = 0;
    union lScores tmp = { xmm:
        index
    };
    printf("NCODE  %10i %10i %10i %10i\n" , tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]
          );
#endif
    const unsigned int JD = _mm_cvtsi128_si32(index);
    __m128i lindex = _mm_srli_si128(index, 8);
    const unsigned int JM = _mm_cvtsi128_si32(lindex);
    lindex = _mm_srli_si128(lindex, 4);
    const unsigned int JI = _mm_cvtsi128_si32(lindex);
    register const unsigned char ctmp1 = ((unsigned char) JM ) << 4;
    register const unsigned char ctmp2 = ((unsigned char) JI ) << 2;
    register const unsigned char ctmp3 = ((unsigned char) JD ) & 3;
    unsigned char res = (unsigned char) (ctmp1 | ctmp2 | ctmp3);
#ifdef XALIT_DEBUG
    printf("XALIT CPMA %10i %8u %8u %8u\n", count++, JI, JM, JD);
#endif
    return res;
}

static inline void Dcode(const unsigned char Data, unsigned int * const JM,
                         unsigned int * const JI, unsigned int * const JD)
{
    register const unsigned char ctmp = Data >> 2;
    *JD = (unsigned int) (Data & 3);
    *JI = (unsigned int) (ctmp & 3);
    *JM = (unsigned int) ((Data >> 4) & 3);
}

static __m128i MaxP(__m128i * const Data, const __m128i Offset,
                    const TransitionScores * const restrict Transition)
{
    __m128i __KM           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(MATCH,MATCH,MATCH,MATCH));
    __m128i __TransitionsM = LoadStoredIntegerVector(&(Transition->From[MATCH]));
#ifdef XALIT_DEBUG_PREVIOUS
union lScores tmp = { xmm:
        __TransitionsM
    };
    printf("XALIT M %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsM         = _mm_add_epi32(__TransitionsM, __KM);

    __m128i __KI           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(INSERTION,INSERTION,INSERTION,INSERTION));
    __m128i __TransitionsI = LoadStoredIntegerVector(&(Transition->From[INSERTION]));
#ifdef XALIT_DEBUG_PREVIOUS
    tmp.xmm = __TransitionsI;
    printf("XALIT I %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsI         = _mm_add_epi32(__TransitionsI, __KI);

    __m128i __KD           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(DELETION,DELETION,DELETION,DELETION));
    __m128i __TransitionsD = LoadStoredIntegerVector(&(Transition->From[DELETION]));
#ifdef XALIT_DEBUG_PREVIOUS
    tmp.xmm = __TransitionsD;
    printf("XALIT D %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsD         = _mm_add_epi32(__TransitionsD, __KD);

    __m128i __TransitionsX = LoadStoredIntegerVector(&(Transition->From[EXTRA]));

#ifdef XALIT_DEBUG
#ifndef XALIT_DEBUG_PREVIOUS
    union lScores tmp;
#endif
    tmp.xmm = Offset;
    printf("XALIT KMID %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsM;
    printf("XALIT M %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsI;
    printf("XALIT I %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsD;
    printf("XALIT D %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsX;
    printf("XALIT X %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif

    __m128i __Index        = (__m128i) _mm_setzero_si128();
		const __m128i __One    = _mm_set1_epi32(1);
	  const __m128i __Three  = _mm_set1_epi32(3);
    const __m128i __Two    = _mm_add_epi32(__One, __One);

    __m128i __Mask         = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
    __Index                = _my_blendv_epi8(__Index, __Three, __Mask);

    __Mask                 = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
    __Index                = _my_blendv_epi8(__Index, __Two, __Mask);

    __Mask                 = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
    __Index                = _my_blendv_epi8(__Index, __One, __Mask);

    *Data                  = __TransitionsX;
#ifdef XALIT_DEBUG
    tmp.xmm = __TransitionsX;
    printf("XALIT MAXVAL %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);

    tmp.xmm = __Index;
    printf("XALIT INDEX %10i %10i %10i %10i\n",  tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    return __Index;
}

static __m128i MaxP_border(__m128i * const Data, const __m128i Offset,
                           const ScoreTuple * const Border, const TransitionScores * const restrict Transition)
{
    __m128i __KM           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(MATCH,MATCH,MATCH,MATCH));
    __m128i __TransitionsM = LoadStoredIntegerVector(&(Transition->From[MATCH]));
//    __TransitionsM         = _my_insert_epi32_POS1(__TransitionsM, (int) Border->To[MATCH]);
#ifdef XALIT_DEBUG_PREVIOUS
union lScores tmp = { xmm:
        __TransitionsM
    };
    printf("XALIT M BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsM         = _mm_add_epi32(__TransitionsM, __KM);

    __m128i __KI           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(INSERTION,INSERTION,INSERTION,INSERTION));
    __m128i __TransitionsI = LoadStoredIntegerVector(&(Transition->From[INSERTION]));
//    __TransitionsI         = _my_insert_epi32_POS1(__TransitionsI, (int) Border->To[INSERTION]);
#ifdef XALIT_DEBUG_PREVIOUS
    tmp.xmm = __TransitionsI;
    printf("XALIT I BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsI         = _mm_add_epi32(__TransitionsI, __KI);

    __m128i __KD           = _mm_shuffle_epi32(Offset, _MM_SHUFFLE(DELETION,DELETION,DELETION,DELETION));
    __m128i __TransitionsD = LoadStoredIntegerVector(&(Transition->From[DELETION]));
//    __TransitionsD         = _my_insert_epi32_POS1(__TransitionsD, (int) Border->To[DELETION]);
#ifdef XALIT_DEBUG_PREVIOUS
    tmp.xmm = __TransitionsD;
    printf("XALIT D BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    __TransitionsD         = _mm_add_epi32(__TransitionsD, __KD);

//    __m128i __TransitionsX = _mm_loadl_epi64((__m128i*) &(Transition->From[EXTRA]));
    __m128i __TransitionsX = LoadStoredIntegerVector(Border);

#ifdef XALIT_DEBUG
#ifndef XALIT_DEBUG_PREVIOUS
    union lScores tmp;
#endif
    tmp.xmm = Offset;
    printf("XALIT KMID BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsM;
    printf("XALIT M BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsI;
    printf("XALIT I BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsD;
    printf("XALIT D BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __TransitionsX;
    printf("XALIT X BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif

    __m128i __Index        = (__m128i) _mm_setzero_si128();
		const __m128i __One    = _mm_set1_epi32(1);
	  const __m128i __Three  = _mm_set1_epi32(3);
    const __m128i __Two    = _mm_add_epi32(__One, __One);

    __m128i __Mask         = _mm_cmpgt_epi32(__TransitionsD, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsD, __Mask);
    __Index                = _my_blendv_epi8(__Index, __Three, __Mask);

    __Mask                 = _mm_cmpgt_epi32(__TransitionsM, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsM, __Mask);
    __Index                = _my_blendv_epi8(__Index, __Two, __Mask);

    __Mask                 = _mm_cmpgt_epi32(__TransitionsI, __TransitionsX);
    __TransitionsX         = _my_blendv_epi8(__TransitionsX, __TransitionsI, __Mask);
    __Index                = _my_blendv_epi8(__Index, __One, __Mask);

    *Data                  = __TransitionsX;

#ifdef XALIT_DEBUG
    tmp.xmm = __TransitionsX;
    printf("XALIT MAXVAL BORDER %10i %10i %10i %10i\n", tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
    tmp.xmm = __Index;
    printf("XALIT INDEX BORDER %10i %10i %10i %10i\n",  tmp.Element[MATCH],
           tmp.Element[INSERTION],
           tmp.Element[DELETION],
           tmp.Element[DUMMY]);
#endif
    return __Index;
}

// This functions fills in IPMB, IPME and CALI that are then required !
int xalit_sse2(const struct Profile * const restrict prf, const size_t N1, const size_t N2, const size_t bseq,
               const PFSequence * const restrict PFSeq, char * const restrict CALI, union lScores * const restrict iop,
               struct Alignment * const restrict alignment, const _Bool * const restrict Lock)
{
    int IPM[2];
    unsigned int K3 = (unsigned int) prf->Length;
    unsigned int JS = 0;
    const size_t lseq = PFSeq->Length;
    const unsigned char * const restrict Sequence = PFSeq->ProfileIndex;
    //////////////////////////////////////////////////////////////////////////////////////////////
    // Prologue
    //////////////////////////////////////////////////////////////////////////////////////////////
    const unsigned int SequenceBegin = alignment->Region.Sequence.Begin;
    const unsigned int SequenceEnd   = alignment->Region.Sequence.End;
    const int CutOff                 = alignment->Score;

#ifdef XALIT_DEBUG
    fprintf(stdout,"XALIT ALIGN %i %i %i | %u %u %i %lu\n",
            alignment->Score,
            alignment->Region.Sequence.Begin, alignment->Region.Sequence.End, SequenceBegin, SequenceEnd, CutOff, prf->Length );
#endif

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Allocate memory
    const size_t Memory = (prf->Length+1) * ((size_t)(SequenceEnd-SequenceBegin+2));
    const unsigned int prfLength = prf->Length;

#ifdef XALIT_DEBUG
    fprintf(stdout, "XALIT Memory is %lu\n", Memory);
#endif
    unsigned char * const restrict cpma = (unsigned char*) malloc(Memory*sizeof(unsigned char));
    if (cpma == NULL) {
        fputs("Unable to allocate memory to store coding/decoding path\n", stderr);
        return -1;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Loop through the path
    //////////////////////////////////////////////////////////////////////////////////////////////
    unsigned int K1 = 0;
    register const TransitionScores * restrict pTransitions = prf->Scores.Insertion.Transitions;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Beginning of sequence
    register __m128i __KMID = _mm_set1_epi32((int) NLOW);
    register __m128i __Index;
    register const StoredIntegerFormat * pMatch = prf->Scores.Match.Alphabet;
    register const size_t AlignStep   = prf->Scores.Match.AlignStep;

    if (SequenceBegin == 1) {
        const register ScoreTuple * const restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
        __Index = MaxP_border(&iop[0].xmm, __KMID, &FirstSequenceProtein[0], pTransitions);

        // Store temporary data packed as much as possible into one byte
        cpma[K1] = Ncode(__Index);

        // Loop through the rest of the profile
        for (unsigned int iprf=1; iprf<=prfLength; ++iprf) {
            pTransitions++;
            const int KD = iop[iprf-1].Element[DELETION] + (int) pMatch[_D];
            pMatch += AlignStep;
            __KMID  = _my_insert_epi32_POS0(__KMID, KD);
            __Index = MaxP_border(&(iop[iprf].xmm), __KMID, &FirstSequenceProtein[iprf], pTransitions);
            cpma[++K1] = Ncode(__Index);
        }
    } else {
        __Index = MaxP(&(iop[0].xmm), __KMID, pTransitions);

        // Store temporary data packed as much as possible into one byte
        cpma[K1] = Ncode(__Index);

        // Loop through the rest of the profile
        for (unsigned int iprf=1; iprf<=prfLength; ++iprf) {
            pTransitions++;
            const int KD = iop[iprf-1].Element[DELETION] + (int) pMatch[_D];
            pMatch += AlignStep;
            __KMID = /*_mm_set_epi32(NLOW, NLOW, NLOW, KD);*/ _my_insert_epi32_POS0(__KMID, KD);
            __Index = MaxP(&(iop[iprf].xmm), __KMID, pTransitions);
            cpma[++K1] = Ncode(__Index);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Loop through the internal sequence indices
    for (unsigned int iseq = SequenceBegin; iseq<SequenceEnd; ++iseq) {
        // Protected region
        __KMID = _mm_set1_epi32(NLOW);
        if (Lock[iseq]) {
            iop[N1-1].Element[MATCH] = NLOW;
            for ( size_t i=N1; i<N2; ++i) {
                // _mm_storel_pi((__m64*) &iop[i].xmm, (__m128) __KMID);
                iop[i].xmm = __KMID;
            }
        }

        const unsigned int SequenceIndex = Sequence[iseq-1]; // Fortran one based index
#ifdef XALIT_DEBUG
        fprintf(stdout, "XALIT SEQUENCE %lu LETTER %u\n",iseq, SequenceIndex);
#endif
        // Pointers to Score data
        pTransitions = prf->Scores.Insertion.Transitions;
        register const StoredIntegerFormat * restrict pInsertion = prf->Scores.Insertion.Alphabet;
        pMatch = prf->Scores.Match.Alphabet;

        register int KOPM = iop[0].Element[MATCH];
        register int KI   = iop[0].Element[INSERTION] + (int) pInsertion[SequenceIndex];
        __KMID = _mm_set_epi32(KI,NLOW,NLOW,NLOW);
        __Index    = MaxP(&iop[0].xmm, __KMID, pTransitions);
        cpma[++K1] = Ncode(__Index);

        // Loop through the rest of the profile
        for (unsigned int iprf=1; iprf<=prfLength; ++iprf) {
            pTransitions++;
            pInsertion   += AlignStep;
            const int KM  = KOPM + (int) pMatch[SequenceIndex];
//         __KMID        = _mm_insert_epi32(__KMID, KM, MATCH);
            KI            = iop[iprf].Element[INSERTION] + (int) pInsertion[SequenceIndex];
//         __KMID        = _mm_insert_epi32(__KMID, KI, INSERTION);
            const int KD  = iop[iprf-1].Element[DELETION] + (int) pMatch[_D];
//         __KMID        = _mm_insert_epi32(__KMID, KD, DELETION);
            __KMID        = _mm_set_epi32(KI,KM,NLOW,KD);
            pMatch       += AlignStep;
            KOPM          = iop[iprf].Element[MATCH];
            __Index       = MaxP(&(iop[iprf].xmm), __KMID, pTransitions);
            cpma[++K1]    = Ncode(__Index);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Last sequence index

    // Protected region
    __KMID = _mm_set1_epi32(NLOW);
    if (Lock[SequenceEnd]) {
        iop[N1-1].Element[MATCH] = NLOW;
        for ( size_t i=N1; i<N2; ++i) {
            //_mm_storel_pi((__m64*) &iop[i].xmm, (__m128) __KMID);
            iop[i].xmm = __KMID;
        }
    }

    const unsigned int SequenceIndex = Sequence[SequenceEnd-1]; // Fortran one based index

    // Pointers to Score data
    pTransitions = prf->Scores.Insertion.Transitions;
    register const StoredIntegerFormat * restrict pInsertion = prf->Scores.Insertion.Alphabet;
    pMatch = prf->Scores.Match.Alphabet;

    register int KOPM = iop[0].Element[MATCH];
    register int KI   = iop[0].Element[INSERTION] + (int) pInsertion[SequenceIndex];
    __KMID = _mm_set_epi32(KI,NLOW,NLOW,NLOW);

    // Last Sequence element or not
    if (SequenceEnd == lseq) {
        const ScoreTuple * const restrict LastSequenceProtein = prf->Scores.Insertion.LastSequenceProtein;
        __Index    = MaxP_border(&iop[0].xmm, __KMID, &LastSequenceProtein[0], pTransitions);
        cpma[++K1] = Ncode(__Index);

        // Check if first profile insertion is enough to reach cutoff
#ifdef XALIT_DEBUG
        fprintf(stdout,"XALIT BORDER LAST SEQ %10i %10i %10i\n",0, KI + (int) pTransitions->From[INSERTION].To[EXTRA], KOPM);
#endif
        if ((KI + (int) pTransitions->From[INSERTION].To[EXTRA]) >= CutOff) {
            K3 = 0;
            JS = 1;
        } else {
            // Loop through the rest of the profile
            for (unsigned int iprf=1; iprf<=prfLength; ++iprf) {
                pTransitions++;
                pInsertion   += AlignStep;
                const int KM  = KOPM + (int) pMatch[SequenceIndex];
//            __KMID        = _mm_insert_epi32(__KMID, KM, MATCH);
                KI            = iop[iprf].Element[INSERTION] + (int) pInsertion[SequenceIndex];
//            __KMID        = _mm_insert_epi32(__KMID, KM, INSERTION);
                const int KD  = iop[iprf-1].Element[DELETION] + (int) pMatch[_D];
//            __KMID        = _mm_insert_epi32(__KMID, KM, DELETION);
                __KMID        = _mm_set_epi32(KI,KM,NLOW,KD);
                pMatch       += AlignStep;
                KOPM          = iop[iprf].Element[MATCH];
//            __Index       = MaxP_border(&iop[iprf].xmm, __KMID, &LastSequenceProtein[iprf], pTransitions);
                __Index       = MaxP(&iop[iprf].xmm, __KMID, pTransitions);
                cpma[++K1]    = Ncode(__Index);
#ifdef XALIT_DEBUG
                fprintf(stdout,"XALIT BORDER LAST SEQ %10u %10i %10i %10i\n",iprf,
                        KI + (int) LastSequenceProtein[iprf].To[INSERTION],
                        KM + (int) LastSequenceProtein[iprf].To[MATCH],
                        KD + (int) LastSequenceProtein[iprf].To[DELETION]
                       );
#endif

                if ((KI + (int) LastSequenceProtein[iprf].To[INSERTION]) >= CutOff) {
                    JS = 1;
                    K3 = iprf;
                    break;
                } else if ((KM + (int) LastSequenceProtein[iprf].To[MATCH]) >= CutOff) {
                    JS = 2;
                    K3 = iprf;
                    break;
                } else if ((KD + (int) LastSequenceProtein[iprf].To[DELETION]) >= CutOff) {
                    JS = 3;
                    K3 = iprf;
                    break;
                }
            }
        }

    } else {
        __Index    = MaxP(&iop[0].xmm, __KMID, pTransitions);
        cpma[++K1] = Ncode(__Index);

        // Check if first profile insertion is enough to reach cutoff
#ifdef XALIT_DEBUG
        fprintf(stdout,"XALIT LAST SEQ %10i %10i %10i\n",0, KI + (int) pTransitions->From[INSERTION].To[EXTRA],KOPM);
#endif
        if ((KI + (int) pTransitions->From[INSERTION].To[EXTRA]) >= CutOff) {
            K3 = 0;
            JS = 1;
        } else {
            // Loop through the rest of the profile
            for (unsigned int iprf=1; iprf<=prfLength; ++iprf) {
                pTransitions++;
                pInsertion   += AlignStep;
                const int KM  = KOPM + (int) pMatch[SequenceIndex];
//            __KMID        = _mm_insert_epi32(__KMID, KM, MATCH);
                KI            = iop[iprf].Element[INSERTION] + (int) pInsertion[SequenceIndex];
//            __KMID        = _mm_insert_epi32(__KMID, KM, INSERTION);
                const int KD  = iop[iprf-1].Element[DELETION] + (int) pMatch[_D];
//            __KMID        = _mm_insert_epi32(__KMID, KM, DELETION);
                __KMID        = _mm_set_epi32(KI,KM,NLOW,KD);
                pMatch       += AlignStep;
                KOPM          = iop[iprf].Element[MATCH];
                __Index       = MaxP(&iop[iprf].xmm, __KMID, pTransitions);
                cpma[++K1]    = Ncode(__Index);
#ifdef XALIT_DEBUG
                fprintf(stdout,"XALIT LAST SEQ %10u %10i %10i %10i\n",iprf,
                        KI + (int) pTransitions->From[INSERTION].To[EXTRA],
                        KM + (int) pTransitions->From[MATCH].To[EXTRA],
                        KD + (int) pTransitions->From[DELETION].To[EXTRA]
                       );
#endif
                if ((KI + (int) pTransitions->From[INSERTION].To[EXTRA]) >= CutOff) {
                    JS = 1;
                    K3 = iprf;
                    break;
                } else if ((KM + (int) pTransitions->From[MATCH].To[EXTRA]) >= CutOff) {
                    JS = 2;
                    K3 = iprf;
                    break;
                } else if ((KD + (int) pTransitions->From[DELETION].To[EXTRA]) >= CutOff) {
                    JS = 3;
                    K3 = iprf;
                    break;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Epilogue
    //////////////////////////////////////////////////////////////////////////////////////////////
    //Epilogue:

    int K2          = SequenceEnd-1;
    IPM[1]          = K3 - prf->Length - 1;
    unsigned int J0 = 0;
    int J1          = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Fill in unused profile indices with character '-'
    for (unsigned int iprf = prfLength; iprf>K3; iprf--) {
        CALI[++J1] = (unsigned char) '-';
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    // Continue backward analysis of profile up to JS = 0
    const register char * const restrict CABC = prf->CABC;
    do {
#ifdef XALIT_DEBUG
        printf("XALIT BACKTRACE %10i %10i %8u %8u %10i %u %u %c %c\n", K1, K2, K3, JS, J1+1, Sequence[K2], CABC[Sequence[K2]], CABC[Sequence[K2]], PFSeq->OriginalSequence[K2]);
#endif
        // Beware that ordering is important here
        if ( JS == 1 ) {
            ///////////////////////////////////////////////
            // Insertion treated first
            CALI[++J1]  = (unsigned char) 32 + PFSeq->OriginalSequence[K2];
            K1         -= (int) prfLength + 1;
            --K2;
        } else if ( JS == 2 ) {
            ///////////////////////////////////////////////
            // Match treated second
            CALI[++J1]  = PFSeq->OriginalSequence[K2];
            K1         -= (int) prfLength + 2;
            --K2;
            --K3;
        } else if (JS == 3 ) {
            ///////////////////////////////////////////////
            // Deletion treated last
            CALI[++J1] = '-';
            --K1;
            --K3;
        }

        // Decode next traceback
        unsigned int JM, JI, JD;
        Dcode(cpma[K1], &JM, &JI, &JD);

        J0 = JS;
        // Beware that ordering is important here
        if ( JS == 1 ) {
            ///////////////////////////////////////////////
            // Insertion treated first
            JS = JI;
        } else if ( JS == 2) {
            ///////////////////////////////////////////////
            // Match treated second
            JS = JM;
        } else if (JS == 3 ) {
            ///////////////////////////////////////////////
            // Deletion treated last
            JS = JD;
        }
    } while (JS > 0);

    // Free no more needed memory
    free(cpma);

    K3 = (K3 >= prfLength) ? 0 : K3;

    for (int iprf=K3; iprf>=1; iprf--) {
        CALI[++J1] = '-';
    }

#ifdef XALIT_DEBUG
    fputs("\nXALIT SEQUENCE ", stderr);
    for (int i=1; i<=J1; ++i) fputc(CALI[i],stderr);
    fputs("\n", stderr);
    fflush(stderr);
#endif

    const int LALI = J1;
    J1 = (J1 + 1)/2;
    for (int iprf = LALI/2+1; iprf<=LALI; ++iprf) {
        register const char c = CALI[iprf];
        CALI[iprf] = CALI[J1];
        CALI[J1]   = c;
        --J1;
    }

    CALI[LALI+1] = '\0';

    IPM[0] = K3 + 1;

    alignment->IPMB = IPM[0];
    alignment->IPME = IPM[1];

    return 0;
}
