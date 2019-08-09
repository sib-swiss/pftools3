/*******************************************************
                        PFTOOLS
 *******************************************************
  Jan 16, 2011 Normalization.c
 *******************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdint.h>
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include "../include/pfProfile.h"

//   Subroutine CPAve(IMPP,IDMP,LPRF,CABC,NABC,PAVE)
// 
//   Integer           IMPP(0:27,0:IDMP)
//   Character         CABC(0:26)
//   Real              PAVE(0:26)
// 
//   Do  I1=0,NABC
//       PAVE(I1)=0
//   End Do
//   Do  I1=1,LPRF
//       Do  I2=0,NABC
// 	       PAVE(I2)=PAVE(I2)+IMPP(I2,I1) 
//       End do  
//   End Do 
// 
//   Return
//   End

int InitAverage(struct Profile * const restrict prf) //const union Scores * Matrices, const size_t prfLength, const size_t AlphabetLength, SAverage * const Average)
{
 const size_t AlignedStep          = prf->Scores.Match.AlignStep;
 const StoredIntegerFormat * restrict lMatch = prf->Scores.Match.Alphabet;
 const size_t AlphabetLength       = prf->Alphabet_Length;
 SAverage * const restrict Average = &prf->Average;
 
 /* Allocate memory */
 float * const restrict Weights = _mm_malloc(AlphabetLength*sizeof(float),16);
 if (Weights == 0) {
    Average->Weights = 0;
    Average->size    = 0;
    return -1;
 } else {
    memset(Weights, 0, AlphabetLength*sizeof(float));
    const size_t prfLength = prf->Length;
    for (size_t iprf=0; iprf<prfLength; ++iprf) {
      for (size_t alpha=0; alpha<AlphabetLength; ++alpha) Weights[alpha] += (float) lMatch[alpha];
      lMatch += AlignedStep;
    }
    Average->Weights = Weights;
    Average->size    = AlphabetLength;
 }
 return 0;
}

// Subroutine CFAve(ISEQ,IDMS,LSEQ,CABC,NABC,FAVE)
// 
//   Integer*2         ISEQ(IDMS)
//   Character         CABC(0:26)
//   Real              FAVE(0:26)
// 
//   Do  I1=0,NABC
//       FAVE(I1)=0
//   End Do
//   Do  I1=1,LSEQ
//       J1=ISEQ(I1)
//       FAVE(J1)=FAVE(J1)+1
//   End Do 
//   Do  I1=0,NABC
//       FAVE(I1)=FAVE(I1)/LSEQ
//   End Do
// 
//   Return
//   End


float ComputeAverageFrequencies(const PFSequence * const restrict Sequence, const SAverage * const restrict Average)
{
    const size_t AlphabetLength = Average->size;
    const uintptr_t stack = (uintptr_t) alloca(AlphabetLength*sizeof(float)+15);
    float * const restrict data = (float*) ( stack & ~15 );
    memset(data, 0, AlphabetLength*sizeof(float));
    
    const size_t SeqLength = Sequence->Length;
    for (size_t iseq=0; iseq<SeqLength; ++iseq) {
      const size_t index = (size_t) Sequence->ProfileIndex[iseq];
      data[index] += 1.0f;
    }
    const float scale = 1.0f/( (float)SeqLength );
    float value = 0.0f;
    const float * const restrict Weights = Average->Weights;
    for (size_t alpha=0; alpha<AlphabetLength; ++alpha) value += scale*Weights[alpha]*data[alpha];
    return value;
}

// Subroutine NtoR(R,N,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
// 
//         Real              RNOP(KNPM,MAXN)
// 
//         If     (IFUN.EQ.1) then 
//            X=(R-RNOP(1,INOR)) / RNOP(2,INOR) 
//         Else if(IFUN.EQ.2) then 
//            X=( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
//      *      *( RNOP(5,INOR) * R + RNOP(4,INOR) )
//         Else if(IFUN.EQ.3) then 
//            X=( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
//      *      *( RNOP(5,INOR) * R + RNOP(4,INOR) ) + RAVE
//         End if 
//            N=INT(X)
//            If(Real(N).LT.X) N=N+1
//         Return
//         End 

int N2R_1(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
 const float tmp = (( R - Coefs[0] ) / Coefs[1] ) + 0.5f;
 return (int) tmp;
}

int N2R_2(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
 const float tmp0 = 1.0f - expf(Coefs[1]*(float) SeqLength - Coefs[2]);
 const float tmp1 = Coefs[0]*tmp0*(Coefs[4]*R + Coefs[3]);
 return (int) (tmp1 + 0.5f);
}

int N2R_3(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
 const float tmp0 = 1.0f - expf(Coefs[1]*(float) SeqLength - Coefs[2]);
 const float tmp1 = Coefs[0]*tmp0*(Coefs[4]*R + Coefs[3]);
 return (int) (tmp1 + AverageValue + 0.5f);
}

// Subroutine RtoN(N,R,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
// 
// Real              RNOP(KNPM,MAXN)
// 
// If     (IFUN.EQ.1) then 
//     R=RNOP(1,INOR)+RNOP(2,INOR)*Real(N)
// Else if(IFUN.EQ.2) then 
//     R=( Real(N) /
// *       ( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
// *       - RNOP(4,INOR) ) / RNOP(5,INOR)
// Else if(IFUN.EQ.3) then 
//     R=( (N-RAVE) /
// *       ( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
// *       - RNOP(4,INOR) ) / RNOP(5,INOR)
// End if 
// Return
// End

float R2N_1(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
  return Coefs[0] + Coefs[1]*((float) N);
}

float R2N_2(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
  const float tmp0 = 1.0f - expf(Coefs[1]*(float) SeqLength - Coefs[2]);
  const float tmp1 = ((float) N )/tmp0 - Coefs[3];
  
  return tmp1/Coefs[4];
}

float R2N_3(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength)
{
  const float tmp0 = 1.0f - expf(Coefs[1]*(float) SeqLength - Coefs[2]);
  const float tmp1 = ((float) N - AverageValue)/tmp0 - Coefs[3];
  
  return tmp1/Coefs[4];
}
