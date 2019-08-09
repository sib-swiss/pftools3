/*******************************************************
                        PFTOOLS
 *******************************************************
  May 29, 2013 random.h
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#ifndef _RANDOM_H_
#define _RANDOM_H_
#include <stdlib.h>
#include <stdbool.h>
#include "config.h"
#include "pfConfig.h"

#define MBIG            1000000000L
#define MSEED           161803398L
#define MZ              0L
#define FAC             (1.0f/1000000000.0f)
#define OVER            0.999999999f

struct Random {
  long int ma[56] __attribute__((aligned(16)));
  size_t inext;
  size_t inextp;
  long int Seed;
};

void InitializeGenerator(const long int Seed, struct Random * const sr);

extern inline bool __ALWAYS_INLINE YesOrNo(struct Random * const sr) 
{
    if (++sr->inext == 56) sr->inext=1;
    if (++sr->inextp == 56) sr->inextp=1;
    register long mj = sr->ma[sr->inext] - sr->ma[sr->inextp];
    if (mj < MZ) mj += MBIG;
    sr->ma[sr->inext]=mj;
    return (mj < MBIG/2) ? true : false;
}

extern inline long __ALWAYS_INLINE FlatDistributionValue(struct Random * const sr) 
{
  if (++sr->inext == 56) sr->inext=1;
  if (++sr->inextp == 56) sr->inextp=1;
  register long mj =sr->ma[sr->inext]-sr->ma[sr->inextp];
  if (mj < MZ) mj += MBIG;
  sr->ma[sr->inext]=mj;

  mj = ( mj > MBIG ) ? MBIG : mj;
  return mj;
}

extern inline float __ALWAYS_INLINE NormalizedFlatDistributionValue(struct Random * const sr)
{
  register long int * const ma = sr->ma;
  if (++(sr->inext) == 56)  sr->inext=1;
  if (++(sr->inextp) == 56) sr->inextp=1;
  register long int mj = ma[sr->inext] - ma[sr->inextp];
  if (mj < MZ) mj += MBIG;
  ma[sr->inext]=mj;

  if( mj >= MBIG) {
    return OVER;
  } else {
    return ((float) mj)*FAC;
  }
}
#endif /* _RANDOM_H_ */