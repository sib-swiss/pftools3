/*******************************************************
                        PFTOOLS
 *******************************************************
  May 29, 2013 random.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/

#include "random.h"
#include <string.h>
/* Taken from the Numerical receipe */

void InitializeGenerator(const long int Seed, struct Random * const sr)
{
  long mj,mk;
  size_t i,ii,k;

  /*
   * THIS IS NECESSARY AS IT SEEMS ALL VALUES ARE NOT SET WITHIN THE INITIALIZATION
   * SOME OLD ONE REMAINS !!!
   */
  memset(sr->ma, 0, 56*sizeof(long int));

  long int lSeed = Seed;
  long int * const ma = (long int *) sr->ma;

  mj = MSEED -(lSeed < 0 ? -lSeed : lSeed);
  mj %= MBIG;
  ma[55]=mj;
  mk=1;
  for (i=1;i<=54;++i) {
    ii=(21*i) % 55;
    ma[ii]=mk;
    mk=mj-mk;
    if (mk < MZ) mk += MBIG;
    mj=ma[ii];
  }
  for (k=1;k<=4;++k) {
    for (i=1;i<=55;++i) {
      ma[i] -= ma[1+(i+30) % 55];
      if (ma[i] < MZ) ma[i] += MBIG;
    }
  }

  sr->Seed    = Seed;
  sr->inext   = 0;
  sr->inextp  = 31;
}

