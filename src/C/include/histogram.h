/************************************************************
 * HMMER - Biological sequence analysis with profile HMMs
 * Copyright (C) 1992-2003 Washington University School of Medicine
 * All Rights Reserved
 *
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* histogram.c
 * SRE, Sat Jan 20 16:16:17 1996
 *
 * Accumulation, printing, and fitting of score histograms
 * from database searches.
 *
 * CVS $Id: histogram.c,v 1.18 2003/04/14 16:00:16 eddy Exp $
 ************************************************************
 * Basic API:
 *
 * struct histogram_s *h;
 *
 * h = AllocHistogram(min_hint, max_hint, lumpsize);
 *
 * while (getting scores x) AddToHistogram(h, x);
 *
 * ExtremeValueFitHistogram(h, high_hint);
 * PrintASCIIHistogram(fp, h);
 * FreeHistogram(h);
 */
#ifndef _HISTOGRAM_H_
#define _HISTOGRAM_H_
#include <stdio.h>
/* Structure: histogram_s
 *
 * Keep a score histogram.
 *
 * The main implementation issue here is that the range of
 * scores is unknown, and will go negative. histogram is
 * a 0..max-min array that represents the range min..max.
 * A given score is indexed in histogram array as score-min.
 * The AddToHistogram() function deals with dynamically
 * resizing the histogram array when necessary.
 */
struct histogram_s {
  int *histogram;               /* counts of hits                     */
  int  min;                     /* elem 0 of histogram == min         */
  int  max;                     /* last elem of histogram == max      */
  int  highscore;               /* highest active elem has this score */
  int  lowscore;                /* lowest active elem has this score  */
  int  lumpsize;                /* when resizing, overalloc by this   */
  int  total;                   /* total # of hits counted            */

  float *expect;                /* expected counts of hits            */
  int    fit_type;              /* flag indicating distribution type  */
  float  param[3];              /* parameters used for fits           */
  float  chisq;                 /* chi-squared val for goodness of fit*/
  float  chip;                  /* P value for chisquared             */
};
#define HISTFIT_NONE     0      /* no fit done yet               */
#define HISTFIT_EVD      1      /* fit type = extreme value dist */
#define HISTFIT_GAUSSIAN 2      /* fit type = Gaussian           */
#define EVD_MU           0      /* EVD fit parameter mu          */
#define EVD_LAMBDA       1      /* EVD fit parameter lambda      */
#define EVD_WONKA        2      /* EVD fit fudge factor          */
#define GAUSS_MEAN       0      /* Gaussian parameter mean       */
#define GAUSS_SD         1      /* Gaussian parameter std. dev.  */


extern struct histogram_s *AllocHistogram(int min, int max, int lumpsize);
extern void FreeHistogram(struct histogram_s *h);
extern void UnfitHistogram(struct histogram_s *h);
extern void AddToHistogram(struct histogram_s *h, float sc);
extern void PrintASCIIHistogram(FILE *fp, struct histogram_s *h);
extern void PrintXMGRHistogram(FILE *fp, struct histogram_s *h);
extern void PrintXMGRDistribution(FILE *fp, struct histogram_s *h);
extern void PrintXMGRRegressionLine(FILE *fp, struct histogram_s *h);
extern void EVDBasicFit(struct histogram_s *h);
extern int  ExtremeValueFitHistogram(struct histogram_s *h, int censor, int lowscore, float high_hint);
extern void ExtremeValueSetHistogram(struct histogram_s *h, float mu, float lambda, float low, float high, int ndegrees);
extern int  GaussianFitHistogram(struct histogram_s *h, float high_hint);
extern void GaussianSetHistogram(struct histogram_s *h, float mean, float sd);
extern double EVDDensity(float x, float mu, float lambda);
extern double EVDDistribution(float x, float mu, float lambda);
extern double ExtremeValueP (float x, float mu, float lambda);
extern double ExtremeValueP2(float x, float mu, float lambda, int N);
extern double ExtremeValueE (float x, float mu, float lambda, int N);
extern float  EVDrandom(float mu, float lambda);
extern int    EVDMaxLikelyFit(float *x, int *y, int n,float *ret_mu, float *ret_lambda);
extern int    EVDCensoredFit(float *x, int *y, int n, int z, float c, float *ret_mu, float *ret_lambda);
extern void   Lawless416(float *x, int *y, int n, float lambda, float *ret_f, float *ret_df);
extern void   Lawless422(float *x, int *y, int n, int z, float c, float lambda, float *ret_f, float *ret_df);


#endif /* _HISTOGRAM_H_ */
