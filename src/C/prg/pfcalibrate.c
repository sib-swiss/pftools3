/*******************************************************
 *                        PFTOOLS
 *******************************************************
 *  Oct 8, 2015 pfcalibrate.c
 *******************************************************
 * (C) 2013 SIB Swiss Institute of Bioinformatics
 *     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#ifdef USE_AFFINITY
#include <sched.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <locale.h>

#if !defined(__USE_WINAPI__)
#include <pthread.h>
#else
#include <windows.h>
#endif

#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif

#include <inttypes.h>
#include <getopt.h>
#include <limits.h>

#ifdef USE_AFFINITY
# include <unistd.h>
#endif

#ifdef USE_PDF
#ifdef USE_PLPLOT
#define __WITH_REPORT__
#include <float.h>
#include "plplot/plplot.h"
#include "hpdf.h"
#include <setjmp.h>
#endif
#endif

#define HEADER "%--------------------------------------------------------------------%\n"\
"|                         PFCALIBRATE v"PF_VERSION "                         |\n"\
"%--------------------------------------------------------------------%\n"\
"| Built on " __DATE__ " at " __TIME__ ".                                  |\n"

#define BASE_E 2.71828182845904524
#define NUMBER_OF_METHODS 9

#include "pfProfile.h"
#include "pfSequence.h"
#include "system.h"
#include "histogram.h"

#define _NEEDS_HEURISTIC_CALIBRATION_
#define _NEEDS_HEURISTIC_
#define _NEEDS_FILTER_
#include "threads.h"

struct Method {
	char name[16];
	char NormalizationText[32];
	char LogText[32];
	int (*fct)(const int * const restrict FilterScores, const size_t SequenceCount, void * restrict const DATA);
};

enum StatisticMethod { AVERAGE, MINIMUM, MAXIMUM };

static SystemInfo System;
static int indexDatabase = 0;
#ifdef USE_AFFINITY
# ifdef __WITH_REPORT__
static const char opt_to_test[] = "hsVd:D:F:H:m:B:!:?:z:0:1:2:3:4:5:7:8:9:Q:E:N:R:L:M:o:t:T:r:fbq:k+z";
# else
static const char opt_to_test[] = "hsVd:D:F:H:m:B:!:?:z:0:1:2:3:4:5:7:8:9:Q:E:N:R:L:M:o:t:T:fbq:k+z";
# endif
#else
# ifdef __WITH_REPORT__
static const char opt_to_test[] = "hsVd:D:F:H:m:B:!:?:z:0:1:2:3:4:5:7:8:9:Q:E:N:R:L:M:o:t:T:r:+z";
# else
static const char opt_to_test[] = "hsVd:D:F:H:m:B:!:?:z:0:1:2:3:4:5:7:8:9:Q:E:N:R:L:M:o:t:T:+z";
# endif
#endif

static int From_perl_script = 0;
static const struct option long_options[] =
{
	/*
	 * These options set a flag.
	 */

	/*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	{"help",               			no_argument,       	0,	'h'},
	{"sse2",										no_argument,				0,	's'},
	{"verbose",									no_argument,				0,	'V'},
	/* Print ouptut methods */
	{"dump-filter-scores",			required_argument,	0,	'D'},
	{"dump-heuristic-scores",		required_argument,	0,	'd'},
	#ifdef __WITH_REPORT__
	{"report",									required_argument,	0,	'r'},
	#endif
	/* Database indexing options */
	{"filter-db",								required_argument,	0,	'F'},
	{"heuristic-db",						required_argument,	0,	'H'},
	/* Others */
	{"method",									required_argument,	0,	'm'},
	{"pfscale-logarithm-base",	required_argument,	0,	'B'},
	{"pfscale-upper-threshold", required_argument,	0,	'!'},
	{"pfscale-lower-threshold", required_argument,	0,	'?'},

	{"evd-bin-center",          required_argument,	0,	'3'},
	{"evd-bin-width",						required_argument,	0,	'2'},
	{"evd-average-length",			required_argument,	0,	'4'},
	{"evd-tail-area",						required_argument,	0,	'5'},

	{"res_count",								required_argument,	0,	'z'},
	{"first_rank",							required_argument,	0,	'0'},
	{"last_rank",								required_argument,	0,	'1'},

	{"seed",										required_argument,	0,	'7'},
	{"pam-distance",						required_argument,	0,	'8'},
	{"profile-sampling",				required_argument,	0,	'9'},
	{"profile-only-on-match",		no_argument,				0,	'z'},
	{"quantile",								required_argument,	0,	'Q'},
	{"max-exponential",					required_argument,	0,	'E'},

	{"normalized_score",				required_argument,	0,	'N'},
	{"raw_score",								required_argument,	0,	'R'},
	{"mode",										required_argument,	0,	'M'},
	{"output",									required_argument,	0,	'o'},
	{"perl",										no_argument,				0,	'+'},
	/* SMP options */
	{"nthreads",								required_argument,	0,	't'},
	{"max-heuristic-nthreads",	required_argument,	0,	'T'},
#ifdef USE_AFFINITY
	{"no-affinity",							no_argument,				0,	'f'},
	{"split", 									no_argument,				0,	'b'},
	{"thread-affinity",					required_argument,	0,	'q'},
	{"no-shared-core",					no_argument,				0,	'k'},
#endif
	{0, 0, 0, 0}
};

static void* (*thread_heuristic_ptr)(void*);

static _Bool OutputVerbose = false;
static bool PerformHeuristic = false;

#ifdef __WITH_REPORT__
static char * HeuristicTempName;	/* Temporary file name for heuristic plot */
static char * FilterTempName;		/* Temporary file name for filter plot */
static _Bool CreatePDFReport=false;	/* Create or not a PDF report */
static PLINT red[]   = { 255, 0, 230, 255,   0,   0,  39, 125,   8,   0};
static PLINT green[] = { 255, 0, 230,   0, 255,   0,  80,   0,   0,   0 };
static PLINT blue[]  = { 255, 0, 230,   0,   0, 255, 204, 125,   0, 255 };
jmp_buf env;

#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler  (HPDF_STATUS   error_no,
								HPDF_STATUS   detail_no,
								void         *user_data)
{
	printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT)error_no,
					(HPDF_UINT)detail_no);
	longjmp(env, 1);
}
#endif

/* PFScale */
double logarithmic_base = 10.0; 	/* logarithmic base of parameters */
double upperProbRange = 0.0001;	 	/* upper threshold of probability range */
double lowerProbRange = 0.000001;	/* lower threshold of probability range */
float first_rank = 0.0f;
float last_rank  = 0.0f;

/* EVD */
float binCenter = 0.0f;
float binWidth  = 1.0f;
float avgLength = -1.0f;
float resCount = 0.0f;
float tailArea = 0.2f;

/* Heuristic */


/* Statistics */
union __32bitData { float FloatScores; int SignedScores;};


extern int QuantileRegression(const unsigned int * const restrict HeuristicScores, const int * const restrict FilterScore,
															float coefs[], const double quantile, const int size);

static void __attribute__((noreturn)) Usage(FILE * stream)
{
	fprintf(stream,
					"Calibrate a profile:\n"
					"pfcalibrate [options] --method <method> --mode <mode> profile\n"
					"E.g. pfcalibrate --method pfscale --filter-db window20.seq --heuristic-db profile --profile-sampling 25 <profile>\n"
					"     pfcalibrate --filter-db window20.seq --heuristic-db profile --profile-sampling 25 <profile>\n\n"
					"Options:\n"
					"  Profile\n"
					"    --mode                              [-M] : mode to use for normalization (default 1)\n\n"
					"  Filter calibration\n"
					"    --filter-db <file>                  [-F] : FASTA sequence database for filter calibration\n\n"
					"    --method <string>                   [-m] : specify the method to use (default evd_tail)\n"
					"      == pfscale\n"
					"        --pfscale-logarithm-base <string>    : logarithmic base of parameters\n"
					"                                               you may use 'e' or 'E' for natural base\n"
					"                                               default is base %2.0lf.\n"
					"        --pfscale-upper-threshold <float>    : upper threshold of probability range\n"
					"                                               default is %lg.\n"
					"        --pfscale-lower-threshold <float>    : lower threshold of probability range\n"
					"                                               default is %lg.\n"
					"      == evd_full\n"
					"        --evd-bin-center <float>             : default is %.1f.\n"
					"        --evd-bin-width <float>              : default is %.1f.\n"
					"        --evd-average-length <float>         : default is %.1f.\n"
					"      == evd_peak\n"
					"      == evd_tail\n"
					"        --evd-tail-area <float>              : default is %.1f.\n\n"
					"  UNDESCRIPTIVE OPTIONS:\n"
					"    --res_count <float>                      : default is %.1f\n"
					"    --first_rank <float>                     : default is %.1f\n"
					"    --last_rank <float>                      : default is %.1f\n\n"
					"  Heuristic calibration\n"
					"    --heuristic-db <file>               [-H] : FASTA sequence database file name for heuristic calibration\n"
					"                                               if no such file exists, you may specify 'profile'\n"
					"                                               as database file name to trigger automatic generation\n"
					"                                               of sequences based on the profile concensus and scoring.\n"
					"    --seed <uint>                            : Random generator initial seed (default 1234)\n"
					"    --profile-sampling <uint>                : specifies the number of seed sequences that will be\n"
					"                                               generated from the profile. Only required when a profile\n"
					"                                               is used for the generation of mutated sequences;\n"
					"                                               default is 50.\n"
					"    --profile-only-on-match                  : profile sequence generation using only Match scores\n"
					"    --pam-distance [start,stop,step,sampling]: set PAM distance and interval for mutation of seed sequences,\n"
					"                                               default is [%zu,%zu,%zu,%zu]\n"
					"    --quantile <float>                       : quantile value (default 0.05)\n"
					"    --max-exponential <float>                : maximum exponent used in exponential when\n"
					"                                               generating sequences from the profile.\n"
					"                                               default is 8.0.\n\n"
					"  Database\n"
					"    --database-index <file>             [-i] : use indices stored in given file\n\n"
					"  Optimizations\n"
					"    --sse2                                   : enforce SSE 2 only instruction set\n"
					"    --nthreads <uint>                   [-t] : max number of threads to use\n"
					"                                               default to all available cores\n"
					"    --max-heuristic-nthreads <uint>     [-T] : max number of threads to use for\n"
					"                                               heuristic phase only. (IO bounds)\n"
					"                                               default to all available cores\n"
#ifdef USE_AFFINITY
					"    --no-affinity                            : disable CPU affinity file\n"
					"    --thread-affinity                        : file containing thread mask,\n"
					"                                               one row for one thread\n"
					"    --no-shared-core                         : Prevent core resource sharing\n"
					"    --split                                  : if both SSE 2 & 4.1 are available,\n"
					"                                               split half-half using linked resources\n\n"
#else
					"\n"
#endif
					"  Printing output\n"
					"    --dump-filter-scores <file>         [-D] : output filter scores to given file\n"
					"    --dump-heuristic-scores <file>      [-d] : output filter - heuristic scores to given file\n"
					"    --output <file>                     [-o] : output profile to given file (default stdout)\n"
#ifdef __WITH_REPORT__
					"    --report <file>                     [-r] : generate PDF report to given file\n"
#endif
					"\n"
					"  Other\n"
					"    --verbose                           [-V] : verbose on stderr\n"
					"    --help                              [-h] : output command help\n\n"
					" Version " PF_VERSION " built on " __DATE__ " at " __TIME__ ".\n",
				 logarithmic_base, upperProbRange, lowerProbRange,
				 binCenter, binWidth, avgLength, tailArea,
				 resCount, first_rank, last_rank,
				 PamDistanceStart, PamDistanceStop, PamDistanceStep, PamSampling);
	exit(0);
}

static union __32bitData Statistics(const int * restrict FilterScores, const size_t Count, enum StatisticMethod M)
{
	register int value = 0;
	register size_t N = Count;
	union __32bitData result;
	switch (M) {
		case AVERAGE:
		{
			for (size_t i=0; i<N; ++i) value += FilterScores[i];
			result.FloatScores = (float) value / (float) N;
			return result;
		}
		break;
		case MINIMUM:
		{
			value = INT_MAX;
			for (size_t i=0; i<N; ++i) value = FilterScores[i] < value ? FilterScores[i] : value;
			result.SignedScores = value;
			return result;
		}
		break;
		case MAXIMUM:
		{
			value = INT_MIN;
			for (size_t i=0; i<N; ++i) value = FilterScores[i] > value ? FilterScores[i] : value;
			result.SignedScores = value;
			return result;
		}
		break;
		default:
			fputs("Call to Statistics with unknown method\n", stderr);
			exit(1);
	}
}

static int pfscale(const int * const restrict FilterScores, const size_t SequenceCount, void * restrict const DATA)
{
	const size_t N = SequenceCount;
	{
		const float ftmp = (float) SequenceCount;
		if (upperProbRange == 0.0f) {
			if (last_rank == 0.0) last_rank  = roundf(ftmp/20.0f);
			upperProbRange = last_rank/resCount;
		}
		if (lowerProbRange == 0.0f) {
			if (first_rank == 0.0f) first_rank = roundf(ftmp/500.0f);
			lowerProbRange = first_rank/resCount;
		}
	}
	const double RL   = 1.0/log(logarithmic_base);
	const double RDBS = RL*log((double) resCount);
	const double Emax = lowerProbRange != 0.0 ? -RL*log(lowerProbRange) : 100.0;
	const double Emin = -RL*log(upperProbRange);
	double * const Weights     = (double*) _mm_malloc(SequenceCount*sizeof(double), 16);
	double * const Frequencies = (double*) _mm_malloc(SequenceCount*sizeof(double), 16);
	if (Weights == NULL || Frequencies == NULL) {
		fputs("Unable to allocate memory for weights and Frequencies\n", stderr);
		return 1;
	}

	for (size_t iseq=0; iseq<N; ++iseq) {
		const double dtmp = (double) (1+iseq);
		Frequencies[iseq] = RDBS-RL*log(dtmp - 0.5);
		Weights[iseq]     = dtmp*0.5;
	}

	double TotalWeight      = 0.0;
	double AverageFrequency = 0.0;
	double AverageS         = 0.0;

	for (size_t iseq=0; iseq<N; ++iseq) {
		if (Frequencies[iseq] >= Emin && Frequencies[iseq] <= Emax ) {
			AverageS           += Weights[iseq]*(double)FilterScores[iseq];
			AverageFrequency   += Weights[iseq]*Frequencies[iseq];
			TotalWeight        += Weights[iseq];
		}
	}
	AverageFrequency /= TotalWeight;
	AverageS         /= TotalWeight;

	double XSV = 0.0, XFV = 0.0, XCO = 0.0;
	for (size_t iseq=0; iseq<N; ++iseq) {
		if (Frequencies[iseq] >= Emin && Frequencies[iseq] <= Emax ) {
			const double dtmp1 = (double) FilterScores[iseq] - AverageS;
			const double dtmp2 = Frequencies[iseq] - AverageFrequency;
			XSV += Weights[iseq]*(dtmp1*dtmp1);
			XFV += Weights[iseq]*(dtmp2*dtmp2);
			XCO += Weights[iseq]*dtmp1*dtmp2;
		}
	}
	XSV = sqrt(XSV/TotalWeight);
	XFV = sqrt(XFV/TotalWeight);
	XCO = (XCO/TotalWeight)/(XSV*XFV);


	const float XB = (float) (XCO*XFV/XSV);
	const float XA = (float) (AverageFrequency-XB*AverageS);

	//   fprintf(stdout, "# -LogP = %8.4f + %12.8f * raw-score\n#\n"
	// 		  "#   rank  raw-score   -logFreq   -logProb\n"
	// 		  "#\n", XA, XB);
	//   for (size_t iseq=0; iseq<N; ++iseq) {
	//     const float dtmp = XA + XB*(float)FilterScores[iseq];
	//     fprintf(stdout,"%8lu %10i %10.4f %10.4f %10.4f\n", 1+iseq, FilterScores[iseq], Frequencies[iseq], dtmp, Weights[iseq]);
	//   }

	float * const restrict results = DATA;
	results[0] = XA;
	results[1] = XB;

#ifdef __WITH_REPORT__
	/* Generate plot to be included in report */
	if (CreatePDFReport) {
		FilterTempName = tempnam(0,"F");
		if (! FilterTempName) {
			fputs("Unable to create a temporary file name\n",stderr);
		}
		else {
#ifndef PL_DOUBLE
#error pfcalibrate does not support plplot in float format yet
#else
			double * const X = (double*) _mm_malloc(SequenceCount*sizeof(double),16);
			if (X == NULL) {
				fputs("Cannot allocate memory for graph.\n", stderr);
				exit(1);
			}
			size_t i = 4;
			__m128 __maxX = _mm_set1_ps(FLT_MIN);
			__m128 __minX = _mm_set1_ps(FLT_MAX);
			__m128d __maxY = _mm_set1_pd(DBL_MIN);
			// Crappy code SSE 4.1 and AVX would be much better !!!
			while (i < SequenceCount) {
				__m128d __fY1 = _mm_load_pd(&Frequencies[i-4]);
				__m128d __fY2 = _mm_load_pd(&Frequencies[i-2]);
				__m128 __fX = _mm_cvtepi32_ps(*(__m128i*) &FilterScores[i-4]);
				i += 4;
				__maxX = _mm_max_ps(__maxX, __fX);
				__minX = _mm_min_ps(__minX, __fX);
				__fY1  = _mm_max_pd(__fY1, __fY2);
				__maxY = _mm_max_pd(__maxY, __fY1);

				__m128d __X1 = _mm_cvtps_pd(__fX);
				__m128d __X2 = _mm_cvtps_pd(_mm_movehl_ps(__fX, __fX));

				_mm_store_pd(&X[i-8], __X1);
				_mm_store_pd(&X[i-6], __X2);
			}
			i -= 4;
			__m128 __zero = _mm_setzero_ps();
			while (i < SequenceCount) {
				__m128 __fX  = _mm_cvtsi32_ss(__zero, FilterScores[i]);
				__maxY = _mm_max_sd(__maxY, *(__m128d*) &Frequencies[i]);
				i += 1;
				__maxX = _mm_max_ps(__maxX, __fX);
				__minX = _mm_min_ps(__minX, __fX);
				__m128d __X1 = _mm_cvtss_sd(_mm_castps_pd(__zero), __fX);

				_mm_store_sd( &X[i-1], __X1);

			}
#endif

			// movhlps xmm1, xmm0	Move top two floats to lower part of xmm1
			__m128 __tMaxX = _mm_movehl_ps(__maxX, __maxX);
			__m128d __tMaxY = _mm_unpackhi_pd(__maxY, __maxY);
			__m128 __tMinX = _mm_movehl_ps(__minX, __minX);
			// maxps   xmm0,xmm1	Get maximum/minimum of the two sets of floats
			__maxX = _mm_max_ps(__maxX, __tMaxX);
			__maxY = _mm_max_pd(__maxY, __tMaxY);
			__minX = _mm_min_ps(__minX, __tMinX);
			// shufps  xmm1,xmm0,$55	Move second float to lower part of xmm1 0101 0101 = 55
			__tMaxX = _mm_shuffle_ps(__maxX, __maxX, 0b01010101);
			__tMinX = _mm_shuffle_ps(__minX, __minX, 0b01010101);
			// maxps   xmm0,xmm1	Get minimum of the two remaining floats
			__maxX = _mm_max_ps(__maxX, __tMaxX);
			__minX = _mm_min_ps(__minX, __tMinX);

			float maxX, minX;
			double maxY;
			_mm_store_ss(&maxX, __maxX);
			_mm_store_sd(&maxY, __maxY);
			_mm_store_ss(&minX, __minX);

			PLFLT x_max = (PLFLT) (1.1*maxX);
			PLFLT x_min = (PLFLT) (1.1*minX);
			PLFLT y_max = (PLFLT) (1.1*maxY);
			PLFLT y_min = (PLFLT) 0.0;

			plsdev("jpeg");
			plsfnam(FilterTempName);
			plscmap0(red, green, blue, 10);
			plspage(0,0,1024,800,0,0);

			plinit();
			plfont(2);

			plenv(x_min, x_max, y_min, y_max, 0, 0);
			if (logarithmic_base == BASE_E)
				pllab("Scores", "-Ln(#gn)", "");
			else
				pllab("Scores", "-Log#d10#u(#gn)", "");
			plcol0(2);

			plcol0(5);plssym(0,.75);
			plpoin((PLINT) SequenceCount, X, Frequencies, 17);

			plcol0(3);
			const PLFLT ymin = (PLFLT) (XA + x_min * XB);
			const PLFLT ymax = (PLFLT) (XA + x_max * XB);
			pljoin(x_min, ymin, x_max, ymax);

			plend();
	#ifdef PL_DOUBLE
			_mm_free(X);
	#endif


		}
	}
	#endif
	_mm_free(Weights);
	_mm_free(Frequencies);
	return 0;
}

static int evd_full(const int * const restrict FilterScores, const size_t SequenceCount, void * restrict const DATA)
{
	const size_t N = SequenceCount;
	int min = INT_MAX;
	int max = INT_MIN;
	for (size_t i=0; i<N; ++i) {
		const int value = FilterScores[i];
		min = (min > value) ? value : min;
		max = (max < value) ? value : max;
	}
	//   min = (int) ((((float)min - binCenter)/binWidth + (min < 0 ? -0.5f : 0.5f)) * binWidth + binCenter);
	//   max = (int) ((((float)max - binCenter)/binWidth + (max < 0 ? -0.5f : 0.5f)) * binWidth + binCenter);
	if (OutputVerbose) fprintf(stderr, "EVD Full: min(%i) max(%i) bin center (%lf) bin width (%lf)\n", min, max, binCenter, binWidth);
	struct histogram_s * H = AllocHistogram(min, max, 16);
	for (size_t i=0; i<N; ++i) {
		//     float value = (((float) FilterScores[i] - binCenter)/binWidth + (FilterScores[i] < 0 ? -0.5f : 0.5f))* binWidth + binCenter;
		const float value = (float) FilterScores[i];
		AddToHistogram(H, value);
	}

	ExtremeValueFitHistogram(H, 0, min, 99999.);

	float * restrict const result = DATA;
	result[0] = (float) ((log(avgLength) - H->param[EVD_MU]*H->param[EVD_LAMBDA])/log(logarithmic_base));
	result[1] = (float) (H->param[EVD_LAMBDA]/log(logarithmic_base));

	//   if (OutputVerbose) {
	//     PrintASCIIHistogram(stdout, H);
	//     fprintf(stdout, "# -LogP = %8.4f + %12.8f * raw-score\n#\n", result[0], result[1]);
	//   }
	FreeHistogram(H);
}

static int evd_peak(const int * const restrict FilterScores, const size_t SequenceCount, void * restrict const DATA)
{
	const size_t N = SequenceCount;
	int min = INT_MAX;
	int max = INT_MIN;
	for (size_t i=0; i<N; ++i) {
		const int value = FilterScores[i];
		min = (min > value) ? value : min;
		max = (max < value) ? value : max;
	}
	//   min = (int) ((((float)min - binCenter)/binWidth + (min < 0 ? -0.5f : 0.5f)) * binWidth + binCenter);
	//   max = (int) ((((float)max - binCenter)/binWidth + (max < 0 ? -0.5f : 0.5f)) * binWidth + binCenter);
	if (OutputVerbose) fprintf(stderr, "EVD Full: min(%i) max(%i) bin center (%lf) bin width (%lf)\n", min, max, binCenter, binWidth);
	struct histogram_s * H = AllocHistogram(min, max, 16);
	for (size_t i=0; i<N; ++i) {
		//     float value = (((float) FilterScores[i] - binCenter)/binWidth + (FilterScores[i] < 0 ? -0.5f : 0.5f))* binWidth + binCenter;
		const float value = (float) FilterScores[i];
		AddToHistogram(H, value);
	}

	ExtremeValueFitHistogram(H, 2, (min+max)/2, 99999.);

	float * restrict const result = DATA;
	result[0] = (float) ((log(avgLength) - H->param[EVD_MU]*H->param[EVD_LAMBDA])/log(logarithmic_base));
	result[1] = (float) (H->param[EVD_LAMBDA]/log(logarithmic_base));

	//   if (OutputVerbose) {
	//     PrintASCIIHistogram(stdout, H);
	//     fprintf(stdout, "# -LogP = %8.4f + %12.8f * raw-score\n#\n", result[0], result[1]);
	//   }
	FreeHistogram(H);
}

static int evd_tail(const int * const restrict FilterScores, const size_t SequenceCount, void * restrict const DATA)
{
	const size_t N = SequenceCount;
	int min = INT_MAX;
	int max = INT_MIN;
	for (size_t i=0; i<N; ++i) {
		const int value = FilterScores[i];
		min = (min > value) ? value : min;
		max = (max < value) ? value : max;
	}

	if (OutputVerbose) fprintf(stderr, "EVD Full: min(%i) max(%i) bin center (%lf) bin width (%lf)\n", min, max, binCenter, binWidth);
	struct histogram_s * H = AllocHistogram(min, max, 100);
	for (size_t i=0; i<N; ++i) {
		const float value = (float) FilterScores[i];
		AddToHistogram(H, value);
	}

	const int TargetArea = (int) (tailArea*(float)SequenceCount);
	int sum = 0;
	int sc = H->highscore;
	const int hmin = H->min;
	while (sc > H->lowscore) {
		sum += H->histogram[sc-- - hmin];
		if (sum >= TargetArea) break;
	}

	ExtremeValueFitHistogram(H, 1, sc, 99999.);
	float * restrict const result = DATA;
	result[0] = (logf(avgLength) - H->param[EVD_MU]*H->param[EVD_LAMBDA])/logf(logarithmic_base);
	result[1] = H->param[EVD_LAMBDA]/logf(logarithmic_base);

	if (OutputVerbose) {
		PrintASCIIHistogram(stdout, H);
		fprintf(stdout, "# -LogP = %8.4f + %12.8f * raw-score\n#\n", result[0], result[1]);
	}

	FreeHistogram(H);
}

static struct Method Methods[NUMBER_OF_METHODS] = {
	{ "consensus",    "", "*", NULL },
	{ "minimum",      "", "*", NULL },
	{ "maximum",      "", "*", NULL },
	{ "average",      "", "*", NULL },
	{ "remove_level", "",       "",  NULL },
	{ "pfscale",      "-Log10", "*", pfscale },
	{ "evd_full",     "-Log10", "*", evd_full },
	{ "evd_peak",     "-Log10", "*", evd_peak },
	{ "evd_tail",     "-Log10", "*", evd_tail }
};

static void iqsort(int * const restrict data, const int N)
{
	int i,j;
	int itmp;
	int t,v;

	if (N<=1) return;
	v = data[0];
	i = 0;
	j = N;
	for (;;)
	{
		while(data[++i] > v && i <  N) {}
		while(data[--j] < v) {}
		if (i>=j)
		{ break; }
		else
		{
			t = data[i];
			data[i] = data[j];
			data[j] = t;
		}
	}
	t = data[i-1];
	data[i-1] = data[0];
	data[0] = t;
	iqsort(data,i-1);
	iqsort(data+i,N-i);
}

static int dumpScores(const int FilterScores[], const int HeuristicScores[],
                      const char * const restrict SequenceFile, FILE * const restrict stream,
											const struct Profile * const restrict prf, const float * const restrict Coeffs,
											const DBSequence_t * restrict const DB)
{
	Sequence SeqData;

	/* Allocate memory to hold sequence */
	SeqData.Data.Memory = malloc(DB->MaxSequenceSize*sizeof(unsigned char));
	if (SeqData.Data.Memory == NULL) {
		fputs("Thread Cannot allocate memory for sequence.\n", stderr);
		return 1;
	}

	/* Open sequence file*/
	SETUP_DATABASE_ACCESS(SequenceFile);
	register const size_t N = DB->SequenceCount;

	fputs("# Sequence Normalized_score Filter_score Heuristic_score\n", stream);
	for (size_t i=0; i<N; ++i) {
		PFSequence * const PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB->DataPtr[i].Sequence));
		const float NormalizedScore = prf->RawToNormalized(FilterScores[i], Coeffs, 0, 0);
		fprintf(stream, "%s\t%f\t%i\t%i\n", &SeqData.Data.Header[1], NormalizedScore, FilterScores[i], HeuristicScores[i]);
	}

	/* close sequence file */
	UNSET_DATABASE_ACCESS();
	free(SeqData.Data.Memory);

	return 0;
}


int main(int argc, char *argv[])
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL STRUCTURES
	////////////////////////////////////////////////////////////////////////////////////////////////

	struct Profile * restrict prf;	/* Profile */
	DBSequence_t DB;    			      /* Sequence Database File */
	Sequence SeqData;			          /* Sequence data to work on */
	struct timeval _t0, _t1;		    /* Timing structures */

	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL DATA
	////////////////////////////////////////////////////////////////////////////////////////////////

	PFSequence * PFSeq;															/* Pointer to translated alphabet sequence */
	size_t nCPUs=0;																	/* number of threads */
	size_t nCPUsHeuristic=0;              					/* maximum number of threads for heuristic phase */
	#ifdef __NUMA__
	size_t nNodes=0;																/* number of NUMA nodes (max is NUMA API dictated) */
	size_t nThreadsPerNodes=0;											/* maximum number of threads per NUMA node */
	#endif
	size_t HeuristicCounter = 0;										/* number of sequences passing heuristic */
	size_t FilterCounter = 0;												/* number of sequences passing filter */
	int res, Score;
	int HeuristicCutOff = -1;												/* Default heuristic cutoff from command line, if not zero then enforces that value*/
	_Bool modeGiven = false;												/* True if --mode is given */
	int Mode = 1;																		/* Default Mode for normalization from command line, if not zero then enforces that value */
	int * restrict FilterScores = NULL;   					/* Array of filter scores */
	unsigned int * restrict HeuristicScores = NULL;	/* Array of heuristic scores used when dumping both filter and heuristic scores */
	unsigned char * restrict PamTypes = NULL;				/* Array of pam types used for colr setting in graphs */
	char * ProfileFile = NULL;											/* Profile file */
	char *DBFileName = NULL;												/* Sequence file containing random sequences */
	char *SeqDB = NULL;															/* Sequence file containing sequences used to define profile*/
	int method = 8;																	/* Index to the method to apply, default is evd_tail */
	_Bool HasFilterData = false;										/* Has there been a filter FASTA database file */
	_Bool HasHeuristicData = false;									/* Has there been an heuristic FASTA database file */
	_Bool ProfileChanged = false;										/* Indicates modification is to be applied */
	SCutOffItem NewCutItem;													/* Holds the new cutoff */
	SNormalizationItem NewFilterNormItem;						/* Holds the new normalization for filter*/
	SNormalizationItem NewHeuristicNormItem;				/* Holds the new normalization for heuristic*/
	char * outputFile = NULL;												/* Output profile to given file */
	char * HeuristicScoresOutputFile = NULL;				/* Output filter - heuristic scores to given file */
	char * FilterScoresOutputFile = NULL;						/* Output filter scores to given file */
	double quantile = 0.05;													/* Quantile value used in heuristic calibration */
	float MaxExponent = 8.0f;												/* Maximum exponent to the exponential in generating sequences from the profile */
	unsigned int ProfileSampling = 50;							/* When heuristic-db = profile, this gives the number of generated sequences */
	enum Version HeuristicVersion = SSE2;						/* Trigger which sse version to use in the heuristic */
	enum Version OtherVersion = SSE2;								/* Trigger which sse version to use in the filter and alignment */
	enum GeneratedSequenceOptions GenOpts = GENERATE_MATCH | GENERATE_INSERTION | GENERATE_DELETION;

	size_t * shares = 0;
	struct ThreadData *threads_arg = NULL;					/* Allocate stack memory for posix thread structures */
	#if !defined(__USE_WINAPI__)
	pthread_t *threads = 0;
	#else
	HANDLE * threads = 0;
	#endif
#ifdef USE_AFFINITY
	Affinity_Mask_t * Thread_masks[2] = {0,0};						/* Define variables to hold thread affinity mask */
	unsigned int Thread_count[2] = {0,0};
	pthread_attr_t * restrict threads_attr = NULL;
	char buffer[128] __attribute__((aligned(16)));	/* buffer to read affinity file mask */
	_Bool noAffinity = false;												/* disable use of cpu affinity */
	_Bool split = false;
	_Bool noSharedCore = false;											/* Prevent hyperthreading or AMD compute unit to share resources */
	_Bool GivenAffinityFile = false;								/* File holding a mask for each thread */
	char * AffinityMaskFileName;										/* Name of affinity mask file provided by option m */
#endif
#ifdef __WITH_REPORT__
	char * ReportFileName=0;												/* Base file name for the PDF report */
	size_t FilterDBSequenceCount=0;									/* Used in PDF creation */
	size_t FilterDBSize = 0;												/* Used in PDF creation */
	size_t HeuristicDBSequenceCount=0;							/* Used in PDF creation */
	size_t HeuristicDBSize = 0;											/* Used in PDF creation */
#endif
	_Bool HeuristicSequencesFromProfile = false;		/* Used to clean generated profile sequence file */


	////////////////////////////////////////////////////////////////////////////////////////////////
	// SYSTEM ARCHITECTURE ANALYSIS
	////////////////////////////////////////////////////////////////////////////////////////////////
	getSystemInfo(&System);

	/* Check for minimum requirement */
	if (!(System.Extensions & MM_SSE2)) {
		fputs("pfsearch requires at least a CPU capable of SSE 2.\n", stderr);
		exit(1);
	}

	/* Allow fast SSE 4.1 extensions ? */
	if (System.Extensions & MM_SSE41) {
		xali1_ptr = xali1_sse41;
		HeuristicVersion = SSE41;
		OtherVersion = SSE41;
	}
	else {
		xali1_ptr = xali1_sse2;
		HeuristicVersion = SSE2;
		OtherVersion = SSE2;
	}

	/* Thousands separator requires this */
	setlocale(LC_ALL, "");

	////////////////////////////////////////////////////////////////////////////////////////////////
	// OPTIONS
	////////////////////////////////////////////////////////////////////////////////////////////////
	while (1) {
		/* getopt_long stores the option index here. */
		int option_index = 0;

		const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) break;
		switch (c) {
#ifdef USE_AFFINITY
			case 'f':
				noAffinity = true;
				break;
			case 'k':
				noSharedCore = true;
				break;
			case 'b':
				if (HeuristicVersion == SSE41) {
					split = true;
				} else {
					fputs("Split not possible without SSE 4.1\n", stderr);
					exit(0);
				}
				break;
			case 'q':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
#endif
			case 'F':
				HasFilterData = true;
				DBFileName = optarg;
				break;
			case 'H':
				HasHeuristicData = true;
				SeqDB = optarg;
				break;
			case 'V':
				OutputVerbose = true;
				break;
			case 'M':
				modeGiven = true;
				Mode = atoi(optarg);
				break;
			case 'D':
				FilterScoresOutputFile = optarg;
				break;
			case 'd':
				HeuristicScoresOutputFile = optarg;
				break;
			case 's':
				xali1_ptr = xali1_sse2;
				HeuristicVersion = SSE2;
				OtherVersion = SSE2;
				break;
			case 't':
				nCPUs = (int) atoi(optarg);
				break;
			case 'T':
				nCPUsHeuristic = (int) atoi(optarg);
				break;
			case 'C':
				HeuristicCutOff = atoi(optarg);
				break;
			case 'B':
				if (optarg[0] == 'e' || optarg[0] == 'E')
					logarithmic_base = BASE_E;
				else
					logarithmic_base = atof(optarg);
				break;
			case '2':
				binWidth = (float) atof(optarg);
				break;
			case '3':
				binCenter = (float) atof(optarg);
				break;
			case '4':
				avgLength = atof(optarg);
				break;
			case '5': /* Tail area */
				tailArea = atof(optarg);
				break;
			case '7': /* Seed */
				Seed = (long int) atol(optarg);
				break;
			case '8': /* pam distance */
				if (sscanf(optarg,"[%lu,%lu,%lu,%lu]", &PamDistanceStart, &PamDistanceStop, &PamDistanceStep, &PamSampling) != 4) {
					fprintf(stderr,"Failed to read pam distances from %s\n", optarg);
					exit(1);
				}
				break;
			case '9': /* sampling */
			{
				const int itmp = atoi(optarg);
				if (itmp > 0) ProfileSampling = (size_t) itmp;
			}
			break;
			case 'z': /* only-on-match */
				GenOpts = GENERATE_MATCH;
				break;
			case 'Q':
				quantile = atof(optarg);
				if (quantile < 0.0 || quantile > 1.0) {
					fputs("Invalid argument\n", stderr);
					exit(1);
				}
				break;
			case 'E':
				MaxExponent = (float) atof(optarg);
				break;
			case 'm':
			{
				int index=0;
				bool Found = false;
				do {
					if (strncmp(optarg, Methods[index].name, 16) == 0) {
						method = (int) index;
						Found = true;
						break;
					}
				} while (++index<NUMBER_OF_METHODS);
				if (! Found) {
					fputs(" The given method has not been found, please check among the correct possibilities\n", stderr);
					exit(1);
				}
			}
			break;
			case 'o':
				outputFile = optarg;
				break;
			case '+':
				From_perl_script = 1;
				break;
#ifdef __WITH_REPORT__
			case 'r':
				CreatePDFReport = true;
				ReportFileName = optarg;
				break;
#endif
			case 'h':
			default:
				Usage(stdout);
		}
	}

	if (optind >= argc) {
		fputs("Expecting arguments after options\n", stderr);
		Usage(stderr);
	}
	if (optind < argc) ProfileFile = argv[optind];

	if ( method < 0 ) method = 5;

	if (OutputVerbose) {
		fputs(HEADER
#ifdef __USE_MMAP__
		"| Using Linux kernel MMAP function.                                  |\n"
#endif
		,stderr);
		printSystemInfo(stderr, &System);
		if ( HeuristicVersion == SSE2 && (System.Extensions & MM_SSE41)) fputs("Enforcing SSE 2...\n", stderr);
		fprintf(stderr, "Using method %i: %s\n", method, Methods[method].name);

	}

	if (HasFilterData && method < 0 ) {
		fputs("Filter calibration or statictics requires at least a method.\n", stderr);
		exit(1);
	}

	if (From_perl_script) {
		upperProbRange = 0.0;
		lowerProbRange = 0.0;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// INPUT ANALYSIS
	////////////////////////////////////////////////////////////////////////////////////////////////

	/* allocates memory for the profile structure */
	prf = (struct Profile *) _mm_malloc(sizeof(struct Profile), 16);
	if (prf == NULL) {
		fputs("Unable to allocate memory for the profile structure\n", stderr);
		exit(1);
	}

	/*
	 * Read the profile and output some infos
	 */
	gettimeofday(&_t0,0);
	res = ReadProfile(ProfileFile, prf, true);
	gettimeofday(&_t1,0);
	if (res < 0) {
		fputs("Error found reading profile.\n", stderr);
		exit(1);
	}
	if (OutputVerbose) {
		const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
		fprintf(stderr, "Profile reading took %lf seconds.\n", T);

		fprintf(stderr,"Profile %s has length %lu and alphabet size of %lu\n",
						ProfileFile, prf->Length, prf->Alphabet_Length);

		fputs("Alphabet Mapping\n",stderr);
		for (size_t i=0; i<ALPHABET_SIZE; ++i) {
			fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
			if ((i+1) % 8 == 0 ) fputs("\n",stderr);
		}
		fputs("\n\n",stderr);

		fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
	}

	/*
	 * test existance of both level and mode
	 */
	{
		register const int res = SetProfileMode(prf, Mode);
		if (res < 0) {
			fprintf(stderr, "Unable to find cutoff with mode set to %i, error %i\n", Mode, res);
			exit(1);
		}
	}

	/*
	 * Retrieve number of cores
	 */
	nCPUs = (nCPUs == 0) ? (size_t) System.nOverallCores : nCPUs;

#ifdef USE_AFFINITY
  if (noAffinity) {
    // -----------------------------------------------------------------------------
    //                        ***  NO AFFINITY ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose) fputs("Thread affinity disabled\n", stderr);
    Thread_count[0] = System.nOverallCores;
    Thread_masks[0] = (Affinity_Mask_t*) malloc(System.nOverallCores*sizeof(Affinity_Mask_t));
    for (size_t thread=0; thread<System.nOverallCores; ++thread) {
      CPU_ZERO_S(sizeof(Affinity_Mask_t),&Thread_masks[0][thread]);
      for (int i=0; i<(int) System.nOverallCores; ++i)
				CPU_SET_S(i, sizeof(Affinity_Mask_t), (cpu_set_t*) &Thread_masks[0][thread]);
    }
  }
  else if (GivenAffinityFile) {
    // -----------------------------------------------------------------------------
    //                     ***  INPUT FILE HOLDING NASKS ***
    // -----------------------------------------------------------------------------
    if (OutputVerbose)
      fprintf(stderr,"Parsing file %s for affinity mask and number of threads\n", AffinityMaskFileName);
			FILE* in = fopen(AffinityMaskFileName, "r");
			if (in == NULL) {
				fprintf(stderr, "Cannot open thread affinity file %s.\n", optarg);
				exit(1);
			}
    size_t lines = 0;
    while (!feof(in)) {
			int num = fread(buffer, sizeof(char), 64, in);
			for (unsigned int i=0; i<num; i++) if (buffer[i] == '\n') lines++;
		}
		rewind(in);
		if (lines != 0) {
			if (lines > System.nOverallCores) lines = System.nOverallCores;
			Thread_masks[0] = (Affinity_Mask_t*) malloc(lines*sizeof(Affinity_Mask_t));
			for (size_t i=0; i<lines; i++) {
					const int nItems = fscanf(in, "%s\n", buffer);
					const size_t tmp_size = strlen(buffer) - 1;
					CPU_ZERO_S(sizeof(Affinity_Mask_t),&Thread_masks[0][i]);
					for (int j=tmp_size; j>=0; j--) {
						if (buffer[j] != '0') CPU_SET_S(j, sizeof(Affinity_Mask_t), (cpu_set_t*) &Thread_masks[0][i]);
					}
			}
			Thread_count[0] = lines;
			if (OutputVerbose) fprintf(stderr,"Found %2lu threads affinity masks.",nCPUs);
		}
		else {
			if (OutputVerbose) printf("Cannot understand cpu mask, keep on normally\n");
    }
    fclose(in);
  }
  else if (noSharedCore) {
    if (OutputVerbose)
      fputs("No sharing of core resources will be used: Intel Hyperthreading or AMD Compute Unit\n", stderr);
		Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
		if (Thread_count[0] == 0) {
			fputs("No potential affinity mask found !!!\n", stderr);
			exit(0);
		}
  }
  else {
    // -----------------------------------------------------------------------------
    //                        *** OPERATING SYSTEM CHOICE ***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, -1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
  }

  {
    register size_t total = (size_t) (Thread_count[0] + Thread_count[1]);
    if (nCPUs > total) nCPUs = total;
  }

  threads_attr = (pthread_attr_t*) alloca(nCPUs*sizeof(pthread_attr_t));
  {
    register const Affinity_Mask_t * current = &Thread_masks[0][0];
    for (size_t i=0; i<nCPUs; ++i) {
      pthread_attr_init(&threads_attr[i]);
      if (i == (size_t) Thread_count[0]) current = &Thread_masks[1][0];
      pthread_attr_setaffinity_np(&threads_attr[i], sizeof(Affinity_Mask_t), (cpu_set_t*) current);
      ++current;
    }
  }
#endif
	if (OutputVerbose) fprintf(stderr, "Job will be dispatched over %lu cores.\n", nCPUs);

	if (HasHeuristicData && PamSampling>0) {
		/* Correct PamDistanceStart*/
		if (PamDistanceStart == 0) PamDistanceStart = 1;

		/* Compute the number of pam distances to execute */
		size_t count = 0;
		for (size_t i=PamDistanceStart; i<=PamDistanceStop; i+=PamDistanceStep) ++count;
		PamDistanceCounter = count;
		if (OutputVerbose) fprintf(stderr, "Pam distances from %lu to %lu by %lu, overall %lu distances.\n",
			PamDistanceStart, PamDistanceStop, PamDistanceStep, PamDistanceCounter);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// SPECIAL METHOD REQUIRING NO DATABASES
	////////////////////////////////////////////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////
	// COMNON THREAD ALLOCATION
	////////////////////////////////////////////////////////////////////////////////////////////////
	/* Prepare structure common to filter and alignment */
	shares = alloca((nCPUs+1)*sizeof(size_t));

	/* Allocate stack memory for posix thread structures */
	threads_arg = alloca(nCPUs*sizeof(struct ThreadData));
#if !defined(__USE_WINAPI__)
	threads = (pthread_t*) alloca(nCPUs*sizeof(pthread_t));
#else
	threads = (HANDLE*) alloca(nCPUs*sizeof(HANDLE));
#endif

	////////////////////////////////////////////////////////////////////////////////////////////////
	// FILTER
	////////////////////////////////////////////////////////////////////////////////////////////////
	if (HasFilterData) {
		/*
		 * Read the FASTA file
		 */
		gettimeofday(&_t0,0);
		res = AnalyzeFASTAStructure(DBFileName, &DB);
		gettimeofday(&_t1,0);

#ifdef __WITH_REPORT__
		FilterDBSequenceCount = DB.SequenceCount;
		FilterDBSize = (size_t) DB.FileSize;
#endif

		if (OutputVerbose) {
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Sequence file indexing took %lf seconds.\n", T);
		}
		if (res != 0) {
			fputs("Error found.\n", stderr);
			return 1;
		}
		if (OutputVerbose) {
			fprintf(stderr,
							"FASTA file %s analyzed\n"
							"\tFound %lu sequences within %lu bytes\n"
							"\tBiggest sequence entry is %lu bytes\n",
					 DBFileName, DB.SequenceCount, DB.FileSize, DB.MaxSequenceSize);
		}

		/* Allocate memory for the filter scores */
		FilterScores = _mm_malloc(DB.SequenceCount*sizeof(int), 16);
		if (FilterScores == NULL) {
			fputs("Unable to allocate memory for the filter scores\n",stderr);
			exit(1);
		}

		/* Compute the new share for each thread */
		size_t SequenceShare = DB.SequenceCount / nCPUs;
		SequenceShare += (DB.SequenceCount % nCPUs) > (nCPUs-1) ? 1 : 0;
		shares[0] = 0;
		for (size_t i=1; i<nCPUs; ++i) shares[i] = i*SequenceShare;

		shares[nCPUs] = DB.SequenceCount;

#ifdef __USE_MMAP__
		const int fd = open(DBFileName, O_RDONLY );
		const size_t length = (size_t) DB.FileSize;
		const char * const SequenceFileMap = (char *) mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
		if (SequenceFileMap == NULL) {
			fputs("Unable to map sequence file to memory\n", stderr);
			exit(1);
		}
		const char * const restrict SequenceFile = SequenceFileMap;
#else
		const char * const restrict SequenceFile = DB;
#endif

		/* Dispatch to threads */
		{
			gettimeofday(&_t0,0);
			for (size_t i=0; i<nCPUs; ++i) {
				threads_arg[i].prf                       = prf;
				threads_arg[i].DB                        = &DB;
				threads_arg[i].SequenceID                = NULL;
				threads_arg[i].SequenceFile              = SequenceFile;
				threads_arg[i].threadId                  = i;
				threads_arg[i].start                     = shares[i];
				threads_arg[i].stop                      = shares[i+1];
				threads_arg[i].FilterScores              = FilterScores;
				/* dumps real filter score, does not quit if above threshold */
				threads_arg[i].counter                   = 1L;

				if (pthread_create (&threads[i],
#ifdef USE_AFFINITY
					&threads_attr[i],
#else
					NULL,
#endif
					thread_xali1,
					(void*) &threads_arg[i]) != 0)
				{
					return 1;
				}
			}
		}

		for (size_t i=0; i<nCPUs; i++) {
			pthread_join(threads[i], NULL);
		}
		gettimeofday(&_t1,0);

		if ( avgLength <= 0.0 ) {
			unsigned long TotalSequenceLength = 0;
			for (size_t i=0; i<nCPUs; i++) TotalSequenceLength += threads_arg[i].counter;
			resCount = (float)TotalSequenceLength;
			avgLength = resCount/(float)DB.SequenceCount;
			if (OutputVerbose)
				fprintf(stderr, "Database file %s contains %lu residues, average is %lf\n", DBFileName, TotalSequenceLength,avgLength);

		}
		else {
			if (OutputVerbose) fprintf(stderr, "Database average sequence residue length is set to %lf\n", avgLength);
			resCount = (int) ((float)DB.SequenceCount*avgLength + 0.5f);
		}

		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Filter took %lf seconds to treat on %li cores.\n", t, nCPUs);
		}

		/* Call the correct method */
		if (Methods[method].fct) {
			memset(&NewFilterNormItem, 0, sizeof(SNormalizationItem));
			memset(&NewCutItem, 0, sizeof(SCutOffItem));
			/* Sort the score from higher to lower */
			iqsort(FilterScores, (int) DB.SequenceCount);

			if (Methods[method].fct(FilterScores, DB.SequenceCount, (void*) &(NewFilterNormItem.RNOP[0])) >= 0) {
				ProfileChanged = true;
			}
		}
		else {
			switch (method) {
				case 1: /* minimum */
				{
					const union __32bitData res = Statistics(FilterScores, DB.SequenceCount, MINIMUM);
					fprintf(stdout,"Minimum score is %i\n", res.SignedScores);
				}
				break;
				case 2: /* maximum */
				{
					const union __32bitData res = Statistics(FilterScores, DB.SequenceCount, MAXIMUM);
					fprintf(stdout,"Maximum score is %i\n", res.SignedScores);
				}
				break;
				case 3: /* average */
				{
					const union __32bitData res = Statistics(FilterScores, DB.SequenceCount, AVERAGE);
					fprintf(stdout,"Average score is %f\n", res.FloatScores);
				}
				break;
				default:
					fputs("Unknown method after filter performed!\n", stderr);
			}
		}

		if (FilterScoresOutputFile) {
			FILE * const outscore = fopen(FilterScoresOutputFile,"w");
			if (outscore == NULL) {
				fprintf(stderr, "Unable to open file %s to dump scores\n", FilterScoresOutputFile);
			}
			else {
				unsigned int * const restrict TempHeurisiticScores = malloc(DB.SequenceCount*sizeof(unsigned int));
				if (TempHeurisiticScores == NULL) {
					fprintf(stderr, "Unable to allocate memory for heuristic scores\n");
					exit(1);
				}

				TransposeMatrix TIMatch;
				if (HeuristicVersion == SSE41) {
					TIMatch.i = TransposeAndConvertMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
				} else {
					TIMatch.f = TransposeAndConvertToFloatMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
				}

				for (size_t i=0; i<nCPUs; ++i) {
					threads_arg[i].prf              = prf;
					threads_arg[i].DB               = &DB;
					threads_arg[i].SequenceID       = NULL;
					threads_arg[i].SequenceFile     = SequenceFile;
					threads_arg[i].threadId         = i;
					threads_arg[i].start            = shares[i];
					threads_arg[i].stop             = shares[i+1];
					threads_arg[i].HeuristicScores  = TempHeurisiticScores;
					threads_arg[i].version          = HeuristicVersion;
					threads_arg[i].TransposeMatch   = TIMatch;
					if (pthread_create (&threads[i],
	#ifdef USE_AFFINITY
						&threads_attr[i],
	#else
						NULL,
	#endif
						thread_heuristic,
						(void*) &threads_arg[i]) != 0)
					{
						return 1;
					}
				}

				for (size_t i=0; i<nCPUs; i++) {
					pthread_join(threads[i], NULL);
				}

				if (dumpScores(FilterScores, TempHeurisiticScores, SequenceFile, outscore, prf, &(NewFilterNormItem.RNOP[0]), &DB) < 0) {
					fputs("Error within dumpScores\n", stderr);
					exit(1);
				}

				free(TempHeurisiticScores);
				fclose(outscore);
			}
		}



		/* Close the database file */
#ifdef __USE_MMAP__
		munmap((void*)SequenceFileMap, length);
		close(fd);
#endif
		FreeDBStructure(&DB);

		_mm_free(FilterScores);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// HEURISTIC
	////////////////////////////////////////////////////////////////////////////////////////////////
	if ( HasHeuristicData ) {
		/*
		 * Do we need to create on the fly some sequences from the profile
		 * to calibrate the heuristic
		 */
		if (strncmp(SeqDB, "profile", 7) == 0) {
			char * const GeneratedSequencesFile = tempnam("", "GS");
			if (GeneratedSequencesFile != NULL) {
				const int res = GenerateSequences(prf, MaxExponent, GeneratedSequencesFile, Seed,
																					GenOpts, ProfileSampling);
				if (res == 0) {
					SeqDB = GeneratedSequencesFile;
					if (OutputVerbose) fprintf(stderr, "Heuristic DB generated sequences placed in %s\n", GeneratedSequencesFile);
				}
				else {
					fputs("Error within the creation of sequences from the profile\n", stderr);
					switch(res) {
						case (-1):
							perror("GenerateSequences");
							break;
						case (-2):
							fputs("Unable to allocate sufficient memory\n", stderr);
							break;
						case (-3):
							fputs("Unable to extend memory to accomodate generated sequence\n", stderr);
							break;
						case (-4):
							fputs("Impossible to generate profile sequence, too many sequential insertions possible\n", stderr);
							break;
						default:
							fputs("Undefined error code returned\n", stderr);
					}
					exit(1);
				}
			}
			else {
				fputs("Unable to get a temporary file name to store generated sequences from profile\n", stderr);
				exit(1);
			}
			HeuristicSequencesFromProfile = true;
		}


		/* Set initial data */
		memset(&NewHeuristicNormItem, 0, sizeof(SNormalizationItem));

		gettimeofday(&_t0,0);
		res = AnalyzeFASTAStructure(SeqDB, &DB);
		gettimeofday(&_t1,0);
#ifdef __WITH_REPORT__
		HeuristicDBSequenceCount = DB.SequenceCount;
		HeuristicDBSize = (size_t) DB.FileSize;
#endif
		if (OutputVerbose) {
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Heuristic sequence file indexing took %lf seconds.\n", T);
		}
		if (res != 0) {
			fputs("Error found.\n", stderr);
			return 1;
		}
		if (OutputVerbose) {
			fprintf(stderr,
							"FASTA file %s analyzed\n"
							"\tFound %lu sequences within %lu bytes\n"
							"\tBiggest sequence entry is %lu bytes\n",
					 SeqDB, DB.SequenceCount, DB.FileSize, DB.MaxSequenceSize);
		}

		gettimeofday(&_t0,0);
		TransposeMatrix TIMatch;
#ifdef USE_AFFINITY
		TransposeMatrix TIMatch1;
		if (split) {
			TIMatch.i  = TransposeAndConvertMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
			TIMatch1.f = TransposeAndConvertToFloatMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
		} else
#endif
			if (HeuristicVersion == SSE41) {
				TIMatch.i = TransposeAndConvertMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
			} else {
				TIMatch.f = TransposeAndConvertToFloatMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
			}
			gettimeofday(&_t1,0);
		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Transposing Match matrix took %lf seconds.\n", t);
		}

		/* Limit number of threads for heuristic */
		if ( nCPUsHeuristic == 0) nCPUsHeuristic = nCPUs;

		/* Share according to file size */
		{
			size_t FileShare = (size_t) DB.FileSize / nCPUsHeuristic;
			FileShare += ((size_t) DB.FileSize % nCPUsHeuristic) > (nCPUsHeuristic-1) ? 1 : 0;
			const s_Data * DataPtr = DB.DataPtr;
			register size_t counter = 0;
			shares[0] = 0;
			for (size_t i=1; i<nCPUsHeuristic; ++i) {
				register size_t tmp = i*FileShare;
				while ( (size_t) DataPtr->Sequence.Offset < tmp) { ++DataPtr; ++counter; }
				shares[i] = counter;
			}
			shares[nCPUsHeuristic] = DB.SequenceCount;
		}

		const size_t PamMemory = (1 + PamSampling*PamDistanceCounter);
		/* Allocate memory for sequence heuristic scores */
		HeuristicScores = _mm_malloc( DB.SequenceCount*PamMemory*sizeof(unsigned int), 16);
		if (HeuristicScores == NULL) {
			fputs("Unable to allocate memory for weights and Frequencies\n", stderr);
			exit(1);
		}

		/* Allocate memory for the filter scores */
		FilterScores = _mm_malloc(DB.SequenceCount*PamMemory*sizeof(int), 16);
		if (FilterScores == NULL) {
			fputs("Unable to allocate memory for the filter scores\n",stderr);
			exit(1);
		}

		/* Allocate memory for pam types */
		PamTypes = (unsigned char *) malloc(DB.SequenceCount*PamMemory*sizeof(unsigned char));
		if (PamTypes == NULL) {
			fputs("Unable to allocate memory for the pam types\n",stderr);
			exit(1);
		}

#ifdef __USE_MMAP__
		const int fd = open(SeqDB, O_RDONLY );
		const size_t length = (size_t) DB.FileSize;
		const char * const SequenceFileMap = (char *) mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
		if (SequenceFileMap == NULL) {
			fputs("Unable to map sequence file to memory\n", stderr);
			exit(1);
		}
		const char * const restrict SequenceFile = SequenceFileMap;
#else
		const char * const restrict SequenceFile = SeqDB;
#endif

		gettimeofday(&_t0,0);
		for (size_t i=0; i<nCPUsHeuristic; ++i) {
			threads_arg[i].prf             = prf;
			threads_arg[i].DB              = &DB;
			threads_arg[i].SequenceID      = NULL;
			threads_arg[i].SequenceFile    = SequenceFile;
			threads_arg[i].start           = shares[i];
			threads_arg[i].stop            = shares[i+1];
			threads_arg[i].threadId        = i;
			threads_arg[i].TransposeMatch  = TIMatch;
			threads_arg[i].HeuristicScores = HeuristicScores;
			threads_arg[i].FilterScores    = FilterScores;
			threads_arg[i].PamTypes        = PamTypes;
			threads_arg[i].version         = HeuristicVersion;
			threads_arg[i].strand          = FORWARD;
#ifdef USE_AFFINITY
			if (pthread_create (&threads[i],  &threads_attr[i], thread_heuristic_calibration,  (void*) &threads_arg[i]) != 0)
#else
			if (pthread_create (&threads[i],  NULL, thread_heuristic_calibration,  (void*) &threads_arg[i]) != 0)
#endif
			{
				fputs("Fail to create thread.\n", stderr);
				exit(0);
			}
		}

		for (size_t i=0; i<nCPUsHeuristic; ++i) {
			pthread_join(threads[i], NULL);
		}
		gettimeofday(&_t1,0);
		double t;
		if (OutputVerbose) {
			t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Heuristic calibration took %lf seconds to treat on %li cores.\n", t, nCPUsHeuristic);
		}

		/* Free some memory already */
		_mm_free(TIMatch.f);
#ifdef USE_AFFINITY
		if (split) _mm_free(TIMatch1.f);
#endif



		/* Calibrate with quantile regression */
		if (quantile > 0.0 ) {
			int result = QuantileRegression(HeuristicScores, FilterScores, &NewHeuristicNormItem.RNOP[0], quantile,
																			DB.SequenceCount*PamMemory);
			if (result != 0)
				fputs("quantile regression error\n",stderr);
			else {
				if (OutputVerbose)
					fprintf(stderr, "Heuristic %.1f%% quantile regression : Rx=(%lf, %lf)\n", 100.0*quantile, NewHeuristicNormItem.RNOP[0], NewHeuristicNormItem.RNOP[1]);
				ProfileChanged = true;
			}
		}
		else
			ProfileChanged = true;

		/* Now output the data */
		if (HeuristicScoresOutputFile) {
			FILE * const outscore = fopen(HeuristicScoresOutputFile,"w");
			if (outscore == NULL) {
				fprintf(stderr, "Unable to open file %s to dump heuristic - filter scores\n", HeuristicScoresOutputFile);
			}
			else {
				/* Allocate memory to hold sequence */
				SeqData.Data.Memory = malloc(DB.MaxSequenceSize*sizeof(unsigned char));
				if (SeqData.Data.Memory == NULL) {
					fputs("Thread Cannot allocate memory for sequence.\n", stderr);
					return 1;
				}

				/* Open sequence file*/
				SETUP_DATABASE_ACCESS(SequenceFile);
				fputs("# Sequence Normalized_score Filter_score Heuristic_score\n", outscore);

				for (size_t d=0; d<DB.SequenceCount; d++) {
					PFSequence * const PFSeq = GET_DATABASE_SEQUENCE(&SeqData, &(DB.DataPtr[d].Sequence));
					for (size_t e=0; e<PamMemory; ++e) {
						const float NormalizedScore = prf->RawToNormalized(FilterScores[PamMemory*d+e], &(NewFilterNormItem.RNOP[0]), 0, 0);
						fprintf(outscore, "%s_%u\t%f\t%i\t%u\n", &SeqData.Data.Header[1], (unsigned int) PamTypes[PamMemory*d+e], NormalizedScore, FilterScores[PamMemory*d+e], HeuristicScores[PamMemory*d+e]);
					}
				}
				/* close sequence file */
				UNSET_DATABASE_ACCESS();
				free(SeqData.Data.Memory);

				fclose(outscore);
			}
		}

		/* Close the database file */
#ifdef __USE_MMAP__
		munmap((void*)SequenceFileMap, length);
		close(fd);
#endif
		FreeDBStructure(&DB);

		if (HeuristicSequencesFromProfile) {
			if (unlink(SeqDB) != 0) {
				perror("Error in removing generated sequence file");
			}
			else {
				if (OutputVerbose) fprintf(stderr, "Cleaning temporary file %s\n", SeqDB);
			}
		}

#ifdef __WITH_REPORT__
		/* Generate plot to be included in report */
		if (CreatePDFReport) {
			HeuristicTempName = tempnam(0,"H");
			if (! HeuristicTempName) {
				fputs("Unable to create a temporary file name\n",stderr);
			}
			else {
				size_t SequenceCount = (size_t) (DB.SequenceCount*PamMemory);
#ifndef PL_DOUBLE
				/* Transform Heuristic and Filter scores into float and get the maximum, minimum */
				float * const X = (float) &FilterScores[0];
				float * const Y = (float) &HeuristicScores[0];
				size_t i = 4;
				__m128 __maxX = _mm_set1_ps(FLT_MIN);
				__m128 __minX = _mm_set1_ps(FLT_MAX);
				__m128 __maxY = __maxX;

				while (i < SequenceCount) {
					__m128 __fH = _mm_cvtepi32_ps(*(__m128i*) &HeuristicScores[i-4]);
					__m128 __fF = _mm_cvtepi32_ps(*(__m128i*) &FilterScores[i-4]);
					i += 4;
					__maxX = _mm_max_ps(__maxX, __fF);
					__minX = _mm_min_ps(__minX, __fF);
					__maxY = _mm_max_ps(__maxY, __fH);
					_mm_store_ps((float*) &HeuristicScores[i-8], __fH);
					_mm_store_ps((float*) &FilterScores[i-8], __fF);
				}
				i -= 4;
				__m128 __zero = _mm_setzero_ps();
				while (i < SequenceCount) {
					__m128 __fH = _mm_cvtsi32_ss(__zero, HeuristicScores[i]);
					__m128 __fF = _mm_cvtsi32_ss(__zero, FilterScores[i]);
					i += 1;
					__maxX = _mm_max_ss(__maxX, __fF);
					__minX = _mm_min_ss(__minX, __fF);
					__maxY = _mm_max_ss(__maxY, __fH);
					_mm_store_ss((float*) &HeuristicScores[i-1], __fH);
					_mm_store_ss((float*) &FilterScores[i-1], __fF);
				}
#else
				double * const X = (double*) _mm_malloc(SequenceCount*sizeof(double),16);
				double * const Y = (double*) _mm_malloc(SequenceCount*sizeof(double),16);
				if (X == NULL || Y == NULL) {
					fputs("Cannot allocate memory for graph.\n", stderr);
					_mm_free(FilterScores);
					_mm_free(HeuristicScores);
					exit(1);
				}
				size_t i = 4;
				__m128 __maxX = _mm_set1_ps(FLT_MIN);
				__m128 __minX = _mm_set1_ps(FLT_MAX);
				__m128 __maxY = __maxX;
				// Crappy code SSE 4.1 and AVX would be much better !!!
				while (i < SequenceCount) {
					__m128 __fH = _mm_cvtepi32_ps(*(__m128i*) &HeuristicScores[i-4]);
					__m128 __fF = _mm_cvtepi32_ps(*(__m128i*) &FilterScores[i-4]);
					i += 4;
					__maxX = _mm_max_ps(__maxX, __fF);
					__minX = _mm_min_ps(__minX, __fF);
					__maxY = _mm_max_ps(__maxY, __fH);
					__m128d __X1 = _mm_cvtps_pd(__fF);
					__m128d __Y1 = _mm_cvtps_pd(__fH);
					__m128d __X2 = _mm_cvtps_pd(_mm_movehl_ps(__fF, __fF));
					__m128d __Y2 = _mm_cvtps_pd(_mm_movehl_ps(__fH, __fH));

					_mm_store_pd(&X[i-8], __X1);
					_mm_store_pd(&X[i-6], __X2);
					_mm_store_pd(&Y[i-8], __Y1);
					_mm_store_pd(&Y[i-6], __Y2);
				}
				i -= 4;
				__m128 __zero = _mm_setzero_ps();
				while (i < SequenceCount) {
					__m128 __fH = _mm_cvtsi32_ss(__zero, HeuristicScores[i]);
					__m128 __fF = _mm_cvtsi32_ss(__zero, FilterScores[i]);
					i += 1;
					__maxX = _mm_max_ss(__maxX, __fF);
					__minX = _mm_min_ss(__minX, __fF);
					__maxY = _mm_max_ss(__maxY, __fH);
					__m128d __X1 = _mm_cvtss_sd(_mm_castps_pd(__zero), __fF);
					__m128d __Y1 = _mm_cvtss_sd(_mm_castps_pd(__zero), __fH);

					_mm_store_sd( &X[i-1], __X1);
					_mm_store_sd( &Y[i-1], __Y1);
				}
#endif

				// movhlps xmm1, xmm0	Move top two floats to lower part of xmm1
				__m128 __tMaxX = _mm_movehl_ps(__maxX, __maxX);
				__m128 __tMaxY = _mm_movehl_ps(__maxY, __maxY);
				__m128 __tMinX = _mm_movehl_ps(__minX, __minX);
				// maxps   xmm0,xmm1	Get maximum/minimum of the two sets of floats
				__maxX = _mm_max_ps(__maxX, __tMaxX);
				__maxY = _mm_max_ps(__maxY, __tMaxY);
				__minX = _mm_min_ps(__minX, __tMinX);
				// shufps  xmm1,xmm0,$55	Move second float to lower part of xmm1 0101 0101 = 55
				__tMaxX = _mm_shuffle_ps(__maxX, __maxX, 0b01010101);
				__tMinX = _mm_shuffle_ps(__minX, __minX, 0b01010101);
				__tMaxY = _mm_shuffle_ps(__maxY, __maxY, 0b01010101);
				// maxps   xmm0,xmm1	Get minimum of the two remaining floats
				__maxX = _mm_max_ps(__maxX, __tMaxX);
				__maxY = _mm_max_ps(__maxY, __tMaxY);
				__minX = _mm_min_ps(__minX, __tMinX);

				float maxX, maxY, minX;
				_mm_store_ss(&maxX, __maxX);
				_mm_store_ss(&maxY, __maxY);
				_mm_store_ss(&minX, __minX);

				PLFLT x_max = (PLFLT) maxX;
				PLFLT x_min = (PLFLT) minX;
				PLFLT y_max = (PLFLT) (1.1*maxY);
				PLFLT y_min = (PLFLT) 0.0;

				plsdev("jpeg");
				plsfnam(HeuristicTempName);
				plscmap0(red, green, blue, 10);
				plspage(0,0,1024,800,0,0);

				plinit();
				plfont(2);

				x_min -= 10;
				x_max += 10;
				y_min -= 10;
				y_max += 10;

				plenv(x_min, x_max, y_min, y_max, 0, 0);
				pllab("Filter scores", "Heuristic scores", "");
				plcol0(2);

				// 	plcol0(5);plssym(0,.75);
				// 	plpoin((PLINT) SequenceCount, X, Y, 17);
				plscmap1n((PLINT) (1+PamDistanceCounter));
				const PLFLT invMaxType = (PLFLT) (1.0f/(float)(PamDistanceStop));
				const PLFLT Red[]   = {1, 0, 0, 0, 1};
				const PLFLT Green[] = {0, 0, 0, 1, 1};
				const PLFLT Blue[]  = {0, 1, 1, 0, 0};
				const PLFLT Intensity[] = {0,0.1,0.9*invMaxType, 1-0.45*invMaxType, 1};

				plscmap1l ((PLBOOL) 1, 5, Intensity, Red, Green, Blue, NULL);
				plssym(0,.75);

				for (size_t ival=0; ival<DB.SequenceCount*PamMemory; ival++) {
					PLFLT color;
					if (PamTypes[ival] != 0) {
						color = (PLFLT) 0.3 + (((float)PamTypes[ival])*invMaxType);
						color = color > 0.99999 ? 0.99999 : color;
					}
					else
						color = (PLFLT) 0.0;

					plcol1(color);
					plpoin((PLINT) 1, &X[ival], &Y[ival], 17);
				}

				plcol0(3);
				const PLFLT ymin = (PLFLT) (NewHeuristicNormItem.RNOP[0] + x_min * NewHeuristicNormItem.RNOP[1]);
				const PLFLT ymax = (PLFLT) (NewHeuristicNormItem.RNOP[0] + x_max * NewHeuristicNormItem.RNOP[1]);
				pljoin(x_min, ymin, x_max, ymax);

				plend();
#ifdef PL_DOUBLE
				_mm_free(X);
				_mm_free(Y);
#endif
			}
		}
#endif

		/* Free last memory */
		_mm_free(FilterScores);
		_mm_free(HeuristicScores);
		free(PamTypes);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// APPLY PROFILE MODIFICATION
	////////////////////////////////////////////////////////////////////////////////////////////////
	if (ProfileChanged) {
		/*
		 * Modifications in the filter calibration
		 */
		if (HasFilterData) {
			if (OutputVerbose) {
				fprintf(stderr,"New normalization coefficients from method %s are Rx=(%lf, %lf)\n",
								Methods[method].name, NewFilterNormItem.RNOP[0], NewFilterNormItem.RNOP[1]);
			}

			/* Add remaining information to the normalization item */
			memcpy(NewFilterNormItem.CNTX, Methods[method].NormalizationText, 32*sizeof(char)); // setup TEXT
			NewFilterNormItem.MNOR = 0; 		// setup FUNCTION to LINEAR, refer to io.c
			NewFilterNormItem.NNOR = Mode;	// setup MODE
			NewFilterNormItem.NNPR = 0;			// setup PRIORITY to 0

			/* Update the normalization item */
			memcpy(&prf->NormalizationData.Values[prf->ModeIndex], &NewFilterNormItem, sizeof(SNormalizationItem));
		}

		/*
		 * Modifications in the heuristic calibration
		 */
		int HeuristicModeIndex;
		if (HasHeuristicData) {
			/* Test if the heuristic mode already exists */
			HeuristicModeIndex = SeekProfileMode(prf, -Mode);
			if (HeuristicModeIndex < 0) {
				if (OutputVerbose) {
					fprintf(stderr,"Creation of heuristic normalization mode data using id %i\n", -Mode);
				}
				/* Do we have enough space left */
				if ( prf->NormalizationData.JNOR < MAXN) {
					HeuristicModeIndex = prf->NormalizationData.JNOR++;
				}
				else {
					fputs("Not enough space left to store a new normalization mode, skipping.\n", stderr);
					goto NO_WRITE;
				}
				snprintf(NewHeuristicNormItem.CNTX, 32, "Heuristic %.1f%%", 100.0*quantile);	// description as read in profile TEXT
				NewHeuristicNormItem.MNOR = 0; 				// setup FUNCTION to LINEAR, refer to io.c
				NewHeuristicNormItem.NNOR = -Mode;		// Mode as read in profile MODE
				NewHeuristicNormItem.NNPR = 0;				// Priority as read in profile PRIORITY
			}
			else {
				const SNormalizationItem * const Modeptr = &prf->NormalizationData.Values[HeuristicModeIndex];
				strncpy(NewHeuristicNormItem.CNTX, Modeptr->CNTX, 32); 	// setup TEXT
				NewHeuristicNormItem.MNOR = 0; 			     	 // setup FUNCTION to LINEAR, refer to io.c
				NewHeuristicNormItem.NNOR = Modeptr->NNOR; // setup MODE
				NewHeuristicNormItem.NNPR = Modeptr->NNPR; // setup PRIORITY to 0
			}

			/* Update the normalization item */
			memcpy(&prf->NormalizationData.Values[HeuristicModeIndex], &NewHeuristicNormItem, sizeof(SNormalizationItem));

			NO_WRITE:
			if (OutputVerbose) {
				fprintf(stderr,"New heuristic normalization coefficients are Rx=(%lf, %lf)\n",
								NewHeuristicNormItem.RNOP[0], NewHeuristicNormItem.RNOP[1]);
			}

		}

		/* Update the all the cutoff item using that mode */
		/* WARNING: ASSUMPTION ON THE FACT THAT THERE IS ONLY 1 MODE PER CUTOFF LEVEL */
		const int cutCount = prf->CutOffData.JCUT;
		for( int icut=0; icut<cutCount; icut++) {
			SCutOffItem * const cutitem = &prf->CutOffData.Values[icut];
			if (cutitem->MCUT[0] == Mode) {
				const int RawCutOff = prf->NormalizedToRaw(cutitem->RCUT[0], prf->NormalizationData.Values[prf->ModeIndex].RNOP, 0.0, 0);
				cutitem->ICUT = RawCutOff;
				if (HasHeuristicData) {
					const float ftmp =  ((float) RawCutOff) * NewHeuristicNormItem.RNOP[1] + NewHeuristicNormItem.RNOP[0] + 0.5f;
					if ( ftmp <= 0.0f ) {
						cutitem->HCUT = 0U;
						fprintf(stderr, "Heuristic setting not possible for normalized cutoff %lf, raw score %i.\n", cutitem->RCUT[0], RawCutOff);
						exit(1);
					}
					else {
						cutitem->HCUT = (unsigned int) ftmp;
					}
				}
			}
		}

		/* Write the new profile */
		char line[512];
		snprintf(line, 512, "CC   %s", argv[0]);
		for (int i=1; i<argc-1; ++i) {
			char * pos = &line[strlen(line)];
			sprintf(pos, " %s", argv[i]);
		}
		{
			char * pos = &line[strlen(line)];
			sprintf(pos, " %s\n", argv[argc-1]);
		}

		if (outputFile != NULL) {
			FILE * const out = fopen(outputFile,"w");
			if (out == NULL) {
				fprintf(stderr,"Unable to open %s for writing.\n", outputFile);
				exit(1);
			}
			WriteProfile(ProfileFile, prf, line, out);
			fclose(out);
		}
		else {
			WriteProfile(ProfileFile, prf, line, stdout);
		}

		////////////////////////////////////////////////////////////////////////////////////////////////
		// CREATE PDF REPORT
		////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef __WITH_REPORT__
#define BORDER_SIZE 10.0f
#define RIGHT_MARGIN 20.0f
#define LEFT_MARGIN 20.0f
#define TOP_MARGIN 20.0f
#define BOTTOM_MARGIN 20.0f
		if (CreatePDFReport) {
			HPDF_Doc  pdf;
			HPDF_Font font;
			HPDF_Page page;
			char fname[256] __attribute__((aligned(16)));

			/* Create file name */
			{
				const size_t NameLength = strlen(ReportFileName);
				if (ReportFileName[NameLength-1] == 'f' &&
						ReportFileName[NameLength-2] == 'd' &&
						ReportFileName[NameLength-3] == 'p' &&
						ReportFileName[NameLength-4] == '.') {
						ReportFileName[NameLength-4] = '\0';
					}
					snprintf(fname, 255, "%s.pdf", ReportFileName);
			}
			if (OutputVerbose) fprintf(stderr, "Creating PDF report in %s\n", fname);

			pdf = HPDF_New (error_handler, NULL);
			if (!pdf) {
				printf ("Error: cannot create PdfDoc object\n");
				exit(1);
			}

			if (setjmp(env)) {
				HPDF_Free (pdf);
				exit(1);
			}

			/* create default-font */
			font = HPDF_GetFont (pdf, "Helvetica", NULL);

			/* add a new page object. */
			page = HPDF_AddPage (pdf);

			/* set to A4 format */
			HPDF_Page_SetSize(page, HPDF_PAGE_SIZE_A4, HPDF_PAGE_PORTRAIT);
			const HPDF_REAL PageWidth  = HPDF_Page_GetWidth(page) - 2.0f*BORDER_SIZE - LEFT_MARGIN - RIGHT_MARGIN;
			const HPDF_REAL PageHeight = HPDF_Page_GetHeight(page) - 2.0f*BORDER_SIZE - TOP_MARGIN - BOTTOM_MARGIN;

			/* print the lines of the page. */
			HPDF_Page_SetLineWidth (page, 0.0);

			/*
			 * Header -------------------------------------------------------------------------
			 */
			const char * cptr;
			char Text[128];
			HPDF_Page_SetGrayFill(page, 0.8);
			HPDF_Page_Rectangle(page,
													BORDER_SIZE+LEFT_MARGIN, BORDER_SIZE+BOTTOM_MARGIN + PageHeight-20,
											 PageWidth, 20 );
			HPDF_Page_Fill(page);
			HPDF_Page_SetGrayFill(page, 0);

			cptr = "Profile data";
			HPDF_Page_SetFontAndSize (page, font, 16);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10, BORDER_SIZE+BOTTOM_MARGIN + PageHeight - 16);
			HPDF_Page_ShowText (page, cptr);
			HPDF_Page_EndText (page);

			/*HPDF_Page_TextRect  (HPDF_Page            page,
			 *                     HPDF_REAL            left,
			 *                     HPDF_REAL            top,
			 *                     HPDF_REAL            right,
			 *                     HPDF_REAL            bottom,
			 *                     const char          *text,
			 *                     HPDF_TextAlignment   align,
			 *                     HPDF_UINT           *len);   */
#define WRITE_ALIGN(TEXT1,TEXT2, ROW) {\
		const HPDF_REAL center_x = BORDER_SIZE+LEFT_MARGIN + 10 + center;\
		const HPDF_REAL center_y = BORDER_SIZE+BOTTOM_MARGIN + PageHeight - 20 -5 - ROW*(0 + 10);\
		const HPDF_REAL tw = HPDF_Page_TextWidth (page, TEXT1);\
		HPDF_Page_BeginText (page);\
		HPDF_Page_MoveTextPos (page, center_x-5-tw, center_y);\
		HPDF_Page_ShowText (page, TEXT1);\
		HPDF_Page_EndText (page);\
		HPDF_Page_BeginText (page);\
		HPDF_Page_MoveTextPos (page, center_x+5, center_y);\
		HPDF_Page_ShowText (page, TEXT2);\
		HPDF_Page_EndText (page);\
		}
		HPDF_Page_SetFontAndSize (page, font, 10);
		HPDF_REAL center = 65;
		WRITE_ALIGN("Identification", prf->Identification, 1);
		WRITE_ALIGN("Description", prf->Description,2);
		WRITE_ALIGN("AC number", prf->AC_Number,3);
		WRITE_ALIGN("Date",prf->Date,4);
		WRITE_ALIGN("File name", ProfileFile,5);
		snprintf(Text,128,"%lu", prf->Length);
		WRITE_ALIGN("Length", Text,6);
		size_t i=1;
		while (prf->CABC[i] != 'X' && i<ALPHABET_MEMORY_SIZE) Text[i-1] = prf->CABC[i++];
		Text[i-1] = '\0';
		WRITE_ALIGN("Alphabet", Text, 7);
		cptr = (prf->isCircular) ? "Yes" : "No";
		WRITE_ALIGN("Circular", cptr,8);
		HPDF_REAL lastYpos = BORDER_SIZE+BOTTOM_MARGIN + PageHeight - 20 -5 - 8*(0 + 10);
		{
			const HPDF_REAL ROW = 9;
			const HPDF_REAL center_x = BORDER_SIZE+LEFT_MARGIN + 10 + center;
			HPDF_REAL center_y = BORDER_SIZE+BOTTOM_MARGIN + PageHeight - 20 -5 - ROW*(0 + 10);
			const HPDF_REAL tw = HPDF_Page_TextWidth (page, "Consensus");
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, center_x-5-tw, center_y);
			HPDF_Page_ShowText (page, "Concensus");
			HPDF_Page_EndText (page);
			HPDF_UINT len = 0;
			const HPDF_UINT Total = strlen(prf->Sequence);
			HPDF_Page_SetFontAndSize (page, font, 8);


			while (len < Total) {
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, center_x+5, center_y);
				snprintf(Text,128,"%.*s", 80, &prf->Sequence[len]);
				HPDF_Page_ShowText (page, Text);
				center_y -= 10;
				len += 80;
				lastYpos -= 10;
				HPDF_Page_EndText (page);
			}

		}
#undef WRITE_ALIGN

		lastYpos -= 5; //BORDER_SIZE+BOTTOM_MARGIN + PageHeight - 20 -5 - 9*(0 + 10);
		HPDF_Page_SetGrayFill(page, 0.8);
		HPDF_Page_Rectangle(page,
												BORDER_SIZE+LEFT_MARGIN+10, lastYpos - 5 - 15,
											PageWidth-10, 15 );
		HPDF_Page_Fill(page);
		HPDF_Page_SetGrayFill(page, 0);
		cptr = "Normalization modes";
		HPDF_Page_SetFontAndSize (page, font, 12);
		HPDF_Page_BeginText (page);
		HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10 + 10, lastYpos - 5 - 12);
		HPDF_Page_ShowText (page, cptr);
		HPDF_Page_EndText (page);

		const char (*CNOR)[16] = (const char (*)[16]) prf->NormalizationData.CNOR;
		center = 75;
		lastYpos -= 5 + 15 + 10;
		HPDF_Page_SetFontAndSize (page, font, 10);
		for (int iNorm=0;iNorm<prf->NormalizationData.JNOR; ++iNorm) {
			const SNormalizationItem * const normItem = &(prf->NormalizationData.Values[iNorm]);
			const HPDF_REAL center_x = BORDER_SIZE+LEFT_MARGIN + 10 + center;
			const HPDF_REAL center_y = lastYpos - 5 - iNorm*(0 + 10);
			snprintf(Text, 128,"Mode %i", normItem->NNOR);
			const HPDF_REAL tw = HPDF_Page_TextWidth (page, Text);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, center_x-5-tw, center_y);
			HPDF_Page_ShowText(page, Text);
			HPDF_Page_EndText(page);
			HPDF_Page_BeginText (page);
			snprintf(Text,128, "function: \'%s\' text: \'%s\' priority: %i", CNOR[normItem->MNOR],
							 normItem->CNTX, normItem->NNPR);
			for (int i=0; i<prf->NormalizationData.JNOP[normItem->MNOR]; ++i) {
				const size_t l = strlen(Text);
				snprintf(&Text[l], 128-l, " R%1.1i: %lf", 1+i, normItem->RNOP[i]);
			}

			HPDF_Page_MoveTextPos (page, center_x+5, center_y);
			HPDF_Page_ShowText(page, Text);
			HPDF_Page_EndText(page);
		}

		lastYpos -=  prf->NormalizationData.JNOR*(0 + 10);HPDF_Page_SetGrayFill(page, 0.8);
		HPDF_Page_Rectangle(page,
												BORDER_SIZE+LEFT_MARGIN+10, lastYpos - 5 - 15,
											PageWidth-10, 15 );
		HPDF_Page_Fill(page);HPDF_Page_SetGrayFill(page, 0.);
		cptr = "Cutoff modes";
		HPDF_Page_SetFontAndSize (page, font, 12);
		HPDF_Page_BeginText (page);
		HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10 + 10, lastYpos - 5 - 12);
		HPDF_Page_ShowText (page, cptr);
		HPDF_Page_EndText (page);

		lastYpos -= 5 + 15 + 10;
		HPDF_Page_SetFontAndSize (page, font, 10);
		for (int iCut=0;iCut<prf->CutOffData.JCUT; ++iCut) {
			const SCutOffItem * const cutItem = &(prf->CutOffData.Values[iCut]);
			const HPDF_REAL center_x = BORDER_SIZE+LEFT_MARGIN + 10 + center;
			const HPDF_REAL center_y = lastYpos - 5 - iCut*(0 + 10);
			snprintf(Text, 128,"Level %i", cutItem->MCLE);
			const HPDF_REAL tw = HPDF_Page_TextWidth (page, Text);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, center_x-5-tw, center_y);
			HPDF_Page_ShowText(page, Text);
			HPDF_Page_EndText(page);
			HPDF_Page_BeginText (page);
			if (cutItem->HCUT > 0) {
				snprintf(Text,128, "Mode %i score: %i normalized score: %.1f heuristic score: %i text: \'%s\'",
								 cutItem->MCUT[0],cutItem->ICUT, cutItem->RCUT[0], cutItem->HCUT, cutItem->CCUT);
			}
			else {
				snprintf(Text,128, "Mode %i score: %i normalized score: %.1f text: \'%s\'",
								 cutItem->MCUT[0], cutItem->ICUT, cutItem->RCUT[0], cutItem->CCUT);
			}
			HPDF_Page_MoveTextPos (page, center_x+5, center_y);
			HPDF_Page_ShowText(page, Text);
			HPDF_Page_EndText(page);
		}

		lastYpos -= 5 + prf->CutOffData.JCUT*(0 + 10);
		HPDF_Page_SetGrayFill(page, 0.8);
		HPDF_Page_Rectangle(page,
												BORDER_SIZE+LEFT_MARGIN, lastYpos-20,
												PageWidth, 20 );
		HPDF_Page_Fill(page);HPDF_Page_SetGrayFill(page, 0.);

		snprintf(Text, 128, "Calibration of mode %i", Mode);
		HPDF_Page_SetFontAndSize (page, font, 16);
		HPDF_Page_BeginText (page);
		HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10, lastYpos - 16);
		HPDF_Page_ShowText (page, Text);
		HPDF_Page_EndText (page);
		lastYpos -= 5 + 20;
		//       HPDF_REAL tw = HPDF_Page_TextWidth (page, cptr);
		//       HPDF_Page_BeginText (page);
		//       HPDF_Page_MoveTextPos (page, (HPDF_Page_GetWidth(page) - tw) / 2,
		// 		  HPDF_Page_GetHeight (page) - BORDER_SIZE - TOP_MARGIN);
		//       HPDF_Page_ShowText (page, cptr);
		//       HPDF_Page_EndText (page);

		//       HPDF_Page_BeginText (page);
		//       HPDF_Page_MoveTextPos (page, BORDER_SIZE, HPDF_Page_GetHeight (page) - BORDER_SIZE);
		//       HPDF_Page_ShowText (page, prf->Sequence);
		//       HPDF_Page_EndText (page);
		/*
		 * Score calibration --------------------------------------------------------------
		 */
		const HPDF_REAL Graph_width = 0.5; /* ration related to the page width */
		if (FilterTempName) {
			HPDF_Page_SetGrayFill(page, 0.8);
			HPDF_Page_Rectangle(page,
													BORDER_SIZE+LEFT_MARGIN+10, lastYpos - 5 - 15,
											 PageWidth-10, 15 );
			HPDF_Page_Fill(page);HPDF_Page_SetGrayFill(page, 0.);
			cptr = "Standard score cutoff";
			HPDF_Page_SetFontAndSize (page, font, 12);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10 + 10, lastYpos - 5 - 12);
			HPDF_Page_ShowText (page, cptr);
			HPDF_Page_EndText (page);

			lastYpos -= 5 + 15;
#define WRITE_ALIGN(TEXT1,TEXT2, ROW) {\
		const HPDF_REAL center_x = BORDER_SIZE+LEFT_MARGIN + 10 + center;\
		const HPDF_REAL center_y = lastYpos -35 - (ROW)*(0 + 10);\
		const HPDF_REAL tw = HPDF_Page_TextWidth (page, TEXT1);\
		HPDF_Page_BeginText (page);\
		HPDF_Page_MoveTextPos (page, center_x-5-tw, center_y);\
		HPDF_Page_ShowText (page, TEXT1);\
		HPDF_Page_EndText (page);\
		HPDF_Page_BeginText (page);\
		HPDF_Page_MoveTextPos (page, center_x+5, center_y);\
		HPDF_Page_ShowText (page, TEXT2);\
		HPDF_Page_EndText (page);\
	}

		HPDF_Page_SetFontAndSize (page, font, 10);
		center = 75;
		{
			char * cptr;
			for (size_t i=0; i<strlen(DBFileName);++i) if (DBFileName[i] == '/') cptr = &DBFileName[i+1];
			WRITE_ALIGN("Database", cptr, 1);
		}
		snprintf(Text, 128, "%'lu bytes\n", FilterDBSize);
		WRITE_ALIGN("Database size", Text, 2);
		snprintf(Text, 128, "%'lu", FilterDBSequenceCount);
		WRITE_ALIGN("# sequences", Text, 3);
		snprintf(Text,128, "%'lu", (unsigned long) resCount);
		WRITE_ALIGN("# residues", Text,4);
		snprintf(Text, 128, "%lf", avgLength);
		WRITE_ALIGN("Mean length", Text, 5);
		snprintf(Text, 128, "%lf",logarithmic_base);
		WRITE_ALIGN("Logarithmic base", Text, 6);
		WRITE_ALIGN("Method", Methods[method].name, 7);
		if (method == 5) /* pfscale */ {
			snprintf(Text, 128, "%.2lf", first_rank);
			WRITE_ALIGN("First rank", Text, 8);
			snprintf(Text, 128, "%.2lf", last_rank);
			WRITE_ALIGN("Last rank", Text, 9);
			snprintf(Text, 128, "%lg", upperProbRange);
			WRITE_ALIGN("Upper probability", Text, 10);
			snprintf(Text, 128, "%lg", lowerProbRange);
			WRITE_ALIGN("Lower probability", Text, 11);

			const SNormalizationItem * const normItem = &(prf->NormalizationData.Values[prf->ModeIndex]);
			WRITE_ALIGN("Function", CNOR[normItem->MNOR], 13);
			char template[] = "R1";
			for (int i=0; i<prf->NormalizationData.JNOP[normItem->MNOR]; ++i) {
				template[1] += (unsigned char) i;
				snprintf(Text, 128, "%lf", normItem->RNOP[i]);
				WRITE_ALIGN(template, Text, 14+i);
			}
		}

		HPDF_Image JPGScores = HPDF_LoadJpegImageFromFile(pdf, FilterTempName);
		HPDF_REAL iw = HPDF_Image_GetWidth (JPGScores);
		HPDF_REAL ih = HPDF_Image_GetHeight (JPGScores);
		if (iw > Graph_width*PageWidth) {
			ih *= Graph_width*PageWidth/iw;
			iw = Graph_width*PageWidth;
		}

		HPDF_Page_DrawImage (page, JPGScores, BORDER_SIZE+LEFT_MARGIN+(1.- Graph_width)*PageWidth, lastYpos - ih - 5, iw, ih);
		lastYpos -= ih;
		}
		/*
		 * Heuristic score calibration ----------------------------------------------------
		 */
		if (HeuristicTempName) {
			HPDF_Image JPGScores = HPDF_LoadJpegImageFromFile(pdf, HeuristicTempName);
			HPDF_REAL iw = HPDF_Image_GetWidth (JPGScores);
			HPDF_REAL ih = HPDF_Image_GetHeight (JPGScores);
			if (iw > Graph_width*PageWidth) {
				ih *= Graph_width*PageWidth/iw;
				iw = Graph_width*PageWidth;
			}
			if (lastYpos < ih+BORDER_SIZE+BOTTOM_MARGIN) {
				HPDF_Page_SetFontAndSize (page, font, 6);
				HPDF_Page_BeginText (page);
				HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN, BORDER_SIZE+BOTTOM_MARGIN - 5 - 12);
				HPDF_Page_ShowText (page, "(C) SIB Swiss Institute of Bioinformatics, PFCalibrate V " PF_VERSION);
				HPDF_Page_EndText (page);
				{
					HPDF_Page_BeginText (page);
					time_t t = time(NULL);
					snprintf(Text, 128, "%s", ctime(&t));
					const HPDF_REAL tw = HPDF_Page_TextWidth (page, Text);\
					HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN+PageWidth-tw, BORDER_SIZE+BOTTOM_MARGIN - 5 - 12);
					HPDF_Page_ShowText (page, Text);
					HPDF_Page_EndText (page);
				}
				page = HPDF_AddPage(pdf);
				lastYpos = BORDER_SIZE+BOTTOM_MARGIN+PageHeight;
			}
			// 	lastYpos -= 5 + 20;
			HPDF_Page_SetGrayFill(page, 0.8);
			HPDF_Page_Rectangle(page,
													BORDER_SIZE+LEFT_MARGIN+10, lastYpos - 5 - 15,
											 PageWidth-10, 15 );
			HPDF_Page_Fill(page);HPDF_Page_SetGrayFill(page, 0.);
			cptr = "Heuristic score cutoff";
			HPDF_Page_SetFontAndSize (page, font, 12);
			HPDF_Page_BeginText (page);
			HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN + 10 + 10, lastYpos - 5 - 12);
			HPDF_Page_ShowText (page, cptr);
			HPDF_Page_EndText (page);

			lastYpos -= 5 + 15;
			HPDF_Page_SetFontAndSize (page, font, 10);
			int ROW;
			if (!HeuristicSequencesFromProfile) {
				char * cptr = SeqDB;
				for (size_t i=0; i<strlen(SeqDB);++i) if (SeqDB[i] == '/') cptr = &SeqDB[i+1];
				WRITE_ALIGN("Database", cptr, 1);
				snprintf(Text, 128, "%'lu bytes\n", HeuristicDBSize);
				WRITE_ALIGN("Database size", Text, 2);
				ROW = 3;
			}
			else {
				WRITE_ALIGN("Database", "generated based on the profile", 1);
				ROW = 2;
			}
			snprintf(Text, 128, "%'lu", HeuristicDBSequenceCount);
			WRITE_ALIGN("# sequences", Text, ROW);
			snprintf(Text, 128, "%lu", Seed);
			WRITE_ALIGN("Generator seed", Text, ROW+1);
			snprintf(Text, 128, "%lu", nCPUsHeuristic);
			WRITE_ALIGN("# threads", Text, ROW+2);
			if (PamSampling*PamDistanceCounter) {
				snprintf(Text ,128, "%lu", PamSampling);
				WRITE_ALIGN("PAM sampling", Text, ROW+3);
				snprintf(Text ,128, "%lu", PamDistanceStart);
				WRITE_ALIGN("PAM start", Text, ROW+4);
				snprintf(Text ,128, "%lu", PamDistanceStop);
				WRITE_ALIGN("PAM stop", Text, ROW+5);
				snprintf(Text ,128, "%lu", PamDistanceStep);
				WRITE_ALIGN("PAM step", Text, ROW+6);
				snprintf(Text ,128, "%lu", PamDistanceCounter);
				WRITE_ALIGN("# PAM", Text, ROW+7);
				ROW += 5;
			}
			const size_t npoints = HeuristicDBSequenceCount*(1+PamDistanceCounter*PamSampling);
			snprintf(Text ,128, "%'lu", npoints);
			WRITE_ALIGN("# points", Text, ROW+3);
			snprintf(Text, 128, "%.2f%%", 100.0*quantile);
			WRITE_ALIGN("Quantile", Text, ROW+4);

			const SNormalizationItem * const normItem = &(prf->NormalizationData.Values[prf->HeuristicModeIndex]);
			WRITE_ALIGN("Function", CNOR[NewHeuristicNormItem.MNOR], ROW+6);
			char template[] = "R1";
			for (int i=0; i<prf->NormalizationData.JNOP[NewHeuristicNormItem.MNOR]; ++i) {
				template[1] += (unsigned char) i;
				snprintf(Text, 128, "%lf", NewHeuristicNormItem.RNOP[i]);
				WRITE_ALIGN(template, Text, ROW+7+i);
			}

#undef WRITE_ALIGN

			HPDF_Page_DrawImage (page, JPGScores, BORDER_SIZE+LEFT_MARGIN+(1.- Graph_width)*PageWidth, lastYpos - ih - 5, iw, ih);

		}

		/*
		 * Footer -------------------------------------------------------------------------
		 */
		HPDF_Page_SetFontAndSize (page, font, 6);
		HPDF_Page_BeginText (page);
		HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN, BORDER_SIZE+BOTTOM_MARGIN - 5 - 12);
		HPDF_Page_ShowText (page,  "(C) SIB Swiss Institute of Bioinformatics, Thierry Schuepbach, PFCalibrate V " PF_VERSION);
		HPDF_Page_EndText (page);
		{
			HPDF_Page_BeginText (page);
			time_t t = time(NULL);
			snprintf(Text, 128, "%s", ctime(&t));
			const HPDF_REAL tw = HPDF_Page_TextWidth (page, Text);\
			HPDF_Page_MoveTextPos (page, BORDER_SIZE+LEFT_MARGIN+PageWidth-tw, BORDER_SIZE+BOTTOM_MARGIN - 5 - 12);
			HPDF_Page_ShowText (page, Text);
			HPDF_Page_EndText (page);
		}
		/* save the document to a file */
		HPDF_SaveToFile (pdf, fname);
		HPDF_Free(pdf);
		}
#endif
	}
	exit(0);
}
