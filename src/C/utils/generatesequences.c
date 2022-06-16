/*******************************************************
                        PFTOOLS
 *******************************************************
  Mai 6, 2013 output.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdint.h>
#include <float.h>
#ifdef __INTEL_COMPILER
#include <mathimf.h>
#else
#include <math.h>
#endif
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif

#include "../include/pfProfile.h"
#include "../include/random.h"

static void fcqsort(float * const restrict Scores, char * const restrict Alphabet, const size_t N)
{
  int i,j;
  int itmp;
  float t,v;
  char ct, cv;

  if (N<=1) return;
  v = Scores[0]; cv = Alphabet[0];
  i = 0;
  j = N;
  for (;;)
    {
      while(Scores[++i] < v && i <  N) {}
      while(Scores[--j] > v) {}
      if (i>=j)
	{ break; }
      else
	{
	  t = Scores[i]; ct = Alphabet[i];
	  Scores[i] = Scores[j]; Alphabet[i] = Alphabet[j];
	  Scores[j] = t; Alphabet[j] = ct;
	}
    }
  t = Scores[i-1]; ct = Alphabet[i-1];
  Scores[i-1] = Scores[0]; Alphabet[i-1] = Alphabet[0];
  Scores[0] = t; Alphabet[0] = ct;
  fcqsort(Scores, Alphabet, i-1);
  fcqsort(Scores+i, Alphabet+i, N-i);
}


#ifndef _PFEMIT
int GenerateSequences(const struct Profile * const prf, const float MaxExponent, const char * const FileName,
		      const long int Seed, const enum GeneratedSequenceOptions Options, const unsigned int N)
#else
static int localGenerateSequences(const struct Profile * const prf, const float MaxExponent, FILE * fd,
			     const long int Seed, const enum GeneratedSequenceOptions Options, const unsigned int N)
#endif
{
  struct Random sr;
  int res = 0;

  /* Open file */
#ifndef _PFEMIT
  FILE * const fd = fopen(FileName, "w");
#endif
  if (fd == NULL) return -1;

  /* Localy store profile informations */
  const size_t Length = prf->Length;
  const char * const ProfileAlphabet = &(prf->CABC[1]); // WARNING: alphabet first value is for unknown
  size_t AlphabetLength = 1;
  while ( prf->CABC[AlphabetLength] != 'X' )  ++AlphabetLength;
  AlphabetLength -= 1;

  /* What should we take into account ? (MATCH/INSERTION/DELETION) */
  size_t TotalLength = 0;
  if (Options & GENERATE_MATCH) TotalLength += AlphabetLength;
  if (Options & GENERATE_INSERTION) TotalLength += AlphabetLength;
  if (Options & GENERATE_DELETION) TotalLength += 1;

  const size_t LDA = (TotalLength + 3) & ~3;
//    fprintf(stderr, "Total generated length is %lu rounded to %lu\n", TotalLength, LDA);

  /* Allocate memory to hold probability */
  float * const Probabilities = (float*) _mm_malloc(LDA*Length*sizeof(float), 16);
  if (Probabilities == NULL) {
    res = -2;
    goto FIN;
  }
  /* Allocate memory to hold the ordered alphabet */
  char * const Alphabets = (char*) malloc(Length*LDA*sizeof(char));
  if (Alphabets == NULL) {
    res = -2;
    goto FIN;
  }
  memset(Alphabets, '?', Length*LDA*sizeof(char));

  /* Compute the Probabilities of each profile position */
  const size_t MatchInsertionLDA = prf->Scores.Match.AlignStep;
  const StoredIntegerFormat * restrict Match = prf->Scores.Match.Alphabet; 	  // WARNING: alphabet first value is for unknown
  const StoredIntegerFormat * restrict Insertion = prf->Scores.Insertion.Alphabet;  // WARNING: alphabet[_D] is deletion
  Insertion += MatchInsertionLDA;

  float * restrict ProbPtr = Probabilities;
  char * restrict AlphaPtr = Alphabets;
  for (size_t iprf=0; iprf<Length; ++iprf) {
    register float maximum = (float) NLOW;
    size_t i=0;
    if (Options & GENERATE_MATCH) {
      while (i<AlphabetLength) {
				const register float x = (float) Match[1+i];
				ProbPtr[i] = x;
				maximum = (maximum < x) ? x : maximum;
				AlphaPtr[i] = ProfileAlphabet[i];
				++i;
      }
    }
    if (Options & GENERATE_INSERTION) {
      size_t k = 0;
      while (k<AlphabetLength) {
				const register float x = (float) Insertion[1+k];
				ProbPtr[i] = x;
				maximum = (maximum < x) ? x : maximum;
				// Set lower case to insertion to recognize it
				AlphaPtr[i++] = ProfileAlphabet[k++] + ((unsigned char) 'a' - (unsigned char) 'A');
      }
    }
    if (Options & GENERATE_DELETION) {
      const register float x = (float) Match[_D];
      ProbPtr[i] = x;
      maximum = (maximum < x) ? x : maximum;
      AlphaPtr[i++] = '-';
    }

    while (i<LDA) ProbPtr[i++] = 0.0f;

    if (maximum > 0.0f) {
      const register float scale = MaxExponent/(maximum);
      i = 0;
      while (i<TotalLength) {
				const float x = ProbPtr[i]*scale;
				ProbPtr[i++] = expf(x);
      }
    }
    else if (maximum < 0.0f) {
      const register float scale = -MaxExponent/(maximum);
      const float RealMaxExponent = 2.0f*MaxExponent;
      i = 0;
      while (i<TotalLength) {
				const float x = RealMaxExponent + ProbPtr[i]*scale;
				ProbPtr[i++] = expf(x);
      }
    }
    else {
      i = 0;
      while (i<TotalLength) {
				const float x = MaxExponent + ProbPtr[i]*MaxExponent;
				ProbPtr[i++] = expf(x);
      }
    }

    /* Sort the alphabets and the scores according to the score */
    fcqsort(ProbPtr, AlphaPtr, TotalLength);
#ifdef _PFEMIT
     fprintf(stderr, "%4zu '%.*s'\n", iprf, (int)TotalLength, AlphaPtr);
#endif

    /* Sum up the smaller one */
    float sum = ProbPtr[0];
    i = 1;
    while (i<TotalLength) {
      const register float x = ProbPtr[i];
      sum += x;
      ProbPtr[i++] = sum;
    }

    const float invSum = 1.0f/sum;
    i=0;
    while (i<(TotalLength-1)) {
      ProbPtr[i++] *= invSum;
    }
    ProbPtr[TotalLength-1] = 1.0f;

    Match     += MatchInsertionLDA;
    Insertion += MatchInsertionLDA;
    ProbPtr   += LDA;
    AlphaPtr  += LDA;
  }

  /* Initialize the random generator */
  InitializeGenerator(Seed, &sr);
  size_t SequenceBufferSize = 2*Length;
  char * SequenceBuffer = malloc(SequenceBufferSize*sizeof(char));
  if (SequenceBuffer == NULL) {
      res = -2;
      goto FIN2;
  }
  char * restrict LastSPtr = &SequenceBuffer[SequenceBufferSize-2];

  for (unsigned int iseq=0; iseq<N; ++iseq) {
    ProbPtr = Probabilities;
    AlphaPtr = Alphabets;
    const int count = snprintf(SequenceBuffer, SequenceBufferSize, ">Rnd_Seq_%u\n", 1+iseq);
    char * restrict SPtr = SequenceBuffer + count;
    size_t iprf=0;
    unsigned int CR = 0;
    unsigned int InsertCount = 0;
    while (iprf<Length) {
      /* Get a random value in [0;1] */
      register const float RndVal = NormalizedFlatDistributionValue(&sr);
      size_t ii=0;
      while (RndVal > ProbPtr[ii]) { ++ii; }
      const unsigned char * restrict cptr = &AlphaPtr[ii];
      if (SPtr >= LastSPtr) {
				/* Increase Buffer size */
				SequenceBufferSize += 256;
				LastSPtr += 32;
				char * const newSequenceBuffer = realloc(SequenceBuffer,SequenceBufferSize*sizeof(char));
				if (newSequenceBuffer == NULL) {
						res = -3;
						goto FIN2;
				}
				const size_t currentLength = (uintptr_t) SPtr - (uintptr_t) SequenceBuffer;
				SequenceBuffer = newSequenceBuffer;
				SPtr = newSequenceBuffer + currentLength;
      }
#ifndef _PFEMIT
      if (*cptr != (unsigned char) '-') {
				*SPtr++ = *cptr;
				if (++CR == 80) {
					*SPtr++ = '\n';
					CR = 0;
				}
				if ( *cptr >= (unsigned char) 'a') {
					if (++InsertCount > 100) {
							res = -4;
							goto FIN3;
					}

					continue;
				}
				InsertCount = 0;
      }
#else
      *SPtr++ = *cptr;
      if (++CR == 80) {
				*SPtr++ = '\n';
				CR = 0;
      }
      if ((*cptr != (unsigned char) '-') && (*cptr >= (unsigned char) 'a')) continue;
#endif
      ProbPtr  += LDA;
      AlphaPtr += LDA;
      ++iprf;
    }
    SPtr[0] = '\n';
    SPtr[1] = '\0';
    if (fputs(SequenceBuffer, fd) < 0) {
      res = -1;
      break;
    }
  }

  FIN3:
  free(SequenceBuffer);
  FIN2:
  free(Alphabets);
  FIN1:
  _mm_free(Probabilities);
  FIN:
#ifndef _PFEMIT
  fclose(fd);
#endif
  return res;
}

#ifdef _PFEMIT
#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <sys/time.h>

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
    "Generate sequences based on probabilities from a profile\n"
    "Used in heuristic calibration when sequences used to build the profile are missing\n"
		" pfemit [options] profile\n"
		" Options:\n"
		"   -N <uint>   : number of sequences (default 10)\n"
		"   -s <ulong>  : random generator seed (default 1234)\n"
		"   -E <float>  : maximum exponent value\n"
		"   -v          : verbose mode\n"
		"   -h          : prints this help\n\n",
	stream);
  exit(0);
}

int main(int argc, char *argv[])
{
  struct Profile prf;
  struct timeval _t0, _t1;
  char * ProfileFile;
  const char opt_to_test[]     = "N:vs:E:h";
  size_t N=10;
  bool OutputVerbose = false;
  unsigned long Seed = 1234L;
  float MaxExponent = 8.0f;

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////

  while (1) {
    const int c = getopt(argc, argv, opt_to_test);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
      case 'N':
      {
				int val = (int) atoi(optarg);
				N = (size_t) (val > 0) ? val : -val;
      }
      break;
      case 'v':
				OutputVerbose = true;
				break;
      case 's':
				Seed = atol(optarg);
				break;
      case 'E':
				MaxExponent = (float) atof(optarg);
				break;
      case 'h':
      default:
				Usage(stdout);
    }
  }

  if (optind == argc) {
    fputs("Error in given options\n", stderr);
    Usage(stderr);
  } else {
    ProfileFile = argv[optind];
  }

  /*
   * Read the profile and output some infos
   */
  gettimeofday(&_t0,0);
  const int res = ReadProfile(ProfileFile, &prf, false);
  gettimeofday(&_t1,0);
  if (res < 0) {
    fputs("Error found reading profile.\n", stderr);
    exit(1);
  }
  if (OutputVerbose) {
    const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
    fprintf(stderr, "Profile reading took %lf seconds.\n", T);

    fprintf(stderr,"Profile %s has length %lu and alphabet size of %lu\n",
	  ProfileFile, prf.Length, prf.Alphabet_Length);

    fputs("Alphabet Mapping\n",stderr);
    for (size_t i=0; i<ALPHABET_SIZE; ++i) {
      fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf.Alphabet_Mapping[i]);
      if ((i+1) % 8 == 0 ) fputs("\n",stderr);
    }
    fputs("\n\n",stderr);

    fprintf(stderr,"Disjoint set: %i to %i\n", prf.DisjointData.NDIP[0], prf.DisjointData.NDIP[1]);
  }

  if (localGenerateSequences(&prf, MaxExponent, stdout, Seed, GENERATE_MATCH | GENERATE_INSERTION | GENERATE_DELETION, N) != 0) {
    fputs("Error in generation of sequences\n", stderr);
  }

  FreeProfile(&prf, false);
}
#endif /* _PFEMIT*/
