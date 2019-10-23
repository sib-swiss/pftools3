/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfscan.c
 ***************************************************************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <pthread.h>
#include <alloca.h>
#include <stdint.h>
#include <getopt.h>
#ifdef USE_AFFINITY
# include <unistd.h>
# include <sched.h>
#endif


#define HEADER "%----------------------------------------------------------------------------%\n"\
	       "|                               PFSCAN  v" PF_VERSION "                                 |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "| Built on " __DATE__ " at " __TIME__ ".                                          |\n"
#include "../include/pfProfile.h"
#include "../include/pfRegexp.h"
#include "../include/pfSequence.h"
#include "../include/system.h"

#define _NEEDS_HEURISTIC_
#define _NEEDS_FILTER_
#define _NEEDS_ALIGNMENT_
#define _NEEDS_REGEX_
#include "threads_array.h"

static SystemInfo System;

static const char opt_to_test[] = "i:t:T:hsL:M:n0:W:o:VN:R:45fFqU:"
#ifdef USE_AFFINITY
 "02:3"
#endif
;

static const struct option long_options[] =
{
        /*
	 * These options set a flag.
	 */


        /*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	{"help",               		no_argument,       	0,	'h'},
	{"sse2",			            no_argument,        0,	's'},
	{"verbose",			          no_argument,				0,	'V'},
	/* Database indexing options */
	{"fasta",                 no_argument,				0,	'f'},
	{"fastq",                 no_argument,				0,	'q'},
	{"embl",                  no_argument,  			0,  'F'},
	{"database-index",        required_argument,	0,	'i'},
	/* Others */
	{"level" ,			          required_argument,	0,	'L'},
	{"mode",                  required_argument,	0,	'M'},
	{"no-heuristic",          no_argument,				0,	'n'},
	{"filter-normalized-cutoff", required_argument,	0,	'N'},
	{"matrix-only",	          no_argument,				0,	'4'},
	{"pattern-only",          no_argument,				0,	'5'},
	{"max-regex-match",       required_argument,	0,	'R'},
	{"unknown-symbol",				required_argument,	0,	'U'},
	/* SMP options*/
	{"nthreads",              required_argument,	0,	't'},
	{"max-heuristic-nthreads",required_argument,	0,	'T'},
#ifdef USE_AFFINITY
	{"no-affinity",           no_argument,        0,  '0'},
	{"thread-affinity",	      required_argument,	0,	'2'},
	{"no-shared-core",        no_argument,				0,	'3'},
#endif
	/* Print ouptut methods*/
	{"output-method",		required_argument,	0,	'o'},
	{"output-length",		required_argument,	0,	'W'},
	{0, 0, 0, 0}
};

unsigned int OutputPrintWidth = 60;
bool OutputVerbose = false;

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
    "Scan a protein sequence with a profile library:\n"
		" pfscanV3 [options] profiles sequences\n\n"
		" Options:\n"
		"  Profile\n"
		"   --level <int>                      [-L] : level to use for cutoff (default 0)\n"
		"   --mode <int>                            : mode to use for normalization (default 1)\n"
		"   --pattern-only                          : only pattern profile will be treated\n"
		"   --matrix-only                           : only matrix profile will be treated\n"
		"                                             default is to treat all\n"
		"   --unknown-synbol <character>            : change unknown symbol to given character\n\n"
		"  Regular expressions\n"
		"   --max-regex-match <uint>                : maximum number of returned matches per\n"
		"                                             sequence (default 16)\n\n"
		" Heuristic\n"
		"   --no-heuristic                     [-n] : bypass heuristic\n\n"
		" Filter\n"
		"   --filter-normalized-cutoff <float> [-N] : filter normalized cutoff value\n"
		"                                             heuritic cutoff will be adjusted\n"
		"                                             provided profiles data allows\n\n"
		"  Database\n"
		"   --fasta                            [-f] : FASTA file database as input\n"
		"   --fastq                            [-q] : FASTQ file database as input\n"
		"   --embl                             [-F] : EMBL SwissProt file database as input\n"
		"                                             Note that AC and ID will not be printed out!\n"
		"   --database-index <file>            [-i] : use indices stored in given file\n\n"
		" Optimizations\n"
		"   --sse2                                  : enforces SSE 2 only instruction set,\n"
		"                                             default to using SSE 4.1\n"
		"   --nthreads <uint>                  [-t] : max number of threads to use,\n"
		"                                             default to all available cores\n"
		"   --max-heuristic-nthreads <uint>         : max number of threads to use for\n"
		"                                             heuristic phase only. (IO bounds)\n"
		"                                             default to all available cores\n"
#ifdef USE_AFFINITY
		"   --no-affinity                           : disable CPU affinity file\n"
		"   --thread-affinity                       : file containing thread mask,\n"
		"                                             one row per thread\n"
		"   --no-shared-core                        : Prevent core resource sharing\n\n"
#else
		"\n"
#endif
    "  Printing output\n"
    "   --output-method <uint>     [-o] : printing output method (default 5)\n"
		"                                     == 0 replicates the pfsearch output without options\n"
    "                                     == 1 InterPro\n"
    "                                     == 2 IncMatch\n"
    "                                     == 3 PSMaker\n"
		"                                     == 4 Pfscan\n"
		"                                     == 5 Pfscan long\n"
    "                                     == 6 xPSA output\n"
		"                                     == 7 tsv output (single line tab delimited)\n"
		"                                     == 8 SAM output\n"
		"                                     == 9 Print a classification\n"
		"                                     == 10 Turtle/RDF output\n",
		stream);
	fprintf(stream,
    "   --output-length <uint>     [-W] : maximum number of column for sequence\n"
    "                                     output printing (default %u)\n"
    "  Other\n"
    "   --verbose                  [-V] : verbose on stderr\n"
    "   --help                     [-h] : output command help\n\n"
    " Version " PF_VERSION " built on " __DATE__ " at " __TIME__ ".\n",
    OutputPrintWidth);
    exit(0);
}

int main (int argc, char *argv[])
{
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL STRUCTURES
  ////////////////////////////////////////////////////////////////////////////////////////////////
  struct Profile Rootprf;		/* Root profile structure */
  struct Profile * * MatrixPrfs=NULL;	/* array of pointer to matrix profiles */
  const struct Profile * * PatternPrfs=NULL;	/* array of pointer to pattern profiles */
  struct RegEx regex;			/* Regex structure */
  DBSequence_t DB;			/* Sequence Database File */
  Sequence SeqData;			/* Sequence data to work on */
  struct timeval _t0, _t1;		/* Timing structures */

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // LOCAL DATA
  ////////////////////////////////////////////////////////////////////////////////////////////////

  PFSequence * PFSeq;			           /* Pointer to translated alphabet sequence */
  size_t nCPUs=0UL;		               /* number of threads */
  size_t nCPUsHeuristic=0UL;         /* maximum number of threads for heuristic phase */
  _Bool NoHeuristic=false;		       /* Bypass the heuristic */
  size_t HeuristicCounter = 0UL;		 /* number of sequences passing heuristic */
  size_t FilterCounter = 0UL;		     /* number of sequences passing filter */
  int res, Score;
  int Level=0;			/* Default level used from command line, if not zero then enforces that value */
  int Mode=1;				/* Default Mode for normalization from command line, if not zero then enforces that value */
  float FilterNormalizedCutoff = 0.0f;	/* Used when supplying a floating point cutoff -> normalized value */
  _Bool ImportIndices = false;		      /* Does import indices from file */
  char *ImportFileName = NULL;		      /* If so this is the file name */
  int * restrict FilterScores = NULL;   /* Array of filter scores */
  unsigned int * restrict HeuristicScores = NULL; /* Array of heuristic scores used when dumping both filter and heuristic scores */
  char * ProfileFile;			              /* Profile file */
  enum ProfileType prfType = PF_PATTERN | PF_MATRIX; /* Type of profile to allow */
  char *DBFileName;				            /* FASTA sequence file */
	_Bool isFASTA = false;							/* FASTA sequence file */
	_Bool isFASTQ = false;              /* FASTQ sequence file */
	_Bool isEMBL = false;               /* SwissProt EMBL sequence file */

  size_t * shares = 0;
  struct ID * restrict IDs = NULL;	    /* Allocate memory for profile/sequence index */
  struct ThreadData *threads_arg = NULL;/* Allocate stack memory for posix thread structures */
  pthread_t *threads = 0;
  enum Version HeuristicVersion = SSE2; /* Trigger SSE version to use for heuristic */
  enum Version OtherVersion = SSE2;	    /* Trigger SSE version to use for filter and alignment */
#ifdef USE_AFFINITY
  Affinity_Mask_t * Thread_masks[2] = {0,0};	/* Define variables to hold thread affinity mask */
  unsigned int Thread_count[2] = {0,0};
  pthread_attr_t * restrict threads_attr = NULL;
  char buffer[128] __attribute__((aligned(16))); /* buffer to read affinity file mask */
  _Bool noAffinity = false;		                   /* disable use of cpu affinity */
  _Bool split = false;
  _Bool noSharedCore = false;		    /* Prevent hyperthreading or AMD compute unit to share resources */
  _Bool GivenAffinityFile = false;	/* File holding a mask for each thread */
  char * AffinityMaskFileName;		  /* Name of affinity mask file provided by option m */
#endif
  _Bool HasGivenPrintMethod=false;
	char UnknownSymbol = '\0';

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // SYSTEM ARCHITECTURE ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  getSystemInfo(&System);

  /* Check for minimum requirement */
  if (!(System.Extensions & MM_SSE2)) {
      fputs("pfscanV3 requires at least a CPU capable of SSE 2.\n", stderr);
      exit(1);
  }

  /* Allow fast SSE 4.1 extensions ? */
  if (System.Extensions & MM_SSE41) {
      xali1_ptr = xali1_sse41;
      xalit_ptr = xalit_sse41;
      xalip_ptr = xalip_sse41;
      HeuristicVersion = SSE41;
      OtherVersion = SSE41;
  }
  else {
      xali1_ptr = xali1_sse2;
      xalit_ptr = xalit_sse2;
      xalip_ptr = xalip_sse2;
      HeuristicVersion = SSE2;
      OtherVersion = SSE2;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////
  // OPTIONS
  ////////////////////////////////////////////////////////////////////////////////////////////////
  PrintFunction = &PrintPfscanLOpt;
  while (1) {
    /* getopt_long stores the option index here. */
    int option_index = 0;

    const int c = getopt_long (argc, argv, opt_to_test, long_options, &option_index);

    /* Detect the end of the options. */
    if (c == -1) break;
    switch (c) {
#ifdef USE_AFFINITY
			case '0':
				noAffinity = true;
				break;
			case '3':
				noSharedCore = true;
				break;
			case '2':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
#endif
			case '4':
				prfType ^= PF_PATTERN;
				break;
      case '5':
				prfType ^= PF_MATRIX;
				break;
      case 'V':
				OutputVerbose = true;
				break;
      case 'o':
				{
					const int method = atoi(optarg);
					switch (method) {
						case 0: PrintFunction = &PrintDefault; break;
						case 1: PrintFunction = &PrintInterpro; break;
						case 2: PrintFunction = &PrintIncmatch; break;
						case 3: PrintFunction = &PrintPSMaker; break;
						case 4: PrintFunction = &PrintPfscan; break;
						case 5: PrintFunction = &PrintPfscanLOpt; break;
						case 6: PrintFunction = &PrintxPSA; break;
						case 7: PrintFunction = &PrintTSV; break;
						case 8: PrintFunction = &PrintSAM; break;
						case 9: PrintFunction = (PrintFunctionPtr) &PrintClassification; break;
                        case 10: PrintFunction = &PrintTurtle; break;
						default:
							fputs("Unrecognized ouput method.\n", stderr);
							exit(1);
					}
					HasGivenPrintMethod = true;
				}
				break;
      case 'W':
				OutputPrintWidth = (unsigned int) atoi(optarg);
				break;
      case 'M':
				Mode = atoi(optarg);
				break;
      case 'L':
				Level = atoi(optarg);
				break;
      case 'n':
				NoHeuristic = true;
				break;
      case 't':
				nCPUs = (size_t) atoi(optarg);
				break;
      case 'T':
				nCPUsHeuristic = (size_t) atoi(optarg);
				break;
      case 's':
				fputs("Enforcing SSE 2 as requested.\n" ,stderr);
				xali1_ptr = xali1_sse2;
				xalit_ptr = xalit_sse2;
				xalip_ptr = xalip_sse2;
				HeuristicVersion = SSE2;
				OtherVersion = SSE2;
				break;
			case 'f':
				isFASTA = true;
				isEMBL = false;
				break;
			case 'q':
				isFASTQ = true;
				isEMBL = false;
				break;
			case 'F':
				isEMBL = true;
				break;
      case 'i':
				ImportIndices = true;
				ImportFileName = optarg;
				break;
      case 'N':
				FilterNormalizedCutoff = (float) atof(optarg);
				break;
      case 'R':
				MaxRegexNumber = (size_t) atoi(optarg);
				break;
			case 'U':
				UnknownSymbol = optarg[0];
				break;
      case 'h':
      default:
				Usage(stdout);
    }
  }

  if (optind >= argc) {
    fputs("Expected arguments after options\n", stderr);
    Usage(stderr);
  } else {
    ProfileFile = argv[optind];
    DBFileName = argv[optind+1];
  }

  if (OutputVerbose) {
   fputs(HEADER
#ifdef __USE_MMAP__
        "| Using Linux kernel MMAP function.                                          |\n"
#endif
      ,stderr);
    printSystemInfo(stderr, &System);
    if (OtherVersion == SSE2 && (System.Extensions & MM_SSE41)) {
	fputs("Enforcing SSE 2...\n", stderr);
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // INPUT ANALYSIS
  ////////////////////////////////////////////////////////////////////////////////////////////////

  /* Check profile type argument */
  if ((int) prfType == 0) {
      fputs("Try to keep at least one type of profile!\n", stderr);
      exit(1);
  }

  /*
   * Read the profile and output some infos ------------------------------------------------------
   */
  gettimeofday(&_t0,0);
  const int ProfileCount = ReadProfile(ProfileFile, &Rootprf, true);
  gettimeofday(&_t1,0);
  if (ProfileCount < 0) {
    fputs("Error found reading profile.\n", stderr);
    exit(1);
  }
  else if (ProfileCount == 1) {
    fputs("pfscanV3 is not meant to be used with a single profile, use pfsearchV3 to get better performance in such case.\n",stderr);
  }

  /* Create an array of pointer to Profiles separating pattern from matrix */
  int MatrixProfileCount = 0, PatternProfileCount = 0;
  if (ProfileCount > 0) {
    const struct Profile * tmpPrf = &Rootprf;
    do {
      if (tmpPrf->Type == PF_MATRIX)
				++MatrixProfileCount;
      else
			++PatternProfileCount;
      tmpPrf = tmpPrf->next;
    } while (tmpPrf);

    if (MatrixProfileCount) MatrixPrfs = (struct Profile **) alloca(MatrixProfileCount*sizeof(struct Profile*));
    if (PatternProfileCount) PatternPrfs = (const struct Profile **) alloca(PatternProfileCount*sizeof(struct Profile*));

    size_t iMatrixPrf = 0, iPatternPrf = 0;
    tmpPrf = &Rootprf;
    do {
      if (tmpPrf->Type == PF_MATRIX) {
				MatrixPrfs[iMatrixPrf++] = (struct Profile *) tmpPrf;
      }
      else {
				PatternPrfs[iPatternPrf++] = (struct Profile *) tmpPrf;
      }
      tmpPrf = tmpPrf->next;
    } while (tmpPrf);
  }

  if (OutputVerbose) {
    const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
    fprintf(stderr, "Profiles reading took %lf seconds: %i pattern, %i matrix found\n", T, PatternProfileCount, MatrixProfileCount);
  }

  if ((prfType & PF_MATRIX) == 0) {
      MatrixProfileCount = 0;
      if (OutputVerbose) fputs("Matrix Profile will not be computed.\n", stderr);
  }
  if ((prfType & PF_PATTERN) == 0) {
      PatternProfileCount = 0;
      if (OutputVerbose) fputs("Pattern Profile will not be computed.\n", stderr);
  }

  if (OutputVerbose) {
    for (int i=0; i<MatrixProfileCount;++i) {
      const struct Profile * const prf = MatrixPrfs[i];
      fprintf(stderr,"Profile %s has length %lu and alphabet size of %lu\n",
	    prf->Identification, prf->Length, prf->Alphabet_Length);

      fputs("Alphabet Mapping\n",stderr);
      for (size_t i=0; i<ALPHABET_SIZE; ++i) {
				fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
				if ((i+1) % 8 == 0 ) fputs("\n",stderr);
      }
      fputs("\n",stderr);

      fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
    }
  }

  /* Set the Cutoff Level and the Normalization */
  if (MatrixProfileCount) {
    int iMatrixPrf=0;
    do {
      struct Profile * const prf = MatrixPrfs[iMatrixPrf];
      register const int res = SetProfileLevelAndMode(prf, Level, Mode);
      if (res < 0) {
				fprintf(stderr, "Profile %s does not contain cutoff with level set to %i and mode set to %i, error %i\n",
				        prf->Identification, Level, Mode, res);
				exit(1);
      }
      if (UnknownSymbol) prf->CABC[0] = UnknownSymbol;
    } while (++iMatrixPrf < MatrixProfileCount);
    if (OutputVerbose) {
      fprintf(stderr, "Using level %i, mode %i\n", Level, Mode);
      for (int i=0; i<MatrixProfileCount; ++i) {
				fprintf(stderr,"Profile %s filter cutoff (%i), heuristic cutoff(%u)\n",
								MatrixPrfs[i]->Identification, MatrixPrfs[i]->CutOff, MatrixPrfs[i]->HeuristicCutOff);
      }
    }
  }

  /*
   * Read the FASTA file -------------------------------------------------------------------------
   */
	if (isFASTA) {
		gettimeofday(&_t0,0);
		if (!ImportIndices) {
			res = AnalyzeFASTAStructure(DBFileName, &DB);
		}
		else {
			FILE* inIndex = fopen(ImportFileName, "rb");
			if (inIndex != NULL) {
				res = ImportDBStructure(inIndex, DBFileName, &DB);
				fclose(inIndex);
			}
			else {
				fprintf(stderr,"Unable to open index file %s, will analyze database instead.\n",ImportFileName);
				res = 1;
			}
		}
		gettimeofday(&_t1,0);
		if (OutputVerbose) {
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Sequence file %s indexing took %lf seconds.\n", DBFileName, T);
		}
		if (res != 0) {
			fputs("Error found.\n", stderr);
			exit(1);
		}
		if (OutputVerbose) {
			fprintf(stderr,
							"FASTA file %s analyzed\n"
							"\tFound %lu sequences within %lu bytes\n"
							"\tBiggest sequence entry is %lu bytes\n",
						  DBFileName, DB.SequenceCount, DB.FileSize, DB.MaxSequenceSize);
		}
	}
	else if (isFASTQ) {
		gettimeofday(&_t0,0);
		if (!ImportIndices) {
			res = AnalyzeFASTQStructure(DBFileName, &DB);
		}
		else {
			FILE* inIndex = fopen(ImportFileName, "rb");
			if (inIndex != NULL) {
				res = ImportDBStructure(inIndex, DBFileName, &DB);
				fclose(inIndex);
			}
			else {
				fprintf(stderr,"Unable to open index file %s, will analyze database instead.\n",ImportFileName);
				res = 1;
			}
		}
		gettimeofday(&_t1,0);
		if (OutputVerbose) {
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Sequence file %s indexing took %lf seconds.\n", DBFileName, T);
		}
		if (res != 0) {
			fputs("Error found.\n", stderr);
			exit(1);
		}
		if (OutputVerbose) {
			fprintf(stderr,
							"FASTQ file %s analyzed\n"
							"\tFound %lu sequences within %lu bytes\n"
							"\tBiggest sequence entry is %lu bytes\n",
						  DBFileName, DB.SequenceCount, DB.FileSize, DB.MaxSequenceSize);
		}
	}
	else if (isEMBL) {
		gettimeofday(&_t0,0);
		if (!ImportIndices) {
			res = AnalyzeEMBLStructure(DBFileName, &DB);
		}
		else {
			FILE* inIndex = fopen(ImportFileName, "rb");
			if (inIndex != NULL) {
				res = ImportDBStructure(inIndex, DBFileName, &DB);
				fclose(inIndex);
			}
			else {
				fprintf(stderr,"Unable to open index file %s, will analyze database instead.\n",ImportFileName);
				res = 1;
			}
		}
		gettimeofday(&_t1,0);
		if (OutputVerbose) {
			const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr, "Sequence file %s indexing took %lf seconds.\n", DBFileName, T);
		}
		if (res != 0) {
			fputs("Error found.\n", stderr);
			exit(1);
		}
		if (OutputVerbose) {
			fprintf(stderr,
							"SwissProt/EMBL file %s analyzed\n"
							"\tFound %lu sequences within %lu bytes\n"
							"\tBiggest sequence entry is %lu bytes\n",
						  DBFileName, DB.SequenceCount, DB.FileSize, DB.MaxSequenceSize);
		}
	}

	/*
   * Retrieve number of cores --------------------------------------------------------------------
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
    if (OutputVerbose) fprintf(stderr,"Parsing file %s for affinity mask and number of threads\n", AffinityMaskFileName);
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
  else if ( split ) {
    // -----------------------------------------------------------------------------
    //                 ***  HALF SSE 2 HALF SSE 4.1 HYPERTHREADING***
    // -----------------------------------------------------------------------------
    Thread_count[0] = getMasks(&System, -1, -1, 1, &Thread_masks[0]);
    if (Thread_count[0] == 0) {
      fputs("No potential affinity mask found !!!\n", stderr);
      exit(0);
    }
    Thread_count[1] = getMasks(&System, -1, -1, 2, &Thread_masks[1]);
    if (Thread_count[1] == 0) {
      fputs("No potential affinity mask found with hyperthreading !!!\n", stderr);
      exit(0);
    }
    if (OutputVerbose) fprintf(stderr, "%u threads will use SSE 4.1 and %u SSE 2\n", Thread_count[0], Thread_count[1]);
  }
  else if (noSharedCore) {
    if (OutputVerbose) fputs("No sharing of core resources will be used: Intel Hyperthreading or AMD Compute Unit\n", stderr);
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
  if (OutputVerbose) fprintf(stderr, "Job dispatched over %lu cores.\n", nCPUs);

#ifdef __USE_MMAP__
  const int fd = open(DBFileName, O_RDONLY );
  const size_t length = (size_t) DB.FileSize;
  const char * const SequenceFileMap = (const char * const restrict) mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
  if (SequenceFileMap == NULL) {
    fputs("Unable to map sequence file to memory\n", stderr);
    exit(1);
  }
  const char * const restrict SequenceFile = SequenceFileMap;
	close(fd);
#else
  const char * const restrict SequenceFile = DBFileName;
#endif

  /* Prepare structure common to filter and alignment */
  shares = alloca((nCPUs+1)*sizeof(size_t));

  /* Allocate stack memory for posix thread structures */
  threads_arg = alloca(nCPUs*sizeof(struct ThreadData));
  threads = (pthread_t*) alloca(nCPUs*sizeof(pthread_t));

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // PATTERN MODE
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (PatternProfileCount) {
    /* Initialize regex structure */
    const int res = InitRegExFromProfile(&regex, PatternPrfs, PatternProfileCount, MaxRegexNumber);
    if (res<0) {
      fprintf(stderr, "Regex initialization return error code %i\n", res);
      exit(1);
    }

    /* Share according to file size */
    {
      size_t FileShare = (size_t) DB.FileSize / nCPUs;
      FileShare += ((size_t) DB.FileSize % nCPUs) > (nCPUs-1) ? 1 : 0;
      const s_Data * DataPtr = DB.DataPtr;
      register size_t counter = 0;
      shares[0] = 0;
      for (size_t i=1; i<nCPUs; ++i) {
				register size_t tmp = i*FileShare;
				while ( (size_t) DataPtr->Sequence.Offset < tmp) { ++DataPtr; ++counter; }
				shares[i] = counter;
      }
      shares[nCPUs] = DB.SequenceCount;
    }

    gettimeofday(&_t0,0);
    /* Dispatch to threads */
    for (size_t i=0; i<nCPUs; ++i) {
      threads_arg[i].prf           = PatternPrfs;
      threads_arg[i].regex	       = &regex;
      threads_arg[i].FASTA         = &DB;
      threads_arg[i].SequenceFile  = SequenceFile;
      threads_arg[i].profileCount  = (size_t )PatternProfileCount;
      threads_arg[i].threadId      = i;
      threads_arg[i].start         = shares[i];
      threads_arg[i].stop          = shares[i+1];
#ifdef USE_AFFINITY
      if (pthread_create (&threads[i],  &threads_attr[i], thread_regex,  (void*) &threads_arg[i]) != 0)
#else
      if (pthread_create (&threads[i],  NULL, thread_regex,  (void*) &threads_arg[i]) != 0)
#endif
      {
				fputs("Fail to create thread.\n", stderr);
				exit(0);
      }
    }

    for (size_t i=0; i<nCPUs; ++i) {
      pthread_join(threads[i], NULL);
    }
    gettimeofday(&_t1,0);
    double t;
    if (OutputVerbose) {
      t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
      fprintf(stderr,"Pattern lookup took %lf seconds to treat on %li cores.\n", t, nCPUs);
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // MATRIX MODE
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (MatrixProfileCount) {
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // HEURISTIC
    ////////////////////////////////////////////////////////////////////////////////////////////////

    /* Allocate memory for sequence YesNo or ID to be done */
    IDs = _mm_malloc(DB.SequenceCount*MatrixProfileCount*sizeof(struct ID), 16);
    if (IDs == NULL) {
			fputs("Cannot allocate memory.\n", stderr);
			goto END;
    }

     /* Common base informations */
    for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].prf              = (const struct Profile **) MatrixPrfs;
			threads_arg[i].FASTA            = &DB;
			threads_arg[i].Array            = IDs;
			threads_arg[i].SequenceFile     = SequenceFile;
			threads_arg[i].profileCount     = (size_t )MatrixProfileCount;
			threads_arg[i].threadId         = i;
			threads_arg[i].NormalizedCutoff = FilterNormalizedCutoff;
		}

    if (!NoHeuristic) {
      /* Compute Match Score Matrix transpose */
      gettimeofday(&_t0,0);
      TransposeMatrix * TIMatch = (TransposeMatrix*) alloca(MatrixProfileCount*sizeof(TransposeMatrix));
#ifdef USE_AFFINITY
      TransposeMatrix * TIMatch1 = (TransposeMatrix*) alloca(MatrixProfileCount*sizeof(TransposeMatrix));
      if (split) {
				for (int i=0;i<MatrixProfileCount;++i) {
					TIMatch[i].i  = TransposeAndConvertMatchMatrix(&(MatrixPrfs[i]->Scores), MatrixPrfs[i]->Alphabet_Length, MatrixPrfs[i]->Length);
					TIMatch1[i].f = TransposeAndConvertToFloatMatchMatrix(&(MatrixPrfs[i]->Scores), MatrixPrfs[i]->Alphabet_Length, MatrixPrfs[i]->Length);
				}
      } else
#endif
      if (HeuristicVersion == SSE41) {
				for (int i=0;i<MatrixProfileCount;++i)
					TIMatch[i].i = TransposeAndConvertMatchMatrix(&(MatrixPrfs[i]->Scores), MatrixPrfs[i]->Alphabet_Length, MatrixPrfs[i]->Length);
			}
			else {
				for (int i=0;i<MatrixProfileCount;++i)
					TIMatch[i].f = TransposeAndConvertToFloatMatchMatrix(&(MatrixPrfs[i]->Scores), MatrixPrfs[i]->Alphabet_Length, MatrixPrfs[i]->Length);
      }
      gettimeofday(&_t1,0);
      if (OutputVerbose) {
				const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
				fprintf(stderr,"Transposing Match matrices took %lf seconds.\n", t);
      }

      /* Limit number of threads for heuristic */
      if ( nCPUsHeuristic == 0) nCPUsHeuristic = nCPUs;

      /* Share according to file size */
      {
				size_t FileShare = (size_t)DB.FileSize / nCPUsHeuristic;
				FileShare += ((size_t)DB.FileSize % nCPUsHeuristic) > (nCPUsHeuristic-1) ? 1 : 0;
				const s_Data * DataPtr =DB.DataPtr;
				register size_t counter = 0;
				shares[0] = 0;
				for (size_t i=1; i<nCPUsHeuristic; ++i) {
					register size_t tmp = i*FileShare;
					while ( (size_t) DataPtr->Sequence.Offset < tmp) { ++DataPtr; ++counter; }
					shares[i] = counter;
				}
				shares[nCPUsHeuristic] =DB.SequenceCount;
      }
      gettimeofday(&_t0,0);
#ifdef USE_AFFINITY
      if (split) {
				for (size_t i=0; i<Thread_count[0]; ++i) {
					threads_arg[i].start          = shares[i];
					threads_arg[i].stop           = shares[i+1];
					threads_arg[i].TransposeMatch = TIMatch;
					threads_arg[i].version        = SSE41;
					if (pthread_create (&threads[i],  &threads_attr[i], thread_heuristic,  (void*) &threads_arg[i]) != 0)
					{
						fputs("Fail to create thread.\n", stderr);
						exit(0);
					}
				}
				for (size_t i=0; i<Thread_count[1]; ++i) {
					threads_arg[Thread_count[0]+i].start          = shares[Thread_count[0]+i];
					threads_arg[Thread_count[0]+i].stop           = shares[Thread_count[0]+i+1];
					threads_arg[Thread_count[0]+i].TransposeMatch = TIMatch1;
					threads_arg[i].version                        = SSE2;
					if (pthread_create (&threads[Thread_count[0]+i],  &threads_attr[Thread_count[0]+i], thread_heuristic,  (void*) &threads_arg[Thread_count[0]+i]) != 0)
					{
						fputs("Fail to create thread.\n", stderr);
						exit(0);
					}
				}
      } else {
#endif
				for (size_t i=0; i<nCPUsHeuristic; ++i) {
					threads_arg[i].start          = shares[i];
					threads_arg[i].stop           = shares[i+1];
					threads_arg[i].TransposeMatch = TIMatch;
					threads_arg[i].version        = HeuristicVersion;
#ifdef USE_AFFINITY
					if (pthread_create (&threads[i],  &threads_attr[i], thread_heuristic,  (void*) &threads_arg[i]) != 0)
#else
					if (pthread_create (&threads[i],  NULL, thread_heuristic,  (void*) &threads_arg[i]) != 0)
#endif
					{
						fputs("Fail to create thread.\n", stderr);
						exit(0);
					}
				}
#ifdef USE_AFFINITY
      }
#endif
      for (size_t i=0; i<nCPUsHeuristic; ++i) {
				pthread_join(threads[i], NULL);
      }
      gettimeofday(&_t1,0);
      double t;
      if (OutputVerbose) {
				t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
				fprintf(stderr,"Heuristic took %lf seconds to treat on %li cores.\n", t, nCPUsHeuristic);
      }
      for (int i=0;i<MatrixProfileCount;++i)  _mm_free(TIMatch[i].f);
  #ifdef USE_AFFINITY
      if (split) for (int i=0;i<MatrixProfileCount;++i) _mm_free(TIMatch1[i].f);
  #endif

      /* Gather the one that passed the heuristic */
      HeuristicCounter = 0;
      const size_t limit =DB.SequenceCount*MatrixProfileCount;
      for (size_t iseq=0; iseq<limit; ++iseq) {
				if (IDs[iseq].PrfId >= 0) {
					IDs[HeuristicCounter].PrfId = IDs[iseq].PrfId;
					IDs[HeuristicCounter].SeqId = IDs[iseq].SeqId;
					++HeuristicCounter;
				}
      }
      if (OutputVerbose)
				fprintf(stderr,"Overall there are %lu/%lu sequences passing heuristic. These took %lf seconds to treat on %lu cores.\n",
				        HeuristicCounter,DB.SequenceCount*MatrixProfileCount, t, nCPUsHeuristic);
    }
    else {
      HeuristicCounter = 0;
      for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
				const register unsigned int lseq = (unsigned int) iseq;
				for (int iprf=0; iprf<MatrixProfileCount; iprf++) {
					IDs[HeuristicCounter].PrfId = iprf;
					IDs[HeuristicCounter].SeqId = lseq;
					HeuristicCounter++;
				}
      }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // FILTER
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (HeuristicCounter) {
      size_t maxProfileSize = 0L;
      for (int i=0; i<MatrixProfileCount; ++i) {
				const struct Profile * const prf = MatrixPrfs[i];
				maxProfileSize = maxProfileSize < prf->Length ? prf->Length : maxProfileSize;
      }

      /* Compute the new share for each thread */
      size_t SequenceShare = HeuristicCounter / nCPUs;
      SequenceShare += (HeuristicCounter % nCPUs) > (nCPUs-1) ? 1 : 0;
      shares[0] = 0;
      for (size_t i=1; i<nCPUs; ++i) shares[i] = i*SequenceShare;

      shares[nCPUs] = HeuristicCounter;

      /* Dispatch to threads */
      {
				gettimeofday(&_t0,0);
				for (size_t i=0; i<nCPUs; ++i) {
					threads_arg[i].start            = shares[i];
					threads_arg[i].stop             = shares[i+1];
					threads_arg[i].Array            = IDs;
					threads_arg[i].counter          = 0;
					threads_arg[i].MaxProfileSize   = maxProfileSize;
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

      /* Gather the one that passed xali1 */
      FilterCounter = 0;
      for (size_t iseq=0; iseq<HeuristicCounter; ++iseq) {
				if (IDs[iseq].PrfId >= 0 ) {
					IDs[FilterCounter].PrfId = IDs[iseq].PrfId;
					IDs[FilterCounter].SeqId = IDs[iseq].SeqId;
					++FilterCounter;
				}
      }
      if (OutputVerbose) {
				const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
				fprintf(stderr,"Overall there are %lu/%lu sequences passing filter. These took %lf seconds to treat on %li cores.\n",
				        FilterCounter, HeuristicCounter, t, nCPUs);
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // ALIGNMENT
    ////////////////////////////////////////////////////////////////////////////////////////////////
    if (FilterCounter > 0) {
			/* Output header for some printing function */
			if (PrintFunction == &PrintSAM) {

				printf("@HD\tVN:1.6\tSO:coordinate\n");
				for (int i=0; i<MatrixProfileCount;++i) {
					const struct Profile * const prf = MatrixPrfs[i];
					/* RNAME */
					const char * cptr = prf->AC_Number;
					while(*cptr != ';' && *cptr != '\0') { cptr++; }
					const int AClen = (int) ((uintptr_t) cptr - (uintptr_t) prf->AC_Number);
					printf("@SQ\tSN:%.*s|%s\tLN:%zu\tDS:%s\n",
							 	 AClen, prf->AC_Number, prf->Identification, prf->Length, prf->Description);
				}
			} else if (PrintFunction == &PrintTurtle) {
                printf("PREFIX ys:<http://example.org/yoursequence/>\n");
                printf("PREFIX yr:<http://example.org/yourrecord/>\n");
                printf("PREFIX up:<http://purl.uniprot.org/core/>\n");
                printf("PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns#>\n");
                printf("PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>\n");
                printf("PREFIX faldo:<http://biohackathon.org/resource/faldo#>\n");
                printf("PREFIX edam:<http://edamontology.org/>\n");
                if (strstr(ProfileFile, "hamap") > 0) {
                    printf("PREFIX profile:<http://purl.uniprot.org/hamap/>\n");
                } else {
                    printf("PREFIX profile:<http://example.org/yourprofiledb>\n");
                }
                
            }

      /* Initialize the print mutex */
      pthread_mutex_init(&PrintLock, NULL);

      /* Compute the new share for each thread */
      {
				size_t SequenceShare = FilterCounter / nCPUs;
				SequenceShare += (FilterCounter % nCPUs) > (nCPUs-1) ? 1 : 0;
				shares[0] = 0;
				for (size_t i=1; i<nCPUs; ++i) {
					shares[i] = i*SequenceShare;
				}
				shares[nCPUs] = FilterCounter;
      }

      /* Dispatch to threads */
      gettimeofday(&_t0,0);
      for (size_t i=0; i<nCPUs; ++i) {
				threads_arg[i].start       = shares[i];
				threads_arg[i].stop        = shares[i+1];
#ifdef USE_AFFINITY
				if (pthread_create (&threads[i],  &threads_attr[i], thread_xaliPT,  (void*) &threads_arg[i]) != 0)
#else
				if (pthread_create (&threads[i],  NULL, thread_xaliPT,  (void*) &threads_arg[i]) != 0)
#endif
				{
					return 1;
				}
      }

      for (size_t i=0; i<nCPUs; i++) {
				pthread_join(threads[i], NULL);
      }
      gettimeofday(&_t1,0);

      unsigned int AlignedSequencesCounter = threads_arg[0].counter;
      for (size_t i=1; i<nCPUs; i++) AlignedSequencesCounter += threads_arg[i].counter;

      if (OutputVerbose) {
				const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
				fprintf(stderr,"Overall there are %u aligned sequences found. These took %lf seconds to align on %li cores.\n", AlignedSequencesCounter, t, nCPUs);
      }

      /* Free the print mutex */
      pthread_mutex_destroy(&PrintLock);
    }
    _mm_free(IDs);
  }


  /* Free Memory */
END:
  FreeProfile(&Rootprf, false);

#ifdef __USE_MMAP__
  munmap((void*)SequenceFileMap, length);
#endif
#ifdef USE_AFFINITY
  if (Thread_masks[0]) free(Thread_masks[0]);
  if (Thread_masks[1]) free(Thread_masks[1]);
#endif

  FreeDBStructure(&DB);
  freeSystemInfo(&System);

  exit(0);
}
