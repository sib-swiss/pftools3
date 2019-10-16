/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
 *  Oct 3, 2011 pfsearch.c
 ***************************************************************************************************
 * (C) 2011 SIB Swiss Institute of Bioinformatics
 *     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#ifdef __USE_WINAPI__
#include <windows.h>
#else
#include <pthread.h>
#endif
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include <stdint.h>
#include <getopt.h>
#ifdef USE_AFFINITY
#include <unistd.h>
#include <sched.h>
#endif


#define HEADER \
"%--------------------------------------------------------------------%\n"\
"|                           PFSEARCH v" PF_VERSION "                          |\n"\
"%--------------------------------------------------------------------%\n"\
"| Built on " __DATE__ " at " __TIME__ ".                                  |\n"
#define __USE_INLINE_FUNCTIONS__
#include "../include/pfProfile.h"
#include "../include/pfRegexp.h"
#include "../include/pfSequence.h"
#include "../include/system.h"

#define _NEEDS_HEURISTIC_
#define _NEEDS_FILTER_
#define _NEEDS_ALIGNMENT_
#define _NEEDS_REGEX_
#define _NEEDS_OPTIMAL_ALIGNMENT_
// #define _NEEDS_SCORES_ONLY_
#include "threads.h"

static SystemInfo System;

static const char opt_to_test[] = "C:H:u:QI:i:t:T:hPpasdDL:M:Rrn0:W:o:VxbyN:fqFU:"
#ifdef USE_AFFINITY
"01:2"
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
	/* Profile */
	{"level",                     required_argument,	0,	'L'},
	{"mode",                      required_argument,	0,	'M'},
	{"reverse-profile",           no_argument,				0,	'R'},
	{"unknown-symbol",						required_argument,	0,	'U'},
	/* Sequence */
	{"optimal",                   no_argument,				0,	'a'},
	{"complement",                no_argument,				0,	'x'},
	{"reverse-sequence",          no_argument,				0,	'r'},
	{"satisfy-cutoff",            no_argument,				0,	'y'},
	{"both",                      no_argument,				0,	'b'},
	/* Heuristic */
	{"dump-heuristic-scores",     no_argument,				0,	'P'},
	{"dump-filter-sequences",     no_argument,				0,	'd'},
	{"heuristic-cutoff",          required_argument,	0,	'H'},
	{"no-heuristic",              no_argument,				0,	'n'},
	/* Filter */
	{"dump-filter-scores",        no_argument,				0,	'p'},
	{"dump-alignment-sequences",  no_argument,				0,	'D'},
	{"filter-cutoff",	            required_argument,	0,	'C'},
	{"filter-normalized-cutoff",  required_argument,  0,	'N'},
	/* Regex */
	{"max-regex-match",	          required_argument,	0,	'u'},
	{"header-only",	              no_argument,				0,	'Q'},
	/* Database indexing options */
	{"fasta",                     no_argument,				0,	'f'},
	{"fastq",                     no_argument,				0,	'q'},
	{"embl",                      no_argument,  			0,  'F'},
	{"database-index",            required_argument,	0,	'i'},
	/* Others */
	{"help",                      no_argument,				0,	'h'},
	{"sse2",                      no_argument,				0,	's'},
	{"verbose",                   no_argument,				0,	'V'},
	/* SMP options*/
	{"nthreads",                  required_argument,	0,	't'},
	{"max-heuristic-nthreads",    required_argument,	0,	'T'},
#ifdef USE_AFFINITY
	{"no-affinity",               no_argument,				0,	'0'},
	{"thread-affinity",           required_argument,	0,	'1'},
	{"no-shared-core",            no_argument,				0,	'2'},
#endif
	/* Print ouptut methods*/
	{"output-method",	            required_argument,	0,	'o'},
	{"output-length",	            required_argument,	0,	'W'},
	{0, 0, 0, 0}
};

unsigned int OutputPrintWidth = 60;
_Bool OutputVerbose = false;
size_t MaxRegexNumber = 16;

static void __attribute__((noreturn)) Usage(FILE * stream)
{
	fputs(
		"Scan a protein sequence library for profile matches:\n"
		" pfsearchV3 [options] [profile|regex{...}|pattern{...}] database\n\n"
		" Options:\n"
		"  Profile\n"
		"   --level <int>                      [-L] : level to use for cutoff (default 0)\n"
		"   --mode <int>                            : mode to use for normalization (default 1)\n"
		"   --reverse-profile                       : reverse the profile\n"
		"   --unknown-symbol <character>            : change unknown symbol to given character\n\n"
		"  Sequence\n"
		"   --optimal                          [-a] : report optimal alignment for all sequences\n"
		"     --satisfy-cutoff                      : but satisfying cutoff\n"
		"   --reverse-sequence                      : read sequence backward\n"
		"   --complement                            : does a complementary sequence alignment.\n"
		"                                             note that --reverse is NOT automatic\n"
		"   --both                             [-b] : compute both forward and reverse complemented\n\n"
		"  Regular expressions\n"
		"   --max-regex-match <uint>                : maximum number of returned matches per sequence\n"
		"   --header-only                           : search in headers only\n\n"
		"  Database\n"
		"   --fasta                            [-f] : FASTA file database as input\n"
		"   --fastq                            [-q] : FASTQ file database as input\n"
		"   --embl                             [-F] : EMBL SwissProt file database as input\n"
		"                                             Note that AC and ID will not be printed out!\n"
		"   --database-index <file>            [-i] : use indices stored in given file (optional)\n\n"
		"  Heuristic\n"
		"   --no-heuristic                     [-n] : bypass heuristic\n"
		"   --heuristic-cutoff <uint>          [-H] : heuristic cutoff value\n"
		"   --dump-heuristic-scores                 : only print heuristic scores to stdout\n\n"
		"  Filter\n"
		"   --filter-cutoff <int>              [-C] : filter raw cutoff value\n"
		"   --filter-normalized-cutoff <float> [-N] : filter normalized cutoff value\n"
		"   --dump-filter-scores                    : only print filter scores to stdout\n"
		"   --dump-filter-sequences            [-d] : dump passed heuristic sequences\n\n"
		"  Alignment\n"
		"   --dump-alignment-sequences         [-D] : dump passed heuristic and filter\n"
		"                                             sequences\n\n"
		"  Optimizations\n"
		"   --sse2                                  : enforces SSE 2 only instruction set,\n"
		"                                             default to using SSE 4.1\n"
		"   --nthreads <uint>                  [-t] : max number of threads to use\n"
		"                                             default to all available cores\n"
		"   --max-heuristic-nthreads <uint>         : max number of threads to use for\n"
		"                                             heuristic phase only. (IO bounds)\n"
		"                                             default to all available cores\n"
#ifdef USE_AFFINITY
		"   --no-affinity                           : disable CPU affinity file\n"
		"   --thread-affinity                       : file containing thread mask, one row per thread\n"
		"   --no-shared-core                        : Prevent core resource sharing\n\n"
#else
		"\n"
#endif
		"  Printing output\n"
		"   --output-method <uint>  [-o] : printing output method (default 0)\n"
		"                                     == 0 replicates the pfsearch output without options\n"
    "                                     == 1 InterPro\n"
    "                                     == 2 IncMatch\n"
    "                                     == 3 PSMaker\n"
		"                                     == 4 Pfscan\n"
		"                                     == 5 Pfscan long\n"
    "                                     == 6 xPSA output\n"
		"                                     == 7 tsv output (single line tab delimited)\n"
		"                                     == 8 SAM output\n"
		"                                     == 9 Family classification (ONLY for pfsearch)\n"
		"                                     == 10 Turtle/RDF output for HAMAP as SPARQL style rules (pfsearch)\n",
		stream);
	fprintf(stream,
		"   --output-length <uint>  [-W] : maximum number of column for sequence\n"
		"                                  output printing (default %u)\n"
		"  Other\n"
		"   --verbose               [-V] : verbose on stderr\n"
		"   --help                  [-h] : output command help\n\n"
		" Version " PF_VERSION " built on " __DATE__ " at " __TIME__ ".\n", OutputPrintWidth);
	exit(0);
}

int main (int argc, char *argv[])
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL STRUCTURES
	////////////////////////////////////////////////////////////////////////////////////////////////
	struct Profile * restrict prf;    /* Profile */
	struct RegEx regex;               /* Regex structure */
	DBSequence_t DB;                  /* Sequence Database File */
	Sequence SeqData;                 /* Sequence data to work on */
	struct timeval _t0, _t1;          /* Timing structures */

	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL DATA
	////////////////////////////////////////////////////////////////////////////////////////////////

	PFSequence * PFSeq;	                /* Pointer to translated alphabet sequence */
	size_t nCPUs=0;                     /* number of threads */
	size_t nCPUsHeuristic=0;            /* maximum number of threads for heuristic phase */
	_Bool NoHeuristic=false;            /* Bypass heuristic */
	size_t HeuristicCounter = 0;        /* number of sequences passing heuristic */
	size_t FilterCounter = 0;           /* number of sequences passing filter */
	int res, Score;
	int HeuristicCutOff = -1;           /* Default heuristic cutoff from command line, if not zero then enforces that value*/
	unsigned int FilterCutoff = 0;      /* Default filter cutoff from command line, if not zero then enforces that value */
	float FloatingFilterCutoff = 0.0f;  /* Used when supplying a floating point cutoff rather than integer -> normalized value */
	_Bool NormalizedCutoff = false;     /* Used to trigger the above */
	int Level=0;                        /* Default level used from command line, if not zero then enforces that value */
	int Mode=1;	                        /* Default Mode for normalization from command line, if not zero then enforces that value */
	_Bool ImportIndices = false;        /* Does import indices from file */
	char *ImportFileName = NULL;        /* If so this is the file name */
	_Bool DumpHeuristicOnly = false;    /* If set only dump heuritic scores for each sequence */
	_Bool DumpFilterOnly = false;       /* If set only dump filter scores for each sequence */
	_Bool DumpFilterSequences = false;  /* Dump sequence that passed the heuristic */
	_Bool DumpAlignmentSequences = false;		/* Dump sequence that passed the heuristic and the filter */
	_Bool DumpBoth = false;             /* Dump both heuristic and filter scores */
	_Bool IsReversed = false;           /* Will reverse profile before analyzing, optimal search only */
	_Bool ReverseSequence = false;      /* Will reverse all sequences before analysing */
	_Bool Optimal = false;              /* Only report best alignment */
	_Bool DoBothFwdAndRev = false;      /* Perform both forward and reverse complement */
	enum Strand strand = FORWARD;       /* Default strand */
	int * restrict FilterScores = NULL; /* Array of filter scores */
	unsigned int * restrict HeuristicScores=NULL;	/* Array of heuristic scores used when dumping both filter and heuristic scores */
	unsigned int * restrict SequenceID = NULL;		/* Allocate memory for sequence ID to be done */
	char * ProfileFile = NULL;          /* Profile file */
	_Bool isFASTA = false;							/* FASTA sequence file */
	_Bool isFASTQ = false;              /* FASTQ sequence file */
	_Bool isEMBL = false;               /* SwissProt EMBL sequence file */
	const char * DBFileName = NULL;     /* Database pointer to use */
	_Bool UsingRegex = false;           /* Are we in regex mode or not? */
	_Bool IsAPattern = false;	          /* Is it a pattern rather than a regex */
	_Bool RegexOnlyHeader=false;        /* Trigger search within header only */
	_Bool RegexInFile = false;          /* Are regex provided on stdin or in file */
	char * RegexSource = NULL;          /* Pointer to the regex source argument */
	_Bool DoComplement = false;		 			/* Complementary mode to use*/
	size_t * shares = 0;

	struct ThreadData *threads_arg = NULL; /* Allocate stack memory for posix thread structures */
	enum Version HeuristicVersion = SSE2;  /* Trigger heuristic SSE version to use for heuristic */
	enum Version OtherVersion = SSE2;      /* Trigger heuristic SSE version to use for filter and alignment */
#if !defined(__USE_WINAPI__)
	pthread_t *threads = NULL;
#else
	HANDLE * threads = NULL;
#endif
#ifdef USE_AFFINITY
	Affinity_Mask_t * Thread_masks[2] = {0,0}; /* Define variables to hold thread affinity mask */
	unsigned int Thread_count[2] = {0,0};
	pthread_attr_t * restrict threads_attr = NULL;
	char buffer[128] __attribute__((aligned(16))); /* buffer to read affinity file mask */
	_Bool noAffinity = false;                      /* disable use of cpu affinity */
	_Bool noSharedCore = false;             /* Prevent hyperthreading or AMD compute unit to share resources */
	_Bool GivenAffinityFile = false;        /* File holding a mask for each thread */
	char * AffinityMaskFileName;            /* Name of affinity mask file provided by option m */
#endif
	_Bool HasGivenPrintMethod=false;

	void* (*thread_alignment)(void*) = thread_xaliPT;
	char UnknownSymbol = '\0';

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
		xalit_ptr = xalit_sse41;
		xalip_ptr = xalip_sse41;
		HeuristicVersion = SSE41;
		OtherVersion = SSE41;
	} else {
		xali1_ptr = xali1_sse2;
		xalit_ptr = xalit_sse2;
		xalip_ptr = xalip_sse2;
		HeuristicVersion = SSE2;
		OtherVersion = SSE2;
	}
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
			case '0':
				noAffinity = true;
				break;
			case '1':
				GivenAffinityFile = true;
				AffinityMaskFileName = optarg;
				break;
			case '2':
				noSharedCore = true;
				break;
#endif
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
				if (HeuristicCutOff != -1) {
					fprintf(stderr, "Choose either no heuristic or setting heuristic value\n");
					exit(1);
				}
				NoHeuristic = true;
				HeuristicCutOff = 0;
				break;
			case 'a':
				NoHeuristic = true;
				Optimal = true;
				thread_alignment = thread_optimal;
				break;
			case 'y':
				SatisfyCutoff = true;
				break;
			case 'p':
				DumpFilterOnly = true;
				break;
			case 'P':
				DumpHeuristicOnly = true;
				break;
			case 't':
				nCPUs = (size_t) atoi(optarg);
				break;
			case 'T':
				nCPUsHeuristic = (size_t) atoi(optarg);
				break;
			case 'd':
				DumpFilterSequences = true;
				break;
			case 'D':
				DumpAlignmentSequences = true;
				break;
			case 's':
				xali1_ptr = xali1_sse2;
				xalit_ptr = xalit_sse2;
				xalip_ptr = xalip_sse2;
				HeuristicVersion = SSE2;
				OtherVersion = SSE2;
				break;
			case 'i':
				ImportIndices = true;
				ImportFileName = optarg;
				break;
			case 'H':
				if (NoHeuristic) {
					fprintf(stderr, "Choose either no heuristic or setting heuristic value\n");
					exit(1);
				}
				{
					const char * ptr = optarg;
					_Bool isFloat = false;
					while(*ptr != '\0') if (*ptr == '.') {isFloat = true; break;} else ptr++;
					if (isFloat) {
						fputs("Heuristic cutoff is meant to be an unsigned integer!", stderr);
						exit(1);
					}
				}
				HeuristicCutOff = atoi(optarg);
				break;
			case 'C':
				NormalizedCutoff = false;
				FilterCutoff = (unsigned int) atoi(optarg);
				break;
			case 'N':
				NormalizedCutoff = true;
				FloatingFilterCutoff = (float) atof(optarg);
				break;
			case 'R':
				IsReversed = true;
				if (OutputVerbose) fputs("Working on reversed profile\n", stderr);
				break;
			case 'r':
				ReverseSequence = true;
				if (OutputVerbose) fputs("Working on reversed sequences\n", stderr);
				break;
			case 'u':
				MaxRegexNumber = (size_t) atol(optarg);
				break;
			case 'Q':
				RegexOnlyHeader = true;
				break;
			case 'x':
				DoComplement = true;
				break;
			case 'b':
				DoBothFwdAndRev = true;
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
			case 'U':
				UnknownSymbol = optarg[0];
				break;
			default:
				fprintf(stderr,"Option %c is unknown\n", c);
			case 'h':
				Usage(stdout);
		}
	}

	if (optind >= argc) {
		fputs("Expected arguments after options\n", stderr);
		Usage(stderr);
	}
	else {
		/* Test if Regex or profile */
		if (strncmp(argv[optind], "regex{", 6) == 0) {
			UsingRegex  = true;
			RegexSource = argv[optind];
		}
		else if (strncmp(argv[optind], "pattern{", 8) == 0) {
			IsAPattern  = true;
			UsingRegex  = true;
			RegexSource = argv[optind];
		}
		else {
			ProfileFile = argv[optind];
		}
		DBFileName = argv[optind+1];
	}

	if (OutputVerbose) {
		fputs(HEADER
#ifdef __USE_MMAP__
		"| Using Linux kernel MMAP function.                                  |\n"
#endif
		,stderr);
		printSystemInfo(stderr, &System);
#ifdef USE_32BIT_FORMAT
		fputs("Using 32 bit format integer for scores\n", stderr);
#endif
		if (OtherVersion == SSE2 && (System.Extensions & MM_SSE41)) {
			fputs("Enforcing SSE 2...\n", stderr);
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// COMMAND LINE ARGUMENTS COMPLIANCE
	////////////////////////////////////////////////////////////////////////////////////////////////
	{
		int count = 0;
		if (isFASTA) count++;
		if (isEMBL) count++;
		if (isFASTQ) count++;

		if (count > 1 || count <= 0) {
			fputs("Please provide a single sequence database, FASTA or FASTQ or ENA\n", stderr);
			exit(1);
		}
	}

	if (DoBothFwdAndRev) {
		if (ReverseSequence || DoComplement) {
			fputs("Do not specify --both and (--reverse-sequence and/or --complement\n", stderr);
			exit(1);
		}
		strand = BOTH;
		DoComplement = true;
	}
	else if (ReverseSequence && DoComplement) strand = REVERSE_COMPLEMENT;
	else if (ReverseSequence) strand = REVERSE;
	else if (DoComplement) strand = COMPLEMENT;
	else strand = FORWARD;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// INPUT ANALYSIS
	////////////////////////////////////////////////////////////////////////////////////////////////

	if (!UsingRegex) {
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
		const int ProfileCount = ReadProfile(ProfileFile, prf, !IsReversed);
		gettimeofday(&_t1,0);
		if (ProfileCount < 0) {
			fputs("Error found reading profile.\n", stderr);
			exit(1);
		}
		const double T = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;

		if (ProfileCount > 1 && (void*)PrintFunction != (void*)&PrintClassification) {
			fputs("pfsearchV3 is not meant to be used with multiple profiles other than in the case of classification\n"
			"Use pfscanV3 for that purpose.\n",stderr);
			exit(1);
		}
		else if (ProfileCount > 1) {
			if (OutputVerbose) fprintf(stderr, "Profiles reading took %lf seconds: %i found\n", T, ProfileCount);
			if (prf->Type == PF_MATRIX) {
				if (OutputVerbose) {
					fprintf(stderr,"Master classification profile %s has length %lu and alphabet size of %lu\n",
									ProfileFile, prf->Length, prf->Alphabet_Length);
					fputs("Alphabet Mapping\n",stderr);
					for (size_t i=0; i<ALPHABET_SIZE; ++i) {
						fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
						if ((i+1) % 8 == 0 ) fputs("\n",stderr);
					}
					fputs("\n\n",stderr);

					fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
					fprintf(stderr,"Description:\n"
					               "     Master: %s\n", prf->Description);
				}
				struct Profile * tmpPrf = prf;
				int i=1;
				while (tmpPrf->next) {
					tmpPrf = tmpPrf->next;
					if (tmpPrf->Type == PF_PATTERN) {
						fputs("There is some pattern profile within the given list of profiles, classification cannot be performed.\n", stderr);
						exit(1);
					}
					if (OutputVerbose) fprintf(stderr, "   Slave %2i: %s\n", i++, tmpPrf->Description);
				}
			}
			else {
				fputs("There is some pattern profile within the given list of profiles, classification cannot be performed.\n",	stderr);
				exit(1);
			}
		}
		else {
			if ((void*)PrintFunction == (void*)&PrintClassification) {
				fputs("Classification in not meant to be used with a single profile, the profile file should contain several.\n", stderr);
				exit(1);
			}

			/* Do we have a single pattern profile to treat */
			if (prf->Type == PF_PATTERN) {
				UsingRegex = true;
				IsAPattern = true;
				RegexSource = prf->Pattern;
				if (OutputVerbose) {
					fprintf(stderr, "Profile reading took %lf seconds: %i found\n"
					                "Profile pattern is %s\n", T, ProfileCount, prf->Pattern);
				}
			}
			else if (OutputVerbose) {
				fprintf(stderr, "Profile reading took %lf seconds: %i found\n", T, ProfileCount);
				if ((void*)PrintFunction == (void*)&PrintClassification ) {
					fputs("You asked for classification but the given file only contains one profile\n", stderr);
					return 1;
				}
				fprintf(stderr,"Profile %s has length %lu and alphabet size of %lu\n",
								ProfileFile, prf->Length, prf->Alphabet_Length);

				fputs("Alphabet Mapping\n",stderr);
				for (size_t i=0; i<ALPHABET_SIZE; ++i) {
					fprintf(stderr,"Map %c=%2u  ", (char) ((unsigned char) 'A' + (unsigned char) i), (unsigned int) prf->Alphabet_Mapping[i]);
					if ((i+1) % 8 == 0 ) fputs("\n",stderr);
				}
				fputs("\n",stderr);

				fprintf(stderr,"Disjoint set: %i to %i\n", prf->DisjointData.NDIP[0], prf->DisjointData.NDIP[1]);
			}
		}

		/* Set the Mode and Level provided we are on a matrix profile */
		if (prf->Type == PF_MATRIX) {
			if (IsReversed) {
				if (OutputVerbose) fputs("Reversing profile ...\n", stderr);
				struct Profile * const restrict rprf = ReverseProfile(prf);
				if (rprf == NULL) {
					fputs("Error while reversing profile\n", stderr);
					exit(1);
				} else {
					FreeProfile(prf, true);
					PrepareExtraTable(rprf);
					prf = rprf;
				}
			}

			/* Set the Cutoff Level and the Normalization */
			{
				register const int res = SetProfileLevelAndMode(prf, Level, Mode);
				if (res < 0) {
					fprintf(stderr, "Unable to find cutoff with level set to %i and mode set to %i, error %i\n", Level, Mode, res);
					exit(1);
				}
				if (OutputVerbose) {
					fprintf(stderr, "Using level %i, mode %i : %s with coefficients Rx=(", Level, Mode, NormalizationModeName[prf->ModeIndex]);
					register const int nCoefs = prf->NormalizationData.JNOP[prf->ModeIndex]-1;
					for (int c=0; c<nCoefs; ++c) {
						fprintf(stderr,"%lf, ", prf->NormalizationCoefs[c]);
					}
					fprintf(stderr,"%lf)\n", prf->NormalizationCoefs[nCoefs]);
				}
			}

			/* In case of given normalized cutoff, apply the normalization */
			if (NormalizedCutoff) {
				prf->NormalizedCutOff = FloatingFilterCutoff;

				/* Check the mode type since some cannot be global for all sequences */
				if (prf->NormalizationType != GLE_ZSCAVE) {
					FilterCutoff = prf->NormalizedToRaw(FloatingFilterCutoff, prf->NormalizationCoefs, 0.0f, 0);
					if (OutputVerbose)
						fprintf(stderr,"Translating normalized cutoff %lf to raw cutoff %i\n",FloatingFilterCutoff, FilterCutoff);
					const int hcut = ComputeHeuristicCutoff(prf, Mode, FilterCutoff);
					/* We do not want to impose the new heuristic cutoff but rather keep the command line option first */
					if (hcut > 0 && HeuristicCutOff == -1) {
						HeuristicCutOff = hcut;
						if (OutputVerbose)
							fprintf(stderr,"Translating filter cutoff %i to heuristic cutoff %i\n", FilterCutoff, HeuristicCutOff);
					}
				}
			}

			if (UnknownSymbol) {
				prf->CABC[0] = UnknownSymbol;
			}
		}
		else {
			if (IsReversed) {
				fputs("For the time being Pattern profile cannot be reversed.\n", stderr);
				exit(1);
			}
		}
	}
	else {
		/* Parse the regex argument to get the real regular expression */
		/* Create a copy of the source on the stack */
		const size_t regexSourceSize = 1+strlen(RegexSource);
		char * ctmp = (char *) alloca(regexSourceSize*sizeof(char));
		for (size_t i=0; i<regexSourceSize; ++i) ctmp[i] = RegexSource[i];
		ctmp[regexSourceSize] = '\0';

		/* Check format extract '{' '}' border */
		char * Start;
		char * ptr = ctmp;
		char * const End = &ctmp[regexSourceSize-1];
		_Bool ErrorinParsing = false;

		while ( *ptr != '{' && ptr <= End) ++ptr;
		if (ptr != End) {
			Start = ++ptr;
			ptr = End;
			while ( *ptr != '}' && ptr <= End) --ptr;
			if (ptr != Start)
				*ptr = '\0';
			else
				ErrorinParsing = true;
		}
		else
			ErrorinParsing = true;

		if (! ErrorinParsing)
			RegexSource = Start;
		else {
			fprintf(stderr, "Unable to retrieve correctly the regular expression or pattern from %s\n", RegexSource);
			exit(1);
		}
	}

	if (UsingRegex) {
		if (ReverseSequence) {
			fputs("For the time being Pattern profile cannot have reversed sequences.\n", stderr);
			exit(1);
		}
		if (DoComplement) {
			fputs("Sequences complementation is not possible on pattern\n", stderr);
			exit(1);
		}

		/* Initialize regex structure */
		const int res = InitRegExFromString(&regex, RegexSource, IsAPattern, MaxRegexNumber);
		if (res<0) {
			fprintf(stderr, "Regex initialization return error code %i\n", res);
			exit(1);
		}
	}

	/*
	 * Read the FASTA file
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

	if (OutputVerbose) fprintf(stderr, "Job dispatched over %lu cores.\n", nCPUs);

#ifdef __USE_MMAP__
	const int fd = open(DBFileName, O_RDONLY );
	const size_t length = (size_t) DB.FileSize;
	const char * const SequenceFileMap = (const char * const restrict) mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0);
	if (SequenceFileMap == MAP_FAILED) {
		fputs("Unable to map sequence file to memory\n", stderr);
		exit(1);
	}
	const char * const restrict SequenceFile = SequenceFileMap;
	close(fd);
#else
	const char * const restrict SequenceFile = DB;
#endif

	/* Prepare structure common to filter and alignment */
	shares = alloca((nCPUs+1)*sizeof(size_t));

	/* Allocate stack memory for posix thread structures */
	threads_arg = alloca(nCPUs*sizeof(struct ThreadData));
#if !defined(__USE_WINAPI__)
	threads = (pthread_t*) alloca(nCPUs*sizeof(pthread_t));
#else
	threads = (HANDLE*) alloca(nCPUs*sizeof(HANDLE));
#endif

	/* Dispatch common information to threads */
	for (size_t i=0; i<nCPUs; ++i) {
		threads_arg[i].prf          = prf;
		threads_arg[i].DB           = &DB;
		threads_arg[i].SequenceID   = NULL;
		threads_arg[i].strand       = strand;
		threads_arg[i].SequenceFile = SequenceFile;
		threads_arg[i].threadId     = i;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// REGEX MODE
	////////////////////////////////////////////////////////////////////////////////////////////////

	if (UsingRegex) {
		/* Initialize the print mutex */
#if !defined(__USE_WINAPI__)
		pthread_mutex_init(&PrintLock, NULL);
#else
		InitializeCriticalSection(&PrintLock);
#endif

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
		const unsigned long SearchInHeader = RegexOnlyHeader ? 1 : 0;

		gettimeofday(&_t0,0);
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].regex        = &regex;
			threads_arg[i].start        = shares[i];
			threads_arg[i].stop         = shares[i+1];
			threads_arg[i].counter      = SearchInHeader;
#if !defined(__USE_WINAPI__)
#ifdef USE_AFFINITY
			if (pthread_create (&threads[i],  &threads_attr[i], thread_regex,  (void*) &threads_arg[i]) != 0)
#else
				if (pthread_create (&threads[i],  NULL, thread_regex,  (void*) &threads_arg[i]) != 0)
#endif
#else
					threads[i] = CreateThread(NULL, 0, thread_regex, (void*) &threads_arg[i], 0, NULL);
				if (threads[i] == NULL)
#endif
				{
					fputs("Fail to create thread.\n", stderr);
					exit(0);
				}
		}

		size_t MatchCounter = 0;
#if !defined(__USE_WINAPI__)
		for (size_t i=0; i<nCPUs; ++i) {
			pthread_join(threads[i], NULL);
			MatchCounter += threads_arg[i].counter;
		}
#else
		WaitForMultipleObjects(nCPUs, threads, TRUE, INFINITE);
		for (size_t i=0; i<nCPUs; ++i) {
			MatchCounter += threads_arg[i].counter;
			CloseHandle(threads[i]);
		}
#endif
		gettimeofday(&_t1,0);
		double t;
		if (OutputVerbose) {
			t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Overall there are %lu matching string in %lu sequences. These took %lf seconds to treat on %lu cores.\n",
							MatchCounter, DB.SequenceCount, t, nCPUs);
		}

		/* Free the print mutex */
#if !defined(__USE_WINAPI__)
		pthread_mutex_destroy(&PrintLock);
#else
		DeleteCriticalSection(&PrintLock);
#endif

		FreeRegEx(&regex);
		goto END_NON_PRF;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// DUMP BOTH HEURISTIC AND CUTOFF
	////////////////////////////////////////////////////////////////////////////////////////////////

	/* Get Heuristic cutofstderr);f from command line */
	if (HeuristicCutOff > -1) prf->HeuristicCutOff = HeuristicCutOff;
	if (FilterCutoff) prf->CutOff = (int) FilterCutoff;

	if (OutputVerbose) fprintf(stderr,"Heuristic cutoff set to %u\nFilter cutoff set to %i\n",
		prf->HeuristicCutOff, prf->CutOff);

	DumpBoth = DumpHeuristicOnly && DumpFilterOnly;

	////////////////////////////////////////////////////////////////////////////////////////////////
	// HEURISTIC
	////////////////////////////////////////////////////////////////////////////////////////////////

	if (prf->HeuristicCutOff > 0 || DumpHeuristicOnly && !DumpFilterOnly) {
		/* Allocate memory for sequence heuristic to be done */
		const size_t HeuristicScoresSize = (strand == BOTH) ? 2*DB.SequenceCount : DB.SequenceCount;
		HeuristicScores = _mm_malloc( HeuristicScoresSize*sizeof(unsigned int), 16);
		if (HeuristicScores == NULL) {
			fputs("Cannot allocate memory.\n", stderr);
			goto END;
		}

		/* Compute Match Score Matrix transpose */
		gettimeofday(&_t0,0);
		TransposeMatrix TIMatch;
		if (HeuristicVersion == SSE41) {
			TIMatch.i = TransposeAndConvertMatchMatrix(&(prf->Scores), prf->Alphabet_Length, prf->Length);
		}
		else {
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

		gettimeofday(&_t0,0);
		for (size_t i=0; i<nCPUsHeuristic; ++i) {
			threads_arg[i].start           = shares[i];
			threads_arg[i].stop            = shares[i+1];
			threads_arg[i].HeuristicScores = HeuristicScores;
			threads_arg[i].TransposeMatch  = TIMatch;
			threads_arg[i].version         = HeuristicVersion;
#if !defined(__USE_WINAPI__)
#ifdef USE_AFFINITY
			if (pthread_create (&threads[i],  &threads_attr[i], thread_heuristic,  (void*) &threads_arg[i]) != 0)
#else
				if (pthread_create (&threads[i],  NULL, thread_heuristic,  (void*) &threads_arg[i]) != 0)
#endif
#else
					threads[i] = CreateThread(NULL, 0, thread_heuristic, (void*) &threads_arg[i], 0, NULL);
				if (threads[i] == NULL)
#endif
				{
					fputs("Fail to create thread.\n", stderr);
					exit(0);
				}
		}

#if !defined(__USE_WINAPI__)
		for (size_t i=0; i<nCPUsHeuristic; ++i) {
			pthread_join(threads[i], NULL);
		}
#else
		WaitForMultipleObjects(nCPUsHeuristic, threads, TRUE, INFINITE);
		for (size_t i=0; i<nCPUsHeuristic; ++i) {
			CloseHandle(threads[i]);
		}
#endif
		gettimeofday(&_t1,0);
		double t;
		if (OutputVerbose) {
			t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Heuristic took %lf seconds to treat on %li cores.\n", t, nCPUsHeuristic);
		}
		_mm_free(TIMatch.f);

		/* Do we go for dump only, then output and quit */
		if (DumpHeuristicOnly && !DumpBoth) {
			/* Allocate memory to hold sequence */
			SeqData.Data.Memory = malloc(DB.MaxSequenceSize*sizeof(unsigned char));
			if (SeqData.Data.Memory == NULL) {
				fputs("Thread Cannot allocate memory for sequence.\n", stderr);
				_mm_free(HeuristicScores); HeuristicScores = NULL;
				goto END;
			}

			/* Open sequence file*/

#ifndef __USE_MMAP__
			FILE* inSequence = fopen(DB, "r");
#endif

			if (strand != BOTH) {
				char StrandChar;
				switch(strand) {
					case FORWARD: StrandChar = '+'; break;
					case REVERSE: StrandChar = 'r'; break;
					case REVERSE_COMPLEMENT: StrandChar = '-'; break;
					case COMPLEMENT: StrandChar = 'c'; break;
					case BOTH: StrandChar = '?'; break;
				}
				for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
					/* Read sequence */
#ifndef __USE_MMAP__
					PFSeq = ReadSequenceIndex(&SeqData, inSequence, DB.DataPtr[iseq]);
#else
					PFSeq = MMAP_ReadSequenceIndex(&SeqData, SequenceFileMap, &(DB.DataPtr[iseq].Sequence), 0
# ifdef MMAP_DEBUG
					, 0, 0, length
# endif
					);
#endif
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping, 0);

					char * ptr = SeqData.Data.Header;
					while (*ptr != ' ' && *ptr != '\n') ptr++;
					if (*ptr == ' ') *ptr = '\0';

					/* Ouput results */
					fprintf(stdout, "%u\t%s\t%c\n", HeuristicScores[iseq], SeqData.Data.Header, StrandChar);
				}
			}
			else {
				for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
					/* Read sequence */
#ifndef __USE_MMAP__
					PFSeq = ReadSequenceIndex(&SeqData, inSequence, DB.DataPtr[iseq]);
#else
					PFSeq = MMAP_ReadSequenceIndex(&SeqData, SequenceFileMap, &(DB.DataPtr[iseq].Sequence), 0
# ifdef MMAP_DEBUG
					, 0, 0, length
# endif
					);
#endif
					/* Translate first sequence */
					PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping, 0);

					char * ptr = SeqData.Data.Header;
					while (*ptr != ' ' && *ptr != '\n') ptr++;
					if (*ptr == ' ') *ptr = '\0';

					/* Ouput results */
					fprintf(stdout, "%u\t%s\t%c\n", HeuristicScores[2*iseq], SeqData.Data.Header, '+');
					fprintf(stdout, "%u\t%s\t%c\n", HeuristicScores[2*iseq+1], SeqData.Data.Header, '-');

				}
			}
#ifndef __USE_MMAP__
			fclose(inSequence);
#endif
			_mm_free(HeuristicScores); HeuristicScores = NULL;
			free(SeqData.Data.Memory);
			goto END;
		}

		/* Allocate memory to hold the sequences indices */
		SequenceID = (unsigned int*) _mm_malloc(HeuristicScoresSize*sizeof(unsigned int), 16);
		if (SequenceID == NULL) {
			fputs("Cannot allocate memory for sequence indices.\n", stderr);
			goto END;
		}

		/* Gather the one that passed the heuristic */
		HeuristicCounter = 0UL;
		register const unsigned int lHeuristicCutOff = prf->HeuristicCutOff;
		if (strand != BOTH) {
			unsigned int SequenceModificationMask = 0U;
			if (ReverseSequence) SequenceModificationMask |= 0x40000000;
			if (DoComplement) SequenceModificationMask |= 0x80000000;
			for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
				if (HeuristicScores[iseq] >= lHeuristicCutOff) {
					SequenceID[HeuristicCounter] = SequenceModificationMask | ((unsigned int) iseq & 0x3FFFFFFF);
					++HeuristicCounter;
				}
			}
		}
		else {
			for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
				if (HeuristicScores[iseq] >= lHeuristicCutOff) {
					SequenceID[HeuristicCounter] = ((unsigned int) iseq & 0x3FFFFFFF);
					SequenceID[HeuristicCounter+1] = 0x70000000 | ((unsigned int) iseq & 0x3FFFFFFF);
					HeuristicCounter += 2;
				}
			}
		}
		if (OutputVerbose) {
			fprintf(stderr,"Overall there are %lu/%lu sequences passing heuritic. These took %lf seconds to treat on %lu cores.\n",
							HeuristicCounter, DB.SequenceCount, t, nCPUsHeuristic);
		}
		/* Print out the sequences passing heuristic cutoff */
		if (DumpFilterSequences) {
			/* Allocate memory to hold sequence */
			SeqData.Data.Memory = malloc(DB.MaxSequenceSize*sizeof(unsigned char));
			if (SeqData.Data.Memory == NULL) {
				fputs("Program cannot allocate memory for sequence.\n", stderr);
				goto END;
			}

#ifndef __USE_MMAP__
			FILE* inSequence = fopen(DB, "r");
#endif

			for (size_t iseq=0; iseq<HeuristicCounter; ++iseq) {
				/* Read sequence */
				const size_t sequence_index = SequenceID[iseq];
#ifndef __USE_MMAP__
			PFSeq = ReadSequenceIndex(&SeqData, inSequence,&(DB.DataPtr[sequence_index].Sequence);
#else
			PFSeq = MMAP_ReadSequenceIndex(&SeqData, SequenceFileMap, &(DB.DataPtr[sequence_index].Sequence), 0
#  ifdef MMAP_DEBUG
				,0, 0, length
#  endif
				);
#endif
				/* Ouput results */
				fprintf(stdout, "%s\n%s\n", SeqData.Data.Header, SeqData.ProfileData.ProfileIndex);
			}
#ifndef __USE_MMAP__
			fclose(inSequence);
#endif
			free(SeqData.Data.Memory);
			goto END;
		}
		if (DumpBoth) {
			HeuristicCounter = DB.SequenceCount;
			if (strand != BOTH) {
				for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
					SequenceID[iseq] = ((unsigned int) iseq & 0x3FFFFFFF);
				}
			}
			else {
				for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
					SequenceID[2*iseq  ] = ((unsigned int) iseq & 0x3FFFFFFF);
					SequenceID[2*iseq+1] = 0x70000000 | ((unsigned int) iseq & 0x3FFFFFFF);
				}
			}
		}
		else {
			_mm_free(HeuristicScores);
			HeuristicScores = NULL;
		}
	}
	else if (!Optimal) {
		if (!NoHeuristic) {
			fprintf(stderr, "Profile %s (%s) is not calibrated for heuristic use.\n", prf->AC_Number, prf->Description);
			goto END;
		}
		if (OutputVerbose) fputs("Bypassing heuristic computation...\n",stderr);

		/* Allocate memory to hold the sequences indices */
		HeuristicCounter = (strand == BOTH) ? 2*DB.SequenceCount : DB.SequenceCount;
		SequenceID = (unsigned int*) _mm_malloc(HeuristicCounter*sizeof(unsigned int), 16);
		if (SequenceID == NULL) {
			fputs("Cannot allocate memory for sequence indices.\n", stderr);
			goto END;
		}
		if (strand != BOTH) {
			unsigned int SequenceModificationMask = 0U;
			if (ReverseSequence) SequenceModificationMask |= 0x40000000;
			if (DoComplement) SequenceModificationMask |= 0x80000000;
			for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
				SequenceID[iseq] = SequenceModificationMask | ((unsigned int) iseq & 0x3FFFFFFF);
			}
		}
		else {
			for (size_t iseq=0; iseq<DB.SequenceCount; ++iseq) {
				SequenceID[2*iseq  ] = ((unsigned int) iseq & 0x3FFFFFFF);
				SequenceID[2*iseq+1] = 0xC0000000 | ((unsigned int) iseq & 0x3FFFFFFF);
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// FILTER
	////////////////////////////////////////////////////////////////////////////////////////////////
	if (!Optimal) {
		/* Allocate memory for the filter scores */
		FilterScores = _mm_malloc(HeuristicCounter*sizeof(int), 16);
		if (FilterScores == NULL) {
			fputs("Unable to allocate memory for the filter scores\n",stderr);
			exit(1);
		}

		/* Compute the new share for each thread */
		size_t SequenceShare = HeuristicCounter / nCPUs;
		SequenceShare += (HeuristicCounter % nCPUs) > (nCPUs-1) ? 1 : 0;
		shares[0] = 0;
		for (size_t i=1; i<nCPUs; ++i) shares[i] = i*SequenceShare;

		shares[nCPUs] = HeuristicCounter;

		/* Dispatch to threads */
		{
			const unsigned long realFilterScore = DumpFilterOnly ? 1L : 0L;
			gettimeofday(&_t0,0);
			for (size_t i=0; i<nCPUs; ++i) {
				threads_arg[i].start        = shares[i];
				threads_arg[i].stop         = shares[i+1];
				threads_arg[i].SequenceID   = SequenceID;
				threads_arg[i].FilterScores = FilterScores;
				threads_arg[i].counter      = realFilterScore;
#if !defined(__USE_WINAPI__)
				if (pthread_create (&threads[i],
#ifdef USE_AFFINITY
					&threads_attr[i],
#else
					NULL,
#endif
					thread_xali1,
					(void*) &threads_arg[i]) != 0)
#else
					threads[i] = CreateThread(NULL, 0, thread_xali1, (void*) &threads_arg[i], 0, NULL);
				if (threads[i] == NULL)
#endif
				{
					return 1;
				}
			}
		}
#if !defined(__USE_WINAPI__)
		for (size_t i=0; i<nCPUs; i++) {
			pthread_join(threads[i], NULL);
		}
#else
		WaitForMultipleObjects(nCPUs, threads, TRUE, INFINITE);
		for (size_t i=0; i<nCPUs; ++i) {
			CloseHandle(threads[i]);
		}
#endif

		gettimeofday(&_t1,0);

		if (DumpFilterOnly) {
			/* Allocate memory to hold sequence */
			SeqData.Data.Memory = malloc(DB.MaxSequenceSize*sizeof(unsigned char));
			if (SeqData.Data.Memory == NULL) {
				fputs("Pfsearch cannot allocate memory for sequence.\n", stderr);
				goto END;
			}

			/* Open sequence file*/
#ifndef __USE_MMAP__
			FILE* inSequence = fopen(DB, "r");
#endif

			for (size_t iseq=0; iseq<HeuristicCounter; ++iseq) {
				/* Read sequence */
				const unsigned int index = SequenceID[iseq] & 0x3FFFFFFF;
#  ifndef __USE_MMAP__
				PFSeq = ReadSequenceIndex(&SeqData, inSequence, &(DB.DataPtr[index].Sequence));
#  else
				PFSeq = MMAP_ReadSequenceIndex(&SeqData, SequenceFileMap, &(DB.DataPtr[index].Sequence), 0
#    ifdef MMAP_DEBUG
				, 0, 0, length
#    endif
				);
#  endif
				/* Translate first sequence */
				PFSeq = TranslateSequenceToIndex(PFSeq, prf->Alphabet_Mapping, 0);

				char * ptr = SeqData.Data.Header;
				while (*ptr != ' ' && *ptr != '\n') ptr++;
				if (*ptr == ' ') *ptr = '\0';

				/* Get Strand */
				char StrandChar;
				{
					const unsigned int flags = (SequenceID[iseq] & 0xC0000000);
					switch(flags) {
						case 0x00000000: StrandChar = '+'; break;
						case 0xC0000000: StrandChar = '-'; break;
						case 0x80000000: StrandChar = 'r'; break;
						case 0x40000000: StrandChar = 'c'; break;
					}
				}
				/* Ouput results */
				if (DumpBoth)
					fprintf(stdout, "%u\t%i\t%s\t%c\n", HeuristicScores[iseq], FilterScores[iseq], SeqData.Data.Header, StrandChar);
				else
					fprintf(stdout, "%i\t%s\t%c\n", FilterScores[iseq], SeqData.Data.Header, StrandChar);
			}

#ifndef __USE_MMAP__
			fclose(inSequence);
#endif
			free(SeqData.Data.Memory);
			goto END;
		}


		/* Gather the one that passed xali1 */
		FilterCounter = 0;
		register const int lFilterCutoff = prf->CutOff;
		for (size_t iseq=0; iseq<HeuristicCounter; ++iseq) {
			if ( FilterScores[iseq] >= lFilterCutoff ) {
				SequenceID[FilterCounter++] = SequenceID[iseq];
			}
		}
		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Overall there are %lu/%lu sequences passing filter. These took %lf seconds to treat on %li cores.\n",
							FilterCounter, HeuristicCounter, t, nCPUs);
		}

		/* Print out the sequences passing heuristic and filter cutoff */
		if (DumpAlignmentSequences) {
			/* Allocate memory to hold sequence */
			SeqData.Data.Memory = malloc(DB.MaxSequenceSize*sizeof(unsigned char));
			if (SeqData.Data.Memory == NULL) {
				fputs("Program cannot allocate memory for sequence.\n", stderr);
				goto END;
			}

#ifndef __USE_MMAP__
			FILE* inSequence = fopen(DB, "r");
#endif

			for (size_t iseq=0; iseq<FilterCounter; ++iseq) {
				/* Read sequence */
				const size_t sequence_index = SequenceID[iseq] & 0x3FFFFFFF;
#ifndef __USE_MMAP__
				PFSeq = ReadSequenceIndex(&SeqData, inSequence, &(DB.DataPtr[sequence_index].Sequence));
#else
				PFSeq = MMAP_ReadSequenceIndex(&SeqData, SequenceFileMap, &(DB.DataPtr[sequence_index].Sequence), 0
#  ifdef MMAP_DEBUG
				, 0, 0, length
#  endif
				);
#endif
				/* Ouput results */
				fprintf(stdout, "%s\n%s\n", SeqData.Data.Header, SeqData.ProfileData.ProfileIndex);
			}
#ifndef __USE_MMAP__
			fclose(inSequence);
#endif
			free(SeqData.Data.Memory);
			goto END;
		}
		_mm_free(FilterScores); FilterScores = NULL;
	}
	else {
		FilterCounter = DB.SequenceCount;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////
	// ALIGNMENT
	////////////////////////////////////////////////////////////////////////////////////////////////
	if (FilterCounter > 0) {
		/* Output header for some printing function */
		if (PrintFunction == &PrintSAM) {
			/* RNAME */
			const char * cptr = prf->AC_Number;
			while(*cptr != ';' && *cptr != '\0') { cptr++; }
			const int AClen = (int) ((uintptr_t) cptr - (uintptr_t) prf->AC_Number);
			printf("@HD\tVN:1.6\tSO:coordinate\n"
				     "@SQ\tSN:%.*s|%s\tLN:%zu\tDS:%s\n",
					   AClen, prf->AC_Number, prf->Identification, prf->Length, prf->Description);
		} else if (PrintFunction == &PrintTurtle) {
            printf("PREFIX ys:<http://example.org/yoursequence/>\n");
            printf("PREFIX yr:<http://example.org/yourrecord/>\n)");
            printf("PREFIX up:<http://purl.uniprot.org/core/>\n");
            printf("PREFIX rdf:<http://www.w3.org/1999/02/22-rdf-syntax-ns>\n");
            printf("PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>\n");
            printf("PREFIX faldo:<http://biohackathon.org/resource/faldo#>\n");
            printf("PREFIX signature:<htp://purl.uniprot.org/hamap/>\n");
            printf("PREFIX edam:<http://edamontology.org/>\n");
            printf("PREFIX hamap:<http://hamap.expasy.org/rdf/>\n");
        }


		/* Initialize the print mutex */
#if !defined(__USE_WINAPI__)
		pthread_mutex_init(&PrintLock, NULL);
#else
		InitializeCriticalSection(&PrintLock);
#endif
		/* Compute the new share for each thread */
		{
			size_t SequenceShare = FilterCounter / nCPUs;
			SequenceShare += (FilterCounter % nCPUs) > (nCPUs-1) ? 1 : 0;
			shares[0] = 0;
			for (size_t i=1; i<nCPUs; ++i) {
				shares[i] = i*SequenceShare;
				//  	fprintf(stderr,"share %lu starts at %lu and stops at %lu\n", i, shares[i-1], shares[i]);
			}
			shares[nCPUs] = FilterCounter;
			//       fprintf(stderr,"share %lu starts at %lu and stops at %lu\n", nCPUs, shares[nCPUs-1], shares[nCPUs]);
		}

		/* Dispatch to threads */
		gettimeofday(&_t0,0);
		for (size_t i=0; i<nCPUs; ++i) {
			threads_arg[i].start       = shares[i];
			threads_arg[i].stop        = shares[i+1];
#if !defined(__USE_WINAPI__)
#ifdef USE_AFFINITY
			if (pthread_create (&threads[i],  &threads_attr[i], thread_alignment,  (void*) &threads_arg[i]) != 0)
#else
				if (pthread_create (&threads[i],  NULL, thread_alignment,  (void*) &threads_arg[i]) != 0)
#endif
#else
				threads[i] = CreateThread(NULL, 0, thread_alignment, (void*) &threads_arg[i], 0, NULL);
				if (threads[i] == NULL)
#endif
				{
					return 1;
				}
		}
#if !defined(__USE_WINAPI__)
		for (size_t i=0; i<nCPUs; i++) {
			pthread_join(threads[i], NULL);
		}
#else
		WaitForMultipleObjects(nCPUs, threads, TRUE, INFINITE);
		for (size_t i=0; i<nCPUs; ++i) {
			CloseHandle(threads[i]);
		}
#endif
		gettimeofday(&_t1,0);

		unsigned int AlignedSequencesCounter = threads_arg[0].counter;
		for (size_t i=1; i<nCPUs; i++) AlignedSequencesCounter += threads_arg[i].counter;

		if (OutputVerbose) {
			const double t = (double) (_t1.tv_sec - _t0.tv_sec) + (double) (_t1.tv_usec - _t0.tv_usec) * 0.000001;
			fprintf(stderr,"Overall there are %u aligned sequences found. These took %lf seconds to align on %li cores.\n", AlignedSequencesCounter, t, nCPUs);
		}

		/* Free the print mutex */
#if !defined(__USE_WINAPI__)
		pthread_mutex_destroy(&PrintLock);
#else
		DeleteCriticalSection(&PrintLock);
#endif
	}

	/* Free Memory */
	END:
	if (SequenceID) _mm_free(SequenceID);
	if (FilterScores) _mm_free(FilterScores);
	if (HeuristicScores) _mm_free(HeuristicScores);
	FreeProfile(prf, true);
	END_NON_PRF:
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
