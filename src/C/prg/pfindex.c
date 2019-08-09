// 	if (ExportIndices) {
// 		FILE *io = fopen(ExportFileName, "wb");
// 		if ( io != NULL ) {
// 			const int itmp = ExportStructure(io, &FASTA);
// 			if (OutputVerbose) {
// 				fprintf(stderr, itmp>0 ? "Export of indices failed, check space for %s\n": "Export of indices to file %s\n", ExportFileName);
// 			}
// 			fclose(io);
// 		} else {
// 			if (OutputVerbose)
// 				fprintf(stderr, "Export of indices failed, check write permission for %s\n", ExportFileName);
// 		}
// 	}
// 	

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
#include <stdbool.h>
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

#define HEADER "%----------------------------------------------------------------------------%\n"\
	       "|                              PFSEARCH v" PF_VERSION "                                |\n"\
	       "%----------------------------------------------------------------------------%\n"\
	       "| Built on " __DATE__ " at " __TIME__ ".                                          |\n"
#define __USE_INLINE_FUNCTIONS__
#include "../include/pfSequence.h"

static const char opt_to_test[] = "hVfqFo:";
static const struct option long_options[] =
{
	/*
	 * These options set a flag.
	 */
	
	
	/*
	 * These options don't set a flag. We distinguish them by their indices.
	 */
	/* Database indexing options */
	{"fasta",                     no_argument,				0,	'f'},
	{"fastq",                     no_argument,				0,	'q'},
	{"embl",                      no_argument,  			0,  'F'},
	/* Others */
	{"help",                      no_argument,				0,	'h'},
	{"verbose",                   no_argument,				0,	'V'},
	{"output", 										required_argument,	0,	'o'},
	{0, 0, 0, 0}
};

static void __attribute__((noreturn)) Usage(FILE * stream)
{
	fputs(
		"Index a sequence library for profile matches:\n"
		" pfindex [options] database\n\n"
		" Options:\n"
		"  Database\n"
		"   --fasta                            [-f] : FASTA file database as input\n"
		"   --fastq                            [-q] : FASTQ file database as input\n"
		"   --embl                                  : SwissProt/EMBL file database as input, default\n\n"
		"  Output\n"
		"   --output <file>                    [-o] : output index file name\n\n"
		"  Other\n"
		"   --verbose                          [-V] : verbose on stderr\n"
		"   --help                             [-h] : output command help\n\n"
		" Version " PF_VERSION " built on " __DATE__ " at " __TIME__ ".\n", stream);
	exit(0);
}

int main (int argc, char *argv[])
{
	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL STRUCTURES
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	DBSequence_t DB;                  /* Sequence Database File */
	struct timeval _t0, _t1;          /* Timing structures */
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// LOCAL DATA
	////////////////////////////////////////////////////////////////////////////////////////////////
	const char * DBFileName = NULL;     /* Database pointer to use */
	const char * OutFileName = NULL;    /* Output index file name */
	_Bool isFASTA = false;							/* FASTA sequence file */
	_Bool isFASTQ = false;              /* FASTQ sequence file */
	_Bool isEMBL = true;                  /* SwissProt EMBL sequence file */
	_Bool OutputVerbose = false;
		
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
			case 'V':
				OutputVerbose = true;
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
			case 'o':
				OutFileName = optarg;
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
		DBFileName = argv[optind];
	}
	
	
	////////////////////////////////////////////////////////////////////////////////////////////////
	// COMMAND LINE ARGUMENTS COMPLIANCE
	////////////////////////////////////////////////////////////////////////////////////////////////
	if (OutFileName == NULL) {
		fputs("Please specitfy an output index file name\n", stderr);
		exit(1);
	}
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
	////////////////////////////////////////////////////////////////////////////////////////////////
	// INPUT ANALYSIS
	////////////////////////////////////////////////////////////////////////////////////////////////
	
	/*
	 * Read the FASTA file
	 */
	int res;
	if (isFASTA) {
		gettimeofday(&_t0,0);

		res = AnalyzeFASTAStructure(DBFileName, &DB);
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
		res = AnalyzeFASTQStructure(DBFileName, &DB);
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
		res = AnalyzeEMBLStructure(DBFileName, &DB);
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

	FILE *io = fopen(OutFileName, "wb");
	if ( io != NULL ) {
		const int itmp = ExportDBStructure(io, &DB);
		if (OutputVerbose) {
			fprintf(stderr, itmp>0 ? "Export of indices failed, check space for %s\n": "Export of indices to file %s\n", OutFileName);
		}
		fclose(io);
	}
	else {
		if (OutputVerbose)
			fprintf(stderr, "Export of indices failed, check write permission for %s\n", OutFileName);
	}
	
	FreeDBStructure(&DB);
	
	exit(0);
}
