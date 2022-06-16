/*******************************************************
                        PFTOOLS
 *******************************************************
  Jun 12, 2013 pam.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"
#include <stdlib.h>
#include <string.h>
#include "random.h"
#include "../include/pfProfile.h"
#include "../include/pfSequence.h"

/* WARNING: outSequence has to be allocated using malloc */
void mutate(struct Random * const restrict Generator, const struct Profile * const restrict prf,
	    const PFSequence * const restrict inSequence, Sequence * const restrict outSequence,
	    const size_t MutationNumber)
{
  /*
   * Keep in mind that the corresponding indices used to translate the alphabet have 0 as unknown
   * Furthermore Alphabet length is the dimension size, not the real size, therfore we look for the
   * second occurance of 'X'.
   * At last, we expect outSequence to be already allocated (malloc) to at least the size of inSequence.
   * From then it is adjusted.
   */

	size_t AlphabetLength = 1;
	while ( prf->CABC[AlphabetLength] != 'X' )  ++AlphabetLength;
	AlphabetLength -= 1;

	const float fAlphabetLength = (float) AlphabetLength;
	size_t SequenceLength       = inSequence->Length;
	const size_t N              = MutationNumber;
	size_t MemorySize           = outSequence->Size;
	unsigned char * Memory      = outSequence->Data.Memory;

	memcpy(Memory, inSequence->ProfileIndex, inSequence->Length*sizeof(unsigned char));

	for (size_t i=0; i<N; ++i) {
		/* pick up ramdomly the possible transformations */
		const float NewTransform = (fAlphabetLength+1.0f)*NormalizedFlatDistributionValue(Generator);
		const float f = NormalizedFlatDistributionValue(Generator);
		const size_t position = (size_t) (((float)SequenceLength)*f);

		if (NewTransform >= 1.0f) {
			/* Single mutation */
			Memory[position] = (unsigned char) NewTransform;
		}
		else if (NewTransform >= 0.5f) {
			/* Deletion */
			if (SequenceLength > 1) {
				SequenceLength -= 1;
				for (size_t p=position; p<SequenceLength; ++p) Memory[p] = Memory[p+1];
			}
		}
		else {
			/* Insertion */
			char Buffer[32] __attribute__((aligned(16)));
			size_t count;
			do {
				const float f = NormalizedFlatDistributionValue(Generator);
				float gd_param = f > 0.2f ? f : 0.21f;
				count = 0;
				while (gd_param > 0.2f) { gd_param *= gd_param; ++count; }
			} while ( count > 32);

			/* Check if memory is ok */
			{
				const size_t tmp = SequenceLength + count;
				if (tmp >= MemorySize) {
					MemorySize = 10+tmp;
					unsigned char * const ptr = (unsigned char*) malloc(MemorySize*sizeof(unsigned char));
					if (ptr == NULL) {
						fputs("PAM mutation: unable to extend sequence size\n", stderr);
						exit(1);
					}
					char * restrict ptr2 = ptr;
					register size_t p;
					for (p=0; p<position; ++p) *ptr2++ = (char) Memory[p];
					for (p=0; p<count; ++p) *ptr2++ = 1 + (unsigned char) (fAlphabetLength*NormalizedFlatDistributionValue(Generator));
					for (p=position; p<tmp; ++p) *ptr2++ = (char) Memory[p];
					ptr[tmp] = '\0';
					free(Memory);
					Memory = ptr;
				}
				else {
					register const size_t end = position+count;
					memmove(&Memory[end], &Memory[position], (SequenceLength - position)*sizeof(char));
					for (size_t p=position; p<end; ++p) Memory[p] = 1 + (unsigned char) (fAlphabetLength*NormalizedFlatDistributionValue(Generator));
				}
				SequenceLength += count;
			}
		}
		if (SequenceLength < 2) break;
	}
	outSequence->Size                     = MemorySize;
	outSequence->Data.Memory              = Memory;
	outSequence->ProfileData.ProfileIndex = Memory;
	outSequence->ProfileData.Length       = (size_t) SequenceLength;
}

#ifdef _PFPAM
#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>

static void __attribute__((noreturn)) Usage(FILE * stream)
{
  fputs(
    "Generate PAM distance matrices from a given sequence, used for profile calibration\n"
		" pfpam [options] profile sequence\n"
		" Options:\n"
		"   -N <uint>   : number of mutation (default 10)\n"
		"   -s <ulong>  : random generator seed (default 1234)\n"
		"   -v          : verbose mode\n"
		"   -h          : prints this help\n\n",
	stream);
  exit(0);
}

int main(int argc, char *argv[])
{
  struct Profile prf;
  struct Random Generator;
  char * ProfileFile;
  const char opt_to_test[] = "N:vs:h";
  size_t N=10;
  bool OutputVerbose = false;
  PFSequence pfseq_in;
  Sequence SeqOut;
  unsigned long Seed = 1234L;

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

  char * Sequence = argv[optind+1];

  /*
   * Read the profile and output some infos
   */
  int res = ReadProfile(ProfileFile, &prf, false);
  if (res <= 0) {
    fprintf(stderr, "Error found reading profile %s, code %i.\n", ProfileFile, res);
    return 1;
  }
  if (OutputVerbose) {
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

  size_t length = strlen(Sequence);
  pfseq_in.ProfileIndex = malloc(1+length);
  pfseq_in.Length = length;
  strcpy(pfseq_in.ProfileIndex, Sequence);
  printf(" Original sequence      : %s\n", pfseq_in.ProfileIndex);
  TranslateSequenceToIndex(&pfseq_in, prf.Alphabet_Mapping, 0);

  InitializeGenerator(Seed, &Generator);

  SeqOut.Data.Memory = malloc(1+length);
  SeqOut.Size        = (1+length);
  mutate(&Generator, &prf, &pfseq_in, &SeqOut,N);
  PFSequence * restrict pfseq_out = &SeqOut.ProfileData;

  for (size_t i=0; i<pfseq_out->Length; ++i)
		pfseq_out->ProfileIndex[i] = prf.CABC[pfseq_out->ProfileIndex[i]];
  pfseq_out->ProfileIndex[pfseq_out->Length] = '\0';
  printf(" New sequence           : %s\n", pfseq_out->ProfileIndex);

}
#endif /* _PFPAM*/
