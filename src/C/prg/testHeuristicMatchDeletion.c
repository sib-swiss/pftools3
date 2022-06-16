#include <stdlib.h>
#include <stdio.h>
#include "../../include/pfProfile.h"

int main(int argc, char *argv[])
{
  char LINE[256];
  char *cptr;
  const char cMDDM[] = "MDDM ";
  const char cMIIM[] = "MIIM ";
  struct Profile prf;

  if (argc < 1) return 1;

  ReadProfile(argv[1], &prf);
  struct SMatch Match = prf.Scores.Match;

  const size_t prfLength = prf.Length + 1;
  const StoredIntegerFormat * Matches = Match.Alphabet;
  const size_t AlignedStep  = Match.AlignStep;
  const TransitionScores * restrict InsertionLine = prf.Scores.Insertion.Transitions;

  printf("Match matrix with alignment %lu\n\n",Match.AlignStep );
  printf("    | ");
  for (size_t alpha=0; alpha<prf.Alphabet_Length; ++alpha) {
    printf("%4lu ", alpha);
  }
  fputs("\n", stdout);

  printf("    | ");
  for (size_t alpha=0; alpha<prf.Alphabet_Length; ++alpha) {
    fputs("-----",stdout);
  }
  fputs("\n", stdout);

  for (size_t iprf=0; iprf<prf.Length; ++iprf) {
    const StoredIntegerFormat MDDM = InsertionLine[iprf].From[MATCH].To[DELETION]+InsertionLine[iprf+1].From[DELETION].To[MATCH];
    const StoredIntegerFormat MIIM = InsertionLine[iprf].From[MATCH].To[INSERTION]+InsertionLine[iprf+1].From[INSERTION].To[MATCH];
    StoredIntegerFormat Minimum;
    const char * STR;
    if (MIIM < MDDM) {
      Minimum = MDDM;
      STR = cMDDM;
    } else {
      Minimum = MIIM;
      STR = cMIIM;
    }
    const StoredIntegerFormat * MatchLine = &Matches[iprf*AlignedStep];
    memset(LINE, 0, 256*sizeof(char));
    cptr = LINE + strlen(LINE);
    printf("%3lu | ", iprf+1);
    for (size_t alpha=0; alpha<prf.Alphabet_Length; ++alpha) {
      if (MatchLine[alpha] < Minimum) {
	fputs(STR, stdout);
	sprintf(cptr, "%i,%i ", MatchLine[alpha], Minimum);
	cptr = LINE + strlen(LINE);
      } else if (MatchLine[alpha] == NLOW) {
        printf("NLOW ");
      } else {
        printf("%4i ", MatchLine[alpha]);
      }
    }
    printf("\t%s\n", LINE);
  }
  FreeProfile(&prf);
  return 0;
}
