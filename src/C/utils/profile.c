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
#include "../include/pfProfile.h"
#include "../include/pfSequence.h"

int RescoreAlignment(const struct Profile * const prf, char * const AlignedSequence,
		     const struct Alignment * const alignment, const size_t SequenceLength)
{
  int Score = 0;
  int iprf=1;
  enum VectorPosition PreviousState;
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Prologue
  //////////////////////////////////////////////////////////////////////////////////////////////
  if (alignment->IPMB != 1 && alignment->IPME != -1) {
      fputs("Currently we do not provide classification for other than semiglobal alignment.\n", stderr);
      exit(1);
  }
  
  const size_t AlignmentLength            = strlen(AlignedSequence);
  const char * restrict Seq               = AlignedSequence;
  const char * restrict const MaxSequence = &AlignedSequence[AlignmentLength];
  unsigned int SequenceIndex              = alignment->Region.Sequence.Begin;
  const unsigned int SequenceBegin        = alignment->Region.Sequence.Begin;
  const unsigned int SequenceEnd          = alignment->Region.Sequence.End;
  
  /* Get first alignment state */
  if (*Seq == '-') 
    PreviousState = DELETION;
  else
    PreviousState = ( *Seq < 'a' ) ? MATCH : INSERTION;
  
  const TransitionScores * restrict pTransitions = prf->Scores.Insertion.Transitions;
  const StoredIntegerFormat * pMatch             = prf->Scores.Match.Alphabet;
  const StoredIntegerFormat * pInsertion         = prf->Scores.Insertion.Alphabet;
  const size_t AlignStep                         = prf->Scores.Match.AlignStep;
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // First Sequence Line
  ////////////////////////////////////////////////////////////////////////////////////////////// 
  if (SequenceBegin == 1) {
    register const ScoreTuple * restrict FirstSequenceProtein = prf->Scores.Insertion.FirstSequenceProtein;
    int lScore = (int) FirstSequenceProtein->To[PreviousState];
//     printf("Alignment starts at the beginning of the sequence\n");    
    /* We need to further keep looking that we are not in the case of a multiple deletion entrance */
    if (PreviousState == DELETION) {
      pInsertion += AlignStep;
      ++pTransitions;
      ++iprf;
      while (*(++Seq) == '-') {
// 	fputc(Seq[-1], stdout);
	lScore += ((int) pMatch[_D]) + ((int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(DELETION, DELETION)]);
	++pTransitions;
	pInsertion += AlignStep;
	pMatch     += AlignStep;
	++iprf;
      }
    }
    else if (PreviousState == MATCH) {
//        pInsertion += AlignStep;
//        ++pTransitions;
    }
    Score += lScore;
//     printf(" Profile pos: %4i Score: %i\n", iprf, Score);
  }
  else {
    Score += (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(EXTRA, PreviousState)];
//     ++pTransitions;
//     pInsertion += AlignStep;
//     printf(" Profile pos: %4i Score: %i\n", iprf, Score);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Loop through the internal sequence indices
  //////////////////////////////////////////////////////////////////////////////////////////////
  
  if (Seq < MaxSequence) { 
    do {
//       fputc(*Seq, stdout);
      if (*Seq == '-') {
	 Score += (int) pMatch[_D] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, DELETION)]; ;
	 PreviousState = DELETION;
	 ++pTransitions;
	 pInsertion += AlignStep;
	 pMatch     += AlignStep;
	 ++iprf;
      }
      else {
	const size_t index = (size_t) TranslateCharToIndex(*Seq, prf->Alphabet_Mapping);
	if ( *Seq < 'a' ) {
	  Score   += (int) pMatch[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, MATCH)];
	  PreviousState = MATCH;
	  ++pTransitions;
	  pInsertion += AlignStep;
	  pMatch     += AlignStep;
	  ++iprf;
	}
	else {
	  Score    += (int) pInsertion[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, INSERTION)];
// 	  printf(" Profile pos: %4i ?I: %i i: %i\n", iprf, pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, INSERTION)], pInsertion[index]);
	  PreviousState  = INSERTION;
	}
	++SequenceIndex;
      }
//       printf(" Profile pos: %4i Score: %i\n", iprf, Score);
    } while ( ++Seq < MaxSequence && SequenceIndex < SequenceLength);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  // Epilogue
  //////////////////////////////////////////////////////////////////////////////////////////////
  if (SequenceEnd == SequenceLength) {
    fputs(": End on sequence end", stdout);
  }
  else {
//     fputc(*Seq, stdout);
    if (*Seq == '-') {
      Score += (int) pMatch[_D] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)]; ;
      PreviousState = DELETION;
      ++pTransitions;
      pInsertion += AlignStep;
      pMatch     += AlignStep;
    }
    else {
      const size_t index = (size_t) TranslateCharToIndex(*Seq, prf->Alphabet_Mapping);
      if ( *Seq < 'a' ) {
	Score   += (int) pMatch[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)];
	PreviousState = MATCH;
	++pTransitions;
	pInsertion += AlignStep;
	pMatch     += AlignStep;
      }
      else {
	Score    += (int) pInsertion[index] + (int) pTransitions->Element[TRANSITION_INDEX_FROM_TO(PreviousState, EXTRA)];
	PreviousState  = INSERTION;
      }
    }
  }
  return Score;
}

/*
 * Treats the level of cutoff by seeking the index corresponding to the wanted level number,
 * then storing it into the profile Level member.
 * BE CAREFUL THIS IS THE INDEX AND NOT THE VALUE OF THE LEVEL !!!
 * 
 * The function shall the index to the level if found, -1 if not found and -2 upon multiple
 * existance of the same level.
 */
int SeekProfileLevel(const struct Profile * const prf, const int Level)
{
  int index = -1;
  const SCutOffItem * const restrict cutItems = prf->CutOffData.Values;
  _Bool levelfound = false;
  const size_t N = prf->CutOffData.JCUT;
  for (size_t icut=0; icut<N; ++icut) {
    if ( cutItems[icut].MCLE == Level ) {
      index = (int) icut;
      if (levelfound) return -2;
      levelfound = true;
    }
  }
  return index;
}
/*
 * This function returns the index of the highest priority normalization 
 * matching the given mode or -1 if not found.
 */
int SeekProfileMode(const struct Profile * const prf, const int Mode)
{
  int index = -1;
  register const int SeekMode = Mode;
  const SNormalizationItem * const restrict NormItems = &(prf->NormalizationData.Values[0]);
  
  int Priority =  -1;
  
  for (int iNormalizationMode=0; iNormalizationMode<prf->NormalizationData.JNOR; ++iNormalizationMode) {
    if ( SeekMode == NormItems[iNormalizationMode].NNOR ) {
      if (index < 0) {
	Priority = NormItems[iNormalizationMode].NNPR;
	index    = iNormalizationMode;
      } 
      else if (NormItems[iNormalizationMode].NNPR < Priority) {
	Priority = NormItems[iNormalizationMode].NNPR;
	index    = iNormalizationMode;
      }
    }
  }
  return index;
}

int SetProfileMode(struct Profile * const prf, const int Mode)
{
    const int index = SeekProfileMode(prf, Mode);
    if (index >= 0) {
      prf->ModeIndex = (short int) index;
      prf->NormalizationCoefs              = prf->NormalizationData.Values[index].RNOP;
      register const int NormalizationType = prf->NormalizationData.Values[index].MNOR;
      prf->NormalizationType               = NormalizationType;
      switch (NormalizationType) {
	case LINEAR:
	  prf->NormalizedToRaw = &N2R_1;
	  prf->RawToNormalized = &R2N_1;
	  break;
	case GLE_ZSCORE:
	  prf->NormalizedToRaw = &N2R_2;
	  prf->RawToNormalized = &R2N_2;
	  break;
	case GLE_ZSCAVE:
	  prf->NormalizedToRaw = &N2R_3;
	  prf->RawToNormalized = &R2N_3;
	  if (InitAverage(prf)<0) return -5;
	  break;
	default:
	  return -4;
      }
    }
    else {
      prf->ModeIndex          = (short int) -1;
      prf->NormalizationCoefs = NULL;
      prf->NormalizationType  = 0;;
      prf->NormalizedToRaw    = NULL;
      prf->RawToNormalized    = NULL;
    }
    const int index2 = SeekProfileMode(prf, -Mode);
    prf->HeuristicModeIndex = (index2 >= 0) ? (short int) index : -1;
    
    return index;
}

void SeekProfileLevelAndMode(int * const restrict LevelIndex, int * const restrict ModeIndexWithinLevel, int * const restrict ModeIndex,
			     int * const restrict HeuristicModeIndex, const struct Profile * const prf, const int Level, const int Mode)
{
  int localLevelIndex = -1;
  int localModeIndexWithinLevel = -1;
  int localHeuristicModeIndex = -1;
  
  /* First test if Mode exists */
  int localModeIndex = SeekProfileMode(prf, Mode);
  if (localModeIndex < 0 ) goto FIN;
  localHeuristicModeIndex = SeekProfileMode(prf, -Mode);
  
  /* Now scan the level to see if there exists one satisfying Level and Mode */
  if ( Level < MAXC ) {
    {
      const SCutOffItem * const restrict cutItem = prf->CutOffData.Values;
      register const int M = prf->CutOffData.JCUT;
      _Bool Levelfound = false;
      for (int iCutoff=0; iCutoff<M; ++iCutoff) {
				if (cutItem[iCutoff].MCLE == Level ) {
					localLevelIndex = iCutoff;
					if (Levelfound) {
						localLevelIndex = -2;
						goto FIN;
					}
					Levelfound = true;
				}
      }
      if (!Levelfound) {
				localLevelIndex = -1;
				goto FIN;
      }
    }
    
    const SCutOffItem * const restrict cutItem = &prf->CutOffData.Values[localLevelIndex];
    register const int N = cutItem->JCNM;
    for (int iCutoffMode=0; iCutoffMode<N; ++iCutoffMode) {
      const int CutoffMode = cutItem->MCUT[iCutoffMode];
      if (CutoffMode == Mode) {
				if (localModeIndexWithinLevel < 0) 
					localModeIndexWithinLevel = iCutoffMode;
				else {
					localModeIndexWithinLevel = -2;
					goto FIN;
				}
      }
    }
  }
  FIN:
  *LevelIndex           = localLevelIndex;
  *ModeIndexWithinLevel = localModeIndexWithinLevel;
  *ModeIndex            = localModeIndex;
  *HeuristicModeIndex   = localHeuristicModeIndex;
  return;
}
/*
 * SetProfileLevelAndMod returns:
 * --------------------------------
 *  0 : all OK
 * -1 : Level not found
 * -2 : Mutliple levels satisfy conditions
 * -3 : Mode not found
 * -4 : Unknown Mode in Level (might never exist, need to check)
 * -5 : Average initialization failed
 */
int SetProfileLevelAndMode(struct Profile * const prf, const int Level, const int Mode)
{
  int lLevel=-1, lModeWithinLevel=-1, lMode=-1, lHeuristicMode=-1;
  SeekProfileLevelAndMode(&lLevel, &lModeWithinLevel, &lMode, &lHeuristicMode, prf, Level, Mode);
  
  prf->HeuristicCutOff      = 0;
  prf->ModeIndex            = -1;
  prf->ModeIndexWithinLevel = -1;
  prf->HeuristicModeIndex   = -1;
  prf->LevelIndex           = -1;
  prf->NormalizedCutOff     = 0.0f;
  prf->NormalizationCoefs   = NULL;
  prf->NormalizationType    = 0;
  prf->NormalizedToRaw      = NULL;
  prf->RawToNormalized      = NULL;
  
  if(lMode < 0) return -3;
  
  if (lLevel >= 0) {
    if (lMode >= 0) {
      prf->LevelIndex                      = (short int) lLevel;
      prf->CutOff                          = prf->CutOffData.Values[lLevel].ICUT;
      prf->HeuristicCutOff                 = prf->CutOffData.Values[lLevel].HCUT;
      prf->ModeIndex                       = (short int) lMode;
      prf->ModeIndexWithinLevel            = (short int) lModeWithinLevel;
      prf->NormalizedCutOff                = prf->CutOffData.Values[lLevel].RCUT[lMode];
      prf->NormalizationCoefs              = prf->NormalizationData.Values[lMode].RNOP;
      register const int NormalizationType = prf->NormalizationData.Values[lMode].MNOR;
      prf->NormalizationType               = NormalizationType;
      switch (NormalizationType) {
				case LINEAR:
					prf->NormalizedToRaw = &N2R_1;
					prf->RawToNormalized = &R2N_1;
					break;
				case GLE_ZSCORE:
					prf->NormalizedToRaw = &N2R_2;
					prf->RawToNormalized = &R2N_2;
					break;
				case GLE_ZSCAVE:
					prf->NormalizedToRaw = &N2R_3;
					prf->RawToNormalized = &R2N_3;
					if (InitAverage(prf)<0) return -5;
					break;
				default:
					return -4;
      }
      if (lHeuristicMode>=0) prf->HeuristicModeIndex = lHeuristicMode;
    }
    else {
      return lMode - 2; 
    }
  }
  else {
    return lLevel;
  }
  return 0;
}
