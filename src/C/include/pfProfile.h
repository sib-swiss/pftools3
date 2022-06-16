/***************************************************************************************************
                        PFTOOLS
 ***************************************************************************************************
  Oct 3, 2011 pfProfile.h
 ***************************************************************************************************
 (C) 2011 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 ***************************************************************************************************/
#ifndef _PROFILE_H
#define _PROFILE_H
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <sys/types.h>
#include <string.h>
#include <stdio.h>
#include <emmintrin.h>
#include "pfConfig.h"

#define ALPHABET       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
#define ALPHABET_SIZE  26
/* Extra space for :
 *	undefined letter -> index 0
 * 	deletion cost    -> ALPHABET_SIZE+1
 *	stop codon       -> ALPHABET_SIZE+2
 */
#define ALPHABET_EXTRA_LETTERS  3
#define ALPHABET_MEMORY_SIZE   (ALPHABET_SIZE+ALPHABET_EXTRA_LETTERS)
#define TRANSITION_SIZE 	      9
#define PROFILE_MAX_LINE_SIZE   512

#define KDIS                    2
#define KDPM                    2
/* Maximum number of normalization modes, half is for score and the other for heuristic score */
#define MAXN                    8
/* Maximum number of function type in normalization (LINEAR,...) */
#define KNOR                    3
/* Maximum number of coefficients used in normalization function */
#define KNPM                    5
/* Maximum number of cutoffs */
#define MAXC                    8
/* Maximum number of alignments per sequence */
#define NALI 1024

#ifndef _BEST_IS_NEGATIVE_
/*
 * Lowest score value used when forbidding paths in the computation
 * Remember that we need to be able to go underneath without underflow
 */
# ifndef __USE_32BIT_INTEGER__
#  define NLOW              -16383
#  define STORED_INT_MIN	SHRT_MIN
#  define STORED_INT_MAX	SHRT_MAX
# else
// Division by 4 is required to prevent the pfplot tagging to underflow
#  define NLOW             -536870912/4
#  define STORED_INT_MIN        INT_MIN
#  define STORED_INT_MAX	      INT_MAX
# endif
# define NLOW_16                 -16383
#else
/*
 * Highest score value used when forbidding paths in the computation
 * Remember that we need to be able to go over without underflow
 */
# ifndef __USE_32BIT_INTEGER__
#  define NLOW               16383
#  define STORED_INT_MIN	SHRT_MIN
#  define STORED_INT_MAX	SHRT_MAX
# else
// Division by 4 is required to prevent the pfplot tagging to underflow
#  define NLOW              536870912/4
#  define STORED_INT_MIN	      INT_MIN
#  define STORED_INT_MAX	      INT_MAX
# endif
# define NLOW_16                  16383
#endif

/************************** Insertion Matrices ************************
 * There are 3 different matrices for:
 *    - alphabet
 *    - boundaries
 *    - transitions
 */
#define _UNKNOWN 	              0
#define _D   	 	  ALPHABET_SIZE+1
#define _STOP 	 	ALPHABET_SIZE+2
/*
 ************************************************************************************************
 *                                    BOUNDARIES                                                *
 ************************************************************************************************
 */
/* Price to pay to enter the profile alignment in First Sequence Protein */
#define _B0  0
/* Price to pay to enter the profile alignment elsewhere in the sequence */
#define _B1  1
/* Price to pay to leave the profile alignment in Last Sequence Protein */
#define _E0  2
/* Price to pay to leave the profile alignment elsewhere in the sequence */
#define _E1  3

/* Price to pay to enter the profile alignment used in both First Sequence Protein and Extra*/
/* From entrance to Match */
#define _BM  4
/* From entrance to Insertion */
#define _BI  5
/* From entrance to Deletion */
#define _BD  6
/* From entrance to ??? THIS SEEMS TO BE UNUSED*/
#define _BE  7

/* Extensions to pay to end the profile alignment used in both LastSequence and Extra*/
/* From Match to exit */
#define _ME  8
/* From Insertion to exit */
#define _IE  9
/* From Deletion to exit */
#define _DE 10

#define INSERTION_BOUNDARIES_SIZE 11

/*
 ************************************************************************************************
 *                                   TRANSITIONS                                                *
 ************************************************************************************************
 */
 /* TABLE OF TRANSITIONS
  *
  * CODE    | x < NDIP1  | NDIP1 <= i < NDIP2  | NDIP2 <= i  | part of
  * =============================================================================================
  * _XM     | B1   + BM  |       B1 + BM       |  NLOW + BM  |
  * _XI     | B1   + BI  |       B1 + BI       |  NLOW + BI  | Insertion score EXTRA
  * _XD     | B1   + BD  |       B1 + BD       |  NLOW + BD  |
  *         |            |                     |             |
  * _MX     | NLOW + ME  |       E1 + ME       |   E1  + ME  | Insertion score MATCH
  * _IX     | NLOW + IE  |       E1 + IE       |   E1  + IE  | Insertion score INSERTION
  * _DX     | NLOW + DE  |       E1 + DE       |   E1  + DE  | Insertion score DELETION
  *---------+------------+---------------------+-------------+-----------------------------------
  * _YM     | B0 + BM    |       B0 + BM       |  NLOW + BM  |
  * _YI     | B0 + BI    |       B0 + BI       |  NLOW + BM  | FIRST SEQUENCE PROTEIN
  * _YD     | B0 + BD    |       B0 + BD       |  NLOW + BM  |
  *---------+------------+---------------------+-------------+-----------------------------------
  * _MY     | NLOW + ME  |       E0 + ME       |   E0  + ME  |
  * _IY     | NLOW + IE  |       E0 + IE       |   E0  + IE  | LAST SEQUENCE PROTEIN
  * _DY     | NLOW + DE  |       E0 + DE       |   E0  + DE  |
  *---------+------------+---------------------+-------------+-----------------------------------
  */

////////////////////////////// FOR SSE2 to WORK //////////////////////////////////////////////
// WARNING: MATCH=2,INSERTION=3,DELETION=0,DUMMY=1
enum VectorPosition {
	/* Positions within both 4-tuple and array of 4-tuple */
	MATCH=2,
	INSERTION=3,
	DELETION=0,
	EXTRA=1,
	/* Position of empty space within 4-tuple */
	DUMMY=1,
	/* Positions of transition from ? to ? within array of 4-tuples */
	/* MATCH VECTOR */
	_MM = 4*MATCH+MATCH,     _MI = 4*MATCH+INSERTION,     _MD = 4*MATCH+DELETION,     _MX    = 4*MATCH+EXTRA,
	/* INSERTION VECTOR */
	_IM = 4*INSERTION+MATCH, _II = 4*INSERTION+INSERTION, _ID = 4*INSERTION+DELETION, _IX    = 4*INSERTION+EXTRA,
	/* DELETION VECTOR */
	_DM = 4*DELETION+MATCH,  _DI = 4*DELETION+INSERTION,  _DD = 4*DELETION+DELETION,  _DX    = 4*DELETION+EXTRA,
	/* EXTRA VECTOR */
	_XM = 4*EXTRA+MATCH,     _XI = 4*EXTRA+INSERTION,     _XD = 4*EXTRA+DELETION,     _DUMMY = 4*EXTRA+EXTRA,
	/* Overall Size of transtion structure */
	INSERTION_TRANSITIONS_SIZE = 16 /* 4 times sizeof(ScoreTuple) */
};
#define PREPROCESSOR_MATCH_PLACE 2

enum ProfileType { PF_MATRIX=0b1, PF_PATTERN=0b10};

#define TRANSITION_INDEX_FROM_TO(A,B) (4*(A)+(B))

////////////////////////////// FOR SSE to WORK //////////////////////////////////////////////
// WARNING: MATCH and INSERTION should be next to each other (in this ordering)
//          and either in the upper part of xmm or the lower part.
//          Do no mix them and correct the storing function according to your above choice.
extern inline void __ALWAYS_INLINE
StoreMatchInsertion( __m64 * const _address, const __m128 _reg)
{
	//_mm_storel_pi(_address, _reg);
	_mm_storeh_pi(_address, _reg);
}
// For vertical xali1 we have to define a store match and deletion, hence shuffle will be
// necessary to accomodate the above.
extern inline void __ALWAYS_INLINE
StoreMatchDeletion( __m64 * const _address, const __m128i _reg)
{
	const __m128i tmp = _mm_shuffle_epi32 (_reg, _MM_SHUFFLE(DUMMY,DUMMY,DELETION,MATCH));
	_mm_storel_epi64((__m128i*) _address, tmp);
}

// COVERAGE FUNCTION TO STORE SCORE AND LENGTH in 64 bits
union ScoreLength {
	struct {
		int Score;
		int Length;
	} Element;
	__m64 mm;
};
static inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
StoreScoreLength( union ScoreLength * const _address, __m128i _score, __m128i _length)
{
	_score = _mm_unpacklo_epi32(_score, _length);
	//_mm_storel_pi(_address, _reg);
	_mm_storeh_pi(&_address->mm, (__m128)_score);
}

/*
 ************************************************************************************************
 *                              FUNCTION POINTERS PART 1                                        *
 ************************************************************************************************
 */
typedef int (*NormalizedToRawFunctionPtr)(const float, const float * const restrict, const float, const size_t);
typedef float (*RawToNormalizedFunctionPtr)(const int, const float * const restrict, const float, const size_t);

/*
 ************************************************************************************************
 *                                    DEFINITIONS                                               *
 ************************************************************************************************
 */
#ifdef __USE_32BIT_INTEGER__
typedef int StoredIntegerFormat;
typedef union ScoreTuple {
	StoredIntegerFormat To[4];
	StoredIntegerFormat From[4];
	__m128i vector;
} ScoreTuple;
#else
typedef short int StoredIntegerFormat;
typedef union ScoreTuple {
	StoredIntegerFormat To[4];
	StoredIntegerFormat From[4];
	__m64 vector;
} ScoreTuple;
#endif

typedef union {
	ScoreTuple From[4];
	StoredIntegerFormat Element[INSERTION_TRANSITIONS_SIZE];
} TransitionScores;

typedef struct InsertionScores {
	StoredIntegerFormat * Alphabet;
	StoredIntegerFormat * Boundaries;
	TransitionScores * Transitions;
	ScoreTuple * FirstSequenceProtein;
	ScoreTuple * LastSequenceProtein;
} struct_InsertionScores;

/************************** Match Matrix   ************************/
typedef struct MatchScores {
	StoredIntegerFormat * alphabet;
} struct_MatchScores;

typedef struct Disjoint {
	char CDIS[KDIS][16];
	int  JDIP[KDIS];
	int  NDIP[KDPM];
	int  MDIS;
} SDisjoint;

typedef struct NormalizationItem {
	float RNOP[KNPM];		// Normalization coefficients
	char  CNTX[32];			// description as read in profile TEXT
	int   MNOR;					// index of the function type within CNOR (MNOR<KNOR)
	int   NNOR;					// Mode as read in profile MODE
	int   NNPR;					// Priority as read in profile PRIORITY
} SNormalizationItem;

typedef struct Normalization {
	SNormalizationItem Values[MAXN];
	char  CNOR[KNOR][16];		// description of the mode function
	int   JNOP[KNOR];       // number of coefficients
	int   JNOR;
} SNormalization;

enum NormalizationMode { LINEAR=0, GLE_ZSCORE=1, GLE_ZSCAVE=2 };
extern const char NormalizationModeName[3][16];

typedef struct Average {
	float * Weights;
	size_t size;
} SAverage;

typedef struct CutOffItem {
	char  CCUT[32];    // description as read in profile TEXT
	float RCUT[MAXN];  // Mode normalized cutoff as read in profile N_SCORES
	int   MCUT[MAXN];  // Modes as read in profile MODE
	int   MCLE;        // Level as read in profile LEVEL
	int   ICUT;        // Filter cutoff SCORE
	unsigned int HCUT; // Heuristic cutoff H_SCORE
	int   JCNM;	       // Number of modes defined, used dimension of mode arrays
} SCutOffItem;

typedef struct CutOff {
	SCutOffItem Values[MAXC];
	int JCUT;
} SCutOff;

struct Profile {
	char Identification[64] __attribute__((aligned(16)));
	char AC_Number[64] __attribute__((aligned(16)));
	char Date[128] __attribute__((aligned(16)));
	char Description[256] __attribute__((aligned(16)));
	unsigned char Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));
	unsigned char Complement_Alphabet_Mapping[ALPHABET_SIZE+2] __attribute__((aligned(16)));
	char CABC[ALPHABET_SIZE+2];
	enum ProfileType Type;
	char * Pattern;
	char * Sequence;
	size_t Length;
	size_t Alphabet_Length;

	union Scores {
		struct SInsertion {
			StoredIntegerFormat * Alphabet;
			StoredIntegerFormat * Boundaries;
			TransitionScores * Transitions;
			ScoreTuple * FirstSequenceProtein;
			ScoreTuple * LastSequenceProtein;
			size_t AlignStep;
			struct_MatchScores _dummy;
		} Insertion;
		struct SMatch {
			struct_InsertionScores _dummy;
			size_t AlignStep;
			StoredIntegerFormat * Alphabet;
		} Match;
	} Scores;

	SNormalization NormalizationData;
	SDisjoint DisjointData;
	SCutOff CutOffData;

	NormalizedToRawFunctionPtr NormalizedToRaw;
	RawToNormalizedFunctionPtr RawToNormalized;
	SAverage Average;

	short int LevelIndex;						/* WARNING: This is not the real level but the index corresponding within the cutoff array */
	short int ModeIndex;          	/* WARNING: This is not the real mode but the index corresponding within the normalization array */
	short int HeuristicModeIndex;		/* WARNING: This is not the real mode but the index corresponding within the normalization array */
	short int ModeIndexWithinLevel;	/* WARNING: This is not the real Mode but the index within the cutoff array of modes. */

	unsigned int HeuristicCutOff;
	int CutOff;
	float NormalizedCutOff;
	enum NormalizationMode NormalizationType;
	float * restrict NormalizationCoefs;

	_Bool isCircular;
	_Bool isReversed;

	struct Profile * next;
	struct Profile * previous;
};

#ifndef PFSEQ
#define PFSEQ
typedef struct PFSequence {
	unsigned char * ProfileIndex;
	unsigned char * OriginalSequence;
	size_t Length;
} PFSequence;
#endif

struct UniProtMatch {
	char (*UniqueIdentifier)[16];
	char (*EntryName)[16];
	char *State;
	char (*DB)[2];
	int *DBIndex;
	unsigned int * HeuristicScore;
	int * FilterScore;
	size_t Count;
	unsigned int True_positive;
	unsigned int Unknown;
	unsigned int Partial;
	unsigned int False_negative;
	unsigned int False_posistive;
};

/**************************** Strand ****************************/
enum Strand {FORWARD=0, REVERSE_COMPLEMENT=3, BOTH=4, REVERSE=1, COMPLEMENT=2};

/******************* Function choice enumerator *****************/
enum Version { SSE2=0, SSE41=1};

/************** Sequence generated from profile options *********/
enum GeneratedSequenceOptions {
	GENERATE_MATCH=1,
	GENERATE_INSERTION=2,
	GENERATE_DELETION=4
};

/********************* Heuristic structures *********************/
typedef union TransposeMatrix { const int * i; float * f;} TransposeMatrix;

/******************* XALIP & XALIT structures *******************/
union URegion {
	struct {
		int Begin;
		int End;
	};
	__m64 mm;
};

union Positions {
	struct {
		union URegion Protected;
		union URegion Sequence;
	} Region;
	__m128i xmm;
} __attribute__((aligned(16))) ;

union lScores {
	/* SSE 4.1 can work on integer */
	int Element[4];
	__m128i xmm;
	/* Other have to rely upon float */
	float Elementf[4];
	__m128 xmmf;
};

struct Alignment {
//     union Positions Region __attribute__((aligned(16)));
	union {
	/* WARNING: fasearch outputs use pointer to union URegion, hence the ordering is important
		*          keep Protected and then Sequence
		*/
		struct {
			union URegion Protected;
			union URegion Sequence;
		};
		union URegion array[2];
		__m128i xmm;
	} Region __attribute__((aligned(16)));
	/* The following order has to be kept due to sse zeroing, 128 bit alignment is required too*/
	//    int ProtectedRegionBegin;  // JAL1
	//    int ProtectedRegionEnd;    // JAL2;
	//    int SequenceBegin;         // JALB;
	//    int SequenceEnd;           // JALE;
	/* ---- */
	int Score;                    // JALS
	int IPMB;			 // Filled during xalit, used in coverage as length of found sequence
	int IPME;			 // Filled during xalit
};

struct DBSequence;
typedef struct PrintInput {
	const char * SequenceFile;
	const struct DBSequence * DB;
	unsigned int SeqId;
} PrintInput_t;

/*
 ************************************************************************************************
 *                           FUNCTION POINTERS PART 2                                           *
 ************************************************************************************************
 */
typedef unsigned int (*HeuristicFunctionPtr)(const TransposeMatrix TransposeMatch,
																						 const size_t Alphabet_Length,
																						 const size_t Profile_Length,
																						 const PFSequence * const restrict Sequence);
typedef void (*PrintFunctionPtr)(const struct Profile * const, const char * * const restrict,
																 const struct Alignment * const, char * const, const size_t,
																 const float, const int, const PrintInput_t * const);

/*
 ************************************************************************************************
 *                                   GLOBAL VARIABLES                                           *
 ************************************************************************************************
 */
extern unsigned int OutputPrintWidth;

/*
 ************************************************************************************************
 *                              PROFILE FUNCTION DECLARATIONS                                   *
 ************************************************************************************************
 */
PFIMPEXP int PrepareExtraTable(struct Profile * const prf);
PFIMPEXP struct Profile * ReversePprofilerofile(const struct Profile * const restrict inprf);
PFIMPEXP int ReadProfile(const char * const restrict FileName, struct Profile * const prf, const _Bool SetExtraTable);
PFIMPEXP int WriteProfile(const char * const restrict FileNameIn, struct Profile * const prf,
		          const char * const restrict AdditionalComments, FILE * const prfStreamOut);
PFIMPEXP int GenerateDummyProfile(struct Profile * const prf, const size_t Length);
PFIMPEXP void FreeProfile(struct Profile * const prf, const _Bool IsPointer);
PFIMPEXP int SeekProfileLevel(const struct Profile * const prf, const int Level);
PFIMPEXP int SeekProfileMode(const struct Profile * const prf, const int Mode);
PFIMPEXP int SetProfileMode(struct Profile * const prf, const int Mode);
PFIMPEXP void SeekProfileLevelAndMode(int * const restrict LevelIndex, int * const restrict ModeIndexWithinLevel,
				      int * const restrict ModeIndex, int * const restrict HeuristicModeIndex,
				      const struct Profile * const prf, const int Level, const int Mode);
PFIMPEXP int SetProfileLevelAndMode(struct Profile * const prf, const int Level, const int Mode);
PFIMPEXP int RescoreAlignment(const struct Profile * const prf, char * const AlignedSequence,
			      const struct Alignment * const alignment, const size_t SequenceLength);
PFIMPEXP int GenerateSequences(const struct Profile * const prf, const float MaxEponent, const char * const FileName,
		               const long int Seed, const enum GeneratedSequenceOptions Options, const unsigned int N);
PFIMPEXP int ReadProfileMatch(const char * const restrict FileName, struct UniProtMatch * Data);
PFIMPEXP void FreeProfileMatch(struct UniProtMatch * Data);

PFIMPEXP struct Profile * ReverseProfile(const struct Profile * const restrict inprf);

#ifndef __USE_INLINE_FUNCTIONS__
PFIMPEXP void InitializeDefault(union Scores * const matrices, char * const MatchSymbol, char * const InsertionSymbol);
PFIMPEXP void FreeScores(union Scores * const matrices);
PFIMPEXP int AllocateScores(union Scores * const matrices, const size_t Alphabet_Length, const size_t Profile_Length);
PFIMPEXP int ComputeHeuristicCutoff(struct Profile * const prf, const int Mode, const int CutOff);
#else
#include "pfProfileInline.h"
#endif

/*
 ************************************************************************************************
 *                              FILTER FUNCTION DECLARATIONS                                    *
 ************************************************************************************************
 */
PFIMPEXP int xali1_sse2 (const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                         int * const WORK, const size_t BSEQ, const size_t LSEQ, const int CutOff, const _Bool LOPT);
PFIMPEXP int xali1_sse41(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                         int * const WORK, const size_t BSEQ, const size_t LSEQ, const int CutOff, const _Bool LOPT);

/*
 ************************************************************************************************
 *                              ALIGNMENT FUNCTION DECLARATIONS                                 *
 ************************************************************************************************
 */
PFIMPEXP int xalip_sse2(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                        union lScores * const restrict iop, union Positions * const restrict iom,
                        union Positions * const restrict ioi, const int bseq, const int lseq,
                        struct Alignment * const restrict alignment,
                        _Bool * const restrict Lock, const size_t N1, const size_t N2, const _Bool Lopt,
                        const int Cutoff, const size_t MaxNumberOfAlignment);
PFIMPEXP int xalip_sse41(const struct Profile * const restrict prf, const unsigned char * const restrict Sequence,
                         union lScores * const restrict iop, union Positions * const restrict iom,
                         union Positions * const restrict ioi, const int bseq, const int lseq,
                         struct Alignment * const restrict alignment,
                         _Bool * const restrict Lock, const size_t N1, const size_t N2, const _Bool Lopt,
                         const int Cutoff, const size_t MaxNumberOfAlignment);
PFIMPEXP int xalit_sse2(const struct Profile * const restrict prf, const size_t N1, const size_t N2, const size_t bseq,
                        const PFSequence * const restrict PFSeq, char * const restrict CALI,
                        union lScores * const restrict iop,
                        struct Alignment * const restrict alignment, const _Bool * const restrict Lock);
PFIMPEXP int xalit_sse41(const struct Profile * const restrict prf, const size_t N1, const size_t N2, const size_t bseq,
                         const PFSequence * const restrict PFSeq, char * const restrict CALI,
                         union lScores * const restrict iop,
                         struct Alignment * const restrict alignment, const _Bool * const restrict Lock);

/*
 ************************************************************************************************
 *                             HEURISTIC FUNCTION DECLARATIONS                                  *
 ************************************************************************************************
 */
PFIMPEXP unsigned int heuristic(const struct Profile * const restrict prf, const PFSequence * const restrict Sequence,
                                const unsigned int CutOff);
PFIMPEXP unsigned int TransposeHeuristic_sse2(const TransposeMatrix TransposeMatch, const size_t Alphabet_Length,
                                              const size_t Profile_Length, const PFSequence * const restrict Sequence);
PFIMPEXP unsigned int TransposeHeuristic_sse41(const TransposeMatrix TransposeMatch, const size_t Alphabet_Length,
                                               const size_t Profile_Length, const PFSequence * const restrict Sequence);
PFIMPEXP unsigned int TransposeHeuristicGivenMemory_sse2(const float * const restrict TransposeMatch, float * const Memory,
                                                         const size_t Alphabet_Length, const size_t Profile_Length,
                                                         const PFSequence * const restrict Sequence);
PFIMPEXP unsigned int TransposeHeuristicGivenMemory_sse41(const int * const restrict TransposeMatch, int * const Memory,
                                                          const size_t Alphabet_Length, const size_t Profile_Length,
                                                          const PFSequence * const restrict Sequence);
#ifndef __USE_INLINE_FUNCTIONS__
PFIMPEXP const int * TransposeAndConvertMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
                                                    const size_t Profile_Length);
PFIMPEXP void TransposeAndConvertMatchMatrixGivenMemory(int * const restrict TIMatch, const union Scores * const Matrices,
                                                        const size_t Alphabet_Length, const size_t Profile_Length,
                                                        const size_t Aligned_Profile_Length);
PFIMPEXP float * TransposeAndConvertToFloatMatchMatrix(const union Scores * const Matrices, const size_t Alphabet_Length,
                                                       const size_t Profile_Length);
PFIMPEXP void TransposeAndConvertToFloatMatchMatrixGivenMemory(float * const restrict TIMatch, const union Scores * const Matrices,
                                                               const size_t Alphabet_Length, const size_t Profile_Length,
                                                               const size_t Aligned_Profile_Length);
#else
#  include "pfHeuristicInline.h"
#endif
/*
 ************************************************************************************************
 *                                NORMALIZATION FUNCTION DECLARATIONS                           *
 ************************************************************************************************
 */
PFIMPEXP int InitAverage(struct Profile * const restrict prf); //const union Scores * Matrices, const size_t prfLength, const size_t AlphabetLength, SAverage * const Average);
PFIMPEXP float ComputeAverageFrequencies(const PFSequence * const restrict Sequence, const SAverage * const restrict Average);

PFIMPEXP int N2R_1(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);
PFIMPEXP int N2R_2(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);
PFIMPEXP int N2R_3(const float R, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);
PFIMPEXP float R2N_1(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);
PFIMPEXP float R2N_2(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);
PFIMPEXP float R2N_3(const int N, const float * const restrict Coefs, const float AverageValue, const size_t SeqLength);


/*
 ************************************************************************************************
 *                             OUTPUT FUNCTION DECLARATIONS                                     *
 ************************************************************************************************
 */
PFIMPEXP void PrintSimple(const struct Profile * const prf, const char * * const AlignedSequence,
													const struct Alignment * const alignment, char * const Header,
													const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintDefault(const struct Profile * const prf, const char * * const AlignedSequence,
													 const struct Alignment * const alignment, char * const Header,
													 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintPfscanLOpt(const struct Profile * const prf, const char * * const AlignedSequence,
															const struct Alignment * const alignment, char * const Header,
															const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintInterpro(const struct Profile * const prf, const char * * const AlignedSequence,
														const struct Alignment * const alignment, char * const Header,
														const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintPfscan(const struct Profile * const prf, const char * * const AlignedSequence,
			  const struct Alignment * const alignment, char * const Header,
			  const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintTSV(const struct Profile * const prf, const char * * const AlignedSequence,
		       const struct Alignment * const alignment, char * const Header,
		       const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintClassification(const struct Profile * const prf, char * * const AlignedSequence,
																	const struct Alignment * const alignment, char * const Header,
																	const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintIncmatch(const struct Profile * const prf, const char * * const AlignedSequence,
			    const struct Alignment * const alignment, char * const Header,
			    const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintPSMaker(const struct Profile * const prf, const char * * const AlignedSequence,
			   const struct Alignment * const alignment, char * const Header,
			   const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintxPSA( const struct Profile * const prf, const char * * const AlignedSequence,
												 const struct Alignment * const alignment, char * const Header,
												 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintPsScan(const struct Profile * const prf, const char * * const AlignedSequence,
													const struct Alignment * const alignment, char * const Header,
													const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintOneLine( const struct Profile * const prf, const char * * const AlignedSequence,
														const struct Alignment * const alignment, char * const Header,
														const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintTurtle( const struct Profile * const prf, const char * * const AlignedSequence,
														const struct Alignment * const alignment, char * const Header,
														const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

PFIMPEXP void PrintSAM(const struct Profile * const prf, const char * * const AlignedSequence,
											 const struct Alignment * const alignment, char * const Header,
											 const size_t SequenceLength, const float RAVE, const int N, const PrintInput_t * const extra);

/*
 ************************************************************************************************
 *                                   INLINE FUNCTIONS                                           *
 ************************************************************************************************
 */
extern inline void __ALWAYS_INLINE NextInsertionProfile( struct SInsertion * Insertion)
{
  const size_t step       = Insertion->AlignStep;
  Insertion->Alphabet    += step;
  Insertion->Boundaries  += INSERTION_BOUNDARIES_SIZE;
  Insertion->Transitions++;
}

extern inline void __ALWAYS_INLINE PreviousInsertionProfile( struct SInsertion * Insertion)
{
  const size_t step       = Insertion->AlignStep;
  Insertion->Alphabet    -= step;
  Insertion->Boundaries  -= INSERTION_BOUNDARIES_SIZE;
  Insertion->Transitions--;
}

extern inline void __ALWAYS_INLINE CopyPreviousInsertionProfile(struct SInsertion * Insertion)
{
  const size_t step = Insertion->AlignStep;
  memcpy(Insertion->Alphabet, Insertion->Alphabet - step, step*sizeof(StoredIntegerFormat));
  memcpy(Insertion->Boundaries, Insertion->Boundaries - INSERTION_BOUNDARIES_SIZE, INSERTION_BOUNDARIES_SIZE*sizeof(StoredIntegerFormat));
  memcpy(Insertion->Transitions, Insertion->Transitions - 1, sizeof(TransitionScores));
}

extern inline void __ALWAYS_INLINE NextMatchProfile( struct SMatch * Match)
{
  const size_t step = Match->AlignStep;
  Match->Alphabet  += step;
}

extern inline void __ALWAYS_INLINE PreviousMatchProfile( struct SMatch * Match)
{
  const size_t step = Match->AlignStep;
  Match->Alphabet  -= step;
}

extern inline size_t __ALWAYS_INLINE GetInsertionMemory(const char * const key, struct SInsertion * const Insertion, StoredIntegerFormat ** pointer)
{
	/* return 0 if single value, 1 if vector and 2 if error */

	StoredIntegerFormat * Alphabet    = Insertion->Alphabet;
	StoredIntegerFormat * Boundaries  = Insertion->Boundaries;
	StoredIntegerFormat * Transitions = Insertion->Transitions->Element;

	switch(key[0]) {
		case 'I':
			switch(key[1]) {
				case '\0': *pointer = Alphabet    +   1; return 1; break;
				case '0' : *pointer = Alphabet    +   0; return 0; break;
				case 'M' : *pointer = Transitions + _IM; return 0; break;
				case 'I' : *pointer = Transitions + _II; return 0; break;
				case 'D' : *pointer = Transitions + _ID; return 0; break;
				case 'E' : *pointer = Boundaries  + _IE; return 0; break;
				default  : return 2;
			};
			break;
		case 'B':
			switch(key[1]) {
				case '0' : *pointer = Boundaries + _B0; return 0; break;
				case '1' : *pointer = Boundaries + _B1; return 0; break;
				case 'M' : *pointer = Boundaries + _BM; return 0; break;
				case 'I' : *pointer = Boundaries + _BI; return 0; break;
				case 'D' : *pointer = Boundaries + _BD; return 0; break;
				case 'E' : *pointer = Boundaries + _BE; return 0; break;
				default  : return 2;
			};
			break;
		case 'E':
			switch(key[1]) {
				case '0' : *pointer = Boundaries + _E0; return 0; break;
				case '1' : *pointer = Boundaries + _E1; return 0; break;
				default  : return 2;
			};
			break;
		case 'M':
			switch(key[1]) {
				case 'M' : *pointer = Transitions + _MM; return 0; break;
				case 'I' : *pointer = Transitions + _MI; return 0; break;
				case 'D' : *pointer = Transitions + _MD; return 0; break;
				case 'E' : *pointer = Boundaries  + _ME; return 0; break;
				default  : return 2;
			};
			break;
		case 'D':
			switch(key[1]) {
				case 'M' : *pointer = Transitions + _DM; return 0; break;
				case 'I' : *pointer = Transitions + _DI; return 0; break;
				case 'D' : *pointer = Transitions + _DD; return 0; break;
				case 'E' : *pointer = Boundaries  + _DE; return 0; break;
				default  : return 2;
			};
			break;
		default:
			return 2;
	};
}

static inline void __ALWAYS_INLINE FreeAverage(SAverage * const Average)
{
	if (Average->Weights) _mm_free(Average->Weights);
}

extern inline HeuristicFunctionPtr __ALWAYS_INLINE GetHeuristicVersion(const enum Version version)
{
	if (version == SSE2)
		return &TransposeHeuristic_sse2;
	else
		return &TransposeHeuristic_sse41;
}

#endif
