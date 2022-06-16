/*******************************************************
                        PFTOOLS
 *******************************************************
  May 30, 2013 QuantileRegression.c
 *******************************************************
 (C) 2013 SIB Swiss Institute of Bioinformatics
     Thierry Schuepbach (thierry.schuepbach@sib.swiss)
 *******************************************************/
#include "config.h"

/*      SUBROUTINE RQ0(M,N,M5,N2,A,B,T,TOLER,IFT,X,E,S,WA,WB)
C
C     Modified to remove SOL and related vars -- only good for single tau
C     M Number of Observations
C     N Number of Parameters
C     M5 = M+5  row dimension for WA
C     N2 = N+2  col dimension for WA
C     A is the X matrix
C     B is the y vector
C     T, the desired quantile
C     TOLER, smallest detectable |x-y|/x machine precision to the 2/3
C     IFT exit code:
C		0-ok
C		else dimensions inconsistent or T not in (0,1)
C     X the parameter estimate betahat
C     E is the residual vector
C     S is an integer work array (M)
C     WA is a real work array (M5,N2)
C     WB is another real work array (M)
C     Utilization:  If you just want a solution at a single quantile you
C     The algorithm  is a slightly modified version of Algorithm AS 229
C     described in Koenker and D'Orey, "Computing Regression Quantiles,
C     Applied Statistics, pp. 383-393.
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER I,J,K,KL,KOUNT,KR,M,M1,M2,M3,M4,M5,IFT
      INTEGER N,N1,N2,OUT,S(M)
      LOGICAL STAGE,TEST,INIT,IEND
      DOUBLE PRECISION MIN,MAX
      DOUBLE PRECISION B(M),A(M,N),X(N),WA(M5,N2),WB(M),E(M) */

extern void rq0_(const int * M,
		 const int * N,
		 const int * M5,
		 const int * N2,
		 double * A,
		 double * B,
		 const double * T,
		 const double * TOLER,
		 int * IFT,
		 double * X,
		 double * E,
		 int * S,
		 double * WA,
		 double * WB
		);

int QuantileRegression(const unsigned int * const restrict HeuristicScores, const int * const restrict FilterScore,
		       float coefs[], const double quantile, const int size)
{
  int result;
  const int N = 2;
  const int N2 = N + 2;
  const double TOLER = 1.E-8;

  /* computes M5 in order to have cache aligned columns */
  const int M5 = (size + 5 + (64/sizeof(double)-1)) & ~(64/sizeof(double)-1);

  /* allocates memory (N2(WA) + 1(WB) + 1(S) + 1(E) + 1(S)*/
  double * const mem         = (double*) _mm_malloc((N2+4)*M5*sizeof(double), 64);
  const size_t ldaCacheSize  = (size_t) ((size + 0x7) & ~0x7);
  double * const data        = (double*) _mm_malloc(3*ldaCacheSize*sizeof(double), 64);
  double * const restrict A0 = data;
  double * const restrict A1 = &data[ldaCacheSize];
  double * const restrict B  = &data[2*ldaCacheSize];
  if (mem != NULL & data != NULL) {
    /* Set input data */
    register int maxFilter = 0;
    for (int i=0; i<size; ++i) {
	A0[i] = 1.0;
	maxFilter = maxFilter < FilterScore[i] ? FilterScore[i] : maxFilter;
	A1[i] = (double) FilterScore[i];
	B[i]  = (double) HeuristicScores[i];
    }

    /* Prepare a guess */
    double guess[2] = { 0.0, 1.0/(float) maxFilter};
    /* Compute quantile regression */
    rq0_(&size, &N, &M5, &N2, A0, B, &quantile, &TOLER, &result,
	 guess /*X*/, &mem[M5] /*E*/, (int*) &mem[2*M5] /*S*/,
	 &mem[4*M5] /*WA*/, &mem[3*M5] /*WB*/
	);
    _mm_free(mem);
    _mm_free(data);
    if (result == 0) {
//       printf(" y = %lf + %lf * x\n", guess[0], guess[1]);
      coefs[0] = guess[0];
      coefs[1] = guess[1];
      return 0;
    }
    else {
      return -1;
    }
  }
  else
    return -1;
}