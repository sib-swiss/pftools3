*----------------------------------------------------------------------*     
* $Id: RtoN.f,v 2.4 2003/04/10 11:58:37 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RtoN(N,R,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)

      Real              RNOP(KNPM,MAXN)

      If     (IFUN.EQ.1) then 
         R=RNOP(1,INOR)+RNOP(2,INOR)*Real(N)
      Else if(IFUN.EQ.2) then 
         R=( Real(N) /
     *      ( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
     *      - RNOP(4,INOR) ) / RNOP(5,INOR)
      Else if(IFUN.EQ.3) then 
         R=( (N-RAVE) /
     *      ( RNOP(1,INOR)*(1.0-EXP(RNOP(2,INOR)*LSEQ-RNOP(3,INOR) )))
     *      - RNOP(4,INOR) ) / RNOP(5,INOR)
      End if 
      Return
      End
