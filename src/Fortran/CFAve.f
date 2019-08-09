*----------------------------------------------------------------------*     
* $Id: CFAve.f,v 2.5 2003/04/10 11:58:37 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine CFAve(ISEQ,IDMS,BSEQ,LSEQ,CABC,NABC,FAVE)

      Integer*2         ISEQ(IDMS)
      Integer           BSEQ
      Character         CABC(0:26)
      Real              FAVE(0:26)

      Do  I1=0,NABC
         FAVE(I1)=0
      End Do
      Do  I1=BSEQ,LSEQ
         J1=ISEQ(I1)
         FAVE(J1)=FAVE(J1)+1
      End Do 
      Do  I1=0,NABC
         FAVE(I1)=FAVE(I1)/(LSEQ-BSEQ)
      End Do

      Return
      End
