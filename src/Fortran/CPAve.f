*----------------------------------------------------------------------*     
* $Id: CPAve.f,v 2.4 2003/04/10 11:58:37 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine CPAve(IMPP,IDMP,LPRF,CABC,NABC,PAVE)

      Integer           IMPP(0:27,0:IDMP)
      Character         CABC(0:26)
      Real              PAVE(0:26)

      Do  I1=0,NABC
         PAVE(I1)=0
      End Do
      Do  I1=1,LPRF
         Do  I2=0,NABC
            PAVE(I2)=PAVE(I2)+IMPP(I2,I1) 
         End do  
      End Do 

      Return
      End
