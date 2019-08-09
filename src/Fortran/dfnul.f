*----------------------------------------------------------------------*     
* $Id: dfnul.f,v 2.4 2003/03/26 14:40:29 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine DFNUL(FABC,P0,DABC,MABC)

      Real              FABC(0:26)
      Character*(*)     DABC

      If(MABC.LT.20) then
         MABC=4
         DABC='ACGT'
         Do I1=1,4
            FABC(I1)=0.25
         End do
      Else
         MABC=20
         DABC='ACDEFGHIKLMNPQRSTVWY'
         FABC( 1)=.08713
         FABC( 2)=.03347
         FABC( 3)=.04687
         FABC( 4)=.04953
         FABC( 5)=.03977
         FABC( 6)=.08861
         FABC( 7)=.03362
         FABC( 8)=.03689
         FABC( 9)=.08048
         FABC(10)=.08536
         FABC(11)=.01475
         FABC(12)=.04043
         FABC(13)=.05068
         FABC(14)=.03826
         FABC(15)=.04090
         FABC(16)=.06958
         FABC(17)=.05854
         FABC(18)=.06472
         FABC(19)=.01049
         FABC(20)=.02992
      End if 

      If(P0.LE.0.0) P0=1.0
      Do I1=1,MABC
         FABC(I1)=FABC(I1)*P0
      End do

      Return
      End
