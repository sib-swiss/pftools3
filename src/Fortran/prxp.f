*----------------------------------------------------------------------*     
* $Id: prxp.f,v 1.2 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*        Print psa alignment for multiple matches
*        Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine PRXP
     *   (CALI,IABG,IAED,IPBG,IPED,LPRF,NW)
      
      Character*01      CALI(*)
      Integer           IABG
      Integer           IAED
      Integer           LPRF
      Integer           IPBG
      Integer           IPED

      Character*512     RCOUT

      IIDE=LPRF-IPED

      INL=IAED-IABG+IPBG+IIDE
      INB=INL/NW
      INR=INL-INB*NW
C      IN1=IABG-1
      ii2=IABG

      IPOS=1

* write psa style alignment

      Do IN2=1,INB
         Do ii1=1,NW
            If(IPOS.LT.IPBG.OR.IPOS.GT.IPED) then
               RCOUT(ii1:ii1)='-'
               IPOS=IPOS+1
            Else
               RCOUT(ii1:ii1)=CALI(ii2)
               ii2=ii2+1
               K1=Ichar(RCOUT(ii1:ii1))
               If(K1.GE.65.AND.K1.LE. 90) IPOS=IPOS+1
            End if
         End do
C         IN1=IN1+NW
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,NW)
      End do
      If(INR.GT.0) then
         Do ii1=1,INR
            If(IPOS.LT.IPBG.OR.IPOS.GT.IPED) then
               RCOUT(ii1:ii1)='-'
               IPOS=IPOS+1
            Else
               RCOUT(ii1:ii1)=CALI(ii2)
               ii2=ii2+1
               K1=Ichar(RCOUT(ii1:ii1))
               If(K1.GE.65.AND.K1.LE. 90) IPOS=IPOS+1
            End if
         End do
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,INR)
      End if            
      
      Return

      End
