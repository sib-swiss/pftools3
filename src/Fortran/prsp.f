*----------------------------------------------------------------------*     
* $Id: prsp.f,v 1.4 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*        Print sequence or psa alignment
*        Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine PRSP
     *   (CABC,ISEQ,CALI,IABG,IAED,NW,OPTS,OPTX)
      
      Character         CABC(0:26)
      Integer*2         ISEQ(*)
      Character*01      CALI(*)
      Integer           IABG
      Integer           IAED
      Logical           OPTS
      Logical           OPTX

      Character*512     RCOUT


      INL=IAED-IABG+1
      INB=INL/NW
      INR=INL-INB*NW
      IN1=IABG-1

* write sequence of aligned region

      If(OPTS) then
         Do IN2=1,INB
            Do ii1=1,NW
               RCOUT(ii1:ii1)=CABC(ISEQ(IN1+ii1))
            End do
            IN1=IN1+NW
            Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,NW)
         End do
         If(INR.GT.0) then
            Do ii1=1,INR
               RCOUT(ii1:ii1)=CABC(ISEQ(IN1+ii1))
            End do
            Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,INR)
         End if            

* write psa style alignment
* Do NOT use to display individual matches of a circular profile
* use prxp() instead.

      Else if(OPTX) then
         Do IN2=1,INB
            Do ii1=1,NW
               RCOUT(ii1:ii1)=CALI(IN1+ii1)
            End do
            IN1=IN1+NW
            Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,NW)
         End do
         If(INR.GT.0) then
            Do ii1=1,INR
               RCOUT(ii1:ii1)=CALI(IN1+ii1)
            End do
            Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,INR)
         End if            
      End if

      Return

      End
