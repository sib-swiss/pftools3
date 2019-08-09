*----------------------------------------------------------------------*     
* $Id: prali.f,v 2.8 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine PRALI
     *   (LPRF,CHIP,CHMP,IDMP,LSEQ,LREV,
     *   CALI,LALI,JALB,JALE,NW)

      Character*01      CHIP(0:IDMP)
      Character*01      CHMP(IDMP)
      Character*01      CALI(*)
      Logical           LREV

      Character*512     RCPR
      Character*512     RCSQ
      Character*534     RCOUT
      

* Find beginning of alignment, sequence, profile                

      JB=0
      KB=1 
      Do I1=1,LALI
         If(CALI(I1).EQ.'-') then 
            JB=JB+1
            If(JB.GT.LPRF) JB=JB-LPRF
            KB=KB+1
         Else 
            Go to   2
         End if 
      End do 
 2    LB=JALB 
  
C      Write(NERR,*) 'JB: ',JB,' KB: ',KB,' LB: ',LB

* Find end of alignment, sequence, profile                

      JE=LPRF
      KE=LALI 
      Do I1=LALI,KB+1,-1
         If(CALI(I1).EQ.'-') then 
            JE=JE-1
            KE=KE-1
         Else 
            Go to   3
         End if 
      End do 
 3    Continue
      
C      Write(NERR,*) 'JE: ',JE,' KE: ',KE,' LALI: ',LALI

* Write alignment  

      IX=0

      Do I1=KB,KE

         If(IX.GE.NW) then 
            NPE=JB
            If(.NOT.LREV) then 
               NSE=LB
C               If(RCSQ(1:1).NE.'-') NSE=NSE-1
               NSE=NSE-1
            Else
               NSE=LSEQ-LB+1 
C               If(RCSQ(1:1).NE.'-') NSE=NSE+1
               NSE=NSE+1
            End if
            NPE=NPE-LPRF-1
            NSE=NSE-LSEQ-1
            Write(6,'(''#'')')
            Write(RCOUT(1:),'(''# P'',I9,'' '')') NPB
            RCOUT(14:)=RCPR(1:NW)
            Write(RCOUT(14+NW:),'(I9)') NPE
            Write(6,'(534A)')(RCOUT(ii1:ii1),ii1=1,NW+22)
            Write(RCOUT(1:),'(''# S'',I9,'' '')') NSB
            RCOUT(14:)=RCSQ(1:NW)
            Write(RCOUT(14+NW:),'(I9)') NSE
            Write(6,'(534A)')(RCOUT(ii1:ii1),ii1=1,NW+22)
            IX=0
         End if 

         IX=IX+1
         K1=Ichar(CALI(I1))
         If     (K1.GE.65.AND.K1.LE. 90) then 
            JB=JB+1
            If(JB.GT.LPRF) JB=JB-LPRF
            LB=LB+1
            RCPR(IX:IX)=CHMP(JB)
            RCSQ(IX:IX)=CALI(I1)
         Else if(K1.GE.97.AND.K1.LE.122) then 
            LB=LB+1
            RCPR(IX:IX)=CHIP(JB)
            RCSQ(IX:IX)=CALI(I1)
         Else if(CALI(I1).EQ.'-') then 
            JB=JB+1
            If(JB.GT.LPRF) JB=JB-LPRF
            RCPR(IX:IX)=CHMP(JB)
            RCSQ(IX:IX)=CALI(I1)
         End if 
         
         If(IX.EQ.1) then
            NPB=JB
            If(K1.GE.97.AND.K1.LE.122) NPB=NPB+1 
            If(.NOT.LREV) then 
               NSB=LB
               If(RCSQ(1:1).NE.'-') NSB=NSB-1
            Else
               NSB=LSEQ-LB+1 
               If(RCSQ(1:1).NE.'-') NSB=NSB+1
            End if
         End if
         
      End do  

      If(IX.GT.0) then 
         NPE=JB
         If(.NOT.LREV) then 
            NSE=LB
C            If(RCSQ(1:1).NE.'-') NSE=NSE-1
            NSE=NSE-1
         Else
            NSE=LSEQ-LB+1 
C            If(RCSQ(1:1).NE.'-') NSE=NSE+1
            NSE=NSE+1
         End if
         NPE=NPE-LPRF-1
         NSE=NSE-LSEQ-1
         Write(6,'(''#'')')    
         RCPR(IX+1:)=' '
         RCSQ(IX+1:)=' '
         Write(RCPR(IX+1:IX+10),'(I9)') NPE
         Write(RCSQ(IX+1:IX+10),'(I9)') NSE
         IX=IX+10
         Write(6,'(''# P'',I9,'' '',512A)')
     *      NPB,(RCPR(ii1:ii1),ii1=1,IX)
         Write(6,'(''# S'',I9,'' '',512A)')
     *      NSB,(RCSQ(ii1:ii1),ii1=1,IX)
      End if
      Write(6,'(''#'')')    

      Return 
      End
