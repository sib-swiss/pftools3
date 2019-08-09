*----------------------------------------------------------------------*     
* $Id: wrsam.f,v 2.4 2003/03/26 14:40:29 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine WRSAM
     *   (NOUT,IDMP,RIHM,RMHM,LHMM,NABC,FABC,FLOW,RSCA,BLOG) 

      Include         'pfind.f' 
      Include         'hmdat.f'

      Real             FABC(0:26)
      Real*8           BLOG
      Real             RPEX(64)

* Log to Prob

      Do  10 I1=0,LHMM
         Do I2=1,46
            If     (RIHM(I2,I1).LE.FLOW) then
               RIHM(I2,I1)=0.0
            Else If(RIHM(I2,I1).LE.-6.0) then
               RIHM(I2,I1)=0.000001
            Else
               RIHM(I2,I1)=EXP(RIHM(I2,I1))
            End if 
         End do
         If(I1.GT.0) then  
            Do I2=1,27
               If     (RMHM(I2,I1).LE.FLOW) then
                  RMHM(I2,I1)=0.0
               Else If(RIHM(I2,I1).LE.-6.0) then
                  RMHM(I2,I1)=0.000001
               Else
                  RMHM(I2,I1)=EXP(RMHM(I2,I1))
               End if 
            End do
         End if 
 10   Continue

* header 

      RSCA=RSCA/BLOG
      Write(NOUT,'(
     *   ''% HMM converted from generalized profile'',
     *   '' with pftools program htop.'',
     *   /,''% Logarithic base of profile: '',F12.8,''.'',
     *   /,''% Offset: profile score = '',E10.4,'' + Log(Prob).'',
     *   /,''MODEL:'')')
     *   exp(BLOG),RSCA

      If     (NABC.EQ.20) then 
         Write(6,'(''alphabet protein'')')    
         NPAR=49
      Else if(NABC.EQ. 4) then
         Write(6,'(''alphabet DNA'')')    
         NPAR=17
      Else
         Go to 100
      End if 

* Write position 0

      RPEX( 1)=0.0  
      RPEX( 2)=0.0  
      RPEX( 3)=0.0  

      RPEX( 4)=0.0  
      RPEX( 5)=0.0  
      RPEX( 6)=0.0  

      RPEX( 7)=0.0  
      RPEX( 8)=RIHM(MI, 0)
      RPEX( 9)=RIHM(II, 0)

      J2=9
      Do I2=1,NABC
         J2=J2+1
         RPEX(J2)=0
      End do 
      Do I2=1,NABC
         J2=J2+1
         RPEX(J2)=RIHM(I2, 0)
      End do 

      Write(NOUT,'(''Begin'',3(/,F8.6,2F9.6))')
     *   (RPEX(ii1),ii1=1,9)
      If(NABC.EQ.20) then
         Write(NOUT,'(F8.6,4F9.6,7(/,F8.6,4F9.6))')
     *      (RPEX(ii1),ii1=10,49)
      Else
         Write(NOUT,'(F8.6,3F9.6,/,F8.6,3F9.6)')
     *      (RPEX(ii1),ii1=10,17)
      End if

* Write position 1 to LHMM 

      Do  I1=1,LHMM

         RPEX( 1)=RIHM(DD,I1-1)
         RPEX( 2)=RIHM(MD,I1-1)
         RPEX( 3)=RIHM(ID,I1-1)

         RPEX( 4)=RIHM(DM,I1-1)
         RPEX( 5)=RIHM(MM,I1-1)
         RPEX( 6)=RIHM(IM,I1-1)

         RPEX( 7)=RIHM(DI,I1  )
         RPEX( 8)=RIHM(MI,I1  )
         RPEX( 9)=RIHM(II,I1  )

         J2=9
         Do I2=1,NABC
            J2=J2+1
            RPEX(J2)=RMHM(I2,I1)
         End do 
         Do I2=1,NABC
            J2=J2+1
            RPEX(J2)=RIHM(I2,I1)
         End do 

         If(I1.LT.LHMM) then
            Write(NOUT,'(/,I4,3(/,F8.6,2F9.6))')
     *         I1,(RPEX(ii1),ii1=1,9)
         Else 
            Write(NOUT,'(/,''  -1'',3(/,F8.6,2F9.6))')
     *         (RPEX(ii1),ii1=1,9)
         End if
         If(NABC.EQ.20) then
            Write(NOUT,'(F8.6,4F9.6,7(/,F8.6,4F9.6))')
     *         (RPEX(ii1),ii1=10,49)
         Else
            Write(NOUT,'(F8.6,3F9.6,/,F8.6,3F9.6)')
     *         (RPEX(ii1),ii1=10,17)
         End if

      End do 

* End 

      RPEX( 1)=0.0
      RPEX( 2)=0.0
      RPEX( 3)=0.0

      RPEX( 4)=RIHM(DM,LHMM)
      RPEX( 5)=RIHM(MM,LHMM)
      RPEX( 6)=RIHM(IM,LHMM)

      RPEX( 7)=0.0
      RPEX( 8)=0.0
      RPEX( 9)=0.0

      Do I2=10,9+2*NABC
         RPEX(I2)=0.0
      End do

      Write(NOUT,'(/,''End'',3(/,F8.6,2F9.6))')
     *   (RPEX(ii1),ii1=1,9)
      If(NABC.EQ.20) then
         Write(NOUT,'(F8.6,4F9.6,7(/,F8.6,4F9.6))')
     *      (RPEX(ii1),ii1=10,49)
      Else
         Write(NOUT,'(F8.6,3F9.6,/,F8.6,3F9.6)')
     *      (RPEX(ii1),ii1=10,17)
      End if

      Write(6,'(''ENDMODEL'')')    

 100  Return
      End 
