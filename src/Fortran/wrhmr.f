*----------------------------------------------------------------------*     
* $Id: wrhmr.f,v 2.4 2003/03/26 14:40:29 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine WRHMR
     *   (NOUT,NERR,IDMP,RIHM,RMHM,LHMM,NABC,CABC,FLOW,RSCA,BLOG) 

      Include         'pfind.f' 
      Include         'hmdat.f'

      Real*8           BLOG
      Character        CABC(0:26)

      Character*64     RCEX        
      

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

      Write(NOUT,'(''# HMM v1.7'')')
      RSCA=RSCA/BLOG
      Write(NERR,'(
     *   ''# HMM converted from generalized profile'',
     *   '' with pftools program htop.'',
     *   /,''# Logarithic base of profile: '',F12.8,''.'',
     *   /,''# Offset: profile score = '',E10.4,'' + Log(Prob).'')')
     *   exp(BLOG),RSCA

      Write(RCEX(1:8),'(I8)') LHMM 
      Do I1=1,8
         If(RCEX(I1:I1).NE.' ') Go to  12 
      End do
 12   RCEX(1:8)=RCEX(I1:8) // RCEX(1:I1-1)
      Write(NOUT,'(A8,''# M -- length of model'')') RCEX(1:8) 
      If(NABC.EQ.4) then  
         Write(NOUT,'(''4       # alphabet length'')') 
         Write(NOUT,'(''1       # alphabet type'')') 
         Write(NOUT,'( 4A,''    # alphabet'')')(CABC(ii1),ii1=1,NABC) 
      Else
         Write(NOUT,'(''20      # alphabet length'')') 
         Write(NOUT,'(''3       # alphabet type'')') 
         Write(NOUT,'(20A,''    # alphabet'')')(CABC(ii1),ii1=1,NABC) 
      End if
      Write(NOUT,
     *   '(''no      # Optional reference line annotation?'')')
      Write(NOUT,
     *   '(''no      # Optional consensus structure annotation?'')')

* make sure that meaningless parameter at pos. 0 equal 0 

      Do I1=1,27
         RMHM(I1,0)=0.0
      End do

* write HMM probabilites

      Do  I1=0,LHMM

* match state

         Write(RCEX(1:8),'(I8)') I1
         Do I2=1,8
            If(RCEX(I2:I2).NE.' ') Go to  21 
         End do
 21      RCEX='###MATCH_STATE ' // RCEX(I2:8)
         RCEX=RCEX(1:Lblnk(RCEX)) // ' ( ) ( )'
         Write(NOUT,'(64A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

         Write(NOUT,'(F8.6,''        # t_m1'')') RIHM(MM,I1)
         Write(NOUT,'(F8.6,''        # t_d1'')') RIHM(MD,I1)
         Write(NOUT,'(F8.6,''        # t_i0'')') RIHM(MI,I1)

         Do I2=1,NABC
            Write(NOUT,
     *         '(F8.6,''        # Symbol '',A,'' probability'')')
     *         RMHM(I2,I1),CABC(I2)
         End do 

* delete state

         Write(RCEX(1:8),'(I8)') I1
         Do I2=1,8
            If(RCEX(I2:I2).NE.' ') Go to  22 
         End do
 22      RCEX='###DELETE_STATE ' // RCEX(I2:8)
         RCEX=RCEX(1:Lblnk(RCEX)) // ' ( ) ( )'
         Write(NOUT,'(64A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

         Write(NOUT,'(F8.6,''        # t_m1'')') RIHM(DM,I1)
         Write(NOUT,'(F8.6,''        # t_d1'')') RIHM(DD,I1)
         Write(NOUT,'(F8.6,''        # t_i0'')') RIHM(DI,I1)

* insert state 

         Write(RCEX(1:8),'(I8)') I1
         Do I2=1,8
            If(RCEX(I2:I2).NE.' ') Go to  23 
         End do
 23      RCEX='###INSERT_STATE ' // RCEX(I2:8)
         RCEX=RCEX(1:Lblnk(RCEX)) // ' ( ) ( )'
         Write(NOUT,'(64A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

         Write(NOUT,'(F8.6,''        # t_m1'')') RIHM(IM,I1)
         Write(NOUT,'(F8.6,''        # t_d1'')') RIHM(ID,I1)
         Write(NOUT,'(F8.6,''        # t_i0'')') RIHM(II,I1)

         Do I2=1,NABC
            Write(NOUT,
     *         '(F8.6,''        # Symbol '',A,'' probability'')')
     *         RIHM(I2,I1),CABC(I2)
         End do 

      End do 

 100  Return
      End 
