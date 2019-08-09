*----------------------------------------------------------------------*     
* $Id: remsf.f,v 2.7 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine REMSF
     *   (NERR,NMSF,FMSF,
     *   IDMS,CSEQ,NSEQ,LSEQ,
     *   IDMF,RWGT,SQID,
     *   IRC)

      Character*(*)    FMSF

      Character        CSEQ(*) 
      Real*4           RWGT(*)
      Character*(*)    SQID(*) 

      Logical          LBKS
      Character*512    RCIN 
      Character*32     CMIS

      IRC=0
      LBKS=.FALSE.

      If(FMSF.NE.'-') then 
         Open(NMSF,File=FMSF,Status='OLD',Err=900)
      End if 

      CMIS='..'
 1    Read(NMSF,'(A)',Err=901,End=902) RCIN        
      If(Index(RCIN,'..').EQ.0) go to   1

* read length

      CMIS='MSF'
      IX=Index(RCIN,'MSF:')
      If(IX.EQ.0) Go to 904
      CMIS='length'
      Read(RCIN(IX+4:),*,Err=903) LSEQ

* names

      NSEQ=0
 2    Read(NMSF,'(A)',Err=901,End=905) RCIN        
      If(RCIN(1:2).EQ.'//') go to   3
      IX1=Index(RCIN,'Name: ')
      If(IX1.EQ.0) go to   2 
      IX2=Index(RCIN,'Len: ')
      If(IX2.EQ.0) go to   2 
      NSEQ=NSEQ+1
      If(NSEQ.GT.IDMF) go to 910
      SQID(NSEQ)=RCIN(IX1+5:IX2-1)
      IX=Index(RCIN,'Weight: ') 
      CMIS='Weight'
      If(IX.EQ.0) go to 904
      Read(RCIN(IX+9:),*,Err=903)  RWGT(NSEQ)
      Go to   2

* sequence offset -> NOFS

 3    Read(NMSF,'(A)',Err=901,End=905) RCIN        
      If(RCIN.EQ.' ') go to   3
      N1=1
      N2=0

* read sequences 

 10   Read(RCIN,*,Iostat=IOS) K1,K2
      If(IOS.NE.0.OR.K1.LE.0.OR.K2.LE.0) then 
         L1=Lblnk(RCIN) 
         Do I1=1,L1
            If(RCIN(I1:I1).NE.' ') go to  11
         End do
 11      J1=I1
         Do I1=J1,L1
            If(RCIN(I1:I1).EQ.' ') go to  12
         End do
 12      J1=I1
         Do I1=J1+1,L1
            If(RCIN(I1:I1).NE.' ') go to  13
         End do 
 13      NOFS=I1-1
         J1=I1
         Do I1=J1,L1
            If(RCIN(I1:I1).NE.' ') N2=N2+1 
         End do 
C          Backspace(NMSF,Err=907)
         LBKS=.TRUE.
      Else if(NOFS.EQ.0) then
         NOFS=Index(RCIN,' 1')
         N2=K2
         If(NOFS.EQ.0) go to   3
      Else
         N1=K1
         N2=K2
      End if 

      Do  20 I1=0,NSEQ-1

         If(LBKS) then
            LBKS=.FALSE.
            Go to 16
         End if           

 15      Read(NMSF,'(A)',Err=901,End=905) RCIN        
 16      If(RCIN(1:NOFS).EQ.' ') go to  15

         J1=I1*LSEQ+N1
         N5=(I1+1)*LSEQ
         Do  19 I2=NOFS+1,Lblnk(RCIN) 
            If(RCIN(I2:I2).NE.' ') then
               If(J1.GT.IDMS) go to 906
               If(J1.GT.N5) go to 911
               CSEQ(J1)=RCIN(I2:I2)
               J1=J1+1
            End if 
 19      Continue
         If(I1.NE.0.AND.N4.NE.N2-N1+1) go to 909
         
         If(J1.NE.I1*LSEQ+N2+1) then
            N3=J1-I1*LSEQ-1
            If(N3.NE.LSEQ) then
               Go to 909
            End if
         End if 
         N4=N2-N1+1

 20   Continue

      If(N2.GE.LSEQ) go to 100 

 25   Read(NMSF,'(A)',Err=901,End=905) RCIN        
      If(RCIN.EQ.' ') go to  25 
      N1=N2+1
      Go to  10

 100  Return

* errors

 900  Write(NERR,*) 'Error: Unable to open MSF file'//
     *   ' ''',FMSF(1:Lblnk(FMSF)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to read MSF file'//
     *   ' ''',FMSF(1:Lblnk(FMSF)),'''.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Unexpected end of MSF file. '//
     *   'Unable to find ''',CMIS(1:Lblnk(CMIS)),''' keyword.'
      Write(NERR,*) '       No sequences read.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Badly formatted MSF file. '//
     *   'Unable to read ''',CMIS(1:Lblnk(CMIS)),''' parameter.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      Write(NERR,*) '       No sequences read.'
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Badly formatted MSF file. '//
     *   'Unable to find ''',CMIS(1:Lblnk(CMIS)),''' keyword.'
      Write(NERR,*) '       No sequences read.'
      IRC=1
      Go to 100
 905  Write(NERR,*) 'Error: Unexpected end of MSF file. '//
     *   'Unable to read sequences in MSF file.'
      IRC=1
      Go to 100
 906  Write(NERR,*) 'Error: Multiple sequence length exceeds ',
     *   'buffer size (',IDMS,').'
      IRC=1
      Go to 100
 909  Write(NERR,*) 'Error: Sequence length differs from expected ',
     *   'length.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
 910  Write(NERR,*) 'Error: Number of sequences in MSF file exceeds ',
     *   'buffer size (',IDMF,').'
      IRC=1
      Go to 100
 911  Write(NERR,*) 'Error: Sequence length exceeds expected ',
     *   'length.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100

      End
