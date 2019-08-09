*----------------------------------------------------------------------*     
* $Id: remsa.f,v 1.2 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine REMSA
     *   (NERR,NMSF,FMSF,
     *   IDMS,CSEQ,NSEQ,LSEQ,
     *   IDMF,RWGT,SQID,
     *   IRC)

      Character*(*)    FMSF

      Character        CSEQ(*) 
      Real*4           RWGT(*)
      Character*(*)    SQID(*) 

      Character*512    RCIN 

      IRC=0
      NSEQ=0
      LSEQ=0

      If(FMSF.NE.'-') then 
         Open(NMSF,File=FMSF,Status='OLD',Err=900)
      End if 

 1    Read(NMSF,'(A)',Err=901,End=999,Iostat=IOS) RCIN
 2    If(RCIN(1:1).NE.'>') go to   1   

* read sequence weight in xpsa header if available

      If(NSEQ+1.GT.IDMF) go to 906
      IW=Index(RCIN,'weight=')
      If(IW.NE.0) then
         Read(RCIN(IW+7:),*,Err=905) RWGT(NSEQ+1)
      Else
         RWGT(NSEQ+1)=1.0
      End if

* read sequence identifier

      L=Lblnk(RCIN)
      Do I1=2,L
         If(RCIN(I1:I1).NE.' ') go to 3
      End do
 3    Do I2=I1,L
         If(RCIN(I2:I2).EQ.' ') go to 4
      End do
 4    If(I2-I1.GT.64) I2=I1+64
      SQID(NSEQ+1)=RCIN(I1:I2)

* read sequences in multiple alignment

      J1=NSEQ*LSEQ
 10   RCIN=' '
      Read(NMSF,'(A)',Err=901,Iostat=IOS) RCIN
      L=Lblnk(RCIN)
      If((RCIN(1:1).EQ.'>').OR.
     *   (L.EQ.0.AND.IOS.EQ.-1)) go to 20
      Do 15 I1=1,L
         If(RCIN(I1:I1).NE.' ') then
            J1=J1+1
            If(J1.GT.IDMS) go to 902
            CSEQ(J1)=RCIN(I1:I1)
         End if
 15   Continue
      Go to 10

 20   If(J1.EQ.(NSEQ*LSEQ)) go to 903
      If(NSEQ.EQ.0) then
         LSEQ=J1
      Else if(J1.NE.(NSEQ+1)*LSEQ) then
         go to 904
      End if
      NSEQ=NSEQ+1
      
      If(IOS.NE.-1) go to 2

 100  Return

 110  IX=Lblnk(SQID(NSEQ+1))
      If(IX.GT.1) then
         Write(NERR,*) '       While processing sequence ',
     *      SQID(NSEQ+1)(1:IX)
      End if
      Go to 100

* errors

 900  Write(NERR,*) 'Error: Unable to open MSA file'//
     *   ' ''',FMSF(1:Lblnk(FMSF)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to read MSA file'//
     *   ' ''',FMSF(1:Lblnk(FMSF)),'''.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: MSA length exceeds buffer ',
     *   'size (',IDMS,').'
      IRC=1
      Go to 110
 903  Write(NERR,*) 'Error: Unexpected end of sequence.'//
     *   ' Sequence has zero length.'
      IRC=1
      Go to 110
 904  Write(NERR,*) 'Error: Sequence length differs from expected ',
     *   'length.'
      IRC=1
      Go to 110
 905  Write(NERR,*) 'Error: Unable to read weight in sequence ',
     *   'header.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
 906  Write(NERR,*) 'Error: Number of sequences in MSA file exceeds ',
     *   'buffer size (',IDMF,').'
      IRC=1
      Go to 100

 999  If(NSEQ.EQ.0) then
         Write(NERR,*) 'Error: No sequence has been read.'
         IRC=1
      End if
      Go to 100

      End
