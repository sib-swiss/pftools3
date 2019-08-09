*       Program 6ft
*----------------------------------------------------------------------*
* $Id: 6ft.f,v 2.8 2003/11/28 11:53:33 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: 6-frame translation of DNA sequence into protein
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*

      Include          'ardim.f' 
      Parameter        (NOUT=6)

      Character*512     RCIN
      Character*512     CHDE 
      Character*64      CHID
      Character*512     RCEX

      Integer*2         ISEQ(IDMS)
      Character         CSEQ(IDMS)
      
      Character         CABC(0:63) 
      Character         CNBC(0:15) 

      Character*04      CDNA
      Integer           NBSQ
      Integer           IOS
      
      Include          'sterr.f'
      
* options
      
      Integer           NW
      Logical           OPTS
      Logical           OPTR

* Initializations

      Data              CABC( 0)/'K'/
      Data              CABC( 1)/'N'/
      Data              CABC( 2)/'K'/
      Data              CABC( 3)/'N'/
      Data              CABC( 4)/'T'/
      Data              CABC( 5)/'T'/
      Data              CABC( 6)/'T'/
      Data              CABC( 7)/'T'/
      Data              CABC( 8)/'R'/
      Data              CABC( 9)/'S'/
      Data              CABC(10)/'R'/
      Data              CABC(11)/'S'/
      Data              CABC(12)/'I'/
      Data              CABC(13)/'I'/
      Data              CABC(14)/'M'/
      Data              CABC(15)/'I'/
      Data              CABC(16)/'Q'/
      Data              CABC(17)/'H'/
      Data              CABC(18)/'Q'/
      Data              CABC(19)/'H'/
      Data              CABC(20)/'P'/
      Data              CABC(21)/'P'/
      Data              CABC(22)/'P'/
      Data              CABC(23)/'P'/
      Data              CABC(24)/'R'/
      Data              CABC(25)/'R'/
      Data              CABC(26)/'R'/
      Data              CABC(27)/'R'/
      Data              CABC(28)/'L'/
      Data              CABC(29)/'L'/
      Data              CABC(30)/'L'/
      Data              CABC(31)/'L'/
      Data              CABC(32)/'E'/
      Data              CABC(33)/'D'/
      Data              CABC(34)/'E'/
      Data              CABC(35)/'D'/
      Data              CABC(36)/'A'/
      Data              CABC(37)/'A'/
      Data              CABC(38)/'A'/
      Data              CABC(39)/'A'/
      Data              CABC(40)/'G'/
      Data              CABC(41)/'G'/
      Data              CABC(42)/'G'/
      Data              CABC(43)/'G'/
      Data              CABC(44)/'V'/
      Data              CABC(45)/'V'/
      Data              CABC(46)/'V'/
      Data              CABC(47)/'V'/
      Data              CABC(48)/'O'/
      Data              CABC(49)/'Y'/
      Data              CABC(50)/'O'/
      Data              CABC(51)/'Y'/
      Data              CABC(52)/'S'/
      Data              CABC(53)/'S'/
      Data              CABC(54)/'S'/
      Data              CABC(55)/'S'/
      Data              CABC(56)/'O'/
      Data              CABC(57)/'C'/
      Data              CABC(58)/'W'/
      Data              CABC(59)/'C'/
      Data              CABC(60)/'L'/
      Data              CABC(61)/'F'/
      Data              CABC(62)/'L'/
      Data              CABC(63)/'F'/

      Data              CNBC( 0)/'X'/ 
      Data              CNBC( 1)/'T'/ 
      Data              CNBC( 2)/'X'/ 
      Data              CNBC( 3)/'X'/ 
      Data              CNBC( 4)/'X'/ 
      Data              CNBC( 5)/'P'/ 
      Data              CNBC( 6)/'R'/ 
      Data              CNBC( 7)/'L'/ 
      Data              CNBC( 8)/'X'/ 
      Data              CNBC( 9)/'A'/ 
      Data              CNBC(10)/'G'/ 
      Data              CNBC(11)/'V'/ 
      Data              CNBC(12)/'X'/ 
      Data              CNBC(13)/'S'/ 
      Data              CNBC(14)/'X'/ 
      Data              CNBC(15)/'X'/ 

      IRC=0
      NBSQ=0
      NW=60
      IOS=0

      CDNA='ACGT'

      Call Repar(OPTS,OPTR,NW,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''6ft 2.3 revision 5.d'',//
     *      ''Usage: 6ft [-[r|s]hW] < seq-library-file '',/
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -r: translate only reverse (antisense) strand.'',/
     *      ''    -s: translate only sense strand.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -W<value>:'',/
     *      ''        specifies the output width (default: 60).'',/
     *      )')
         Call Exit(IRC)
      End if

 1    Read(5,'(A)',End=920,Err=900) RCIN
      If(RCIN(1:1).NE.'>') go to   1

 2    If(RCIN(1:1).NE.'>') Go to 920
      J1=Index(RCIN,' ')
      CHID='>x|' // RCIN(2:J1)
      LNID=Lblnk(CHID)
      LNDE=512-LNID+3
      CHDE=RCIN(J1+1:J1+LNDE)
      RCEX=CHID(1:LNID) // '   ' // CHDE(1:LNDE)  
      
      LSEQ=0
 3    Read(5,'(A)',Iostat=IOS,Err=900) RCIN
      If(RCIN(1:1).EQ.'>') go to  10
      L1=Lblnk(RCIN) 
      
      Do   5 I1=1,L1
         K1=Ichar(RCIN(I1:I1))
         If(K1.GE.97) then 
            K1=K1-32
            RCIN(I1:I1)=char(K1)
         End if
         If(K1.GT.90.OR.K1.LT.65) go to   5
         LSEQ=LSEQ+1
         If(LSEQ.GT.IDMS) Go to 901
         ISEQ(LSEQ)=Index(CDNA,(RCIN(I1:I1)))-1
         If(ISEQ(LSEQ).EQ.-1) ISEQ(LSEQ)=-64
 5    Continue
      RCIN=' '
      If(IOS.EQ.0) Go to 3

 10   Continue 

      If(LSEQ.GE.3) then 
         NBSQ=NBSQ+1
      Else
         Go to 902
      End if

* convert sequence into codons. 

      J1=LNID+1
      J2=LNID+2
      J3=Lblnk(RCEX)

* - plus strand
      
      If(.NOT.OPTR) then

         Do  12 I1=1,LSEQ-2
            N1=ISEQ(I1)*16+ISEQ(I1+1)*4+ISEQ(I1+2)
            If(N1.LT.0) then 
               N1=ISEQ(I1)*4+ISEQ(I1+1)
               If(N1.GE.0) then
                  CSEQ(I1)=CNBC(N1)
               Else 
                  CSEQ(I1)='X'
               End if
            Else
               CSEQ(I1)=CABC(N1)
            End if 
 12      Continue

         RCEX(J1:J2)='_1'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,1,LSEQ-2,NW)
         RCEX(J1:J2)='_2'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,2,LSEQ-2,NW)
         RCEX(J1:J2)='_3'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,3,LSEQ-2,NW)

      End if

* - minus strand
      
      If(.NOT.OPTS) then

         Do  13 I1=LSEQ,1,-1
            ISEQ(I1)=3-ISEQ(I1)
            If(ISEQ(I1).GT.4) ISEQ(I1)=-64
 13      Continue
         K1=1
         Do  14 I1=LSEQ,3,-1
            N1=ISEQ(I1)*16+ISEQ(I1-1)*4+ISEQ(I1-2)
            If(N1.LT.0) then 
               N1=ISEQ(I1)*4+ISEQ(I1-1)
               If(N1.GE.0) then 
                  CSEQ(K1)=CNBC(N1)
               Else 
                  CSEQ(K1)='X'
               End if
            Else
               CSEQ(K1)=CABC(N1)
            End if 
            K1=K1+1
 14      Continue

         RCEX(J1:J2)='_4'
         If(OPTR) RCEX(J1:J2)='_1'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,1,LSEQ-2,NW)
         RCEX(J1:J2)='_5'
         If(OPTR) RCEX(J1:J2)='_2'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,2,LSEQ-2,NW)
         RCEX(J1:J2)='_6'
         If(OPTR) RCEX(J1:J2)='_3'
         Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,J3)
         Call Prsq(CSEQ,3,LSEQ-2,NW)

      End if

      If(IOS.EQ.0) Go to   2

 100  Call Exit(IRC)

 900  Write(NERR,*) 'Error: Unable to read sequence from standard'//
     *   ' input.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: sequence length exceeds buffer size (',
     *   IDMS,').'
      Write(NERR,*) '       While processing sequence ',
     *      CHID(4:Lblnk(CHID))
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: sequence length too short.'
      Write(NERR,*) '       While processing sequence ',
     *      CHID(4:Lblnk(CHID))
      IRC=1
      Go to 100
 920  If(NBSQ.EQ.0) then
         Write(NERR,*) 'Warning: No FASTA sequence found. '
         IRC=-1
      End if
      Go to 100
      End

*----------------------------------------------------------------------*
      Subroutine Repar(OPTS,OPTR,NW,IRC)

      Character*64      CARG

      Integer           NW
      Logical           OPTS
      Logical           OPTR
      
      IRC=0
      
      OPTS=.FALSE.
      OPTR=.FALSE.
      
      N1=Iargc()
      
      K1=0
      I2=1
      Do I1=1,N1
         Call GetArg(I2,CARG)
         If(CARG(1:1).EQ.'-') then
            If     (Index(CARG,'h').NE.0) then
               go to 900
            Else if(Index(CARG,'s').NE.0) then
               OPTS=.TRUE.
            Else if(Index(CARG,'r').NE.0) then
               OPTR=.TRUE.
            Else if(Index(CARG,'W').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NW
               Else
                  I2=I2+1
                  If(I2.GT.N1) Go to 900
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NW
               End if
            End if
         Else
            K1=K1+1
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
      End do
      
 20   If(K1.GT.0.OR.(OPTS.AND.OPTR)) go to 900
      If(NW.LE.0.OR.NW.GT.512) NW=60
      
 100  Return
 900  IRC=-1
      Go to 100
      End

*----------------------------------------------------------------------*
      Subroutine Prsq(CSEQ,BSEQ,LSEQ,NW)

      Character*512     RCOUT
      Character         CSEQ(*)
      Integer           BSEQ

      If((LSEQ-BSEQ).LT.0) go to 900
      INL=(LSEQ-BSEQ)/3+1
      INB=INL/NW
      INR=INL-INB*NW

      ii2=BSEQ
      Do IN2=1,INB
         Do ii1=1,NW
            RCOUT(ii1:ii1)=CSEQ(ii2)
            ii2=ii2+3
         End do
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,NW)
      End do
      If(INR.GT.0) then
         Do ii1=1,INR
            RCOUT(ii1:ii1)=CSEQ(ii2)
            ii2=ii2+3
         End do
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,INR)
      End if

 100  Return
 
 900  Write(6,'(1A)')' '
      Go to 100

      End
*----------------------------------------------------------------------*
      Include          'lblnk.f'
