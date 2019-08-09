*       Program pfw
*----------------------------------------------------------------------*     
* $Id: pfw.f,v 2.10 2004/01/09 09:27:36 vflegel Exp $
*----------------------------------------------------------------------*     
*       Function: Calculate weights for individual sequences of a 
*                 multiple sequence alignment.      
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* DATA
*----------------------------------------------------------------------*

* array sizes, I/O units

C       Parameter        (IDM1=16777216)
C        Parameter        (IDM1=1048576)
C        Parameter        (IDM2=   9999)
C        Parameter        (IDM3=   2024)

      Include          'ardim.f'

      Parameter        (NMSF=     11)
      Parameter        (ICSL=IDMS/26)

      Include          'sterr.f'

* weights and distances 

      Real*4            RWGT(IDMF)
      Integer           NDIS(IDMF) 
      Integer           NIND(IDMF) 

* multiple sequence alignment:

      Character*512     FMSF

      Character*64      SQID(IDMF)
      Character         CSEQ(IDMS)
      Character         CSQR(IDMP)

      Logical           LXSC(IDMP)

* character set 

      Character         CSET(26,ICSL)
      Integer           NSET(   ICSL)

* work fields 

      Character*512     RCIO
      Character*512     ROUT
      Character*32      RCTM
      Logical           OPTM
      Character         B

* functions

      Integer           Getc
C       Integer           Fputc

* character translation 

      Character*27      ABCU
      Character*27      ABCL

      Data              ABCU/'ABCDEFGHIJKLMNOPQRSTUVWXYZ-'/
      Data              ABCL/'abcdefghijklmnopqrstuvwxyz.'/

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0
      IDUM=-2

*read command line 

      Call Repar
     *   (FMSF,NRAN,RX,RW,IRAN,OPTM,IRC) 
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''pfw 2.3 revision 5.d'',//
     *      ''Usage: pfw [ -hmNXRW ] [ msf-file | - ] '',
     *      ''[ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -h: print usage help text.'',/
     *      ''    -m: input sequences in MSA format.'',/
     *      ''    -N<value>:'',/
     *      ''        number of shuffles per sequence (default:'',
     *      '' 100).'',/
     *      ''    -X<value>:'',/
     *      ''        gap excision threshold (default: 0.5).'',/
     *      ''    -R<value>:'',/
     *      ''        random number seed, negative integer '',
     *      ''(default: -123456789).'',/
     *      ''    -W<value>:'',/
     *      ''        total weight (default: 1).'',/
     *      )')
         Write(NERR,'(
     *      '' valid (but deprecated) parameters are:'',/,
     *      ''  [N=shuffles-per-seq]      use option -N instead'',/
     *      ''  [X=gap-excision]          use option -X instead'',/
     *      ''  [R=random-seed]           use option -R instead'',/
     *      ''  [W=total-weight]          use option -W instead'',/
     *      )')
         Call Exit(IRC)
      End if

* read msf-file

C       If(FMSF.EQ.'-') then 
C   1      Open(NMSF,Status='SCRATCH') 
C          Do I1=1,2*IDM1
C             If(Getc(B).NE.0) go to   2
C             N1=Fputc(NMSF,B) 
C          End do 
C   2      Rewind(NMSF)
C       End if 

      If(FMSF.EQ.'-') then 
 1       Open(NMSF,Status='SCRATCH',Err=901)
 2       Continue 
         Do I1=1,512
            If(Getc(B).NE.0) go to 3
            If(Ichar(B).EQ.10) then
               Write(NMSF,'(512A)',Err=903)(CSEQ(ii1),ii1=1,I1-1)
               Go to   2   
            Else 
               CSEQ(I1)=B
            End if
         End do  
         Go to   2
 3       Rewind(NMSF)

      End if 

      If(OPTM) then
         Call REMSA
     *      (NERR,NMSF,FMSF,
     *      IDMS,CSEQ,NSEQ,LSEQ,
     *      IDMF,RWGT,SQID,
     *      IRC)
      Else
         Call REMSF
     *      (NERR,NMSF,FMSF,
     *      IDMS,CSEQ,NSEQ,LSEQ,
     *      IDMF,RWGT,SQID,
     *      IRC)
      End if
      If(IRC.NE.0) go to 100

* Check needed buffer size

      If((NSEQ*LSEQ).GT.ICSL) go to 904

*----------------------------------------------------------------------*
* DATA PROCESSING
*----------------------------------------------------------------------*

* make sequence uppercase 

      Do I1=1,LSEQ*NSEQ
         If     (Index(ABCU,CSEQ(I1)).NE.0) then 
            Continue
         Else
            IX=Index(ABCL,CSEQ(I1))
            If(IX.GT.0) then
               CSEQ(I1)=ABCU(IX:IX)
            Else 
               CSEQ(I1)='-'
            End if
         End if
      End do

* gap excision

      NXSC=INT(RX*NSEQ)

      Do I1=1,LSEQ
         J1=I1
         K1=0
         Do I2=1,NSEQ
            If(CSEQ(J1).EQ.'-') K1=K1+1
            J1=J1+LSEQ
         End do
         If(K1.LE.NXSC) then 
            LXSC(I1)=.TRUE.
         Else
            LXSC(I1)=.FALSE.
         End if
      End do

* generate character set profile 

      Call GCSET
     *   (IDMS,CSEQ,NSEQ,LSEQ,
     *   ICSL,NSET,CSET) 

* major loop 

* - initialize new weights

*       Do  10 I1=1,NSEQ
      RWGT(I1)=0
 10   Continue 

      Do  20 I1=1,NRAN*NSEQ

* - generate random sequences 

         Call RanSQ(IRAN,ICSL,NSET,CSET,LSEQ,CSQR)

* - compare random sequence to real sequences
*
*      NDIS(i): distance of random sequence to real sequence i  
*      MIND   : minimal distance 
*      KMIN   : # of real sequences with minimal distance
*      NIND(i): indices of real sequences with minimal distance   

         J3=1 
         MIND=LSEQ
         KMIN=0
         Do  15 I2=1,NSEQ
            NDIS(I2)=0
            Do  12 I3=1,LSEQ
               If(LXSC(I3)) then
                  If(CSQR(I3).NE.CSEQ(J3)) NDIS(I2)=NDIS(I2)+1
               End if
               J3=J3+1
 12         Continue
            If     (NDIS(I2).LT.MIND) then
               MIND=NDIS(I2)
               KMIN=1
               NIND(1)=I2
            Else if(NDIS(I2).EQ.MIND) then
               KMIN=KMIN+1
               NIND(KMIN)=I2
            End if 
 15      Continue

         R1=1.0/KMIN
         Do  16 I2=1,KMIN
            J2=NIND(I2)
            RWGT(J2)=RWGT(J2)+R1
 16      Continue

 20   Continue

      Do  30 I1=1,NSEQ
         RWGT(I1)=RWGT(I1)/(NSEQ*NRAN)
C          Write(6,'(A16,''Weight: '',F6.4)')
C    *        SQID(I1)(1:16),RWGT(I1)  
 30   Continue 

*----------------------------------------------------------------------*
* OUTPUT SECTION
*----------------------------------------------------------------------*

      Rewind(NMSF)

* output msa format with xpsa header

      If(OPTM) then
         K1=0
 40      Read(NMSF,'(A)',End=900,Err=902) RCIO
         If(RCIO(1:1).NE.'>') go to 40
         
* build 'weight=value' pair

 41      K1=K1+1
         RCTM=' weight= '
         Write(RCTM(10:),*) RW*RWGT(K1)

* remove extraoneous blanks

         J1=8
         Do I1=9,Lblnk(RCTM)
            If(RCTM(I1:I1).NE.' ') then 
               J1=J1+1
               RCTM(J1:J1)=RCTM(I1:I1)
               RCTM(I1:I1)=' '
            End if
         End do
         ILTM=Lblnk(RCTM)

* find and replace preexisting 'weight=value' in xpsa header

         L=Lblnk(RCIO)
         IX=Index(RCIO(1:L),' weight=')

         If(IX.EQ.0) then
            If((L+ILTM+1).GE.512) L=L-ILTM-1
            RCIO(L+1:L+1+ILTM)=RCTM
            Write(6,'(512A)')(RCIO(ii1:ii1),ii1=1,L+1+ILTM)
         Else
            Do I1=IX+7,L
               If(RCIO(I1:I1).NE.' ') go to 45
            End do
 45         Do I2=I1,L
               If(RCIO(I2:I2).EQ.' ') go to 46
            End do
 46         ROUT(1:IX)=RCIO(1:IX)
            ROUT(IX:IX+ILTM)=RCTM
            ROUT(IX+ILTM:)=RCIO(I2:L)
            L=Lblnk(ROUT)
            Write(6,'(512A)')(ROUT(ii1:ii1),ii1=1,L)
         End if

 48      RCIO=' '
         Read(NMSF,'(A)',Err=902,Iostat=IOS) RCIO
         If(RCIO(1:1).EQ.'>') go to 41
         L=Lblnk(RCIO)
         If(L.EQ.0.AND.IOS.EQ.-1) go to 100
         Write(6,'(512A)')(RCIO(ii1:ii1),ii1=1,L)
         Go to  48 

* output MSF format
         
      Else
 51      Read(NMSF,'(A)',End=900,Err=902) RCIO
         L=Lblnk(RCIO)
         Write(6,'(512A)')(RCIO(ii1:ii1),ii1=1,L)
         If(Index(RCIO(1:L),'..').EQ.0) go to  51
         
         K1=1
 52      Read(NMSF,'(A)',End=900,Err=902) RCIO
         L=Lblnk(RCIO)
         IX=Index(RCIO(1:L),'Weight: ') 
         If(IX.NE.0) then 
            If(RW*RWGT(K1).LT.10) then
               Write(RCIO(IX+8:),'(F6.4)') RW*RWGT(K1)
               L=IX+13 
            Else
               Write(RCIO(IX+8:),*) RW*RWGT(K1)
               L=Lblnk(RCIO) 
            End if  
            K1=K1+1
         End if   
         Write(6,'(512A)')(RCIO(ii1:ii1),ii1=1,L)
         If(RCIO(1:2).NE.'//') go to  52
         
 53      Read(NMSF,'(A)',End=100,Err=902) RCIO
         L=Lblnk(RCIO)
         Write(6,'(512A)')(RCIO(ii1:ii1),ii1=1,L)
         Go to  53 
      End if


 100  Call Exit(IRC)

* errors

 900  Go to 100 

 901  Write(NERR,*) 'Error: Unable to create temporary file.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Unable to read MSF file.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Unable to write to temporary file.'
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Sequence length times sequence number ',
     *   'exceeds buffer size (',ICSL,').'
      IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FMSF,NRAN,RX,RW,IRAN,OPTM,IRC) 

      Character*512     FMSF 
      Character*512     CPAR
      Logical           OPTM

      NRAN=100
      RX=0.5
      IRAN=-123456789
      RW=1.0
      IRC=0
      FMSF=' '
      OPTM=.FALSE.
      
      N1=Iargc()
      If(N1.LT.1) go to 900  
      
      K1=0
      I2=1
      Do  50 I1=1,N1
         Call GetArg(I2,CPAR)
         If     (CPAR(1:1).EQ.'-'.AND.CPAR(2:2).NE.' ') then
            If(Index(CPAR,'h').NE.0) then 
               IRC=1
               Go to 100
            End if
            If(Index(CPAR,'m').NE.0) OPTM=.TRUE.
            If(Index(CPAR,'N').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) NRAN
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) NRAN
               End if
            End if
            If(Index(CPAR,'X').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RX
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RX
               End if
            End if
            If(Index(CPAR,'R').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) IRAN
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) IRAN
               End if
            End if
            If(Index(CPAR,'W').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RW
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RW
               End if
            End if
 
         Else if(CPAR(1:2).EQ.'N=') then
            Read(CPAR(3:),*,Err=900) NRAN 
         Else if(CPAR(1:2).EQ.'X=') then
            Read(CPAR(3:),*,Err=900) RX 
         Else if(CPAR(1:2).EQ.'R=') then
            Read(CPAR(3:),*,Err=900) IRAN  
         Else if(CPAR(1:2).EQ.'W=') then
            Read(CPAR(3:),*,Err=900) RW 

         Else if(FMSF.EQ.' ') then
            FMSF=CPAR
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 60         

 50   Continue

 60   If(FMSF.EQ.' ') go to 900
      If(NRAN.LE.0) go to 900
      If(IRAN.GT.0) IRAN=-1-IRAN
      If(IRAN.EQ.0) IRAN=-123456789

 100  Return 
 900  IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine GCSET
     *   (IDMS,CSEQ,NSEQ,LSEQ,
     *   ICSL,NSET,CSET) 

* CSEQ(IDMS)
* NSET(ICSL)

      Character         CSEQ(*)
      Integer           NSET(*)
      Character         CSET(26,ICSL)
      
      Do  10 I1=1,LSEQ
         NSET(   I1)=0
         J2=I1-LSEQ
         Do   8 I2=1,NSEQ
            J2=J2+LSEQ
            Do   5 I3=1,NSET(I1)
               If(CSEQ(J2).EQ.CSET(I3,I1)) go to   8 
 5          Continue
            NSET(I1)=NSET(I1)+1
            CSET(NSET(I1),I1)=CSEQ(J2)
 8       Continue
 10   Continue

C       Do  20 I1=1,LSEQ
C          Write(6,'(I6,'' '',26A)')
C    *        NSET(I1),(CSET(ii1,I1),ii1=1,NSET(I1))
C  20   Continue

      Return
      End
*----------------------------------------------------------------------*
      Subroutine RanSQ(IRAN,ICSL,NSET,CSET,LSEQ,CSQR)

      Integer         NSET(*)
      Character       CSET(26,ICSL) 
      Character       CSQR(*)
      
      Do  10 I1=1,LSEQ
         If(NSET(I1).EQ.0) then
            CSQR(I1)='-'
         Else
            K1=RAN2(IRAN)*NSET(I1)+1
            CSQR(I1)=CSET(K1,I1)
         End if
 10   Continue

      Return
      End
*----------------------------------------------------------------------*
      FUNCTION RAN2(IDUM)
      PARAMETER (M=714025,IA=1366,IC=150889,RM=1.4005112E-6)
      DIMENSION IR(97)
      DATA IFF /0/
      IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
         IFF=1
         IDUM=MOD(IC-IDUM,M)
         DO 11 J=1,97
            IDUM=MOD(IA*IDUM+IC,M)
            IR(J)=IDUM
 11      CONTINUE
         IDUM=MOD(IA*IDUM+IC,M)
         IY=IDUM
      ENDIF
      J=1+(97*IY)/M
      IF(J.GT.97.OR.J.LT.1) PAUSE
      IY=IR(J)
      RAN2=IY*RM
      IDUM=MOD(IA*IDUM+IC,M)
      IR(J)=IDUM
      RETURN
      END
*----------------------------------------------------------------------*
      Include          'remsf.f'
      Include          'remsa.f'
      Include          'lblnk.f'
