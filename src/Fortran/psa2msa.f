*       Program psa2msa
*----------------------------------------------------------------------*
* $Id: psa2msa.f,v 2.11 2003/11/28 11:53:33 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: Converts a pfsearch -x output file into Pearson/Fasta
*                 multiple sequence alignment format
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*
* DATA
*----------------------------------------------------------------------*

      Include          'ardim.f'

      Parameter        (NSEQ=  11)

      Character*4096    FSEQ
      Character*01      CSEQ(IDMS)
      Integer           IPRF(0:IDMP)

      Character*512     RCIO
      Character         B

      Logical           OPTL
      Logical           OPTU
      Logical           OPTP
      Logical           OPTD
      Logical           LEOF

      Integer           Getc

      Include          'sterr.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      RC=0
      LPRF=0
      RCIO=' '
      IRC=0
      LEOF=.FALSE.
      NBSQ=0

* command line arguments

      Call Repar
     *   (FSEQ,OPTL,OPTU,OPTP,OPTD,NM,NW,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''psa2msa 2.3 revision 5.d'',//
     *      ''Usage: psa2msa [ -dhlpuMW ] [ psa-file | - ] '',
     *      ''[ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -d: replace periods by dashes on output.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -l: replace upper case letters by lower case.'',/
     *      ''    -p: replace dashes by periods on output.'',/
     *      ''    -u: replace lower case letters by upper case.''
     *      )')
         Write(NERR,'(
     *      ''    -M<value>:'',/
     *      ''        maximal insertion length (default: -1).'',/
     *      ''    -W<value>:'',/
     *      ''        specifies the output width (default: 60).'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',/,
     *      ''    [M=max-insert-length]  use option -M instead'',/
     *      ''    [W=output-width]       use option -W instead'',/
     *      )')
         Call Exit(IRC)
      End if

* open input file

      If(FSEQ.EQ.'-') then
 1       Open(NSEQ,Status='SCRATCH',Err=901)
 2       Continue
         Do I1=1,512
            If(Getc(B).NE.0) go to 3
            If(Ichar(B).EQ.10) then
               Write(NSEQ,'(512A)')(CSEQ(ii1),ii1=1,I1-1)
               Go to   2
            Else
               CSEQ(I1)=B
            End if
         End do
         Go to   2
 3       Rewind(NSEQ)
      Else
         Open(NSEQ,File=FSEQ,Status='OLD',Err=900)
      End if

* first sequence:

      Call GetSeq(NSEQ,CSEQ,LSEQ,RCIO,LDES,NBSQ,IDMS,IRC)
      If(IRC.GT.0.OR.NBSQ.EQ.0) go to 100
      If(IRC.LT.0) LEOF=.TRUE.

      Call UpdIPRF(CSEQ,LSEQ,LPRF,IPRF,IDMP,IRC)
      If(IRC.NE.0) Go to 100

* next sequence

 10   If(.NOT.LEOF) then
         Call GetSeq(NSEQ,CSEQ,LSEQ,RCIO,LDES,NBSQ,IDMS,IRC)
         If     (IRC.EQ.-1) then
            LEOF=.TRUE.
         Else if(IRC.NE. 0) then
            Go to 100
         End if

         Call UpdIPRF(CSEQ,LSEQ,LPRF,IPRF,IDMP,IRC)
         If(IRC.NE.0) Go to 100

         Go to  10
      End if

 50   Continue

      Rewind(NSEQ)
      NBSQ=0
      LEOF=.FALSE.

      LMSA=LPRF
      Do I1=0,LPRF
         LMSA=LMSA+IPRF(I1)
      End do
      LMSB=LMSA

 60   Call GetSeq(NSEQ,CSEQ,LSEQ,RCIO,LDES,NBSQ,IDMS,IRC)
      If     (IRC.EQ.-1) then
         LEOF=.TRUE.
         IRC=0
         Go to  70
      Else if(IRC.NE. 0) then
         Go to 100
      End if

 70   Call UpdSeq(CSEQ,LSEQ,LMSA,IPRF,LPRF,IDMP)
      Call EdtSeq(CSEQ,LMSA,OPTL,OPTU,OPTP,OPTD)
      If(NM.GE.0) Call CutSeq(CSEQ,LMSA,LMSB,IPRF,LPRF,NM,IDMP)

      Write(6,'((512A))')(RCIO(ii1:ii1),ii1=1,LDES)
      Call PRSQ(CSEQ,LMSB,NW)

      If(.NOT.LEOF) Go to  60

 100  Call Exit(IRC)

* errors

 900  Write(NERR,*) 'Error: Unable to open PSA file'//
     *   ' ''',FSEQ(1:Lblnk(FSEQ)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to create temporary file.'
      IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FSEQ,OPTL,OPTU,OPTP,OPTD,NM,NW,IRC)

      Character*(*)     FSEQ
      Character*4096    CARG

      Logical           OPTL
      Logical           OPTU
      Logical           OPTP
      Logical           OPTD

      IRC=0

      OPTL=.FALSE.
      OPTU=.FALSE.
      OPTP=.FALSE.
      OPTD=.FALSE.

      NM=-1
      NW=60

      N1=Iargc()
      K1=0
      I2=1
      Do  I1=1,N1
         Call GetArg(I2,CARG)
         If(CARG(1:1).EQ.'-'.AND.CARG(2:2).NE.' ') then
            If(Index(CARG,'h').NE.0) go to 900
            If(Index(CARG,'l').NE.0) OPTL=.TRUE.
            If(Index(CARG,'u').NE.0) OPTU=.TRUE.
            If(Index(CARG,'p').NE.0) OPTP=.TRUE.
            If(Index(CARG,'d').NE.0) OPTD=.TRUE.
            If(Index(CARG,'M').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NM
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NM
               End if
            End if
            If(Index(CARG,'W').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NW
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NW
               End if
            End if
         Else if(CARG(1:2).EQ.'M=') then
            Read(CARG(3:),*,Err=900) NM
         Else if(CARG(1:2).EQ.'W=') then
            Read(CARG(3:),*,Err=900) NW
         Else
            FSEQ=CARG
            K1=K1+1
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
      End do

 20   If(K1.LT.1) then
         FSEQ='-'
      Else if(K1.GT.1) then
         Go to 900
      End if
      If(NM.LT.0) NM=-1
      If(NW.LE.0.OR.NW.GT.512) NW=60

 100  Return
 900  IRC=-1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine GetSeq(NSEQ,CSEQ,LSEQ,RCIO,LDES,NBSQ,IDMS,IRC)

      Include          'sterr.f'

      Character*01      CSEQ(*)
      Character*(*)     RCIO

      Character*512     RCIN
      Save              RCIN

      IRC=0
      LSEQ=0
      LDES=0

 1    Continue
      If(RCIN(1:1).NE.'>') then
         Read(NSEQ,'(A)',Err=901,End=902) RCIN
         Go to   1
      End if

      LDES=Lblnk(RCIN)
      RCIO(1:LDES)=RCIN(1:LDES)

 2    RCIN=' '
      Read(NSEQ,'(A)',Err=901,End=900) RCIN
 5    L1=Lblnk(RCIN)
      If(RCIN(1:1).EQ.'>') then
         Go to 90
      Else
         Do I1=1,L1
            K1=Ichar(RCIN(I1:I1))
            If((K1.GE.65.AND.K1.LE. 90).OR.
     *         (K1.GE.97.AND.K1.LE.122).OR.
     *         (RCIN(I1:I1).EQ.'-').OR.
     *         (RCIN(I1:I1).EQ.'.')) then
               LSEQ=LSEQ+1
               If(LSEQ.GT.IDMS) go to 903
               CSEQ(LSEQ)=RCIN(I1:I1)
            End if
         End do
      End if
      If(IRC.NE.-1) go to   2

 90   If(LSEQ.EQ.0) go to 904
      NBSQ=NBSQ+1

 100  Return

 110  Write(NERR,*) '       While processing sequence ',
     *   RCIO(1:LDES)
      Go to 100

* errors

 900  If(RCIN.EQ.' '.AND.LSEQ.EQ.0) then
         Write(NERR,*) 'Error: Unexpected end of sequence.'//
     *      ' Sequence has zero length.'
         IRC=1
      Else
         IRC=-1
         Go to 5
      End if
      Go to 110

 901  Write(NERR,*) 'Error: Unable to read PSA sequence file.'
      IRC=1
      Go to 100
 902  If(NBSQ.EQ.0) then
         Write(NERR,*) 'Error: Unable to find sequences in PSA file.'
         IRC=1
      Else
         IRC=-1
      End if
      Go to 100
 903  Write(NERR,*) 'Error: Sequence length exceeds buffer ',
     *   'size (',IDMS,').'
      IRC=1
      Go to 110
 904  Write(NERR,*) 'Error: Unexpected end of sequence.'//
     *   ' Sequence has zero length.'
      IRC=1
      Go to 110

      End
*----------------------------------------------------------------------*
      Subroutine UpdIPRF(CSEQ,LSEQ,LPRF,IPRF,IDMP,IRC)

      Include          'sterr.f'

      Character*01      CSEQ(*)
      Integer           IPRF(0:IDMP)

      IRC=0

      If(LPRF.EQ.0) then
         Do I1=1,LSEQ
            K1=Ichar(CSEQ(I1))
            If((K1.GE.65.AND.K1.LE. 90).OR.
     *         (CSEQ(I1).EQ.'-')) then
               LPRF=LPRF+1
            End if
         End do
         If(LPRF.GT.IDMP) go to 901
         Do I1=0,LPRF
            IPRF(I1)=0
         End do
      End if

      J1=0
      M1=0
      Do I1=1,LSEQ
         K1=Ichar(CSEQ(I1))
         If((K1.GE.65.AND.K1.LE. 90).OR.
     *      (CSEQ(I1).EQ.'-')) then
            IPRF(J1)=MAX(IPRF(J1),M1)
            J1=J1+1
            M1=0
         Else
            M1=M1+1
         End if
      End do
      IPRF(J1)=MAX(IPRF(J1),M1)

      If(J1.NE.LPRF) go to 900

 100  Return

* errors

 900  Write(NERR,*) 'Error: Conflicting sequence length after ',
     *   'insertion removal.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Sequence length exceeds buffer ',
     *   'size (',IDMP,').'
      IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine UpdSeq(CSEQ,LSEQ,LMSA,IPRF,LPRF,IDMP)

      Character*01      CSEQ(*)
      Integer           IPRF(0:IDMP)

      J1=LPRF
      M1=LMSA
      N1=0

      Do I1=LSEQ,1,-1
         K1=Ichar(CSEQ(I1))
         If((K1.GE.65.AND.K1.LE. 90).OR.
     *      (CSEQ(I1).EQ.'-')) then
            K1=IPRF(J1)-N1
            M1=M1-K1
            Do I2=M1+1,M1+K1
               CSEQ(I2)='.'
            End do
            If(J1.EQ.LPRF) then
               Do I2=M1+1,M1+N1
                  CSEQ(I2)=CSEQ(I2+K1)
               End do
               Do I2=M1+N1+1,M1+N1+K1
                  CSEQ(I2)='.'
               End do
            Else

* place dots at the center of insert region

               L1=(N1+1)/2
               Do I2=M1+1,M1+L1
                  CSEQ(I2)=CSEQ(I2+K1)
               End do
               Do I2=M1+L1+1,M1+L1+K1
                  CSEQ(I2)='.'
               End do
            End if
            CSEQ(M1)=CSEQ(I1)
            M1=M1-1
            N1=0
            J1=J1-1
         Else
            CSEQ(M1)=CSEQ(I1)
            M1=M1-1
            N1=N1+1
         End if
      End do
      Do I2=N1+1,IPRF(J1)
         CSEQ(M1)='.'
         M1=M1-1
      End do

 100  Return
      End
*----------------------------------------------------------------------*
      Subroutine EdtSeq(CSEQ,LMSA,OPTL,OPTU,OPTP,OPTD)

      Character*01      CSEQ(*)

      Logical           OPTU
      Logical           OPTL
      Logical           OPTP
      Logical           OPTD

      If(OPTU) then
         Do  I1=1,LMSA
            K1=Ichar(CSEQ(I1))
            If(K1.LE.122.AND.K1.GE.97) CSEQ(I1)=Char(K1-32)
         End do
      End if

      If(OPTL) then
         Do  I1=1,LMSA
            K1=Ichar(CSEQ(I1))
            If(K1.LE.90.AND.K1.GE.65) CSEQ(I1)=Char(K1+32)
         End do
      End if

      If(OPTP) then
         Do  I1=1,LMSA
            If(CSEQ(I1).EQ.'-') CSEQ(I1)='.'
         End do
      End if

      If(OPTD) then
         Do  I1=1,LMSA
            If(CSEQ(I1).EQ.'.') CSEQ(I1)='-'
         End do
      End if

      Return
      End
*----------------------------------------------------------------------*
      Subroutine CutSeq(CSEQ,LMSA,LMSB,IPRF,LPRF,NM,IDMP)

      Character*01      CSEQ(*)
      Integer           IPRF(0:IDMP)

      M1=NM/2
      L1=NM-M1
      J1=0
      K1=0
      Do I1=0,LPRF
         If(IPRF(I1).GT.NM) then
            Do I2=1,M1
               J1=J1+1
               K1=K1+1
               CSEQ(J1)=CSEQ(K1)
            End do
            K1=K1+IPRF(I1)-NM
            Do I2=1,L1
               J1=J1+1
               K1=K1+1
               CSEQ(J1)=CSEQ(K1)
            End do
         Else
            Do I2=1,IPRF(I1)
               J1=J1+1
               K1=K1+1
               CSEQ(J1)=CSEQ(K1)
            End do
         End if
         If(I1.NE.LPRF) then
            J1=J1+1
            K1=K1+1
            CSEQ(J1)=CSEQ(K1)
         End if
      End do
      LMSB=J1

      Return
      End
*----------------------------------------------------------------------*
      Subroutine PRSQ(CSEQ,LMSB,NW)

      Character*01      CSEQ(*)
      Character*512     RCOUT

      INB=LMSB/NW
      INR=LMSB-INB*NW
      IN1=0

      Do IN2=1,INB
         Do ii1=1,NW
            RCOUT(ii1:ii1)=CSEQ(IN1+ii1)
         End do
         IN1=IN1+NW
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,NW)
      End do
      If(INR.GT.0) then
         Do ii1=1,INR
            RCOUT(ii1:ii1)=CSEQ(IN1+ii1)
         End do
         Write(6,'(512A)')(RCOUT(ii1:ii1),ii1=1,INR)
      End if

 100  Return
      End

*----------------------------------------------------------------------*
      Include          'lblnk.f'
