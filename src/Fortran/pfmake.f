*       Program pfmake
*----------------------------------------------------------------------*
* $Id: pfmake.f,v 2.11 2003/11/28 11:53:33 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: Constructs a profile from a multiple sequence
*                 alignment
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*
* DATA
*----------------------------------------------------------------------*

* array sizes, I/O units

      Include          'ardim.f'

C      Parameter        (IDM1=1048576)
C      Parameter        (IDM2=   2048)

      Parameter        (NOUT=      6)
      Parameter        (NMSF=      5)
      Parameter        (NCMP=     12)
      Parameter        (NPRF=     13)

* profile and sequence fields :

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'sterr.f'

C      Character*64      FOUT

* weights and distances

      Real*4            RWGT(IDMF)

* multiple sequence alignment:

      Character*4096    FMSF

      Character*64      SQID(IDMF)
      Character         CSEQ(IDMS)

* frequency profile

      Real              RPRF(26,IDMP)
      Real              FUNK(IDMP)
      Real              FRES(IDMP)
      Real              FDEL(IDMP)

* score profile

      Real              SPRF(26,IDMP)
      Real              SDEL(IDMP)
      Real              SINS(0:IDMP)

* symbol comparison table

      Character*4096    FCMP

      Real              RCMP(26,26)

* profile

      Character*4096    FPRF
      Character*30      FDAT


* function return types

      Integer           Xblnk
      External          Xblnk

* options

      Logical           OPT0
      Logical           OPT1
      Logical           OPT2
      Logical           OPT3
      Logical           OPTC
      Logical           LLLT
      Logical           OPTM
      Logical           LRNM
      Logical           LEOF

      Logical           LBLK
      Logical           LSYM
      Logical           LWGE
      Real              RE
      Real              RF
      Real              RG
      Real              RH
      Real              RI
      Real              RL
      Real              RM
      Real              RS
      Real              RT
      Real              RX

* character translation

      Character*27      ABCU
      Character*27      ABCL

      Data              ABCU/'ABCDEFGHIJKLMNOPQRSTUVWXYZ-'/
      Data              ABCL/'abcdefghijklmnopqrstuvwxyz.'/

* initialization of controlled vocabularies

      Include          'cvini.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0
      IDUM=-2
      BLOG=0.0
      LEOF=.FALSE.
      LRNM=.FALSE.

* read command line

      Call Repar
     *   (FMSF,FCMP,FPRF,OPT0,OPT1,OPT2,OPT3,OPTC,LSYM,LWGE,LBLK,
     *   LLLT,OPTM,RE,RF,RG,RH,RI,RL,RM,RS,RT,RX,NLOW,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''pfmake 2.3 revision 5.d'',//
     *      ''Usage: pfmake [ -0123abcehlsEFGHILMSTX ] [ msf-file'',
     *      '' | - ] score-matrix [ profile-file ] [ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -0: global alignment mode.'',/
     *      ''    -1: domain global alignment mode.'',/
     *      ''    -2: semi-global alignment mode.'',/
     *      ''    -3: local alignment mode.'',/
     *      ''    -a: asymmetrical gap weighting.'',/
     *      ''    -b: block profile mode.'',/
     *      ''    -c: circular profile.'',/
     *      ''    -e: enable endgap-weighting mode.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -m: input sequences in MSA format.'',/
     *      ''    -l: do not impose limit on line length.'',/
     *      ''    -s: symmetrical gap weighting.''
     *      )')
         Write(NERR,'(
     *      ''    -E<value>:'',/
     *      ''        gap extension penalty (default: 0.2).'',/
     *      ''    -F<value>:'',/
     *      ''        output score multiplier (default: 100)'',/
     *      ''    -G<value>:'',/
     *      ''        gap opening penalty (default: 2.1)'',/
     *      ''    -H<value>:'',/
     *      ''        high cost initiation/termination score (default'',
     *      '': *)'',/
     *      ''    -I<value>:'',/
     *      ''        gap penalty multiplier increment (default:'',
     *      '' 0.1)''
     *      )')
         Write(NERR,'(
     *      ''    -L<value>:'',/
     *      ''        low cost initiation/termination score (default:'',
     *      ''  0).'',/
     *      ''    -M<value>:'',/
     *      ''        maximum gap penalty multiplier (default:'',
     *      '' 0.333).'',/
     *      ''    -S<value>:'',/
     *      ''        score matrix multiplier (default: 0.1)'',/
     *      ''    -T<value>:'',/
     *      ''        gap region threshold (default: 0.01)'',/
     *      ''    -X<value>:'',/
     *      ''        gap excision threshold (default: 0.5)'',/
     *      )')
         Write(NERR,'(
     *      '' valid (but deprecated) parameters are:'',/,
     *      ''  [E=gap-extension-weight]      use option -E instead'',/
     *      ''  [F=output-score-multiplier]   use option -F instead'',/
     *      ''  [G=gap-weigth]                use option -G instead'',/
     *      ''  [H=high-cost-init-term-score] use option -H instead'',/
     *      ''  [I=Ginc-multiplier]           use option -I instead'',/
     *      ''  [L=low-cost-init-term-score]  use option -L instead'',/
     *      ''  [M=Gmax-multiplier]           use option -M instead'',/
     *      ''  [S=score-matrix-multiplier]   use option -S instead'',/
     *      ''  [T=gap-region-threshhold]     use option -T instead'',/
     *      ''  [X=gap-excision-threshold]    use option -X instead'',/
     *      )')
         Call Exit(IRC)
      End if

* read msf or msa-file

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

* read score-matrix

      Call RECMP(NERR,NCMP,FCMP,NABC,CABC,RCMP,IRC)
      If(IRC.NE.0) go to 100

* read input profile

      If(FPRF.NE.' ')
     *   Call REPRF
     *   (NPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   BLOG,FABC,P0,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,CHIL,IIPL,ILIP,
     *   CHID,IIPD,CHMD,IMPD,
     *   INBP,LEOF,.FALSE.,IRC)

      If(IRC.GT.0) go to 100
      If(IRC.NE.0.OR.FPRF.EQ.' ') then
         IRC=0
         CPID=' '
         CPAC=' '
         CPDT=' '
         CPDE=' '
         MDIS=0
         JNOR=0
         JCUT=0
         LHDR=0
         LFTR=0
      End if

*----------------------------------------------------------------------*
* DATA PROCESSING SECTION
*----------------------------------------------------------------------*

* normalize weights

      R1=0.0
      Do I1=1,NSEQ
         R1=R1+RWGT(I1)
      End do
      If(R1.EQ.0.0) then
         R1=1/Real(NSEQ)
         Do I1=1,NSEQ
            RWGT(I1)=R1
         End do
      Else
         Do I1=1,NSEQ
            RWGT(I1)=RWGT(I1)/R1
         End do
      End if

* make sequences uppercase

      Do I1=1,LSEQ*NSEQ
         If     (Index(ABCU,CSEQ(I1)).EQ.0) then
            IX=Index(ABCL,CSEQ(I1))
            If(IX.GT.0) then
               CSEQ(I1)=ABCU(IX:IX)
            Else
               CSEQ(I1)='-'
            End if
         End if
      End do

* remove endgaps

      If(.NOT.LWGE) then
         J1=0
         Do  10 I1=1,NSEQ
            Do   6 I2=J1+1,J1+LSEQ
               If(CSEQ(I2).NE.'-') go to  7
               CSEQ(I2)=' '
 6          Continue
 7          Continue

            Do   8 I2=J1+LSEQ,J1+1,-1
               If(CSEQ(I2).NE.'-') go to  9
               CSEQ(I2)=' '
 8          Continue
 9          Continue
            J1=J1+LSEQ
 10      Continue
      End if

* make frequency profile
*
*   for residue i and profile position j:
*      RPRF(i,j) = residue frequency(residue,position)
*   for profile position j:
*      FUNK(  j) = fraction of unknown residues
*      FRES(  j) = fraction of residues and gaps
*      RDEL(  l) = fraction of gaps

* - initialize

      Do  12  I1=1,LSEQ
         Do  11 I2=1,NABC
            RPRF(I2,I1)=0.0
 11      Continue
         FUNK(I1)=0.0
         FRES(I1)=0.0
         FDEL(I1)=0.0
 12   Continue

      J1=0
      K1=1
      Do  20 I1=1,NSEQ
         K2=1
         Do  19 I2=J1+1,J1+LSEQ
            IX=0
            Do  15 I3=1,NABC
               If(CSEQ(I2).EQ.CABC(I3)) then
                  IX=I3
                  Go to  16
               End if
 15         Continue
 16         Continue

            FRES(   K2)=FRES(   K2)+RWGT(K1)
            If     (IX.NE.0) then
               RPRF(IX,K2)=RPRF(IX,K2)+RWGT(K1)
            else if(CSEQ(I2).EQ.'-') then
               FDEL(   K2)=FDEL(   K2)+RWGT(K1)
            else if(CSEQ(I2).EQ.' ') then
               FRES(   K2)=FRES(   K2)-RWGT(K1)
            else
               FUNK(   K2)=FUNK(   K2)+RWGT(K1)
            end if
            K2=K2+1
 19      Continue
         J1=J1+LSEQ
         K1=K1+1
 20   Continue

* normalize frequencies:

      Do  25 I1=1,LSEQ
         FDEL(I1)=FDEL(I1)/FRES(I1)
         R1=FRES(I1)-FUNK(I1)
         If(R1.EQ.0) go to  25
         Do  23 I2=1,NABC
            RPRF(I2,I1)=RPRF(I2,I1)/R1
 23      Continue
 25   Continue

* construct score profile:
*
*   for known residue i and profile match position j:
*      SPRF(i,j) = match extension score
*   for profile position j:
*      SDEL(  j) = gap penalty multipliers for deletion gaps
*      SINS(  j) = gap penalty multipliers for insertion gaps

* - gap penalty multipliers:

      K1=0
      J1=0
      SINS(0)=1.0
      Do   30 I1=1,LSEQ
         If(FDEL(I1).LE.RT) then
            If(K1.GT.0) then
               XDEL=RM/(1.0+K1*RI)
               SINS(J1-1)=XDEL
               Do  27 I2=J1,I1-1
                  SINS(I2)=XDEL
                  SDEL(I2)=XDEL
 27            Continue
               K1=0
            End if
            SINS(I1)=1.0
            SDEL(I1)=1.0
         else
            If(K1.EQ.0) J1=I1
            K1=K1+1
         end if
 30   Continue

      K1=0
      J1=0
      Do  50 I1=1,LSEQ
         If(FDEL(I1).GT.1-RX) then

* - gap excision

            J1=J1+1
         else
            K1=K1+1

* - score matrix application

            Do  33 I2=1,NABC
               SPRF(I2,K1)=0.0
               Do  32 I3=1,NABC
                  SPRF(I2,K1)=SPRF(I2,K1)+RCMP(I3,I2)*RPRF(I3,I1)
 32            Continue
               SPRF(I2,K1)=RS*SPRF(I2,K1)
C                Write(6,'(I4,1x,A,F10.4)') K1,CABC(I2),SPRF(I2,K1)
 33         Continue
            If(J1.GT.0) then
               SINS(K1-1)=-SINS(I1-1)
               J1=0
            End if
            SDEL(K1)=SDEL(I1)
            SINS(K1)=SINS(I1)
            FDEL(K1)=FDEL(I1)
         end if
 50   Continue
      If(J1.GT.0) SINS(K1)=-SINS(K1)
      LPRF=K1


* initialize generalized profile

* - header

      If(CPID.EQ.' ') CPID='SEQUENCE_PROFILE'
      If(CPAC.EQ.' ') CPAC='ZZ99999'
      CALL Fdate(FDAT)
      If(CPDT.NE.' ') then
         CPDT=FDAT(1:Lblnk(FDAT)) // ' ! ' // CPDT
      Else
         CPDT=FDAT(1:Lblnk(FDAT))
      End if
      If(FMSF.EQ.'-') FMSF='stdin'
      If(CPDE.EQ.' ') CPDE='Generated from MSF file: '''
     *   // FMSF(1:Lblnk(FMSF))
     *   // '''.'

* - accessories

      If(MDIS.EQ.0) then
         MDIS=2
         N1=MIN(5,LPRF/10)
         NDIP(1)=1   +N1
         NDIP(2)=LPRF-N1
      End If

      If(JNOR.EQ.0.OR.JCUT.EQ.0) then
         JNOR=1
         MNOR(1)=1
         NNOR(1)=1
         NNPR(1)=1
         CNTX(1)='No_units'
         RNOP(1,1)=0.0
         RNOP(2,1)=1/RF
      End if

      If(JCUT.EQ.0) then
         JCUT=2

         MCLE(1)=0
         CCUT(1)='!'
         ICUT(1)=INT(RF*8.5)
         JCNM(1)=1
         RCUT(1,1)=8.5
         MCUT(1,1)=1

         MCLE(2)=-1
         CCUT(2)='?'
         ICUT(2)=INT(RF*6.5)
         JCNM(2)=1
         RCUT(1,2)=6.5
         MCUT(1,2)=1
      End if

* - defaults for match and insert position

      CHID='-'
      Do  65 I1=1,26
         IIPD(I1)=0
 65   Continue

      If(LSYM) then
         N1=-NINT(RF*RG/2)
         N2=-NINT(RF*RG/2)
      Else
         N1=-NINT(RF*RG)
         N2=0
      End if
      N3=-NINT(RF*RE)

      IIPD(B0)=0
      IIPD(B1)=NLOW
      IIPD(E0)=0
      IIPD(E1)=NLOW

      IIPD(BM)=0
      IIPD(BI)=NLOW
      IIPD(BD)=NLOW
      IIPD(BE)=NLOW
      IIPD(MM)=0
      IIPD(MI)=N1
      IIPD(MD)=N1
      IIPD(ME)=0
      IIPD(IM)=N2
      IIPD(II)=0
      IIPD(ID)=NLOW
      IIPD(IE)=NLOW
      IIPD(DM)=N2
      IIPD(DI)=NLOW
      IIPD(DD)=0
      IIPD(DE)=NLOW

      IIPD(I0)=0

      CHMD='X'
      CHID='-'
      Do  66 I1=1,26
         IIPD(I1)=N3
         IMPD(I1)=0
 66   Continue

      IIPD(M0)=0

      IMPD(D )=N3

      Do  68 I1=0,27
         IMPP(I1,0)=NLOW
 68   Continue

* build insert positions

      Do  75 I1=0,LPRF

         CHIP(I1)=CHID

         Do I2=0,46
            IIPP(I2,I1)=IIPD(I2)
         End do

         If(SINS(I1).LT.0) then
            XINS=0.0
            SINS(I1)=-SINS(I1)
         Else
            XINS=SINS(I1)
         End if

         NINS=-NINT(RF*RE*SINS(I1))
         Do  I2=1,NABC
            IIPP(I2,I1)=NINS
         End do

         If     (LSYM) then
            IIPP(MI,I1)=-NINT(RF*RG*XINS      /2)
            IIPP(IM,I1)=-NINT(RF*RG*XINS      /2)
            If(XINS.EQ.0) then
               If(I1.NE.LPRF) IIPP(MD,I1)=-NINT(RF*RG*SINS(I1)  /2)
               If(I1.NE.   0) IIPP(DM,I1)=-NINT(RF*RG*SINS(I1)  /2)
            Else
               If(I1.NE.LPRF) IIPP(MD,I1)=-NINT(RF*RG*SDEL(I1+1)/2)
               If(I1.NE.   0) IIPP(DM,I1)=-NINT(RF*RG*SDEL(I1  )/2)
            End if
         Else if(I1.NE.LPRF) then
            IIPP(MI,I1)=-NINT(RF*RG*SDEL(I1+1))
            IIPP(MD,I1)=-NINT(RF*RG*SDEL(I1+1))
         End if

 75   Continue

* - build match position

      Do  80 I1=1,LPRF
         Do  I2=0,27
            IMPP(I2,I1)=IMPD(I2)
         End do
         IMPP( D,I1)=-NINT(RF*RE*SDEL(I1))
         Do  I2=1,NABC
            IMPP(I2,I1)=NINT(RF*SPRF(I2,I1))
         End do

         J2=0
         K2=0
         Do I2=1,NABC
            If(IMPP(I2,I1).GT.K2) then
               J2=I2
               K2=IMPP(I2,I1)
            End if
         End do
         If(J2.NE.0) then
            CHMP(I1)=CABC(J2)
         Else
            CHMP(I1)=CHMD
         End if
 80   Continue

* block profile

      If(LBLK) then
         K1=0
         J1=0
 85      J1=J1+1
         If((FDEL(J1).LE.RT.AND.K1+1.NE.J1).OR.
     *      (FDEL(J1).GT.RT.AND.J1.EQ.LPRF)) then

* - gap region: define center position

            N1=K1
            N3=J1-1
            If(FDEL(J1).GT.RT) N3=LPRF

* - - beginning, end

            If(N1.EQ.0) then
               N2=N3
               GO to 90
            End if

            If(N3.EQ.LPRF) then
               N2=N1
               Go to 90
            End if

* - - gap excision position ?

            N2=N1
            Do I1=N1,N3
               If(IIPP(MI,I1).EQ.0) then
                  N2=I1
                  Go to  90
               End if
            End do

* - - find center position

            R1=0
            R2=0
            Do I1=N1+1,N3
               R1=R1+FDEL(I1)-FDEL(I1-1)
               If(R1.GT.R2) then
                  N2=I1
                  R2=R1
               End if
            End do

* - edit gap region

 90         Continue
            Do I1=N1,N3
               If(I1.LT.N2) then
                  IIPP(DM,I1)=IIPD(DM)
                  IIPP(IM,I1)=IIPD(IM)
                  IIPP(MI,I1)=IIPD(MI)
               End if
               If(I1.GT.N2) then
                  IIPP(MD,I1)=IIPD(MD)
                  IIPP(IM,I1)=IIPD(IM)
                  IIPP(MI,I1)=IIPD(MI)
               End if
            End do

* - gap region done

            K1=J1

* - no gap region

         Else
            If(FDEL(J1).LE.RT) K1=J1
         End if
         If(J1.LT.LPRF) go to 85

      End if

* begin-to, to-end scores

      IIPP(BM,   0)=IIPP(MM,   0)
      IIPP(BI,   0)=IIPP(MI,   0)
      IIPP(BD,   0)=IIPP(MD,   0)
      IIPP(ME,LPRF)=IIPP(MM,LPRF)
      IIPP(IE,LPRF)=IIPP(IM,LPRF)
      IIPP(DE,LPRF)=IIPP(DM,LPRF)

* alignment mode

* - define high and low init/term scores in integer units

      NH=NLOW
      If(RH.GT.Real(NLOW)) then
         RH=-RH*RF
         If(RH.GT.Real(NLOW)) NH=NINT(RH)
      End if
      NL=NLOW
      If(-RF*RL.GT.Real(NLOW)) NL=NINT(-RL*RF)

* - global alignment mode (also initialization for other modes)

      Do I1=0,LPRF
         IIPP(B0,  I1)=NH
         IIPP(B1,  I1)=NH
         IIPP(E0,  I1)=NH
         IIPP(E1,  I1)=NH
         IIPP(BE,  I1)=IIPD(BE)
      End do
      IIPP(B0,   0)=NL
      IIPP(E0,LPRF)=NL

      IIPD(B0)=NH
      IIPD(B1)=NH
      IIPD(E0)=NH
      IIPD(E1)=NH


* - domain-global alignment mode

      If(OPT1.OR.OPT2.OR.OPT3) then
         IIPP(B1,   0)=NL
         IIPP(E1,LPRF)=NL
      End if

* - semiglobal alignment mode

      If(OPT2.OR.OPT3) then
         Do I1=0,LPRF
            IIPP(B0,I1)=NL
            IIPP(E0,I1)=NL
         End do
         IIPD(B0)=NL
         IIPD(E0)=NL
      End if

* - local alignment mode

      If(OPT3) then
         Do I1=1,LPRF-1
            IIPP(B1,I1)=NL
            IIPP(E1,I1)=NL
         End do
         IIPD(B1)=NL
         IIPD(E1)=NL
      End if

* - default M0 value:

      R=0.0
      Do I1=1,LPRF
         Do I2=1,NABC
            R=R+IMPP(I2,I1)
         End do
      End do
      IMPD(M0)=INT(R/(LPRF*NABC))

      Do I1=1,LPRF
         IMPP(M0,I1)=IMPD(M0)
      End do

* linear or circular profile ?

      If(OPTC) then
         Do I1=0,46
            K1=MAX(IIPP(I1,   0),IIPP(I1,LPRF))
            IIPP(I1,   0)=K1
            IIPP(I1,LPRF)=K1
         End do
         IIPP(MM,0)=MAX(IIPP(MM,0),IIPP(BM,0),IIPP(ME,0))
         IIPP(MI,0)=MAX(IIPP(MI,0),IIPP(BI,0))
         IIPP(MD,0)=MAX(IIPP(MD,0),IIPP(BD,0))
         IIPP(IM,0)=MAX(IIPP(IM,0),IIPP(IE,0))
         IIPP(DM,0)=MAX(IIPP(DM,0),IIPP(DE,0))
         LPCI=.TRUE.
      Else
         LPCI=.FALSE.
      End if

*----------------------------------------------------------------------*
* OUTPUT SECTION
*----------------------------------------------------------------------*

C      FOUT='stdout'

* add command-line to footer lines

      LFTR=LFTR+1
      If(LFTR.GT.1024) LFTR=1024
      Do I1=LFTR,2,-1
         CFTR(I1)=CFTR(I1-1)
      End do

      CFTR(1)='CC   /GENERATED_BY="'
      Call Recmd(CFTR(1)(21:130))
      IC=Lblnk(CFTR(1))
      CFTR(1)(IC+1:)='";'

      Call WRPRF
     *   (NOUT,LLLT,LRNM,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   BLOG,FABC,P0,
     *   CHID,IIPD,CHMD,IMPD,
     *   IRC)

 100  Call Exit(IRC)
      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FMSF,FCMP,FPRF,OPT0,OPT1,OPT2,OPT3,OPTC,LSYM,LWGE,LBLK,
     *   LLLT,OPTM,RE,RF,RG,RH,RI,RL,RM,RS,RT,RX,NLOW,IRC)

      Character*4096    CPAR
      Character*(*)     FMSF
      Character*(*)     FCMP
      Character*(*)     FPRF

      Logical           OPT0
      Logical           OPT1
      Logical           OPT2
      Logical           OPT3
      Logical           OPTC
      Logical           LLLT
      Logical           OPTM

      Logical           LBLK
      Logical           LSYM
      Logical           LWGE

      OPT0=.FALSE.
      OPT1=.FALSE.
      OPT2=.TRUE.
      OPT3=.FALSE.
      OPTC=.FALSE.
      LLLT=.TRUE.
      OPTM=.FALSE.

      LBLK=.FALSE.
      LSYM=.TRUE.
      LWGE=.FALSE.

      RE=0.2
      RF=100
      RG=2.1
      RH=Real(NLOW)
      RI=0.1
      RL=0.0
      RM=0.333
      RS=0.1
      RT=0.01
      RX=0.5

      IRC=0
      FMSF=' '
      FCMP=' '
      FPRF=' '

      N1=Iargc()

      K1=0
      I2=1
      Do  10 I1=1,N1
         Call GetArg(I2,CPAR)
         If     (CPAR(1:1).EQ.'-'.AND.CPAR(2:2).NE.' ') then
            If(Index(CPAR,'0').NE.0) then
               OPT0=.TRUE.
               OPT2=.FALSE.
            else if(Index(CPAR,'1').NE.0) then
               OPT1=.TRUE.
               OPT2=.FALSE.
            else if(Index(CPAR,'2').NE.0) then
               OPT2=.TRUE.
            else if(Index(CPAR,'3').NE.0) then
               OPT3=.TRUE.
               OPT2=.FALSE.
            End if
            If(Index(CPAR,'h').NE.0) then
               IRC=1
               Goto 100
            End if
            If(Index(CPAR,'a').NE.0) LSYM=.FALSE.
            If(Index(CPAR,'b').NE.0) LBLK=.TRUE.
            If(Index(CPAR,'c').NE.0) OPTC=.TRUE.
            If(Index(CPAR,'e').NE.0) LWGE=.TRUE.
            If(Index(CPAR,'l').NE.0) LLLT=.FALSE.
            If(Index(CPAR,'m').NE.0) OPTM=.TRUE.
            If(Index(CPAR,'s').NE.0) LSYM=.TRUE.
            If(Index(CPAR,'E').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RE
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RE
               End if
            End if
            If(Index(CPAR,'F').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RF
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RF
               End if
            End if
            If(Index(CPAR,'G').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RG
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RG
               End if
            End if
            If(Index(CPAR,'H').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RH
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RH
               End if
            End if
            If(Index(CPAR,'I').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RI
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RI
               End if
            End if
            If(Index(CPAR,'L').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RL
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RL
               End if
            End if
            If(Index(CPAR,'M').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RM
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RM
               End if
            End if
            If(Index(CPAR,'S').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RS
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RS
               End if
            End if
            If(Index(CPAR,'T').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RT
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RT
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

         Else if(CPAR(1:2).EQ.'E=') then
            Read(CPAR(3:),*,Err=900) RE
         Else if(CPAR(1:2).EQ.'F=') then
            Read(CPAR(3:),*,Err=900) RF
         Else if(CPAR(1:2).EQ.'G=') then
            Read(CPAR(3:),*,Err=900) RG
         Else if(CPAR(1:2).EQ.'H=') then
            Read(CPAR(3:),*,Err=900) RH
         Else if(CPAR(1:2).EQ.'I=') then
            Read(CPAR(3:),*,Err=900) RI
         Else if(CPAR(1:2).EQ.'L=') then
            Read(CPAR(3:),*,Err=900) RL
         Else if(CPAR(1:2).EQ.'M=') then
            Read(CPAR(3:),*,Err=900) RM
         Else if(CPAR(1:2).EQ.'S=') then
            Read(CPAR(3:),*,Err=900) RS
         Else if(CPAR(1:2).EQ.'T=') then
            Read(CPAR(3:),*,Err=900) RT
         Else if(CPAR(1:2).EQ.'X=') then
            Read(CPAR(3:),*,Err=900) RX

         Else if(FMSF.EQ.' ') then
            FMSF=CPAR
         Else if(FCMP.EQ.' ') then
            FCMP=CPAR
         Else if(FPRF.EQ.' ') then
            FPRF=CPAR
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
 10   Continue

 20   If(FMSF.EQ.' '.OR.FCMP.EQ.' ') go to 900

      If(LBLK) then
         LSYM=.TRUE.
         LWGE=.FALSE.
      End if

 100  Return
 900  IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine RECMP(NERR,NCMP,FCMP,NABC,CABC,RCMP,IRC)

      Character*(*)    FCMP

      Integer          NABC
      Character        CABC(0:26)

      Real             RCMP(26,26)

      Character*512    RCIN
      Character*32     CMIS


      Open(NCMP,File=FCMP,Status='OLD',Err=900)

      CMIS='..'
 1    Read(NCMP,'(A)',Err=901,End=902) RCIN
      IX=Index(RCIN,'..')
      If(IX.EQ.0) go to   1

* read alphabet

      NABC=0
      Do  10 I1=1,IX-1
         If(RCIN(I1:I1).NE.' ') then
            NABC=NABC+1
            CABC(NABC)=RCIN(I1:I1)
         End if
 10   Continue

      If(NABC.LE.0.OR.NABC.GT.26) go to 905
      Do  20  I1=1,NABC
 11      Read(NCMP,'(A)',Err=901,End=903) RCIN
         If(RCIN.EQ.' ') go to  11
         Read(RCIN,*,Err=904)(RCMP(I1,ii1),ii1=I1,NABC)
         Read(RCIN,*,Err=904)(RCMP(ii1,I1),ii1=I1,NABC)
 20   Continue

 100  Return

* errors

 900  Write(NERR,*) 'Error: Unable to open matrix file'//
     *   ' ''',FCMP(1:Lblnk(FCMP)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to read matrix file'//
     *   ' ''',FCMP(1:Lblnk(FCMP)),'''.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Unexpected end of matrix file. '//
     *   'Unable to find ''',CMIS,''' keyword.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Unexpected end of matrix file. '//
     *   'Discrepancies between alphabet length and matrix lines.'
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Unable to read matrix values.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
 905  Write(NERR,*) 'Error: Unable to handle matrix alphabet size.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Include          'remsf.f'
      Include          'remsa.f'
      Include          'recmd.f'
      Include          'reprf.f'
      Include          'wrprf.f'
      Include          'lblnk.f'
      Include          'Xblnk.f'

