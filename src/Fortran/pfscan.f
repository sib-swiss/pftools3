*       Program pfscan
*----------------------------------------------------------------------*
* $Id: pfscan.f,v 2.17 2003/12/09 13:42:42 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: Scan a DNA or protein sequences with a profile library
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*
* DATA
*----------------------------------------------------------------------*

* array dimensions and I/O units

      Include          'ardim.f'

      Parameter        (NOUT=   6)

      Parameter        (NPRF=  11)
      Parameter        (NSEQ=  12)

* profile

      Character*4096    FPRF

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'pxdat.f'
      Include          'avdat.f'
      Include          'sterr.f'

      Logical           LUNI
      Logical           LNOR
      Logical           LREV
      Logical           LTRA
      Logical           LPFA
      Logical           LEOF
      Logical           LOUT

      Character         DABC(0:26)

* sequence

      Character*4096    FSEQ

      Character*64      CSID
      Character*64      CSAC
      Character*64      CSFH
      Character*512     CSDE

      Integer           LSEQ
      Integer           BSEQ

      Integer*2         ISEQ(IDMS)
      Character         CSEQ(IDMS)

      Logical           LCKS(IDMS)

      Integer*2         IS

* number of sequences, profiles read

      Integer           INBS
      Integer           INBP


* function return types

      Integer           Xblnk
      External          Xblnk

* options and command line parameters

      Logical           OPTA
      Logical           OPTB
      Logical           OPTF
      Logical           OPTL
      Logical           OPTR
      Logical           OPTS
      Logical           OPTU
      Logical           OPTX
      Logical           OPTY
      Logical           OPTZ
      Logical           OPLU
      Logical           OPTM
      Logical           OPTK
      Logical           OPTO
      Logical           OPTD
      Logical           OPTV

      Logical           LCMM
      Logical           OOPR
      Logical           LDRS

* alignments

      Integer           NALI
      Integer           IALS(IDMN)
      Integer           IALB(IDMN)
      Integer           IAL1(IDMN)
      Integer           IAL2(IDMN)
      Integer           IALE(IDMN)

      Integer           LALI
      Character         CALI(IDMA)
      Character         CPMA(IDMM)

* multiple matches of circular profile

      Integer           NMAT(IDMN)
C      Integer           IMMB(IDMN,IDML)
C      Integer           IMME(IDMN,IDML)
      Integer           PK1E(IDML)
      Integer           PK1B(IDML)
      Integer           PK2E(IDML)
      Integer           PK2B(IDML)
      Integer           PK3E(IDML)
      Integer           PK3B(IDML)
      Integer           PJSE(IDML)
      Integer           PJ1E(IDML)
      Integer           PJ1B(IDML)
      Integer           PLAL(IDML)
      Integer           IMSC(IDML)

* path matrix fields

      Integer           IOPM(0:IDMP)
      Integer           IOPI(0:IDMP)
      Integer           IOPD(0:IDMP)

      Integer           IOMB(0:IDMP)
      Integer           IOM1(0:IDMP)
      Integer           IOM2(0:IDMP)

      Integer           IOIB(0:IDMP)
      Integer           IOI1(0:IDMP)
      Integer           IOI2(0:IDMP)

      Integer           IODB(0:IDMP)
      Integer           IOD1(0:IDMP)
      Integer           IOD2(0:IDMP)

* work fields
C      Character*256     RCIN
      Logical           SOPM

* initialization of controlled vocabularies

      Include          'cvini.f'
      Include          'abini.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0

      LUNI=.FALSE.
      LNOR=.FALSE.
      LREV=.FALSE.
      LTRA=.FALSE.
      LPFA=.FALSE.
      LEOF=.FALSE.
      LDRS=.FALSE.

      LEOF=.FALSE.
      CABC(0)='-'

      INBS=0
      INBP=0
      BSEQ=1

* read command line arguments

      Call Repar(
     *   OPTA,OPTB,OPTF,OPTL,OPTR,OPTS,OPTU,OPLU,OPTX,OPTY,OPTZ,OPTM,
     *   OPTK,FPRF,FSEQ,LCUC,NW,NMOD,OPTO,OPTD,OPTV,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''pfscan 2.3 revision 5.d'',//
     *      ''Usage: pfscan[ -abCdfhlLmMkrsuvWxyz ] [ seq-file'',
     *      '' | - ] [ profile-library-file | - ] [ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -a: report optimal alignment for all profiles.'',/
     *      ''    -b: search complementary strand of DNA sequences.'',/
     *      ''    -f: input sequence file is in FASTA format.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -l: indicate highest cut-off level (number).'',/
     *      ''    -L: indicate highest cut-off level (text).'',/
     *      ''    -m: report individual matches for circular'',
     *      '' profiles.'',/
     *      ''    -r: use raw score.'',/
     *      ''    -u: force profile disjointness to UNIQUE.''
     *      )')
         Write(NERR,'(
     *      ''    -C<value>:'',/
     *      ''        cut-off level to be used for match selection.'',
     *      '' Same as parameter L.'',/
     *      ''    -M<value>:'',/
     *      ''        set the normalization mode to use for the'',
     *      '' score computation.'',/
     *      ''        Overrides the profile PRIORITY parameter.'',/
     *      )')
         Write(NERR,'(
     *      ''   output modifiers:'',/,
     *      ''    -d: impose length limit on profile description.'',/
     *      ''    -k: output using the xPSA header (using keyword'',
     *      ''=value pairs).'',/
     *      ''    -s: list sequences of the matched regions.'',/
     *      ''    -v: suppress warnings on stderr.'',/
     *      ''    -x: list alignments in PSA format.'',/
     *      ''    -y: list alignments in human readable form.'',/
     *      ''    -z: indicate profile start and stop positions.'',/
     *      ''    -W<value>:'',/
     *      ''        specifies the output width. Same as parameter'',
     *      '' W.'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',/,
     *      ''    [L=cut-off-level]  use option -C instead'',/
     *      ''    [W=output-width]   use option -W instead'',/
     *      )')
         Call Exit(IRC)
      End if


      If(FSEQ.EQ.'-') then
         MSEQ=5
      Else
         MSEQ=NSEQ
      End if

      If(FPRF.EQ.'-') then
         MPRF=5
      Else
         MPRF=NPRF
      End if

      If(OPTR) then
         LDRS=.TRUE.
      End if

      SOPM=OPTM

* read sequence

      If(OPTF) then
         Call RFSEQ
     *      (MSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,CSFH,LSEQ,ISEQ,LEOF,
     *      INBS,OPTV,IRC)
      Else
         Call RESEQ
     *      (MSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,LSEQ,ISEQ,LEOF,
     *      INBS,OPTV,IRC)
      End if
      If(IRC.GT.0) go to 100

      If(LEOF.AND.IRC.LT.0) then
         If(INBS.GE.1) IRC=0
         go to 100
      End if

* check parameter consistency

      If(.NOT.OPTV.AND.OPTA.AND.LCUC.NE.0) then
         Write(NERR,*) 'Warning: Option -a is set. Ignoring command',
     *      ' line cut-off level (option -C).'
         LCUC=0
      End if
      OOPR=OPTR

* determine amino acid composition

      Call CFAve(ISEQ,IDMS,BSEQ,LSEQ,CABC,NABC,FAVE)

* backtranslate sequence into characters

      Do  I1=1,LSEQ
         CSEQ(I1)=CABC(ISEQ(I1))
      End do

* alignment and ouptut format switches

      If(OPTX.OR.OPTY.OR.OPTZ.OR.OPTM) then
         LTRA=.TRUE.
      Else
         LTRA=.FALSE.
      End if
      If(OPTS.OR.OPTX) then
         LPFA=.TRUE.
      Else
         LPFA=.FALSE.
      End if

* compute alignments if any output format switch is specified

      LOUT=(LTRA.OR.LPFA)

*----------------------------------------------------------------------*
* major loop over profiles
*----------------------------------------------------------------------*

 1    Continue

* save alphabet

      OPTR=OOPR
      MABC=NABC
      Do  I1=1,NABC
         DABC(I1)=CABC(I1)
      End do

* read profile

      Call REPRF
     *   (MPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   BLOG,FABC,P0,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,CHIL,IIPL,ILIP,
     *   CHID,IIPD,CHMD,IMPD,
     *   INBP,LEOF,OPTV,IRC)

      If(IRC.NE.0.OR.LEOF) Go to 100

* check parameter consistency

      If(.NOT.OPTV.AND.OPTB.AND.NABC.GT.4) then
         Write(NERR,*) 'Warning: Not a DNA sequence. ',
     *      'Ignoring option -b.'
         OPTB=.FALSE.
      End if
      OPTM=SOPM
      If(.NOT.OPTV.AND.OPTM.AND..NOT.LPCI) then
         Write(NERR,*) 'Warning: Profile not circular. ',
     *      'Ignoring option -m.'
         Write(NERR,*) '         While processing profile ',
     *      CPID(1:Lblnk(CPID))
         OPTM=.FALSE.
      End if

* cut-off value

      LCUT=0
      Do   3 I1=1,JCUT
         If(   (MCLE(I1).GE.LCUC.AND.MCLE(I1).LT.LCUT.AND.LCUC.LT.LCUT).
     *      OR.(MCLE(I1).LE.LCUC.AND.MCLE(I1).GT.LCUT.AND.LCUC.GT.LCUT))
     *      LCUT=MCLE(I1)
 3    Continue

      If (.NOT.OPTV.AND.LCUT.NE.LCUC) then
         Write(NERR,*) 'Warning: Cut-off level',LCUC,' is not '//
     *      'defined. Using nearest level',LCUT,'.'
         Write(NERR,*) '         While processing profile ',
     *      CPID(1:Lblnk(CPID))
      End if

      KCUT=0

* find normalisation mode and cut-off level

      LNOR=.FALSE.
      If(OPTO) then
         J1=0
*   find the specified normalisation mode
         Do I1=1,JNOR
            If(NNOR(I1).EQ.NMOD) J1=I1
         End do
*   exit if mode number not specified in profile
         If(J1.EQ.0) Go to 901

*   find mode number in list of cut-off levels
         INOR=0
         Do I1=1,JCUT
            If(MCLE(I1).EQ.LCUT) then
               If(JCNM(I1).NE.0) then
                  J2=0
                  Do I2=1,JCNM(I1)
                     If(MCUT(I2,I1).EQ.NMOD) J2=I2
                  End do
                  If(J2.EQ.0) Go to 903
                  LNOR=.TRUE.
                  INOR=J1
                  MNUM=NNOR(J1)
                  IFUN=MNOR(J1)
                  KCUT=ICUT(I1)
                  XCUT=RCUT(J2,I1)
               Else
                  Go to 903
               End if
            End if
         End do
         If(INOR.EQ.0) Go to 901
      Else

* search cut-off level in profile if normalisation was not specified

         Do   6 I1=1,JCUT
            If(MCLE(I1).EQ.LCUT) then
               INOR=0
               If(JCNM(I1).NE.0) then
                  LNOR=.TRUE.
C                  J2=1
                  Do  5 I2=1,JCNM(I1)
                     J3=0
                     Do  4 I3=1,JNOR
                        If(MCUT(I2,I1).EQ.NNOR(I3)) J3=I3
 4                   Continue
                     If(J3.EQ.0) Go to 904

                     If     (I2.EQ.1) then
                        INOR=J3
                        MNUM=NNOR(J3)
                        NPRI=NNPR(J3)
                        IFUN=MNOR(J3)
                        KCUT=ICUT(I1)
                        XCUT=RCUT(I2,I1)
                     Else if(NNPR(J3).LT.NPRI) then
                        INOR=J3
                        MNUM=NNOR(J3)
                        NPRI=NNPR(J3)
                        IFUN=MNOR(J3)
                        KCUT=ICUT(I1)
                        XCUT=RCUT(I2,I1)
                     End if
 5                Continue
               End if

               If(JCNM(I1).EQ.0.OR.INOR.EQ.0) then
                  KCUT=ICUT(I1)
                  LNOR=.FALSE.
               End if
            End if
 6       Continue
      End if

      If(.NOT.LNOR) OPTR=.TRUE.
*      If(OPTR) LNOR=.FALSE.

* disjoint definition
* set disjointness position to begin and end of profile if necessary

      If(MDIS.EQ.1.OR.OPTU.OR.OPTA) then
         LUNI=.TRUE.
         If(MDIS.EQ.1) then
            NDIP(1)=1
            NDIP(2)=LPRF
         End if
      Else
         LUNI=.FALSE.
      End if

* - initialize profile lock

      If(.NOT.LPCI) then
         Do  8 I1=0,NDIP(1)-1
            IIPP(E0,I1)=NLOW
            IIPP(E1,I1)=NLOW
 8       Continue

         Do  9 I1=NDIP(2),LPRF
            IIPP(B0,I1)=NLOW
            IIPP(B1,I1)=NLOW
 9       Continue
      End if

C      End if

* profile extra parameters

      MLOW=NLOW/4*3
      Do  10 I1=0,LPRF
         IIPX( XM,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BM,I1))
         IIPX( XI,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BI,I1))
         IIPX( XD,I1) = MAX(MLOW,IIPP( B1,I1) + IIPP( BD,I1))
         IIPX( YM,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BM,I1))
         IIPX( YI,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BI,I1))
         IIPX( YD,I1) = MAX(MLOW,IIPP( B0,I1) + IIPP( BD,I1))
         IIPX( MX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( ME,I1))
         IIPX( IX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( IE,I1))
         IIPX( DX,I1) = MAX(MLOW,IIPP( E1,I1) + IIPP( DE,I1))
         IIPX( MY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( ME,I1))
         IIPX( IY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( IE,I1))
         IIPX( DY,I1) = MAX(MLOW,IIPP( E0,I1) + IIPP( DE,I1))
 10   Continue

* check alphabet

      If(NABC.EQ.MABC) then
         Do  I1=1,NABC
            If(CABC(I1).NE.DABC(I1)) go to  15
         End do
         Go to  21
      End if
 15   Continue

* reconvert sequence into numbers (if necessary)

      Do  20 I1=1,LSEQ
         ISEQ(I1)=0
         Do  19 I2=1,NABC
            If(CSEQ(I1).EQ.CABC(I2)) then
               ISEQ(I1)=I2
               Go to  20
            End if
 19      Continue
 20   Continue
      If(.NOT.OPTR.AND.IFUN.EQ.3)
     *   Call CFAve(ISEQ,IDMS,BSEQ,LSEQ,CABC,NABC,FAVE)

 21   Continue

* compute cut-off in raw score units

      If(.NOT.OPTR) then
         If(IFUN.EQ.3) then
            Call CPAve(IMPP,IDMP,LPRF,CABC,NABC,PAVE)
            RAVE=0
            Do  I1=0,NABC
               RAVE=RAVE+FAVE(I1)*PAVE(I1)
            End do
         End if
         Call NtoR(XCUT,KCUT,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
      End if

      JSEQ=0

 25   Continue

* should alignments of repeats be computed

      If(OPTM) then
         LCMM=.TRUE.
      Else
         LCMM=.FALSE.
      End if


* compute optimal alignment score

      Call XALI1
     *   (LPRF,LPCI,
     *   KCUT,IIPP,IMPP,IIPX,
     *   BSEQ,LSEQ,ISEQ,
     *   IOPM,IOPI,IOPD,
     *   IOPT,LUNI,
     *   IRC)

      If(OPTA) then
         Continue
      Else if(IOPT.LT.KCUT) then
         go to  50
      End if

* do not compute alignments if no format modifier specified
* and the disjointness is 'UNIQUE'

      If(.NOT.LOUT.AND.LUNI) go to  30

* initialize sequence lock

      Do  I1=1,LSEQ
         LCKS(I1)=.FALSE.
      End do

* find optimal match

      If(LUNI) then
         Call XALIP
     *      (NABC,LPRF,LPCI,NDIP(1),NDIP(2),
     *      IOPT,IIPP,IMPP,CHIP,CHMP,IIPX,
     *      BSEQ,LSEQ,ISEQ,LCKS,
     *      IOPM,IOPI,IOPD,
     *      IOMB,IOM1,IOM2,IOIB,IOI1,IOI2,IODB,IOD1,IOD2,
     *      NALI,IALS,IALB,IAL1,IAL2,IALE,
     *      LUNI,
     *      IRC)
      Else

* find all matches

         Call XALIP
     *      (NABC,LPRF,LPCI,NDIP(1),NDIP(2),
     *      KCUT,IIPP,IMPP,CHIP,CHMP,IIPX,
     *      BSEQ,LSEQ,ISEQ,LCKS,
     *      IOPM,IOPI,IOPD,
     *      IOMB,IOM1,IOM2,IOIB,IOI1,IOI2,IODB,IOD1,IOD2,
     *      NALI,IALS,IALB,IAL1,IAL2,IALE,
     *      LUNI,
     *      IRC)
      End if
      If(IRC.NE.0) then
         Write(NERR,*) '       While processing sequence ',
     *      CSID(1:Lblnk(CSID))
         IRC=0
         Go to 50
      End if

* remove sequence lock if alignments are to be generated

      If(LTRA) then
         Do I1=1,LSEQ
            LCKS(I1)=.FALSE.
         End do
      End if

* OUTPUT

 30   Continue

      If(LUNI) then
         NALI=1
         IALS(1)=IOPT
      End if
      Do  40 I1=1,NALI
         JSEQ=JSEQ+1

         If(LTRA) then
            Call XALIT
     *         (CABC,LPRF,LPCI,NDIP(1),NDIP(2),
     *         IIPP,IMPP,IIPX,
     *         LSEQ,ISEQ,LCKS,
     *         IOPM,IOPI,IOPD,
     *         LALI,CALI,CPMA,IMSC,
     *         IALS(I1),IALB(I1),IALE(I1),
     *         NMAT(I1),
     *         PK2E,PK2B,PK3E,PK3B,
     *         PJ1E,PJ1B,PLAL,
     *         OPTM,
     *         IPMB,IPME,
     *         IRC)
            If(IRC.GT.0) Go to 100
            Do  I2=IAL1(I1),IAL2(I1)
               LCKS(I2)=.FALSE.
            End do
         End if

         If(.NOT.OPTK) then
            Call WPRSM(JSEQ,NMAT(I1),.FALSE.,
     *         LUNI,LOUT,LNOR,LREV,LPFA,OPTZ,OPTL,OPLU,NW,
     *         CPID,CPAC,CPDE,OPTD,OPTR,LDRS,
     *         IALS(I1),IALB(I1),IALE(I1),NALI,IPMB,IPME,
     *         JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *         RNOP,KNPM,MAXN,INOR,IFUN,MNUM,LSEQ,RAVE)
         Else
            Call XPRSM(JSEQ,NMAT(I1),.FALSE.,
     *         LUNI,LOUT,LNOR,LREV,OPTZ,OPTL,OPLU,OPTB,
     *         CSID,CSAC,CSFH,CPID,CPAC,OPTR,OPTF,LDRS,
     *         IALS(I1),IALB(I1),IALE(I1),NALI,IPMB,IPME,
     *         JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *         RNOP,KNPM,MAXN,INOR,IFUN,MNUM,LSEQ,RAVE)
         End if

         If     (OPTS) then
            Call PRSP(CABC,ISEQ,CALI,IALB(I1),IALE(I1),NW,OPTS,OPTX)
         Else if(OPTX) then
            Call PRSP(CABC,ISEQ,CALI,1,LALI,NW,OPTS,OPTX)
         Else if(OPTY) then
            Call PRALI
     *         (LPRF,CHIP,CHMP,IDMP,LSEQ,LREV,
     *         CALI,LALI,IALB(I1),IALE(I1),NW)
         End if

*----------------------------------------------------------------------*
* Compute single repeats for circular profiles
*----------------------------------------------------------------------*

         If(OPTM) then
            If(NMAT(I1).GT.1) then
               Do I2=NMAT(I1),1,-1
                  If(.NOT.OPTK) then
                     Call WPRSM(JSEQ,NMAT(I1)-I2+1,.TRUE.,
     *                  LUNI,LOUT,LNOR,LREV,LPFA,OPTZ,OPTL,OPLU,NW,
     *                  CPID,CPAC,CPDE,OPTD,OPTR,LDRS,
     *                  IMSC(I2),PK2B(I2),PK2E(I2),NALI,
     *                  PK3B(I2),PK3E(I2)-LPRF-1,
     *                  JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *                  RNOP,KNPM,MAXN,INOR,IFUN,MNUM,LSEQ,RAVE)
                  Else
                     Call XPRSM(JSEQ,NMAT(I1)-I2+1,.TRUE.,
     *                  LUNI,LOUT,LNOR,LREV,OPTZ,OPTL,OPLU,OPTB,
     *                  CSID,CSAC,CSFH,CPID,CPAC,OPTR,OPTF,LDRS,
     *                  IMSC(I2),PK2B(I2),PK2E(I2),NALI,
     *                  PK3B(I2),PK3E(I2)-LPRF-1,
     *                  JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *                  RNOP,KNPM,MAXN,INOR,IFUN,MNUM,LSEQ,RAVE)
                  End if

                  If     (OPTS) then
                     Call PRSP(CABC,ISEQ,CALI,PK2B(I2),PK2E(I2),
     *                  NW,OPTS,OPTX)
                  Else if(OPTX) then
                     Call PRXP(CALI,PJ1B(I2),PJ1B(I2)+PLAL(I2)-1,
     *                  PK3B(I2),PK3E(I2),LPRF,NW)
                  Else if(OPTY) then
                     Call PMALI
     *                  (LPRF,CHIP,CHMP,IDMP,LSEQ,LREV,CALI,
     *                  PK2B(I2),PK3E(I2),PK3B(I2),PJ1E(I2),PJ1B(I2),NW)
                  End if

               End do
* End loop over repeats
            End if
         End if
* End loop over matches

 40   Continue

 50   If(OPTB) then

         If(LREV) then

* regenerate original strand

            J1=LSEQ
            Do  I1=1,LSEQ/2
               IS=ISEQ(I1)
               ISEQ(I1)=ISEQ(J1)
               ISEQ(J1)=IS
               J1=J1-1
            End do

            Do  I1=1,LSEQ
               If(ISEQ(I1).NE.0) ISEQ(I1)=NABC-ISEQ(I1)+1
            End do

            LREV=.FALSE.
            Go to   1
         End if
      Else
         Go to   1
      End if

*----------------------------------------------------------------------*
* Complementary strand
*----------------------------------------------------------------------*

* generate complementary sequence

      J1=LSEQ
      Do  I1=1,LSEQ/2
         IS=ISEQ(I1)
         ISEQ(I1)=ISEQ(J1)
         ISEQ(J1)=IS
         J1=J1-1
      End do

      Do  I1=1,LSEQ
         If(ISEQ(I1).NE.0) ISEQ(I1)=NABC-ISEQ(I1)+1
      End do

      LREV=.TRUE.

      Go to  25

 100  Call Exit(IRC)

* errors

 901  Write(NERR,*) 'Error: Normalisation mode',NMOD,' is not'//
     *   ' defined.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Normalisation mode',NMOD,' is not'//
     *   ' defined for level',LCUT,' cut-off.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Normalisation mode(s) of level',LCUT,
     *   ' is not defined in the profile.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100

      End
*----------------------------------------------------------------------*
      Subroutine Repar(
     *   OPTA,OPTB,OPTF,OPTL,OPTR,OPTS,OPTU,OPLU,OPTX,OPTY,OPTZ,OPTM,
     *   OPTK,FPRF,FSEQ,LCUC,NW,NMOD,OPTO,OPTD,OPTV,IRC)

      Logical           OPTA
      Logical           OPTB
      Logical           OPTF
      Logical           OPTL
      Logical           OPLU
      Logical           OPTR
      Logical           OPTS
      Logical           OPTU
      Logical           OPTX
      Logical           OPTY
      Logical           OPTZ
      Logical           OPTM
      Logical           OPTK
      Logical           OPTO
      Logical           OPTD
      Logical           OPTV

      Character*(*)     FPRF
      Character*(*)     FSEQ
      Character*4096   CARG

      IRC=0

      OPTA=.FALSE.
      OPTB=.FALSE.
      OPTF=.FALSE.
      OPTL=.FALSE.
      OPLU=.FALSE.
      OPTR=.FALSE.
      OPTS=.FALSE.
      OPTU=.FALSE.
      OPTX=.FALSE.
      OPTY=.FALSE.
      OPTZ=.FALSE.
      OPTM=.FALSE.
      OPTO=.FALSE.
      OPTK=.FALSE.
      OPTD=.FALSE.
      OPTV=.FALSE.

      LCUC=0
      NW=60
      N1=Iargc()

      K1=0
      I2=1
      Do I1=1,N1
         Call GetArg(I2,CARG)
         If     (CARG(1:1).EQ.'-'.
     *      AND.CARG(2:2).NE.' '.AND.K1.LT.1) then
            If(Index(CARG,'h').NE.0) go to 900
            If(Index(CARG,'a').NE.0) OPTA=.TRUE.
            If(Index(CARG,'b').NE.0) OPTB=.TRUE.
            If(Index(CARG,'d').NE.0) OPTD=.TRUE.
            If(Index(CARG,'f').NE.0) OPTF=.TRUE.
            If(Index(CARG,'l').NE.0) OPTL=.TRUE.
            If(Index(CARG,'L').NE.0) OPLU=.TRUE.
            If(Index(CARG,'r').NE.0) OPTR=.TRUE.
            If(Index(CARG,'s').NE.0) OPTS=.TRUE.
            If(Index(CARG,'u').NE.0) OPTU=.TRUE.
            If(Index(CARG,'x').NE.0) OPTX=.TRUE.
            If(Index(CARG,'y').NE.0) OPTY=.TRUE.
            If(Index(CARG,'z').NE.0) OPTZ=.TRUE.
            If(Index(CARG,'m').NE.0) OPTM=.TRUE.
            If(Index(CARG,'k').NE.0) OPTK=.TRUE.
            If(Index(CARG,'v').NE.0) OPTV=.TRUE.
            If(Index(CARG,'W').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NW
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NW
               End if
            End if
            If(Index(CARG,'C').NE.0) then
               If(CARG(3:3).NE.' ') then
                  CARG(1:2)='  '
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
               End if
               Read(CARG,*,Err=900) LCUC
            End if
            If(Index(CARG,'M').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NMOD
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NMOD
               End if
               OPTO=.TRUE.
            End if
         Else if(K1.LE.1) then
            K1=K1+1
            If     (K1.EQ.1) then
               FSEQ=CARG
            Else if(K1.EQ.2) then
               FPRF=CARG
            End if
         Else

* - cut-off level on command line

            If(CARG(1:2).EQ.'L=') then
               CARG(1:2)=' '
               LCUC=0
               Read(CARG,*,Err=900) LCUC

* - output width on command line

            Else if(CARG(1:2).EQ.'W=') then
               Read(CARG(3:),*,Err=900) NW
            End if
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
      End do

 20   If(K1.NE.2) IRC=1

 100  If(NW.LE.0.OR.NW.GT.512) NW=60
      Return
 900  IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Include          'reprf.f'
      Include          'reseq.f'
      Include          'rfseq.f'
      Include          'xali1.f'
      Include          'xalip.f'
      Include          'RtoN.f'
      Include          'NtoR.f'
      Include          'CFAve.f'
      Include          'CPAve.f'
      Include          'wprsm.f'
      Include          'xprsm.f'
      Include          'xalit.f'
      Include          'lblnk.f'
      Include          'prali.f'
      Include          'pmali.f'
      Include          'prsp.f'
      Include          'prxp.f'
      Include          'Xblnk.f'

