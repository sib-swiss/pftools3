*       Program pfsearch
*----------------------------------------------------------------------*
* $Id: pfsearch.f,v 2.30 2003/12/09 13:42:42 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: Scan a protein or DNA sequence library for profile
*                 matches
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

* sequence

      Character*4096    FSEQ

      Character*64      CSID
      Character*64      CSAC
      Character*64      CSFH
      Character*512     CSDE

      Integer           LSEQ
      Integer           BSEQ

      Integer*2         ISEQ(IDMS)

      Logical           LCKS(IDMS)

      Integer*2         IS

* number of sequences read

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
      Logical           OPLU
      Logical           OPTR
      Logical           OPTS
      Logical           OPTU
      Logical           OPTX
      Logical           OPTY
      Logical           OPTZ
      Logical           OPTM
      Logical           OPTK
      Logical           OPTJ
      Logical           OPTO
      Logical           OPTD
      Logical           OPTV

      Logical           LCMM

      Integer           NCUC
      Integer           KCUC
      real              XCUC
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

C     Character*1024    RCIN
C     Character*1024    RCOUT

* initialization of controlled vocabularies

      Include          'cvini.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0

      LUNI=.FALSE.
      LNOR=.FALSE.
      LREV=.FALSE.
      LTRA=.FALSE.
      LPFA=.FALSE.
      LDRS=.FALSE.

      LEOF=.FALSE.
      CABC(0)='-'

      INBS=0
      INBP=0

* read command line arguments

      Call Repar(
     *   OPTA,OPTB,OPTF,OPTL,OPLU,OPTR,OPTS,OPTU,OPTX,OPTY,OPTZ,OPTM,
     *   OPTK,OPTJ,FPRF,FSEQ,NCUC,KCUC,XCUC,NW,NMOD,OPTO,OPTD,OPTV,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''pfsearch 2.3 revision 5.d'',//
     *      ''Usage: pfsearch [ -abCdfhlLmMkrsuvWxyz ] [ profile-file'',
     *      '' | - ] [ seq-library-file | - ] [ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -a: report optimal alignment for all sequences.'',/
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
     *      ''        cut-off value. An integer value forces -r.'',
     *      '' Same as parameter C.'',/
     *      ''    -M<value>:'',/
     *      ''        set the normalization mode to use for the'',
     *      '' score computation.'',/
     *      ''        Overrides the profile PRIORITY parameter.'',/
     *      )')
         Write(NERR,'(
     *      ''   output modifiers:'',/,
     *      ''    -d: impose length limit on sequence description.'',/
     *      ''    -k: output using the xPSA header (using keyword'',
     *      ''=value pairs).'',/
     *      ''    -j: output using the xPSA header adding the sequence'',
     *      '' matched by itself.'',/
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
     *      ''    [C=cut-off-value]  use option -C instead'',/
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

      If(OPTJ) then
          OPTK=.TRUE.
      End if

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

      If(IRC.GT.0) go to 100

* check parameter consistency

      If(.NOT.OPTV.AND.OPTM.AND..NOT.LPCI) then
         Write(NERR,*) 'Warning: Profile not circular. ',
     *      'Ignoring option -m.'
         OPTM=.FALSE.
      End if
      If(.NOT.OPTV.AND.OPTB.AND.NABC.GT.4) then
         Write(NERR,*) 'Warning: Not a DNA profile. ',
     *      'Ignoring option -b.'
         OPTB=.FALSE.
      End if
      If(.NOT.OPTV.AND.OPTA.AND.NCUC.GT.0) then
         Write(NERR,*) 'Warning: Option -a is set. Ignoring command',
     *      ' line cut-off (option -C).'
         NCUC=0
      End if
      If(.NOT.OPTV.AND.OPTR.AND.NCUC.GT.1) then
         Write(NERR,*) 'Warning: Option -r is set. Please use only',
     *      ' integer command line cut-off (option -C) values.'
         NCUC=0
      End if

* get cut-off for mode nb specified on command line

      KCUT=0
      If(OPTO) then
         J1=0
*   find the specified normalisation mode
         Do 1 I1=1,JNOR
            If(NNOR(I1).EQ.NMOD) J1=I1
 1       Continue
*   exit if mode number not specified in profile
         If(J1.EQ.0) Go to 901

*   if cut-off value was specified on command line do not search levels
         If(NCUC.EQ.2) then
            LNOR=.TRUE.
            INOR=J1
            MNUM=NNOR(J1)
            IFUN=MNOR(J1)
*   exit if an integer cut-off value was specified with option M
         Else if(NCUC.EQ.1) then
            Go to 902
*   find mode number in list of level 0 cut-off values
         Else
            INOR=0
            Do 3 I1=1,JCUT
               If(MCLE(I1).EQ.0) then
                  If(JCNM(I1).NE.0) then
                     J2=0
                     Do 2 I2=1,JCNM(I1)
                        If(MCUT(I2,I1).EQ.NMOD) J2=I2
 2                   Continue
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
 3          Continue
            If(INOR.EQ.0) Go to 901
         End if
      Else

* get cut-off and normalisation modes from profile

         Do   6 I1=1,JCUT
            If(MCLE(I1).EQ.0) then
               INOR=0
               If(JCNM(I1).NE.0) then
                  LNOR=.TRUE.
*               J2=1
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

* cut-off from command line

      If(NCUC.EQ.1) then
         KCUT=KCUC
         OPTR=.TRUE.
      Else if(NCUC.EQ.2.AND.LNOR) then
         XCUT=XCUC
      Else if(.NOT.OPTV.AND.NCUC.EQ.2.AND..NOT.LNOR) then
         Write(NERR,*) 'Warning: Profile does not provide ',
     *      'normalization. Ignoring command line cut-off.'
      End if
      If(.NOT.LNOR) OPTR=.TRUE.

* disjointness definition
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

* average match score for average amino acid composition

      If(.NOT.OPTR.AND.IFUN.EQ.3)
     *   Call CPAve(IMPP,IDMP,LPRF,CABC,NABC,PAVE)

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
* major loop over sequences
*----------------------------------------------------------------------*

* read sequence

 20   Continue
      If(LEOF) go to 100
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

      JSEQ=0

 25   Continue

* should alignments of repeats be computed

      If(OPTM) then
         LCMM=.TRUE.
      Else
         LCMM=.FALSE.
      End if

      BSEQ=1

* compute cut-off in raw score units

      If(.NOT.OPTR) then
         If(IFUN.EQ.3) then
            Call CFAve(ISEQ,IDMS,BSEQ,LSEQ,CABC,NABC,FAVE)
            RAVE=0
            Do  I1=0,NABC
               RAVE=RAVE+FAVE(I1)*PAVE(I1)
            End do
         End if
         Call NtoR(XCUT,KCUT,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
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

      Do I1=1,LSEQ
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
         Write(NERR,*) '       While processing profile ',
     *      CPID(1:Lblnk(CPID))
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

      if(OPTJ .AND. (NALI > 0)) then
            Do I1=1,LSEQ
               CALI(I1)=Char(Ichar(CABC(ISEQ(I1))))
            End do
            Call JPRSM(CSID,CSAC,CSFH,OPTF,1,LSEQ)
            Call PRSP(CABC,ISEQ,CALI,1,LSEQ,NW,OPTS,OPTX)
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
               LCKS(I2)=.TRUE.
            End do
         End if

         If(.NOT.OPTK) then
            Call WPRSM(JSEQ,NMAT(I1),.FALSE.,
     *         LUNI,LOUT,LNOR,LREV,LPFA,OPTZ,OPTL,OPLU,NW,
     *         CSID,CSAC,CSDE,OPTD,OPTR,LDRS,
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
     *                  CSID,CSAC,CSDE,OPTD,OPTR,LDRS,
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


 50   Continue

      If(OPTB) then
         If(LREV) then
            LREV=.FALSE.
            Go to  20
         Else
            Continue
         End if
      Else
         Go to  20
      End if

*----------------------------------------------------------------------*
* Complementary strand
*----------------------------------------------------------------------*

      LREV=.TRUE.

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

      Go to  25

 100  Call Exit(IRC)

* errors

 901  Write(NERR,*) 'Error: Normalisation mode ',NMOD,' is not'//
     *   ' defined.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Cut-off must be a real number when option'//
     *   ' -M is used.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Normalisation mode ',NMOD,' is not'//
     *   'defined for level 0 cut-off.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Normalisation mode(s)',MCUT(I2,I1),
     *   ' of level',MCLE(I1),' is not defined in the profile.'
      Write(NERR,*) '       While processing profile ',
     *   CPID(1:Lblnk(CPID))
      IRC=1
      Go to 100

      End
*----------------------------------------------------------------------*
      Subroutine Repar(
     *   OPTA,OPTB,OPTF,OPTL,OPLU,OPTR,OPTS,OPTU,OPTX,OPTY,OPTZ,OPTM,
     *   OPTK,OPTJ,FPRF,FSEQ,NCUC,KCUC,XCUC,NW,NMOD,OPTO,OPTD,OPTV,IRC)

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
      Logical           OPTJ
      Logical           OPTO
      Logical           OPTD
      Logical           OPTV

      Character*(*)     FPRF
      Character*(*)     FSEQ
      Character*4096    CARG

      IRC=0
      NCUC=0
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
      OPTJ=.FALSE.
      OPTD=.FALSE.
      OPTV=.FALSE.

      NW=60

      N1=Iargc()

      K1=0
      I2=1
      Do  10 I1=1,N1
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
            If(Index(CARG,'j').NE.0) OPTJ=.TRUE.
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
               If(Index(CARG,'.').EQ.0) then
                  NCUC=1
                  Read(CARG,*,Err=900) KCUC
               Else
                  NCUC=2
                  Read(CARG,*,Err=900) XCUC
               End if
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
               FPRF=CARG
            Else if(K1.EQ.2) then
               FSEQ=CARG
            End if
         Else

* - cut-off value on command line

            If     (CARG(1:2).EQ.'C=') then
               CARG(1:2)='  '
               If(Index(CARG,'.').EQ.0) then
                  NCUC=1
                  Read(CARG,*,Err=900) KCUC
               Else
                  NCUC=2
                  Read(CARG,*,Err=900) XCUC
               End if

* - output width on command line

            Else if(CARG(1:2).EQ.'W=') then
               Read(CARG(3:),*,Err=900) NW
            End if
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
 10   Continue

 20   If (K1.NE.2) IRC=1

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
      Include          'jprsm.f'
      Include          'xalit.f'
      Include          'lblnk.f'
      Include          'prali.f'
      Include          'pmali.f'
      Include          'prsp.f'
      Include          'prxp.f'
      Include          'Xblnk.f'

