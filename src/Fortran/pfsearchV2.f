*       Program pfsearch
*----------------------------------------------------------------------*     
*       Function: Scan a protein or DNA sequence library for profile 
*                 matches 
*       Author:   Philipp Bucher
*       Version:  This file is part of pftools release 2.2 June 1999
*----------------------------------------------------------------------*     
* DATA
*----------------------------------------------------------------------*     

* array dimensions and I/O units

        Include         'ardim.f' 

        Parameter        (NOUT=   6)    

        Parameter        (NPRF=  11)    
        Parameter        (NSEQ=  12)    

* profile 

        Character*64      FPRF

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

* sequence

        Character*64      FSEQ

        Character*64      CSID
        Character*64      CSAC
        Character*256     CSDE

        Integer           LSEQ
        Integer*2         ISEQ(IDMS)

        Logical           LCKS(IDMS)

        Integer*2         IS

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
        
        Integer           NCUC
        Integer           KCUC
        real              XCUC

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

        Character*256     RCIN

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
        LEOF=.FALSE.

        LEOF=.FALSE.
        CABC(0)='-'
 
* read command line arguments

        Call Repar(
     *     OPTA,OPTB,OPTF,OPTL,OPLU,OPTR,OPTS,OPTU,OPTX,OPTY,OPTZ,
     *     FPRF,FSEQ,NCUC,KCUC,XCUC,NW,IRC)
        If(IRC.NE.0) then 
           Write(NERR,'(
     *      ''Usage: pfsearch [ -abflLrsuxyz ] [ profile-file | - ] '',
     *      ''[ seq-library-file | - ] [ parameters ]'',//,
     *      ''   valid parameters are:'',//,
     *      ''                 [C=cut-off-value]          '',/
     *      ''                 [W=output-width]           '',/
     *        )')
           Stop
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

* read profile

        Call REPRF
     *    (MPRF,FPRF,
     *     CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *     BLOG,FABC,P0,
     *     CDIS,JDIP,MDIS,NDIP,
     *     CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP, 
     *     JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT, 
     *     IDMP,CHIP,IIPP,CHMP,IMPP,
     *     CHID,IIPD,CHMD,IMPD,
     *     IRC)

        If(IRC.GT.0) go to 100

* cut-off form profile 

           KCUT=0
        Do   6 I1=1,JCUT
           If(MCLE(I1).EQ.0) then
                 INOR=0
              If(JCNM(I1).NE.0) then 
                 LNOR=.TRUE.
                       J2=1
                 Do  5 I2=1,JCNM(I1)
                          J3=0
                    Do  4 I3=1,JNOR
                       If(MCUT(I2,I1).EQ.NNOR(I3)) J3=I3
    4               Continue
                 
                    If     (I2.EQ.1) then 
                       INOR=J3
                       NPRI=NNPR(J3)    
                       IFUN=MNOR(J3)
                       KCUT=ICUT(I1)
                       XCUT=RCUT(I2,I1)
                    Else if(NNPR(J3).LT.NPRI) then
                       INOR=J3
                       NPRI=NNPR(J3)    
                       IFUN=MNOR(J3)
                       KCUT=ICUT(I1)
                       XCUT=RCUT(I2,I1)
                    End if
    5            Continue
              End if 

              If(JCNM(I1).EQ.0.OR.INOR.EQ.0) then
                 KCUT=ICUT(I1)
                 LNOR=.FALSE.
              End if
           End if
    6   Continue

* cut-off from command line 

        If     (NCUC.EQ.1) then 
           KCUT=KCUC 
           LNOR=.FALSE.
        Else if(NCUC.EQ.2.AND.LNOR) then
           XCUT=XCUC 
        End if

        If(OPTR) LNOR=.FALSE.

* disjointness definition

        If(MDIS.EQ.1.OR.OPTU.OR.OPTA) then
           LUNI=.TRUE.
        Else
           LUNI=.FALSE.

* - initialize profile lock

           If(.NOT.LPCI) then
              Do  8 I1=0,NDIP(1)-1 
                 IIPP(E0,I1)=NLOW
                 IIPP(E1,I1)=NLOW
    8         Continue

              Do  9 I1=NDIP(2),LPRF 
                 IIPP(B0,I1)=NLOW
                 IIPP(B1,I1)=NLOW
    9         Continue
           End if

        End if

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
	   print *, "DBG", i1, "B1", IIPP( B1,I1), "BM", IIPP( BM,I1)
     *             , IIPX( XM,I1)
   10   Continue

* average match score for average amino acid composition 

        If(LNOR.AND.IFUN.EQ.3) 
     *     Call CPAve(IMPP,IDMP,LPRF,CABC,NABC,PAVE)

* alignment and ouptut format switches 

        If(OPTX.OR.OPTY.OR.OPTZ) then 
           LTRA=.TRUE.
        Else
           LTRA=.FALSE.
        End if
        If(OPTS.OR.OPTX) then
           LPFA=.TRUE.
        Else
           LPFA=.FALSE.
        End if
  
*----------------------------------------------------------------------*
* major loop over sequences
*----------------------------------------------------------------------*

* read sequence  

   20   Continue
        If(LEOF) go to 100
        If(OPTF) then 
           Call RFSEQ
     *      (MSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,LSEQ,ISEQ,LEOF,RCIN,IRC)
        Else 
           Call RESEQ
     *      (MSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,LSEQ,ISEQ,LEOF,RCIN,IRC)
        End if 
        If(IRC.EQ.-1) go to 100

        If(IRC.NE.0) then 
           Write(NERR,'(
     *      ''Sequence file unreadable or in wrong format.''
     *        )')
           Stop
        End if 

        JSEQ=0

   25   Continue

* compute cut-off in raw score units

        If(LNOR) then
           If(IFUN.EQ.3) then
              Call CFAve(ISEQ,IDMS,LSEQ,CABC,NABC,FAVE)
                 RAVE=0
              Do  I1=0,NABC
                 RAVE=RAVE+FAVE(I1)*PAVE(I1)
              End do
           End if
           Call NtoR(XCUT,KCUT,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
        End if

* compute optimal alignment score

        Call XALI1
     *    (NABC,CABC,LPRF,LPCI,
     *     KCUT,IDMP,IIPP,IMPP,CHIP,CHMP,IIPX,
     *     IDMS,LSEQ,ISEQ,
     *     IOPM,IOPI,IOPD,
     *     IOPT,LUNI,  
     *     IRC)

        If(OPTA) then
           Continue         
        Else if(IOPT.LT.KCUT) then
           go to  50
        End if 

        If(LUNI) go to  30 

* initialize sequence lock

        Do I1=1,LSEQ
           LCKS(I1)=.FALSE.
        End do  

* find multiple matches 

        Call XALIP
     *    (NABC,CABC,LPRF,LPCI,NDIP(1),NDIP(2),
     *     KCUT,IDMP,IIPP,IMPP,CHIP,CHMP,IIPX,
     *     IDMS,LSEQ,ISEQ,LCKS,
     *     IOPM,IOPI,IOPD,
     *     IOMB,IOM1,IOM2,IOIB,IOI1,IOI2,IODB,IOD1,IOD2,
     *     IDMN,NALI,IALS,IALB,IAL1,IAL2,IALE, 
     *     IRC)

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
     *          (NABC,CABC,LPRF,LPCI,NDIP(1),NDIP(2),
     *           KCUT,IDMP,IIPP,IMPP,CHIP,CHMP,IIPX,
     *           IDMS,LSEQ,ISEQ,LCKS,
     *           IOPM,IOPI,IOPD,
     *           LALI,IDMA,CALI,IDMM,CPMA,
     *           IALS(I1),IALB(I1),IALE(I1), 
     *           IPMB,IPME,
     *           IRC)
              Do  I2=IAL1(I1),IAL2(I1)
                 LCKS(I2)=.TRUE.
              End do
           End if 

           Call WPRSM(JSEQ,
     *       LUNI,LNOR,LREV,LPFA,OPTZ,OPTL,OPLU,NW,
     *       CSID,CSAC,CSDE,
     *       IALS(I1),IALB(I1),IALE(I1),NALI,IPMB,IPME,
     *       JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *       RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)

           If     (OPTS) then
              Write(6,'((60A))')(CABC(ISEQ(ii1)),ii1=IALB(I1),IALE(I1))
           Else if(OPTX) then
              Write(6,'((60A))')(CALI(ii1),ii1=1,LALI)
           Else if(OPTY) then 
              Call PRALI
     *          (LPRF,CHIP,CHMP,IDMP,LSEQ,LREV,
     *           CALI,LALI,IALB(I1),IALE(I1))
           End if 

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

  100   Stop
        End
*----------------------------------------------------------------------*     
        Subroutine Repar(
     *     OPTA,OPTB,OPTF,OPTL,OPLU,OPTR,OPTS,OPTU,OPTX,OPTY,OPTZ,
     *     FPRF,FSEQ,NCUC,KCUC,XCUC,NW,IRC)

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

        Character*64      FPRF
        Character*64      FSEQ
        Character*64      CARG

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

        NW=132

        N1=Iargc()

           K1=0
        Do  10 I1=1,N1
           Call GetArg(I1,CARG)
           If     (CARG(1:1).EQ.'-'.
     *                AND.CARG(2:2).NE.' '.AND.K1.LT.1) then
              If(Index(CARG,'a').NE.0) OPTA=.TRUE.
              If(Index(CARG,'b').NE.0) OPTB=.TRUE.
              If(Index(CARG,'f').NE.0) OPTF=.TRUE.
              If(Index(CARG,'l').NE.0) OPTL=.TRUE.
              If(Index(CARG,'L').NE.0) OPLU=.TRUE.
              If(Index(CARG,'r').NE.0) OPTR=.TRUE.
              If(Index(CARG,'s').NE.0) OPTS=.TRUE.
              If(Index(CARG,'u').NE.0) OPTU=.TRUE.
              If(Index(CARG,'x').NE.0) OPTX=.TRUE.
              If(Index(CARG,'y').NE.0) OPTY=.TRUE.
              If(Index(CARG,'z').NE.0) OPTZ=.TRUE.
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
                    Read(CARG,*) KCUC
                 Else 
                    NCUC=2
                    Read(CARG,*) XCUC
                 End if 
              Else if(CARG(1:2).EQ.'W=') then
                 Read(CARG(3:64),*,Err=900) NW
              End if
           End if
   10   Continue 

        If (K1.NE.2) IRC=1 
           
  100   Return
  900   IRC=1
        Go to 100
        End 
*----------------------------------------------------------------------*     
        Include          'reprf.f'
        Include          'reseq.f'
        Include          'rfseq.f'
#include "xali1.f"
        Include          'xalip.f'
        Include          'RtoN.f'
        Include          'NtoR.f'
        Include          'CFAve.f'
        Include          'CPAve.f'
        Include          'wprsm.f'
        Include          'xalit.f'
        Include          'lblnk.f'
        Include          'prali.f'
