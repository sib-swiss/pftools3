*       Program ptof
*----------------------------------------------------------------------*
* $Id: ptof.f,v 2.12 2003/12/01 13:33:04 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: converts a protein profile into a framesearch profile
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
      Parameter        (NNUL=  12)

* profile

      Character*4096    FPRF

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'

      Include          'sterr.f'


* function return types

      Integer           Xblnk
      External          Xblnk

* parameters and options

      Logical           OPTR
      Logical           LLLT

      Real              RB
      Real              RF
      Real              RI
      Real              RX
      Real              RY
      Real              RZ

      Integer           INBP
      Logical           LEOF
      Logical           LRNM

* initialization of controlled vocabularies

      Include          'cvini.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0
      FPRF='stdout'
      LRNM=.FALSE.
      INBP=0

* command line arguments

      Call Repar(FPRF,OPTR,LLLT,RB,RF,RI,RX,RY,RZ,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''ptof 2.3 revision 5.d'',//
     *      ''Usage: ptof [ -hlrBFIXYZ ] [ profile | - ] '',
     *      ''[ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -h: print usage help text.'',/
     *      ''    -l: do not impose limit on line length.'',/
     *      ''    -r: parameters given as normalized score units.''
     *      )')
         Write(NERR,'(
     *      ''    -B<value>:'',/
     *      ''        minimal initiation/termination score (default: '',
     *      ''-50 or -0.5 with option -r).'',/
     *      ''    -F<value>:'',/
     *      ''        frameshift error penalty (default: '',
     *      ''-100 or -1.0 with option -r).'',/
     *      ''    -I<value>:'',/
     *      ''        insert score multiplier (default: 1/3).'',/
     *      ''    -X<value>:'',/
     *      ''        stop codon penalty (default: '',
     *      ''-100 or -1.0 with option -r).'',/
     *      ''    -Y<value>:'',/
     *      ''        intron opening penalty (default: '',
     *      ''-300 or -3.0 with option -r).'',/
     *      ''    -Z<value>:'',/
     *      ''        intron extension penalty (default: '',
     *      ''-1 or -0.01 with option -r).'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',/,
     *      ''    [B=min-beg-end-score]   use option -B instead'',/
     *      ''    [F=deletion-frameshift] use option -F instead'',/
     *      ''    [I=insert-score-mult]   use option -I instead'',/
     *      ''    [X=stop-codon-score]    use option -X instead'',/
     *      ''    [Y=intron-open-penalty] use option -Y instead'',/
     *      ''    [Z=intron-ext-penalty]  use option -Z instead'',/
     *      )')
         Call Exit(IRC)
      End if

* read profile

      If(FPRF.EQ.'-') then
         MPRF=5
      Else
         MPRF=NPRF
      End if

      Call REPRF
     *   (MPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   BLOG,FABC,P0,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,CHIL,IIPL,ILIP,
     *   CHID,IIPD,CHMD,IMPD,
     *   INBP,LEOF,.FALSE.,IRC)

      If(IRC.NE.0) go to 100

* expanded profile size check

      If(LPRF*3.GT.IDMP) go to 901
      If(NABC.LE.4) go to 920

* add command-line to footer lines

 10   Continue
      LFTR=LFTR+1
      If(LFTR.GT.1024) LFTR=1024
      Do I1=LFTR,2,-1
         CFTR(I1)=CFTR(I1-1)
      End do

      CFTR(1)='CC   /GENERATED_BY="'
      Call Recmd(CFTR(1)(21:130))
      IC=Lblnk(CFTR(1))
      CFTR(1)(IC+1:)='";'

*----------------------------------------------------------------------*
* CONVERSION SECTION
*----------------------------------------------------------------------*

* Adjust frame search parameters to normalization function

      If(OPTR) then
         R2=0.01
         K1=1
         J1=NNPR( 1)
         Do I1=2,JNOR
            If(NNPR(I1).LT.J1) then
               K1=I1
               JN=NNPR(I1)
            End if
         End do
         If(MNOR(K1).EQ.1) R2=RNOP(2,K1)

         RB=RB / R2
         RF=RF / R2
         RX=RX / R2
         RY=RY / R2
         RZ=RZ / R2
      End if

* Expand alphabet

      NABC=NABC+1
      CABC(NABC)='O'
      NX=NINT(RX)
      J2=0
      Do I2=1,NABC-1
         J2=J2+IIPD(I2)
      End do
      IIPD(NABC)=J2/(NABC-1)
      J2=0
      Do I2=1,NABC-1
         J2=J2+IIPP(I2,0)
      End do
      IIPP(NABC,0)=J2/(NABC-1)
      Do I1=1,LPRF
         IMPP(NABC,I1)=NX
         J2=0
         Do I2=1,NABC-1
            J2=J2+IIPP(I2,I1)
         End do
         IIPP(NABC,I1)=J2/(NABC-1)
      End do

* - default begin, end score

      IIPD(B0)=MAX(IIPD(B0),NINT(RB))
      IIPD(B1)=MAX(IIPD(B1),NINT(RB))
      IIPD(E0)=MAX(IIPD(E0),NINT(RB))
      IIPD(E1)=MAX(IIPD(E1),NINT(RB))

* - eliminate II transitions

      Do I1=0,LPRF
         If(IIPP(BI,I1).NE.NLOW) IIPP(BI,I1)=IIPP(BI,I1)-IIPP(II,I1)
         If(IIPP(MI,I1).NE.NLOW) IIPP(MI,I1)=IIPP(MI,I1)-IIPP(II,I1)
         If(IIPP(DI,I1).NE.NLOW) IIPP(DI,I1)=IIPP(DI,I1)-IIPP(II,I1)
         Do I2=0,26
            If(IIPP(I2,I1).NE.NLOW)
     *         IIPP(I2,I1)=IIPP(I2,I1)+IIPP(II,I1)
         End do
         IIPP(II,I1)=0
      End do

* Expand profile

* -define architecture of profile according to topology

      If(LPCI) then
         NBB=1
         NB=2
         NE=3*LPRF-1
         NEE=3*LPRF
      Else
         NBB=1
         NB=1
         NE=3*LPRF-2
         NEE=3*LPRF-2
      End if

* - move last insert position to new end of profile

      Do I1=0,46
         IIPP(I1,NEE)=IIPP(I1,LPRF)
      End do
      CHIP(NEE)=CHIP(LPRF)

* - move other positions

      K1=NE
      Do I1=LPRF,1,-1
         Do I2=0,27
            IMPP(I2,K1)=IMPP(I2,I1)
         End do
         CHMP(K1)=CHMP(I1)
         K1=K1-3
      End do

      K1=NE-2
      Do I1=LPRF-1,1,-1
         Do I2= 0,26
            IIPP(I2,K1)=NINT(IIPP(I2,I1)*RI)
         End do
         Do I2=27,46
            IIPP(I2,K1)=IIPP(I2,I1)
         End do
         CHIP(K1)=CHIP(I1)
         K1=K1-3
      End do

* - fill in m-1, m+1 positions

      Do I1=NB+1,NEE,3
         Do I2=0,27
            IMPP(I2,I1)=0
         End do
         CHMP(I1)='>'
      End do

      J1=NB-1
      If(J1.LE.0) J1=J1+3
      Do I1=J1,NEE,3
         Do I2=0,27
            IMPP(I2,I1)=0
         End do
         CHMP(I1)='<'
      End do

*- fill i-1 position

      Do I2= 0,26
         IIPP(I2,NB)=0
      End do
      Do I2=27,46
         IIPP(I2,NB)=NLOW
      End do
      IIPP(MM,NB)=0
      IIPP(DD,NB)=0
      IIPP(MI,NB)=NINT(RF/2)
      IIPP(II,NB)=NINT(RF)
      IIPP(IM,NB)=NINT(RF/2)
      CHIP(NB)='*'
      Do I1=NB+3,NEE-1,3
         Do I2= 0,47
            IIPP(I2,I1)=IIPP(I2,NB)
         End do
         CHIP(I1)='*'
      End do

* - fill i+1 position

      J1=3
      If(LPCI) J1=1
      Do I2= 0,26
         IIPP(I2,J1)=0
      End do
      Do I2=27,46
         IIPP(I2,J1)=NLOW
      End do
      IIPP(MM,J1)=0
      IIPP(DD,J1)=0
      IIPP(MI,J1)=NINT(RY/2)
      IIPP(II,J1)=MIN(-1,NINT(RZ))
      IIPP(IM,J1)=NINT(RY/2)
      CHIP(J1)=':'
      Do I1=J1+3,NEE-1,3
         Do I2=0,46
            IIPP(I2,I1)=IIPP(I2,J1)
         End do
         CHIP(I1)=':'
      End do

* frameshift deletetion scores:

      J1=2
      If(LPCI) J1=1
      K1=3
      If(LPCI) K1=1
      Do I1=J1,LPRF-1
         IIPP(MD,K1)=(NINT(RF-IMPP( D,K1+1))/2)
         IIPP(DM,K1+1)=NINT(RF-IIPP(MD,K1)-IMPP( D,K1+1))
         K1=K1+3
      End do

* begin and end scores

      If(.NOT.LPCI) then
         Do I1=2,NEE,3
            IIPP(B0,I1-1)=IIPP(B0,I1)
            IIPP(B1,I1-1)=IIPP(B1,I1)
            IIPP(BM,I1-1)=IIPP(BM,I1)
            IIPP(BI,I1-1)=IIPP(BI,I1)
            IIPP(BD,I1-1)=IIPP(BD,I1)
            IIPP(E0,I1-1)=IIPP(E0,I1)
            IIPP(E1,I1-1)=IIPP(E1,I1)
            IIPP(ME,I1-1)=IIPP(ME,I1)
            IIPP(IE,I1-1)=IIPP(IE,I1)
            IIPP(DE,I1-1)=IIPP(DE,I1)
            IIPP(B0,I1+1)=IIPP(B0,I1)
            IIPP(B1,I1+1)=IIPP(B1,I1)
            IIPP(BM,I1+1)=IIPP(BM,I1)
            IIPP(BI,I1+1)=IIPP(BI,I1)
            IIPP(BD,I1+1)=IIPP(BD,I1)
            IIPP(E0,I1+1)=IIPP(E0,I1)
            IIPP(E1,I1+1)=IIPP(E1,I1)
            IIPP(ME,I1+1)=IIPP(ME,I1)
            IIPP(IE,I1+1)=IIPP(IE,I1)
            IIPP(DE,I1+1)=IIPP(DE,I1)
         End do
      Else
         Do I1=2,NEE,3
            IIPP(B0,I1-1)=IIPP(B0,I1-2)
            IIPP(B1,I1-1)=IIPP(B1,I1-2)
            IIPP(BM,I1-1)=IIPP(BM,I1-2)
            IIPP(BI,I1-1)=IIPP(BI,I1-2)
            IIPP(BD,I1-1)=IIPP(BD,I1-2)
            IIPP(E0,I1-1)=IIPP(E0,I1-2)
            IIPP(E1,I1-1)=IIPP(E1,I1-2)
            IIPP(ME,I1-1)=IIPP(ME,I1-2)
            IIPP(IE,I1-1)=IIPP(IE,I1-2)
            IIPP(DE,I1-1)=IIPP(DE,I1-2)
            IIPP(B0,I1  )=IIPP(B0,I1+1)
            IIPP(B1,I1  )=IIPP(B1,I1+1)
            IIPP(BM,I1  )=IIPP(BM,I1+1)
            IIPP(BI,I1  )=IIPP(BI,I1+1)
            IIPP(BD,I1  )=IIPP(BD,I1+1)
            IIPP(E0,I1  )=IIPP(E0,I1+1)
            IIPP(E1,I1  )=IIPP(E1,I1+1)
            IIPP(ME,I1  )=IIPP(ME,I1+1)
            IIPP(IE,I1  )=IIPP(IE,I1+1)
            IIPP(DE,I1  )=IIPP(DE,I1+1)
         End do
      End if

      Do I1=NBB,NEE-1
         IIPP(B0,I1)=MAX(IIPD(B0),IIPP(B0,I1))
         IIPP(B1,I1)=MAX(IIPD(B1),IIPP(B1,I1))
         IIPP(E0,I1)=MAX(IIPD(E0),IIPP(E0,I1))
         IIPP(E1,I1)=MAX(IIPD(E1),IIPP(E1,I1))
      End do

* disjointness definition

      If(MDIS.EQ.2) then
         NDIP(1)=NDIP(1)*3-2
         NDIP(2)=NDIP(2)*3-2
      End if

* adjust length

      LPRF=NEE

*----------------------------------------------------------------------*
* OUTPUT SECTION
*----------------------------------------------------------------------*

* write profile

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

* errors

 901  Write(NERR,*) 'Error: Expanded profile length exceeds buffer ',
     *   'size (',IDMP,').'
      IRC=1
      Go to 100

* warnings

 920  Write(NERR,*) 'Warning: Profile seems to be DNA. Check input '//
     *   'file.'
      Go to 10

      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FPRF,OPTR,LLLT,RB,RF,RI,RX,RY,RZ,IRC)

      Character*(*)     FPRF
      Character*4096    CARG

      Logical           OPTR
      Logical           LLLT
      Real              RB
      Real              RF
      Real              RI
      Real              RX
      Real              RY
      Real              RZ

      IRC=0

* initializations

      FPRF='-'

      OPTR=.FALSE.
      LLLT=.TRUE.
      RB=-50
      RF=-100
      RI=1/3.0
      RX=-100
      RY=-300
      RZ=-1

      N1=Iargc()

      K1=0
      I2=1
      Do I1=1,N1
         Call GetArg(I2,CARG)
         If     (CARG(1:1).EQ.'-'.AND.CARG(2:2).NE.' ') then
            If(Index(CARG,'h').NE.0) go to 900
            If(Index(CARG,'l').NE.0) LLLT=.FALSE.
            If(Index(CARG,'r').NE.0) then
               OPTR=.TRUE.
               RB=-0.5
               RF=-1.0
               RX=-1.0
               RY=-3.0
               RZ=-0.01
            End if
            If(Index(CARG,'B').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RB
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RB
               End if
            End if
            If(Index(CARG,'F').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RF
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RF
               End if
            End if
            If(Index(CARG,'I').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RI
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RI
               End if
            End if
            If(Index(CARG,'X').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RX
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RX
               End if
            End if
            If(Index(CARG,'Y').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RY
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RY
               End if
            End if
            If(Index(CARG,'Z').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RZ
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RZ
               End if
            End if
         Else if(CARG(1:2).EQ.'B=') then
            Read(CARG(3:),*,Err=900) RB
         Else if(CARG(1:2).EQ.'F=') then
            Read(CARG(3:),*,Err=900) RF
         Else if(CARG(1:2).EQ.'I=') then
            Read(CARG(3:),*,Err=900) RI
         Else if(CARG(1:2).EQ.'X=') then
            Read(CARG(3:),*,Err=900) RX
         Else if(CARG(1:2).EQ.'Y=') then
            Read(CARG(3:),*,Err=900) RY
         Else if(CARG(1:2).EQ.'Z=') then
            Read(CARG(3:),*,Err=900) RZ
         Else if(K1.LE.0) then
            K1=K1+1
            If     (K1.EQ.1) then
               FPRF=CARG
            End if
         Else
            Go to 900
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
      End do

 20   Continue

 100  Return
 900  IRC=-1
      Go to 100
      End
*----------------------------------------------------------------------*
      Include          'reprf.f'
      Include          'recmd.f'
      Include          'wrprf.f'
      Include          'lblnk.f'
      Include          'Xblnk.f'
