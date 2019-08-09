*       Program ptoh 
*----------------------------------------------------------------------*     
* $Id: ptoh.f,v 2.11 2003/12/01 13:33:05 vflegel Exp $
*----------------------------------------------------------------------*     
*       Function: Reformats profile -> hmm: in-fmt=PROSITE / out-fmt=SAM    
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* DATA
*----------------------------------------------------------------------*     

      Include          'ardim.f' 

* array dimensions and I/O units

      Parameter        (NOUT=   6)    

      Parameter        (NPRF=  11)
      Parameter        (NNUL=  12)

* profile 

      Character*512     FPRF

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'

      Include          'sterr.f'

* null model 

      Character*512     FNUL

      Integer           IABC(26)
      Integer           JABC(26)
      Character*20      DABC

* HMM 

      Include          'hmdat.f' 


* function return types

      Integer           Xblnk
      External          Xblnk

* parameters and options

      Logical           OPTF
      Logical           OPFF
      Logical           OPTH
      Logical           OPTS
      Logical           LEOF

      Real              RD 
      Real              RI
      Real*8            DL
      Integer           INBP

* initialization of controlled vocabularies

      Include          'cvini.f' 

*----------------------------------------------------------------------*     
* INPUT SECTION 
*----------------------------------------------------------------------*     

      IRC=0
      FPRF='stdout'
      INBP=0

* command line arguments
      
      Call Repar(FPRF,FNUL,OPTF,OPFF,OPTH,OPTS,RD,RI,DL,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''ptoh 2.3 revision 5.d'',//
     *      ''Usage: ptoh [ -fhsFDIL ] [ profile-file | - ] '',
     *      ''[ random-model-file ] [ parameters ]'',//
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -f: emulate domain- or semi-global alignment mode.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -s: output in SAM format (if not set: HMMER1 '',
     *      ''format).'',/
     *      ''    -F: emulate local alignment mode.''
     *      )')
         Write(NERR,'(
     *      ''    -D<value>:'',/
     *      ''        delete-to-delete transition probabilities '',
     *      ''(default: 0.9).'',/
     *      ''    -I<value>:'',/
     *      ''        insert-to-insert transition probabilities '',
     *      ''(default: 0.99).'',/
     *      ''    -L<value>:'',/
     *      ''        logarithmic base (default: 1.0233739).'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',/,
     *      ''    [D=del-to-del-prob]  use option -D instead'',/
     *      ''    [I=ins-to-ins-prob]  use option -I instead'',/
     *      ''    [L=log-base]         use option -L instead'',/
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

      If(NABC.GE.20) then 
         DABC='ACDEFGHIKLMNPQRSTVWY'
         MABC=20
      Else
         DABC='ACGT'
         MABC=4
      End if

* read null model

      If     (FNUL.NE.' ') then 
         Call RHNUL(NNUL,FNUL,FABC,MABC,IRC)
         If(IRC.NE.0) go to 100
      Else if(P0.GT.0.0) then
         Call NRNUL(NABC,CABC,FABC,P0,DABC,MABC) 
      Else
         Call DFNUL(FABC,P0,DABC,MABC)
      End if

* reduce / rearrange alphabet

      K1=0
      Do I1=1,NABC
         N1=Index(DABC(1:MABC),CABC(I1)) 
         If(N1.NE.0) then
            IABC(N1)=I1
            K1=K1+1
         End if
      End do 

      If(K1.NE.MABC) go to 901

      Do I1=0,LPRF
         Do I2=1,MABC
            JABC(I2)=IIPP(IABC(I2),I1)
         End do
         Do I2=1,MABC
            IIPP(I2,I1)=JABC(I2)
         End do
         If(I1.NE.0) then  
            Do I2=1,MABC
               JABC(I2)=IMPP(IABC(I2),I1)
            End do
            Do I2=1,MABC
               IMPP(I2,I1)=JABC(I2)
            End do
         End if
      End do 
      NABC=MABC
      Do I1=1,NABC
         CABC(I1)=DABC(I1:I1) 
      End do 

* convert profile into Log(Prob) 

      If     (DL  .NE.0.0) then 
         DL=LOG(DL)
      Else if(BLOG.NE.0.0) then 
         DL=LOG(BLOG)
      Else
         DL=(LOG(2.0)/30.0)
      End if

      Do I1=0,LPRF
         Do I2=0,46 
            If(I1.EQ.0.OR.I2.GT.27) then
               RIHM(I2,I1)=Real(IIPP(I2,I1))*DL
            Else
               RMHM(I2,I1)=Real(IMPP(I2,I1))*DL
               RIHM(I2,I1)=Real(IIPP(I2,I1))*DL
            End if
         End do 
      End do 
      FLOW=Real(NLOW)*DL

* subtract null model

      Do I1=1,NABC
         R1=LOG(FABC(I1))
         RIHM(I1, 0)=RIHM(I1, 0)+R1
         Do I2=1,LPRF
            RIHM(I1,I2)=RIHM(I1,I2)+R1
            RMHM(I1,I2)=RMHM(I1,I2)+R1
         End do
      End do
      
* modify begin state

      RIHM(MM,   0)=RIHM(BM,   0)
      RIHM(MI,   0)=RIHM(BI,   0)
      RIHM(MD,   0)=RIHM(BD,   0)
      RIHM(DM,   0)=FLOW
      RIHM(DI,   0)=FLOW
      RIHM(DD,   0)=FLOW

* modify end state 

      RIHM(MM,LPRF)=RIHM(ME,LPRF)
      RIHM(IM,LPRF)=RIHM(IE,LPRF)
      RIHM(DM,LPRF)=RIHM(DE,LPRF)
      RIHM(MD,   0)=FLOW
      RIHM(ID,   0)=FLOW
      RIHM(DD,   0)=FLOW
      
* scale HMM

      If(OPTS) then 
         Call SCHMM(NOUT,
     *      IDMP,RIHM,RMHM,LPRF,NABC,FLOW,FSCA,OPTF,OPFF,RD,RI)  
      Else 
         Call SCHMM(NERR,
     *      IDMP,RIHM,RMHM,LPRF,NABC,FLOW,FSCA,OPTF,OPFF,RD,RI)  
      End if

* print HMM  

      If(OPTS) then
         If(NABC.EQ.4) then
            R         =RIHM(3, 0) 
            RIHM(3, 0)=RIHM(2, 0) 
            RIHM(2, 0)=R
            Do I1=1,LPRF
               R         =RIHM(3,I1) 
               RIHM(3,I1)=RIHM(2,I1) 
               RIHM(2,I1)=R
               R         =RMHM(3,I1) 
               RMHM(3,I1)=RMHM(2,I1) 
               RMHM(2,I1)=R
            End do
         End if
         Call WRSAM(NOUT,
     *      IDMP,RIHM,RMHM,LPRF,NABC,FABC,FLOW,FSCA,DL) 
      Else
         Call WRHMR(NOUT,NERR,
     *      IDMP,RIHM,RMHM,LPRF,NABC,CABC,FLOW,FSCA,DL) 
      End if

 100  Call Exit(IRC)

* errors
      
 901  Write(NERR,*) 'Error: Incompatible alphabets between profile and',
     *   ' null model.'
      IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FPRF,FNUL,OPTF,OPFF,OPTH,OPTS,RD,RI,DL,IRC)

      Character*(*)     FPRF
      Character*(*)     FNUL
      Character*512     CARG 

      Logical           OPTF
      Logical           OPFF
      Logical           OPTH
      Logical           OPTS

      Real*8            DL

      IRC=0

* initializations

      FPRF='-'
      FNUL=' '

      OPTF=.FALSE.
      OPFF=.FALSE.
      OPTH=.TRUE.
      OPTS=.FALSE.

      DL=1.0233739
      RI=0.99
      RD=0.9
      
      N1=Iargc()

      K1=0
      I2=1
      Do I1=1,N1
         Call GetArg(I2,CARG)
         If     (CARG(1:1).EQ.'-'.AND.CARG(2:2).NE.' ') then
            If(Index(CARG,'h').NE.0) go to 900
            If(Index(CARG,'f').NE.0) OPTF=.TRUE.
            If(Index(CARG,'F').NE.0) OPFF=.TRUE.
            If(Index(CARG,'s').NE.0) then
               OPTS=.TRUE.
               OPTH=.FALSE.
            End if
            If(Index(CARG,'D').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RD
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RD
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
            If(Index(CARG,'L').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) DL
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) DL
               End if
            End if
         Else if(CARG(1:2).EQ.'D=') then
            Read(CARG(3:),*,Err=900) RD 
         Else if(CARG(1:2).EQ.'I=') then
            Read(CARG(3:),*,Err=900) RI
         Else if(CARG(1:2).EQ.'L=') then
            Read(CARG(3:),*,Err=900) DL
         Else if(K1.LE.1) then
            K1=K1+1
            If     (K1.EQ.1) then
               FPRF=CARG
            Else if(K1.EQ.2) then
               FNUL=CARG
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
      Subroutine NRNUL(NABC,CABC,FABC,P0,DABC,MABC) 

      Character         CABC(0:26)
      Real              FABC(0:26)
      Character*(*)     DABC
      Real              RABC(20)

      X=0
      Do I1=1,NABC
         N1=Index(DABC(1:MABC),CABC(I1))
         If(N1.NE.0) then
            RABC(N1)=FABC(I1)
            X=X+RABC(N1)
         End if
      End do

      X=P0/X 

      Do I1=1,MABC
         FABC(I1)=RABC(I1)*X
      End do

      Return
      End
*----------------------------------------------------------------------*
      Include          'reprf.f'
      Include          'rhnul.f'
      Include          'wrsam.f' 
      Include          'wrhmr.f' 
      Include          'schmm.f' 
      Include          'dfnul.f'
      Include          'lblnk.f' 
      Include          'Xblnk.f'
