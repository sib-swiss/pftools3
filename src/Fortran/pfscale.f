*       Program pfscale 
*----------------------------------------------------------------------*     
* $Id: pfscale.f,v 2.11 2003/11/28 11:53:33 vflegel Exp $
*----------------------------------------------------------------------*     
*       Function: Fits paramters of an extreme value distribution to a
*                 profile score distribution. Input: sorted score list.
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* DATA
*----------------------------------------------------------------------*

* array sizes, I/O units

      Include          'ardim.f' 

      Parameter        (NOUT=   6)
      Parameter        (NSCL=  10)
      Parameter        (NPRF=  11)    


* function return types

      Integer           Xblnk
      External          Xblnk

* Parameters and options

      Character*512     FSCL  
      Character*512     FPRF  

      Real*8            DL
      Integer           NN
      Integer           IMNB
      Real              RP
      Real              RQ

* Profile fields

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'sterr.f'
      Include          'cvini.f'

* score statistics  

      Real              XSCO(IDMC)
      Real              XFRQ(IDMC)
      Real              XWGT(IDMC)
      Logical           LMNB
      Logical           LLLT
      Logical           LRNM
      Logical           LEOF

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0
      LRNM=.FALSE.
      LEOF=.FALSE.

* read command line 

      Call Repar(FSCL,FPRF,LLLT,DL,NN,RP,RQ,LMNB,IMNB,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''pfscale 2.3 revision 5.d'',//
     *      ''Usage: pfscale [ -lhLMNPQ ] [ score-list | - ] '' 
     *      ''[ profile-file ] [ parameters ]'',/
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -h: print usage help text.'',/
     *      ''    -l: do not impose limit on line length.'',/
     *      ''    -L<value>:'',/
     *      ''        logarithmic base of parameters (default: 10).'',/
     *      ''    -M<value>:'',/
     *      ''        profile mode number to scale.'',/
     *      ''    -N<value>:'',/
     *      ''        database size (default: 14147368).'',/
     *      ''    -P<value>:'',/
     *      ''        upper threshold of probability range (default:'',
     *      '' 0.0001).'',/
     *      ''    -Q<value>:'',/
     *      ''        lower threshold of probability range (default:'',
     *      '' 0.000001).'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',/,
     *      ''    [L=log-base]        use option -L instead'',/,
     *      ''    [M=mode-nb]         use option -M instead'',/,
     *      ''    [N=db-size]         use option -N instead'',/,
     *      ''    [P=upper-threshold] use option -P instead'',/
     *      ''    [Q=lower-threshold] use option -Q instead'',/
     *      )')
         Call Exit(IRC)
      End if
         
      If(FSCL.NE.' '.AND.FSCL.NE.'-') then
         MSCL=NSCL 
         Open(MSCL,File=FSCL,Status='OLD',Err=900)
      Else
         MSCL=5
      End if

      RL=1/LOG(DL)

      RDBS=RL*LOG(Real(NN))
      If(RQ.NE.0) then
         EMAX=-RL*(LOG(RQ))
      Else
         EMAX=100
      End if
      EMIN=-RL*(LOG(RP))

      Do   5 I1=1,IDMC
         Read(MSCL,*,End= 10,Err=901) XSCO(I1)
         XFRQ(I1)=RDBS-RL*LOG(I1-0.5)
         XWGT(I1)=I1*0.5
 5    Continue

 10   NSCO=I1-1       
      If(NSCO.LT.1) Go to 904

      XSM=0
      XFM=0
      XWM=0
      Do  20 I1=1,NSCO
         If(XFRQ(I1).GE.EMIN.AND.XFRQ(I1).LE.EMAX) then 
            XSM=XSM+XWGT(I1)*XSCO(I1)
            XFM=XFM+XWGT(I1)*XFRQ(I1)
            XWM=XWM+XWGT(I1)
         End if
 20   Continue
      XSM=XSM/XWM
      XFM=XFM/XWM


      XSV=0
      XFV=0
      XCO=0
      Do  30 I1=1,NSCO
         If(XFRQ(I1).GE.EMIN.AND.XFRQ(I1).LE.EMAX) then 
            XSV=XSV+XWGT(I1)*( (XSCO(I1)-XSM)**2 ) 
            XFV=XFV+XWGT(I1)*( (XFRQ(I1)-XFM)**2 )
            XCO=XCO+XWGT(I1)*(XFRQ(I1)-XFM)*(XSCO(I1)-XSM)
         End if
 30   Continue
      XSV=(XSV/XWM)**0.5
      XFV=(XFV/XWM)**0.5
      XCO=(XCO/XWM)/(XFV*XSV)

      XB=XCO*XFV/XSV
      XA=XFM-XB*XSM

      If(FPRF.NE.' ') Go to  50 

* Case 1: no profile input file - print list

      Write(NOUT,
     *   '(''# -LogP ='',F8.4,'' + '',F12.8,'' * raw-score'')') XA,XB 
      
      Write(NOUT,'(''#'')')
      Write(NOUT,'(''#   rank raw-score  -logFreq  -logProb'')')
      Write(NOUT,'(''#'')')

      Do  40 I1=1,NSCO
         XPRE=XA+XB*XSCO(I1) 
         Write(NOUT,'(I8,F10.2,F10.4,F10.4)')
     *      I1,XSCO(I1),XFRQ(I1),XPRE
 40   Continue
      Go to 100

* Case 2: Modify profile input file

 50   Continue 

      Call REPRF
     *   (NPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   BLOG,FABC,P0,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,CHIL,IIPL,ILIP,
     *   CHID,IIPD,CHMD,IMPD,
     *   INBP,LEOF,.FALSE.,IRC)

      If(IRC.NE.0) go to 100

* Check normalization modes (Lowest priority mode will be updated)

      If(LMNB) then
         Do I1=1,JNOR
            If(NNOR(I1).EQ.IMNB) then
               J1=I1
               Go to 55
            End if
         End do
         Go to 903
      Else
         J1=1
         K1=NNPR(J1)
         Do I1=2,JNOR
            If(NNPR(I1).LT.K1) then
               K1=NNPR(I1)
               J1=I1
            End if 
         End do 
      End if

 55   CNTX(J1)='-LogE'
      If(MNOR(J1).NE.1) go to 902
 
* add normalisation parameters: 

      RNOP(1,J1)=XA
      RNOP(2,J1)=XB

* define cut-offs:

      Do I1=1,JCUT
         Do I2=1,JCNM(I1)
            If(MCUT(I2,I1).EQ.NNOR(J1)) then
               Call NtoR
     *            (RCUT(I2,I1),ICUT(I1),RNOP,KNPM,MAXN,J1,1,LSEQ,RAVE)
               Go to 60
            End if
         End do
 60      Continue
      End do 

* rescaling command

      LFTR=LFTR+1
      If(LFTR.GT.1024) LFTR=1024
      Do I1=LFTR,2,-1
         CFTR(I1)=CFTR(I1-1)
      End do

      CFTR(1)='CC   /RESCALED_BY="'
      Call Recmd(CFTR(1)(21:130))
      IC=Lblnk(CFTR(1))
      CFTR(1)(IC+1:)='";'

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

 900  Write(NERR,*) 'Error: Unable to open score file'//
     *   ' ''',FSCL(1:Lblnk(FSCL)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to read score value.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Normalization mode number',NNOR(J1),
     *   ' is not linear. Scaling impossible.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Normalization mode number',IMNB,
     *   ' is not defined in profile.'
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Score list does not contain any ',
     *   'valid scores.'
      IRC=1
      Go to 100

      End
*----------------------------------------------------------------------*
      Subroutine Repar(FSCL,FPRF,LLLT,DL,NN,RP,RQ,LMNB,IMNB,IRC) 

      Include          'sterr.f'

      Character*(*)     FSCL  
      Character*(*)     FPRF

      Real*8            DL
      Integer           NN
      Real              RP
      Real              RQ
      Logical           LLLT
      Logical           LMNB
      Integer           IMNB

      Character*512     CARG 

* initializations

      LLLT=.TRUE.
      LLLT=.FALSE.
      IMNB=1
      FSCL=' '
      FPRF=' '
      DL=10.0
      NN=14147368
      RP=0.0001
      RQ=0.000001
      IRC=0

* interpret command line arguments 

      N1=Iargc()

      K1=0
      I2=1
      Do I1=1,N1
         Call GetArg(I2,CARG)
         If     (CARG(1:1).EQ.'-'.AND.CARG(2:2).NE.' ') then
            If(Index(CARG,'l').NE.0) LLLT=.FALSE.
            If(Index(CARG,'h').NE.0) go to 900
            If(Index(CARG,'L').NE.0) then
               If(CARG(3:3).NE.' ') then
                  If(CARG(3:3).EQ.'e'.OR.CARG(3:3).EQ.'E') then
                     DL=EXP(1.0) 
                  Else
                     Read(CARG(3:),*,Err=900) DL
                  End if
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  If(CARG(1:1).EQ.'e'.OR.CARG(1:1).EQ.'E') then
                     DL=EXP(1.0) 
                  Else
                     Read(CARG,*,Err=900) DL
                  End if
               End if
            End if
            If(Index(CARG,'N').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) NN
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) NN
               End if
            End if
            If(Index(CARG,'P').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RP
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RP
               End if
            End if
            If(Index(CARG,'Q').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) RQ
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) RQ
               End if
            End if
            If(Index(CARG,'M').NE.0) then
               If(CARG(3:3).NE.' ') then
                  Read(CARG(3:),*,Err=900) IMNB
                  LMNB=.TRUE.
               Else
                  I2=I2+1
                  Call GetArg(I2,CARG)
                  Read(CARG,*,Err=900) IMNB
                  LMNB=.TRUE.
               End if
            End if

         Else if(CARG(1:3).EQ.'L=e'.OR.CARG(1:3).EQ.'L=E') then
            DL=EXP(1.0) 
         Else if(CARG(1:2).EQ.'L=') then
            Read(CARG(3:),*,Err=900) DL 
         Else if(CARG(1:2).EQ.'N=') then
            Read(CARG(3:),*,Err=900) NN 
         Else if(CARG(1:2).EQ.'P=') then
            Read(CARG(3:),*,Err=900) RP
         Else if(CARG(1:2).EQ.'Q=') then
            Read(CARG(3:),*,Err=900) RQ
         Else if(CARG(1:2).EQ.'M=') then
            Read(CARG(3:),*,Err=900) IMNB
            LMNB=.TRUE.
         Else if(K1.LE.1) then
            K1=K1+1
            If     (K1.EQ.1) then
               FSCL=CARG
            Else if(K1.EQ.2) then
               FPRF=CARG
            End if
C         Else
C            Go to 900
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 20
      End do
      
 20   If(DL.LE.1) Go to 901
      If(NN.LT.1) Go to 902
      If(RP.LT.0.OR.RP.GT.1) Go to 903
      If(RQ.LT.0.OR.RQ.GT.1) Go to 904
      If(LMNB.AND.IMNB.LT.1) Go to 905
 100  Return
 900  IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Value of option -L is out of bound.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Value of option -N is out of bound.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Value of option -P is out of bound.'
      IRC=1
      Go to 100
 904  Write(NERR,*) 'Error: Value of option -Q is out of bound.'
      IRC=1
      Go to 100
 905  Write(NERR,*) 'Warning: Mode numbers should be greater or ',
     *   'equal to 1.'
      Go to 100
      End
*----------------------------------------------------------------------*
      Include          'reprf.f'
      Include          'wrprf.f'
      Include          'recmd.f'
      Include          'NtoR.f'
      Include          'lblnk.f'
      Include          'Xblnk.f'
