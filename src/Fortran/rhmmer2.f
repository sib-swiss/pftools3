*----------------------------------------------------------------------*     
* $Id: rhmmer2.f,v 2.7 2003/11/28 11:57:32 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RHMMER2 
     *   (NPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   BLOG,FABC,P0,
     *   CHID,IIPD,CHMD,IMPD, 
     *   INBP,LEOF,IRC)

      Character*(*)     FPRF

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'sterr.f'

* work fields

      Character*512     RCIN
      Logical           LOPN
      Logical           LIDF
      Logical           LEOF
      Logical           LRPF
      Character*8       CHIN(0:26)
      Character*32      CMIS
      Integer           LREA

* initializatons 
      
      IRC=0
      BL=LOG(2.0)
      LEOF=.FALSE.
      LRPF=.FALSE.

* open hmm file 
      
      LIDF=.FALSE.
      CMIS='HMMER2.0'
      Inquire(File=FPRF,OPENED=LOPN)
      If(LOPN) go to   1
      If(FPRF.NE.'-'.OR.NPRF.NE.5)
     *   Open(NPRF,File=FPRF,Status='OLD',Err=900)

* read first line

 1    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:8).NE.'HMMER2.0') go to   1

* name (Id)
      
      CMIS='NAME'
 2    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:4).NE.'NAME') go to   2

      K2=Lblnk(RCIN)
      Do I1=K2,5,-1
         If(RCIN(I1:I1).EQ.' ') then
            K1=I1+1
*** To be modified into goto ?
            go to 3
***
         End if
      End do 

* - lower case to upper, dash to underscore

 3    CPID=' '
      J1=0
      Do I1=K1,K2
         J1=J1+1
         M1=Ichar(RCIN(I1:I1))
         If(M1.LE.122.AND.M1.GE.97) then
            CPID(J1:J1)=Char(M1-32)
         Else
            CPID(J1:J1)=RCIN(I1:I1)
            If(CPID(J1:J1).EQ.'-') CPID(J1:J1)='_'
         End if
      End do 
      LIDF=.TRUE.

* AC, DE

      CPAC=' '
      CPDT=' '
      CPDE=' '
      CMIS='LENG'

 4    Read(NPRF,'(A)',Err=901,End=902) RCIN

      If(RCIN(1:2).EQ.'//') then 
         go to 920
      Else if(RCIN(1:4).EQ.'ACC ') then
         CPAC=RCIN(7:)
         Go to   4
      Else if(RCIN(1:4).EQ.'LENG') then
         Go to   6
      Else if(RCIN(1:4).NE.'DESC') then
         Go to   4
      Else 
         CPDE=RCIN(7:)
      End if

* length

 5    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:4).NE.'LENG') go to   5

 6    Read(RCIN(5:),*,Err=910) LPRF
      If(LPRF.GE.IDMP.OR.LPRF.LE.0) go to 904

* alphabet

      CMIS='ALPH'
 7    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:4).NE.'ALPH') go to   7

      If(Index(RCIN,'Amino').NE.0) then 
         NABC=20
         RCIN='ACDEFGHIKLMNPQRSTVWY'
         Read(RCIN,'(26A)')(CABC(ii1),ii1=1,NABC)
      Else if(Index(RCIN,'Nucleic').NE.0) then 
         NABC=4
         RCIN='ACGT'
         Read(RCIN,'(26A)')(CABC(ii1),ii1=1,NABC)
      Else
         Go to 911
      End if 

* comment lines 

      LHDR=0
      CMIS='XT'
 8    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:2).NE.'XT') then 
         LHDR=LHDR+1
         If(LHDR.LT.512)
     *      CHDR(LHDR)='CC   ' // RCIN
         Go to   8
      End if 
      
* XT (only NB and EC are relevant)

      Read(RCIN(3:),*,Err=910) NXNB,NXNN,NXEC,NXEJ,NXCT,NXCC,NXJB,NXJJ

* Null model (NI is null model residue extension score).

      CMIS='NULT'
 9    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:4).NE.'NULT') Go to 9

      Read(RCIN(5:),*,Err=910) NI
      BLOG=EXP(BL/1000) 
      P0=EXP(Real(NI)/1000.0*BL)
      
      CMIS='NULE'
 10   Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:4).NE.'NULE') Go to  10

      Read(RCIN(5:),*,Err=910)(IIPP(ii1,0),ii1=1,NABC)
      Do I1=1,NABC
         FABC(I1)=(EXP(Real(IIPP(I1,0))/1000*BL))/20.0
      End do 

* Normalization parameters, cut-off values        

      RNOP(1,1)=-LOG10(350.0)
      RNOP(2,1)=1.0
      CMIS='HMM'
 11   Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:3).EQ.'EVD') then 
         CMIS='EVD'
         Read(RCIN(4:),*,Err=910) RNOP(1,1),RNOP(2,1)
      Else if(RCIN(1:3).NE.'HMM') then
         Go to  11 
      End if 

      MDIS=2
      NDIP(1)=1
      NDIP(2)=LPRF

      JNOR=1
      MNOR(1)=1
      NNOR(1)=1
      NNPR(1)=1
      CNTX(1)='bitscore'

      RNOP(1,1)=-RNOP(1,1)*RNOP(2,1)/LOG(10.0)
      RNOP(2,1)= RNOP(2,1)*0.001/LOG(10.0)

* cut-off value

      JCUT=1
      MCLE(1)=0
      CCUT(1)=' '
      JCNM(1)=1
      MCUT(1,1)=1

*** This is done in htop.f:
C     RCUT(1,1)=RC
C     ICUT(1)=NINT((RCUT(1,1)-RNOP(1,1))/RNOP(2,1))+1
***
      
* defaults

      CHID='-'
      Do I1=0,26
         IIPD(I1)=0
      End do   
      
      IIPD(B0)=NLOW
      IIPD(B1)=NLOW
      IIPD(E0)=NLOW
      IIPD(E1)=NLOW
      
      IIPD(BM)=0
      IIPD(BI)=NLOW
      IIPD(BD)=NLOW
      IIPD(BE)=NLOW
      IIPD(MM)=NLOW
      IIPD(MI)=NLOW
      IIPD(MD)=NLOW
      IIPD(ME)=0
      IIPD(IM)=NLOW
      IIPD(II)=NLOW
      IIPD(ID)=NLOW
      IIPD(IE)=NLOW
      IIPD(DM)=NLOW
      IIPD(DI)=NLOW
      
      IIPD(DD)=0
      IIPD(DE)=NLOW

      CHMD='X'
      Do I1=0,26
         IMPD(I1)=0
      End do   
      IMPD(D )=0

* initialize main model

      LPCI=.FALSE.

      Do I1=0,LPRF
         CHIP(I1)='-' 
         Do I2=0,47
            IIPP(I2,I1)=IIPD(I2)
         End do 
      End do 
      Do I2=0,NABC
         IIPP(I2, 0)=NLOW
      End do
      Do I1=1,LPRF
         CHMP(I1)='X'
         Do I2=0,47
            IMPP(I2,I1)=IMPD(I2)
         End do 
      End do 

* main model parameters

      CMIS='HMM'
      Go to  13
 12   Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
 13   If(RCIN(1:3).NE.'HMM') Go to  12 
      Read(NPRF,'(A)',Err=901,End=903) RCIN
      If(RCIN(1:2).EQ.'//') go to 903

* first line (position zero) 

      Read(NPRF,'(A)',Err=901,End=903) RCIN
      If(RCIN(1:2).EQ.'//') go to 903
      Call c2i(RCIN,CHIN,1,3)
      If(CHIN(1).EQ.'*') then
         IIPP(BM, 0)=NLOW
      Else
         Read(CHIN(1),*,Err=912) IIPP(BM, 0)
      End if 
      If(CHIN(2).EQ.'*') then
         IIPP(BI, 0)=NLOW
      Else
         Read(CHIN(2),*,Err=912) IIPP(BI, 0)
      End if 
      If(CHIN(3).EQ.'*') then
         IIPP(BD, 0)=NLOW
      Else
         Read(CHIN(3),*,Err=912) IIPP(BD, 0)
      End if 

* subsequent lines (positions 1 to LPRF) 

      LREA=0
      Do I1=1,LPRF
         Read(NPRF,'(A)',Err=901,End=903) RCIN
         If(RCIN(1:2).EQ.'//') go to 903
         Call c2i(RCIN,CHIN,0,NABC)
         Do I2=1,NABC
            If(CHIN(I2).EQ.'*') then
               IMPP(I2,I1)=NLOW
            Else
               Read(CHIN(I2),*,Err=912) IMPP(I2,I1)
            End if 
         End do
         Read(NPRF,'(A)',Err=901,End=903) RCIN
         Call c2i(RCIN,CHIN,0,NABC)
         Do I2=1,NABC
            If(CHIN(I2).EQ.'*') then
               IIPP(I2,I1)=NLOW
            Else
               Read(CHIN(I2),*,Err=912) IIPP(I2,I1)
            End if 
         End do
         Read(NPRF,'(A)',Err=901,End=903) RCIN
         Call c2i(RCIN,CHIN,0,9)
         Do I2=1,NABC
            If(CHIN(I2).EQ.'*') then
               IIPP(I2,I1+1)=NLOW
            Else
               Read(CHIN(I2),*,Err=912) IIPP(I2,I1+1)
            End if 
         End do
         IIPP(MM,I1)=IIPP( 1,I1+1)
         IIPP(MI,I1)=IIPP( 2,I1+1)
         IIPP(MD,I1)=IIPP( 3,I1+1)
         IIPP(IM,I1)=IIPP( 4,I1+1)
         IIPP(II,I1)=IIPP( 5,I1+1)
         IIPP(DM,I1)=IIPP( 6,I1+1)
         IIPP(DD,I1)=IIPP( 7,I1+1)
         IIPP(B1,I1-1)=IIPP( 8,I1+1)
         IIPP(B0,I1-1)=IIPP( 8,I1+1)
         IIPP(E1,I1)=IIPP( 9,I1+1)
         IIPP(E0,I1)=IIPP( 9,I1+1)
         LREA=LREA+1
      End do     

      IIPP(BM,   0)=IIPP(B0,   0)
      IIPP(B0,   0)=0
      IIPP(B1,   0)=0

      IIPP(BM,LPRF)=IIPP(B0,LPRF)
      IIPP(E0,LPRF)=0
      IIPP(E1,LPRF)=0

* null model corrections, other modifications:

      Do I1=0,LPRF
         If(IIPP(B0,I1).NE.NLOW)
     *      IIPP(B0,I1)=IIPP(B0,I1)+NXNB
         If(IIPP(B1,I1).NE.NLOW)
     *      IIPP(B1,I1)=IIPP(B1,I1)+NXNB
         If(IIPP(E0,I1).NE.NLOW)
     *      IIPP(E0,I1)=IIPP(E0,I1)+NXEC
         If(IIPP(E1,I1).NE.NLOW) 
     *      IIPP(E1,I1)=IIPP(E1,I1)+NXEC
         X=0
         Do I2=1,NABC
            If(IIPP(I2,I1).NE.NLOW) then
               IIPP(I2,I1)=IIPP(I2,I1)-NI
            End if 
            X=X+FABC(I2)*IIPP(I2,I1) 
         End do 
         IIPP(I0,I1)=NINT(X)
         If(IIPP(I0,I1).LE.NLOW+27) IIPP(I0,I1)=NLOW
      End do 
      Do I1=1,LPRF
         X=0
         Do I2=1,NABC
            If(IMPP(I2,I1).NE.NLOW) then
               IMPP(I2,I1)=IMPP(I2,I1)-NI
            End if 
            X=X+FABC(I2)*IMPP(I2,I1) 
         End do 
         IMPP(I0,I1)=NINT(X)
         If(IMPP(I0,I1).LE.NLOW+27) IMPP(I0,I1)=NLOW
      End do 

* add consensus characters: 

      Do I1=1,LPRF
         K1=IMPP( 0,I1) 
         Do I2=1,NABC
            If(IMPP(I2,I1).GE.K1) then
               K1=IMPP(I2,I1)
               CHMP(I1)=CABC(I2)
            End if
         End do
      End do  
      
      INBP=INBP+1
      LIDF=.FALSE.
      LRPF=.TRUE.
      If(LREA.NE.LPRF) Go to 921

 15   Read(NPRF,'(A)',Err=901,End=905) RCIN
      If(RCIN(1:2).NE.'//') go to  15

 100  If(LRPF) IRC=0
      Return

 110  If(LIDF) then
         Write(NERR,*) '       While processing profile ',
     *      CPID(1:Lblnk(CPID))
      End if
      Go to 100

* errors

 900  Write(NERR,*) 'Error: Unable to open profile file'//
     *   ' ''',FPRF(1:Lblnk(FPRF)),'''.'
      IRC=1
      Go to 110
 901  Write(NERR,*) 'Error: Unable to read profile file'//
     *   ' ''',FPRF(1:Lblnk(FPRF)),'''.'
      IRC=1
      Go to 110
 902  LEOF=.TRUE.
      If(INBP.LT.1.OR.CMIS.NE.'HMMER2.0') then
         Write(NERR,*) 'Error: Unexpected end of file. '//
     *      'Unable to find profile ',CMIS(1:Lblnk(CMIS)),' line.'
         Write(NERR,*) '       No profile read.'
         IRC=1
      Else
C         LEOF=.TRUE.
         IRC=-1
      End if
      Go to 110
 903  LEOF=.TRUE.
      If(INBP.LT.1) then
         Write(NERR,*) 'Error: Unexpected end of file. '
         Write(NERR,*) '       Profile length does not match ',
     *      'expected length (',LPRF,').'
         Write(NERR,*) '       No profile read.'
         IRC=1
      Else
C         LEOF=.TRUE.
         IRC=0
      End if
      Go to 110
 904  Write(NERR,*) 'Error: Profile length exceeds buffer size (',
     *   IDMP,').'
      IRC=1
      Go to 110
 905  LEOF=.TRUE.
      Go to 100

 910  Write(NERR,*) 'Error: Unable to read value of parameter(s) '//
     *   CMIS(1:Lblnk(CMIS))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 110
 911  Write(NERR,*) 'Error: Unknown alphabet type defined.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 110
 912  Write(NERR,*) 'Error: Unable to read profile values.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 110

* warnings

 920  Write(NERR,*) 'Warning: Unexpected profile separator (''//'')',
     *   ' found.'
      Write(NERR,*) '         Unable to find ',
     *   CMIS(1:Lblnk(CMIS)),' line.'
      If(LIDF) then
         Write(NERR,*) '         While processing profile ',
     *      CPID(1:Lblnk(CPID))
      End if
      Write(NERR,*) '         Reading next profile.'
      IRC=0
      LIDF=.FALSE.
      CMIS='HMMER2.0'
      Go to 1
 921  Write(NERR,*) 'Warning: Expected (',LPRF,' and ',
     *   'actual (',LREA,') profile sizes do not match.'
      Write(NERR,*) '         Check profile syntax.'
      If(LIDF) then
         Write(NERR,*) '         While processing profile ',
     *      CPID(1:Lblnk(CPID))
      End if
      Go to 15
      
      End
*----------------------------------------------------------------------*
      Subroutine c2i(RCIN,CHIN,N1,N2)
      
      Character*(*)     RCIN
      Character*8       CHIN(0:26)
      Logical           LSTR

      K1=N1-1
      J1=1
      LSTR=.FALSE.
      Do I1=1,Lblnk(RCIN)+1
         If(LSTR) then
            If(RCIN(I1:I1).NE.' ') then
               Continue
            Else
               CHIN(K1)=RCIN(J1:I1-1)
               LSTR=.FALSE.
            End if 
         Else
            If(RCIN(I1:I1).NE.' ') then
               K1=K1+1
*** To be modified into Go to ?
               If(K1.GT.N2) Return
***
               J1=I1
               LSTR=.TRUE.
            Else
               Continue
            End if
         End if 
      End do 
      Return
      End    

