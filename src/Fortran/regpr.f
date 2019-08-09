*----------------------------------------------------------------------*     
* $Id: regpr.f,v 2.5 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine REGPR
     *   (NGPR,FGPR,
     *   RG,RE,RF,RO,LSYM,
     *   CPID,CPAC,CPDE,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP, 
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT, 
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   CHID,IIPD,CHMD,IMPD,
     *   IRC)

      Character*(*)     FGPR
      Character*512     RCIN  
      Character         B

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'

      Include          'sterr.f'
      
      Integer           IPRF(32)
      Character         CPRF

      Logical           LSYM

      Integer           Getc

      IRC=0

* open input file
      
      If(FGPR.EQ.'-') then
 1       Open(NGPR,Status='SCRATCH',Err=901)
 2       Continue 
         Do I1=1,512
            If(Getc(B).NE.0) go to 3
            If(Ichar(B).EQ.10) then
               Write(NGPR,'(512A)',Err=903)(RCIN(ii1:ii1),ii1=1,I1-1)
               Go to   2  
            Else
               RCIN(I1:I1)=B
            End if
         End do 
         Go to 905
 3       Rewind(NGPR)
      Else
         Open(NGPR,File=FGPR,Status='OLD',Err=900)
      End if
      
* initialize 

* - profile header

      CPID='GCG_PROFILE'
      CPAC='GC99999'
      CPDE='Automatically reformatted from file ''' 
     *   // FGPR(1:Lblnk(FGPR))
     *   // '''.' 

* - accessories

      LPCI=.FALSE.

      JNOR=1
      MNOR(1)=1
      NNOR(1)=1
      NNPR(1)=1   
      CNTX(1)='OrigScore'
      RNOP(1,1)=0.0
      RNOP(2,1)=1/RF

      JCUT=1
      MCLE(1)=0
      CCUT(1)=' '
      ICUT(1)=0
      JCNM(1)=1
      RCUT(1,1)=0.0
      MCUT(1,1)=1 

* - defaults for match and insert position  

      CHID='-'
      Do  15 I1=1,26 
         IIPD(I1)=0
 15   Continue
      
      IIPD(B0)=0
      IIPD(B1)=NLOW
      IIPD(E0)=0
      IIPD(E1)=NLOW

      IIPD(BM)=0
      IIPD(BI)=NLOW
      IIPD(BD)=NLOW
      IIPD(BE)=NLOW
      IIPD(MM)=0
      IIPD(MI)=NLOW
      IIPD(MD)=NLOW
      IIPD(ME)=0
      IIPD(IM)=0
      IIPD(II)=0
      IIPD(ID)=NLOW
      IIPD(IE)=NLOW
      IIPD(DM)=0
      IIPD(DI)=NLOW
      IIPD(DD)=0
      IIPD(DE)=NLOW

      IIPD(I0)=0

      CHMD='X'
      Do  16 I1=1,26 
         IMPD(I1)=0
 16   Continue

      IIPD(M0)=0 
      IMPD(D )=0

      Do  18 I1=0,27
         IMPP(I1,0)=NLOW
 18   Continue 

* read alphabet

 25   Read(NGPR,'(A)',End=905,Err=902) RCIN
      If(RCIN( 1: 4).NE.'Cons') go to  25
      
      IC1=Index(RCIN,'Gap')
      K1=0
      Do  29 I1=5,IC1-1
         If(RCIN(I1:I1).NE.' ') then
            K1=K1+1
            CABC(K1)=RCIN(I1:I1)
         End if
 29   Continue
      NABC=K1

* read numbers

      K1=0
 30   Read(NGPR,'(A)',End= 50,Err=902) RCIN
      If(RCIN( 1: 1).EQ.'!') go to  30

* - input line

C      RCIN(1024:1024)='@'
      CPRF=RCIN(2:2)
      Read(RCIN(3:512),*,Err=910,End= 50)
     *   (IPRF(ii1),ii1=1,NABC+2)
      Do  34 I2=1,NABC
         IPRF(I2)=NINT(Real(IPRF(I2))/100*RF+RO)
 34   Continue
      If(LSYM) then 
         NGO=-NINT(Real(IPRF(NABC+1))/200*RF*RG)
      Else
         NGO=-NINT(Real(IPRF(NABC+1))/100*RF*RG)
      End if 
      NGE=-NINT(Real(IPRF(NABC+2))/100*RF*RE)

* - build insert position 

      CHIP(K1)=CHID 
      Do  36 I1=0,46
         IIPP(I1,K1)=IIPD(I1)
 36   Continue
      Do  37 I1=1,NABC
         IIPP(I1,K1)=NGE
 37   Continue
      IIPP(MI,K1)=NGO
      IIPP(MD,K1)=NGO
      If(LSYM) then
         IIPP(IM,K1)=NGO
         IIPP(DM,K1)=NGO
      End if 

* - build match position

      K1=K1+1
      If(K1.GT.IDMP) go to 915
      CHMP(K1)=CPRF
      Do  43 I1=1,NABC
         IMPP(I1,K1)=IPRF(I1)
 43   Continue 
      IMPP( 0,K1)=0
      IMPP(D ,K1)=NGE

      Go to  30

 50   LPRF=K1
      If(LPRF.LE.0) go to 920

* - disjointness definition

      MDIS=2
      NDIP(1)=1+LPRF/10
      NDIP(2)=LPRF-LPRF/10

* - defaults for gap weights

      NGO=IIPP(MI,0) 
      NGE=IMPP( D,1) 
      Do  53 I1=1,LPRF-1
         NGO=MIN(NGO,IIPP(MI,I1))
         NGE=MIN(NGE,IMPP( D,I1+1))
 53   Continue

      IIPD(MI)=NGO
      IIPD(MD)=NGO
      If(LSYM) then 
         IIPD(DM)=NGO
         IIPD(IM)=NGO
      End if
      
      Do  54 I1=1,NABC
         IIPD(I1)=NGE
 54   Continue 
      IMPD( D)=NGE

* - last insert position 

      CHIP(K1)=CHID 
      Do  60 I1=0,46
         IIPP(I1,K1)=IIPD(I1)
 60   Continue

* - domain global mode:

      IIPP(B1,   0)=0
      IIPP(E1,LPRF)=0

* - move DM scores one position forward: 

      If(LSYM) then
         IIPP(DM, 0)=IIPD(DM)
         Do  65 I1=LPRF, 1,-1
            IIPP(DM,I1)=IIPP(DM,I1-1)
 65      Continue
      End if  

 100  Return 

* errors
 
 900  Write(NERR,*) 'Error: Unable to open profile file'//
     *   ' ''',FGPR(1:Lblnk(FGPR)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to create temporary file.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Unable to read profile file.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Unable to write to temporary file.'
      IRC=1
      Go to 100
 905  Write(NERR,*) 'Error: Unexpected end of file. '//
     *   'Unable to find profile alphabet.'
      IRC=1
      Go to 100
 910  Write(NERR,*) 'Error: Unable to read profile values.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
 915  Write(NERR,*) 'Error: Profile length exceeds buffer size (',
     *   IDMP,').'
      IRC=1
      Go to 100
 920  Write(NERR,*) 'Error: Unexpected end of profile. Profile has '//
     *   'zero length. Check profile syntax.'
      IRC=1
      Go to 100

      End
