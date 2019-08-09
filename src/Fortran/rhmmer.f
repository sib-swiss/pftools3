*----------------------------------------------------------------------*     
* $Id: rhmmer.f,v 2.5 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RHMMER(NPRF,FPRF,LPRF,RIHM,RMHM,IDMP,NABC,CABC,
     *   INBP,IRC)

      Include          'pfind.f' 
      Include          'sterr.f'

      Character*(*)     FPRF

      Include          'hmdat.f'

      Character         CABC(0:26)

      Character*512     RCIN
      Character*32      CMIS

* initializatons 
      
      IRC=0

* open HMM file

      If(FPRF.NE.'-') then
         Open(NPRF,File=FPRF,Status='OLD',Err=900)
      End if

* read header lines

 1    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:8).EQ.'HMMER2.0') Go to 910
      If(RCIN(1:1).EQ.'#') Go to   1
      CMIS='profile length'
      Read(RCIN,*,Err=903,End=903) LPRF 
      If(LPRF.GT.IDMP) Go to 904

 2    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   2
      CMIS='alphabet length'
      Read(RCIN,*,Err=903,End=903) NABC
      If(NABC.GT.26) Go to 904

 3    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   3

 4    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   4
      CMIS='alphabet'
      Read(RCIN,'(26A)',Err=903,End=903)
     *   (CABC(ii1),ii1=1,NABC)

 5    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   5

 6    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   6

* read model parameters

      CMIS=' '
      Do I1=0,LPRF

* - match state

         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(MM,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(MD,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(MI,I1)

         Do I2=1,NABC
            If(next(NPRF,RCIN).NE.0) Go to 902
            Read(RCIN,*,Err=903) RMHM(I2,I1)
         End do 

* - delete state

         If(next(NPRF,RCIN).NE.0) Go to 902
         Read(RCIN,*,Err=903) RIHM(DM,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(DD,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(DI,I1)

         RMHM( D,I1)=1.0

* - insert state

         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(IM,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(ID,I1)
         If(next(NPRF,RCIN).NE.0) Go to 902 
         Read(RCIN,*,Err=903) RIHM(II,I1)

         Do I2=1,NABC
            If(next(NPRF,RCIN).NE.0) Go to 902 
            Read(RCIN,*,Err=903) RIHM(I2,I1)
         End do 

      End do 

 100  Return

 900  Write(NERR,*) 'Error: Unable to open profile file'//
     *   ' ''',FPRF(1:Lblnk(FPRF)),'''.'
      IRC=1
      Go to 100

 901  Write(NERR,*) 'Error: Unable to read profile file'//
     *   ' ''',FPRF(1:Lblnk(FPRF)),'''.'
      IRC=1
      Go to 100

 902  Write(NERR,*) 'Error: Unexpected end of file.'
      Write(NERR,*) '       No profile read.'
      IRC=1
      Go to 100

 903  Write(NERR,*) 'Error: Unable to read value of parameter(s) '//
     *   CMIS(1:Lblnk(CMIS))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100

 904  Write(NERR,*) 'Error: Value of parameter '''//
     *   CMIS(1:Lblnk(CMIS)),''' exceeds buffer size.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100

 910  Write(NERR,*) 'Error: HMM profile file in wrong format.'//
     *   ' Try without option -o.'
      IRC=1
      Go to 100

      End
*----------------------------------------------------------------------*
      Integer function next(NPRF,RCIN)

      Character*(*)     RCIN       

 1    Read(NPRF,'(A)',Iostat=NEXT) RCIN
      If(RCIN(1:1).EQ.'#') go to   1
      Return

      End
