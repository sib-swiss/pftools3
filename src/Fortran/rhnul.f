*----------------------------------------------------------------------*     
* $Id: rhnul.f,v 2.5 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RHNUL(NNUL,FNUL,FABC,NABC,IRC)   

      Include          'sterr.f'

      Character*(*)     FNUL 

      Real              FABC(0:26)

      Character*512     RCIN
      Character*512     CBUF     
      Integer           IOST

      IOST=0

      Open(NNUL,File=FNUL,Status='OLD',Err=900)

 1    Read(NNUL,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:1).EQ.'#') Go to   1

      If(RCIN(1:5).EQ.'Amino'.AND.NABC.LT.20) then
         Go to 903
      Else if(RCIN(1:5).EQ.'Amino'.AND.NABC.GE.20) then
         Continue
      Else if(RCIN(1:7).EQ.'Nucleic'.AND.NABC.NE.4) then
         Go to 903
      Else if(RCIN(1:7).EQ.'Nucleic'.AND.NABC.EQ.4) then
         Continue
      Else 
         Go to 1
      End if 

      M=0
 2    RCIN=' '
      Read(NNUL,'(A)',Err=901,Iostat=IOST) RCIN
      If(IOST.NE.0.AND.RCIN.EQ.' ') Go to 5
      If(RCIN(1:1).EQ.'#'.OR.Lblnk(RCIN).EQ.0) Go to   2
      
      L=Index(RCIN,'#')
      If(L.EQ.0) then 
         L=Lblnk(RCIN)
      Else
         L=Lblnk(RCIN(1:L-1))
      End if 
      CBUF(M+1:)=RCIN(1:L)
      M=M+L+1
      If(M.GT.512) Go to 907
      Go to 2
      
 5    Read(CBUF(1:M+1),*,End=905,Err=906)(FABC(ii1),ii1=1,NABC)
 
 100  Return

 900  Write(NERR,*) 'Error: Unable to open random model file'//
     *   ' ''',FNUL(1:Lblnk(FNUL)),'''.'
      IRC=1
      Go to 100
 901  Write(NERR,*) 'Error: Unable to read random model file'//
     *   ' ''',FNUL(1:Lblnk(FNUL)),'''.'
      IRC=1
      Go to 100
 902  Write(NERR,*) 'Error: Unexpected end of file. '//
     *      'Unable to read random model.'
      IRC=1
      Go to 100
 903  Write(NERR,*) 'Error: Random model alphabet does not match '//
     *      'profile alphabet.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Lblnk(RCIN))
      IRC=1
      Go to 100
 905  Write(NERR,*) 'Error: Not enough values in random model.'
      Write(NERR,*) '       values read: ',
     *   CBUF(1:Lblnk(CBUF))
      IRC=1
      Go to 100
 906  Write(NERR,*) 'Error: Unable to read random model values.'
      IRC=1
      Go to 100
 907  Write(NERR,*) 'Error: Random model exceeds buffer size.'
      Write(NERR,*) '       at line: ',
     *   CBUF(1:Lblnk(CBUF))
      IRC=1
      Go to 100

      End
