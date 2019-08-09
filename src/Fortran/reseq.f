*----------------------------------------------------------------------*     
* $Id: reseq.f,v 2.12 2003/11/18 10:50:07 vflegel Exp $
*----------------------------------------------------------------------*     
*        Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RESEQ
     *   (NSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,LSEQ,ISEQ,LEOF,
     *   INBS,OPTV,IRC)

      Include          'ardim.f'
      Include          'sterr.f'
     

* function return types

      Integer           Xblnk
      External          Xblnk

* reads sequence file in Swiss-Prot/EMBL format 

      Character*(*)     FSEQ
      Character         CABC(0:26)
      Character*(*)     CSID
      Character*(*)     CSAC
      Character*(*)     CSDE
      Integer*2         ISEQ(*)
      Logical           LEOF 
      Logical           LOPN
      Logical           LIDF
      Logical           OPTV
      Character*2       CMIS
      Character*256     RCIN
      
      IRC=0
      LIDF=.FALSE.
      LEOF=.FALSE.
      CMIS='ID'
      Inquire(File=FSEQ,OPENED=LOPN)
      If(LOPN) go to   1
      If(FSEQ.NE.'-') Open(NSEQ,File=FSEQ,Status='OLD',Err=900)
 1    Read(NSEQ,'(A)',Err=901,End=999,Iostat=IOS) RCIN
      If(RCIN(1:1).EQ.'>') go to 903
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:2).NE.'ID') go to   1   
      IC=Index(RCIN(6:17),' ')+5
      CSID=RCIN( 6:IC)
      LIDF=.TRUE.

      CMIS='AC'
 2    Read(NSEQ,'(A)',Err=901,End=902,Iostat=IOS) RCIN
      If(RCIN(1:2).NE.'AC'.AND.RCIN(1:2).NE.'DE') go to   2   
      If(RCIN(1:2).EQ.'AC') then 
         IX=Index(RCIN,';')-1
         If(IX.LE.0) go to 916
         CSAC=RCIN(6:IX) // '|'
      Else
         CSDE=RCIN( 6:Xblnk(RCIN,64))
         Go to   4
      End if 

      CMIS='DE'
 3    Read(NSEQ,'(A)',Err=901,End=902,Iostat=IOS) RCIN
      If(RCIN(1:2).NE.'DE') go to   3   
      CSDE=RCIN( 6:Xblnk(RCIN,64))

      CMIS='SQ'
 4    Read(NSEQ,'(A)',Err=901,End=902,Iostat=IOS) RCIN
      If(RCIN(1:2).NE.'DE') go to   5   
      LR=Xblnk(CSDE,64)
      If(CSDE(LR:LR).EQ.Char(13)) LR=LR-1
      JN=Xblnk(RCIN,64)
      If((LR+JN-4).LT.512)
     *   CSDE=CSDE(1:LR) // ' ' // RCIN( 6:JN)

 5    Continue
      If(RCIN(1:2).EQ.'SQ') go to   7

 6    Read(NSEQ,'(A)',Err=901,End=902,Iostat=IOS) RCIN
      If(RCIN(1:2).NE.'SQ') go to   6   

 7    Continue

      J1=0
 10   Read(NSEQ,'(A)',Err=901,End= 20,Iostat=IOS) RCIN
      If(RCIN(1:2).EQ.'//') go to  20

      Do 15 I1=1,Xblnk(RCIN,64)
         N1=0
         K1=Ichar(RCIN(I1:I1))
         If(K1.GE.97) then 
            K1=K1-32
            RCIN(I1:I1)=Char(K1)
         End if
         If(K1.GT.90.OR.K1.LT.65) go to  15 

         Do 13 I2=1,NABC
            If(CABC(I2).EQ.RCIN(I1:I1)) then
               N1=I2
               Go to   14
            End if
 13      Continue
 14      J1=J1+1
         If(J1.GT.IDMS) go to 915
         ISEQ(J1)=N1
 15   Continue
      Go to  10 

 20   LSEQ=J1
      If(LSEQ.LE.0) go to 917
      INBS=INBS+1

 100  If(IOS.EQ.-1.AND.IRC.LT.1) then
         LEOF=.TRUE.
         If(INBS.LT.1) go to 921
      End if

 105  Return

 110  If(LIDF) then
         Write(NERR,*) '       While processing sequence ',
     *      CSID(1:Lblnk(CSID))
      End if
      Go to 100


* errors

 900  Write(NERR,*) 'Error: Unable to open sequence file'//
     *   ' ''',FSEQ(1:Xblnk(FSEQ,64)),'''.'
      IRC=1
      Go to 110
 901  Write(NERR,*) 'Error: Unable to read sequence file'//
     *   ' ''',FSEQ(1:Xblnk(FSEQ,64)),'''.'
      IRC=1
      Go to 110
 902  Write(NERR,*) 'Error: Unexpected end of file. '//
     *   'Unable to find sequence ',CMIS,' line.'
      Write(NERR,*) '       No sequence read.'
      IRC=1
      Go to 110
 903  Write(NERR,*) 'Error: Sequence file in wrong format.'//
     *   ' Try option -f.'
      IRC=1
      Go to 110
 915  Write(NERR,*) 'Error: Sequence length exceeds buffer ',
     *   'size (',IDMS,').'
      IRC=1
      Go to 110
 916  Write(NERR,*) 'Error: Missing semicolon.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 917  Write(NERR,*) 'Error: Unexpected end of sequence.'//
     *   ' Sequence has zero length.'
      IRC=1
      Go to 110

* warnings

 920  If(.NOT.OPTV) then
         Write(NERR,*)'Warning: Unexpected sequence separator (''//'')',
     *      ' found. Unable to find sequence '//CMIS//' line.'
         If(LIDF) then
            Write(NERR,*) '       While processing sequence ',
     *         CSID(1:Lblnk(CSID))
         End if
         Write(NERR,*) '         Reading next sequence.'
      End if
      IRC=0
      LIDF=.FALSE.
      CMIS='ID'
      Go to 1
 921  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: No sequence has been read.'
         Write(NERR,*) '         No alignment produced.'
      End if
      IRC=1
      Go to 105

 999  IRC=-1
      Go to 100

      End

