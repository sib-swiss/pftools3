*----------------------------------------------------------------------*     
* $Id: rfseq.f,v 2.10 2003/11/18 10:50:07 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RFSEQ
     *   (NSEQ,FSEQ,NABC,CABC,CSID,CSAC,CSDE,CSFH,LSEQ,ISEQ,LEOF,
     *   INBS,OPTV,IRC)
      
      Include          'ardim.f'
      Include          'sterr.f'


* function return types

      Integer           Xblnk
      External          Xblnk

* reads sequence file in Pearson Fasta format 

      Character*(*)     FSEQ
      Character         CABC(0:26)
      Character*(*)     CSID
      Character*(*)     CSAC
      Character*(*)     CSFH
      Character*(*)     CSDE
      Integer*2         ISEQ(*)
      Logical           LEOF
      Logical           LOPN
      Logical           OPTV
      Character*512     RCIN
      
      IRC=0
      LEOF=.FALSE.
      CSID=' '
      CSAC=' '
      CSDE=' '

      If(RCIN(1:1).EQ.'>') Go to  2

      Inquire(File=FSEQ,OPENED=LOPN)
      If(LOPN) go to   1
      If(FSEQ.NE.'-') Open(NSEQ,File=FSEQ,Status='OLD',Err=900)
 1    Read(NSEQ,'(A)',Err=901,End=999,Iostat=IOS) RCIN
      If(RCIN(1:1).NE.'>') go to   1   

 2    L=Xblnk(RCIN,64)
      IX2=Index(RCIN(1:L),' ')-1
      If(IX2.LE.0) IX2=L
      Do I1=IX2,2,-1
         If(Index(':;|',RCIN(I1:I1)).NE.0) go to   3
      End do
 3    IX1=I1
      Do I1=IX2+1,L
         If(RCIN(I1:I1).NE.' ') go to  4
      End do
 4    IX3=I1 
      CSAC=RCIN(2:IX1)
*      FIXME
*      Write(NERR,*) 'SEQUENCE_READ=',CSAC
      CSID=RCIN(IX1+1:IX2)
      CSDE=RCIN(IX3:L)
      CSFH=RCIN(2:IX2)

      J1=0
 10   Read(NSEQ,'(A)',Err=901,End= 20,Iostat=IOS) RCIN
      L=Xblnk(RCIN,64)
      If(RCIN(1:1).EQ.'>') go to  20

      Do 15 I1=1,L
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
      If(LSEQ.LE.0) go to 916
      INBS=INBS+1
 100  If(IOS.EQ.-1.AND.IRC.LT.1) then
         LEOF=.TRUE.
         If(INBS.LT.1) go to 920
      End if

 105  Return

 110  IX=Lblnk(CSID)
      If(IX.GT.1) then
         Write(NERR,*) '       While processing sequence ',
     *      CSID(1:IX)
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
 915  Write(NERR,*) 'Error: Sequence length exceeds buffer ',
     *   'size (',IDMS,').'
      IRC=1
      Go to 110
 916  Write(NERR,*) 'Error: Unexpected end of sequence.'//
     *   ' Sequence has zero length.'
      IRC=1
      Go to 110

 920  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: No sequence has been read.'
         Write(NERR,*) '         No alignment produced.'
      End if
      IRC=-1
      Go to 105

 999  IRC=-1
      Go to 100 
      End

