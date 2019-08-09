*----------------------------------------------------------------------*     
* $Id: reprf.f,v 2.17 2003/11/18 10:50:07 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine REPRF
     *   (NPRF,FPRF,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   BLOG,FABC,P0,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,CHIL,IIPL,ILIP,
     *   CHID,IIPD,CHMD,IMPD,
     *   INBP,LEOF,OPTV,IRC)

* profile data fields

      Character*(*)     FPRF

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'dfdat.f'
      Include          'pfdat.f'
      Include          'sterr.f'

      Character*512     RCIN


* function return types

      Integer           Xblnk
      External          Xblnk

* work fields

      Integer           NKEY 
      Character*64      CPAR
      Character*512     CVAL 
      Logical           LNEW
      Logical           LIDF
      Logical           LZCO
      Integer           ISCO(26)

      Character*64      CH64
      Character*2       CMIS
      Logical           LOPN
      Logical           LEOF
      Logical           OPTV
      
      IRC=0
      RCIN=' '
      LZCO=.FALSE.
      LEOF=.FALSE.

* open profile file 

      Inquire(File=FPRF,OPENED=LOPN)
      If(LOPN) go to   1
      If(FPRF.NE.'-'.OR.NPRF.NE.5)
     *   Open(NPRF,File=FPRF,Status='OLD',Err=900)

* profile-id 
      
 1    LIDF=.FALSE.
      CMIS='ID'
      Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:2).NE.'ID') go to   1   
      LR=Xblnk(RCIN,64)
      IC=Index(RCIN,';')
      If(IC.LE.0) go to 916
      CPID=RCIN( 6:IC-1)
      If(Index(RCIN(IC:LR),'MATRIX').EQ.0) go to 924 
      LIDF=.TRUE.

* ac-number 
      
      CMIS='AC'
 2    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:2).NE.'AC') go to   2   
      IX=Index(RCIN,';')-1
      If(IX.LE.0) go to 916
C      CPAC=RCIN( 6:IC-1)
      CPAC=RCIN(6:IX) // '|'

* date, description 

      CMIS='DE'
      CPDT=' '
 3    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920
      If(RCIN(1:2).EQ.'DT') then 
         LR=Xblnk(RCIN,64)
         CPDT=RCIN( 6:LR)
         Go to   3
      End if 
      If(RCIN(1:2).NE.'DE') go to   3   
      LR=Xblnk(RCIN,64)
      CPDE=RCIN( 6:LR)

* go to first MA line, read additional header lines 
      
      LHDR=0
      CMIS='MA'
 5    Read(NPRF,'(A)',Err=901,End=902) RCIN
      If(RCIN(1:2).EQ.'//') go to 920 
      If(RCIN(1:2).NE.'MA') then
         If(LHDR.LT.512) then
            LHDR=LHDR+1
            CHDR(LHDR)=RCIN
         End if
         go to   5   
      End if
*      LR=Lblnk(RCIN)

* initialize position-independent profile parameters

* - general specifications

      CVAL='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      NABC=26
      Read(CVAL,'(26A)')(CABC(ii1),ii1=1,NABC)

      LPRF=0
      LENP=0

      LPCI=.FALSE.
      P0=0.0

* - disjoint mode

      MDIS=1

* - normalization modes

      JNOR=0 
      Do  11 I1=1,MAXN
         NNOR(I1)=I1
         NNPR(I1)=I1
 11   Continue

* - cut-off

      JCUT=0

* - defaults for match and insert position  

      JI=-1
      JM=0
      
      Call INITD(IIPD,IMPD,CHMD,CHID,NLOW)

      Do  18 I1=0,27
         IMPP(I1,0)=NLOW
 18   Continue 

* intialize profile input buffers.

      LR=Xblnk(RCIN,64)
      If(Ichar(RCIN(LR:LR)).EQ.13)LR=LR-1 
      JR=5

* read next parameter 

 20   Continue
      LFTR=0
      Call NEXTP(NPRF,RCIN,LR,JR,NKEY,LNEW,CPAR,CVAL,LFTR,CFTR,IRC) 
      If(IRC.GT.0) go to 940
      If(IRC.LT.0) go to 930

* interpret parameters

* - GENERAL_SPEC:

      If     (NKEY.EQ.1) then 

         If     (CPAR.EQ.'ALPHABET') then
            Read(CVAL,*,Err=910) CH64
            NABC=Lblnk(CH64)
            If(NABC.GT.26) Go to 908
            Read(CH64,'(64A)',Err=910)(CABC(ii1),ii1=1,NABC)

            If(NABC.LT.20) then
               CABC(0)='N'
            Else
               CABC(0)='X'
            End if
         Else if(CPAR.EQ.'LENGTH') then
            Read(CVAL,*,Err=910) LENP
         Else if(CPAR.EQ.'TOPOLOGY'.AND.CVAL.EQ.'CIRCULAR') then 
            LPCI=.TRUE.
         Else if(CPAR.EQ.'LOG_BASE') then
            Read(CVAL,*,Err=910) BLOG
            If(P0.EQ.0.0) P0=1.0
         Else if(CPAR.EQ.'P0') then
            Read(CVAL,*,Err=910) P0 
         Else if(CPAR.EQ.'P') then
            L1=Xblnk(CVAL,64) 
            J1=1
            Do  I1=1,L1
               If(CVAL(I1:I1).EQ.',') J1=J1+1
            End do 
            If(J1.GT.26) Go to 908
            Read(CVAL,*,Err=910)(FABC(ii1),ii1=1,J1)
         End if
         
* - DISJOINT:

      Else if(NKEY.EQ.2) then 

         If     (CPAR.EQ.'DEFINITION') then 
            Do  31 I1=1,KDIS
               If(CVAL.EQ.CDIS(I1)) MDIS=I1
 31         Continue
         Else
            If     (CPAR.EQ.'N1') then
               Read(CVAL,*,Err=910) NDIP(1) 
            Else if(CPAR.EQ.'N2') then
               Read(CVAL,*,Err=910) NDIP(2) 
            End if
         End if           

* - NORMALIZATION:

      Else if(NKEY.EQ.3) then 

         If(LNEW) then
            If(JNOR.GE.MAXN) Go to 912
            JNOR=JNOR+1
            CNTX(JNOR)=' '
         End if

         If     (CPAR.EQ.'FUNCTION') then
            If(CVAL.EQ.'GRIBSKOV') CVAL='GLE_ZSCORE'
            Do  41 I1=1,KNOR
               If(CVAL.EQ.CNOR(I1)) MNOR(JNOR)=I1
 41         Continue
            If (MNOR(JNOR).EQ.0) Go to 911
         Else if(CPAR.EQ.'MODE') then
            Read(CVAL,*,Err=910) NNOR(JNOR)
         Else if(CPAR.EQ.'PRIORITY') then
            Read(CVAL,*,Err=910) NNPR(JNOR)
         Else if(CPAR.EQ.'TEXT') then
            Read(CVAL,*,Err=910) CNTX(JNOR)
         Else
            If     (CPAR.EQ.'R1') then
               Read(CVAL,*,Err=910) RNOP(1,JNOR) 
            Else if(CPAR.EQ.'R2') then
               Read(CVAL,*,Err=910) RNOP(2,JNOR) 
            Else if(CPAR.EQ.'R3') then
               Read(CVAL,*,Err=910) RNOP(3,JNOR) 
            Else if(CPAR.EQ.'R4') then
               Read(CVAL,*,Err=910) RNOP(4,JNOR) 
            Else if(CPAR.EQ.'R5') then
               Read(CVAL,*,Err=910) RNOP(5,JNOR) 
            End if
         End if

* - CUT_OFF:

      Else if(NKEY.EQ.4) then 

         If(LNEW) then 
            If(JCUT.GE.MAXC) Go to 913
            JCUT=JCUT+1
            MCLE(JCUT)=0
            CCUT(JCUT)=' '
            JCNM(JCUT)=0
         End if

         If     (CPAR.EQ.'LEVEL') then
            Read(CVAL,*,Err=910) MCLE(JCUT)
            If(MCLE(JCUT).EQ.0) LZCO=.TRUE.
         Else if(CPAR.EQ.'SCORE') then
            Read(CVAL,*,Err=910) ICUT(JCUT)
         Else if(CPAR.EQ.'TEXT') then
            Read(CVAL,*,Err=910) CCUT(JCUT)
         Else if(CPAR.EQ.'N_SCORE') then
            L1=Xblnk(CVAL,64) 
            J1=1
            Do  51 I1=1,L1
               If(CVAL(I1:I1).EQ.',') J1=J1+1
 51         Continue
            Read(CVAL,*,Err=910)(RCUT(ii1,JCUT),ii1=1,J1)
            JCNM(JCUT)=J1
         Else if(CPAR.EQ.'MODE') then
            L1=Xblnk(CVAL,64) 
            J1=1
            Do  52 I1=1,L1
               If(CVAL(I1:I1).EQ.',') J1=J1+1
 52         Continue
            Read(CVAL,*,Err=910)(MCUT(ii1,JCUT),ii1=1,J1)
         End if

* - DEFAULT:

      Else if(NKEY.EQ.5) then 

         If     (CPAR.EQ.' ') then
            Call INITD(IIPD,IMPD,CHMD,CHID,NLOW)
         Else if(CPAR.EQ.'SY_I') then
            Read(CVAL,*,Err=910) CHID 
         Else if(CPAR.EQ.'SY_M') then
            Read(CVAL,*,Err=910) CHMD
         Else if(CPAR.EQ.'M0'.OR.CPAR.EQ.'M'.OR.CPAR.EQ.'D') then
            Call GETPM(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)
            If(IRC.NE.0) Go to 950
            If(INDX.EQ.1) then
               Do  61 I1=1,NABC
                  IMPD(I1)=ISCO(I1)
 61            Continue
            Else
               IMPD(INDX)=ISCO(1)
            End if
         Else
            Call GETPI(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)
            If(IRC.NE.0) Go to 950
            If(INDX.EQ.1) then
               Do  62 I1=1,NABC
                  IIPD(I1)=ISCO(I1)
 62            Continue
            Else
               IIPD(INDX)=ISCO(1)
            End if
         End if
         
* - I:

      Else if(NKEY.EQ.6) then 

         If(LNEW) then
            JI=JI+1

            If (JI.GE.IDMP.OR.JM.GT.IDMP) then
               Go to 915
            End if

            CHIP(JI)=CHID 
            Do  71 I1=0,46
               IIPP(I1,JI)=IIPD(I1)
 71         Continue
            If(JI.GT.JM) then
               JM=JM+1
               CHMP(JM)=CHMD 
               Do  72 I1=0,27
                  IMPP(I1,JM)=IMPD(I1)
 72            Continue
            End if
         End if 

         If     (CPAR.EQ.' ') then
            Go to  20
         Else if(CPAR.EQ.'SY') then
            Read(CVAL,*,Err=910) CHIP(JI) 
         Else
            Call GETPI(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)
            If(IRC.NE.0) Go to 950
            If(INDX.EQ.1) then
               Do  73 I1=1,NABC
                  IIPP(I1,JI)=ISCO(I1)
 73            Continue
            Else
               IIPP(INDX,JI)=ISCO(1)
            End if
         End if

* - M:

      Else if(NKEY.EQ.7) then 

         If(LNEW) then
            JM=JM+1

            If (JI.GE.IDMP.OR.JM.GT.IDMP) then
               Go to 915
            End if

            CHMP(JM)=CHMD 
            Do  81 I1=0,27
               IMPP(I1,JM)=IMPD(I1)
 81         Continue
            If(JM-1.GT.JI) then
               JI=JI+1
               CHIP(JI)=CHID 
               Do  82 I1=0,46
                  IIPP(I1,JI)=IIPD(I1)
 82            Continue
            End if
         End if 

         If     (CPAR.EQ.' ') then
            Go to  20
         Else if(CPAR.EQ.'SY') then
            Read(CVAL,*,Err=910) CHMP(JM) 
         Else
            Call GETPM(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)
            If(IRC.NE.0) Go to 950
            If(INDX.EQ.1) then
               Do  83 I1=1,NABC
                  IMPP(I1,JM)=ISCO(I1)
 83            Continue
            Else
               IMPP(INDX,JM)=ISCO(1)
            End if
         End if
      Else 
         Continue
      End if  

C      Print *,JI,JM,NKEY,' ',LNEW,' ',
C     *   CPAR(1:Lblnk(CPAR)),'=',
C     *   CVAL(1:Lblnk(CVAL))

      Go to  20

 90   Continue
      IRC=0

* last insert position 
* store values for multiple matches of circular profiles

      If(LPCI) then
         If(JI.EQ.JM) then 
            ILIP=0
            CHIL=CHIP(0)
            CHIP(    0)=CHIP(   JI)
            Do  91 I1=0,46
               IIPL(I1)=IIPP(I1, 0)
               IIPP(I1, 0)=IIPP(I1,JI)
 91         Continue
         Else
            JI=JI+1
            ILIP=JI
            CHIL=CHID
            CHIP(   JI)=CHIP(    0)
            Do  92 I1=0,46
               IIPL(I1)=IIPD(I1)
               IIPP(I1,JI)=IIPP(I1, 0)
 92         Continue
         End if
      Else if(JI.LT.JM) then 
         JI=JI+1
         CHIP(   JI)=CHID
         Do  93 I1=0,46
            IIPP(I1,JI)=IIPD(I1)
 93      Continue
      End if


* set length of profile
 
      LPRF=JI

* check consistency of some parameters

 94   If(.NOT.LZCO.OR.JCUT.EQ.0) go to 919
      If(LPRF.LT.1) then
         Go to 918
      Else
         INBP=INBP+1
      End if
      If(NDIP(1).LE.0.OR.NDIP(1).GT.LPRF) then
         NDIP(1)=1
         Go to 923
      End if
      If(NDIP(2).LE.0.OR.NDIP(2).GT.LPRF) then
         NDIP(2)=LPRF
         Go to 923
      End if
      If(NDIP(2).LT.NDIP(1)) then
         I1=NDIP(2)
         NDIP(2)=NDIP(1)
         NDIP(1)=I1
      End if
C      Do 95 I1=0,JI
C         Write(NERR,*) 'I[',I1,']: '
C         Write(NERR,*)(IIPP(ii1,I1),ii1=0,46)
C 95   Continue
C
C      Do 96 I1=0,JM
C         Write(NERR,*) 'M[',I1,']: '
C         Write(NERR,*)(IMPP(ii1,I1),ii1=0,27)
C 96   Continue

 100  Return

 110  If(LIDF.AND.IRC.NE.0) then
         Write(NERR,*) '       While processing profile ',
     *      CPID(1:Lblnk(CPID))
      End if
      Go to 100

* errors

 900  Write(NERR,*) 'Error: Unable to open profile file'//
     *   ' ''',FPRF(1:Xblnk(FPRF,64)),'''.'
      IRC=1
      Go to 110
 901  Write(NERR,*) 'Error: Unable to read profile file'//
     *   ' ''',FPRF(1:Xblnk(FPRF,64)),'''.'
      IRC=1
      Go to 110
 902  If(INBP.LT.1.OR.LIDF) then
         Write(NERR,*) 'Error: Unexpected end of file. '//
     *      'Unable to find profile ',CMIS,' line.'
         If(INBP.LT.1)
     *      Write(NERR,*) '       No profile read.'
         IRC=1
      Else
         LEOF=.TRUE.
         IRC=0
      End if
      Go to 110
 906  Write(NERR,*) 'Error: Unknown data block at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 908  Write(NERR,*) 'Error: Too many values for parameter '//
     *   CPAR(1:Lblnk(CPAR))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 910  Write(NERR,*) 'Error: Unable to read value of parameter '//
     *   CPAR(1:Lblnk(CPAR))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 911  Write(NERR,*) 'Error: Unknown value '//
     *   CVAL(1:Xblnk(CVAL,64)),' for parameter '//
     *   CPAR(1:Lblnk(CPAR))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 912  Write(NERR,*) 'Error: Maximum number of NORMALIZATION '//
     *   'data blocks exceeded '
      IRC=1
      Go to 110
 913  Write(NERR,*) 'Error: Maximum number of CUT_OFF '//
     *   'data blocks exceeded '
      IRC=1
      Go to 110
 914  Write(NERR,*) 'Error: Unknown parameter '//
     *   CPAR(1:Lblnk(CPAR))
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 915  Write(NERR,*) 'Error: Profile length exceeds buffer size (',
     *   IDMP,').'
      IRC=1
      Go to 110
 916  Write(NERR,*) 'Error: Missing semicolon.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 917  Write(NERR,*) 'Error: Missing ''MATRIX'' keyword in ID line.'
      Write(NERR,*) '       at line: ',
     *   RCIN(1:Xblnk(RCIN,64))
      IRC=1
      Go to 110
 918  Write(NERR,*) 'Error: Unexpected end of profile. Profile has '//
     *   'zero length. Check profile syntax.'
      IRC=1
      Go to 110

 919  Write(NERR,*) 'Error: No level 0 ''CUT-OFF'' data-block'//
     *   ' in profile.'
      Write(NERR,*) '       Check profile syntax.'
      IRC=1
      Go to 110

* warnings

 920  If(.NOT.OPTV) then
         Write(NERR,*)'Warning: Unexpected profile separator (''//'')'//
     *      ' found. Unable to find ',CMIS,' line.'
         If(LIDF) then
            Write(NERR,*) '         While processing profile ',
     *         CPID(1:Lblnk(CPID))
         End if
         Write(NERR,*) '         Reading next profile.'
      End if
      IRC=0
      Go to 1

 921  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: Parameter name too long at line: ',
     *      RCIN(1:Xblnk(RCIN,64))
         If(LIDF) then
            Write(NERR,*) '         While processing profile ',
     *         CPID(1:Lblnk(CPID))
         End if
         Write(NERR,*) '         Skipping parameter.'
      End if
      IRC=0
      Go to 20
      
 922  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: Parameter value too long at line: ',
     *      RCIN(1:Xblnk(RCIN,64))
         If(LIDF) then
            Write(NERR,*) '         While processing profile ',
     *         CPID(1:Lblnk(CPID))
         End if
         Write(NERR,*) '         Skipping parameter.'
      End if
      IRC=0
      Go to 20

 923  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: Disjointness parameter out of bound. '
         Write(NERR,*) '         Setting parameter to acceptable value.'
         If(LIDF) then
            Write(NERR,*) '         While processing profile ',
     *         CPID(1:Lblnk(CPID))
         End if
      End if
      Go to 94

 924  If(.NOT.OPTV) then
         Write(NERR,*) 'Warning: Missing ''MATRIX'' keyword in ID line.'
         Write(NERR,*) '         at line: ',
     *      RCIN(1:Xblnk(RCIN,64))
         IRC=0
         If(LIDF) then
            Write(NERR,*) '         While processing profile ',
     *         CPID(1:Lblnk(CPID))
         End if
         Write(NERR,*) '         Reading next profile.'
      End if
      IRC=0
      Go to 1


* process warnings occuring in NEXTP()

 930  If(IRC.EQ.-3) Go to 921
      If(IRC.EQ.-4) Go to 922
      Go to 90

* process errors occuring in NEXTP()

 940  If(IRC.EQ.2) Go to 906
      Go to 901

* process errors occuring in GETPI() or GETPM()

 950  If(IRC.GT.0) Go to 910
      If(IRC.EQ.-1) Go to 914
      If(IRC.EQ.-2) Go to 908
      Go to 90

      End


*----------------------------------------------------------------------*
* Initializes all parameters to their implicit defaults
*----------------------------------------------------------------------*

      Subroutine INITD(IIPD,IMPD,CHMD,CHID,NLOW)

      Integer           B0 
      Integer           B1 
      Integer           E0 
      Integer           E1 
      Integer           BM 
      Integer           BI 
      Integer           BD 
      Integer           BE 
      Integer           MM 
      Integer           MI 
      Integer           MD 
      Integer           ME 
      Integer           DM 
      Integer           DI 
      Integer           DD 
      Integer           DE 
      Integer           IM 
      Integer           II 
      Integer           ID 
      Integer           IE 

      Integer           M0
      Integer           I0
      Integer           D 

      Parameter        (B0=  27) 
      Parameter        (B1=  28)
      Parameter        (E0=  29)
      Parameter        (E1=  30)
      Parameter        (BM=  31)
      Parameter        (BI=  32)
      Parameter        (BD=  33)
      Parameter        (BE=  34)
      Parameter        (MM=  35)
      Parameter        (MI=  36)
      Parameter        (MD=  37)
      Parameter        (ME=  38)
      Parameter        (IM=  39)
      Parameter        (II=  40)
      Parameter        (ID=  41)
      Parameter        (IE=  42)
      Parameter        (DM=  43)
      Parameter        (DI=  44)
      Parameter        (DD=  45)
      Parameter        (DE=  46)

      Parameter        (M0=   0)
      Parameter        (I0=   0)
      Parameter        (D =  27)

      Integer           IIPD(0:46)
      Integer           IMPD(0:27)
      Character         CHMD
      Character         CHID

      CHID='-'
      Do  10 I1=1,26 
         IIPD(I1)=0
 10   Continue

      IIPD(B0)=0
      IIPD(B1)=0
      IIPD(E0)=0
      IIPD(E1)=0

      IIPD(BM)=0
      IIPD(BI)=NLOW
      IIPD(BD)=NLOW
      IIPD(BE)=NLOW
      IIPD(MM)=0
      IIPD(MI)=NLOW
      IIPD(MD)=NLOW
      IIPD(ME)=0
      IIPD(IM)=NLOW
      IIPD(II)=0
      IIPD(ID)=NLOW
      IIPD(IE)=NLOW
      IIPD(DM)=NLOW
      IIPD(DI)=NLOW
      IIPD(DD)=0
      IIPD(DE)=NLOW

      IIPD(I0)=0

      CHMD='X'
      Do  15 I1=1,26 
         IMPD(I1)=0
 15   Continue

      IMPD(M0)=0 
      IMPD(D )=0

 100  Return

      End

*----------------------------------------------------------------------*
      Subroutine NEXTP(NPRF,RCIN,LR,JR,NKEY,LNEW,CPAR,CVAL,LFTR,
     *   CFTR,IRC) 

* function return types

      Integer           Xblnk
      External          Xblnk

      Character*(*)     RCIN
      Character*(*)     CPAR 
      Character*(*)     CVAL 
      Character*(*)     CFTR(1024)
      Logical           LNEW 

      LNEW=.FALSE.
 3    Do  4 I1=JR+1,LR
         JR=JR+1
         If(Index(' ;,.        ',RCIN(JR:JR)).EQ.0) go to 10
 4    Continue

 5    Read(NPRF,'(A)',Err=999,End=901) RCIN
      If     (RCIN(1:2).EQ.'MA') then
         LFTR=0
      Else if(RCIN(1:2).EQ.'CC') then
         If(LFTR.LT.1024) then 
            LFTR=LFTR+1
            CFTR(LFTR)=RCIN
         End if
         Go to   5 
      Else if(RCIN(1:2).EQ.'//') then
         Go to 902
      Else
         Go to 102
      End if 

      LR=Xblnk(RCIN,64)
      If(Ichar(RCIN(LR:LR)).EQ.13) LR=LR-1 
      JR=5
      Go to   3

* identify keyword 

 10   If(RCIN(JR:JR).EQ.'/') then 
         LNEW=.TRUE.
         If     (RCIN(JR:JR+13).EQ.'/GENERAL_SPEC:') then
            NKEY=1
            JR=JR+13
         Else if(RCIN(JR:JR+ 9).EQ.'/DISJOINT:') then
            NKEY=2
            JR=JR+ 9
         Else if(RCIN(JR:JR+14).EQ.'/NORMALIZATION:') then
            NKEY=3
            JR=JR+14
         Else if(RCIN(JR:JR+ 8).EQ.'/CUT_OFF:') then
            NKEY=4
            JR=JR+ 8
         Else if(RCIN(JR:JR+ 8).EQ.'/DEFAULT:') then
            NKEY=5
            JR=JR+ 8
         Else if(RCIN(JR:JR+ 2).EQ.'/I:') then
            NKEY=6
            JR=JR+ 2
         Else if(RCIN(JR:JR+ 2).EQ.'/M:') then
            NKEY=7
            JR=JR+ 2
         Else
            Go to 900
         End if
      Else 
         JR=JR-1
      End if 

* read parameter 

      LPAR=0
      CPAR=' '  

 20   JR=JR+1  

      If(JR.GT.LR) then 
 21      Read(NPRF,'(A)',Err=999,End=901) RCIN
         If(RCIN(1:2).EQ.'CC') then
            If(LFTR.LT.1024) then
               LFTR=LFTR+1
               CFTR(LFTR)=RCIN
            End if
            Go to  21
         End if
         LR=Xblnk(RCIN,64)
         If(Ichar(RCIN(LR:LR)).EQ.13) LR=LR-1 
         If(RCIN(1:2).NE.'MA') then 
            JR=LR
            CVAL=' '
            Go to  100
         End if
         LFTR=0
         JR=6
      End if  

      If(RCIN(JR:JR).EQ.'/') then
         JR=JR-1
         CVAL=' '
         Go to  100
      End if   

      If(RCIN(JR:JR).NE.'=') then
         If(LPAR.NE.0.OR.RCIN(JR:JR).NE.' ') then 
            LPAR=LPAR+1
            If(LPAR.LE.64) then
               CPAR(LPAR:LPAR)=RCIN(JR:JR)
            End if
         End if 
         Go to  20
      End if

* read value 

      LVAL=0
      CVAL=' '  

 30   JR=JR+1  

      If(JR.GT.LR) then 
 31      Read(NPRF,'(A)',Err=999,End=901) RCIN
         If(RCIN(1:2).EQ.'CC') then 
            If(LFTR.LT.1024) then
               LFTR=LFTR+1 
               CFTR(LFTR)=RCIN 
            End if
            Go to  31
         End if
         LR=Xblnk(RCIN,64)
         If(Ichar(RCIN(LR:LR)).EQ.13) LR=LR-1
         If(RCIN(1:2).NE.'MA') go to 102
         LFTR=0
         JR=6
      End if  

      If(Index('/;',RCIN(JR:JR)).EQ.0) then
         If(LVAL.NE.0.OR.RCIN(JR:JR).NE.' ') then 
            LVAL=LVAL+1
            If(LVAL.LE.512) then
               CVAL(LVAL:LVAL)=RCIN(JR:JR)
            End if
         End if
         Go to  30
      End if

* check parameter length against buffer size and exit
      
 100  If(LPAR.GT.64) Go to 903
      If(LVAL.GT.512) Go to 904
 101  Return

 102  If(LFTR.LT.1024) then 
         LFTR=LFTR+1
         CFTR(LFTR)=RCIN
      End if
 103  Read(NPRF,'(A)',Err=999,End=901) RCIN
      If(RCIN(1:2).NE.'//') then
         Go to 102
      Else 
         Go to 902
      End if 

* Error unknown data block

 900  IRC=2
      Go to 101

* Warn end of file

 901  IRC=-1
      Go to 101

* Warn profile separator found

 902  IRC=-2
      Go to 101

* Warn parameter name too long

 903  IRC=-3
      Go to 101

* Warn parameter value too long

 904  IRC=-4
      Go to 101

* Error reading file

 999  IRC=1
      Go to 101

      End
*----------------------------------------------------------------------*
      Integer Function PINTR(CVAL,NLOW,IRC)

      Integer           NLOW
      Character*(*)     CVAL 

      If(Index(CVAL,'*').NE.0) then
         PINTR=NLOW
      Else     
         Read(CVAL,*,Err= 101) PINTR
      End if
 100  Return
 101  IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine PINTS(NSCO,ISCO,CVAL,NLOW,NABC,IRC)

* function return types

      Integer           Xblnk
      External          Xblnk

      Integer           ISCO(*)
      Character*(*)     CVAL 
      Integer           NLOW

      L1=Xblnk(CVAL,64) 
      NSCO=0
      J1=1
      Do  10 I1=1,L1
         If(CVAL(I1:I1).EQ.',') then
            J2=I1-1
            NSCO=NSCO+1
            If(NSCO.GT.NABC) Go to 102
            If(J2.LT.J1) then 
               ISCO(NSCO)=0 
            Else
               If(Index(CVAL(J1:J2),'*').NE.0) then
                  ISCO(NSCO)=NLOW
               Else
                  Read(CVAL(J1:J2),*,Err=101) ISCO(NSCO) 
               End if
            End if 
            J1=I1+1
         End if
 10   Continue
      J2=L1
      NSCO=NSCO+1
      If(NSCO.GT.NABC) Go to 102
      If(J2.LT.J1) then 
         ISCO(NSCO)=0 
      Else
         If(Index(CVAL(J1:J2),'*').NE.0) then
            ISCO(NSCO)=NLOW
         Else
            Read(CVAL(J1:J2),*,Err=101) ISCO(NSCO) 
         End if
      End if 
      
 100  Return
      
 101  IRC=1
      Go to 100
 102  IRC=-2
      Go to 100
      End
*----------------------------------------------------------------------*
      Subroutine GETPI(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)

      Character*(*)     CPAR
      Character*(*)     CVAL
      Integer           ISCO(*)
      Integer           NLOW
      Integer           PINTR

      If     (CPAR.EQ.'I0') then
         INDX= 0
      Else if(CPAR.EQ.'I' ) then
         INDX= 1
      Else if(CPAR.EQ.'B0') then
         INDX=27 
      Else if(CPAR.EQ.'B1') then
         INDX=28
      Else if(CPAR.EQ.'E0') then
         INDX=29
      Else if(CPAR.EQ.'E1') then
         INDX=30

      Else if(CPAR.EQ.'BM') then
         INDX=31
      Else if(CPAR.EQ.'BI') then
         INDX=32
      Else if(CPAR.EQ.'BD') then
         INDX=33
      Else if(CPAR.EQ.'BE') then
         INDX=34

      Else if(CPAR.EQ.'MM') then
         INDX=35
      Else if(CPAR.EQ.'MI') then
         INDX=36
      Else if(CPAR.EQ.'MD') then
         INDX=37
      Else if(CPAR.EQ.'ME') then
         INDX=38

      Else if(CPAR.EQ.'IM') then
         INDX=39
      Else if(CPAR.EQ.'II') then
         INDX=40
      Else if(CPAR.EQ.'ID') then
         INDX=41
      Else if(CPAR.EQ.'IE') then
         INDX=42

      Else if(CPAR.EQ.'DM') then
         INDX=43
      Else if(CPAR.EQ.'DI') then
         INDX=44
      Else if(CPAR.EQ.'DD') then
         INDX=45
      Else if(CPAR.EQ.'DE') then
         INDX=46
      Else 
         Go to 900
      End if 

      If(INDX.EQ.1) then
         Call PINTS(NSCO,ISCO,CVAL,NLOW,NABC,IRC)
         If(NSCO.EQ.1) then
            Do  10 I1=2,NABC
               ISCO(I1)=ISCO( 1)
 10         Continue
         End if
      Else
         ISCO(1)=PINTR(CVAL,NLOW,IRC)
      End if 

 100  Return 
 900  IRC=-1
      Go to 100
      End 
*----------------------------------------------------------------------*
      Subroutine GETPM(CPAR,CVAL,INDX,ISCO,NABC,NLOW,IRC)

      Character*(*)     CPAR
      Character*(*)     CVAL
      Integer           ISCO(*)
      Integer           NLOW
      Integer           PINTR

      If     (CPAR.EQ.'M0') then
         INDX= 0
      Else if(CPAR.EQ.'M' ) then
         INDX= 1
      Else if(CPAR.EQ.'D ') then
         INDX=27 
      Else
         Go to 900
      End if

      If(INDX.EQ.1) then
         Call PINTS(NSCO,ISCO,CVAL,NLOW,NABC,IRC)
         If(NSCO.EQ.1) then
            Do  10 I1=2,NABC
               ISCO(I1)=ISCO( 1)
 10         Continue
         End if
      Else
         ISCO(1)=PINTR(CVAL,NLOW,IRC)
      End if 

 100  Return 
 900  IRC=-1
      Go to 100
      End 
