*----------------------------------------------------------------------*     
* $Id: wrprf.f,v 2.9 2003/11/18 10:52:29 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine WRPRF
     *   (NOUT,LLLT,LRNM,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT, 
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   BLOG,FABC,P0,
     *   CHID,IIPD,CHMD,IMPD, 
     *   IRC)

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'sterr.f'

* work fields
      Integer           SLEN
      Integer           LLEN
      Parameter        (SLEN=78)
      Parameter        (LLEN=132)
      Parameter        (LBUF=512)

      Character*512     CBLK
      Character*512     CPAR
      Character*512     RCEX
      Character*16      CHRP(MAXN) 

      Logical           LPRI
      Logical           LLLT
      Logical           LRNM

* initiation, transition and termination score names

      Character*2       CPSN(27:46)
      Data              CPSN(27)/'B0'/
      Data              CPSN(28)/'B1'/
      Data              CPSN(29)/'E0'/
      Data              CPSN(30)/'E1'/
      Data              CPSN(31)/'BM'/
      Data              CPSN(32)/'BI'/
      Data              CPSN(33)/'BD'/
      Data              CPSN(34)/'BE'/
      Data              CPSN(35)/'MM'/
      Data              CPSN(36)/'MI'/
      Data              CPSN(37)/'MD'/
      Data              CPSN(38)/'ME'/
      Data              CPSN(39)/'IM'/
      Data              CPSN(40)/'II'/
      Data              CPSN(41)/'ID'/
      Data              CPSN(42)/'IE'/
      Data              CPSN(43)/'DM'/
      Data              CPSN(44)/'DI'/
      Data              CPSN(45)/'DD'/
      Data              CPSN(46)/'DE'/
      
      IRC=0

* prevent meaningless parameters from being printed

      If(.NOT.LPCI) then 
         IIPP(MM,   0)=IIPD(MM)
         IIPP(MI,   0)=IIPD(MI)
         IIPP(MD,   0)=IIPD(MD)
         IIPP(ME,   0)=IIPD(ME)
         IIPP(DM,   0)=IIPD(DM)
         IIPP(DI,   0)=IIPD(DI)
         IIPP(DD,   0)=IIPD(DD)
         IIPP(DE,   0)=IIPD(DE)
         IIPP(MM,LPRF)=IIPD(MM)
         IIPP(IM,LPRF)=IIPD(IM)
         IIPP(DM,LPRF)=IIPD(DM)
         IIPP(BM,LPRF)=IIPD(BM)
         IIPP(MD,LPRF)=IIPD(MD)
         IIPP(ID,LPRF)=IIPD(ID)
         IIPP(DD,LPRF)=IIPD(DD)
         IIPP(BD,LPRF)=IIPD(BD)
      End if
      
* write header

* - ID line
      
      RCEX='ID   ' 
     *   // CPID(1:Lblnk(CPID))
     *   // '; MATRIX.'
      Call wline(NOUT,LLLT,SLEN,RCEX)
C      Write(NOUT,'(78A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

* - AC line 
      
      IX=Lblnk(CPAC)
      If(CPAC(IX:IX).EQ.'|') IX=IX-1
      RCEX='AC   '  // CPAC(1:IX) // ';'
      Call wline(NOUT,LLLT,SLEN,RCEX)
C      Write(NOUT,'(78A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

* - DT line
      
      If(CPDT.NE.' ') then
         RCEX='DT   ' 
     *      // CPDT(1:Lblnk(CPDT))
         Call wline(NOUT,LLLT,SLEN,RCEX)
C         Write(NOUT,'(132A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))
      End if 

* - DE line
      
      RCEX='DE   ' 
     *   // CPDE(1:Lblnk(CPDE))
      Call wline(NOUT,LLLT,LLEN,RCEX)
C      Write(NOUT,'(132A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))

* - Additional header lines

      Do I1=1,LHDR
         Call wline(NOUT,LLLT,LLEN,CHDR(I1))
C         Write(NOUT,'(132A)')
C     *      (CHDR(I1)(ii1:ii1),ii1=1,Lblnk(CHDR(I1)))
      End do 

* write /GENERAL_SPEC: block 

      CBLK='/GENERAL_SPEC:'   
      JB=14

* - alphapet
      
      Write(CPAR,'(''ALPHABET='''''',26A)')(CABC(ii1),ii1=1,NABC) 
      JP=12+NABC
      CPAR(JP-1:JP)=''';'
      If(JB+JP+1.GT.LBUF) Go to 900
      CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
      JB=JB+JP+1

* - length
C      Write(CPAR,'(''LENGTH='',I6,'';'')') LPRF
      CPAR='LENGTH='
      Write(CPAR(8:),*)LPRF,';'
      JP=Lblnk(CPAR)
      Call slpar(CPAR,JP)
      If(JB+JP+1.GT.LBUF) Go to 900
      CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
      JB=JB+JP+1

* - topology 

      If(LPCI) then
         If(JB+19.GT.LBUF) Go to 900
         CBLK=CBLK(1:JB) // ' TOPOLOGY=CIRCULAR;'
         JB=JB+19
      End if


* - HMM parameters:
*** To be modified ?? Using write(...,*)
      If(BLOG.NE.0) then 
         Write(CPAR,'(''LOG_BASE='',F8.6,'';'')') BLOG
         JP=18
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

         If(P0.NE.1.0) then
            Write(CPAR,'(''P0='',F6.4,'';'')') P0
            JP=10
            Call slpar(CPAR,JP)
            If(JB+JP+1.GT.LBUF) Go to 900
            CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
            JB=JB+JP+1
         End if 
      End if 
***
      Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)
      If(BLOG.NE.0) Call wrnul(NOUT,LRNM,LLLT,SLEN,NABC,FABC)

* write /DISJOINT: block 

      CBLK='/DISJOINT:'
      JB=10

* - definition

      CPAR='DEFINITION=' // CDIS(MDIS)
      JP=Lblnk(CPAR)+1
      CPAR(JP:JP)=';'
      If(JB+JP+1.GT.LBUF) Go to 900
      CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
      JB=JB+JP+1

* - parameters

      KDIP=JDIP(MDIS)
      JP=0
      Do 20 I1=1,KDIP
         Write(CPAR(JP+1:),*)'N',I1,'=',NDIP(I1),';'
         JP=Lblnk(CPAR)+1
 20   Continue
      Call slpar(CPAR,JP)
      If(JB+JP+1.GT.LBUF) Go to 900
      CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
      JB=JB+JP

      Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

* write /NORMALIZATION: block

      LPRI=.FALSE.
      Do  21 I1=1,JNOR
         If(NNOR(I1).NE.NNPR(I1)) LPRI=.TRUE.
 21   Continue

      Do  30 I1=1,JNOR

         CBLK='/NORMALIZATION:'
         JB=15

* - mode 

         Write(CPAR(1:),*)'MODE=',NNOR(I1),';' 
         JP=Lblnk(CPAR)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

* - priority

         If(LPRI) then
            Write(CPAR(1:),*)'PRIORITY=',NNPR(I1),';' 
            JP=Lblnk(CPAR)
            Call slpar(CPAR,JP)
            If(JB+JP+1.GT.LBUF) Go to 900
            CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
            JB=JB+JP+1
         End if         

* - function

         CPAR='FUNCTION=' // CNOR(MNOR(I1))
         JP=Lblnk(CPAR)+1
         CPAR(JP:JP)=';'
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

* - paramaters
         
         KNOP=JNOP(MNOR(I1))
         JP=0
         Do  22 I2=1,KNOP
            If(ABS(RNOP(I2,I1)).LE.10.0
     *         .AND.ABS(RNOP(I2,I1)).GT.0.0001
     *         .OR.RNOP(I2,I1).EQ.0.0) then
               Write(CHRP(I2),'(F10.7)') RNOP(I2,I1)
            Else 
               Write(CHRP(I2),*) RNOP(I2,I1)
            End if
            Write(CPAR(JP+1:),*)'R',I2,'=',CHRP(I2),'; '
            JP=Lblnk(CPAR)+1
 22      Continue
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP

* - text

         If(CNTX(I1).NE.' ') then 
            CPAR='TEXT=''' // CNTX(I1)(1:Lblnk(CNTX(I1))) // ''''
            JP=Lblnk(CPAR)+1
            CPAR(JP:JP)=';'
            If(JB+JP+1.GT.LBUF) Go to 900
            CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
            JB=JB+JP+1
         End if

         Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

 30   Continue

* write /CUT_OFF: block

      Do  40 I1=1,JCUT

         CBLK='/CUT_OFF:'
         JB=9

* - level

         Write(CPAR,*)'LEVEL=',MCLE(I1),';' 
         JP=Lblnk(CPAR)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

* - score

         Write(CPAR,*)'SCORE=',ICUT(I1),';' 
         JP=Lblnk(CPAR)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

* - normalized scores
         
         KCNM=JCNM(I1)
         If(KCNM.GT.0) CPAR='N_SCORE='
         JP=8
         Do 32 I2=1,KCNM
            If(I2.GT.1) then
               JP=JP+2
               CPAR(JP:JP)=','
            End if
            Write(CPAR(JP+1:),*)RCUT(I2,I1)
            JP=Lblnk(CPAR)
C         Write(CPAR,*)'(''N_SCORE='',G12.5,7('','',G12.5))')
C     *      (RCUT(ii1,I1),ii1=1,KCNM)
C         JP=8+13*KCNM
 32      Continue
         JP=JP+2
         CPAR(JP:JP)=';'
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1

* - normalization modes
         
         If(KCNM.GT.0) CPAR='MODE='
         JP=5
         Do 34 I2=1,KCNM
            If(I2.GT.1) then
               JP=JP+2
               CPAR(JP:JP)=','
            End if
            Write(CPAR(JP+1:),*)MCUT(I2,I1)
            JP=Lblnk(CPAR)
 34      Continue
C         Write(CPAR,'(''MODE='',I6,7('','',I6))')
C     *      (MCUT(ii1,I1),ii1=1,KCNM)
C         JP=5+7*KCNM
         JP=JP+2
         CPAR(JP:JP)=';'
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
         
* - text

         If(CCUT(I1).NE.' ') then 
            CPAR='TEXT=''' // CCUT(I1)(1:Lblnk(CCUT(I1))) // ''''
            JP=Lblnk(CPAR)+1
            CPAR(JP:JP)=';'
            If(JB+JP+1.GT.LBUF) Go to 900
            CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
            JB=JB+JP+1
         End if

         Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

 40   Continue

* write /DEFAULT: block 

* - transfer defaults to pos. LPRF+1

      I1=LPRF+1

      Do  51 I2=0,46
         IIPP(I2,I1)=IIPD(I2)
 51   Continue
      CHIP(I1)=CHID
      Do  52 I2=0,27
         IMPP(I2,I1)=IMPD(I2)
 52   Continue
      CHMP(I1)=CHMD

* - reinitialize defaults for match and insert position  

      CHID='-'
      Do  53 I2=1,26 
         IIPD(I2)=0
 53   Continue
      
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
      Do  54 I2=1,26 
         IMPD(I2)=0
 54   Continue

      IMPD(M0)=0 
      IMPD(D )=0

      Do  55 I2=0,27
         IMPP(I2,0)=NLOW
 55   Continue 

* - symbols

      CBLK='/DEFAULT:'
      JB=9

      If(CHMP(I1).NE.CHMD) then
         Write(CPAR,'(''SY_M='''''',A,'''''';'')') CHMP(I1)
         JP=9
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if 

      If(CHIP(I1).NE.CHID) then
         Write(CPAR,'(''SY_I='''''',A,'''''';'')') CHIP(I1)
         JP=9
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if 

* - match extension scores

      K1=0
      Do  56 I2=1,NABC
         If(IMPP(I2,I1).NE.IMPD(I2)) K1=1
 56   Continue 
      Do  57 I2=1,NABC
         If(IMPP(I2,I1).NE.IMPP( 1,I1)) K1=2
 57   Continue 

      If(K1.EQ.0) then
         JP=0 
      Else
         CPAR='M='
         JP=3
         If(K1.EQ.1) then
            If(IMPP( 1,I1).LE.NLOW) then
               CPAR(JP:)='*'
            Else 
               Write(CPAR(JP:),*)IMPP( 1,I1)
            End if 
            JP=Lblnk(CPAR)+1
            CPAR(JP:JP)=';'
         Else
            Do 58 I2=1,NABC
               If(I2.GT.1) then
                  JP=JP+2
                  CPAR(JP:JP)=','
               End if
               If(IMPP(I2,I1).LE.NLOW) then 
                  CPAR(JP+1:)='*'
               Else
                  Write(CPAR(JP+1:),*)IMPP(I2,I1)
               End if
               JP=Lblnk(CPAR)
 58         Continue
            JP=JP+1
            CPAR(JP:JP)=';'


C            Write(CPAR,'(''M='',I6,25('','',I6))')
C     *         (IMPP(ii1,I1),ii1=1,NABC)
C            JP=2+NABC*7
C            CPAR(JP:JP)=';'
C            J2=JP-NABC*7+1
C            Do  59 I2=1,NABC
C               If(IMPP(I2,I1).EQ.NLOW) CPAR(J2:J2+5)='     *'
C               J2=J2+7
C 59         Continue
         End if

         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if

* match extension score for unknown character

      If(IMPP( 0,I1).NE.IMPD( 0)) then
         CPAR='M0='
         JP=4
         If(IMPP( 0,I1).LE.NLOW) then
            CPAR(JP:)='*'
         Else
            Write(CPAR(JP:),*)IMPP( 0,I1)
         End if
         JP=LBlnk(CPAR)+1
         CPAR(JP:JP)=';'
C         Call slpar(CPAR,JP)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if 

* - deletion extension score

      If(IMPP(27,I1).NE.IMPD(27)) then
         CPAR='D='
         JP=3
         If(IMPP(27,I1).LE.NLOW) then
            CPAR(JP:)='*'
         Else
            Write(CPAR(JP:),*)IMPP(27,I1)
         End if
         JP=LBlnk(CPAR)+1
         CPAR(JP:JP)=';'
C         Call slpar(CPAR,JP)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if 

* - insert extension scores

      K1=0
      Do  61 I2=1,NABC
         If(IIPP(I2,I1).NE.IIPD(I2)) K1=1
 61   Continue 
      Do  62 I2=1,NABC
         If(IIPP(I2,I1).NE.IIPP( 1,I1)) K1=2
 62   Continue 

      If(K1.EQ.0) then
         JP=0 
      Else
         CPAR='I='
         JP=3
         If(K1.EQ.1) then
            If(IIPP( 1,I1).LE.NLOW) then 
               CPAR(JP:)='*'
            Else
               Write(CPAR(JP:),*)IIPP( 1,I1)
            End if 
            JP=Lblnk(CPAR)+1
            CPAR(JP:JP)=';'
         Else
            Do 64 I2=1,NABC
               If(I2.GT.1) then
                  JP=JP+2
                  CPAR(JP:JP)=','
               End if
               If(IMPP(I2,I1).LE.NLOW) then 
                  CPAR(JP+1:)='*'
               Else
                  Write(CPAR(JP+1:),*)IIPP(I2,I1)
               End if
               JP=Lblnk(CPAR)
 64         Continue
            JP=JP+1
            CPAR(JP:JP)=';'

C            Write(CPAR,'(''I='',I6,25('','',I6))')
C     *         (IIPP(ii1,I1),ii1=1,NABC)
C            JP=2+NABC*7
C            CPAR(JP:JP)=';'
C            J2=JP-NABC*7+1
C            Do  65 I2=1,NABC
C               If(IIPP(I2,I1).EQ.NLOW) CPAR(J2:J2+5)='     *'
C               J2=J2+7
C 65         Continue
         End if

         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if

      If(IIPP( 0,I1).NE.IIPD( 0)) then
         CPAR='I0='
         JP=4
         Write(CPAR(JP:),*)IIPP( 0,I1)
         JP=Lblnk(CPAR)+1
         CPAR(JP:JP)=';'
C         Call slpar(CPAR,JP)
         Call slpar(CPAR,JP)
         If(JB+JP+1.GT.LBUF) Go to 900
         CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
         JB=JB+JP+1
      End if 

* - initiation, transition and termination scores

      Do  67 I2=27,46
         If(IIPP(I2,I1).NE.IIPD(I2)) then
            If(IIPP(I2,I1).LE.NLOW) then
               Write(CPAR,*) CPSN(I2),'=*'
            Else
               Write(CPAR,*)CPSN(I2),'=',IIPP(I2,I1)
            End if 
            JP=Lblnk(CPAR)+1
            CPAR(JP:JP)=';'
            Call slpar(CPAR,JP)
            If(JB+JP+1.GT.LBUF) Go to 900
            CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
            JB=JB+JP+1
         End if 
 67   Continue  

      Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

* - restore previous defaults from pos. LPRF+1 

      I1=LPRF+1

      Do  68 I2=0,46
         IIPD(I2)=IIPP(I2,I1)
 68   Continue
      CHID=CHIP(I1)
      Do  69 I2=0,27
         IMPD(I2)=IMPP(I2,I1)
 69   Continue
      CHMD=CHMP(I1)

      Do  90 I1=0,LPRF

* write /M: block 

         If(I1.NE.0) then

* - symbol

            CBLK='/M:'
            JB=3

            If(CHMP(I1).NE.CHMD) then
               Write(CPAR,'(''SY='''''',A,'''''';'')') CHMP(I1)
               JP=7
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if 

* - match extension scores

            K1=0
            Do  71 I2=1,NABC
               If(IMPP(I2,I1).NE.IMPD(I2)) K1=1
 71         Continue 
            Do  72 I2=1,NABC
               If(IMPP(I2,I1).NE.IMPP( 1,I1)) K1=2
 72         Continue 

            If(K1.EQ.0) then
               JP=0 
            Else
               CPAR='M='
               JP=3
               If(K1.EQ.1) then
                  If(IMPP( 1,I1).LE.NLOW) then
                     CPAR(JP:)='*'
                  Else 
                     Write(CPAR(JP:),*)IMPP( 1,I1)
                  End if 
                  JP=Lblnk(CPAR)+1
                  CPAR(JP:JP)=';'
               Else
                  Do 73 I2=1,NABC
                     If(I2.GT.1) then
                        JP=JP+2
                        CPAR(JP:JP)=','
                     End if
                     If(IMPP(I2,I1).LE.NLOW) then 
                        CPAR(JP+1:)='*'
                     Else
                        Write(CPAR(JP+1:),*)IMPP(I2,I1)
                     End if
                     JP=Lblnk(CPAR)
 73               Continue
                  JP=JP+1
                  CPAR(JP:JP)=';'



C                  Write(CPAR,'(''M='',I6,25('','',I6))')
C     *               (IMPP(ii1,I1),ii1=1,NABC)
C                  JP=2+NABC*7
C                  CPAR(JP:JP)=';'
C                  J2=JP-NABC*7+1
C                  Do  75 I2=1,NABC
C                     If(IMPP(I2,I1).EQ.NLOW) CPAR(J2:J2+5)='     *'
C                     J2=J2+7
C 75               Continue
               End if

               Call slpar(CPAR,JP)
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if

* match extension score for unknown character

            If(IMPP( 0,I1).NE.IMPD( 0)) then
               CPAR='M0='
               JP=4
               If(IMPP( 0,I1).LE.NLOW) then
                  CPAR(JP:)='*'
               Else
                  Write(CPAR(JP:),*)IMPP( 0,I1)
               End if
               JP=LBlnk(CPAR)+1
               CPAR(JP:JP)=';'
C               Call slpar(CPAR,JP)
               Call slpar(CPAR,JP)
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if 

* - deletion extension score

            If(IMPP(27,I1).NE.IMPD(27)) then
               CPAR='D='
               JP=3
               If(IMPP(27,I1).EQ.NLOW) then
                  CPAR(JP:)='*'
               Else
                  Write(CPAR(JP:),*)IMPP(27,I1)
               End if
               JP=LBlnk(CPAR)+1
               CPAR(JP:JP)=';'
C               Call slpar(CPAR,JP)
               Call slpar(CPAR,JP)
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if 

            Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

         End if

* write /I: block 

         If(I1.NE.LPRF.OR..NOT.LPCI) then 

            CBLK='/I:'
            JB=3

* - symbol

            If(CHIP(I1).NE.CHID) then
               Write(CPAR,'(''SY='''''',A,'''''';'')') CHIP(I1)
               JP=7
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if 

* - insert extension scores

            K1=0
            Do  81 I2=1,NABC
               If(IIPP(I2,I1).NE.IIPD(I2)) K1=1
 81         Continue 
            Do  82 I2=1,NABC
               If(IIPP(I2,I1).NE.IIPP( 1,I1)) K1=2
 82         Continue 

            If(K1.EQ.0) then
               JP=0 
            Else
               CPAR='I='
               JP=3
               If(K1.EQ.1) then
                  If(IIPP( 1,I1).LE.NLOW) then
                     CPAR(JP:)='*'
                  Else 
                     Write(CPAR(JP:),*)IIPP( 1,I1)
                  End if 
                  JP=Lblnk(CPAR)+1
                  CPAR(JP:JP)=';'
               Else
                  Do 83 I2=1,NABC
                     If(I2.GT.1) then
                        JP=JP+2
                        CPAR(JP:JP)=','
                     End if
                     If(IIPP(I2,I1).LE.NLOW) then 
                        CPAR(JP+1:)='*'
                     Else
                        Write(CPAR(JP+1:),*)IIPP(I2,I1)
                     End if
                     JP=Lblnk(CPAR)
 83               Continue
                  JP=JP+1
                  CPAR(JP:JP)=';'



C                  Write(CPAR,'(''I='',I6,25('','',I6))')
C     *               (IIPP(ii1,I1),ii1=1,NABC)
C                  JP=2+NABC*7
C                  CPAR(JP:JP)=';'
C                  J2=JP-NABC*7+1
C                  Do  85 I2=1,NABC
C                     If(IIPP(I2,I1).EQ.NLOW) CPAR(J2:J2+5)='     *'
C                     J2=J2+7
C 85               Continue
               End if

               Call slpar(CPAR,JP)
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if

* insert extension score for unknown character

            If(IIPP( 0,I1).NE.IIPD( 0)) then
               CPAR='I0='
               JP=4
               If(IIPP( 0,I1).LE.NLOW) then
                  CPAR(JP:)='*'
               Else
                  Write(CPAR(JP:),*)IIPP( 0,I1)
               End if
               JP=LBlnk(CPAR)+1
               CPAR(JP:JP)=';'
C               Call slpar(CPAR,JP)
               Call slpar(CPAR,JP)
               If(JB+JP+1.GT.LBUF) Go to 900
               CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
               JB=JB+JP+1
            End if 

* - initiation, transition and termination scores

            Do  87 I2=27,46
               If(IIPP(I2,I1).NE.IIPD(I2)) then
                  If(IIPP(I2,I1).LE.NLOW) then
                     Write(CPAR,*) CPSN(I2),'=*'
                  Else
                     Write(CPAR,*)CPSN(I2),'=',IIPP(I2,I1)
                  End if 
                  JP=Lblnk(CPAR)+1
                  CPAR(JP:JP)=';'
                  Call slpar(CPAR,JP)
                  If(JB+JP+1.GT.LBUF) Go to 900
                  CBLK(JB+2:JB+JP+1)=CPAR(1:JP)
                  JB=JB+JP+1
               End if 
 87         Continue  

            If(CBLK.NE.'/I:') Call wrblk(NOUT,LLLT,LLEN,CBLK,JB)

         End if 

 90   Continue

* - Footer lines

      Do I1=1,LFTR
         Call wline(NOUT,LLLT,LLEN,CFTR(I1))
C         Write(NOUT,'(132A)')
C     *      (CFTR(I1)(ii1:ii1),ii1=1,Lblnk(CFTR(I1)))
      End do 

      Write(NOUT,'(''//'')')

 100  Return

 900  Write(NERR,*) 'Error: Line length exceeds buffer size (',
     *   LBUF,').'
      IRC=1
      Go to 100

      End
*----------------------------------------------------------------------*
      Subroutine slpar(CPAR,JP)

      Character*(*)     CPAR

      K1=0
      Do  10 I1=1,JP
         If(CPAR(I1:I1).NE.' '.OR.CPAR(I1-1:I1-1).EQ.';') then 
            K1=K1+1
            CPAR(K1:K1)=CPAR(I1:I1)
         Else if
     *         (CPAR(I1:I1+1).EQ.' ;'.AND.CPAR(K1:K1).EQ.'.') then
            K1=K1+1
            CPAR(K1:K1)='0' 
         Else if
     *         (CPAR(I1:I1+1).EQ.' ,'.AND.CPAR(K1:K1).EQ.'.') then
            K1=K1+1
            CPAR(K1:K1)='0' 
         End if
 10   Continue
      JP=K1
      
      Return
      End
*----------------------------------------------------------------------*
      Subroutine wrblk(NOUT,LLLT,LEN,CBLK,JB)

      Parameter        (LBUF=512)
      Character*(*)     CBLK
      Character*512     RCEX
      Logical           LLLT
      Integer           LEN
      Integer           NW

      RCEX='MA'


      If(LLLT) then
         NW=LEN
      Else
         NW=LBUF
      End if

      If(JB.LE.NW-5) then
         RCEX(6:)=CBLK(1:JB)
         Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,JB+5)
      Else
         Do  14 I1=NW-6,1,-1
            If(CBLK(I1:I1).EQ.';'
     *         .OR.CBLK(I1:I1).EQ.',') go to 15 
 14      Continue
 15      RCEX(6:)=CBLK(1:I1)
         Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,I1+5)
         RCEX='MA'
         KB=I1+1
         
 20      Continue
         
         Do  22 I1=KB,JB
            If(CBLK(I1:I1).NE.' ') go to 23 
 22      Continue
 23      KB=I1
         Do  24 I1=MIN(KB+NW-2,JB),KB,-1
            If(   CBLK(I1:I1).EQ.';'
     *         .OR.CBLK(I1:I1).EQ.',') go to  25
 24      Continue
 25      RCEX(9:)=CBLK(KB:I1)
         Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,Lblnk(RCEX))
         RCEX='MA'
         KB=I1+1
         If(KB.LT.JB) go to  20 
         
      End if
      
      Return
      End
*----------------------------------------------------------------------*
      Subroutine  wrnul(NOUT,LRNM,LLLT,LEN,NABC,FABC)

      Parameter        (LBUF=512)
      Real              FABC(0:26)
      Character*512     RCEX
      Logical           LRNM
      Integer           NMRS
      Logical           LLLT
      Integer           LEN
      Integer           NW

      If(LLLT) then
         NW=LEN
      Else
         NW=LBUF
      End if
      If(LRNM) then
         NMRS=100
      Else
         NMRS=1
      End if

      RCEX='MA      P='
      JP=11
      Do 5, I1=1,NABC
         Write(RCEX(JP:),*)NMRS*FABC(I1)
         JN=Lblnk(RCEX)
         If(RCEX(JN:JN).EQ.'.') then
            JN=JN+1
            RCEX(JN:JN)='0'
         End if
         If(I1.NE.NABC) then
            JN=JN+1
            RCEX(JN:JN)=','
         End if

         If(JN.GE.NW) then
            Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,JP)
            RCEX(1:10)='MA        '
 3          If(RCEX(JP+1:JP+1).EQ.' ') then
               JP=JP+1
               Go to 3
            End if
    
            RCEX(12:)=RCEX(JP:JN)
            JP=13+JN-JP
         Else
            JP=JN+1
         End if
 5    Continue
      
      If(JP.GT.10) then
         RCEX(JP:JP)=';'
         Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,JP)
      End if


 100  Return

      End
*----------------------------------------------------------------------*
      Subroutine  wline(NOUT,LLLT,LEN,RCEX)
      
      Character*(*)     RCEX
      Logical           LLLT
      Integer           LEN
      Integer           NW

      NW=Lblnk(RCEX)
      If(LLLT.AND.NW.GT.LEN) NW=LEN
      Write(NOUT,'(512A)')(RCEX(ii1:ii1),ii1=1,NW)
      
 100  Return
      End



