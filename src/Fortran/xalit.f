*----------------------------------------------------------------------*     
* $Id: xalit.f,v 2.9 2003/04/11 13:43:53 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine XALIT
     *   (CABC,LPRF,LPCI,N1,N2,
     *   IIPP,IMPP,IIPX, 
     *   LSEQ,ISEQ,LCKS, 
     *   IOPM,IOPI,IOPD, 
     *   LALI,CALI,CPMA,IMSC,
     *   IOPT,JALB,JALE, 
     *   NMAT,
     *   PK2E,PK2B,PK3E,PK3B,
     *   PJ1E,PJ1B,PLAL,
     *   OPTM,
     *   IPMB,IPME,
     *   IRC) 

* profile and sequence fields :

      Include           'ardim.f'
      Include           'gsdat.f'
      Include           'pfdat.f'
      Include           'pxdat.f'
      Include           'sterr.f'

* sequence

      Integer           LSEQ
      Integer*2         ISEQ(IDMS)
      Logical           LCKS(IDMS)

* alignment
      
      Integer           LALI
      Character         CALI(IDMA)
      Character         CPMA(IDMM)
      Character         CH

* multiple matches of circular profile

      Integer           NMAT
C      Integer           IMMB(IDMN,IDML)
C      Integer           IMME(IDMN,IDML)
      Integer           PK2E(*)
      Integer           PK2B(*)
      Integer           PK3E(*)
      Integer           PK3B(*)
      Integer           PJ1E(*)
      Integer           PJ1B(*)
      Integer           PLAL(*)
      Integer           IMSC(*)
C      Integer           IMMS(512)

      Logical           OPTM

* prune delete positions at beginning and end of multiple matches
      Logical           LDFE
* Logical Delete To Begin
      Logical           LDTB
* Position (in alignment) Delete To Begin
      Integer           PDTB
* Score Delete to Begin
      Integer           SDTB
* Temp Delete to Begin score
      Integer           TDTB

* work fields

      Integer           TK1B
      Integer           TK2B
      Integer           TK3B
      Integer           TJ1B
      Integer           TLAL

      Integer           IOPM(0:IDMP)
      Integer           IOPI(0:IDMP)
      Integer           IOPD(0:IDMP)

      Integer           ZM
      Integer           ZI
      Integer           ZD

      Integer           MZ
      Integer           IZ
      Integer           DZ

      Integer           KOPM

      Logical           LIMM


      IRC=0
      MLOW=NLOW/4*3

C       If(NABC.LT.20) then 
C          CABC(0)='N'
C       Else 
C          CABC(0)='X'
C       End if    

      If((LPRF+1)*(JALE-JALB+1).GT.IDMM) then
         Write(NERR,*) ' Error: Alignment path matrix exceeds buffer ',
     *      'size (',IDMM,').'
         IRC=1
         Go to 100
      End if  
      If((JALE-JALB+1).GT.IDMA) then
         Write(NERR,*) ' Error: Alignment exceeds buffer ',
     *      'size (',IDMA,').'
         IRC=1
         Go to 100
      End if  

      K1=0

* beginning of sequence segment  

      If(JALB.EQ.1) then
         ZM=YM
         ZI=YI
         ZD=YD
      Else
         ZM=XM
         ZI=XI
         ZD=XD
      End if
      MZ=MX
      IZ=IX
      DZ=DX

      K1=K1+1
      I2=0

      KM=NLOW
      KI=NLOW
      KD=NLOW

      Call MaxP
     *   (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *   IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *   IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *   IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *   IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
      Call Ncode(CPMA(K1),JM,JI,JD)
      
      Do   9, I2=1,LPRF

         K1=K1+1

         KM=NLOW
         KI=NLOW 
         KD=IOPD(I2-1)+IMPP( D,I2)

         Call MaxP
     *      (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *      IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *      IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *      IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *      IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
         Call Ncode(CPMA(K1),JM,JI,JD)

 9    Continue

* - circular extensions

      If(LPCI) then
         L1=K1-LPRF
         Call Dcode(CPMA(L1),LM,LI,LD)
         If(IOPM( 0).LT.IOPM(LPRF)) then
            IOPM( 0)=IOPM(LPRF)
            LM=JM 
         End if 
         If(IOPI( 0).LT.IOPI(LPRF)) then
            IOPI( 0)=IOPI(LPRF)
            LI=JI 
         End if 
         If(IOPD( 0).LT.IOPD(LPRF)) then
            IOPD( 0)=IOPD(LPRF)
            LD=JD 
         End if 
         Call Ncode(CPMA(L1),LM,LI,LD)

         Do I2=1,LPRF
            L1=L1+1
            KD=IOPD(I2-1)+IMPP( D,I2)
            If(IOPD(I2).GE.KD+IIPP(DD,I2)) then
               Go to  10
            Else
               Call Dcode(CPMA(L1),LM,LI,LD)
               If(IOPM(I2).LT.KD+IIPP(DM,I2)) then
                  IOPM(I2)=KD+IIPP(DM,I2)
                  LM=3
               End if 
               If(IOPI(I2).LT.KD+IIPP(DI,I2)) then
                  IOPI(I2)=KD+IIPP(DI,I2)
                  LI=3
               End if 
               IOPD(I2)=KD+IIPP(DD,I2)
               LD=3
               Call Ncode(CPMA(L1),LM,LI,LD)
            End if
         End do
 10      Continue
      End if
C          Call wtrace(IDMP,IDMM,LPRF,CPMA,K1-LPRF,IOPM,IOPI,IOPD)
* -----------------------------------------------------------


* internal sequence positions

      ZM=XM
      ZI=XI
      ZD=XD

      Do  25 I1=JALB,JALE-1

         If(LCKS(I1)) then
            IOPM(N1-1)=NLOW
            Do I2=N1,N2-1
               IOPM(I2)=NLOW
               IOPI(I2)=NLOW
            End do 
         End if

         J1=ISEQ(I1)
         K1=K1+1
         I2=0

         KM=NLOW
         KI=IOPI( 0)+IIPP(J1, 0)
         KD=NLOW

         KOPM=IOPM( 0)

         Call MaxP
     *      (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *      IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *      IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *      IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *      IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
         Call Ncode(CPMA(K1),JM,JI,JD)

         Do  19 I2=1,LPRF

            K1=K1+1

            KM=KOPM      +IMPP(J1,I2)
            KI=IOPI(I2  )+IIPP(J1,I2)
            KD=IOPD(I2-1)+IMPP( D,I2)

            KOPM=IOPM(I2)

            Call MaxP
     *         (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *         IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *         IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *         IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *         IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
            Call Ncode(CPMA(K1),JM,JI,JD)
 19      Continue

* - circular extensions

         If(LPCI) then
            L1=K1-LPRF
            Call Dcode(CPMA(L1),LM,LI,LD)
            If(IOPM( 0).LT.IOPM(LPRF)) then
               IOPM( 0)=IOPM(LPRF)
               LM=JM 
            End if 
            If(IOPI( 0).LT.IOPI(LPRF)) then
               IOPI( 0)=IOPI(LPRF)
               LI=JI 
            End if 
            If(IOPD( 0).LT.IOPD(LPRF)) then
               IOPD( 0)=IOPD(LPRF)
               LD=JD 
            End if 
            Call Ncode(CPMA(L1),LM,LI,LD)

            Do I2=1,LPRF
               L1=L1+1
               KD=IOPD(I2-1)+IMPP( D,I2)
               If(IOPD(I2).GE.KD+IIPP(DD,I2)) then
                  Go to  20
               Else
                  Call Dcode(CPMA(L1),LM,LI,LD)
                  If(IOPM(I2).LT.KD+IIPP(DM,I2)) then
                     IOPM(I2)=KD+IIPP(DM,I2)
                     LM=3
                  End if 
                  If(IOPI(I2).LT.KD+IIPP(DI,I2)) then
                     IOPI(I2)=KD+IIPP(DI,I2)
                     LI=3
                  End if 
                  IOPD(I2)=KD+IIPP(DD,I2)
                  LD=3
                  Call Ncode(CPMA(L1),LM,LI,LD)
               End if
            End do
 20         Continue
         End if
* -----------------------------------------------------------
C          Call wtrace(IDMP,IDMM,LPRF,CPMA,K1-LPRF,IOPM,IOPI,IOPD)

 25   Continue

* end of sequence 

      If(JALE.EQ.LSEQ) then
         MZ=MY
         IZ=IY
         DZ=DY
      End if  

      If(LCKS(JALE)) then
         IOPM(N1-1)=NLOW
         Do I2=N1,N2-1
            IOPM(I2)=NLOW
            IOPI(I2)=NLOW
         End do 
      End if

      J1=ISEQ(JALE)
      K1=K1+1
      I2=0

      KM=NLOW
      KI=IOPI( 0)+IIPP(J1, 0)
      KD=NLOW 

      KOPM=IOPM( 0)

      Call MaxP
     *   (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *   IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *   IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *   IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *   IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
      Call Ncode(CPMA(K1),JM,JI,JD)

      I2=0
      If(KI+IIPX(IZ,I2).GE.IOPT) then
         JS=1
         Go to  50
      End if

      Do  29 I2=1,LPRF

         K1=K1+1

         KM=KOPM      +IMPP(J1,I2)
         KI=IOPI(I2  )+IIPP(J1,I2)
         KD=IOPD(I2-1)+IMPP( D,I2)

         KOPM=IOPM(I2)

         Call MaxP
     *      (IOPM(I2),IOPI(I2),IOPD(I2),JM,JI,JD,KM,KI,KD,
     *      IIPP(MM,I2),IIPP(MI,I2),IIPP(MD,I2),
     *      IIPP(IM,I2),IIPP(II,I2),IIPP(ID,I2),
     *      IIPP(DM,I2),IIPP(DI,I2),IIPP(DD,I2),
     *      IIPX(ZM,I2),IIPX(ZI,I2),IIPX(ZD,I2))
         Call Ncode(CPMA(K1),JM,JI,JD)

         If     (KI+IIPX(IZ,I2).GE.IOPT) then 
            JS=1
            Go to  50
         Else if(KM+IIPX(MZ,I2).GE.IOPT) then 
            JS=2
            Go to  50
         Else if(KD+IIPX(DZ,I2).GE.IOPT) then 
            JS=3
            Go to  50
         End if 

 29   Continue

* - circular extensions

      If(LPCI) then
         L1=K1-LPRF
         Call Dcode(CPMA(L1),LM,LI,LD)
         If(IOPD( 0).LT.IOPD(LPRF)) then
            IOPD( 0)=IOPD(LPRF)
            LD=JD 
         End if 
         Call Ncode(CPMA(L1),LM,LI,LD)

         Do I2=1,LPRF
            L1=L1+1
            KD=IOPD(I2-1)+IMPP( D,I2)
            If(KD+IIPX(DZ,I2).GE.IOPT) then
               JS=3
               K1=L1
            End if 
            If(IOPD(I2).GE.KD+IIPP(DD,I2)) then
               Go to  30
            Else
               Call Dcode(CPMA(L1),LM,LI,LD)
               IOPD(I2)=KD+IIPP(DD,I2)
               LD=3
               Call Ncode(CPMA(L1),LM,LI,LD)
            End if
         End do
 30      Continue
      End if
* -----------------------------------------------------------

 50   Continue

C          Call wtrace(IDMP,IDMM,LPRF,CPMA,K1-LPRF,IOPM,IOPI,IOPD)

      K2=JALE
      K3=I2
      IPME=K3-LPRF-1
C      Write(NERR,*) 'K1: ',K1,' K2: ',K2,' K3: ',K3,
C     *   ' IPME: ',IPME,' OPTM: ',OPTM


* trace back
      
      JO=0
      K4=0
      K5=IOPT
      J1=0
      Do  I1=LPRF,K3+1,-1
         J1=J1+1
         CALI(J1)='-'
      End do

      LIMM=.FALSE.
      SDTB=MLOW
      TDTB=MLOW
      LDTB=.FALSE.
 60   Continue
C       Write(6,'(5I5)') K1,K2,K3,JS

* insert position

C      Write(NERR,*)'JO: ',JO,' JS: ',JS

      If     (JS.EQ.1) then 
         LDTB=.FALSE.
         If(LIMM.AND.OPTM) then
            LDFE=.FALSE.
            If(JO.EQ.1) then
               IMSC(K4)=IMSC(K4)+IIPP(II,K3)
C               Write(NERR,*)'Score: ',IIPP(II,K3)
            Else if(JO.EQ.2) then
               IMSC(K4)=IMSC(K4)+IIPP(MI,K3)
C               Write(NERR,*)'Score: ',IIPP(IM,K3)
            Else if(JO.EQ.3) then
               IMSC(K4)=IMSC(K4)+IIPP(DI,K3)
C               Write(NERR,*)'Score: ',IIPP(ID,K3)
            End if
            IMSC(K4)=IMSC(K4)+IIPP(ISEQ(K2),K3)
C            Write(NERR,*)'Score: ',IIPP(ISEQ(K2),K3)
         End if
         J1=J1+1
         CALI(J1)=Char(Ichar(CABC(ISEQ(K2)))+32)
C         Write(NERR,*) 'CALI(',J1,'): ',CALI(J1),' No Match',
C     *      ' JS: ',JS
         K1=K1-LPRF-1
         K2=K2-1

* match position

      Else if(JS.EQ.2) then
         LDTB=.FALSE.
         If(LIMM.AND.OPTM) then
            If(JO.EQ.1) then
               IMSC(K4)=IMSC(K4)+IMPP(ISEQ(K2),K3)
C               Write(NERR,*)'Score: ',IMPP(ISEQ(K2),K3)
               IMSC(K4)=IMSC(K4)+IIPP(MI,K3)
C               Write(NERR,*)'Score: ',IIPP(MI,K3)
            Else if(JO.EQ.2) then
               LDFE=.FALSE.
               IMSC(K4)=IMSC(K4)+IMPP(ISEQ(K2),K3)
C               Write(NERR,*)'Score: ',IMPP(ISEQ(K2),K3)
               IMSC(K4)=IMSC(K4)+IIPP(MM,K3)
C               Write(NERR,*)'Score: ',IIPP(MM,K3)
            Else if(JO.EQ.3) then
               If((LDFE.EQV..TRUE.).AND.
     *            (IIPX(MX,K3).GE.(IMSC(K4)+IIPP(MD,K3)))) then
                  PK2E(K4)=K2
                  PK3E(K4)=K3
                  PJ1E(K4)=J1
                  IMSC(K4)=IIPX(MX,K3)+IMPP(ISEQ(K2),K3)
               Else
                  IMSC(K4)=IMSC(K4)+IMPP(ISEQ(K2),K3)
C               Write(NERR,*)'Score: ',IMPP(ISEQ(K2),K3)
                  IMSC(K4)=IMSC(K4)+IIPP(MD,K3)
C               Write(NERR,*)'Score: ',IIPP(MD,K3)
               End if
            End if
            LDFE=.FALSE.
         Else if(.NOT.LIMM.AND.OPTM) then
            K4=K4+1
            PK2E(K4)=K2
            PK3E(K4)=K3
            PJ1E(K4)=J1
            If(K2.GE.LSEQ.AND.IIPX(MY,K3).GT.MLOW) then
               IMSC(K4)=IIPX(MY,K3)
C               Write(NERR,*)'Score: ',IIPX(MY,K3)
            Else if(IIPX(MX,K3).GT.MLOW) then
               IMSC(K4)=IIPX(MX,K3)
C               Write(NERR,*)'Score: ',IIPX(MX,K3)
            Else
               IMSC(K4)=0
C               Write(NERR,*)'Score: 0'
            End if
            IMSC(K4)=IMSC(K4)+IMPP(ISEQ(K2),K3)
C            Write(NERR,*)'Score: ',IMPP(ISEQ(K2),K3)
C            Write(NERR,*) 'End mult match nb: ',K4,' in match'
C            Write(NERR,*) 'PK2E: ',PK2E(K4)
C            Write(NERR,*) 'PK3E: ',PK3E(K4)
C            Write(NERR,*) 'PJ1E: ',PJ1E(K4)
            LIMM=.TRUE.
            LDFE=.FALSE.
            SDTB=MLOW
            TDTB=MLOW
         End if
         J1=J1+1
         CALI(J1)=CABC(ISEQ(K2))
C         Write(NERR,*) 'CALI(',J1,'): ',CALI(J1),' Match',
C     *      ' JS: ',JS
         K1=K1-LPRF-2
         K2=K2-1
         K3=K3-1  
         If(LPCI.AND.K3.EQ.0) then 
            If(LIMM.AND.OPTM) then
               PK2B(K4)=K2+1
               PK3B(K4)=K3+1
               PJ1B(K4)=J1
               PLAL(K4)=PJ1B(K4)-PJ1E(K4)
               If(K2.LE.1.AND.IIPX(YM,K3).GT.MLOW) then
                  IMSC(K4)=IMSC(K4)+IIPX(YM,K3)
C                  Write(NERR,*)'Score: ',IIPX(YM,K3)
               Else if(IIPX(XM,K3).GT.MLOW) then
                  IMSC(K4)=IMSC(K4)+IIPX(XM,K3)
C                  Write(NERR,*)'Score: ',IIPX(XM,K3)
*              Else
*                 Do nothing
               End if
C               Write(NERR,*) 'Begin mult match nb: ',K4,' in match'
C               Write(NERR,*) 'PK2B: ',PK2B(K4)
C               Write(NERR,*) 'PK3B: ',PK3B(K4)
C               Write(NERR,*) 'PJ1B: ',PJ1B(K4)
C               Write(NERR,*) 'PLAL: ',PLAL(K4)
               LIMM=.FALSE.
            End if
            K3=K3+LPRF
            K1=K1+LPRF
         End if

* delete position

      Else if(JS.EQ.3) then 
         LDTB=.TRUE.
         If(LIMM.AND.OPTM) then
            If(JO.EQ.1) then
               TDTB=IMSC(K4)+IIPX(XI,K3)
               IMSC(K4)=IMSC(K4)+IMPP(D,K3)
C               Write(NERR,*)'Score: ',IMPP(D,K3)
               IMSC(K4)=IMSC(K4)+IIPP(DI,K3)
C               Write(NERR,*)'Score: ',IIPP(DI,K3)
            Else if(JO.EQ.2) then
               TDTB=IMSC(K4)+IIPX(XM,K3)
               IMSC(K4)=IMSC(K4)+IMPP(D,K3)
C               Write(NERR,*)'Score: ',IMPP(D,K3)
               IMSC(K4)=IMSC(K4)+IIPP(DM,K3)
C               Write(NERR,*)'Score: ',IIPP(DM,K3)
            Else if(JO.EQ.3) then
               TDTB=IMSC(K4)+IIPX(XD,K3)
               IMSC(K4)=IMSC(K4)+IMPP(D,K3)
C               Write(NERR,*)'Score: ',IMPP(D,K3)
               IMSC(K4)=IMSC(K4)+IIPP(DD,K3)
C               Write(NERR,*)'Score: ',IIPP(DD,K3)
            End if
            If(TDTB.LT.MLOW) TDTB=MLOW
            If(TDTB.GT.SDTB) then
               SDTB=TDTB
               TK1B=K1
               TK2B=K2+1
               TK3B=K3
               TJ1B=J1
               TLAL=TJ1B-PJ1E(K4)
            End if
         Else if(.NOT.LIMM.AND.OPTM) then
            K4=K4+1
            PK2E(K4)=K2
            PK3E(K4)=K3
            PJ1E(K4)=J1
            If(K2.GE.LSEQ.AND.IIPX(DY,K3).GT.MLOW) then
               IMSC(K4)=IIPX(DY,K3)
C               Write(NERR,*)'Score: ',IIPX(DY,K3)
            Else if(IIPX(DX,K3).GT.MLOW) then
               IMSC(K4)=IIPX(DX,K3)
C               Write(NERR,*)'Score: ',IIPX(DX,K3)
            Else
               IMSC(K4)=0
            End if
            IMSC(K4)=IMSC(K4)+IMPP(D,K3)
C            Write(NERR,*)'Score: ',IMPP(D,K3)
C            Write(NERR,*) 'End mult match nb: ',K4,' in del'
C            Write(NERR,*) 'PK2E: ',PK2E(K4)
C            Write(NERR,*) 'PK3E: ',PK3E(K4)
C            Write(NERR,*) 'PJ1E: ',PJ1E(K4)
            LIMM=.TRUE.
            LDFE=.TRUE.
            SDTB=NLOW
            TDTB=NLOW
         End if
         J1=J1+1
         CALI(J1)='-'
C         Write(NERR,*) 'CALI(',J1,'): ',CALI(J1),' Deletion',
C     *      ' JS: ',JS
         K1=K1-1
         K3=K3-1
         If(LPCI.AND.K3.EQ.0) then 
            If(LIMM.AND.OPTM) then
               If(K2.LE.1.AND.IIPX(YD,K3).GT.MLOW) then
                  IMSC(K4)=IMSC(K4)+IIPX(YD,K3)
C                  Write(NERR,*)'Score: ',IIPX(YD,K3)
               Else if(IIPX(XD,K3).GT.MLOW) then
                  IMSC(K4)=IMSC(K4)+IIPX(XD,K3)
C                  Write(NERR,*)'Score: ',IIPX(XD,K3)
*              Else
*                 Do nothing
               End if
               If(.NOT.LDTB.OR.IMSC(K4).LE.-99999
     *            .OR.SDTB.LT.IMSC(K4)) then
                  PK2B(K4)=K2+1
                  PK3B(K4)=K3+1
                  PJ1B(K4)=J1
                  PLAL(K4)=PJ1B(K4)-PJ1E(K4)
C                  Write(NERR,*) 'Begin mult match nb: ',K4,' in del'
C                  Write(NERR,*) 'LDTB: ',LDTB
C                  Write(NERR,*) 'IMSC: ',IMSC(K4)
C                  Write(NERR,*) 'PK2B: ',PK2B(K4)
C                  Write(NERR,*) 'PK3B: ',PK3B(K4)
C                  Write(NERR,*) 'PJ1B: ',PJ1B(K4)
C                  Write(NERR,*) 'PLAL: ',PLAL(K4)
               Else
                  IMSC(K4)=SDTB
                  PK2B(K4)=TK2B
                  PK3B(K4)=TK3B+1
                  PJ1B(K4)=TJ1B
                  PLAL(K4)=TLAL
C                  Write(NERR,*) 'Begin mult match nb: ',K4,' in del_t_b'
C                  Write(NERR,*) 'LDTB: ',LDTB
C                  Write(NERR,*) 'IMSC: ',IMSC(K4)
C                  Write(NERR,*) 'PK2B: ',PK2B(K4)
C                  Write(NERR,*) 'PK3B: ',PK3B(K4)
C                  Write(NERR,*) 'PJ1B: ',PJ1B(K4)
C                  Write(NERR,*) 'PLAL: ',PLAL(K4)
               End if
               LIMM=.FALSE.
            End if
            K3=K3+LPRF
            K1=K1+LPRF
         End if
      End if 

* next traceback step

      Call Dcode(CPMA(K1),JM,JI,JD)
C      Write(NERR,*) 'K1: ',K1,' JM: ',JM,' JI: ',JI,
C     *   ' JD: ',JD,' JS: ',JS

      JO=JS

      If     (JS.EQ.1) then 
         JS=JI
      Else if(JS.EQ.2) then
         JS=JM
      Else if(JS.EQ.3) then
         JS=JD
      End if

      If(JS.NE.0) Go to 60

C      Write(NERR,*) 'JS eq 0 for: ',K2
      
      If(LIMM.AND.OPTM) then
         If(K2.LT.1) then
            If(JO.EQ.1.AND.IIPX(YI,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(YI,K3)
C               Write(NERR,*)'Score: ',IIPX(YI,K3)
            Else if(JO.EQ.2.AND.IIPX(YM,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(YM,K3)
C               Write(NERR,*)'Score: ',IIPX(YM,K3)
            Else if(JO.EQ.3.AND.IIPX(YD,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(YD,K3)
C               Write(NERR,*)'Score: ',IIPX(YD,K3)
            End if
         Else
            If(JO.EQ.1.AND.IIPX(XI,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(XI,K3)
C               Write(NERR,*)'Score: ',IIPX(XI,K3)
            Else if(JO.EQ.2.AND.IIPX(XM,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(XM,K3)
C               Write(NERR,*)'Score: ',IIPX(XM,K3)
            Else if(JO.EQ.3.AND.IIPX(XD,K3).GT.MLOW) then
               IMSC(K4)=IMSC(K4)+IIPX(XD,K3)
C               Write(NERR,*)'Score: ',IIPX(XD,K3)
            End if
         End if
         If(.NOT.LDTB.OR.IMSC(K4).LE.-99999
     *      .OR.SDTB.LT.IMSC(K4)) then 
            PK2B(K4)=K2+1
            PK3B(K4)=K3+1
            PJ1B(K4)=J1
            PLAL(K4)=PJ1B(K4)-PJ1E(K4)
         Else
            IMSC(K4)=SDTB
            PK2B(K4)=TK2B
            PK3B(K4)=TK3B+1
            PJ1B(K4)=TJ1B
            PLAL(K4)=TLAL
         End if
C         Write(NERR,*) 'Begin mult match: ',K4
C         Write(NERR,*) 'PK2B: ',PK2B(K4)
C         Write(NERR,*) 'PK3B: ',PK3B(K4)
C         Write(NERR,*) 'PJ1B: ',PJ1B(K4)
C         Write(NERR,*) 'PLAL: ',PLAL(K4)
      End if

      If(K3.GE.LPRF) K3=0

      Do I1=K3,1,-1
C         Write(NERR,*) 'Add deletions'
         J1=J1+1
         CALI(J1)='-'
      End do

      LALI=J1
      If(OPTM) then
         NMAT=K4
C         Write(NERR,*)'*** Scores:'
C         Write(NERR,*)(IMSC(ii1),ii1=1,NMAT)
      End if

* reverse begin and end positions in alignment

      If(OPTM) then
         Do I1=1,NMAT
            PJ1B(I1)=LALI-PJ1B(I1)+1
            PJ1E(I1)=LALI-PJ1E(I1)
         End do
      End if

* reverse alignment

      J1=(LALI+1)/2
      Do I1=LALI/2+1,LALI 
         CH=CALI(I1)
         CALI(I1)=CALI(J1)
         CALI(J1)=CH         
         J1=J1-1
      End do

      IPMB=K3+1

 100  Return

      End 
*----------------------------------------------------------------------*
      Subroutine MaxP
     *   (IOPM,IOPI,IOPD,JM,JI,JD,KM,KI,KD,
     *   IMM,IMI,IMD,IIM,III,IID,
     *   IDM,IDI,IDD,IZM,IZI,IZD)

      JMM=KM+IMM
      JMI=KM+IMI
      JMD=KM+IMD
      JIM=KI+IIM
      JII=KI+III
      JID=KI+IID
      JDM=KD+IDM
      JDI=KD+IDI
      JDD=KD+IDD

      IOPM=IZM
      JM=0
      IOPI=IZI
      JI=0
      IOPD=IZD
      JD=0

      If(JIM.GT.IOPM) then
         IOPM=JIM
         JM=1
      End if
      If(JMM.GT.IOPM) then
         IOPM=JMM
         JM=2
      End if
      If(JDM.GT.IOPM) then
         IOPM=JDM
         JM=3
      End if

      If(JII.GT.IOPI) then
         IOPI=JII
         JI=1
      End if
      If(JMI.GT.IOPI) then
         IOPI=JMI
         JI=2
      End if
      If(JDI.GT.IOPI) then
         IOPI=JDI
         JI=3
      End if

      If(JID.GT.IOPD) then
         IOPD=JID
         JD=1
      End if
      If(JMD.GT.IOPD) then
         IOPD=JMD
         JD=2
      End if
      If(JDD.GT.IOPD) then
         IOPD=JDD
         JD=3
      End if

      Return
      End
*----------------------------------------------------------------------*
      Subroutine Ncode(CPMA,JM,JI,JD)

      Character         CPMA

      CPMA=Char(JM*16+JI*4+JD) 
      Return 
      End 
*----------------------------------------------------------------------*
      Subroutine Dcode(CPMA,JM,JI,JD)

      Character         CPMA

      JD=Ichar(CPMA)
      JM=JD/16
      JD=JD-JM*16
      JI=JD/ 4
      JD=JD-JI* 4
      Return 
      End 
*----------------------------------------------------------------------*
      Subroutine wtrace(IDMP,IDMM,LPRF,CPMA,K1,IOPM,IOPI,IOPD)

      Character CPMA(IDMM)
      Integer   IOPM(0:IDMP)
      Integer   IOPI(0:IDMP)
      Integer   IOPD(0:IDMP)

      Integer   IM(20)
      Integer   II(20)
      Integer   ID(20)

      J1=1
      Do I1=K1,K1+LPRF
         Call Dcode(CPMA(I1),IM(J1),II(J1),ID(J1))
         J1=J1+1
      End do

      Write(6,'(I4)') K1/(LPRF+1)
      Write(6,'(20(I4,I2))')(IOPM(ii1-1),IM(ii1),ii1=1,LPRF+1)
      Write(6,'(20(I4,I2))')(IOPI(ii1-1),II(ii1),ii1=1,LPRF+1)
      Write(6,'(20(I4,I2))')(IOPD(ii1-1),ID(ii1),ii1=1,LPRF+1)
      Return
      End

