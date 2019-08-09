*----------------------------------------------------------------------*     
* $Id: xali1.f,v 2.6 2003/01/15 14:49:05 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine XALI1
     *   (LPRF,LPCI,
     *   KCUT,IIPP,IMPP,IIPX,
     *   BSEQ,LSEQ,ISEQ,
     *   IOPM,IOPI,IOPD,
     *   IOPT,LOPT,  
     *   IRC)

* profile and sequence fields :

      Include           'ardim.f'
      Include           'gsdat.f'
      Include           'pfdat.f'
      Include           'pxdat.f'

* sequence

      Integer           LSEQ
      Integer           BSEQ
      Integer*2         ISEQ(IDMS)

* work fields

      Integer           IOPM(0:IDMP)
      Integer           IOPI(0:IDMP)
      Integer           IOPD(0:IDMP)

      Integer           KOPM

      Logical           LOPT

      IRC=0

      IOPT=NLOW

* beginning of sequence 

      IOPM(0)=IIPX(YM, 0)
      IOPI(0)=IIPX(YI, 0)
      IOPD(0)=IIPX(YD, 0)

      Do   8, I2=1,LPRF

         KD=IOPD(I2-1)+IMPP( D,I2)

         IOPM(I2)=MAX(KD+IIPP(DM,I2), 
     *      IIPX(YM,I2))  
         IOPI(I2)=MAX(KD+IIPP(DI,I2), 
     *      IIPX(YI,I2))  
         IOPD(I2)=MAX(KD+IIPP(DD,I2), 
     *      IIPX(YD,I2))  
 8    Continue

* - circular extensions:

      If(LPCI) then
         IOPM( 0)=MAX(IOPM( 0),IOPM(LPRF))
         IOPI( 0)=MAX(IOPI( 0),IOPI(LPRF))
         IOPD( 0)=MAX(IOPD( 0),IOPD(LPRF))

         Do   9, I2=1, LPRF
            KD=IOPD(I2-1)+IMPP( D,I2)
            If(IOPD(I2).GE.KD+IIPP(DD,I2)) then
               Go to 10
            Else  
               IOPM(I2)=MAX(KD+IIPP(DM,I2),IOPM(I2))
               IOPI(I2)=MAX(KD+IIPP(DI,I2),IOPI(I2))
               IOPD(I2)=    KD+IIPP(DD,I2)
            End if
 9       Continue    
 10      Continue
      End if 

* -----------------------------------------------------------

* internal sequence positions

      Do  50 I1=BSEQ,LSEQ-1

         J1=ISEQ(I1)
         KI=IOPI( 0)+IIPP(J1, 0)

         KOPM=IOPM( 0)

         IOPM( 0)=MAX(KI+IIPP(IM, 0),IIPX(XM, 0))  
         IOPI( 0)=MAX(KI+IIPP(II, 0),IIPX(XI, 0))  
         IOPD( 0)=MAX(KI+IIPP(ID, 0),IIPX(XD, 0))  
         IOPT    =MAX(KI+IIPX(IX, 0),IOPT) 

* - circular match extensions (presumably not necessary)

C             If(LPCI) then 
C                KM=IOPM(LPRF-1)+IMPP(J1,LPRF)
C                IOPM( 0)=MAX(IOPM( 0),KM+IIPP(MM, 0))                   
C                IOPI( 0)=MAX(IOPI( 0),KM+IIPP(MI, 0))                   
C                IOPD( 0)=MAX(IOPD( 0),KM+IIPP(MD, 0))                   
C                IOPT    =MAX(IOPT    ,KM+IIPX(MX, 0)) 
C             End if
* -----------------------------------------------------------

 20      Continue


         Do  38 I2=1,LPRF

            KM=KOPM      +IMPP(J1,I2)
            KI=IOPI(I2  )+IIPP(J1,I2)
            KD=IOPD(I2-1)+IMPP( D,I2)

            KOPM=IOPM(I2)

            IOPM(I2)=MAX(KM+IIPP(MM,I2), 
     *         KI+IIPP(IM,I2), 
     *         KD+IIPP(DM,I2), 
     *         IIPX(XM,I2))  
            IOPI(I2)=MAX(KM+IIPP(MI,I2), 
     *         KI+IIPP(II,I2), 
     *         KD+IIPP(DI,I2), 
     *         IIPX(XI,I2))  
            IOPD(I2)=MAX(KM+IIPP(MD,I2), 
     *         KI+IIPP(ID,I2), 
     *         KD+IIPP(DD,I2), 
     *         IIPX(XD,I2))  

            IOPT    =MAX(IOPT,
     *         KM+IIPX(MX,I2),
     *         KI+IIPX(IX,I2),
     *         KD+IIPX(DX,I2))

 38      Continue

* - circular extensions

         If(LPCI) then
            
            IOPM( 0)=MAX(IOPM( 0),IOPM(LPRF))
            IOPI( 0)=MAX(IOPI( 0),IOPI(LPRF))
            IOPD( 0)=MAX(IOPD( 0),IOPD(LPRF))
            
            Do  39, I2=1,LPRF
               KD=IOPD(I2-1)+IMPP( D,I2)
               If(IOPD(I2).GE.KD+IIPP(DD,I2)) then 
                  Go to  40
               Else
                  IOPM(I2)=MAX(KD+IIPP(DM,I2),IOPM(I2))
                  IOPI(I2)=MAX(KD+IIPP(DI,I2),IOPI(I2)) 
                  IOPD(I2)=MAX(KD+IIPP(DD,I2),IOPD(I2)) 
                  IOPT    =MAX(KD+IIPX(DX,I2),IOPT)
               End if
 39         Continue
 40         Continue
* -----------------------------------------------------------

         End if 

         If(.NOT.LOPT.AND.IOPT.GE.KCUT) go to 100

 50   Continue

* end of sequence 

      J1=ISEQ(LSEQ)

      KI=IOPI( 0)+IIPP(J1, 0)

      KOPM=IOPM( 0)

      IOPD( 0)=MAX(KI+IIPP(ID, 0),IIPX(XD, 0))  
      IOPT    =MAX(KI+IIPX(IY, 0),IOPT       ) 

* - circular match extension (presumably not necessary)

C             If(LPCI) then
C                KM=IOPM(LPRF-1)+IIPP(J1,LPRF)
C                IOPD( 0)=MAX(IOPD( 0),KM+IIPP(MD, 0))
C                IOPT    =MAX(IOPT,    KM+IIPP(MY, 0))
C             End if
* -----------------------------------------------------------

 60   Continue

      Do  68 I2=1,LPRF

         KM=KOPM      +IMPP(J1,I2)
         KI=IOPI(I2  )+IIPP(J1,I2)
         KD=IOPD(I2-1)+IMPP( D,I2)

         KOPM=IOPM(I2)

         IOPD(I2)=MAX(KM+IIPP(MD,I2), 
     *      KI+IIPP(ID,I2), 
     *      KD+IIPP(DD,I2), 
     *      IIPX(XD,I2))  

         IOPT    =MAX(IOPT,
     *      KM+IIPX(MY,I2),
     *      KI+IIPX(IY,I2),
     *      KD+IIPX(DY,I2))

 68   Continue

* - circular extensions

      If(LPCI) then

         IOPD( 0)=MAX(IOPD( 0),IOPD(LPRF))

         Do  69, I2=1,LPRF
            KD=IOPD(I2-1)+IMPP( D,I2)
            If(IOPD(I2).GE.KD+IIPP(DD,I2)) then 
               Go to 100 
            Else  
               IOPD(I2)=MAX(KD+IIPP(DD,I2),IOPD(I2))
               IOPT    =MAX(KD+IIPP(DY,I2),IOPT    )
            End if
 69      Continue

      End if 
* -----------------------------------------------------------

 100  Return

      End 
