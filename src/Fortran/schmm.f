*----------------------------------------------------------------------*     
* $Id: schmm.f,v 2.3 2003/03/24 14:46:10 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine SCHMM 
     *   (NOUT,IDMP,RIHM,RMHM,LHMM,NABC,FLOW,FSCA,OPTF,OPFF,RD,RI)
      

      Real             FNOR(26)

      Include         'pfind.f'
      Include         'hmdat.f'
      Include         'sterr.f'

      Logical          OPTF
      Logical          OPFF

* rescale emission probabilities

      Do  20 I1=0,LHMM
         Do I2=1,NABC
            FNOR(I2)=RIHM(I2,I1)
         End do 
         Call GetNP(FNOR,NABC,RNOR,FLOW)
         If(RNOR.NE.0.0) then 
            Do I2=1,NABC
               If(RIHM(I2,I1).GT.FLOW) RIHM(I2,I1)=RIHM(I2,I1)+RNOR
            End do 
            If(RIHM(IM,I1).GT.FLOW) RIHM(IM,I1)=RIHM(IM,I1)-RNOR
            If(RIHM(II,I1).GT.FLOW) RIHM(II,I1)=RIHM(II,I1)-RNOR
            If(RIHM(ID,I1).GT.FLOW) RIHM(ID,I1)=RIHM(ID,I1)-RNOR
         End if 

         If(RIHM(II,I1).GE.0) then 
            Write(NOUT,'(''# Insert position '',I4,
     *         '': Log Prob(I->I)='',F8.4,'' reset to -0.01'')')
     *         I1,RIHM(II,I1)
            RIHM(II,I1)=-0.01
         End if 

         If(I1.EQ.0) then
            Do I2=1,NABC
               RMHM(I2,I1)=FLOW
            End do
            Go to  20 
         End if

         Do I2=1,NABC
            FNOR(I2)=RMHM(I2,I1)
         End do 
         Call GetNP(FNOR,NABC,RNOR,FLOW)
         If(RNOR.NE.0.0) then 
            Do I2=1,NABC
               If(RMHM(I2,I1).GT.FLOW) RMHM(I2,I1)=RMHM(I2,I1)+RNOR
            End do 
            If(RIHM(MM,I1).GT.FLOW) RIHM(MM,I1)=RIHM(MM,I1)-RNOR
            If(RIHM(MI,I1).GT.FLOW) RIHM(MI,I1)=RIHM(MI,I1)-RNOR
            If(RIHM(MD,I1).GT.FLOW) RIHM(MD,I1)=RIHM(MD,I1)-RNOR
         End if

         If(RMHM( D,I1).GT.FLOW) then
            RNOR=-RMHM( D,I1)
            RMHM( D,I1)=0.0
            RIHM(DM,I1)=RIHM(DM,I1)-RNOR
            RIHM(DI,I1)=RIHM(DI,I1)-RNOR
            RIHM(DD,I1)=RIHM(DD,I1)-RNOR
         End if 

 20   Continue

* transitions 

      Do  50 I1=LHMM,0,-1

* - insert positions

         FNOR( 1)=RIHM(IM,I1)
         FNOR( 2)=RIHM(ID,I1)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(II,  I1).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(II,I1)))
         If(RIHM(IM,I1).GT.FLOW) RIHM(IM,I1)=RIHM(IM,I1)+RNOR
         If(RIHM(ID,I1).GT.FLOW) RIHM(ID,I1)=RIHM(ID,I1)+RNOR
         If(RIHM(MI,I1).GT.FLOW) RIHM(MI,I1)=RIHM(MI,I1)-RNOR
         If(RIHM(DI,I1).GT.FLOW) RIHM(DI,I1)=RIHM(DI,I1)-RNOR

         If(I1.EQ.0) go to 50

* -match positions 

         FNOR( 1)=RIHM(MM,I1)
         FNOR( 2)=RIHM(MI,I1)
         FNOR( 3)=RIHM(MD,I1)
         Call GetNP(FNOR, 3,RNOR,FLOW)
         If(RIHM(MM,I1).GT.FLOW) RIHM(MM,I1)=RIHM(MM,I1)+RNOR
         If(RIHM(MI,I1).GT.FLOW) RIHM(MI,I1)=RIHM(MI,I1)+RNOR
         If(RIHM(MD,I1).GT.FLOW) RIHM(MD,I1)=RIHM(MD,I1)+RNOR
         J1=I1-1
         If(RIHM(MM,J1).GT.FLOW) RIHM(MM,J1)=RIHM(MM,J1)-RNOR
         If(RIHM(IM,J1).GT.FLOW) RIHM(IM,J1)=RIHM(IM,J1)-RNOR
         If(RIHM(DM,J1).GT.FLOW) RIHM(DM,J1)=RIHM(DM,J1)-RNOR

* - delete positions 

         FNOR( 1)=RIHM(DM,I1)
         FNOR( 2)=RIHM(DI,I1)
         FNOR( 3)=RIHM(DD,I1)
         Call GetNP(FNOR, 3,RNOR,FLOW)
         If(RIHM(DM,I1).GT.FLOW) RIHM(DM,I1)=RIHM(DM,I1)+RNOR
         If(RIHM(DI,I1).GT.FLOW) RIHM(DI,I1)=RIHM(DI,I1)+RNOR
         If(RIHM(DD,I1).GT.FLOW) RIHM(DD,I1)=RIHM(DD,I1)+RNOR
         J1=I1-1
         If(RIHM(MD,J1).GT.FLOW) RIHM(MD,J1)=RIHM(MD,J1)-RNOR
         If(RIHM(ID,J1).GT.FLOW) RIHM(ID,J1)=RIHM(ID,J1)-RNOR
         If(RIHM(DD,J1).GT.FLOW) RIHM(DD,J1)=RIHM(DD,J1)-RNOR

 50   Continue 

      FNOR( 1)=RIHM(MM, 0)
      FNOR( 2)=RIHM(MI, 0)
      FNOR( 3)=RIHM(MD, 0)
      Call GetNP(FNOR, 3,RNOR,FLOW)
      If(RIHM(MM, 0).GT.FLOW) RIHM(MM, 0)=RIHM(MM, 0)+RNOR
      If(RIHM(MI, 0).GT.FLOW) RIHM(MI, 0)=RIHM(MI, 0)+RNOR
      If(RIHM(MD, 0).GT.FLOW) RIHM(MD, 0)=RIHM(MD, 0)+RNOR

      FSCA=-RNOR

* add FIMs

      If(OPTF.OR.OPFF) then

         If(RD.GT.0.99999) RD=0.99999
         If(RI.GT.0.99999) RI=0.99999

         If(RD.GT.0.0) then
            RD=LOG(RD)
         Else
            RD=FLOW
         End if 
         If(RI.GT.0.0) then
            RI=LOG(RI)
         Else
            RI=FLOW
         End if 

* - state M(   0)

         RIHM(MD,   0)=RD
         FNOR( 1)=RIHM(MI,   0)
         FNOR( 2)=RIHM(MM,   0)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(MD,   0).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(MD,   0)))
         RIHM(MI,   0)=RIHM(MI,   0)+RNOR
         RIHM(MM,   0)=RIHM(MM,   0)+RNOR

         RIHM(MI,   0)=RD
         FNOR( 1)=RIHM(MD,   0)
         FNOR( 2)=RIHM(MM,   0)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(MI,   0).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(MI,   0)))
         RIHM(MD,   0)=RIHM(MD,   0)+RNOR
         RIHM(MM,   0)=RIHM(MM,   0)+RNOR

* - state I(   0)

         RIHM(ID,   0)=RD
         FNOR( 1)=RIHM(II,   0)
         FNOR( 2)=RIHM(IM,   0)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(ID,   0).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(ID,   0)))
         RIHM(II,   0)=RIHM(II,   0)+RNOR
         RIHM(IM,   0)=RIHM(IM,   0)+RNOR

         RIHM(II,   0)=RI
         FNOR( 1)=RIHM(ID,   0)
         FNOR( 2)=RIHM(IM,   0)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(II,   0).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(II,   0)))
         RIHM(ID,   0)=RIHM(ID,   0)+RNOR
         RIHM(IM,   0)=RIHM(IM,   0)+RNOR

* - state M(LHMM)

         RIHM(ID,LHMM)=RD
         FNOR( 1)=RIHM(II,LHMM)
         FNOR( 2)=RIHM(IM,LHMM)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(ID,LHMM).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(ID,LHMM)))
         RIHM(II,LHMM)=RIHM(II,LHMM)+RNOR
         RIHM(IM,LHMM)=RIHM(IM,LHMM)+RNOR

         RIHM(II,LHMM)=RI
         FNOR( 1)=RIHM(ID,LHMM)
         FNOR( 2)=RIHM(IM,LHMM)
         Call GetNP(FNOR, 2,RNOR,FLOW)
         If(RIHM(II,LHMM).GT.FLOW)
     *      RNOR=RNOR+LOG(1.0-EXP(RIHM(II,LHMM)))
         RIHM(ID,LHMM)=RIHM(ID,LHMM)+RNOR
         RIHM(IM,LHMM)=RIHM(IM,LHMM)+RNOR

* - state D(LHMM)

         RIHM(DI,LHMM)=RI
         If(RI.GT.FLOW) then
            RIHM(DM,LHMM)=1-EXP(RI) 
         Else
            RIHM(DM,LHMM)=1.0
         End if
         

* - state I(LHMM)

         RIHM(II,LHMM)=RI
         If(RI.GT.FLOW) then
            RIHM(IM,LHMM)=1-EXP(RI) 
         Else
            RIHM(IM,LHMM)=1.0
         End if

      End if

* - internal delete states

      If(OPFF.AND.RD.GT.FLOW) then
         RD=LOG(EXP(RD)/(1-EXP(RD)))
         Do I1=1,LHMM-1
            RIHM(DD,I1)=RD
            FNOR( 1)=RIHM(DD,I1)
            FNOR( 2)=RIHM(DI,I1)
            FNOR( 3)=RIHM(DM,I1)
            Call GetNP(FNOR, 3,RNOR,FLOW)
            RIHM(DD,I1)=RIHM(DD,I1)+RNOR
            RIHM(DI,I1)=RIHM(DI,I1)+RNOR
            RIHM(DM,I1)=RIHM(DM,I1)+RNOR
         End do
      End if 

 100  Return
      End
*----------------------------------------------------------------------*
      Subroutine GetNP(FNOR,INOR,RNOR,FLOW)

      Real              FNOR(*)

      RNOR=0.0

      RMAX=FLOW
      Do I1=1,INOR
         RMAX=MAX(RMAX,FNOR(I1))
      End do 

      If(RMAX.LE.FLOW) go to 100

      RSUM=0.0
      RCUT=MAX(FLOW,RMAX-20.0)
      Do I1=1,INOR
         If(FNOR(I1).GT.RCUT) RSUM=RSUM+EXP(FNOR(I1)-RMAX)
      End do 
      RNOR=-RMAX-LOG(RSUM) 

 100  Return
      End
*----------------------------------------------------------------------*
