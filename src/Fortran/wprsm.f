*----------------------------------------------------------------------*     
* $Id: wprsm.f,v 2.13 2003/12/09 13:42:42 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine WPRSM(JSEQ,JMMA,LPMM,
     *   LUNI,LOUT,LNOR,LREV,LPFA,OPTZ,OPTL,OPLU,NW,
     *   CHID,CHAC,CHDE,OPTD,OPTR,LDRS,
     *   IOPT,JALB,JALE,NALI,IPMB,IPME,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   RNOP,KNPM,MAXN,INOR,IFUN,MNUM,LSEQ,RAVE)

* profile parameters

      Include          'codat.f'

      Real              RNOP(KNPM,MAXN)

* logical switches 

      Logical           LUNI
      Logical           LOUT
      Logical           LNOR
      Logical           LREV
      Logical           LPFA
      Logical           OPTL
      Logical           OPLU
      Logical           OPTZ 
      Logical           OPTD
      Logical           OPTR

      Logical           LPMM
      Logical           LFCL
      Logical           LDRS

* function return types

      Integer           Xblnk
      External          Xblnk

* sequence header 
      
      Character*(*)     CHID
      Character*(*)     CHAC
      Character*(*)     CHDE

* Output fields 

      Character*04      CHLE
      Character*08      CHNS
      Character*07      CHRS 
      Character*24      CHLO 
      Character*18      CHPP
      Character*64      CHER 
      Character*64      CHMI
      Character*64      CHRP

      Character*512     RCEX

* prepare output fields 

* - norm-score

      If(LNOR) then
         Call RtoN
     *      (IOPT,XOPT,RNOP,KNPM,MAXN,INOR,IFUN,LSEQ,RAVE)
         R=XOPT
         If(R.GT.9999.999) R=9999.999 
         If(R.LT.-999.999) R=-999.999
         Write(CHNS,'(F8.3)') R
      End if 

* - cut-off level

      If(OPTL.OR.OPLU) then
         LFCL=.FALSE.
         If(.NOT.OPTR) then
            KI=0
            KJ=0
            Do I1=1,JCUT
               Do I2=1,JCNM(I1)
                  If(MCUT(I2,I1).EQ.MNUM) then
                     If(KI.EQ.0) then
                        KI=I1
                        KJ=I2
                     Else if(RCUT(I2,I1).LT.RCUT(KJ,KI)) then
                        KI=I1
                        KJ=I2
                     End if
                  End if
               End do
            End do
            Do I1=1,JCUT
               Do I2=1,JCNM(I1)
                  If(MCUT(I2,I1).EQ.MNUM) then 
                     If(XOPT.GE.RCUT(I2,I1)
     *                  .AND.RCUT(I2,I1).GE.RCUT(KJ,KI)
     *                  .AND.MCLE(I1).GE.MCLE(KI)) then
                        LFCL=.TRUE.
                        KI=I1
                        KJ=I2
                     End if
                  End if
               End do
            End do
         Else
            KI=1
            Do I1=2,JCUT
               If(ICUT(I1).LT.ICUT(KI)) KI=I1
            End do 
            Do I1=1,JCUT
               If(IOPT.GE.ICUT(I1).AND.ICUT(I1).GE.ICUT(KI)
     *            .AND.MCLE(I1).GE.MCLE(KI)) then
                  LFCL=.TRUE.
                  KI=I1
               End if
            End do 
         End if 
         K=MCLE(KI)
         If(.NOT.LFCL) then
            Write(CHLE,'(''NA'')')
         Else if(OPLU) then
            CHLE=CCUT(KI)
         Else if(K.LT.0.OR.K.GT.9) then 
            Write(CHLE,'(''L='',I2)') K 
         Else
            Write(CHLE,'(''L='',I1,'' '')') K 
         End if
      End if  

* - raw-score

      I=IOPT
      If(I.GT.999999) I=999999
      If(I.LE.-99999) I=-99999
      Write(CHRS,'(I7)') I
      
* - location in sequence

      If(.NOT.LUNI.OR.LOUT) then 
         If(LREV) then 
            JALB=LSEQ-JALB+1
            JALE=LSEQ-JALE+1
         End if 
         Write(CHLO,'('' pos. '',I8,'' -'',I8)') JALB,JALE
         If(LREV) then 
            JALB=LSEQ-JALB+1
            JALE=LSEQ-JALE+1
         End if 

* - location in profile 

         If(OPTZ) then 
            Write(CHPP,'('' ['',I5,'','',I6,'']'')') IPMB,IPME
         End if
      End if 

* - entry-ref 

      L1=Lblnk(CHAC)
      L2=Lblnk(CHID)
      CHER=CHAC(1:L1) // CHID(1:L2)
      LNER=MIN(61,L1+L2)

      If(LREV) then 
         CHER=CHER(1:LNER) // '(-)' 
         LNER=LNER+3
      End if 

* - match-id

      If(LPFA.AND.NALI.GT.1) then 
         CHMI=CHID(1:L2) // '_'
         Write(CHMI(L2+2:),*) JSEQ
         J1=L2+1
         Do I1=L2+2,Lblnk(CHMI)
            If(CHMI(I1:I1).NE.' ') then 
               J1=J1+1
               CHMI(J1:J1)=CHMI(I1:I1)
            End if
         End do
         CHMI(J1+1:J1+1)=' '
         LNMI=J1
         L2=J1
      Else
         CHMI=CHID(1:L2)
         LNMI=L2
      End if

* - multiple repeat ID

*      If(LPFA.AND.LPMM) then 
*         CHMI=CHMI(1:L2) // '_REP'
*         Write(CHMI(L2+5:),*) JMMA
*         J1=L2+4
*         Do I1=L2+5,Lblnk(CHMI) 
*            If(CHMI(I1:I1).NE.' ') then 
*               J1=J1+1
*               CHMI(J1:J1)=CHMI(I1:I1)
*            End if
*         End do
*         CHMI(J1+1:J1+1)=' '
*         LNMI=J1
*      End if

* assemble ouput record

      LNEX=0 

      If(LPFA) then  
         RCEX(LNEX+1:)='>' // CHMI // ' ' 
         LNEX=LNEX+LNMI+2
      End if 

      If(OPTL.OR.OPLU) then
         RCEX(LNEX+1:)=CHLE
         LNEX=LNEX+4
         If(OPLU) LNEX=LNEX-2
      End if

      If(LNOR.AND..NOT.LDRS) then 
         RCEX(LNEX+1:)=CHNS
         LNEX=LNEX+8
      End if

      RCEX(LNEX+1:)=CHRS
      LNEX=LNEX+7

      If(.NOT.LUNI.OR.LOUT) then 
         RCEX(LNEX+1:)=CHLO
         LNEX=LNEX+24
         
         If(OPTZ) then
            RCEX(LNEX+1:)=CHPP
            LNEX=LNEX+15
         End if
      End if   

* - match ID, region nb and repeat nb for standard output and option -x
      
      If(.NOT.LPMM.AND.JMMA.GT.1) then
         Write(CHRP,*) 'REGION',JSEQ
         J1=6
         Do I1=7,Lblnk(CHRP)
            If(CHRP(I1:I1).NE.' ') then
               J1=J1+1
               CHRP(J1:J1)=CHRP(I1:I1)
            End if
         End do
         RCEX(LNEX+1:)=CHRP(1:J1)
         LNEX=LNEX+J1
      End if

      If(LPMM) then
         Write(CHRP,*) JSEQ,'REP'
         L2=Lblnk(CHRP)
         Write(CHRP(L2+1:),*) JMMA
         J1=L2
         Do I1=L2+1,Lblnk(CHRP)
            If(CHRP(I1:I1).NE.' ') then 
               J1=J1+1
               CHRP(J1:J1)=CHRP(I1:I1)
            End if
         End do
         RCEX(LNEX+1:)=CHRP(1:J1)
         LNEX=LNEX+J1
      End if

      RCEX(LNEX+2:)=CHER(1:LNER)
      LNEX=LNEX+1+LNER

      If(OPTD) then
         L3=MAX(0,MIN(NW-LNEX-1,Xblnk(CHDE,64)))
      Else
         L3=MAX(0,MIN(510-LNEX-1,Xblnk(CHDE,64)))
      End if
      RCEX(LNEX+2:)=CHDE(1:L3)
      LNEX=LNEX+1+L3

* Write output record
      Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,LNEX)

 100  Return 
      End 

