*----------------------------------------------------------------------*     
* $Id: xprsm.f,v 2.12 2003/12/09 13:42:42 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine XPRSM(JSEQ,JMMA,LPMM,
     *   LUNI,LOUT,LNOR,LREV,OPTZ,OPTL,OPLU,OPTB,
     *   CHID,CHAC,CHFH,CPID,CPAC,OPTR,OPTF,LDRS,
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
      Logical           OPTL
      Logical           OPLU
      Logical           OPTZ 
      Logical           OPTB 
      Logical           OPTR
      Logical           OPTF

      Logical           LPMM
      Logical           LFCL
      Logical           LDRS


* sequence and profile header 
      
      Character*(*)     CHID
      Character*(*)     CHAC
      Character*(*)     CHFH
      Character*(*)     CPID
      Character*(*)     CPAC

* Output fields 

      Character*13      CHLE
      Character*20      CHNS
      Character*17      CHRS 
      Character*64      CHLO
      Character*21      CHSE
      Character*17      CHTY
      
      Character*21      CHPP

      Character*24      CHPS
      Character*24      CHPE
      Character*09      CHST 
      Character*64      CHER 
      Character*64      CHMN 
      Character*22      CHRN 
      Character*26      CHPN 
      
      Character*64      CHMO
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
         CHNS='norm_score= '
         Write(CHNS(13:),'(F8.3)') R

         J1=11
         Do I1=12,Lblnk(CHNS)
            If(CHNS(I1:I1).NE.' ') then 
               J1=J1+1
               CHNS(J1:J1)=CHNS(I1:I1)
               CHNS(I1:I1)=' '
            End if
         End do
         ILNS=Lblnk(CHNS)

*         Write(NERR,*) 'CHNS: ', CHNS,' ILNS: ',ILNS

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
         If(OPLU) then
            If(.NOT.LFCL.OR.CCUT(KI)(1:1).EQ.' ') then
               Write(CHLE,'(''level_tag=NA'')')
            Else
               Write(CHLE,'(''level_tag='',(A2))') CCUT(KI)
            End if
         Else
            If(.NOT.LFCL) then
               Write(CHLE,'(''level=NA'')')
            Else
               If(K.LT.0.OR.K.GT.9) then 
                  Write(CHLE,'(''level='',I2,'' '')') K 
               Else
                  Write(CHLE,'(''level='',I1,''  '')') K 
               End if
            End if
         End if
         ILLE=Lblnk(CHLE)
*         Write(NERR,*) 'CHLE: ', CHLE,' ILLE: ',ILLE
      End if  
      
* - raw-score

      I=IOPT
      If(I.GT.999999) I=999999
      If(I.LE.-99999) I=-99999
      Write(CHRS,'(''raw_score='',I7)') I
      J1=10
      Do I1=11,Lblnk(CHRS)
         If(CHRS(I1:I1).NE.' ') then 
            J1=J1+1
            CHRS(J1:J1)=CHRS(I1:I1)
            CHRS(I1:I1)=' '
         End if
      End do
      ILRS=Lblnk(CHRS)
*      Write(NERR,*) 'CHRS: ', CHRS,' ILRS: ',ILRS
      
* - location in sequence

      If(.NOT.LUNI.OR.LOUT) then 
         If(LREV) then 
            JALB=LSEQ-JALB+1
            JALE=LSEQ-JALE+1
            CHSE='seq_end= '
            Write(CHSE(10:),*) JALB-LSEQ-1
         Else 
*         Write(CHLO,'('' pos. '',I6,'' -'',I6)') JALB,JALE
            CHSE='seq_end= '
            Write(CHSE(10:),*) JALE-LSEQ-1
         End if
         J1=8
         Do I1=9,Lblnk(CHSE)
            If(CHSE(I1:I1).NE.' ') then 
               J1=J1+1
               CHSE(J1:J1)=CHSE(I1:I1)
               CHSE(I1:I1)=' '
            End if
         End do
         ILSE=Lblnk(CHSE)

         Write(CHLO,*) JALB,'-',JALE
         J1=0
         Do I1=1,Lblnk(CHLO)
            If(CHLO(I1:I1).NE.' ') then 
               J1=J1+1
               CHLO(J1:J1)=CHLO(I1:I1)
               CHLO(I1:I1)=' '
            End if
         End do
         If(LREV) then 
            JALB=LSEQ-JALB+1
            JALE=LSEQ-JALE+1
         End if 
         ILLO=Lblnk(CHLO)
*         Write(NERR,*) 'CHLO: ', CHLO,' ILLO: ',ILLO

* - location in profile

         If(OPTZ) then
*            CHPS='motif_start= '
*            Write(CHPS(14:),*) IPMB
            Write(CHPS,'(''motif_start='',I7)') IPMB
            J1=12
            Do I1=13,Lblnk(CHPS)
               If(CHPS(I1:I1).NE.' ') then 
                  J1=J1+1
                  CHPS(J1:J1)=CHPS(I1:I1)
                  CHPS(I1:I1)=' '
               End if
            End do
*            CHPE='motif_end= '
*            Write(CHPE(12:),*) IPME
            Write(CHPE,'(''motif_end='',I7)') IPME
            J1=10
            Do I1=11,Lblnk(CHPE)
               If(CHPE(I1:I1).NE.' ') then 
                  J1=J1+1
                  CHPE(J1:J1)=CHPE(I1:I1)
                  CHPE(I1:I1)=' '
               End if
            End do
*            Write(CHPP,'('' ['',I5,'','',I6,'']'')') IPMB,IPME
            ILPS=Lblnk(CHPS)
*            Write(NERR,*) 'CHPS: ',CHPS,' ILPS: ',ILPS
            ILPE=Lblnk(CHPE)
*            Write(NERR,*) 'CHPE: ',CHPE,' ILPE: ',ILPE
         End if
      End if 

* - entry-ref 

      If(OPTF) then
         L1=Lblnk(CHFH)
         If (L1.EQ.0) then
            CHMI='unknown'
            ILMI=7
         Else
            CHMI=CHFH(1:L1)
            ILMI=MIN(64,L1)
         Endif
      Else
         L1=Lblnk(CHAC)
         L2=Lblnk(CHID)
         CHMI=CHAC(1:L1) // CHID(1:L2)
         ILMI=MIN(64,L1+L2)
C      CHER=CHAC(1:L1) // CHID(1:L2)
C      LNER=MIN(64,L1+L2)
      End if

      If(OPTB) then
         If(LREV) then 
*         CHER=CHER(1:LNER) // '(-)' 
*         LNER=LNER+3
            CHST='strand=r' 
         Else
            CHST='strand=s'
         End if 
         ILST=8
*        Write(NERR,*) 'CHST: ',CHST
      End if 
      
*      ILER=Lblnk(CHER)
*      Write(NERR,*) 'CHER: ',CHER,' ILER: ',ILER

* - match-number

      If(NALI.GT.1.AND..NOT.LPMM) then 
         CHMN='match_nb= '
         Write(CHMN(11:),*) JSEQ
         J1=9
         Do I1=10,Lblnk(CHMN)
            If(CHMN(I1:I1).NE.' ') then 
               J1=J1+1
               CHMN(J1:J1)=CHMN(I1:I1)
               CHMN(I1:I1)=' '
            End if
         End do
         ILMN=Lblnk(CHMN)
*         Write(NERR,*) 'CHMN: ',CHMN,' ILMN: ',ILMN
      End if
      If(.NOT.LPMM.AND.(JMMA.GT.1.OR.NALI.GT.1)) then
         CHTY='match_type=region'
         ILTY=17
      End if

* - multiple repeat number and type

      If(LPMM) then
         CHTY='match_type=repeat'
         ILTY=17
         CHRN='repeat_nb= '
         Write(CHRN(12:),*) JMMA
         J1=10
         Do I1=11,Lblnk(CHRN)
            If(CHRN(I1:I1).NE.' ') then 
               J1=J1+1
               CHRN(J1:J1)=CHRN(I1:I1)
               CHRN(I1:I1)=' '
            End if
         End do
         ILRN=Lblnk(CHRN)
*         Write(NERR,*) 'CHRN: ',CHRN,' ILRN: ',ILRN

* -- multiple repeat parent match number

         If(NALI.GT.1) then 
            CHPN='match_parent= '
            Write(CHPN(15:),*) JSEQ
            J1=13
            Do I1=14,Lblnk(CHPN)
               If(CHPN(I1:I1).NE.' ') then 
                  J1=J1+1
                  CHPN(J1:J1)=CHPN(I1:I1)
                  CHPN(I1:I1)=' '
               End if
            End do
            ILPN=Lblnk(CHPN)
C            Write(NERR,*) 'CHPN: ',CHPN,' ILPN: ',ILPN
         End if         
      End if

* - profile name (motif)

      CHMO='motif= '
      
      L1=Lblnk(CPID)
      L2=Lblnk(CPAC)
      L3=6
      If(L2-1.GT.0) then
         CHMO=CHMO(1:L3) // CPAC(1:L2)
         L3=L3+L2
      End if
      If(L1.GT.0) then
         CHMO=CHMO(1:L3) // CPID(1:L1)
      End if
      If(L2-1.LE.0.AND.L1.EQ.0) then
         CHMO=CHMO(1:6) // 'unknown'
      End if

      ILMO=Lblnk(CHMO)
*      Write(NERR,*) 'CHMO: ',CHMO,' ILMO: ',ILMO

* assemble ouput record

      LNEX=0 

      RCEX(LNEX+1:)='>' // CHMI 
      LNEX=LNEX+ILMI+2

      If(.NOT.LUNI.OR.LOUT) then 
         RCEX(LNEX:)='/' // CHLO
         LNEX=LNEX+ILLO+1
      End if

      RCEX(LNEX+1:)=CHMO 
      LNEX=LNEX+ILMO+1

      If(LNOR.AND..NOT.LDRS) then 
         RCEX(LNEX+1:)=CHNS
         LNEX=LNEX+ILNS+1
      End if

      RCEX(LNEX+1:)=CHRS
      LNEX=LNEX+ILRS+1

      If(OPTB) then 
         RCEX(LNEX+1:)=CHST
         LNEX=LNEX+ILST+1
      End if

      If(NALI.GT.1) then 
         If(.NOT.LPMM) then
            RCEX(LNEX+1:)=CHMN
            LNEX=LNEX+ILMN+1
         Else
            RCEX(LNEX+1:)=CHPN
            LNEX=LNEX+ILPN+1
         End if
      End if

      If(.NOT.LPMM.AND.(JMMA.GT.1.OR.NALI.GT.1)) then
         RCEX(LNEX+1:)=CHTY
         LNEX=LNEX+ILTY+1
      End if

      If(LPMM) then
         RCEX(LNEX+1:)=CHRN
         LNEX=LNEX+ILRN+1
         RCEX(LNEX+1:)=CHTY
         LNEX=LNEX+ILTY+1
      End if

      If(OPTL.OR.OPLU) then
         RCEX(LNEX+1:)=CHLE
         LNEX=LNEX+ILLE+1
*         If(OPLU) LNEX=LNEX-2
      End if

      If(.NOT.LUNI.OR.LOUT) then 
*         RCEX(LNEX+1:)=CHLO
*         LNEX=LNEX+ILLO
         RCEX(LNEX+1:)=CHSE
         LNEX=LNEX+ILSE+1
         
         If(OPTZ) then
            RCEX(LNEX+1:)=CHPS
            LNEX=LNEX+ILPS+1
            RCEX(LNEX+1:)=CHPE
            LNEX=LNEX+ILPE+1
         End if
      End if   

* - match ID and repeat number for standard output and option -x
      
*      If (.NOT.LPFA.AND.LPMM) then
*         Write(CHRP,*) JSEQ,'_REP'
*         L2=Lblnk(CHRP)
*         Write(CHRP(L2+1:),*) JMMA
*         J1=L2
*         Do I1=L2+1,Lblnk(CHRP)
*            If(CHRP(I1:I1).NE.' ') then 
*               J1=J1+1
*               CHRP(J1:J1)=CHRP(I1:I1)
*            End if
*         End do
*         RCEX(LNEX+1:)=CHRP(1:J1)
*         LNEX=LNEX+J1
*      End if

*      RCEX(LNEX+2:)=CHER(1:LNER)
*      LNEX=LNEX+1+LNER

*      L3=MAX(0,MIN(NW-LNEX-1,Lblnk(CHDE)))
*      RCEX(LNEX+2:)=CHDE(1:L3)
*      LNEX=LNEX+1+L3

* Write output record
      If(LNEX.GT.512) then
         LNEX=512
      End if
      Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,LNEX)

 100  Return 
      End 
