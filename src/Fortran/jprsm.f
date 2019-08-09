*----------------------------------------------------------------------*     
* $Id: xprsm.f,v 2.12 2003/12/09 13:42:42 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine JPRSM(CHID,CHAC,CHFH,OPTF,JALB,JALE)

* logical switches 

      Logical           OPTF

* sequence and profile header 
      
      Character*(*)     CHID
      Character*(*)     CHAC
      Character*(*)     CHFH

* Output fields 
      
      Character*64      CHLO
      Character*64      CHMO
      Character*64      CHMI
      Character*512     RCEX

* prepare output fields 

* - profile name (motif)

      CHMO='motif= '
      
      L1=Lblnk(CHID)
      L2=Lblnk(CHAC)
      L3=6
      If(L2-1.GT.0) then
         CHMO=CHMO(1:L3) // CHAC(1:L2)
         L3=L3+L2
      End if
      If(L1.GT.0) then
         CHMO=CHMO(1:L3) // CHID(1:L1)
      End if
      If(L2-1.LE.0.AND.L1.EQ.0) then
         CHMO=CHMO(1:6) // 'unknown'
      End if

      ILMO=Lblnk(CHMO)

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
      End if

* - location
      Write(CHLO,*) JALB,'-',JALE
      J1=0
      Do I1=1,Lblnk(CHLO)
          If(CHLO(I1:I1).NE.' ') then 
             J1=J1+1
             CHLO(J1:J1)=CHLO(I1:I1)
             CHLO(I1:I1)=' '
          End if
      End do
      ILLO=Lblnk(CHLO)

* assemble ouput record

      LNEX=0 

      RCEX(LNEX+1:)='>' // CHMI 
      LNEX=LNEX+ILMI+2
         
      RCEX(LNEX:)='/' // CHLO
      LNEX=LNEX+ILLO+1

      RCEX(LNEX+1:)=CHMO 
      LNEX=LNEX+ILMO+1

* Write output record
      If(LNEX.GT.512) then
         LNEX=512
      End if
      Write(6,'(512A)')(RCEX(ii1:ii1),ii1=1,LNEX)

 100  Return 
      End 
