*----------------------------------------------------------------------*     
* $Id: Xblnk.f,v 1.3 2003/11/28 11:56:08 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Integer Function  Xblnk(STRING,LIMIT) 
      Character*(*)     STRING
      Integer           LIMIT

      L=Len(STRING)
      Xblnk=0
      J1=0

      Do I1=1,L
         If(STRING(I1:I1).EQ.' ') then
            J1=J1+1
         Else
            J1=0
         End if
         Xblnk=I1-J1
         If(J1.GE.LIMIT) return
      End do 
      Return
      End
