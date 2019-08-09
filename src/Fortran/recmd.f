*----------------------------------------------------------------------*     
* $Id: recmd.f,v 2.4 2003/07/03 13:08:58 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      Subroutine RECMD(RCEX)

      Character*(*)     RCEX

      Character*512     CARG
      Character*512     CMDL

      LCMD=LEN(RCEX)

* Concatenate arguments of command-line

      N=Iargc()

      IC=1
      Do I1=0,N
         Call GetArg(I1,CARG)
         CMDL(IC:)=CARG
         IC=Lblnk(CMDL)+2
      End do

      RCEX=CMDL
      If(IC.GT.LCMD) RCEX(LCMD-1:LCMD)='..'

      Return
      End 
