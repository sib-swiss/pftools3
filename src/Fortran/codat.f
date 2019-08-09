*----------------------------------------------------------------------*     
* $Id: codat.f,v 2.4 2003/02/19 08:38:24 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* CUT_OFF 
      
      Parameter        (MAXC=8)
      
      Integer           JCUT 
      Integer           MCLE(MAXC)
      Character*32      CCUT(MAXC)
      Integer           ICUT(MAXC)
      
      Integer           JCNM(MAXC)
      Real              RCUT(MAXN,MAXC)
      Integer           MCUT(MAXN,MAXC)
      
