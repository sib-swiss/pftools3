*----------------------------------------------------------------------*     
* $Id: nodat.f,v 2.4 2003/02/19 08:38:25 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* NORMALIZATION 
      
      Parameter        (KNOR=3)
      Parameter        (KNPM=5)
      Parameter        (MAXN=8)
      
      Character*16      CNOR(KNOR)
      Integer           JNOP(KNOR)
      
      Integer           JNOR
      Integer           MNOR(MAXN)
      Integer           NNOR(MAXN)
      Integer           NNPR(MAXN)
      Character*32      CNTX(MAXN)
      Real              RNOP(KNPM,MAXN)
      
