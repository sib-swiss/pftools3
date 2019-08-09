*----------------------------------------------------------------------*     
* $Id: pfdat.f,v 2.5 2003/02/19 08:38:25 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
* PROFILE 

      Integer           B0 
      Integer           B1 
      Integer           E0 
      Integer           E1 
      Integer           BM 
      Integer           BI 
      Integer           BD 
      Integer           BE 
      Integer           MM 
      Integer           MI 
      Integer           MD 
      Integer           ME 
      Integer           DM 
      Integer           DI 
      Integer           DD 
      Integer           DE 
      Integer           IM 
      Integer           II 
      Integer           ID 
      Integer           IE 

      Integer           M0
      Integer           I0
      Integer           D 

      Parameter        (B0=  27) 
      Parameter        (B1=  28)
      Parameter        (E0=  29)
      Parameter        (E1=  30)
      Parameter        (BM=  31)
      Parameter        (BI=  32)
      Parameter        (BD=  33)
      Parameter        (BE=  34)
      Parameter        (MM=  35)
      Parameter        (MI=  36)
      Parameter        (MD=  37)
      Parameter        (ME=  38)
      Parameter        (IM=  39)
      Parameter        (II=  40)
      Parameter        (ID=  41)
      Parameter        (IE=  42)
      Parameter        (DM=  43)
      Parameter        (DI=  44)
      Parameter        (DD=  45)
      Parameter        (DE=  46)

      Parameter        (M0=   0)
      Parameter        (I0=   0)
      Parameter        (D =  27)
* 
      Integer           NLOW
      Parameter        (NLOW=-536870912)   
      
      Character*01      CHIP(0:IDMP)
      Integer           IIPP(0:46,0:IDMP)     
      
      Character*01      CHMP(IDMP)
      Integer           IMPP(0:27,0:IDMP)
      
* Storage of circular profile last insert position

      Character*01      CHIL
      Integer           IIPL(0:46)
      Integer           ILIP
      
