*----------------------------------------------------------------------*     
* $Id: ardim.f,v 2.7 2003/03/21 15:55:44 vflegel Exp $
*----------------------------------------------------------------------*     
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*     
      
* Note: for use on machines with less than 32 MB RAM, replace the actual 
*       array Parameter instructions by commented Parameter instructions 
*       or modify the array size definitions manually. 
      
* max. profile length  
      
      Parameter        (IDMP=19999)
C     Parameter        (IDMP= 9999)
C     Parameter        (IDMP= 2499)
      
* max. sequence length 
      
      Parameter        (IDMS=10000000)
C     Parameter        (IDMS= 1000000)
C     Parameter        (IDMS=  250000)
      
* max. number of matches per sequence
      
      Parameter        (IDMN=4000) 
C     Parameter        (IDMN=2000) 
C     Parameter        (IDMN= 500) 

* max. alignment length 

      Parameter        (IDMA=40000)
C     Parameter        (IDMA=20000)
C     Parameter        (IDMA= 5000)

* max. path matrix surface

      Parameter        (IDMM=40000000)
C     Parameter        (IDMM= 4000000)
C     Parameter        (IDMM= 1000000)

* max. repeats per circular profile match

      Parameter        (IDML=512)

* max. nb of sequences in MSF files

      Parameter        (IDMF=2048)

* max. nb of scores in score list

      Parameter        (IDMC=262144)
