*       Program gtop
*----------------------------------------------------------------------*
* $Id: gtop.f,v 2.9 2003/11/28 11:53:33 vflegel Exp $
*----------------------------------------------------------------------*
*       Function: Reformats profiles: in-fmt=GRIBSKOV / out-fmt=PROSITE
*       Author:   Philipp Bucher
*       Contact:  pftools@sib.swiss
*       Version:  File under developpment for release 2.3
*----------------------------------------------------------------------*
*
* DATA
*----------------------------------------------------------------------*

* array dimensions and I/O units

      Include          'ardim.f'

      Parameter        (NOUT=   6)
      Parameter        (NGPR=  11)

* profile and sequence fields

      Character*4096    FGPR

      Include          'psdat.f'
      Include          'gsdat.f'
      Include          'djdat.f'
      Include          'nodat.f'
      Include          'codat.f'
      Include          'pfdat.f'
      Include          'dfdat.f'
      Include          'sterr.f'

* command line options and parameters

      Logical           LSYM
      Logical           LLLT
      Logical           LRNM

* initialization of controlled vocabularies

      Include          'cvini.f'

*----------------------------------------------------------------------*
* INPUT SECTION
*----------------------------------------------------------------------*

      IRC=0
      FGPR='-'
      BLOG=0.0
      LRNM=.FALSE.

* read command line

      RL=NLOW
      Call Repar
     *   (FGPR,LSYM,LLLT,RG,RE,RF,RO,IRC)
      If(IRC.NE.0) then
         Write(NERR,'(/,
     *      ''gtop 2.3 revision 5.d'',//
     *      ''Usage: gtop [ -aslhGEFO ] gcg-profile | - ''
     *      ''[ parameters ]'',/
     *      )')
         Write(NERR,'(
     *      ''   options:'',/,
     *      ''    -a: apply asymmetric gap weighting mode.'',/
     *      ''    -s: apply symmetric gap weighting mode (default).'',/
     *      ''    -l: do not impose limit on line length.'',/
     *      ''    -h: print usage help text.'',/
     *      ''    -G<value>:'',/
     *      ''        gap opening penalty (default: 4.5).'',/
     *      ''    -E<value>:'',/
     *      ''        gap extension penalty (default: 0.05).'',/
     *      ''    -F<value>:'',/
     *      ''        rescaling factor (default: 100).'',/
     *      ''    -O<value>:'',/
     *      ''        output score offset (default: 0). '',/
     *      ''        Added to match scores after multiplication by '',
     *      ''the rescaling factor F.'',/
     *      )')
         Write(NERR,'(
     *      ''   valid (but deprecated) parameters are:'',//,
     *      ''    [G=gap-weigth]              use option -G instead'',/,
     *      ''    [E=gap-extension-weight]    use option -E instead'',/,
     *      ''    [F=output-score-multiplier] use option -F instead'',/
     *      ''    [O=output-score-offset]     use option -O instead'',/
     *      )')
         Call Exit(IRC)
      End if

* read profile

      Call REGPR
     *   (NGPR,FGPR,
     *   RG,RE,RF,RO,LSYM,
     *   CPID,CPAC,CPDE,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   CHID,IIPD,CHMD,IMPD,
     *   IRC)

      If(IRC.NE.0) go to 100

      CPDT=' '
      LHDR=0
      LFTR=0

      Call WRPRF
     *   (NOUT,LLLT,LRNM,
     *   CPID,CPAC,CPDT,CPDE,LHDR,CHDR,LFTR,CFTR,NABC,CABC,LPRF,LPCI,
     *   CDIS,JDIP,MDIS,NDIP,
     *   CNOR,JNOP,JNOR,MNOR,NNOR,NNPR,CNTX,RNOP,
     *   JCUT,MCLE,CCUT,ICUT,JCNM,RCUT,MCUT,
     *   IDMP,CHIP,IIPP,CHMP,IMPP,
     *   BLOG,FABC,P0,
     *   CHID,IIPD,CHMD,IMPD,
     *   IRC)

 100  Call Exit(IRC)
      End
*----------------------------------------------------------------------*
      Subroutine Repar
     *   (FGPR,LSYM,LLLT,RG,RE,RF,RO,IRC)

      Character*4096    CPAR
      Character*(*)     FGPR

      Logical           LSYM
      Logical           LLLT

      LSYM=.TRUE.
      LLLT=.TRUE.
      RG=4.5
      RE=0.05
      RF=100
      RO=0

      FGPR=' '

      N1=Iargc()

      K1=0
      I2=1
      Do  50 I1=1,N1
         Call GetArg(I2,CPAR)
         If     (CPAR(1:1).EQ.'-'
     *      .AND.CPAR(2:2).NE.' '.AND.K1.EQ.0) then
            If(Index(CPAR,'h').NE.0)go to 900
            If(Index(CPAR,'a').NE.0) LSYM=.FALSE.
            If(Index(CPAR,'l').NE.0) LLLT=.FALSE.
            If(Index(CPAR,'G').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RG
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RG
               End if
            End if
            If(Index(CPAR,'E').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RE
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RE
               End if
            End if
            If(Index(CPAR,'F').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RF
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RF
               End if
            End if
            If(Index(CPAR,'O').NE.0) then
               If(CPAR(3:3).NE.' ') then
                  Read(CPAR(3:),*,Err=900) RO
               Else
                  I2=I2+1
                  Call GetArg(I2,CPAR)
                  Read(CPAR,*,Err=900) RO
               End if
            End if
         Else if(K1.EQ.0) then
            K1=K1+1
            FGPR=CPAR
         Else
            If     (CPAR(1:2).EQ.'G=') then
               Read(CPAR(3:),*,Err=900) RG
            Else if(CPAR(1:2).EQ.'E=') then
               Read(CPAR(3:),*,Err=900) RE
            Else if(CPAR(1:2).EQ.'F=') then
               Read(CPAR(3:),*,Err=900) RF
            Else if(CPAR(1:2).EQ.'O=') then
               Read(CPAR(3:),*,Err=900) RO
            End if
         End if
         I2=I2+1
         If(I2.GT.N1) Go to 60
 50   Continue

 60   If(K1.NE.1) Go to 900
 100  Return
 900  IRC=1
      Go to 100
      End
*----------------------------------------------------------------------*
      Include          'lblnk.f'
      Include          'regpr.f'
      Include          'wrprf.f'
