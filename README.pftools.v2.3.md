# Description of *pftools* release 2.3 package and files

Contents:

1. Package description
..* Profile format
..* Programs and manual pages
..* Compilation issues
2. Installation
..* Compiling
..* Installing
3. File list
4. Testing
5. Contact


## Package description
The *pftools* package is a collection of experimental programs supporting the
generalized profile format and search method of PROSITE.
A description of the generalized profile format is given in the file:

   [doc/profile.txt](https://raw.githubusercontent.com/sib-swiss/pftools3/master/doc/profile.txt)

The official web site for the *pftools* package can be found under:

   [pftools](https://web.expasy.org/pftools/)

Further information on PROSITE specific applications of the *pftools* programs
can be found at the following URL:

   [PROSITE](https://prosite.expasy.org/)

Release 2.3 contains FORTRAN 77 source code and manual pages for the following
programs:

```
   pfsearch
   pfscan
   psa2msa
   gtop
   htop
   ptoh
   pfmake
   pfw
   pfscale
   ptof
   6ft
   2ft
```

Some supplementary manual pages describing the *psa* and *xpsa* file formats
have also been included.

These programs have successfully been compiled and tested under Solaris,
IRIX64, HP-UX, Tru64 and Linux operating systems, using the GNU g77 or f77
compilers.

For compilation with native f77 compilers, the Makefile needs to be modified
according to instructions given therein.

For use under HP-UX, the file *sterr.f* should also be modified according to
instructions given therein.

The buffer sizes in release 2.3 of the *pftools* have been significantly
increased to handle large sequences and profiles. Therefore it might be
necessary to reduce the buffer sizes in order to run the programs on machines
with less RAM. This can be achieved by modifying the file *ardim.f* according
to instructions given therein.

The "Integer*2" declarations may cause problems with early versions
of the g77 compiler; solution: replace "Integer*2" by "Integer".


## Installation
The compressed *tar* archive can be extracted to the current directory using
the following command:

```
   tar -xzvf pft2.3.tar.gz
```

This will create a directory containing the source and package files.

In order to compile the programs on Unix like systems with a FORTRAN 77
compiler, change into the source directory and simply type:

```
   make all
```

This will create the 12 binaries listed above. To install the binaries and the
corresponding manual pages, an installation script has been provided. It can
be invoked using:

```
   make install
```

The installation script will ask the user the destination directory for the
*pftools* package. This directory will contain all the supplementary files,
as well as a *bin* directory for the executable programs. The installer will
then create symbolic links to the binary programs in a directory specified by
the user. This directory should be included in the users *PATH* variable for
him to be able to execute the *pftools* programs.
In a last step, the installer will link the man pages in the system wide man
page location, this location can also be specified by the user.
After successful installation, the build directory can safely be deleted.

Note that you should have write permissions on the destination directories
specified during installation. For system wide installation you may need root
privileges, if necessary contact your system administrator.


## File list
Further included in this release are the following demo data files:

```
   CVPBR322
   sh3.seq
   GTPA_HUMAN
   sh3.msf
   ecp.prf
   sh3.gpr
   sh3.prf
   pfam_sh3.hmm
   prosite13.prf
   standard.random
   blosum45.cmp
   score.lis
   coils.prf
   MYSA_HUMAN
   R76849.seq
```

plus the following additional substution matrices in old GCG format:

```
   blosum30.cmp
   blosum50.cmp
   blosum62.cmp
   blosum65.cmp
   blosum80.cmp
   blosum100.cmp
   gonnet.cmp
   pam30.cmp
   pam40.cmp
   pam80.cmp
   pam120.cmp
   pam160.cmp
   pam200.cmp
   pam220.cmp
   pam250.cmp
   pam400.cmp
```

## Testing
A test script is provided. It will simply execute all the *pftools* binaries
with some example files. The test script *test.sh* should be executed in
the same directory as the binaries and the demo files. An example test
output file *test.out* is provided as a reference.
To test the programs, type:

```
   ./test.sh > out
   diff test.out out
```

Some rounding or formatting variation may occur with real number
editing. With g77 there will be lots of small integer rounding
differences in the result of test 7.

Note that the file pfam_sh3.hmm contains a hidden Markov model from the
PFAM A collection release 4.0 (see [Pfam](http://pfam.xfam.org/)).


## Contact
Please send bug reports to:

   pftools@sib.swiss

The pftools package was originaly developed by _Philipp Bucher_ and was
maintained by _Thierry Schuepbach_ at the:

```
SIB Swiss Institute of Bioinformatics
Vital-IT Group
Quartier Sorge, Batiment Amphipole
CH-1015 Lausanne
Switzerland
```

[SIB](https://www.sib.swiss)

