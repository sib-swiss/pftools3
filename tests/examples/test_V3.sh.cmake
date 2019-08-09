#!/bin/sh -ve

#----------------------------------------------------------------------#
# PATHS declaration (set by cmake/make: do not modify them here)
#----------------------------------------------------------------------#
PFSEARCH=$<TARGET_FILE:pfsearch>
PFSCAN=$<TARGET_FILE:pfscan>
PFSEARCHV3=$<TARGET_FILE:pfsearchV3>
PFSCANV3=$<TARGET_FILE:pfscanV3>
PFCALIBRATEV3=$<TARGET_FILE:pfcalibrateV3>
PFINDEX=$<TARGET_FILE:pfindex>

GTOP=$<TARGET_FILE:gtop>
HTOP=$<TARGET_FILE:htop>
PFSCALE=$<TARGET_FILE:pfscale>
PFMAKE=$<TARGET_FILE:pfmake>
PSA2MSA=$<TARGET_FILE:psa2msa>
PFW=$<TARGET_FILE:pfw>
PTOH=$<TARGET_FILE:ptoh>
PTOF=$<TARGET_FILE:ptof>
P2FT=$<TARGET_FILE:2ft> # NB: sh does not allow variable name starting with a digit
 
SORT_PSA=@PERL_SCRIPT_DIR@/sort_fasta.pl # FIXME: use cmake syntax 
MAKE_IUPAC_CMP=@PERL_SCRIPT_DIR@/make_iupac_cmp.pl # FIXME: use cmake syntax 
SCRAMBLE=@PERL_SCRIPT_DIR@/scramble_fasta.pl # FIXME: use cmake syntax

CMPDIR=@DATA_DIR@/Matrices
TMPDIR=/tmp

#----------------------------------------------------------------------#
# The PFTOOLS is a powerful software to align biological sequences. 
# Owing to the 'generalized profile syntax', it allows the fine-tuning 
# of an alignent scoring system, beyond what is feasible by most other 
# software. Despite the PFTOOLS are crippled by a lot of legacy code, 
# they are still extremely useful for precision work . 
#
# Here below are given examples of commands for the latest V3 version. 
#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Searching for the occurence of the SH3 domain within the VAV oncogene,
# using pfsearch V2 ...
#----------------------------------------------------------------------#
$PFSEARCH -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr

#----------------------------------------------------------------------#
# ...and using pfsearch V3:
#----------------------------------------------------------------------#
$PFSEARCHV3 -n -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr

#----------------------------------------------------------------------#
# Create a database of sequences and a database of profiles, each one 
# with two entries. 
#----------------------------------------------------------------------#
cat ./VAV_HUMAN.seq ./VAV_RAT.seq > $TMPDIR/VAV.seq
cat ./sh2.prf       ./sh3.prf     > $TMPDIR/SHX.prf

#----------------------------------------------------------------------#
# The following commands must all produce the same final results
#
# (1) pfsearch (V2) expects ONE profile and MANY sequences:
# (2) pfscan (V2) expects MANY profile and ONE sequences:
# (3) pfsearch (V3) expects ONE profile and MANY sequences:
# (4) pfscan (V3) supports MANY profiles and MANY sequences:
#----------------------------------------------------------------------#
$PFSEARCH -fkxz ./sh2.prf $TMPDIR/VAV.seq > $TMPDIR/SHX.pfsearch2.hit
$PFSEARCH -fkxz ./sh3.prf $TMPDIR/VAV.seq >> $TMPDIR/SHX.pfsearch2.hit

$PFSCAN -fkxz ./VAV_HUMAN.seq $TMPDIR/SHX.prf > $TMPDIR/SHX.pfscan2.hit
$PFSCAN -fkxz ./VAV_RAT.seq $TMPDIR/SHX.prf >> $TMPDIR/SHX.pfscan2.hit

$PFSEARCHV3 -f -n -o 6 ./sh2.prf -f $TMPDIR/VAV.seq > $TMPDIR/SHX.pfsearch3.hit
$PFSEARCHV3 -f -n -o 6 ./sh3.prf -f $TMPDIR/VAV.seq >> $TMPDIR/SHX.pfsearch3.hit

$PFSCANV3 -f -n -o 6 $TMPDIR/SHX.prf $TMPDIR/VAV.seq > $TMPDIR/SHX.pfscan3.hit

#----------------------------------------------------------------------#
# All these commands produces exactly the same list of matched 
# sequences, with the same raw scores and coordinates. 
#
# However the output order is not preserved. 

# If properly calibrated, the V3 heuristic should produce the same 
# results.
#
# Let's verify that the output are comparable (FASTA headers need to be 
# fixed, however):
#----------------------------------------------------------------------#
$SORT_PSA -s $TMPDIR/SHX.pfscan2.hit   > $TMPDIR/SHX.pfscan2.out
$SORT_PSA -s $TMPDIR/SHX.pfscan3.hit   > $TMPDIR/SHX.pfscan3.out
$SORT_PSA -s $TMPDIR/SHX.pfsearch2.hit > $TMPDIR/SHX.pfsearch2.out
$SORT_PSA -s $TMPDIR/SHX.pfsearch3.hit > $TMPDIR/SHX.pfsearch3.out
diff $TMPDIR/SHX.pfscan2.out $TMPDIR/SHX.pfscan3.out    # expecting no difference
diff $TMPDIR/SHX.pfscan2.out $TMPDIR/SHX.pfsearch2.out  # expecting no difference
diff $TMPDIR/SHX.pfscan2.out $TMPDIR/SHX.pfsearch3.out  # expecting no difference

#----------------------------------------------------------------------#
# Pfsearch/pfscan V2 supports the following input formats for sequence: 
# FASTA, Swiss-Prot and EMBL.
#----------------------------------------------------------------------#
$PFSEARCH -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr    # FASTA
$PFSCAN   -f ./VAV_HUMAN.seq ./sh3.prf | sort -nr    # FASTA

$PFSEARCH ./sh3.prf ./GTPA_HUMAN.dat | sort -nr      # SwissProt
$PFSCAN   ./GTPA_HUMAN.dat ./sh3.prf | sort -nr      # SwissProt

$PFSEARCH  ./ecp.prf ./CVPBR322.embl | sort -nr      # EMBL
$PFSCAN    ./CVPBR322.embl ./ecp.prf | sort -nr      # EMBL 

#----------------------------------------------------------------------#
# Pfsearch/pfscan V3 supports the following input formats for sequence: 
# FASTA and FASTQ.
#----------------------------------------------------------------------#
$PFSEARCHV3 -n -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr  # FASTA
$PFSCANV3   -n -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr  # FASTA
cat ./CVPBR322.embl \
| perl -ne 'if(/^ID +(\w+)/){print ">$1\n"}elsif(/^ /){s/[\s\d]+//g;print "$_\n";}' \
> $TMPDIR/CVPBR322.fa
$PFSEARCHV3 -n -f ./ecp.prf $TMPDIR/CVPBR322.fa | sort -nr
$PFSCANV3  -n -f ./ecp.prf $TMPDIR/CVPBR322.fa  | sort -nr

$PFSEARCHV3 -n -q ./hiv.prf ./hiv.fastq | sort -nr    # FASTQ
$PFSCANV3   -n -q ./hiv.prf ./hiv.fastq | sort -nr    # FASTQ

#----------------------------------------------------------------------#
# The following output formats are preserved between V2 and V3
#
# Default output mode:
#----------------------------------------------------------------------#
$PFSEARCH      -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr > $TMPDIR/sh3.VAV_HUMAN.2.hit
$PFSEARCHV3 -n -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr > $TMPDIR/sh3.VAV_HUMAN.3.hit
diff -b  $TMPDIR/sh3.VAV_HUMAN.2.hit $TMPDIR/sh3.VAV_HUMAN.3.hit # expecting no difference # FIXME: could we remove -b ?

$PFSEARCH      -fb ./ecp.prf $TMPDIR/CVPBR322.fa | sort -nr > $TMPDIR/ecp.CVPBR322.2.hit
$PFSEARCHV3 -n -fb ./ecp.prf $TMPDIR/CVPBR322.fa | sort -nr > $TMPDIR/ecp.CVPBR322.3.hit
diff -b $TMPDIR/ecp.CVPBR322.2.hit $TMPDIR/ecp.CVPBR322.3.hit  # expecting no difference # FIXME: it does not work

#----------------------------------------------------------------------#
# The different PSA (i.e FASTA) output formats differ by the content 
# their headers. The usage of XPSA (-o 6) is strongly recommended, the 
# other formats being preserved for legacy.
#----------------------------------------------------------------------#
$PFSEARCHV3 -n -f -o 1 -N 15 ./sh3.prf ./VAV_HUMAN.seq | $SORT_PSA -
$PFSEARCHV3 -n -f -o 2 -N 15 ./sh3.prf ./VAV_HUMAN.seq | $SORT_PSA -
$PFSEARCHV3 -n -f -o 3 -N 15 ./sh3.prf ./VAV_HUMAN.seq | $SORT_PSA -
$PFSEARCHV3 -n -f -o 4 -N 15 ./sh3.prf ./VAV_HUMAN.seq | $SORT_PSA -
$PFSEARCHV3 -n -f -o 6 -N 15 ./sh3.prf ./VAV_HUMAN.seq | $SORT_PSA - # recommended

#----------------------------------------------------------------------#
# Multiple sequence alignment (MSA) can easily be produced from PSA 
# using psa2msa (V2)
#----------------------------------------------------------------------#
$PFSEARCHV3 -n -f -o 6 ./ecp.prf $TMPDIR/CVPBR322.fa | $SORT_PSA - | $PSA2MSA 

#----------------------------------------------------------------------#
# V3 has two new output formats: TSV and SAM
#----------------------------------------------------------------------#

$PFSEARCHV3 -n -f    -o 7 ./ecp.prf $TMPDIR/CVPBR322.fa | sort # FASTA to TSV
$PFSEARCHV3 -n -f -b -o 7 ./ecp.prf $TMPDIR/CVPBR322.fa | sort # FASTA to TSV

$PFSEARCHV3 -n -f    -o 8 ./ecp.prf $TMPDIR/CVPBR322.fa | sort # FASTA to SAM
$PFSEARCHV3 -n -f -b -o 8 ./ecp.prf $TMPDIR/CVPBR322.fa | sort # FASTA to SAM

$PFSEARCHV3 -n -q    -o 7 ./hiv.prf ./hiv.fastq  | sort # FASTQ to TSV
$PFSEARCHV3 -n -q -b -o 7 ./hiv.prf ./hiv.fastq  | sort # FASTQ to TSV

$PFSEARCHV3 -n -q    -o 8 ./hiv.prf ./hiv.fastq  | sort # FASTQ to SAM
$PFSEARCHV3 -n -q -b -o 8 ./hiv.prf ./hiv.fastq  | sort # FASTQ to SAM

#----------------------------------------------------------------------#
# Index can optionaly be used to speed-up database upload for pfsearch. 
# This should be especially usefull for repetitive database searches 
# using the heuristic. 
#----------------------------------------------------------------------#

$PFINDEX -f -o $TMPDIR/CVPBR322.fa.idx $TMPDIR/CVPBR322.fa
$PFSEARCHV3 -n -f                            -o 7 ./ecp.prf $TMPDIR/CVPBR322.fa  | sort > $TMPDIR/A.out    # FASTA to TSV
$PFSEARCHV3 -n -f -i $TMPDIR/CVPBR322.fa.idx -o 7 ./ecp.prf $TMPDIR/CVPBR322.fa  | sort > $TMPDIR/B.out    # FASTA to TSV
diff $TMPDIR/A.out  $TMPDIR/B.out 

$PFINDEX -q -o $TMPDIR/hiv.fastq.idx ./hiv.fastq
$PFSEARCHV3 -n -q -b -o 8 ./hiv.prf ./hiv.fastq | sort > $TMPDIR/A.out
$PFSEARCHV3 -n -q -b -o 8 -i $TMPDIR/hiv.fastq.idx ./hiv.prf ./hiv.fastq | sort > $TMPDIR/B.out
diff $TMPDIR/A.out  $TMPDIR/B.out

#----------------------------------------------------------------------#
# Profile and/or sequence can be reversed on the fly. If both are 
# reversed the match score is the same
#----------------------------------------------------------------------#
$PFSEARCHV3 -f -a       rand_dna.prf rand_dna.seq
$PFSEARCHV3 -f -a -r -R rand_dna.prf rand_dna.seq

#----------------------------------------------------------------------#
# There is a much better support of rev/comp search in DNA sequences, 
# featuring the full IUPAC code
#----------------------------------------------------------------------#
# $PFSEARCH -n ./V4V6.prf V4V6.fastq  > V4V6.hit
# $PFSEARCH -n ./V4V6.recomp.prf V4V6.fastq  >> V4V6.hit

#----------------------------------------------------------------------#
# Or simply for the same result
#----------------------------------------------------------------------#

# $PFSEARCH -n -b ./V4V6.prf V4V6.fastq  > V4V6.b.hit

#----------------------------------------------------------------------#
# To illustrate these capabilities, let build a profile for bacterial
# 16S sequence, starting from the Rfam seed RF00177.msa. First create 
# an adhoc substitution matrix
#----------------------------------------------------------------------#
$MAKE_IUPAC_CMP -M 20 -m -30 > $TMPDIR/iupac20-30.cmp

#----------------------------------------------------------------------#
# ... let build a profile for bacterial 16S sequence, starting from 
# the Rfam seed RF00177.msa.
#----------------------------------------------------------------------#
cat ./RF00177.msa \
| perl -pe 'tr/U/T/ if ! /^>/' \
| $PFW -m -N 25 - \
| $PFMAKE -m -3 -S 0.01 -F 100 - $TMPDIR/iupac20-30.cmp \
| perl -pe 's/^ID   SEQUENCE_PROFILE/ID   16S/' \
> 16S.prf.tmp

#----------------------------------------------------------------------#
# Prepare a randomized databse to calibrate the profiles
#----------------------------------------------------------------------#
cat ./SRR9619541.sample.fastq \
| perl -ne '$x++;if( $x % 4 == 2 ){chomp;s/(.{1,60})/$1\n/g;print ">seq$x\n$_"}' \
| $SCRAMBLE -m permutation -M 10 - \
> $TMPDIR/permut.fa


$PFCALIBRATEV3 -F  $TMPDIR/permut.fa 16S.prf.tmp > 16S.prf

$PFSEARCHV3 -qn 16S.prf SRR9619541.sample.fastq | sort -nr



