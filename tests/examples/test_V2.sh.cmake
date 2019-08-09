#!/bin/sh -ve

PFSEARCH=$<TARGET_FILE:pfsearch>
PFSCAN=$<TARGET_FILE:pfscan>
GTOP=$<TARGET_FILE:gtop>
HTOP=$<TARGET_FILE:htop>
PFSCALE=$<TARGET_FILE:pfscale>
PFMAKE=$<TARGET_FILE:pfmake>
PFW=$<TARGET_FILE:pfw>
PTOH=$<TARGET_FILE:ptoh>
PTOF=$<TARGET_FILE:ptof>
X2FT=$<TARGET_FILE:2ft>
PSA2MSA=$<TARGET_FILE:psa2msa>

CMPDIR=@DATA_DIR@/Matrices

echo "#----------------------------------------------------------------------#"
echo "# PATHS"
echo "#----------------------------------------------------------------------#"
echo "# PFSEARCH=$<TARGET_FILE:pfsearch>"
echo "# PFSCAN=$<TARGET_FILE:pfscan>"
echo "# GTOP=$<TARGET_FILE:gtop>"
echo "# HTOP=$<TARGET_FILE:htop>"
echo "# PFSCALE=$<TARGET_FILE:pfscale>"
echo "# PFMAKE=$<TARGET_FILE:pfmake>"
echo "# PFW=$<TARGET_FILE:pfw>"
echo "# PTOH=$<TARGET_FILE:ptoh>"
echo "# PTOF=$<TARGET_FILE:ptof>"
echo "# 2FT=$<TARGET_FILE:2ft>"
echo "# PSA2MSA=$<TARGET_FILE:psa2msa>"
echo "#----------------------------------------------------------------------#"
echo "# pftools test 1: pfsearch -f -C 6.0 sh3.prf sh3.seq"
echo "#----------------------------------------------------------------------#"
$PFSEARCH -f sh3.prf VAV_HUMAN.seq
echo "#----------------------------------------------------------------------#"
echo "# pftools test 2: pfsearch -bx ecp.prf CVPBR322.embl | psa2msa -du"
echo "#----------------------------------------------------------------------#"
$PFSEARCH -bx ecp.prf CVPBR322.embl | $PSA2MSA -du
echo "#----------------------------------------------------------------------#"
echo "# pftools test 3: pfscan -s GTPA_HUMAN.dat prosite13.prf"
$PFSCAN -s GTPA_HUMAN.dat prosite13.prf
echo "#----------------------------------------------------------------------#"
echo "# pftools test 4: pfscan -by -C 2 CVPBR322.embl ecp.prf"
echo "#----------------------------------------------------------------------#"
$PFSCAN -by -C 2 CVPBR322.embl ecp.prf
echo "#----------------------------------------------------------------------#"
echo "# pftools test 5: gtop -F 50 sh3.gpr | pfsearch -far - sh3.seq"
echo "#                    | sort -nr"
echo "#----------------------------------------------------------------------#"
$GTOP -F 50 sh3.gpr | ${PFSEARCH} -far - sh3.seq | sort -nr
echo "#----------------------------------------------------------------------#"
echo "# pftools test 6: htop pfam_sh3.hmm  | pfsearch -f - sh3.seq"
echo "                     | sort -nr"
echo "#----------------------------------------------------------------------#"
$HTOP pfam_sh3.hmm | ${PFSEARCH} -f - sh3.seq | sort -nr
#echo "#----------------------------------------------------------------------#"
#echo "# pftools test 6b: htop htop_issue.hmm"
#echo "#----------------------------------------------------------------------#"
#FIXME Fortran runtime error: $HTOP htop_issue.hmm
echo "#----------------------------------------------------------------------#"
echo "# pftools test 7: pfw -N 1000 sh3.msf |"
echo "#                  pfmake -b -H 0.6 - blosum45.cmp"
echo "#----------------------------------------------------------------------#"
$PFW -N 1000 sh3.msf | ${PFMAKE} -b -H 0.6 - $CMPDIR/blosum45.cmp
echo "#----------------------------------------------------------------------#"
echo "# pftools test 8: ptoh -L 1.15 ecp.prf"
echo "#----------------------------------------------------------------------#"
$PTOH -L 1.15 ecp.prf
echo "#----------------------------------------------------------------------#"
echo "# pftools test 9: pfscale -N 14147368 -P 0.0001 -Q 0.000001 score.lis"
echo "#                  | sed -n 1,25p"
echo "#----------------------------------------------------------------------#"
$PFSCALE -N 14147368 -P 0.0001 -Q 0.000001 score.lis | sed -n 1,25p
echo "#----------------------------------------------------------------------#"
echo "# pftools test 10: pfsearch -y coils.prf MYSA_HUMAN"
echo "#----------------------------------------------------------------------#"
$PFSEARCH -y coils.prf MYSA_HUMAN

