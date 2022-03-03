#!/bin/sh -ve

#----------------------------------------------------------------------#
# PATHS declaration (set by cmake/make: do not modify them here)
#----------------------------------------------------------------------#
PFSEARCHV3=$<TARGET_FILE:pfsearchV3>

TMPDIR=/tmp/test_pfsearchV3_pattern
mkdir -p $TMPDIR

#----------------------------------------------------------------------#
# The PFTOOLS is a powerful software to align biological sequences.
# Owing to the 'generalized profile syntax', it allows the fine-tuning
# of an alignent scoring system, beyond what is feasible by most other
# software. Despite the PFTOOLS are crippled by a lot of legacy code,
# they are still extremely useful for precision work.
#
# Nota bene to use this script as a testsuite:
# (1) The output order of pfsearch is reproducible, as well as the one of
#     pfsearchV3 with -t 1.
# (2) Refrain using any pipe.
#----------------------------------------------------------------------#

#----------------------------------------------------------------------#
# The following output formats are preserved between V2 and V3
#
# Default output mode:
#----------------------------------------------------------------------#
$PFSEARCHV3 -n -f ./sh3.prf ./VAV_HUMAN.seq | sort -nr > $TMPDIR/sh3.VAV_HUMAN.3.hit
$PFSEARCHV3 -n -fb ./ecp.prf $TMPDIR/CVPBR322.fa | sort -nr > $TMPDIR/ecp.CVPBR322.3.hit
#diff -b $TMPDIR/ecp.CVPBR322.2.hit $TMPDIR/ecp.CVPBR322.3.hit  # expecting no difference


