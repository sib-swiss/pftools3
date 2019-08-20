#!/bin/sh -e

#----------------------------------------------------------------------#
# PATHS declaration (set by cmake/make: do not modify them here)
#----------------------------------------------------------------------#
TMPDIR=/tmp/test_V3
mkdir -p $TMPDIR

#----------------------------------------------------------------------#
# This script verify the correctness of the output of test_V3.sh
#----------------------------------------------------------------------#

sh -ve ./test_V3.sh  > $TMPDIR/test_V3.full 2>&1
perl -ne '$ok=1 if/^### TESTS/;print $_ if $ok' < $TMPDIR/test_V3.full > $TMPDIR/test_V3.out
diff ./test_V3.out $TMPDIR/test_V3.out




