#!/bin/sh -e

#----------------------------------------------------------------------#
# PATHS declaration (set by cmake/make: do not modify them here)
#----------------------------------------------------------------------#
TMPDIR=/tmp/test_V3
mkdir -p $TMPDIR

#----------------------------------------------------------------------#
# This script verify the correctness of the output of test_V3.sh
#----------------------------------------------------------------------#

sh -ve ./test_V3.sh  > $TMPDIR/test_V3.out 2>/dev/null
diff ./test_V3.out $TMPDIR/test_V3.out

# Cleaning
rm -Rf $TMPDIR

