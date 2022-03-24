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
# TESTS FOR PATTERN AND REGEXP SEARCHES
# Run different patterns, and the corresponding regexps, and compare to the expected outputs
#----------------------------------------------------------------------#

# The PROSITE Pattern syntax is described in details in the User Manual (https://prosite.expasy.org/prosuser.html).
# Here is a summary:
# - The standard IUPAC one-letter codes for the amino acids are used.
# - The symbol 'x' is used for a position where any amino acid is accepted.
# - Ambiguities are indicated by listing the acceptable amino acids for a given position, between square parentheses '[ ]'. For example:[ALT] stands for Ala or Leu or Thr.
# - Ambiguities are also indicated by listing between a pair of curly brackets '{ }' the amino acids that are not accepted at a given position. For example: {AM} stands for any amino acid except Ala and Met.
# - Each element in a pattern is separated from its neighbor by a '-'.
# - Repetition of an element of the pattern can be indicated by following that element with a numerical value or a numerical range between parenthesis. Examples: x(3) corresponds to x-x-x, and x(2,4) corresponds to x-x or x-x-x or x-x-x-x.
# - When a pattern is restricted to either the N- or C-terminal of a sequence, that pattern either starts with a '<' symbol or respectively ends with a '>' symbol.
#
# e.g.:
# C-x-C-x(2)-[GP]-[FYW]-x(4,8)-C       (spaces, tabs, and line breaks are ignored)
#     can be translated as
# Cys-any-Cys-any-any-[Gly or Pro]-[Phe or Tyr or Trp]-'four to eight any'-Cys


#----------------------------------------------------------------------#
# Search with fully defined pattern/regexp
#----------------------------------------------------------------------#
$PFSEARCHV3 'pattern{L-N-N-L-L-P-H-A-I-N-L-R-E-V-N-L-R-P-Q-M-S-Q}'  ./VAV_HUMAN.seq >$TMPDIR/fully_defined-pattern.hit
$PFSEARCHV3 'regex{LNNLLPHAINLREVNLRPQMSQ}'                         ./VAV_HUMAN.seq >$TMPDIR/fully_defined-regexp.hit
diff -b $TMPDIR/fully_defined-pattern.hit patterns/fully_defined.hit  # expecting no difference
diff -b $TMPDIR/fully_defined-regexp.hit  patterns/fully_defined.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search with undefined symbol
#----------------------------------------------------------------------#
#NOTE X does not work, need to use x! X is used to search for X character directly used in sequences!
$PFSEARCHV3 'pattern{S-x-C-C-x-K-F}'  ./VAV_HUMAN.seq >$TMPDIR/undefined-pattern.hit
$PFSEARCHV3 'regex{S.CC.KF}'          ./VAV_HUMAN.seq >$TMPDIR/undefined-regexp.hit
diff -b $TMPDIR/undefined-pattern.hit patterns/undefined.hit  # expecting no difference
diff -b $TMPDIR/undefined-regexp.hit  patterns/undefined.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search with ambiguities
#----------------------------------------------------------------------#
$PFSEARCHV3 'pattern{Y-D-C-[VE]-[EV]-N-E-E}'  ./VAV_HUMAN.seq >$TMPDIR/ambiguities-pattern.hit
$PFSEARCHV3 'regex{YDC[VE][EV]NEE}'           ./VAV_HUMAN.seq >$TMPDIR/ambiguities-regexp.hit
diff -b $TMPDIR/ambiguities-pattern.hit patterns/ambiguities.hit  # expecting no difference
diff -b $TMPDIR/ambiguities-regexp.hit  patterns/ambiguities.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search with ambiguities excluded
#----------------------------------------------------------------------#
$PFSEARCHV3 'pattern{M-V-P-M-{PY}-R-V-L-K}'  ./VAV_HUMAN.seq >$TMPDIR/ambiguities_excluded-pattern.hit
$PFSEARCHV3 'regex{MVPM[^PY]RVLK}'           ./VAV_HUMAN.seq >$TMPDIR/ambiguities_excluded-regexp.hit
diff -b $TMPDIR/ambiguities_excluded-pattern.hit patterns/ambiguities_excluded.hit  # expecting no difference
diff -b $TMPDIR/ambiguities_excluded-regexp.hit  patterns/ambiguities_excluded.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search with repetition of element(s)
#----------------------------------------------------------------------#
$PFSEARCHV3 'pattern{R-D(2)-S(2)-G-D}'  ./VAV_HUMAN.seq >$TMPDIR/repetition-pattern.hit
$PFSEARCHV3 'regex{RD{2}S{2}GD}'        ./VAV_HUMAN.seq >$TMPDIR/repetition-regexp.hit
$PFSEARCHV3 'pattern{E-L-K(2,3)-W-M}'   ./VAV_HUMAN.seq >$TMPDIR/repetition_range-pattern.hit
$PFSEARCHV3 'regex{ELK{2,3}WM}'         ./VAV_HUMAN.seq >$TMPDIR/repetition_range-regexp.hit
diff -b $TMPDIR/repetition-pattern.hit       patterns/repetition.hit        # expecting no difference
diff -b $TMPDIR/repetition-regexp.hit        patterns/repetition.hit        # expecting no difference
diff -b $TMPDIR/repetition_range-pattern.hit patterns/repetition_range.hit  # expecting no difference
diff -b $TMPDIR/repetition_range-regexp.hit  patterns/repetition_range.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search with anchor to either the N- or C-terminal of a sequence
#----------------------------------------------------------------------#
#FIXME pattern anchors do not work! So don't know if it has to be < or <-
#FIXME regexp anchor does not look to accept $!
$PFSEARCHV3 'pattern{<M-E-L-x(2)-Q-C}'      ./VAV_HUMAN.seq >$TMPDIR/anchorN-pattern.hit
$PFSEARCHV3 'regex{^MEL.{2}QC}'             ./VAV_HUMAN.seq >$TMPDIR/anchorN-regexp.hit
$PFSEARCHV3 'regex{\AMEL.{2}QC}'            ./VAV_HUMAN.seq >$TMPDIR/anchorNalt-regexp.hit
$PFSEARCHV3 'pattern{V-E(2)-D-Y-S-E-Y-C>}'  ./VAV_HUMAN.seq >$TMPDIR/anchorC-pattern.hit
$PFSEARCHV3 'regex{VE{2}DYSEYC\E}'          ./VAV_HUMAN.seq >$TMPDIR/anchorC-regexp.hit
#diff -b $TMPDIR/anchorN-pattern.hit patterns/anchorN.hit  # expecting no difference
diff -b $TMPDIR/anchorN-regexp.hit  patterns/anchorN.hit  # expecting no difference
#diff -b $TMPDIR/anchorC-pattern.hit patterns/anchorC.hit  # expecting no difference
diff -b $TMPDIR/anchorC-regexp.hit  patterns/anchorC.hit  # expecting no difference

#----------------------------------------------------------------------#
# Search giving multiple matches per sequences
#----------------------------------------------------------------------#
$PFSEARCHV3 'pattern{x(6)-A-x-[YRFQ]-D}'  ./VAV_HUMAN.seq >$TMPDIR/multiple_matches-pattern.hit
$PFSEARCHV3 'regex{.{6}A.[YRFQ]D}'        ./VAV_HUMAN.seq >$TMPDIR/multiple_matches-regexp.hit
diff -b $TMPDIR/multiple_matches-pattern.hit patterns/multiple_matches.hit  # expecting no difference
diff -b $TMPDIR/multiple_matches-regexp.hit  patterns/multiple_matches.hit  # expecting no difference

# Cleaning
rm -Rf $TMPDIR

