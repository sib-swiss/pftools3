#!/usr/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Std;
use Cwd qw( abs_path );
use FindBin qw( $RealBin );

use lib $RealBin;

use Toolbox;
my $tb = Toolbox->new();

my $PFSEARCH       = '$<TARGET_FILE:pfsearch>';
my $PFSEARCH_V3    = '$<TARGET_FILE:pfsearchV3>';
my $PFSCAN_V3      = '$<TARGET_FILE:pfscanV3>';
my $PFINDEX        = '$<TARGET_FILE:pfindex>';

# ------------------------------------------------------------
# Read command line and check arguments
# ------------------------------------------------------------

my $usage = "$0 [options] <profile-file> <sequence-file> <test-dir>

This script test if pfscanV3 and pfsearch are producing the same ouput

options: -h          this help
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'h', \%opt );
if( $opt{h} or @ARGV != 3 ){
    print "$usage\n";
    exit 0;
}
my( $profile_file, $fasta_file, $test_dir ) = @ARGV;
$tb->system( "mkdir -p $test_dir" ) unless -d $test_dir;
$tb->die( "Dir does not exist: $test_dir" ) unless -d $test_dir;

# ------------------------------------------------------------
# Check and prepare datasets
# ------------------------------------------------------------

my @prf_ID = ();

$tb->system( "mkdir -p $test_dir/prf" ) unless -d "$test_dir/prf";
$tb->system( "mkdir -p $test_dir/seq" ) unless -d "$test_dir/seq";
$tb->system( "mkdir -p $test_dir/hit" ) unless -d "$test_dir/hit";
my $all_fh = $tb->open( $profile_file );
my $prf_fh = undef;
my( $id, $type ) = ( '', '' );
while( my $line = <$all_fh> ){
    if( $line =~ /^ID/ ){
        ( $id, $type ) = $line =~ /^ID +(\S+); (\S+)\./;
        $tb->die( "Invalid first line: $line" ) unless $id and $type;
        $tb->die( "Unsupported profile type: $type" ) unless $type eq 'MATRIX';
        push @prf_ID, $id;
        $prf_fh = $tb->open( "> $test_dir/prf/$id.prf" );
    }
    print $prf_fh $line;
    if( $line =~ /^\/\// ){
        close $prf_fh;
    }
}
close $all_fh;

# foreach my $id ( @prf_ID ){
#    my $cmd = join ' ',
#        $PFCALIBRATE_V3,
#        '-F ~/gitlab/PfTools/data/Calibration/window20.seq',
#        '-H profile',
#        "$test_dir/prf/$id.tmp",
#        "> $test_dir/prf/$id.prf";
#    $tb->system( $cmd );
#}
#$tb->system( "rm -f $test_dir/prf/*.tmp $test_dir/prf/all.prf" );

# my $fasta_in  = $tb->open( $fasta_file );
# my $fasta_out = $tb->open( "> $test_dir/seq/all.seq" );
# print $fasta_out $_ while <$fasta_in>;
# close $fasta_in;
# close $fasta_out;
my $cmd = "$PFINDEX -f -o $test_dir/seq/all.idx $fasta_file";
$tb->system( $cmd );

# ------------------------------------------------------------
#
# ------------------------------------------------------------

my @thresh_option    = ( '', '-L -1', '-C 6.0', '-C 8.0', '-C 10.0' ); # Normalized score with pfsearch V2 syntax, i.e. -C <float>

my @sse2_option      = ( '', '--sse2' );
my @nthreads_option  = ( '', '-t 1', '-t 2' ); # none means all !
my @index_option     = ( '', "-i $test_dir/seq/all.idx" ); # order matter
my @heuristic_option = ( '-n', '' );

my @option_V3 = ();
foreach my $sse2 (@sse2_option ){
    foreach my $nthreads ( @nthreads_option ){
        foreach my $index ( @index_option ){
            foreach my $heuristic ( @heuristic_option ){
                push @option_V3, join ' ', $sse2, $nthreads, $index, $heuristic;
            }
        }
    }
}
# ------------------------------------------------------------
# Check that pfsearch (V2) pfscanV3 pfsearchV3 are producing
# the same aligmnet scores
# ------------------------------------------------------------
foreach my $thresh_option ( @thresh_option ){
    $tb->system( "rm -f $test_dir/hit/pfsearch.out" );
    foreach my $id ( @prf_ID ){
        my $cmd = join ' ',
                $PFSEARCH,
                '-f',
                $thresh_option,
                "$test_dir/prf/$id.prf",
                $fasta_file,
            ">> $test_dir/hit/pfsearch.out";
        $tb->system( $cmd );
    }
    my $cmd = q{perl -pe '/ (\d+) /;$_="$1\n"'} . " < $test_dir/hit/pfsearch.out | sort -nr > $test_dir/hit/score.ref.out";
    $tb->system( $cmd );

    my $thresh_option_V3 = $thresh_option;
    $thresh_option_V3 =~ s/\-C /\-N /;
    foreach my $option_V3 ( @option_V3 ){
        my $cmd = join ' ',
            $PFSCAN_V3,
                '-f',
                $thresh_option_V3,
                $option_V3,
                $profile_file,
                $fasta_file,
            "|", q{perl -pe '/ (\d+) /;$_="$1\n"'},
            "| sort -nr > $test_dir/hit/score.scan.out";
        $tb->system( $cmd );
        $tb->system( "diff $test_dir/hit/score.ref.out $test_dir/hit/score.scan.out > $test_dir/hit/diff.out" );
        $tb->die( `cat $test_dir/hit/diff.out` ) if -s "$test_dir/hit/diff.out";
        $tb->system( "rm -f $test_dir/hit/pfsearchV3.out" );
        foreach my $id ( @prf_ID ){
            my $cmd = join ' ',
                $PFSEARCH_V3,
                    '-f',
                    $thresh_option_V3,
                    $option_V3,
                    "$test_dir/prf/$id.prf",
                    $fasta_file,
                ">> $test_dir/hit/pfsearchV3.out";
            $tb->system( $cmd );
        }
        $cmd = q{perl -pe '/ (\d+) /;$_="$1\n"'} . " < $test_dir/hit/pfsearchV3.out | sort -nr > $test_dir/hit/score.search.out";
        $tb->system( $cmd );
        $tb->system( "diff $test_dir/hit/score.ref.out $test_dir/hit/score.scan.out > $test_dir/hit/diff.out" );
        $tb->die( `cat $test_dir/hit/diff.out` ) if -s "$test_dir/hit/diff.out";
        $tb->report( 'test Ok', $option_V3 );
        warn "\n";
    }
}
