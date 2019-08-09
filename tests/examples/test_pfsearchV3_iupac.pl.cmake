#!/bin/env perl

use strict;

use Data::Dumper;
use Getopt::Std;
use Cwd qw( abs_path );
use FindBin qw( $RealBin );

use lib $RealBin;

use Toolbox;
my $tb = Toolbox->new();

my $PFSEARCH_V3    = '$<TARGET_FILE:pfsearchV3>';
my $PFCALIBRATE_V3 = '$<TARGET_FILE:pfcalibrateV3>';
my $PFMAKE         = '$<TARGET_FILE:pfmake>';
my $PFSEARCH       = '$<TARGET_FILE:pfsearch>';
my $SORT_FASTA     = '@PERL_SCRIPT_DIR@/sort_fasta.pl';
my $CMP_DIR        = '@DATA_DIR@/Matrices';

my $TEST_FAILED    = 0; # global

# ------------------------------------------------------------
# Read command line and check arguments
# ------------------------------------------------------------

my $usage = "$0 [options] <test-dir>

This script tests if pfsearchV3 provides full support of DNA
including the full IUPAC codes. A random database is created.

options: -h         this help
         -L <int>   length of a seq in db ( 500 by default)
         -N <int>   db size               (1000 by default)
         -l <int>   profile length        ( 250 by defaul)
         -s <int>   random seed           (1234 by default)
";

my %opt;  # GLOBAL: to store the options
Getopt::Std::getopts( 'hL:N:l:s:', \%opt );
if( $opt{h} or @ARGV != 1 ){
    print "$usage\n";
    exit 0;
}
my( $test_dir ) = @ARGV;
$tb->die( "Dir does not exist: $test_dir" ) unless -d $test_dir;
$test_dir .= '/test.iupac';
$opt{L} = 500  unless exists $opt{L};
$opt{N} = 2000 unless exists $opt{N};
$opt{l} = 250  unless exists $opt{l};
$opt{s} = 1234 unless exists $opt{s};

# ------------------------------------------------------------
# Check and prepare datasets
# ------------------------------------------------------------
$tb->report( '###', 'Prepare a dataset of random sequences (IUPAC)');
$tb->report( 'db size', $opt{N} );
$tb->report( 'seq length', $opt{L} );
$tb->report( 'prf length', $opt{l} );
$tb->report( 'seed', $opt{s} );
srand( $opt{s} );

my @prf_ID = ();
$tb->system( "mkdir -p $test_dir/prf" ) unless -d "$test_dir/prf";
$tb->system( "mkdir -p $test_dir/seq" ) unless -d "$test_dir/seq";
$tb->system( "mkdir -p $test_dir/hit" ) unless -d "$test_dir/hit";

my %nuc_freq = ( # must be integer
    A => 400,
    C => 400,
    G => 300,
    T => 300, # need to consider U ?
    N =>  10, # any base
    R =>   1, # A or G
    Y =>   1, # C or T
    S =>   1, # G or C
    W =>   1, # A or T
    K =>   1, # G or T
    M =>   1, # A or C
    B =>   1, # C or G or T
    D =>   1, # A or G or T
    H =>   1, # A or C or T
    V =>   1, # A or C or G
);

my @master_dna = '';
foreach my $nuc ( sort keys %nuc_freq ){
    foreach( 1 .. $nuc_freq{$nuc} ){
        push @master_dna, $nuc;
    }
}
my $master_len = @master_dna;

sub get_rand_DNA{
    my $L = shift;
    my @dna = ();
    foreach( 1 .. $L ){
        push @dna, $master_dna[ int rand $master_len ] || 'N'; # Mmmh, bad behaviour of rand !!!!
    }
    return join '', @dna;
}

# my $fh = $tb->open( "> $test_dir/seq/calib.fa" );
# foreach my $n ( 1 .. $opt{N} ){
#    my $dna = get_rand_DNA( $opt{L} );
#    my $len = length( $dna );
#    $dna =~ s/(.{1,60})/$1\n/g;
#    print $fh ">dna_$len\_$n\n$dna";
# }
# close $fh;

my $fh = $tb->open( "> $test_dir/seq/db.fa" );
foreach my $n ( 1 .. $opt{N} ){
    my $dna = get_rand_DNA( $opt{L} );
    my $len = length( $dna );
    $dna =~ s/(.{1,60})/$1\n/g;
    print $fh ">dna_$len\_$n\n$dna";
}
close $fh;

$fh = $tb->open( "> $test_dir/prf/sample.fa" );
my $dna = get_rand_DNA( $opt{l} );
$dna =~ s/(.{1,60})/$1\n/g;
print $fh ">sample weight=1.0\n$dna";
close $fh;

warn "\n";

$tb->report( '###', 'Compare default versus reverse/reverse: scores must be always identical, coordinates often identical' );
my $cmd = join ' ',
    $PFMAKE,
    '-m',
    '-3',
    '-S 0.01 -F 100', # -F 1000 causes an oveflow!
    "$test_dir/prf/sample.fa",
    "$CMP_DIR/iupac_50_40.cmp",
    "> $test_dir/prf/sample.prf";

$tb->system( $cmd );
# $cmd = join ' ',
#    $PFCALIBRATE_V3,
#    '--seed 123',
#    "-F $test_dir/seq/calib.fa",
#    "-H $test_dir/prf/sample.fa",
#    '--pam-distance [0,200,10,10]',
##     "--report $test_dir/prf/report.pdf",
#    "$test_dir/prf/sample.tmp",
#    "> $test_dir/prf/sample2.prf";
# $tb->system( $cmd );

# ------------------------------------------------------------
# Check that reversing both the profile and the sequence has
# no effect on the best match scores. And compute a few stats
# about the matched sequences, because of alternative,
# possibly disjoint alignements with  the same score
# ------------------------------------------------------------
$cmd = join ' ',
    $PFSEARCH_V3,
        '-f -a -o 6',
        "$test_dir/prf/sample.prf",
        "$test_dir/seq/db.fa",
    "| $SORT_FASTA -", # just to ease debugging
    "> $test_dir/hit/out_1.psa";
$tb->system( $cmd );
$cmd = join ' ',
    $PFSEARCH_V3,
        '-f -a -o 6',
        '-R -r ',
        "$test_dir/prf/sample.prf",
        "$test_dir/seq/db.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/out_2.psa";
$tb->system( $cmd );
compare_psa_output( "$test_dir/hit/out_1.psa", "$test_dir/hit/out_2.psa" , 1, 0 );

# ------------------------------------------------------------
# Create a profile + its rev/comp version
# Nota Bene: the iupac.cmp matrix is symetrical w.r.t comp
#            local/local alignment mode is symmetrical
# ------------------------------------------------------------

$tb->report( '###', 'Compare scores and match coordinates V3 versus V2' );
$fh = $tb->open( "> $test_dir/prf/forward.fa" );
$dna = get_rand_DNA( $opt{l} );
$dna =~ s/(.{1,60})/$1\n/g;
print $fh ">sample weight=1.0\n$dna";
close $fh;
$cmd = join ' ',
     $SORT_FASTA,
    '-R -C',
    "$test_dir/prf/forward.fa",
    "> $test_dir/prf/revcomp.fa";
$tb->system( $cmd );
$cmd = join ' ',
    $SORT_FASTA,
    '-R -C',
    "$test_dir/seq/db.fa",
    "> $test_dir/seq/db.revcomp.fa";
$tb->system( $cmd );
$cmd = join ' ',
    $PFMAKE,
    '-m',
    '-3',
    '-S 0.01 -F 100',
    "$test_dir/prf/forward.fa",
    "$CMP_DIR/iupac_50_40.cmp",
    "> $test_dir/prf/forward.prf";
$tb->system( $cmd );
$cmd = join ' ',
    $PFMAKE,
    '-m',
    '-3',
    '-S 0.01 -F 100',
    "$test_dir/prf/revcomp.fa",
    "$CMP_DIR/iupac_50_40.cmp",
    "> $test_dir/prf/revcomp.prf";
$tb->system( $cmd );

# $cmd = join ' ',
#    $PFCALIBRATE_V3,
#    "-F $test_dir/seq/calib.fa",
#    "-H $test_dir/prf/forward.fa",
#    '--pam-distance [0,200,10,10]',
##     "--report $test_dir/prf/report.forward.pdf",
#    "$test_dir/prf/forward.tmp",
#    "> $test_dir/prf/forward.prf";
# $tb->system( $cmd );

#$cmd = join ' ',
#    $PFCALIBRATE_V3,
#    "-F $test_dir/seq/calib.fa",
#    "-H $test_dir/prf/revcomp.fa",
#    '--pam-distance [0,200,10,10]',
##     "--report $test_dir/prf/report.revcomp.pdf",
#    "$test_dir/prf/revcomp.tmp",
#    "> $test_dir/prf/revcomp.prf";
# $tb->system( $cmd );

# ------------------------------------------------------------
#
# ------------------------------------------------------------

$cmd = join ' ',
    $PFSEARCH_V3,
        '-f -a -o 6',
        "$test_dir/prf/forward.prf",
        "$test_dir/seq/db.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/forward.psa";
$tb->system( $cmd );
$cmd = join ' ',
    $PFSEARCH,
    '-faxzk',
    "$test_dir/prf/forward.prf",
    "$test_dir/seq/db.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/forward.2.psa";
$tb->system( $cmd );
$tb->report( '###', 'Compare V3 versus V2' );
compare_psa_output( "$test_dir/hit/forward.psa", "$test_dir/hit/forward.2.psa", 0, 0 );

$cmd = join ' ',
    $PFSEARCH_V3,
        '-f -a -o 6',
        "$test_dir/prf/revcomp.prf",
        "$test_dir/seq/db.revcomp.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/revcomp.psa";
$tb->system( $cmd );

$cmd = join ' ',
    $PFSEARCH,
    '-faxzk',
    "$test_dir/prf/revcomp.prf",
    "$test_dir/seq/db.revcomp.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/revcomp.2.psa";
$tb->system( $cmd );
$tb->report( '###', 'Compare V3 versus V2' );
compare_psa_output( "$test_dir/hit/revcomp.psa", "$test_dir/hit/revcomp.2.psa", 0, 0 );

$tb->report( '###', 'Compare default versus revcomp/revcomp (V2): scores must be identical' );
compare_psa_output( "$test_dir/hit/forward.2.psa", "$test_dir/hit/revcomp.2.psa", 1, 1 );
$tb->report( '###', 'Compare default versus revcomp/revcomp (V3): scores must be identical' );
compare_psa_output( "$test_dir/hit/forward.psa", "$test_dir/hit/revcomp.psa", 1, 1 );

exit 0;


$tb->system( "cat $test_dir/hit/forward.psa $test_dir/hit/revcomp.psa  | $SORT_FASTA - > $test_dir/hit/one_one.psa" );
$cmd = join ' ',
    $PFSEARCH_V3,
        '-f -a -o 6',
        '-b',
        "$test_dir/prf/forward.prf",
        "$test_dir/seq/db.fa",
    "| $SORT_FASTA -",
    "> $test_dir/hit/both.psa";
$tb->system( $cmd );
$tb->report( '###', 'Compare V3 versus V2' );
compare_psa_output( "$test_dir/hit/one_one.psa", "$test_dir/hit/both.psa", 1 );


# ------------------------------------------------------------
# Process the results
# ------------------------------------------------------------

sub compare_psa_output{
    my( $psa_file_1, $psa_file_2, $rev, $comp ) = @_; # rev/comp applies only to psa_file_2
    $tb->report( 'compare', "$psa_file_1 versus $psa_file_2" );
    my %match = ();
    my $id    = '';
    my $psa_1 = $tb->open( $psa_file_1 );
    while( <$psa_1> ){
        if( />(\w+)\/(\d+)\-(\d+)/ ){
            $id                    = $1;
            $match{$id}{for}{from} = $2 < $3 ? $2 : $3;
            $match{$id}{for}{to}   = $2 < $3 ? $3 : $2;
            ( $match{$id}{for}{score} ) = /raw_score=(\d+)/;
        }
        else{
            s/\s//g;
            $match{$id}{for}{seq} .= uc $_;
        }
    }
    my $psa_2 = $tb->open( $psa_file_2 );
    while( <$psa_2> ){
        if( />(\w+)\/(\d+)\-(\d+)/ ){
            $id                    = $1;
            $match{$id}{rev}{from} = $2 < $3 ? $2 : $3;
            $match{$id}{rev}{to}   = $2 < $3 ? $3 : $2;
            ( $match{$id}{rev}{score} ) = /raw_score=(\d+)/;
        }
        else{
            s/\s//g;
            $match{$id}{rev}{seq} .= uc $_;
        }
    }
    if( $rev ){
        foreach my $id ( sort keys %match ){
            next unless exists $match{$id}{rev};
#             $match{$id}{rev}{from} = $opt{L} + 1 - $match{$id}{rev}{to};
#             $match{$id}{rev}{to}   = $opt{L} + 1 - $match{$id}{rev}{from};
            $match{$id}{for}{seq}  = join '', reverse split '', $match{$id}{for}{seq};
        }
    }
    if( $comp ){
        $match{$id}{for}{seq} =~ tr/ACGTRYSWKMBDHVN/TGCAYRSWMKVHDBN/;
    }
    my $perfect_count  = 0;
    my $good_count     = 0;
    my $disjoint_count = 0;
    my $other_1_count  = 0;
    my $other_2_count  = 0;
    my $error_count    = 0;
    foreach my $id ( sort keys %match ){
        my $ok = 1;
        $ok = 0 if ! exists $match{$id}{for};
        $ok = 0 if ! exists $match{$id}{rev};
        if( $ok ){
            $ok = 0 unless $match{$id}{for}{score} == $match{$id}{rev}{score};
        };
        if( $ok ){
            if( $match{$id}{for}{seq} eq $match{$id}{rev}{seq} ){
                if( $match{$id}{for}{from} == $match{$id}{rev}{from} and
                    $match{$id}{for}{to}   == $match{$id}{rev}{to} ){
                    $perfect_count++;
                }
                else{
                    $other_1_count++;
                    # $tb->report( "Ambigous case ( seq=$id):" . Dumper $match{$id} );

                }
            }
            else{
                if( $match{$id}{for}{from} == $match{$id}{rev}{from} and
                    $match{$id}{for}{to}   == $match{$id}{rev}{to} ){
                    $good_count++;
                }
                elsif( $match{$id}{for}{to}   < $match{$id}{rev}{from} or
                       $match{$id}{for}{from} > $match{$id}{rev}{to} ){ # the two matches are fully disjoint
                       $disjoint_count++;
                }
                else{
                    $other_2_count++;
                    # $tb->report( "Ambigous case ( seq=$id):" . Dumper $match{$id} );
                }
            }
        }
        else{ # ! $ok
            $error_count++;
#             $tb->warn( "Comparison failed (seq=$id):" . Dumper $match{$id} );
        }
    }
    my $test_ok = $error_count == 0;
    $TEST_FAILED = 1 unless $test_ok;
    $tb->report( 'Ok perfect',  $perfect_count );
    $tb->report( 'Ok good',     $good_count );
    $tb->report( 'Ok disjoint', $disjoint_count );
    $tb->report( 'Ok other 1',  $other_1_count );
    $tb->report( 'Ok other 2',  $other_2_count );
    $tb->report( 'Error',       $error_count );
    $tb->report( 'test',        $test_ok ? 'PASS' : 'FAILED' );
    warn "\n";
}

exit 1 if $TEST_FAILED;
exit 0;


