#!/bin/env perl

use strict;
use Getopt::Std;

# ------------------------------------------------------
# Read command line and check its arguments
# ------------------------------------------------------
my $command = join ' ', $0, @ARGV;
my $usage = "$0 [options] <fasta-input>

options: -h       this help
         -M <float>  match score    (50 by default)
         -m <float> mismatch score (-40 by default)
";
my %opt;  # global to store the options
Getopt::Std::getopts('hM:m:',\%opt);
if( $opt{'h'} or @ARGV ){
    print "$usage\n";
    exit 0;
}
$opt{M} =  50 unless exists $opt{M};
$opt{m} = -40 unless exists $opt{m};

# ------------------------------------------------------
# Define the IUPAC code
# ------------------------------------------------------

my @alpha = (
    A => 'A',
    C => 'C',
    G => 'G',
    T => 'T',    # need to consider U ?

    R => 'AG',   # A or G
    Y => 'CT',   # C or T
    S => 'GC',   # G or C
    W => 'AT',   # A or T
    K => 'GT',   # G or T
    M => 'AC',   # A or C

    B => 'CGT',  # C or G or T
    D => 'AGT',  # A or G or T
    H => 'ACT',  # A or C or T
    V => 'ACG',  # A or C or G

    N => 'ACGT', # any base
);


# ------------------------------------------------------
# Do the job
# ------------------------------------------------------
# print "#\n";
# print "# this file was created with: $command\n";
# print "#\n";
foreach my $i ( 0 .. int( @alpha / 2 ) - 1 ){
    print "     $alpha[ 2 * $i ]";
}
print " ..\n";

foreach my $i ( 0 .. int( @alpha / 2 ) - 1 ){
    my $symb_1 = $alpha[ 2 * $i ];
    my $list_1 = $alpha[ 2 * $i + 1 ];
    foreach my $j ( 0 .. int( @alpha / 2 ) - 1 ){
        my $symb_2 = $alpha[ 2 * $i ];
        my $list_2 = $alpha[ 2 * $j + 1 ];
        my $score = 0;
        my $count = 0;
        foreach my $x ( split //, $list_1 ){
            foreach my $y ( split //, $list_2 ){
                $score += $x eq $y ? $opt{M} : $opt{m};
                $count++;
            }
        }
        if( $j < $i ){
            print ' ' x 6;
        }
        else{
            my $txt = $score > 0 ? int( $score / $count + 0.5 ) : int( $score / $count - 0.5 );
            $txt = ' ' . $txt while 6 > length $txt;
            print $txt;
        }
    }
    print " $symb_1\n";
}


