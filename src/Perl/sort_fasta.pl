#!/bin/env perl

use strict;
use Getopt::Std;

# ------------------------------------------------------
# Read command line and check its arguments
# ------------------------------------------------------
my $usage = "$0 [options] <fasta-input>

options: -h       this help
         -R       reverse sequence
         -G       remove gap
         -U       uppercase
         -C       complement
	     -s       simplify header: keep ID and coordinates
         -S       simplify header: keep ID only
";
my %opt;  # global to store the options
Getopt::Std::getopts('hRGUCsS',\%opt);
if( $opt{'h'} or @ARGV < 1){
    print "$usage\n";
    exit 0;
}

# ------------------------------------------------------
#
# ------------------------------------------------------

my %complement = (
    A => 'T',
    C => 'G',
    G => 'C',
    T => 'A',

    R => 'Y', # A or G => C or T
    Y => 'R', # C or T
    S => 'S', # G or C => G or C
    W => 'W', # A or T => A or T
    K => 'M', # G or T => C or A
    M => 'K', # A or C

    B => 'V', # C or G or T => G or C or A = !T
    D => 'H', # A or G or T => T or C or A = !G
    H => 'D', # A or C or T => T or G or A = !C
    V => 'B', # A or C or G => T or G or C = !A

    N => 'N', # any base

);

# ------------------------------------------------------
# Do the job
# ------------------------------------------------------

my $header   = '';
my %fasta = ();
while( <> ){
	if( /^>/ ){
        $header = $_;
    }
    else{
        $fasta{$header} .= $_;
    }
}
foreach my $header ( sort keys %fasta ){
    my $seq = $fasta{$header};
    $header =~ s/>(\S+).+/>$1/     if $opt{s};
    $header =~ s/>([[^\/]+).+/>$1/ if $opt{S};
    print $header;
    $seq =~ s/\s+//g;
    $seq = join '', reverse split //, $seq if $opt{R};
    $seq =~ s/[^A-z]//g if $opt{G};
    $seq = uc $seq if $opt{U};
    if( $opt{C} ){
        my @seq = ();
        foreach( split //, $seq ){
            if( /[A-Z]/ ){
                push @seq, $complement{$_} || 'N';
            }
            elsif( /[a-z]/ ){
                push @seq, lc( $complement{ uc $_ } ) || 'n';
            }
            else{
                push @seq, $_;
            }
        }
        $seq = join '', @seq;
    }
    $seq =~ s/(.{1,60})/$1\n/g;
    print $seq;
}

