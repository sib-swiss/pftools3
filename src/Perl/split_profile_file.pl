#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

use File::Slurp;

my $prf_file = $ARGV[0]  or die "\n\t$0 hamap.prf\n\n";
open(my $PRF, '<', "$prf_file")  or die "\n\tCannot open [$prf_file]";
my $prf = '';
my $ac  = '';
while(<$PRF>){
    if ( $_ =~ /^\/\/$/ ){
        $prf .= $_;
        write_file("$ac.prf", $prf);
    }
    #AC   MF_00001;
    elsif ( $_ =~ /^ID   / ){
        $prf = $_;
    }
    elsif ( $_ =~ /^AC   ([^;]+);/ ){
        $ac   = $1;
        $prf .= $_;
    }
    else {
        $prf .= $_;
    }
}
close $PRF;

exit 0;

