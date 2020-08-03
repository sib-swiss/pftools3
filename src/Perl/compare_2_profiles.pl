#!/usr/bin/env perl

# Perl embedded modules
use strict;
use warnings;
use diagnostics;

# Command line used:
# pfcalibrate -m pfscale --filter-db /db/shuffled/window20.seq --heuristic-db profile    -M 1 -o $OUT/$shortname $prf -r $OUT/${shortname/.prf/.pdf} --verbose 2>$OUT/${shortname/.prf/.log} --nthreads 4

use File::Slurp;

# Get profiles
my $old_prf = $ARGV[0]  or die "\n\t$0 old_prf new_prf\n\n";
my $new_prf = $ARGV[1]  or die "\n\t$0 old_prf new_prf\n\n";

if ( !-e "$old_prf" || !-r "$old_prf" || -z "$old_prf" ){
    die "\n\t[$old_prf] cannot be read\n\n";
}
if ( !-e "$new_prf" || !-r "$new_prf" || -z "$new_prf" ){
    die "\n\t[$new_prf] cannot be read\n\n";
}


# Read profiles
my $old_param = read_prf("$old_prf");
my $new_param = read_prf("$new_prf");


# Check parameters
if ( $old_param->{'ID'} eq '' || $old_param->{'ID'} ne $new_param->{'ID'} ){
    die "\n\tIDs are different between both profiles [$old_param->{'ID'}] [$new_param->{'ID'}]\n\n";
}
if ( $old_param->{'DT'} eq '' || $old_param->{'DT'} != $new_param->{'DT'} ){
    die "\n\tDates are different between both profiles [$old_param->{'DT'}] [$new_param->{'DT'}]\n\n";
}
CHECK:
for my $p ( keys %{ $old_param } ){
    if ( !defined $old_param->{$p} || !defined $new_param->{$p} ){
        warn "\tParsing problem for [$old_param->{$p}] or [$new_param->{$p}] in $p\n";
    }
    if ( ($old_param->{$p} eq '' && $new_param->{$p} ne '') || ($old_param->{$p} ne '' && $new_param->{$p} eq '') ){
        warn "\tCannot compare parameters $p for $old_param->{'ID'} ($old_param->{'DT'}): [$old_param->{$p}] [$new_param->{$p}]\n";
        exit 1;     # Skip the entry for stat
        last CHECK; # Stop at the 1st failed test
    }
}


# Write result
print join("\t", '#ID', 'DATE', 'OLD:M1_R1', 'OLD:M1_R2', 'OLD:M-1_R1', 'OLD:M-1_R2', 'OLD:L1_score', 'OLD:L0_score', 'OLD:L-1_score',
                                'NEW:M1_R1', 'NEW:M1_R2', 'NEW:M-1_R1', 'NEW:M-1_R2', 'NEW:L1_score', 'NEW:L0_score', 'NEW:L-1_score'), "\n";
print join("\t", $old_param->{'ID'}, $old_param->{'DT'}, $old_param->{'M1_R1'}, $old_param->{'M1_R2'},
                                                         $old_param->{'M-1_R1'}, $old_param->{'M-1_R2'},
                                                         $old_param->{'L1_score'}, $old_param->{'L0_score'}, $old_param->{'L-1_score'},
                                                         $new_param->{'M1_R1'}, $new_param->{'M1_R2'},
                                                         $new_param->{'M-1_R1'}, $new_param->{'M-1_R2'},
                                                         $new_param->{'L1_score'}, $new_param->{'L0_score'}, $new_param->{'L-1_score'},
          ), "\n";

exit 0;


sub read_prf {
    my ($prf_file) = @_;

    my $param;
    $param->{'ID'}        = '';
    $param->{'DT'}        = '';
    $param->{'M1_R1'}     = '';
    $param->{'M1_R2'}     = '';
    $param->{'M-1_R1'}    = '';
    $param->{'M-1_R2'}    = '';
    $param->{'L1_score'}  = '';
    $param->{'L0_score'}  = '';
    $param->{'L-1_score'} = '';
    PRF:
    for my $line ( read_file("$prf_file", chomp =>1) ){
        last PRF  if ( $line =~ /^MA   \/DEFAULT:/ ); # To speed up the parsing

        # ID   AddA; MATRIX.
        if ( $line =~ /^ID   ([^;]+);/ ){
            $param->{'ID'} = $1;
        }
        # DT   DEC-2014 (DATA UPDATE).
        elsif ( $line =~ /^DT   [A-Z]{3}\-(\d{4}) / ){
            $param->{'DT'} = $1;
        }
        # MA   /NORMALIZATION: MODE=1; FUNCTION=LINEAR; R1=7.4199743; R2=0.0003973; TEXT='-LogE';
        elsif ( $line =~ /^MA   \/NORMALIZATION: MODE=1; FUNCTION=LINEAR; R1=([^;]+); R2=([^;]+);/ ){
            $param->{'M1_R1'} = &is_float($1) ? $1 : undef;
            $param->{'M1_R2'} = &is_float($2) ? $2 : undef;
        }
        # MA   /NORMALIZATION: MODE=-1; FUNCTION=LINEAR; R1=-1288983.5; R2=117.584267; TEXT='Heuristic';
        elsif ( $line =~ /^MA   \/NORMALIZATION: MODE=-1; FUNCTION=LINEAR; R1=([^;]+); R2=([^;]+);/ ){
            $param->{'M-1_R1'} = &is_float($1) ? $1 : undef;
            $param->{'M-1_R2'} = &is_float($2) ? $2 : undef;
        }
        # MA   /CUT_OFF: LEVEL=1; SCORE=19660; H_SCORE=1022723; N_SCORE=15.23; MODE=1; TEXT='!';
        elsif ( $line =~ /^MA   \/CUT_OFF: LEVEL=1; SCORE=([^;]+);/ ){
            $param->{'L1_score'} = &is_int($1) ? $1 : undef;
        }
        # MA   /CUT_OFF: LEVEL=0; SCORE=19660; H_SCORE=1022723; N_SCORE=15.23; MODE=1; TEXT='?';
        elsif ( $line =~ /^MA   \/CUT_OFF: LEVEL=0; SCORE=([^;]+);/ ){
            $param->{'L0_score'} = &is_int($1) ? $1 : undef;
        }
        # MA   /CUT_OFF: LEVEL=-1; SCORE=2719; H_SCORE=-969271; N_SCORE=8.5; MODE=1; TEXT='??';
        elsif ( $line =~ /^MA   \/CUT_OFF: LEVEL=\-1; SCORE=([^;]+);/ ){
            $param->{'L-1_score'} = &is_int($1) ? $1 : undef;
        }
    }

    return $param;
}

sub is_float {
    @_ == 1 || croak(q/Usage: is_float(string)/);
    local $_ = $_[0];
    return ( defined && /\A-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?\z/ );
}

sub is_int {
    @_ == 1 || croak(q/Usage: is_int(string)/);
    local $_ = $_[0];
    return ( defined && /\A-?(?:0|[1-9][0-9]*)\z/ );
}
