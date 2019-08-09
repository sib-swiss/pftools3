package Toolbox;

use strict;

use Carp qw( longmess );
use Data::Dumper;
use IO::File;
use IO::Dir;
use Term::ANSIColor;

sub new{
    my( $package, $is_quiet ) = @_;
    my $self = bless {}, $package;
    $self->{verbose} = ! $is_quiet; # verbose by default: report actions on STDERR
    return $self;
}
sub verbose{ # turn verbosity on/off and query it
    my( $self, $verbose ) = @_;
    $self->{verbose} = $verbose  if defined $verbose;
    return $self->{verbose};
}
sub _format_msg{
    my( $self, $title, $arg, $color ) = @_;
    $arg = Dumper $arg if ref $arg;
    $arg = '' unless defined $arg;
    $color = 'black' unless $color;
    chomp $arg;
    $title = ' ' .$title  while 12 > length $title;
    return join'',
        colored( "# $title : ", 'blue' ),
        colored( $arg, $color ),
        "\n";
}
sub report{ # the "two-arguments version of warn"
    my( $self, $title, $arg, $color) = @_; # $title is expected to be a single short word
    warn $self->_format_msg( $title, $arg, $color || 'black' );
}
sub warn{
    my( $self, $arg ) = @_;
    $self->report( 'WARN', $arg, 'red');
}
sub die{
    my( $self, $arg ) = @_;
    die $self->_format_msg( 'ERROR', $arg, 'red');
    exit 1;
}
sub confess{
    my( $self, $arg ) = @_;
    CORE::die $self->_format_msg( 'ERROR', $arg . longmess(), 'red');
}
sub system{
    my( $self, $cmd ) = @_;
    $self->report( 'system', $cmd, 'magenta' ) if $self->{verbose};
    unless( system( $cmd ) == 0 ){
        $self->warn( $cmd )  unless $self->{verbose};
        # the following code was adapted from `perldoc -f system`
        if( $? == -1 ){
            $self->die( "Failed to execute: $!" );
        }
        elsif( $? & 127 ) {
            $self->die( sprintf
            'Child died with signal %d, %s coredump',
            ($? & 127),  ($? & 128) ? 'with' : 'without' );
        }
        else{
            $self->die( sprintf 'Child exited with value %d', $? >> 8 );
        }
    }
}
sub open{
    my( $self, $file_arg ) = @_; # something like 'foo.txt'; '> foot.txt'; '| sort >> foo.txt'
    my $fh = new IO::File;
    $fh->open( $file_arg )  or $self->die( "Cannot open: $file_arg" );
    $self->report( 'open', $file_arg )  if $self->{'verbose'};
    return $fh;
}
sub slurp{
    my( $self, $filename ) = @_;
    my $fh = new IO::File;
    $fh->open( $filename )  or $self->die( "Cannot slurp file: $filename" );
    $self->report( 'slurp', $filename )  if $self->{'verbose'};
    my $buf = do{ local $/; <$fh>; }; # Faster file slurping method I am aware of
    close $fh;
    return $buf;
}
sub scan_tsv{ # somehow a memory greedy method: use only on "small" file
    my( $self, $filename, $delim ) = @_;
    $delim = "\t"  unless $delim;
    my $fh = new IO::File;
    $fh->open( $filename )  or $self->die( "Cannot scan_tsv file: $filename" );
    $self->report( 'scan_tsv', $filename )  if $self->{'verbose'};
    my $buf = do{ local $/; <$fh>; }; # Maybe sysread is even faster!
    close $fh;
    chomp $buf;
    my @tab = ();
    foreach( split /\n/, $buf ){
        next  if /^\#/;
        push @tab, [ split /$delim/o ];
    }
    return \@tab;
}

1;

