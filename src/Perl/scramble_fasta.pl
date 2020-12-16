#!/usr/bin/env perl

use strict;

use List::Util qw( shuffle ); # as a replacement for Math::Random::random_permutation

#--------------------------------------------------------------------#
# How to                                                             #
#--------------------------------------------------------------------#

my $usage='
usage:

   scramble_fasta.pl [options] fasta-file...

This script produces"random" or "shuffled" sequences according to a
model and to a set of parameters. All sequences are read and written
in FASTA format. No distinction is made between DNA and protein
sequences. All letter are converted to uppercase. Gaps are
removed. FASTA headers are modified. The module Math:Random from CPAN
was originally to generate random numbers. This version was rewritten
for portability (i.e. no dependency).

The mandatory option "-m <string>" allows to select a shuffling
method. The meanings and the numerical values of the other parameters
then depend on the chosen method (The default values are given between
brackets).

   scramble_fasta.pl -m synthetic [-P 100 -F sprot34]

produces a sequences of length "-P" with the amino-acid frequencies set
according to "-F". Use "-M" to produce several sequences.

   scramble_fasta.pl -m permutation fasta-file

permutates the amino acids of each sequences in the source file.

   scramble_fasta.pl -m window [-P 20] fasta-file

permutates the amino acids along the sequence in consecutive window of
length "-P".

   scramble_fasta.pl -m codon [-P 20] fasta-file

permutates the codons (three consecutive letters) along the sequence in
consecutive window of length "-P" codons, i.e. 3*P letters.

   scramble_fasta.pl -m scramble fasta-file

behaves like "synthetic" but takes the sequence count and lengths from
fasta-file.

   scramble_fasta.pl -m corrupt [-P 50 -F sprot34] fasta-file

corrupts the amino acids of each sequences in the source
file. Corruption consist in the replacement of randomly sampled
residue with a random amino acid which frequencies obey to "-F". The
elementary corrpution events is performed "-P" times per 100
residues. Note than non integer "-P" are supported.

   scramble_fasta.pl -m bubble [-P 50] fasta-file

swaps randomly chosen pair of adjascent amino acids.  The elementary
swap event is carried on P times every 100 residues. Note that the
sequence is considered periodic: the last amino acid can be swaped with
the first one.

   scramble_fasta.pl -m reverse fasta-file

reverses the sequences.

   scramble_fasta.pl -m sample [-P 50] fasta-file

extracts a random sample of length "-P" from each sequence of fasta-file.

   scramble_fasta.pl -m rubble [-P 10] fasta-file

brokes the sequence in piece of at most "-P" residue long and then
reassamble the sequence at random.

   scramble_fasta.pl -m spoil [-P 10] fasta-file

spoils the amino acids of each sequences in the source file. At each
spoil cycle a residue is removed at a randomly choosen position and a
new residue is inserted at another randomly choosen position. The
elementary spoil events is performed "-P" times per 100 residues. Note
than non integer "-P" are supported.

brokes the sequence in piece of at most "-P" residue long and then
reassamble the sequence at random.

   scramble_fasta.pl -m comb [-P 2 -Q 0] fasta-file

mutates every "-P" residue, with an offset of -Q, with a randomly
selected residue according to frequency "-F". The residue is
however garantied to be changed.

   scramble_fasta.pl -m subset [-P 10] fasta-file

randomly select "-P" sequences off fasta-file.

More options:

    -h This help message.
    -s <int> Seed of the random number generator (current time by default).
    -F = [sprot34|Robinson|Altshul|Dayhoff|Pro-rich|DNA|source], predefined amino acids frequencies, (sprot34 by default).
    -M <integer> multiplier, repeat M times the procedure without reseting the seed of the random number generator.


Example:

   scramble_fasta.pl -m synthetic  -P 50 -s 0 \
   | scramble_fasta.pl -m bubble -P 20 -M 8 -s "an apple a day" - \
   | readseq -p -fMSF

';

####################################################
# sprot34 amino acids frequency                    #
####################################################

my %sprot34_freq=( # B X Z were  removed
    'A' => 0.07554,
    'C' => 0.01698,
    'D' => 0.05305,
    'E' => 0.06322,
    'F' => 0.04077,
    'G' => 0.06846,
    'H' => 0.02241,
    'I' => 0.05730,
    'K' => 0.05941,
    'L' => 0.09342,
    'M' => 0.02357,
    'N' => 0.04530,
    'P' => 0.04927,
    'Q' => 0.04024,
    'R' => 0.05158,
    'S' => 0.07223,
    'T' => 0.05747,
    'V' => 0.06527,
    'W' => 0.01252,
    'Y' => 0.03199
    );

####################################################
# The following amino acids frequency were stolen  #
# in the file "blastkar.c" from the ncbi blast 2.0 #
# distribution                                     #
####################################################

#/*  M. O. Dayhoff amino acid background frequencies   */
my %dayhoff_freq=(
    'A' => 0.08713,
    'C' => 0.03347,
    'D' => 0.04687,
    'E' => 0.04953,
    'F' => 0.03977,
    'G' => 0.08861,
    'H' => 0.03362,
    'I' => 0.03689,
    'K' => 0.08048,
    'L' => 0.08536,
    'M' => 0.01475,
    'N' => 0.04043,
    'P' => 0.05068,
    'Q' => 0.03826,
    'R' => 0.04090,
    'S' => 0.06958,
    'T' => 0.05854,
    'V' => 0.06472,
    'W' => 0.01049,
    'Y' => 0.02992,
    );
# Stephen Altschul amino acid background frequencies #
my %altschul_freq=(
    'A' => 0.08100,
    'C' => 0.01500,
    'D' => 0.05400,
    'E' => 0.06100,
    'F' => 0.04000,
    'G' => 0.06800,
    'H' => 0.02200,
    'I' => 0.05700,
    'K' => 0.05600,
    'L' => 0.09300,
    'M' => 0.02500,
    'N' => 0.04500,
    'P' => 0.04900,
    'Q' => 0.03900,
    'R' => 0.05700,
    'S' => 0.06800,
    'T' => 0.05800,
    'V' => 0.06700,
    'W' => 0.01300,
    'Y' => 0.03200,
    );
# /* amino acid background frequencies from Robinson and Robinson */
my %robinson_freq=(
    'A' => 0.07805,
    'C' => 0.01925,
    'D' => 0.05364,
    'E' => 0.06295,
    'F' => 0.03856,
    'G' => 0.07377,
    'H' => 0.02199,
    'I' => 0.05142,
    'K' => 0.05744,
    'L' => 0.09019,
    'M' => 0.02243,
    'N' => 0.04487,
    'P' => 0.05203,
    'Q' => 0.04264,
    'R' => 0.05129,
    'S' => 0.07120,
    'T' => 0.05841,
    'V' => 0.06441,
    'W' => 0.01330,
    'Y' => 0.03216
    );
my %dna_freq=(
	      'A' => 0.25,
	      'C' => 0.25,
	      'G' => 0.25,
	      'T' => 0.25
);

######################
# "global" variables #
######################

my @id;
my %seq;
my %freq;
my @alphabet;
my %cumul_freq;
my %opt;

########################
# process command line #
########################

require Getopt::Std;
Getopt::Std::getopts('hm:P:Q:F:o:s:M:',\%opt);
die $usage if $opt{'h'} or ! $opt{m};
$opt{'o'}='seq' unless $opt{'o'};
die "wrong output option \"$opt{'o'}\"" unless $opt{'o'}=~/^(seq|freq)$/;
$opt{'M'}=1 unless $opt{'M'};
if( exists $opt{s} ){
    die "Random seed of is not an integer: $opt{s}\n" unless $opt{s} =~ /^\d+$/;
    warn "Seed: $opt{s}\n";
    srand( $opt{s} );
}


################
# Program body #
################

#####################################
# seed  the random number generator #
#####################################

# require Math::Random;
# if(defined $opt{'s'}){
#     Math::Random::random_set_seed_from_phrase($opt{'s'});
# }

# main of main
if($opt{'m'} eq 'synthetic'){
    $opt{'P'}=100 unless $opt{'P'};
    $opt{'F'}='sprot34' unless $opt{'F'};
    doOptF();
}
elsif($opt{'m'} eq 'permutation'){}# do nothing in this case
elsif($opt{'m'} eq 'window'){
    $opt{'P'}=20 unless $opt{'P'};
}
elsif($opt{'m'} eq 'codon'){
    $opt{'P'}=20 unless $opt{'P'};
}
elsif($opt{'m'} eq 'scramble'){
    $opt{'F'}='sprot34' unless $opt{'F'};
    doOptF();
}
elsif($opt{'m'} eq 'corrupt'){
    $opt{'P'}=50 unless $opt{'P'};
    $opt{'F'}='sprot34' unless $opt{'F'};
    doOptF();
}
elsif($opt{'m'} eq 'bubble'){
    $opt{'P'}=50 unless $opt{'P'};
}
elsif($opt{'m'} eq 'reverse'){}
elsif($opt{'m'} eq 'sample'){
    $opt{'P'}=50 unless $opt{'P'};
}
elsif($opt{'m'} eq 'rubble'){
    $opt{'P'}=10 unless $opt{'P'};
}
elsif($opt{'m'} eq 'spoil'){
    $opt{'P'}=10 unless $opt{'P'};
    $opt{'F'}='sprot34' unless $opt{'F'};
    doOptF();
}
elsif($opt{'m'} eq 'comb'){
    $opt{'Q'}=0 unless $opt{'Q'};
    $opt{'P'}=2 unless $opt{'P'};
    $opt{'F'}='sprot34' unless $opt{'F'};
    doOptF();
}
elsif($opt{'m'} eq 'subset'){
  $opt{'P'}=10 unless $opt{'P'};
}
else{
    die "unknown or missing model \"$opt{'m'}\"\n";
}if($opt{'o'}=~/^freq$/){
    foreach(@alphabet){
#	printf "%s %.5f %.5f %.5f %.5f\n",$_,$sprot34_freq{$_},$robinson_freq{$_},$altschul_freq{$_},$dayhoff_freq{$_};
	printf "%s %.5f %.5f\n",$_,$freq{$_},$cumul_freq{$_};
    }
}
else{ # means $opt{'o'} eq 'seq'
    my($header,$seq,$i);
    if($opt{'m'} eq 'synthetic'){
	for($i=0;$i<$opt{'M'};$i++){
	    $header="synthetic_".($i+1);
	    $seq=randseq($opt{'P'},\@alphabet,\%cumul_freq);
	    print toFasta($header,$seq);
	}
    }
    elsif($opt{'m'} eq 'permutation'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_permut_".($i+1);
		# $seq=join('',Math::Random::random_permutation(split(//,$seq{$_})));
		$seq = join( '', shuffle( split( //, $seq{$_} )));
		print toFasta($header,$seq);
	    }
	}
    }
    elsif($opt{'m'} eq 'window'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my @seq;
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_window$opt{'P'}_".($i+1);
		@seq=$seq{$_}=~/(.{1,$opt{'P'}})/g;
		foreach(@seq){
		    # $_=join('',Math::Random::random_permutation(split //)); # cool !
		    $_ = join( '', shuffle split // );
		}
		print toFasta($header,join('',@seq));
	    }
	}
    }
    elsif($opt{'m'} eq 'codon'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my(@seq,@buf);
	my $len=3*$opt{'P'};
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_codon$opt{'P'}_".($i+1);
		@seq=$seq{$_}=~/(.{1,$len})/g;
		foreach(@seq){
		    @buf=$_=~/(.{3})/g;
		    # $_=join('',Math::Random::random_permutation(@buf));
		    $_ = join( '', shuffle @buf );
		}
		print toFasta($header,join('',@seq));
	    }
	}
    }
    elsif($opt{'m'} eq 'scramble'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_scramble_".($i+1);
		$seq=randseq(length $seq{$_},\@alphabet,\%cumul_freq);
		print toFasta($header,$seq);
	    }
	}
    }
    elsif($opt{'m'} eq 'corrupt'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my($mutation_count,$length,$mutation);
	my @seq;
	my @pos;
	my @mutation_site;
	my @rand;
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_corrupt$opt{'P'}_".($i+1);
		$length=length($seq{$_});
		$mutation_count=int(($length*$opt{'P'})/100);
		@seq=split(//,$seq{$_});
		# @pos=Math::Random::random_uniform_integer($mutation_count,0,$length-1);
		@pos = random_uniform_integer( $mutation_count, 0, $length-1 );
		foreach(@pos){
		    $mutation = randalpha( \@alphabet, \%cumul_freq );
		    $seq[$_]=$mutation;
		}
		print toFasta($header,join('',@seq));
	    }
	}
    }
    elsif($opt{'m'} eq 'bubble'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my($mutation_count,$length,$mutation);
	my @seq;
	my @pos;
	my @mutation_site;
	my @rand;
	for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_bubble$opt{'P'}_".($i+1);
		$length=length($seq{$_});
		$mutation_count=int(($length*$opt{'P'})/100);
		@seq=split(//,$seq{$_});
		# @pos=Math::Random::random_uniform_integer($mutation_count,0,$length-1);
		@pos=random_uniform_integer( $mutation_count, 0, $length-1 );
		($seq[$_-1],$seq[$_])=($seq[$_],$seq[$_-1]) foreach @pos;
		print toFasta($header,join('',@seq));
	    }
	}
    }
    elsif($opt{'m'} eq 'reverse'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my @residue;
        foreach(@id){
	   $header=$_."_reverse";
           @residue=split(//,$seq{$_});
           print toFasta($header,join('',reverse(@residue)))
       }
    }
    elsif($opt{'m'} eq 'sample'){
	readSequences(); # read $ARGV[1] into %seq and @id;
        for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		# my $pos=Math::Random::random_uniform_integer(1,0,length($seq{$_})-$opt{'P'});
		my $pos = random_uniform_integer( 1, 0, length($seq{$_}) - $opt{'P'} );
		$header=$_.'/'.($pos+1).'-'.($pos+$opt{'P'});
		print toFasta($header,substr($seq{$_},$pos,$opt{'P'}));
	    }
	}
    }
    elsif($opt{'m'} eq 'rubble'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my @block;
        for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_rubble$opt{'P'}_".($i+1);
		@block=$seq{$_}=~/(.{$opt{'P'}})/g;
		# $seq=join('',Math::Random::random_permutation(@block));
		$seq = join( '', shuffle @block );
		print toFasta($header,$seq);
	    }
	}
    }
    elsif($opt{'m'} eq 'spoil'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my @residue;
        for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_spoil$opt{'P'}_".($i+1);
		@residue=split //,$seq{$_};
		my $len=@residue;
		my $mutation_count=int(($len*$opt{'P'})/100);
		foreach(1..$mutation_count){
                    # my $new_residue=splice @residue,Math::Random::random_uniform_integer(1,0,$len-1),1;
                    # splice @residue,Math::Random::random_uniform_integer(1,0,$len-2),0,$new_residue;
		    my $new_residue = splice @residue, random_uniform_integer( 1, 0, $len - 1 ), 1;
		    splice @residue, random_uniform_integer( 1, 0, $len - 2 ), 0, $new_residue;
		}
		print toFasta($header,join('',@residue));
	    }
	}
    }
    elsif($opt{'m'} eq 'comb'){
	readSequences(); # read $ARGV[1] into %seq and @id;
	my @residue;
        for($i=0;$i<$opt{'M'};$i++){
	    foreach(@id){
		$header=$_."_comb_$opt{'P'}_$opt{'Q'}_".($i+1);
		my $rand;
		@residue=split //,$seq{$_};
		foreach(0..@residue-1){
		    if(($_ % $opt{'P'})==$opt{'Q'}){
			do{
			    $rand=randalpha(\@alphabet,\%cumul_freq);
			}until $rand ne $residue[$_];
			$residue[$_]=lc $rand;
		    }
		}
		print toFasta($header,join('',@residue));
	    }
	}
    }
    elsif($opt{'m'} eq 'subset'){
        readSequences(); # read $ARGV[1] into %seq and @id;
        die "-P <integer> is too large\n" if $opt{'P'}> @id;
        my %rnd;
	my $count=0;
        while($opt{'P'}>$count){
	    # my $x=Math::Random::random_uniform_integer(1,0,@id-1);
	    my $x = random_uniform_integer( 1, 0, @id - 1 );
	    unless($rnd{$x}){
		$rnd{$x}=1;
		$count++;
	    }
        }
        foreach(sort {$a<=>$b} sort keys %rnd){
	    print toFasta($id[$_],$seq{$id[$_]});
        }
    }
}

########################
# Subs are handy tools #
########################

sub doOptF{
    unless($opt{'F' }=~/^(sprot34|Robinson|Altschul|Dayhoff|Pro-rich|DNA|source)$/){
	die "unknow frequency option \"$opt{'F'}\"";
    }
    if($opt{'F'} eq 'sprot34'){%freq=%sprot34_freq}
    elsif($opt{'F'} eq 'Robinson'){%freq=%robinson_freq}
    elsif($opt{'F'} eq 'Altschul'){%freq=%altschul_freq}
    elsif($opt{'F'} eq  'Dayhoff'){%freq=%dayhoff_freq}
    elsif($opt{'F'} eq 'Pro-rich'){
	%freq=%sprot34_freq;
	$freq{$_}*=9/10 foreach sort keys %freq;
	$freq{'P'}+=1/10;
    }
    elsif($opt{'F'} eq 'DNA'){%freq=%dna_freq}

    elsif($opt{'F'} eq 'source'){
	readSequences();
	freqencyFromSequences();
    }
    updateFrequency();
}

########################
# read the source file #
########################

sub readSequences{
    @ARGV=('-') unless @ARGV;
    die "missing filename" unless @ARGV==1;
    foreach my $filename (@ARGV){
	open(FASTA,$filename)or die "cannot open \"$filename\" for reading";
	my($id,$seq);
	while(<FASTA>){
	    if(/^>(\S+)/){
		if($id){
		    push @id,$id;
		    $seq{$id}=$seq;
		}
		$id=$1;
		$seq='';
	    }
	    else{
		s/\W//g;        # remove blanks and gaps
		$seq.=uc $_;    # uppercase only
	    }
	}
	if($id){
	    push @id,$id;
	    $seq{$id}=$seq;
	}
    }
}
##############################
# to save shuffled sequences #
##############################

sub toFasta{
    my($header,$seq)=@_;
    $seq=~s/(.{1,80})/$1\n/g;
    ">$header\n$seq";
}

#####################################
# subs to process letters frequency #
# used by the synthetic model       #
#####################################

sub freqencyFromSequences{
    my $id;
    foreach $id(@id){
	foreach(split(//,$seq{$id})){
	    $freq{$_}++; # unless $_=~/[BXZ]/;
	}
    }
    my $total;
    $total+=$freq{$_} foreach sort keys %freq;
    $freq{$_}/=$total foreach sort keys %freq;
    updateFrequency();
}
sub updateFrequency{
    my $total=0;
    @alphabet = sort keys %freq;
    my $sum=0;
    foreach(@alphabet){
	$sum+=$freq{$_};
	$cumul_freq{$_}=$sum;
    }
}
sub randseq{
    my($length,$alphabet_ref,$cumul_freq_ref)=@_;
    my $seq;
    for(my $i=0;$i<$length;$i++){
	$seq.=randalpha($alphabet_ref,$cumul_freq_ref);
    }
    $seq;
}
sub randalpha{
    my($alphabet_ref,$cumul_freq_ref)=@_;
    # my $random=Math::Random::random_uniform();
    my $random = rand();
    foreach(@{$alphabet_ref}){
	return $_ if $random<$cumul_freq_ref->{$_};
    }
    warn "Random letter generator does not work: random=$random\n";
    return $alphabet_ref->[-1];
}

sub random_uniform_integer{ # this is a standalone replacement for the same method in Math::Random
    my( $n, $low, $high ) = @_;
    my $delta = $high - $low + 1 ;
    my @res = ();
    foreach( 1 .. $n ){
        push @res, $low + rand( $delta );
    }
    return @res;
}


