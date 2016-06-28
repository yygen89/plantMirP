#!/usr/bin/perl
use strict;
use warnings;



my $usage = "$0 

< file in FASTA >
< Number of sequences >
< outfile >

";


my $fasta_file 	= shift or die $usage;
my $number 		= shift or die $usage;
my $outfile 	= shift or die $usage;


my %hash_seq = ();
&readfasta($fasta_file,\%hash_seq);


my @shuffleIDs = keys %hash_seq;
&shuffle(\@shuffleIDs);


my $i = 0;
open OUT,">$outfile" or die;
for my $id ( @shuffleIDs ){
	$i++;
	print OUT ">$id\n$hash_seq{$id}\n";
	last if $i >= $number;
}
close OUT;



########################################################

sub readfasta {
	
	my $usage = "< fasta file > < hash reference >";
	my $infile = shift or die $usage;
	my $hash_seq = shift or die $usage;
	
	unless(-e $infile){
		print STDERR "Error:$infile not exists";
		die;
	}
	
	open IN, $infile || die;
	
	my $c=0;
	my $seqId;
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}


########################################################


sub shuffle {
	
	my $usage = " shuffle < reference of array >";
	my $array = shift or die $usage;
	
	for (my $j = @$array-1; $j > 0; --$j) {
  		my $i = int(rand($j));
   		my $t = $array->[$j]; $array->[$j] = $array->[$i]; $array->[$i] = $t;
	}
}

########################################################
