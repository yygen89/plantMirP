#!/usr/bin/perl
use strict;
use warnings;


my $usage = "$0 < negative FASTA > < ratio of training > ";


my $negative_fasta = shift or die $usage;
my $ratio_training = shift or die $usage;



my %hash_seq = ();
my $total = &readfasta($negative_fasta,\%hash_seq);


my $ratio_test = 1 - $ratio_training;
my @shuffleIDs = keys %hash_seq;
&shuffle(\@shuffleIDs);


open TEST,">negative_test-$ratio_test.fa" or die;
open TRAINING,">negative_training-$ratio_training.fa" or die;

my $n = 0;
my $m = 0;
my $r = 0;	
for my $id ( @shuffleIDs ){
		
	my $seq = $hash_seq{$id};
	$seq = uc($seq);
	$seq =~ s/T/U/ig;
	
	unless(not ($seq =~ /^([A|C|G|U]+)$/)){
			
		if( $n > $ratio_training * $total ){
			$m++;
			print TEST ">$id\n$seq\n";
		}else{
			$r++;
			print TRAINING ">$id\n$seq\n";
		}
		$n++;
	}
}

close TEST;
close TRAINING;
			
			
		





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