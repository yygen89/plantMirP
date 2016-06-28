#!/usr/bin/perl
use strict;
use warnings;


my $usage = "$0 < precusor FASTA > < ratio of training > ";


my $precursor_fasta = shift or die $usage;
my $ratio_training = shift or die $usage;



my %hash_seq = ();
&readfasta($precursor_fasta,\%hash_seq);



my %hash_species = ();
my %hash_species_ID = ();
for my $id ( sort keys %hash_seq ){
	my @tmp = split/-/,$id;
	my $species = shift @tmp;
	$hash_species{$species}++;
	push(@{$hash_species_ID{$species}},$id);
}


my $ratio_test = 1 - $ratio_training;

open OUT,">statistics.txt" or die;
open TEST,">precursor_test-$ratio_test.fa" or die;
open TRAINING,">precursor_training-$ratio_training.fa" or die;

print OUT "Species\tTotal\tTraining\tTest\n";
for my $species ( sort keys %hash_species ){
	
	my @shuffleIDs = @{$hash_species_ID{$species}};
	&shuffle(\@shuffleIDs);
	
	my $n = 0;
	my $m = 0;
	my $r = 0;
	my $total = $hash_species{$species};
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
	
	print OUT "$species\t$total\t$r\t$m\n";
}

close OUT;
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