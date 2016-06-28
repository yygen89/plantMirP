#!/usr/bin/perl
use strict;
use warnings;



my $fasta = shift or die;
my $prefix = shift or die;


my %hash_seq = ();
&readfasta($fasta,\%hash_seq);


my $max_length = 0;
my $n = 1;
my $total = 0;
my $miss = 0;
open SEQ,">valid-$fasta" or die;
for my $id (sort keys %hash_seq ){
	my $seq = $hash_seq{$id};

	$seq = uc($seq);
	$seq =~ s/T/U/ig;
	
	if(length($seq) > $max_length){
		$max_length = length($seq);
	}
	
	unless(not ($seq =~ /^([A|C|G|U]+)$/)){
		print SEQ ">$prefix$n\n$seq\n";
		$n++;
	}else{
		$miss++;
	}
	
	$total++;
}
close SEQ;
		
print "\nMax length:$max_length\n";
print "\nTotal number:$total\n";
print "\nMissed number:$miss\n";


#########################################################


sub readfasta {
	
	my $usage = "< fasta file > < hash reference >";
	my $infile = shift or die $usage;
	my $hash_seq = shift or die $usage;
	
	unless(-e $infile){
		print STDERR "Error:$infile not exists";
		die;
	}
	
	open IN, $infile || die;
	
	my %hash_check_id = ();
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
			
			my $new_id = $seqId;	
			if($new_id =~ /(\S+)/){
				$new_id = $1;
			}
			$seqId = $new_id;
			
			unless(exists $hash_check_id{$seqId}){
				$hash_check_id{$seqId} = 1;
			}else{
				print STDERR "Error: $seqId repeat\n";
				die;
			}
			
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}