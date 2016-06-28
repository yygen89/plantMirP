#!/usr/bin/perl
use strict;
use warnings;


my %hash_negative_seq = ();
&readfasta("negative_training_data.fa",\%hash_negative_seq);



my %hash_segment_struct = ();
&readfasta("CDS_segment_struct.txt",\%hash_segment_struct);



open OUT,">test_CDS_segment_struct.txt" or die;
for my $id ( sort keys %hash_segment_struct ){
	unless(exists $hash_negative_seq{$id}){
		my $s = $hash_segment_struct{$id};
		chomp($s);
		print OUT ">$id\n$s\n";
	}
}
close OUT;




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
	
	my $c=0;
	my $seqId;
	while (defined (my $line=<IN>)) {
		#chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
		} else {
			#$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}
