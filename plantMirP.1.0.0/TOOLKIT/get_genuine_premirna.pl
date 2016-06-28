#!/usr/bin/perl
use strict;
use warnings;




my @related_organisms = qw(ath gma osa ppt mtr sbi aly zma sly);
my $miRBase_hairpin_fasta = "hairpin.fa";



my %hash_seq = ();
&readfasta($miRBase_hairpin_fasta,\%hash_seq);



open OUT,">genuine_premirna.fa" or die;
for my $id ( sort keys %hash_seq ){
	for my $related_organism ( @related_organisms ){
		if($id =~/^$related_organism\b/i){
			my $seq = $hash_seq{$id};
			unless(not ($seq =~ /^([A|C|G|U]+)$/)){
				print OUT ">$id\n$seq\n";
			}
		}
	}
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


########################################################

sub trim {
	my $str = shift or die;
	$str =~ s/^\s*//g;
	$str =~ s/\s*$//g;
	return $str;
}


