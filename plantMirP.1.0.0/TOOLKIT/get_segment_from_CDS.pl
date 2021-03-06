#!/usr/bin/perl
use strict;
use warnings;



my $usage = "$0

<positive sequences in FASTA>
<CDS sequences in FASTA>

";



my $positve_fasta = shift or die $usage;
my $CDS_fasta = shift or die $usage;


&get_segment_from_CDS($positve_fasta,$CDS_fasta);


print "\n\n################# END #################\n\n";



#########################################################

sub get_segment_from_CDS {
	
	my $infile_positive_fasta = shift or die;
	my $infile_CDS_fasta = shift or die;


	# get length distribution of known geuine pre-miRNAs
	my %hash_positive_seq = ();
	&readfasta($infile_positive_fasta,\%hash_positive_seq);
	
	my %hash_positive_length = ();
	for my $id ( keys %hash_positive_seq ){
		my $lng = length($hash_positive_seq{$id});
		$hash_positive_length{$lng}++;
	}
		
	
	# ************************************************* #


	my %hash_seq = ();
	&readfasta($infile_CDS_fasta,\%hash_seq);
	
	
	my %hash_segment_seq = ();
	for my $id ( keys %hash_seq ){

		# non-overlapped segments under a constraint
		# condition that the length distribution of extracted 
		# segments was identical with that of known geuine pre-miRNAs
		
		my $long_negative_seq = $hash_seq{$id};
		
		my $i = 0;
		my $exit_loop = 0;
		while ( $i < length($long_negative_seq) ){
			for my $j ( sort {$a<=>$b} keys %hash_positive_length ){	
				if( ($i + $j) >= length($long_negative_seq) - 1 ){
					$exit_loop = 1;
					last;
				}	
				
				my $subseq = substr($long_negative_seq, $i, $j);		
				$subseq = uc($subseq);
				$subseq =~ tr/T/U/;
				$subseq =~ s/^\s*//g;
				$subseq =~ s/\s*$//g; 
				unless( not ($subseq =~/^([A|C|G|U]+)$/) ){
					unless(exists $hash_segment_seq{$subseq} ){
						$hash_segment_seq{$subseq} = 1;
					}
				}
				
				$i += $j;			
			}
		
			if( $exit_loop ){
				last;
			}
		}
	}
	
	my $n = 1;
	my $segment_fasta = "CDS_segment_seq.fa";
	open OUT,">$segment_fasta" or die;	
	for my $seq ( keys %hash_segment_seq ){
		print OUT ">randomSeqFrom_CDS$n\n$seq\n";
		$n++;
	}	
	close OUT;

}


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
	my %hash_id_check = ();
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
			
			unless(exists $hash_id_check{$seqId}){
				$hash_id_check{$seqId} = 1;
			}else{
				print STDERR "Error:$seqId repeat\n";
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

#########################################################
