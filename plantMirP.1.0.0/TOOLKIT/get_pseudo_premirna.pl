#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);


# ===============================================================================
# 
# To ensure the pseudo pre-miRNAs to be similar with known pre-miRNAs, 
# the pseudo pre-miRNA are srandomly selected from CDS sequences.
# The ratio of single- and multi-stem pseudo pre-miRNAs, length distribution, 
# minimum of base pairings in the hairpins, maximum of free energy of
# secondary structures (including GU wobble pairs) are required to 
# the same that of geuine pre-miRNAs.
# Reference: ZmirP
#
# CDS sequencs of plant can downloaded from EnsemblPlants database:
# http://plants.ensembl.org/info/website/ftp/index.html
#
# ===============================================================================





my $usage = "$0

< positive struct >
< segments struct >
< ratio of negative to positive >

";



my $infile_positive_struct = shift or die $usage;
my $infile_segment_struct = shift or die $usage;
my $ratio_negative_to_positive = shift or die $usage;



&random_select_pseudo_premirna_from_segments(
$ratio_negative_to_positive,$infile_positive_struct,
$infile_segment_struct);



print "\n\n################# END #################\n\n";



#########################################################

sub parse_struct {
	
	my $infile_struct = shift or die;
	my $hash_struct	= shift or die;
	
	unless(-e $infile_struct){
		print STDERR "Error:$infile_struct not exists\n";
		die;
	}
	
	open IN,"<$infile_struct" or die;
	
	my $id = "";
	my %hash_tmp = ();
	
	while( defined(my $line = <IN> ) ){
		chomp( $line );
		next if $line =~ /^\s*$/;	
			
		if( $line =~/>/ ){
			$id = substr($line,1);
		} elsif ( $line =~ /[ACGUTN]+/i ) {
			
			$line = uc( $line );
			$hash_tmp{$id}{seq} = $line;
				
		} elsif ( $line =~ /\.|\(|\)/  ){				
			
			my @tmp = split /\s+\(/,$line;	
			$hash_tmp{$id}{struct} = $tmp[0];
			
			my $mfe = $tmp[1];
			$mfe =~ s/\(|\)//g;			
			$hash_tmp{$id}{mfe}= $mfe;
			
		}
	}
	
	close IN;
	
	
	for my $id ( keys %hash_tmp ){
		
		my $mfe = $hash_tmp{$id}{mfe};
		my $seq = $hash_tmp{$id}{seq};
		my $struct = $hash_tmp{$id}{struct};
		
		unless( $struct =~ /\(/ ) {
			next;
		}
		
		unless( $struct =~ /\)/ ) {
			next;
		}
		
		if (length($struct) != length($seq)){
			next;
		}
		
		unless( looks_like_number( $mfe ) ) {
			next;
		}
		
		$hash_struct->{$id}{mfe} = $mfe;
		$hash_struct->{$id}{seq} = $seq;
		$hash_struct->{$id}{struct} = $struct;
		
	}

}




#########################################################


sub is_bifurcation {
	
    my $struct	= shift or die;
    
    unless ( $struct =~ /^((\.|\()+..(\.|\))+)$/ ){
    	return 1;
    } else {
		return 0; 
	}
 
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


sub shuffle {
	
	my $usage = " shuffle < reference of array >";
	my $array = shift or die $usage;
	
	for (my $j = @$array-1; $j > 0; --$j) {
  		my $i = int(rand($j));
   		my $t = $array->[$j]; $array->[$j] = $array->[$i]; $array->[$i] = $t;
	}
}



#########################################################


sub random_select_pseudo_premirna_from_segments {
	
	my $ratio_negative_to_positive = shift or die;
	my $infile_positive_struct = shift or die;
	my $infile_CDS_struct = shift or die;

	

	# ************************************************* #
	
	
	my %hash_positive_struct = ();
	&parse_struct($infile_positive_struct,\%hash_positive_struct);
	
	
	my $bifurcation_positive = 0;
	my $no_bifurcation_positive = 0;
	my $max_mfe_positive = -10000000;
	my $min_paired_base_positive = 1000000;
	my %hash_positive_single_length = ();
	my %hash_positive_multi_length = ();
	
	for my $id ( keys %hash_positive_struct ){
		
		my $seq = $hash_positive_struct{$id}{seq};
		my $struct = $hash_positive_struct{$id}{struct};
	
		# get valid sequences
		$seq = uc($seq);
		$seq =~ s/^\s*//g;
		$seq =~ s/\s*$//g; 		
		if( not ($seq =~/^([A|C|G|U]+)$/) ){
			print STDERR "Error:$id\n";
			die;
		}
	
	
		# find minimum of base pairings
		my $n_bp = 0;
		while( $struct =~ /\)/g ){
			$n_bp++;
		}
		if( $n_bp < $min_paired_base_positive ){
			$min_paired_base_positive = $n_bp;
		}
		
		# find maximum of free energy of secondary structure
		if( $hash_positive_struct{$id}{mfe} > $max_mfe_positive ){
			$max_mfe_positive = $hash_positive_struct{$id}{mfe};
		}
		
		# The ratio of single- and multi-stem pre-miRNAs
		my $lng = length($seq);
		if( &is_bifurcation($struct) ){			
			$bifurcation_positive++;
			$hash_positive_multi_length{$lng}++;				
		} else {
			$no_bifurcation_positive++;
			$hash_positive_single_length{$lng}++;
		}
		
	}
	
	my $total_positive = $bifurcation_positive + $no_bifurcation_positive;
	
	print "\n\n########### positive samples ###########\n";
	print "minimum of base pairings: $min_paired_base_positive\n";
	print "maximum of free energy:$max_mfe_positive\n";
	print "The number of single-stem:$no_bifurcation_positive\n";
	print "The number of multi-stem:$bifurcation_positive\n";
	print "The number of genuine precursors:$total_positive\n";
	print "###########################################\n\n";
	
	
	
	# ************************************************* #
	
	
	my %hash_segment_struct = ();
	&parse_struct($infile_CDS_struct,\%hash_segment_struct);
		
	
	my $bifurcation_negative = 0;
	my $no_bifurcation_negative = 0;
	my %hash_negative_single = ();
	my %hash_negative_multi = ();

	for my $id ( keys %hash_segment_struct ){
		
		my $mfe = $hash_segment_struct{$id}{mfe};
		my $seq = $hash_segment_struct{$id}{seq};
		my $struct = $hash_segment_struct{$id}{struct};
		
		# get valid sequences
		$seq = uc($seq);
		$seq =~ s/^\s*//g;
		$seq =~ s/\s*$//g; 		
		if( not ($seq =~/^([A|C|G|U]+)$/) ){
			next;
		}

		
		# Criteria 1: minimum of base pairings in the hairpins 
		my $n_bp = 0;
		while( $struct =~ /\)/g ){
			$n_bp++;
		}
		if( $n_bp < $min_paired_base_positive ){
			next;
		}
		
		# Criteria 2: maximum of free energy of secondary structures 
		if($mfe > $max_mfe_positive){
			next;
		}
			
		# The ratio of single- and multi-stem pseudo pre-miRNAs
		my $lng = length($seq);
		if( &is_bifurcation($struct) ){	
			$bifurcation_negative++;	
			$hash_negative_multi{$lng}{$id} = 1;
		} else {			
			$no_bifurcation_negative++;
			$hash_negative_single{$lng}{$id} = 1;
		}	
	}


	print "\n\n########### Candidat segments ###########\n";
	print "The number of single-stem:$no_bifurcation_negative\n";
	print "The number of multi-stem:$bifurcation_negative\n";
	print "###########################################\n\n";
	


	# ************************************************* #
	
	# randomly select multi-stem pseudo pre-miRNAs
	my @miss = ();
	my %hash_multi = ();
	for my $lng ( sort {$a<=>$b} keys %hash_positive_multi_length ){
		# Non this length 
		my $n = $hash_positive_multi_length{$lng};
		unless(exists $hash_negative_multi{$lng}){		
			push(@miss,"multi:$lng\t$n\n");
			next;
		}
		
		#get ID
		my @shuffleIDs = ();
		for my $id ( keys %{$hash_negative_multi{$lng}} ){
			push(@shuffleIDs,$id);
		}
		
		#shuffle ID
		&shuffle(\@shuffleIDs);
		
		#output sequence
		my $i = 0;
		for my $id ( @shuffleIDs ){
			$hash_multi{$id} = $hash_segment_struct{$id}{seq};
			$i++;
			# Criteria 3: length distribution	
			if( $i >= $ratio_negative_to_positive * $n ){
				last;
			}
		}
		
		# Non enough candidate
		if( $i < $ratio_negative_to_positive * $n ){
			my $r = $ratio_negative_to_positive * $n - $i;
			push(@miss,"multi:$lng\t$r\n");
		}
			
	}
	
	
	# randomly select single-stem pseudo pre-miRNAs
	my %hash_single = ();
	for my $lng ( sort {$a<=>$b} keys %hash_positive_single_length ){
		# Non this length 
		my $n = $hash_positive_single_length{$lng};
		unless(exists $hash_negative_single{$lng}){		
			push(@miss,"single:$lng\t$n\n");
			next;
		}
		
		#get ID
		my @shuffleIDs = ();
		for my $id ( keys %{$hash_negative_single{$lng}} ){
			push(@shuffleIDs,$id);
		}
		
		#shuffle ID
		&shuffle(\@shuffleIDs);
		
		#output sequence
		my $i = 0;
		for my $id ( @shuffleIDs ){
			$hash_single{$id} = $hash_segment_struct{$id}{seq};
			$i++;
			# Criteria 3: length distribution	
			if( $i >= $ratio_negative_to_positive * $n ){
				last;
			}
		}
		
		# Non enough candidate
		if( $i < $ratio_negative_to_positive * $n ){
			my $r = $ratio_negative_to_positive * $n - $i;
			push(@miss,"single:$lng\t$r\n");
		}
	}
				
	
	# ************************************ #
	
	my $bifurcation_positive_final = 0;
	my $no_bifurcation_negative_final = 0;
	# ouput selected pseudo pre-miRNAs 
	open OUT,">negative_training_data.fa" or die;
	for my $id (sort keys %hash_multi ){
		print OUT ">$id\n$hash_multi{$id}\n";
		$bifurcation_positive_final++;
	}
	for my $id (sort keys %hash_single ){
		print OUT ">$id\n$hash_single{$id}\n";
		$no_bifurcation_negative_final++;
	}
	close OUT;
	
	
	# ************************************ #
	
	open  OUT,">miss_stat.txt" or die;
	print OUT "@miss\n";
	close OUT;
	
	# ************************************ #
	
	my $total_negative = $no_bifurcation_negative_final + $bifurcation_positive_final;
	print "\n\n########### negative samples ###########\n";
	print "The number of single-stem:$no_bifurcation_negative_final\n";
	print "The number of multi-stem:$bifurcation_positive_final\n";
	print "The number of pseudo precursors:$total_negative\n";
	print "###########################################\n\n";
	
	
}		

#########################################################
