package STRUCT;
#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Basename;

use RNA; # Vienna RNA package perl interface

use TOOLKIT;



#############################################################
####################    Non-features    #####################
#############################################################


sub fold_linux_RNA {
	
	my $infile_fasta 	= shift or die;
	my $outfile_struct	= shift or die;
	
	my %hash_seq = ();
	TOOLKIT::readfasta($infile_fasta,\%hash_seq);
	
	open OUT,">$outfile_struct" or die;
	for my $id ( keys %hash_seq ){
		my $seq = $hash_seq{$id};
		$seq = uc ( $seq );
		$seq =~ s/T/U/ig;
		$seq =~ s/^\s*//g;
		$seq =~ s/\s*$//g;
		unless( not ($seq =~/^([A|C|G|U|N]+)$/) ){
			my ($struct,$mfe) = RNA::fold($seq);
			$mfe = TOOLKIT::round($mfe,6);
			print OUT ">$id\n$seq\n$struct ($mfe)\n";
		}	
	}
	close OUT;
}


#############################################################

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
		
		unless( TOOLKIT::is_numeric($mfe) ) {
			next;
		}
		
		$hash_struct->{$id}{mfe} = $mfe;
		$hash_struct->{$id}{seq} = $seq;
		$hash_struct->{$id}{struct} = $struct;
		
	}

}


#############################################################



sub basePair { 
	
	my $struct  = shift or die;
	my $hash_bp	= shift or die;
	
	my @strs = split('', $struct);
	
	my @leftBrackets = ();
	
	for my $i ( 0 .. $#strs ) {
		if ( $strs[$i] eq "(" ) {
			push @leftBrackets, $i+1;
		} elsif ( $strs[$i] eq ")" ) {
			my $ch = pop @leftBrackets;
			my $rp = $i + 1;
			$hash_bp->{ $ch } = $rp;
			$hash_bp->{ $rp } = $ch;
		}
	}
	
	if( @leftBrackets ){
		warn "unbalanced $struct\n";
	}
}



#############################################################


sub count_n_stems {
	
	my $struct = shift or die;
	
	my $n = 0;
	while( $struct =~ /(\(\(\()+/g ) {
		$n++;
	}
	
	my $m = 0;
	while( $struct =~ /(\)\)\))+/g ) {
		$m++;
	}

	my $n_stem = $n > $m ? $n : $m;
	
	return $n_stem;
	
}


	
#############################################################


sub count_n_loops {
	
	my $struct = shift or die;
	
	my %hash_bp = ();
	&basePair($struct,\%hash_bp);
	
	my @t = sort ( {$a<=>$b} keys %hash_bp );
	
	my $n_loops = 0;
	for my $i ( 1 .. $#t ){	
		if( $hash_bp{ $t[$i] } == $t[$i-1] && $hash_bp{ $t[$i-1] } == $t[$i] ) {
			$n_loops++;
		}
	}
	
	return $n_loops;
	
}



#############################################################

sub Needleman_Wunsch {
	
	my $seqA  = shift or die;
	my $seqB  = shift or die;

	
	my $MATCH     =  1; 
	my $MISMATCH  = -1;
	my $GAP       = -1; 


	my @lattice = ();
	$lattice[0][0]{score}   = 0;
	$lattice[0][0]{pointer} = "none";
	
	for(my $j = 1; $j <= length($seqA); $j++) {
	    $lattice[0][$j]{score}   = $GAP * $j;
	    $lattice[0][$j]{pointer} = "left";
	}
	
	for (my $i = 1; $i <= length($seqB); $i++) {
	    $lattice[$i][0]{score}   = $GAP * $i;
	    $lattice[$i][0]{pointer} = "up";
	}

	
	for(my $i = 1; $i <= length($seqB); $i++) {
	    for(my $j = 1; $j <= length($seqA); $j++) {
	        my ($diagonal_score, $left_score, $up_score);
   
	        my $letterA = substr($seqA, $j-1, 1);
	        my $letterB = substr($seqB, $i-1, 1);     
	                               
	        if ($letterA eq $letterB) {
	            $diagonal_score = $lattice[$i-1][$j-1]{score} + $MATCH;
	         } else {
	            $diagonal_score = $lattice[$i-1][$j-1]{score} + $MISMATCH;
	         }

	        
	        $up_score   = $lattice[$i-1][$j]{score} + $GAP;
	        $left_score = $lattice[$i][$j-1]{score} + $GAP;

	        
	        if ($diagonal_score >= $up_score) {
	            if ($diagonal_score >= $left_score) {
	                $lattice[$i][$j]{score}   = $diagonal_score;
	                $lattice[$i][$j]{pointer} = "diagonal";
	             } else {
	                $lattice[$i][$j]{score}   = $left_score;
	                $lattice[$i][$j]{pointer} = "left";
	             }
	         } else {
	            if ($up_score >= $left_score) {
	                $lattice[$i][$j]{score}   = $up_score;
	                $lattice[$i][$j]{pointer} = "up";
	             } else {
	                $lattice[$i][$j]{score}   = $left_score;
	                $lattice[$i][$j]{pointer} = "left";
	             }
	         }
	     }
	}


	my $alignA = "";
	my $alignB = "";

	my $j = length($seqA);
	my $i = length($seqB);

	while ( 1 ) {
		
	    last if $lattice[$i][$j]{pointer} eq "none"; 

	    if ($lattice[$i][$j]{pointer} eq "diagonal") {
	        $alignA .= substr($seqA, $j-1, 1);
	        $alignB .= substr($seqB, $i-1, 1);
	        $i--;
	        $j--;
	     } elsif ($lattice[$i][$j]{pointer} eq "left") {
	        $alignA .= substr($seqA, $j-1, 1);
	        $alignB .= "-";
	        $j--;
	     } elsif ($lattice[$i][$j]{pointer} eq "up") {
	        $alignA .= "-";
	        $alignB .= substr($seqB, $i-1, 1);
	        $i--;
	     }    
	     
	}

	$alignA = reverse $alignA;
	$alignB = reverse $alignB;
	
	return ($alignA,$alignB);
}






#############################################################
#############   triplet-SVM features    #####################
#############################################################


# =============================
# The local contiguous structure-sequence features
# Reference: triplet-SVM
# =============================

sub triplet {
	
	my $seq = shift or die;
	my $struct = shift or die;
	my $triplet_array = shift or die;
	my $triplet_array_str = shift or die;
	
	my $i = 0;
	my $k_mer = 3;
	my %triplet = ();
	
	if( not ($seq =~/^([A|C|G|T|U|a|c|g|t|u]+)$/) ){
		print "Error: The $seq contains characters others than [acgtuACGTU]\n";
		die;
	}
	
	
	while ( $i < (length $seq) ) {
		
		my $subseq = substr( $seq, $i, $k_mer );
		my $substruct = substr( $struct, $i, $k_mer );
		$i++;
		
		my $ch = substr( $subseq, 1, 1 );
		
		my $str = $ch."$substruct";
		$str =~ s/\)/\(/g;
		
		next if length( $str ) < $k_mer + 1;
		
		unless( exists $triplet{ $str } ){
			$triplet{ $str } = 1;
		} else {
			$triplet{ $str }++;
		}	
	}
	
	# ********************************** #
	
	my $k = 0;
	my @chars = qw( A C G U );
	my @bracket_dot = ( "(((", "((.", "(..", "(.(", ".((", ".(.", "..(", "..." );
	
	my @allTriplet = ();
	for my $base ( @chars ) {
		for my $i ( 0 ..$#bracket_dot ){
			$allTriplet[$k++] = $base.$bracket_dot[$i];
		}
	}
	
	my $ntriplet = 0;
	for my $triplet ( @allTriplet ){
		unless( exists $triplet{ $triplet } ) {
			$triplet{ $triplet } = 0;
		}
		$ntriplet++;
	}
	
	# ********************************** #
	
	my $j = 0;
	for my $triplet ( keys %triplet ) {
		$triplet_array->[$j] = $triplet{$triplet}/$ntriplet;
		$triplet_array_str->[$j] = $triplet;
		$j++;
	}
	
	if ( $j != ( 2 ** $k_mer ) * 4 ){
		print "Error:lack triplet ($seq\t$struct)\n";
		die;
	}
	
}






#############################################################
#################    miRD features    #######################
#############################################################





# =============================
# Calculate the distribution of 
# unpaired base to paired base
# Reference:miRD (Zhang)
# =============================

sub paired_unpaired_base {
	
	my $struct = shift or die;
	
	my $sum = 0;
	while( $struct =~ /\(|\)|\./g ) {
		$sum++;
	}
	
	my $pair = 0;
	while( $struct =~ /\(|\)/g ) {
		$pair++;
	}
	
	my $unpair = 0;
	while( $struct =~ /\./g ) {
		$unpair++;
	}
	
	my $pair2unpairRatio = 0;
	if( $unpair ){
		$pair2unpairRatio = $pair/$unpair;
	}
	
	return ($pair/$sum,$unpair/$sum,$pair2unpairRatio);
	
}




# =============================
# Calculate the distribution of 
# unpaired base to paired base
#Reference:miRD (Zhang)
# =============================

sub unpaired2paired_distribution {
	
	my $struct 	= shift or die;
	my $nbin 	= shift or die;
	my $dis		= shift or die;
	
	
	my $bin = (length $struct) / $nbin;	
	if( $bin - int( $bin ) < 1.e-10 ){
		$bin = int ( $bin );
	} else {
		$bin = int( $bin ) + 1;
	}
	
	my $i = 0;
	my $j = 0;
	
	while ( $i < (length $struct) ) {
		
		my $substruct = substr( $struct, $i, $bin );
				
		my $pair = 0;
		while( $substruct =~ /\(|\)/g ) {
			$pair++;
		}
		
		my $unpair = 0;
		while( $substruct =~ /\./g ) {
			$unpair++;
		}	
		
		my $ratio = 0;
		if( $pair ){
			$ratio = $unpair / $pair;	
		}
		
		$dis->[$j] = $ratio;
		$i += $bin;
		$j++;
		
	}	
	
	while ( $j < $nbin  ){
		$dis->[$j++] = 0;
	}	
	
}




# =============================
# reverse tructure strings and character strings, then 
# aligne the strings with their reversed counterparts.
# the gaps between these strings were filled with "-"

# the aligned strings were divided into ten equal regions,
# and  the ratios of each nucleotide pair type were calculated as
# input featuresaps between these strings were filled with "-"
# 
# Reference:miRD (Zhang)
# =============================

sub reverse_align_bp {
	
	my $seq = shift or die;
	my $struct = shift or die;
	my $n_segment = shift or die;
	my $hash_rabp = shift or die;
	
	
	my $rev_seq = "";
	my $rev_struct = "";
	my @seqs = split ('',$seq);
	my @structs = split ('',$struct);

	# "reverse" sequence and struct.
	for( my $i = $#seqs; $i >= 0; $i-- ){
		$rev_seq .= $seqs[$i];	
		if ( $structs[$i] eq ")" ){
			$rev_struct .= "(";
		} else {
			$rev_struct .= $structs[$i];
		}
	}
	

	# change struct string into sequence string,
	# then excute Needleman Wunsch alignement.
	$struct =~ tr/\)/\(/;
	$struct =~ tr/\(\./AT/;
	$rev_struct =~ tr/\(\./AT/;
	($struct, $rev_struct) = &Needleman_Wunsch( $struct, $rev_struct );
	$struct =~ tr/AT/\(\./;
	$rev_struct =~ tr/AT/\(\./;


	# get the sequence corresponding to struct.
	@seqs = split('',$seq);
	@structs = split('',$struct);
	
	my $j = 0;
	my $result_seq = "";
	
	foreach my $i ( 0 .. $#structs ){
		if ( $structs[$i] eq "-" ){
			$result_seq .= "-";
		} else {
			$result_seq .= $seqs[$j];
			$j++;
		}
	}
	
	
	my @rev_seqs = split('',$rev_seq);
	my @rev_structs = split('',$rev_struct);
	
	$j = 0;
	my $result_rev_seq = "";
	
	foreach my $i ( 0 .. $#rev_structs ){
		if ( $rev_structs[$i] eq "-" ){
			$result_rev_seq .= "-";
		} else {
			$result_rev_seq .= $rev_seqs[$j];
			$j++;
		}
	}
	
	
	
	# ****************************************** #
	
	
	my $strA = uc ( $result_seq );
	my $strB = uc ( $result_rev_seq );
	$strA =~ tr/T/U/;
	$strB =~ tr/T/U/;
	
	
	# all possible combination between A C G U -
	my @bases = qw(A C G U -);	
	my %hash_all_bp = ();
	for my $i ( @bases ) {
		for my $j ( @bases ) {
			if ( $i ne "-" or $j ne "-" ){
				my @tmp = ($i, $j);
				@tmp = sort ( @tmp );
				my $bp = join('',@tmp);
				unless( exists $hash_all_bp{ $bp } ){
					$hash_all_bp{ $bp } = 1;
				} 
			}
		}
	}
	

	if ( length ( $strA ) ne length( $strB ) ){
		print STDERR "Error:the length of two strings is different\n";
		die;
	}
	
	my $string_lng = length ( $strA );
	my $segment_lng = int ( $string_lng / $n_segment ) + 1;
	
	my @arrayA = split ('',$strA);
	my @arrayB = split ('',$strB);
	
	my %hash_total = ();
	for(my $i = 0; $i <= $#arrayA; $i++ ){
		my $t = int ( $i / $segment_lng );

		my $chA = $arrayA[$i];
		my $chB = $arrayB[$i];
		my @tmp = ($chA,$chB);
		@tmp = sort (@tmp);
		my $bp = join('',@tmp);		
		unless ( exists $hash_rabp->{$t}{$bp} ){
			$hash_rabp->{$t}{$bp} = 1;
		} else {
			$hash_rabp->{$t}{$bp}++;
		}
		unless ( exists $hash_total{$bp} ){
			$hash_total{$bp} = 1;
		} else {
			$hash_total{$bp}++;
		}
	}
	

	for my $i (0 ..$n_segment-1){
		for my $j ( sort keys %hash_all_bp ){
			unless ( exists $hash_rabp->{$i}{$j} ){
				$hash_rabp->{$i}{$j} = 0;
			} 
			unless ( exists $hash_total{$j} ){
				$hash_total{$j} = 1;
			}
		}
	}
	
	
	# normalization 
	my $sum = 0;
	for my $i ( sort {$a<=>$b} keys %$hash_rabp ){
		for my $j ( sort  keys %{ $hash_rabp->{ $i } } ){
			$hash_rabp->{ $i }{ $j } /= $hash_total{ $j };
			$sum++;
		}
	}
	
}



#############################################################
#################   miPred features   #######################
#############################################################



# =============================
# %(G+C) = (|C|+|G|)/L * 100
#Reference:miPred
# =============================

sub GC_content {	
	
	my $seq = shift or die;
	
	my @t = split ( '', $seq );	
	
	my $m = 0;
	my $n = 0;
	
	for my $base ( @t ) {
		if ( $base eq "G" or $base eq "g" ) {
			$m++;
		} elsif ( $base eq "C" or $base eq "c" ) {
			$m++;
		}
		$n++;
	}
	
	return ($m/$n) * 100;
	
}



# =============================
# dinucleotide frequencies
# %XY = |XY| * 100 / ( L - 1 )
# where X and Y belong to (A C G U)
#Reference:miPred
# =============================

sub dinucleotide_frequency {
	
	my $seq = shift or die;
	my $dinucleotide = shift or die;
	my $dinucleotide_str = shift or die;
	
	$seq = uc($seq);
	$seq =~ tr/T/U/;
	
	my $n = 0;
	my @tmp = ();
	my @bases = qw ( A C G U );
	for my $i ( @bases ){
		for my $j ( @bases ){
			$tmp[$n++] = $i.$j;
		}
	}
		
	my $i = 0;
	for my $di ( @tmp ){
		my $n = 0;
		while( $seq =~ /$di/g ) {
			$n++;
		}
		$dinucleotide->[$i] = ($n/(length($seq)-1)) * 100;
		$dinucleotide_str->[$i] = $di;
		$i++;
	}
	
	if ( $i != 16 ){
		print STDERR "Error: lack dinucleotide\n";
		die;
	}
	
}




# =============================
# normalized minimum free energy of folding( dG ):
# dG = MFE / L
# Reference:miPred 
# =============================

sub normalized_mfe {
	
	my $struct = shift or die;
	my $mfe = shift or die;
	
	return $mfe/(length($struct));
	
}




# =============================
# MFE Index 1 = dG / %(G+C)
# Reference:miPred
# =============================

sub mfe_index_1 {
	
	my $seq = shift or die;
	my $mfe = shift or die;
	
	my $dG = &normalized_mfe($seq,$mfe);
	
	my $GCc = &GC_content($seq);
	
	if( $GCc ) {
		return $dG / $GCc;
	} else {
		return 0;
	}
	
}



# =============================
# MFE Index 2 = dG / n_stems
# Reference:miPred
# =============================

sub mfe_index_2 {
	
	my $struct = shift or die;
	my $mfe = shift or die;
	
	my $dG = &normalized_mfe($struct,$mfe);
	
	my $n_stems = &count_n_stems($struct);
	
	if( $n_stems ){
		return $dG / $n_stems;
	} else {
		return 0;
	}
	
}



# =============================
# Normalized base-paring propensity (dP)
# dP = tot_bases / L
# Reference:miPred
# =============================

sub normalized_base_paring_propensity {
	
	my $struct = shift or die;
	
	my $tot_bases = 0;
	while( $struct =~ /\)/g ){
		$tot_bases++;
	}
	
	return $tot_bases /(length $struct);
	
}




#############################################################
################## microPred features #######################
#############################################################



# =============================
# MFE Index 3 = dG / n_loops
# Reference:microPred
# =============================

sub mfe_index_3 {
	
	my $struct = shift or die;
	my $mfe = shift or die;
	
	my $dG = &normalized_mfe($struct,$mfe);

	my $n_loops = &count_n_loops($struct);
	
	if( $n_loops ){
		return $dG / $n_loops;
	} else {
		return 0;
	}
}
	



# =============================
# Normalized base pair counts:
# |A-U|/L, |G-C|/L and |G-U|/L
# where |X-Y| is the number of (X-Y) base pairs 
# in the secondary structure
# 
# AND
#
# %(A-U)/n_stems, %(G-C)/n_stems and %(G-U)/n_stems
# where %(X-Y) = |X-Y|/tot_bases. 
# tot_bases is the total number of base pairs
# Reference:microPred
# =============================

sub normalized_base_pair_count {
	
	my $seq = shift or die;
	my $struct = shift or die;
	
	my $ref_array_nbp = shift or die;
	my $ref_array_nbp_stem = shift or die;
	
	my $ref_array_nbp_str = shift or die;
	my $ref_array_nbp_stem_str = shift or die;
	


	$seq = uc ( $seq );
	$seq =~ tr /T/U/;
	
	my %hash_bp = ();
	&basePair($struct,\%hash_bp);
	
	
	my %hash_tmp = ();
	for my $i ( sort {$a<=>$b} keys %hash_bp ){
		my $ch_a = substr($seq,$i-1,1);
		my $ch_b = substr($seq,$hash_bp{$i}-1,1);
		my @tmp = ();
		push(@tmp,$ch_a);
		push(@tmp,$ch_b);
		@tmp = sort {$a cmp $b} @tmp;
		my $ch = join('',@tmp);
		unless( exists $hash_tmp{$ch} ){
			$hash_tmp{$ch} = 1;
		} else {
			$hash_tmp{$ch}++;
		}
	}
	
	my %bp_tmp = ();
	$bp_tmp{AU} = 1;
	$bp_tmp{CG} = 1;
	$bp_tmp{GU} = 1;

	for my $bp ( sort keys %bp_tmp ){
		unless ( exists $hash_tmp{ $bp } ){
			$hash_tmp{ $bp } = 0;
		}
	}
	
	
	my $n = 0;
	for my $bp ( sort keys %hash_tmp ){
		$ref_array_nbp->[$n] = $hash_tmp{ $bp } / (length $seq);
		$ref_array_nbp_str->[$n] = $bp;
		$n++;
	}
	
	
	if ( $n != 3 ){
		print STDERR "Error:lack base pair counts ($n,$seq)\n";
		for my $bp ( sort keys %hash_tmp ){
			print STDERR "$bp\n";
		}
		die;
	}
	
	
	# ************************************** #
	
	my $tot_bases = 0;
	while( $struct =~ /\)/g ){
		$tot_bases++;
	}
	
	my $n_stems = &count_n_stems($struct);
	
	$n = 0;
	for my $bp ( sort keys %hash_tmp ){
		if( $n_stems ){
			$ref_array_nbp_stem->[$n] = ($hash_tmp{ $bp } / $tot_bases) / $n_stems;
		} else {
			$ref_array_nbp_stem->[$n] = 0;
		}
		$ref_array_nbp_stem_str->[$n] = $bp;
		$n++;
	}
	
	# ************************************** #
}




# =============================
# Average base pairs per stem:
# avg_bp_stem = tot_bases / n_stems
# Reference: microPred
# =============================

sub avg_bp_stem {
	
	my $struct = shift or die;
	
	my $tot_bases = 0;
	while( $struct =~ /\)/g ){
		$tot_bases++;
	}
	
	my $n_stems = &count_n_stems($struct);
	
	if( $n_stems ){
		return $tot_bases / $n_stems;
	} else {
		return 0;
	}
	
}




#############################################################
#################### ZmirP features #########################
#############################################################




# =============================
# the size of biggest bulge
# Reference: ZmirP
# =============================

sub biggest_bulge_size {
	
	my $struct = shift or die;
	
	my $max_bulge_size = 0;
	while ( $struct =~ /(\.+)/g ) {
		my $size = length $1;
		if ( $size > $max_bulge_size ) {
			$max_bulge_size = $size;
		}
	}
	
	return $max_bulge_size;
	
}




# =============================
# maximum number of consecutive paired nucleotides
# Reference:ZmirP
# =============================

sub biggest_consecutive_paired_base {
	
	my $struct = shift or die;
	
	my %hash_bp = ();
	&basePair($struct,\%hash_bp);
	
	my $max_consecutive_lng = 0;
	for my $i ( sort {$a<=>$b} keys %hash_bp ){
		my $j = $i;
		while ( $j <= (length $struct) ){
			if ( defined $hash_bp{$j} ){
				$j++;
			} else {
				last;
			}
		}
		if ( $max_consecutive_lng < $j - $i ){
			$max_consecutive_lng = $j - $i;
		}
	}	
	
	return $max_consecutive_lng;
	
}




# =============================
# the number of bulges
# Reference: ZmirP
# =============================

sub count_n_bulges {
	
	my $struct = shift or die;
	
	my $n = 0;
	while( $struct =~ /(\.\.\.)+/g ) {
		$n++;
	}
	
	return $n;
	
}



# =============================
# the normalized number of stems
# normalized_n_stems = n_stems / L
# Reference: ZmirP
# =============================

sub normalized_n_stems {
	
	my $struct = shift or die;
	
	my $n_stems = &count_n_stems($struct);
	
	return $n_stems/(length $struct);
	
}




# =============================
# the normalized number of bulges
# normalized_n_bulges = n_bulges / L
# Reference: ZmirP
# =============================

sub normalized_n_bulges {
	
	my $struct = shift or die;
	
	my $n_bulges = &count_n_bulges($struct);
	
	return $n_bulges/(length $struct);
	
}




# =============================
# the normalized number of loops
# normalized_n_loops = n_loops / L
# Reference: ZmirP
# =============================

sub normalized_n_loops {
	
	my $struct = shift or die;
	
	my $n_loops = &count_n_loops($struct);
	
	return $n_loops/(length $struct);
	
}




# =============================
# Average unpaired bases per bulge:
# avg_unbp_bulge = tot_unpaired_bases / n_bulges
# where tot_unpaired_bases is the total number of unpaired bases
# Reference: ZmirP
# =============================

sub avg_unbp_bulge {
	
	my $struct = shift or die;
	
	my $tot_unpaired_bases = 0;
	while( $struct =~ /\./g ){
		$tot_unpaired_bases++;
	}
	
	my $n_bulges = &count_n_bulges($struct);
	
	if( $n_bulges ){
		return $tot_unpaired_bases / $n_bulges;
	} else {
		return 0;
	}
	
}




# =============================
# MFE Index 4 = dG / tot_bases
# where tot_bases is the total number of base pairs
# Reference:ZmirP
# =============================

sub mfe_index_4 {
	
	my $struct = shift or die;
	my $mfe = shift or die;
	
	my $tot_bases = 0;
	while( $struct =~ /\)/g ){
		$tot_bases++;
	}
	

	my $dG = &normalized_mfe($struct,$mfe);
	
	if( $tot_bases ){
		return $dG / $tot_bases;
	} else {
		return 0;
	}
	
}




# =============================
# MFE Index 5 = dG / n_bulges
# Reference: ZmirP
# =============================

sub mfe_index_5 {
	
	my $struct = shift or die;
	my $mfe = shift or die;
	
	my $dG = &normalized_mfe($struct,$mfe);
	
	my $n_bulges = &count_n_bulges($struct);
	
	if( $n_bulges ){
		return $dG / $n_bulges;
	} else {
		return 0;
	}
	
}


#############################################################
1;