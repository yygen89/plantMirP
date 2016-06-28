package KnowledgeBasedEnergy;
#!/usr/bin/perl
use strict;
use warnings;


use TOOLKIT;


####################################################

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



####################################################

sub reverse_align_bp {
	
	my $seq = shift or die;
	my $struct = shift or die;
	
	
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
	

	if ( length ( $strA ) ne length( $strB ) ){
		print STDERR "Error:the length of two strings is different\n";
		die;
	}
	
	return($strA,$strB);
	
}



####################################################

sub kmer_composition_frequency {
	
	my $bin = shift or die;
	my $r_bin = shift or die;
	my $hash_struct = shift or die;
	my $hash_dist = shift or die;
	my $hash_base_compositions = shift or die;
	
	
	my $k_mer = 3;
	my $hash_potential;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/k_mer/i) {$k_mer=shift @_}
		if ($argument=~/potential/i) {$hash_potential=shift @_}
	}
	
	
	
	# **************************************** #
	
	
	my $score = 0;
	for my $id ( sort keys %{$hash_struct} ){
		
		my $seq = $hash_struct->{$id}{seq};
		my $struct = $hash_struct->{$id}{struct};
		
		if( not ($seq =~/^([A|C|G|T|U|a|c|g|t|u]+)$/) ){
			print "Error: The $seq contains characters others than [acgtuACGTU]\n";
			die;
		}
		
		
		my ($strA,$strB) = &reverse_align_bp($seq,$struct);
		
		
		my @tmp1 = split //,$strA;
		my @tmp2 = split //,$strB;
		
		
		$seq = "";
		for(my $i = 0; $i < length($strA); $i++){
			$seq .= "$tmp1[$i]$tmp2[$i]";
		}
		
		
		# coordinate
		my %lattice = ();
		
		my $iii = 0;
		while ( $iii < (length $seq) ) {
			
			my $subseq = substr( $seq, $iii, $k_mer );
		
			last if length( $subseq ) < $k_mer;

			# coordinate
			$lattice{$iii}{"x"} = $iii;
			$lattice{$iii}{"base"} = $subseq;
			
			$iii++;
		}
		
		
		
		# the number of pair (i,j) at the distance r	
		for(my $i = 0; $i < $iii; $i++) {
			for(my $j = $i + 1; $j < $iii; $j++) {
				
				my $r = abs($lattice{$i}{"x"}-$lattice{$j}{"x"});				
				my $b = int( $r / $r_bin );
				$b = $bin if $b > $bin;
							

				my $base_i = $lattice{$i}{"base"};
				my $base_j = $lattice{$j}{"base"};
				
					
				my $base_ij = "$base_i*$base_j";
				$hash_base_compositions->{$base_ij} = 1;		
				
				
				if( defined $hash_potential ){
					
					if ( exists $hash_potential->{$base_ij}{$b} ){
						my $u = $hash_potential->{$base_ij}{$b};
						$score += $u;
					} 
					
				}else{
					
					unless(exists $hash_dist->{$base_ij}{$b}){
						$hash_dist->{$base_ij}{$b} = 1;
					} else {
						$hash_dist->{$base_ij}{$b}++;
					}	
				}
				
			}
		}

	}
	
	return $score;
}
	

sub statistical_potentials {
	
	my $bin = shift or die;
	my $r_bin = shift or die;
	my $positive_struct = shift or die;
	my $negative_struct = shift or die;
	my $hash_potentials = shift or die;
	
	
	my $k_mer = 3;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/k_mer/i) {$k_mer=shift @_}
	}
	
	
	
	# ================= Positive =========================
	
	
	my %hash_struct = ();
	&parse_struct($positive_struct,\%hash_struct);
	
	my %hash_pos_dist = ();
	my %hash_base_compositions = ();
	&kmer_composition_frequency($bin,$r_bin,\%hash_struct,
	\%hash_pos_dist,\%hash_base_compositions,'k_mer'=>$k_mer);
	
	
	
	# ================= Negative =========================
	
	%hash_struct = ();
	&parse_struct($negative_struct,\%hash_struct);
	
	my %hash_neg_dist = ();
	&kmer_composition_frequency($bin,$r_bin,\%hash_struct,
	\%hash_neg_dist,\%hash_base_compositions,'k_mer'=>$k_mer);
	
	
	# ==========================================
	
	my $epsilon = 0.00001;
	for my $base_ij ( sort keys %hash_base_compositions ){
	for my $b ( 0 .. $bin ){
		
		my $numerator = 1;
		if (exists $hash_pos_dist{$base_ij}{$b}){
			$numerator = $hash_pos_dist{$base_ij}{$b};
		}
		
		my $denominator = 1;
		if (exists $hash_neg_dist{$base_ij}{$b}){
			$denominator = $hash_neg_dist{$base_ij}{$b};
		}
		
		if( ($numerator > 1) or ($denominator > 1) ){
			my $r = -1 * log ( $numerator / $denominator );
			$r = TOOLKIT::round($r,6);
			
			if ( abs($r) >= $epsilon ){
				unless(exists $hash_potentials->{$base_ij}{$b}){
					$hash_potentials->{$base_ij}{$b} = $r;
				} else {
					print STDERR "Error:($base_ij,$b) repeat\n";
					die;
				}
			}
		}
		
	}}
	
	# ==========================================	
}



#############################################################



sub energy_score {
	
	my $bin = shift or die;
	my $r_bin = shift or die;
	my $hash_struct = shift or die;
	my $hash_potential = shift or die;


	my $k_mer = 3;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/k_mer/i) {$k_mer=shift @_}
	}
	

	my %hash_tmp1 = ();
	my %hash_tmp2 = ();
	my $score = &kmer_composition_frequency($bin,$r_bin,$hash_struct,
	\%hash_tmp1,\%hash_tmp2,'k_mer'=>$k_mer,'potential'=>$hash_potential);
	
	
	return $score;

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


sub is_numeric {
	
	use Scalar::Util qw(looks_like_number);
	
    my $v = shift;

    if( looks_like_number( $v ) ){  	
        return 1;
    } else {
        return 0;
    }
    
}


####################################################
1;


