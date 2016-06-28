package FEATURE_SELECTION;
#!/usr/bin/perl
use strict;
use warnings;
use File::Copy;
use File::Basename;

use TOOLKIT;



#############################################################

# Reference: Combining SVMs with various feature selection strategies
# Reference: Efficient feature selection filters for high-dimensional data

sub feature_seletcion {
	
	my $infile_feature = shift or die;
	my $outfile_selected_feature = shift or die;
	my $outfile_keep_feature_name = shift or die;
	my $keep_feature_num = shift or die;
	my $relevance_measure_type 	= shift or die;
	my $max_similarity = shift;
	
	
	unless(-e $infile_feature){
		print STDERR "Error:$infile_feature not exists\n";
		die;
	}
	
	
	my %hash_row = ();
	my @cols = TOOLKIT::toMatrix($infile_feature,\%hash_row);
	
	
	my @rows = ();
	my %hash_type = ();
	my $type_col_index = (scalar @cols) - 1;
	
	for my $row ( sort {$a cmp $b} keys %hash_row ){
		if ( $row !~/^\bName\b/i ){
			push(@rows,$row);
		}
		my $label = $hash_row{$row}{ $cols[$type_col_index] };
		$hash_type{$row} = $label;
	}

	
	my $feature_num = 0;
	my %hash_relevance = ();
	
	
	# ******** Compute the relevance of each feature ********
	
	for my $col ( 1 .. ($type_col_index-1) ){	 #features

		my $whole_sum 	 = 0;
		my $positive_sum = 0;
		my $negative_sum = 0;
		my $positive_num = 0;
		my $negative_num = 0;
		$feature_num++;
		
		
		my $exp_positive_sum = 0;
		my $exp_negative_sum = 0;
	
		my @positives = ();
		my @negatives = ();
		
		for my $row ( @rows ){		
			
			my $v = $hash_row{$row}{ $cols[$col] };
			my $label = $hash_row{$row}{ $cols[$type_col_index] };
			
			unless ( TOOLKIT::is_numeric($v) ){
				print STDERR "Error:$v not a numeric ($row,$cols[$col])\n";
				die;
			}
			
			$whole_sum += $v;
			if ( $label eq "positive" ){		
				$positive_sum += $v;
				$positive_num++;	
				$exp_positive_sum += exp( $v );
				push(@positives,$v);
			} elsif ( $label eq "negative" ){		
				$negative_sum += $v;
				$negative_num++;
				$exp_negative_sum += exp( $v );	
				push(@negatives,$v);
			} else {
				print STDERR "Error:($row,$cols[$col],$v)\n";
				die;
			}
			
		}
		
		
			
		my $av_positive = $positive_sum / $positive_num;
		my $av_negative = $negative_sum / $negative_num;
			
		if ( $relevance_measure_type eq "AMGM" ){
			
			# =============== AMGM ==============
						
			my $AMGM_positive = $exp_positive_sum / ( $positive_num *  exp($av_positive) );
			my $AMGM_negative = $exp_negative_sum / ( $negative_num * exp($av_negative) );
		
			unless(exists $hash_relevance{$cols[$col]}){
				$hash_relevance{$cols[$col]} = 0.5 * ($AMGM_positive + $AMGM_negative);
			} else {
				print STDERR "Error:$cols[$col] repeat\n";
				die;
			}
		
		} elsif ( $relevance_measure_type eq "MM" ){
		
			# ============= MM ================
			
			my $median_positive	= TOOLKIT::median(\@positives);
			my $median_negative	= TOOLKIT::median(\@negatives);
				
			my $MM_positive = abs( $av_positive - $median_positive );
			my $MM_negative = abs( $av_negative - $median_negative );
			
			unless(exists $hash_relevance{$cols[$col]}){
				$hash_relevance{$cols[$col]} = 0.5 * ($MM_positive + $MM_negative);
			} else {
				print STDERR "Error:$cols[$col] repeat\n";
				die;
			}
			
		} else {
			
			# ============= F-score ================
		
			my $whole_av 	= $whole_sum / ($positive_num + $negative_num);
			my $positive_av = $positive_sum / $positive_num;
			my $negative_av = $negative_sum / $negative_num;
			
			my $numerator = ( $positive_av - $whole_av ) * ( $positive_av - $whole_av );
			$numerator   += ( $negative_av - $whole_av ) * ( $negative_av - $whole_av );
			
			my $positive_variance_sum = 0;
			my $negative_variance_sum = 0;
		
			for my $row ( @rows ){
				
				my $v = $hash_row{$row}{ $cols[$col] };
				my $label = $hash_row{$row}{ $cols[$type_col_index] };
				
				if ( $label eq "positive" ){			
					my $t = ( $v - $positive_av ) * ( $v - $positive_av );
					$positive_variance_sum += $t;	
				} elsif ( $label eq "negative" ){			
					my $t = ( $v - $negative_av ) * ( $v - $negative_av );
					$negative_variance_sum += $t;	
				} else {
					print STDERR "Error:($row,$label)\n";
					die;
				}
			}
			
			my $f_score = 0;	
			if ( $positive_num > 1 && $negative_num > 1 ){
				my $denominator = $positive_variance_sum / ( $positive_num - 1 );
				$denominator   += $negative_variance_sum / ( $negative_num - 1 );
				if ( $denominator ){
					$f_score = $numerator / $denominator;
				}
			}
			
			unless(exists $hash_relevance{$cols[$col]}){
				$hash_relevance{$cols[$col]} = $f_score;
			} else {
				print STDERR "Error:$cols[$col] repeat\n";
				die;
			}
			
		}
	}
	
	
	# ******** Output relevance of features ********
	
	my $dest_dir = dirname($infile_feature);
	my $outfile_relevance = "$dest_dir/feature_relevance.txt";
	open SELECTION,">$outfile_relevance" or die;
	
	for my $feature ( sort {$hash_relevance{$b}<=>$hash_relevance{$a}} 
		keys %hash_relevance ){
		print SELECTION "$feature\t$hash_relevance{$feature}\n";
	}
	
	close SELECTION;
	
		
	# ******** Sort the features by decreasing order of relevance ********
	
	my @sorted_features = sort {$hash_relevance{$b}<=>
	$hash_relevance{$a}} keys %hash_relevance;
	
	my $all_feature_num = (scalar @sorted_features) - 1;
	if ( $keep_feature_num > $all_feature_num ){
		$keep_feature_num = $all_feature_num;
	}
	
	my $prev = 0;
	my $Next = 1;

	my %hash_featKeeps = ();
	$hash_featKeeps{ $sorted_features[$prev] } = 0;
	
	for my $f ( 1 .. $keep_feature_num ){
		
		# ============= Compute similarity ================
		
		my $prev_col = $sorted_features[$prev];
		my $current_col = $sorted_features[$f];
		
		my @prev_features_positive 		= ();
		my @prev_features_negative 	  	= ();
		my @current_features_positive 	= ();
		my @current_features_negative 	= ();
		
		for my $row ( @rows ){		# samples
		
			my $label = $hash_row{$row}{ $cols[$type_col_index] };
			
			if ( $label eq "positive" ){	
				
				if (exists $hash_row{$row}{$prev_col} ){
					push(@prev_features_positive,$hash_row{$row}{$prev_col});
				} else {
					print STDERR "Error:($row,$prev_col) not defined\n";
					die;
				}
			
				if ( exists $hash_row{$row}{$current_col}){
					push(@current_features_positive,$hash_row{$row}{$current_col});
				} else {
					print STDERR "Error:($row,$current_col) not defined\n";
					die;
				}
				
			} elsif ( $label eq "negative" ){
				
				if (exists $hash_row{$row}{$prev_col} ){
					push(@prev_features_negative,$hash_row{$row}{$prev_col});
				} else {
					print STDERR "Error:($row,$prev_col) not defined\n";
					die;
				}
			
				if ( exists $hash_row{$row}{$current_col}){
					push(@current_features_negative,$hash_row{$row}{$current_col});
				} else {
					print STDERR "Error:($row,$current_col) not defined\n";
					die;
				}
				
			} else {
				
				print STDERR "Error:($row,$label) not defined\n";
				die;
				
			}
		}
		
		
		my $similarity_positive = TOOLKIT::pearson(\@prev_features_positive,\@current_features_positive);
		my $similarity_negative = TOOLKIT::pearson(\@prev_features_negative,\@current_features_negative);
		my $similarity = 0.5 * ( $similarity_positive + $similarity_negative );
		
			
		# ============= Remove redundant features ================
		
		if ( $similarity <= $max_similarity ){
			unless(exists $hash_featKeeps{$current_col}){
				$hash_featKeeps{$current_col} = $Next;
			} else {
				print STDERR "Error:$current_col repeat\n";
				die;
			}
			
			$prev = $f;
			$Next++;
			
			if ( $Next == $keep_feature_num ){
				last;
			}
		}
	}
		
		
	# ******** Output selected features ********
	
	my $title = "Name";
	for my $col ( sort {$hash_featKeeps{$a} <=> 
		$hash_featKeeps{$b}} keys %hash_featKeeps ){
		$title .= "\t$col";
	}

	
	open  SELECTION,">$outfile_selected_feature" or die;
	print SELECTION "$title\ttype\n";	
		
	for my $row ( @rows ){
		my $str = $row;
		for my $col ( sort {$hash_featKeeps{$a} <=> 
			$hash_featKeeps{$b}} keys %hash_featKeeps ){
			unless( exists $hash_row{$row}{$col}){
				print STDERR "Error:($row,$col) not defined\n";
				die;
			}
			$str .= "\t$hash_row{$row}{$col}";
		}
		my $type = $hash_type{$row};
		print SELECTION "$str\t$type\n";
	}
	
	close SELECTION;
	
	
	# ******** Output selected features ********
		
	open  SELECTION,">$outfile_keep_feature_name" or die;
	for my $col ( sort {$hash_featKeeps{$a} <=> 
		$hash_featKeeps{$b}} keys %hash_featKeeps ){
		my $id = $hash_featKeeps{$col} + 1;
		print SELECTION "$id\t$col\n";
	}
	close SELECTION;
	
	return ($all_feature_num, $Next);
	
}
	
	
#############################################################	
	
sub load_feature {	
	
	my $feature_model = shift or die;
	my $infile_feature = shift or die;
	my $outfile_selected_feature = shift or die;
	
	
	my @array_feature_names = ();
	open LOADIN,"<$feature_model" or die;
	while(<LOADIN>){
		chomp;
		my @tmp = split/\s+/;
		my $id = shift @tmp;
		my $feature = shift @tmp;
		push(@array_feature_names,$feature);
	}
	close LOADIN;
	
	
	my %hash_all_feature = ();
	TOOLKIT::toMatrix($infile_feature,\%hash_all_feature);
	
	
	my $title = "Name";
	for my $feature_name ( @array_feature_names ){
		$title .= "\t$feature_name";
	}
	
	open  SELECTION,">$outfile_selected_feature" or die;
	print SELECTION "$title\ttype\n";
	for my $row ( sort keys %hash_all_feature ){
		my $str = $row;
		for my $feature_name ( @array_feature_names ){
			my $value = $hash_all_feature{$row}{$feature_name};
			$str .= "\t$value";
		}
		my $type = $hash_all_feature{$row}{type};
		print SELECTION "$str\t$type\n";
	}
	close SELECTION;
	
}
	

#############################################################
1;
