package FEATURE;
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;


use TOOLKIT;
use STRUCT;
use KnowledgeBasedEnergy;




####################################################
#################   评估特征   #####################
####################################################

sub F_score {
	
	my $infile_feature = shift or die;
	my $hash_fscore = shift or die;
	

	unless(-e $infile_feature){
		print STDERR "Error:$infile_feature not exists\n";
		die;
	}
	
	
	my %hash_row = ();
	my @cols = TOOLKIT::toMatrix($infile_feature,\%hash_row);
	
	
	my @rows = ();
	for my $row ( sort {$a cmp $b} keys %hash_row ){
		if ( $row !~/^\bname\b/i ){
			push(@rows,$row);
		}
	}
	
	
	my $feature_num = 0;
	for my $col ( 1 .. $#cols - 1 ){
	
		my $whole_sum 	 = 0;
		my $positive_sum = 0;
		my $negative_sum = 0;
		my $positive_num = 0;
		my $negative_num = 0;
		$feature_num++;
		
		for my $row ( @rows ){
			
			my $v = $hash_row{$row}{ $cols[$col] };
			my $label = $hash_row{$row}{ $cols[$#cols] };
			
			$whole_sum += $v;
			if ( $label eq "positive" ){		
				$positive_sum += $v;
				$positive_num++;		
			} elsif ( $label eq "negative" ){		
				$negative_sum += $v;
				$negative_num++;	
			} else {
				print STDERR "Error:($row,$cols[$col],$v)\n";
				die;
			}
		}
		
		my $whole_av 	= $whole_sum / ($positive_num + $negative_num);
		my $positive_av = $positive_sum / $positive_num;
		my $negative_av = $negative_sum / $negative_num;
		
		my $numerator = ( $positive_av - $whole_av ) * ( $positive_av - $whole_av );
		$numerator   += ( $negative_av - $whole_av ) * ( $negative_av - $whole_av );
		
		my $positive_variance_sum = 0;
		my $negative_variance_sum = 0;

		for my $row ( @rows ){
			
			my $v = $hash_row{$row}{ $cols[$col] };
			my $label = $hash_row{$row}{ $cols[$#cols] };
			
			if ( $label eq "positive" ){			
				my $t = ( $v - $positive_av ) * ( $v - $positive_av );
				$positive_variance_sum += $t;	
			} elsif ( $label eq "negative" ){		
				my $t = ( $v - $negative_av ) * ( $v - $negative_av );
				$negative_variance_sum += $t;	
			} else {
				print STDERR "Error:($row,$cols[$col],$v)\n";
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
		
		$hash_fscore->{$cols[$col]} = $f_score;
		
	}
	
	return $feature_num;
	
}



####################################################
#################   提取特征   #####################
####################################################

sub feature_extraction {
	
	my $struct_file = shift or die;
	my $hash_parameter = shift or die;
	my $hash_argument = shift or die;
	local *FEATURE = shift or die;

	
	my $label;
	my $title_flag = 0;
	while (@_) {
		my $argument = shift @_;	
		if ($argument=~/label/i) {$label=shift @_}
		if ($argument=~/title/i) {$title_flag=shift @_}
	}
	
	
	# ************************************** #
	# 获得参数
	# ************************************** #
	

	# parameters for statistic potential features
	my $rbin = $hash_parameter->{rbin};
	my $nbin = $hash_parameter->{nbin};
	
	
	# parameter for feature extraction
	my $is_only_knowledge_based_energy_features = $hash_parameter->{is_only_knowledge_based_energy_features};


	# array & hash reference
	my $adjacents = $hash_argument->{adjacents};
	my $hash_potential = $hash_argument->{hash_potential};


	##########################################


	my %hash_struct = ();
	STRUCT::parse_struct($struct_file,\%hash_struct);


	my $num_featrure  = 0;	
	my %hash_check_id = ();
	for my $id ( sort keys %hash_struct ){
	
		# check ID
		my $new_id = $id;	
		$new_id =~ s/^\s*//g;
		$new_id =~ s/\s*$//g;
		if($new_id =~ /(\S+)/){
			$new_id = $1;
		}
		unless(exists $hash_check_id{$new_id}){
			$hash_check_id{$new_id} = 1;
		} else {
			print STDERR "Error:$new_id repeat\n";
			die;
		}
		
		
		# check sequence
		my $seq = $hash_struct{$id}{seq};
		my $mfe = $hash_struct{$id}{mfe};
		my $struct = $hash_struct{$id}{struct};
		
		$seq =~ tr/[acgtTun\.]/[ACGUUUNN]/;	
		if( not ($seq =~/^([A|C|G|U|N]+)$/) ){
			next;
		}
		
		
		# *************************************** #
		
		
		my $n = 0;
		my $str = $new_id;		
		my %hash_feature = ();
		

		my %hash_seq_tmp = ();
		my %hash_struct_tmp = ();
		$hash_seq_tmp{$id} = $seq;
		$hash_struct_tmp{$id}{seq} = $seq;
		$hash_struct_tmp{$id}{struct} = $struct;
		


		#############################################################
		#############   Novel features    #####################
		#############################################################
			
			
		for my $k ( @$adjacents ){	
			my $energy_score = KnowledgeBasedEnergy::energy_score($nbin,$rbin,
			\%hash_struct_tmp,$hash_potential,'k_mer'=>$k);
			
			$energy_score = TOOLKIT::round($energy_score,6);
			$str .= "\t$energy_score";
			$hash_feature{$n} = "$k-energy_score-3";
			$n++;
		}
		
		
		
		

		########################################
		unless($is_only_knowledge_based_energy_features){
		########################################
		
		
		
		
		# the distribution of unpaired base to paired base in secondary structure
		my $N_bin = 10;
		my @rbp2UnpB_distribution = ();
		STRUCT::unpaired2paired_distribution($struct,$N_bin,\@rbp2UnpB_distribution);
		my $i2 = 0;
		for my $r ( @rbp2UnpB_distribution ){
			$str .= "\t$r";
			$hash_feature{$n} = "rbp2UnpB_dis:$i2";
			$n++;
			$i2++;
		}
			


		# the size of biggest bulge in structure
		my $max_bulge_size = STRUCT::biggest_bulge_size($struct);
		$str .= "\t$max_bulge_size";
		$hash_feature{$n} = "max_bulge_size";
		$n++;


		# normalized stems = n_stems / L
		my $normalized_n_stems = STRUCT::normalized_n_stems($struct);
		$str .= "\t$normalized_n_stems";
		$hash_feature{$n} = "normalized_n_stems";
		$n++;
		

			
		# normalized loops  = n_loops / L
		my $normalized_n_loops = STRUCT::normalized_n_loops($struct);
		$str .= "\t$normalized_n_loops";
		$hash_feature{$n} = "normalized_n_loops";
		$n++;
		
		
	
		
		#############################################################
		#############   triplet-SVM features    #####################
		#############################################################
		

		# triplet elements	
		#my @triplets = ();
		#my @triplets_str = ();
		#STRUCT::triplet($seq,$struct,\@triplets,\@triplets_str);
		#my $i1 = 0;
		#for my $r ( @triplets ){			
		#	$str .= "\t$r";
		#	$hash_feature{$n} = "triplets:$triplets_str[$i1]";
		#	$n++;
		#	$i1++;
		#}
	
		
		
		
		#############################################################
		#################    miRD features    #######################
		#############################################################
		
		
		
		# ratio of paired base to unpaired base (rpb)
		my ($rpb,$rUnpB,$rbp2UnpB) = STRUCT::paired_unpaired_base($struct);
		$str .= "\t$rpb";
		$hash_feature{$n} = "rpb";
		$n++;
		
		
	
		#miRD features
		#my $n_bin = 10;
		#my %hash_rabp = ();
		#STRUCT::reverse_align_bp($seq,$struct,$n_bin,\%hash_rabp);		
		#for my $i ( sort {$a<=>$b} keys %hash_rabp ){
		#	for my $j ( sort  keys %{ $hash_rabp{ $i } } ){
		#		my $r = $hash_rabp{ $i }{ $j };	
		#		$str .= "\t$r";
		#		$hash_feature{$n} = "Zhang:$i$j";
		#		$n++;
		#	}
		#}	
		
		
		
		
		#############################################################
		#################   miPred features   #######################
		#############################################################

		
		# GC content
		# %(G+C) = (|C|+|G|)/L * 100
		my $GCc = STRUCT::GC_content($seq);
		$str .= "\t$GCc";
		$hash_feature{$n} = "GC_content";
		$n++;
		
		
		# dinucleotide frequencies 
		# %XY = |XY| * 100 / ( L - 1 )
		# where X and Y belong to (A C G U)
		my @dinucleotide = ();
		my @dinucleotide_str = ();
		STRUCT::dinucleotide_frequency($seq,\@dinucleotide,\@dinucleotide_str);	
		my $i3 = 0;
		for my $di ( @dinucleotide ){
			$str .= "\t$di";
			$hash_feature{$n} = "dinucleotide:$dinucleotide_str[$i3]";
			$n++;
			$i3++;
		}
		
		
		# dG = MFE / L
		my $normalized_mfe = STRUCT::normalized_mfe($struct,$mfe);
		$str .= "\t$normalized_mfe";
		$hash_feature{$n} = "normalized_mfe";
		$n++;
		
		
		# MFE Index 1 = dG / %(G+C) (where dG = MFE / L)
		my $mfe_index_1 = STRUCT::mfe_index_1($seq,$mfe);
		$str .= "\t$mfe_index_1";
		$hash_feature{$n} = "mfe_index_1";
		$n++;
		
		
		# MFE Index 2 = dG / n_stems
		my $mfe_index_2 = STRUCT::mfe_index_2($struct,$mfe);	
		$str .= "\t$mfe_index_2";
		$hash_feature{$n} = "mfe_index_2";
		$n++;
		
		
			
		# Normalized base-paring propensity (dP = tot_bases / L)
		my $dP = STRUCT::normalized_base_paring_propensity($struct);
		$str .= "\t$dP";
		$hash_feature{$n} = "dP";
		$n++;
		
		
		
		
		
		#############################################################
		################## microPred features #######################
		#############################################################

		
		
		
		# MFE Index 3 = dG / n_loops
		my $mfe_index_3 = STRUCT::mfe_index_3($struct,$mfe);
		$str .= "\t$mfe_index_3";
		$hash_feature{$n} = "mfe_index_3";
		$n++;
		
		
		# |A-U|/L, |G-C|/L and |G-U|/L
		my @nbpc = ();
		my @nbpc_stems = ();
		my @nbpc_str = (); 
		my @nbpc_stems_str = ();
		STRUCT::normalized_base_pair_count($seq,
		$struct,\@nbpc,\@nbpc_stems,\@nbpc_str,\@nbpc_stems_str);
		my $i4 = 0;
		for my $r ( @nbpc ){
			$str .= "\t$r";
			$hash_feature{$n} = $nbpc_str[$i4]."/L";
			$n++;
			$i4++;
		}
		
		
		# %(A-U)/n_stems, %(G-C)/n_stems and %(G-U)/n_stems
		my $i5 = 0;
		for my $r ( @nbpc_stems ){
			$str .= "\t$r";
			$hash_feature{$n} = "%(".$nbpc_stems_str[$i5].")/n_stems";
			$n++;
			$i5++;
		}
		
		
		# Average base pairs per stem
		# avg_bp_stem = tot_bases / n_stems
		my $avg_bp_stem = STRUCT::avg_bp_stem($struct);
		$str .= "\t$avg_bp_stem";
		$hash_feature{$n} = "avg_bp_stem";
		$n++;
		
		
		
		#############################################################
		#################### ZmirP features #########################
		#############################################################
		
		
		
		# the maximum of consecutive paired base
		my $max_consecutive_bp = STRUCT::biggest_consecutive_paired_base($struct);
		$str .= "\t$max_consecutive_bp";
		$hash_feature{$n} = "max_consecutive_bp";
		$n++;
		
		
		
		# normalized bulges = n_bulges / L
		my $normalized_n_bulges = STRUCT::normalized_n_bulges($struct);
		$str .= "\t$normalized_n_bulges";
		$hash_feature{$n} = "normalized_n_bulges";
		$n++;
		
	
		
		# Average unpaired bases per bulge
		# avg_unbp_bulge = tot_unpaired_bases / n_bulges
		# where tot_unpaired_bases is the total number of unpaired bases
		my $avg_unbp_bulge = STRUCT::avg_unbp_bulge($struct);
		$str .= "\t$avg_unbp_bulge";
		$hash_feature{$n} = "avg_unbp_bulge";
		$n++;
		
		
		
		# MFE Index 4 = dG / tot_bases
		# where tot_bases is the total number of base pairs
		my $mfe_index_4 = STRUCT::mfe_index_4($struct,$mfe);
		$str .= "\t$mfe_index_4";
		$hash_feature{$n} = "mfe_index_4";
		$n++;
		
		
		
		# MFE Index 5 = dG / n_bulges
		my $mfe_index_5 = STRUCT::mfe_index_5($struct,$mfe);
		$str .= "\t$mfe_index_5";
		$hash_feature{$n} = "mfe_index_5";
		$n++;
		
		
		
		
		########################################
		}  # End of  unless($is_only_knowledge_based_energy_features)
		########################################
		
		
		
		
		
		# ***************** OUTPUT *************** #
		
		if ( $title_flag ){
			my $title = "Name";
			my %hash_tmp = ();
			for my $f ( sort {$a <=> $b} keys %hash_feature ){
				my $feature = $hash_feature{$f};
				$title .= "\t$feature";
				unless(exists $hash_tmp{$feature}){
					$hash_tmp{$feature} = 1;
				} else {
					print STDERR "Error:($f,$feature) repeat\n";
					die;
				}
			}
			
			if(defined $label){
				print FEATURE "$title\ttype\n";
			}else{
				print FEATURE "$title\n";
			}
			
			$title_flag = 0;
		}
		
		
		# ***************************************** #
		
		if(defined $label){
			print FEATURE "$str\t$label\n";
		}else{
			print FEATURE "$str\n";
		}
		
		unless( $num_featrure ){
			$num_featrure = $n;
		} elsif ( $num_featrure != $n ){
			print STDERR "Error:lack feature ($id\n$seq\n)\n";
			die;
		}
		
	}
	

	return $num_featrure;
	
}


####################################################
1;


