package CROSS_VALIDATION;
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;


use TOOLKIT;
use FEATURE;
use LIBSVM;
use KnowledgeBasedEnergy;




####################################################
#################   交叉检验   #####################
####################################################	

sub cross_validation {
	
	my $positive_fasta 	= shift or die;
	my $negative_fasta 	= shift or die;
	my $outfile_roc 	= shift or die;
	my $outfile_out 	= shift or die;
	my $dest_dir		= shift or die;
	my $hash_parameter	= shift or die;
	my $hash_argument	= shift or die;
	my $libsvm_path		= shift or die;
	
	
	
	
	# ************************************** #
	# 获得参数
	# ************************************** #
	
	
	# parameters for CV
	my $cv = $hash_parameter->{cv};
	my $round_num = $hash_parameter->{round_num};

	# parameters for statistic potential features
	my $rbin = $hash_parameter->{rbin};
	my $nbin = $hash_parameter->{nbin};
	
	# array reference
	my $adjacents = $hash_argument->{adjacents};
	
	
	

	# ****************************************************** #
	# randomly partitioned into k equal size subsamples
	# ****************************************************** #
	
	
	my %hash_pos_seq = ();
	my %hash_neg_seq = ();
	TOOLKIT::readfasta($positive_fasta,\%hash_pos_seq);
	TOOLKIT::readfasta($negative_fasta,\%hash_neg_seq);
	
	my @seqIDs = ();
	for my $id ( sort keys %hash_pos_seq){
		push(@seqIDs,$id);
	}
	for my $id ( sort keys %hash_neg_seq){
		push(@seqIDs,$id);
	}
	
	
	my $sample_size = scalar @seqIDs;	
	$cv = $sample_size if $cv =~ /loo/i;	
	my $subsample_size = int($sample_size/$cv);	
	$subsample_size = 1 if $subsample_size < 1;
	
	
	# ****************************************************** #
	# To reduce variability, multiple rounds of cross-validation 
	# are performed using different partitions, 
	# and the validation results are averaged over the rounds.
	# ****************************************************** #
	
	
	my $round_n = 0;
	my %hash_parse = ();
	for my $round ( 1 .. $round_num ){
		
		$round_n++;
				
		my $round_dir = "$dest_dir/round$round";
		mkdir($round_dir);
		
		my $log_file = "$round_dir/log.txt";
		open CVLOG,">$log_file" or die;
		
		print "\n######## round:$round ########\n";
		print CVLOG "######## round:$round ########\n";
		print CVLOG "cv:$cv\n";

		
		# ****************************************************** #
		# shuffle samples
		# ****************************************************** #
		
	
		
		TOOLKIT::shuffle(\@seqIDs);
		
		
		
		# ****************************************************** #
		# partitize samples
		# ****************************************************** #
		
		my $k = 1;
		my %hash_subsample = ();	
		for my $i ( 0 .. $#seqIDs ){
			$hash_subsample{$k}{ $seqIDs[$i] } = 1;		
			if ( ( ($i+1) % $subsample_size == 0 ) && ( $k < $cv ) ){
				$k++;
			}
		}
		
		print CVLOG "\nSize of subsamples:\n";
		
		my $total_sample = 0;
		for my $k ( sort {$a<=>$b} keys %hash_subsample ){
			my $t = scalar keys %{$hash_subsample{$k}};
			$total_sample += $t;
			print "$k:$t\n";
			print CVLOG "$k:$t\n";
		}
		
		if ( $total_sample != $sample_size ){
			print STDERR "Error:$sample_size\t$total_sample\n";
			die;
		}
		
		
		# ****************************************************** #
		# a single subsample is retained as the validtion data 
		# for testing the model, and the remaining k - 1 subsamples
		# are used as training data.
		# The cross-validation process is then repeated k times (the folds),
		# with each of the k subsamples used exactly once as the validation data.
		# The k results from the folds then can be averaged 
		# (or otherwise combined) to produce a single estimation.
		# ****************************************************** #
		
			
		my %hash_out = ();
		for my $k ( sort {$a<=>$b} keys %hash_subsample ){
			
			
			####################################################
			#################   构造数据集  ####################
			####################################################
			
			
			print "\n********* $k *********\n";
			
			my $dir = "$round_dir/$k";
			mkdir($dir);
			
			my $test_dir = "$dir/test";
			mkdir($test_dir);
			
			my $train_dir = "$dir/train";
			mkdir($train_dir);
			
			
			# ============ Test dataset =============
			
			my $test_fasta = "$test_dir/test.fa";
			open TEST,">$test_fasta" or die;
			for my $id ( keys %{$hash_subsample{$k}} ){
				my $seq;
				if(exists $hash_pos_seq{$id}){
					$seq = $hash_pos_seq{$id};
				}else{
					$seq = $hash_neg_seq{$id};
				}
				print TEST ">$id\n$seq\n";
			}
			close TEST;
			
		
			# ============ Traning dataset =============
			
			my $training_positive_fasta = "$train_dir/positive.fa";
			my $training_negative_fasta = "$train_dir/negative.fa";
			open POS,">$training_positive_fasta" or die;	
			open NEG,">$training_negative_fasta" or die;
			for my $n ( sort {$a<=>$b} keys %hash_subsample ){
				next if $n == $k;
				for my $id ( keys %{$hash_subsample{$n}} ){
					if(exists $hash_pos_seq{$id}){
						print POS ">$id\n$hash_pos_seq{$id}\n";
					}else{
						print NEG ">$id\n$hash_neg_seq{$id}\n";
					}
				}
			}
			close POS;
			close NEG;
			
			
			
			####################################################
			#################   提取特征   #####################
			####################################################
			
			
			
			# ========== Predict secondary structure =============
			
			
			my $training_positive_struct = "$training_positive_fasta.struct";
			STRUCT::fold_linux_RNA($training_positive_fasta,$training_positive_struct);

			my $training_negative_struct = "$training_negative_fasta.struct";
			STRUCT::fold_linux_RNA($training_negative_fasta,$training_negative_struct);

			my $test_struct = "$test_fasta.struct";
			STRUCT::fold_linux_RNA($test_fasta,$test_struct);
			
			
			# ========== Calculate statistical potentials =============
			
			
			my %hash_potential = ();
			
			for my $k ( @$adjacents ){
				KnowledgeBasedEnergy::statistical_potentials($nbin,$rbin,
				$training_positive_struct,$training_negative_struct,
				\%hash_potential,'k_mer'=>$k);
			}
			
	
		
			
			# ==========  Feature extraction ==========
			
			
			# array & hash reference
			my %hash_argument_tmp = ();
			$hash_argument_tmp{adjacents} = $adjacents;
			$hash_argument_tmp{hash_potential} = \%hash_potential;


			my $train_feature = "$train_dir/feature.txt";
			open FEATURETRAIN,">$train_feature" or die;
			
			FEATURE::feature_extraction($training_positive_struct,
			$hash_parameter,\%hash_argument_tmp,\*FEATURETRAIN,
			'title'=>1,'label'=>"positive");

			FEATURE::feature_extraction($training_negative_struct,
			$hash_parameter,\%hash_argument_tmp,\*FEATURETRAIN,
			'label'=>"negative");
			
			close FEATURETRAIN;
			
			
			
		
		
			####################################################
			#################   训练和测试   ###################
			####################################################
		
		
			# ================ Traning =================
				
			
			my $training_svm_model = "$train_dir/model.txt";
			my $training_scaled_model = "$train_dir/scaled_model.txt";
			LIBSVM::training($libsvm_path,$train_feature,$training_svm_model,
			$training_scaled_model,$train_dir);
			
			
	
			# ================ Test =================
			
			
			my $test_feature = "$test_dir/feature.txt";
			open FEATURETEST,">$test_feature" or die;
			
			FEATURE::feature_extraction($test_struct,$hash_parameter,
			\%hash_argument_tmp,\*FEATURETEST,'title'=>1,'label'=>"positive");
			
			close FEATURETEST;
			
			
			my $outfile_out = "$test_dir/$k-out.txt";
			LIBSVM::predict($libsvm_path,$test_feature,$training_svm_model,
			$training_scaled_model,$outfile_out,$test_dir);
		
			
			
			# ==================== out ==================
			
			my %hash_out_tmp = ();
			TOOLKIT::toMatrix($outfile_out,\%hash_out_tmp);
			
			for my $id ( sort keys %hash_out_tmp ){
				for my $col ( keys %{$hash_out_tmp{$id}} ){
					unless(exists $hash_out{$id}{$col}){
						$hash_out{$id}{$col} = $hash_out_tmp{$id}{$col};
					} else {
						print STDERR "Error:($id,$col) repeat\n";
						die;
					}
				}
			}
			
		###################
		} # End k-fold
		###################
				
				
		my $outfile_out = "$round_dir/round$round-out.txt";
		open  CVOUT,">$outfile_out" or die;
		print CVOUT "Name\tpredictLable\tnegative\tpositive\n";
		
		for my $id ( sort keys %hash_out ){
			my $predictLable = $hash_out{$id}{predictType};
			my $proP = $hash_out{$id}{positive};
			my $proN = $hash_out{$id}{negative};
			print CVOUT "$id\t$predictLable\t$proN\t$proP\n";
			
			unless(exists $hash_parse{$id}{positive}){
				$hash_parse{$id}{positive} = $proP;
			} else {
				$hash_parse{$id}{positive}+= $proP;
			}
			unless(exists $hash_parse{$id}{negative}){
				$hash_parse{$id}{negative} = $proN;
			} else {
				$hash_parse{$id}{negative}+= $proN;
			}
		}
		
		close CVOUT;
		close CVLOG;
		
		
		######################
		last if $subsample_size == 1; # LOO
		######################
		
	
	#####################	
	} # End rounds
	#####################
	
	
	# ****************************************************** #
	# the validation results are averaged over the rounds.
	# ****************************************************** #
	
	my %hash_predict = ();
	my $positive_score_threshold = 0.5;
	
	open  CVOUT,">$outfile_out" or die;
	print CVOUT "Name\tpredictType\tnegative\tpositive\n";
	
	for my $id ( sort keys %hash_parse ){		
		
		my $av_proP = &round($hash_parse{$id}{positive}/$round_n,6);
		my $av_proN = &round($hash_parse{$id}{negative}/$round_n,6);
		my $predictLable = -1;
		if ( $av_proP >= $positive_score_threshold ){
			$predictLable = 1;
		}
		print CVOUT "$id\t$predictLable\t$av_proN\t$av_proP\n";
	
		$hash_predict{$id}{proP} = $av_proP;
			
		if(exists $hash_pos_seq{$id}){
			$hash_predict{$id}{real} = 1;
		}elsif ( exists $hash_neg_seq{$id} ){
			$hash_predict{$id}{real} = -1;
		}else {
			print STDERR "Error: label must be 1 or -1 ($id)\n";
			die;
		}	
	}
	close CVOUT;
	
	
	
	
	
	# ================= draw ROC curve ===================
	
	my $cutoff_step = 0.001; 
	&draw_roc_curve(\%hash_predict,$outfile_roc,$cutoff_step);
	
}


####################################################
################## 计算性能  ######################
####################################################

sub draw_roc_curve {
	
	my $hash_predict 	= shift or die;
	my $outfile_roc_data= shift or die;
	my $cutoff_step 	= shift or die;

	unless(defined $cutoff_step){
		print STDERR "Error:cut-off step not defined\n";
		die;
	}
	
	open  OUT,">$outfile_roc_data" or die;
	print OUT "Cutoff\tTOT_P\tTOT_N\tTP\tTN\tFP\tFN\tPr\tAc\tGm\t1-Sp\tSn\tSp\tMcc\n";
	
	my $cutoff = 0;
	while ( $cutoff < 1 ) {	
		
		my %hash_count = ();
		&TP_TN_FP_FN($hash_predict,\%hash_count,$cutoff);
		
		my %hash_perf = ();
		&sn_sp_ac_mcc(\%hash_count,\%hash_perf);
	
		my $TP = $hash_count{TP};
		my $TN = $hash_count{TN};
		my $FP = $hash_count{FP};
		my $FN = $hash_count{FN};
		my $TOT_P = $hash_count{total_pos};
		my $TOT_N = $hash_count{total_neg};
		
		my $sn = $hash_perf{SN};
		my $sp = $hash_perf{SP};
		my $ac = $hash_perf{AC};
		my $pr = $hash_perf{PR};
		my $Gm = $hash_perf{GM};
		my $mcc= $hash_perf{MCC};
			
		my $t = "";
		if( &is_numeric($sp) ){
			$t = &round(1-$sp,6);
		} else {
			$t = "NaN";
		}
		
		$cutoff = &round($cutoff,4);
		print OUT "$cutoff\t$TOT_P\t$TOT_N\t$TP\t$TN\t$FP\t$FN\t$pr\t$ac\t$Gm\t$t\t$sn\t$sp\t$mcc\n";
		
		$cutoff += $cutoff_step;

	}
	
	close OUT;
	
}



####################################################

sub TP_TN_FP_FN {
	
	my $hash_predict= shift or die;
	my $hash_count 	= shift or die;
	my $cutoff		= shift;
	
	my $total_pos = 0;	
	my $total_neg = 0; 	

	foreach my $id ( sort keys %$hash_predict ){		
		if( $hash_predict->{$id}{real} == 1 ){ 		
			$total_pos++;		
			if ( $hash_predict->{$id}{proP} >= $cutoff ){ 
				unless ( exists $hash_count->{TP} ){
					$hash_count->{TP} = 1;
				} else {
					$hash_count->{TP}++;
				}		
			} else { 
				unless ( exists $hash_count->{FN} ){
					$hash_count->{FN} = 1;
				} else {
					$hash_count->{FN}++;
				}
			}				
		} elsif ( $hash_predict->{$id}{real} == -1 ){
			$total_neg++;
			if ( $hash_predict->{$id}{proP} >= $cutoff ){ 
				unless ( exists $hash_count->{FP} ){
					$hash_count->{FP} = 1;
				} else {
					$hash_count->{FP}++;
				}	
			} else { 
				unless ( exists $hash_count->{TN} ){
					$hash_count->{TN} = 1;
				} else {
					$hash_count->{TN}++;
				}
			}
		} else {
			print STDERR "Error: label must be 1 or -1 ($hash_predict->{$id}{real})\n";
			die;
		}				
	}
	
	$hash_count->{total_pos} = $total_pos;
	$hash_count->{total_neg} = $total_neg;
	

	unless ( exists $hash_count->{TP} ){
		$hash_count->{TP} = 0;
	}
	unless ( exists $hash_count->{TN} ){
		$hash_count->{TN} = 0;
	}
	unless ( exists $hash_count->{FP} ){
		$hash_count->{FP} = 0;
	}
	unless ( exists $hash_count->{FN} ){
		$hash_count->{FN} = 0;
	}
	
}



####################################################

sub sn_sp_ac_mcc {
	
	my $hash_count = shift or die;
	my $hash_perf  = shift or die;

	my $TP = $hash_count->{TP};
	my $TN = $hash_count->{TN};
	my $FP = $hash_count->{FP};
	my $FN = $hash_count->{FN};
	
	my $sn = "NaN";
	my $sp = "NaN";
	my $ac = "NaN";
	my $pr = "NaN";
	my $Gm = "NaN";
	my $mcc = "NaN";
	
	if ( $TP + $FN > 0 ) {
		$sn = $TP / ( $TP + $FN );
		$sn = &round($sn,6);
	}

	if ( $TN + $FP > 0  ){
		$sp = $TN / ( $TN + $FP );
		$sp = &round($sp,6);
	} 

	if ( $TP + $FP + $TN + $FN > 0 ) {
		$ac = ($TP + $TN) / ( $TP + $FP + $TN + $FN );
		$ac = &round($ac,6);
	}

 	if ( $TP + $FP > 0 ){
 		$pr = $TP / ( $TP + $FP );
 		$pr = &round($pr,6);
 	}
 	
	if( ($TP+$FN) * ($TN+$FP) * ($TP+$FP) * ($TN+$FN) > 0 ){
		$mcc = ( ($TP*$TN) - ($FN*$FP) ) / sqrt( ($TP+$FN) * ($TN+$FP) * ($TP+$FP) * ($TN+$FN) );
		$mcc = &round($mcc, 6);
	}
	
	if ( (&is_numeric($sn)) && (&is_numeric($sp)) && $sn * $sp >= 0 ){
		$Gm = sqrt( $sn * $sp );
		$Gm = &round($Gm,6);
	}
			
	
	$hash_perf->{SN}  = $sn;
	$hash_perf->{SP}  = $sp;
	$hash_perf->{AC}  = $ac;
	$hash_perf->{PR}  = $pr;
	$hash_perf->{GM}  = $Gm;
	$hash_perf->{MCC} = $mcc;
	
	
	my $t = "";
	if( &is_numeric($sp) ){
		$t = &round(1-$sp,6);
	} else {
		$t = "NaN";
	}
	
	
	return "TP:$TP\tTN:$TN\tFP:$FP\tFN:$FN\tPR:$pr\tAC:$ac\tGM:$Gm\t1-SP:$t\tSN:$sn\tSP:$sp\tMCC:$mcc";
	
}


####################################################

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

sub round {
	
	my $usage = "<  numerical value > < the numbers after the decimal point >";
    my $val = shift;
    my $col = shift;
    
    unless(defined $val){
    	die $usage;
    }
    
    unless(defined $col){
    	die $usage;
    }
    
    unless( &is_numeric($val) ){
    	print STDERR "Error:$val not a numeric";
		die;
	}
	
	unless( &is_numeric($col) ){
    	print STDERR "Error:$col not a numeric";
		die;
	}
    
    my $r = 10 ** $col;
    my $a = ($val > 0) ? 0.5 : -0.5;
    
    return int($val * $r + $a) / $r;
    
}

####################################################
1;
