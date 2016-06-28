package LIBSVM;
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Basename;


####################################################

sub change_label_type {
	
	my $infile_data = shift or die;
	my $feature_file = shift or die;
	
	
	open IN,"<$infile_data" or die;
	open OUT,">$feature_file" or die;
	
	my %hash_name = ();
	my $sample_size = 0;
	while(defined (my $line=<IN>)){
		chomp($line);
		
		$line =~ s/^\s*//g;
		$line =~ s/\s*$//g;
	
		if ($line =~ /^\s*$/) {
			next;
		}
		
		if ($line =~ /^Name\b/i){
			next;
		}
		
		my @tmp = split /\t/,$line;
		my $name = shift @tmp;
		my $type = pop @tmp;
		
		if ( $type =~ /positive/i ){
			$type = 1;
			$line = join("\t",@tmp);
			$line = "$type\t$line";
		} elsif ( $type =~ /negative/i ){
			$type = -1;
			$line = join("\t",@tmp);
			$line = "$type\t$line";
		} else {
			push(@tmp,$type);
			$line = join("\t",@tmp);
		}
			
		print OUT "$line\n";
		unless(exists $hash_name{$name}){
			$hash_name{$name} = 1;
		} else {
			print STDERR "Error:$name repeat in $infile_data\n";
			die;
		}
		$sample_size++;
			
	}
	
	close IN;
	close OUT;
	
}


####################################################

sub toLibsvm {
	
	my $usage = "$0 < infile > < outfile >";
	my $infile = shift or die $usage;
	my $outfile = shift or die $usage;
	
	
	open IN,"<$infile" or die;
	open OUT,">$outfile" or die;  
	while(<IN>){   
		chomp;   
		my @cont=split;   
		my $rlable = $cont[0]; 
		my $str = $rlable;   
		for(my $i=1;$i<@cont;$i++){   
		    if($cont[$i]!=0){   
		        $str .= "\t$i:$cont[$i]";   
		    }   
		}   
		print OUT "$str\n";   
	} 
	close IN;
	close OUT;
	
}



####################################################
################      TRAIN      ###################
####################################################


sub scale_training_data {
	
	my $libsvm_path 		= shift or die;
	my $infile				= shift or die;
	my $outfile_scaled_data = shift or die;	
	my $outfile_scaled_model= shift or die;


	my $current_dir  = getcwd;
	chdir( $libsvm_path );
	
	
	unless(-e "svm-scale"){
		print STDERR "Error:svm-scale not exists in $libsvm_path\n\n\n";
		die;
	}


	my $os = $^O;
	if ( $os =~ /win/i ){
		die if system("svm-scale -s $outfile_scaled_model $infile > $outfile_scaled_data");
	} else {
		die if system("./svm-scale -s $outfile_scaled_model $infile > $outfile_scaled_data");
	}
	
	chdir($current_dir);
	
}


####################################################

sub search_parameter {	
	
	my $libsvm_path 	  = shift or die;
	my $infile 			  = shift or die;
	my $outfile_parameter = shift or die;
	my $cv 			  	  = shift or die;
	my $kernel			  = shift;
	
	
	unless(defined $kernel){
		print STDERR "Error:kernel function not defined\n";
		die;
	}
	
	
	my $os = $^O;
	unless ( $os =~ /win/i ){
		$libsvm_path .= "tools/";
	}
	
	
	my $current_dir = getcwd;
	chdir( $libsvm_path );


	unless(-e "grid.py"){
		print STDERR "Error:grid.py not exists in $libsvm_path\n\n\n";
		die;
	}
	
	
	die if system("python grid.py -v $cv -t $kernel $infile > $outfile_parameter");
	
	
	chdir( $current_dir );
	
	#get the last blank line of file
	my ($c,$g) = ();
	if(-e $outfile_parameter){
		
		open IN,"<$outfile_parameter" or die;	
		my $last_line = "";
		while(<IN>){
			next if /^\s*$/;
			$last_line = $_;
		}
		close IN;	
		
		my @tmp = split /\s+/,$last_line;
		$c = $tmp[0];
		$g = $tmp[1];
		
	} else {
		print STDERR "Error:$outfile_parameter not exists\n";
		die;
	}

	return ($c,$g);
	
}


####################################################


sub do_train {
	
	my $libsvm_path 	= shift or die;
	my $infile			= shift or die;
	my $c 				= shift or die;
	my $g			 	= shift or die;
	my $kernel 			= shift or die;
	my $outfile_model 	= shift or die;


	my $current_dir = getcwd;
	chdir( $libsvm_path );
	
	unless(-e "svm-train"){
		print STDERR "Error:svm-train not exists in $libsvm_path\n";
		die;
	}

	my $os = $^O;
	if ( $os =~ /win/i ){
		die if system("svm-train -c $c -g $g -t $kernel -b 1 $infile $outfile_model");
	} else {	
		die if system("./svm-train -c $c -g $g -t $kernel -b 1 $infile $outfile_model");
	}
	
	chdir( $current_dir );
	
}


####################################################


sub training {
	
	my $libsvm_path = shift or die;
	my $infile_data = shift or die;
	my $svm_model 	= shift or die;
	my $scaled_model= shift or die;
	my $dest_dir	= shift or die;
	
	
	# ****************************************************** #
	
	my $cv = 5;
	my $kernel = 2;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/cv/i) {$cv=shift @_}
		if ($argument=~/kernel/i) {$kernel=shift @_}
	}
	
	print "\nLibsvm::training\n";
	print "cv for grid search:$cv\n";
	print "kernel:$kernel\n\n";
	
	
	# ****************************************************** #


	$libsvm_path =~ s/\/$//g;
	$libsvm_path .= "/";
	
	
	# ****************************************************** #
	
	
	my $feature_file = "$dest_dir/training.feature.txt";
	my $sample_size = &change_label_type($infile_data,$feature_file);
	
	
	
	# ****************************************************** #
	
	my $libsvm_file = "$dest_dir/libsvm.txt";
	&toLibsvm($feature_file,$libsvm_file);
	
	
	# ****************************************************** #
	
	my $scaled_file = "$dest_dir/scaled.txt";
	&scale_training_data($libsvm_path,$libsvm_file,$scaled_file,$scaled_model);


	# ****************************************************** #
	
	
	$cv = $sample_size if $cv =~ /loo/i;
	my $parameter_file = "$dest_dir/parameter.txt";
	my ($c,$g) = &search_parameter($libsvm_path,$scaled_file,$parameter_file,$cv,$kernel);
		
		
	# ****************************************************** #
	
	&do_train($libsvm_path,$scaled_file,$c,$g,$kernel,$svm_model);
	
	
	# ****************************************************** #
	
}
	





####################################################
#################   PREDICTION   ###################
####################################################



sub scale_test_data {

	my $libsvm_path 		= shift or die;
	my $infile_scaled_model	= shift or die;
	my $infile_test_data	= shift or die;
	my $outfile_test_data 	= shift or die;
	
	
	
	my $current_dir = getcwd;
	chdir( $libsvm_path );
	
	
	my $os = $^O;
	if ( $os =~ /win/i ){
		die if system("svm-scale -r $infile_scaled_model $infile_test_data > $outfile_test_data");
	} else {
		die if system("./svm-scale -r $infile_scaled_model $infile_test_data > $outfile_test_data");
	}
	
	chdir($current_dir);
	
}



####################################################

sub parse_result {
	
	my $infile_feature = shift or die;
	my $infile_original_result = shift or die;
	my $outfile_parsed_result = shift or die;
	
	
	#get Name of line
	my %hash_name = ();
	my $Line_number = 1;
	
	open IN,"<$infile_feature" or die;
	while(defined(my $line = <IN> )){
		chomp($line);
		
		$line =~ s/^\s*//g;
		$line =~ s/\s*$//g;
	
		if ($line =~ /^\s*$/) {
			next;
		}
		
		if ($line =~ /^Name\b/i){
			next;
		}
		
		my @tmp = split /\t/,$line;
		my $name = shift( @tmp );
		$hash_name{ $Line_number } = $name;
		$Line_number++;
	}
	close IN;
	
	
	#output paresed results
	$Line_number = 0;
	my $col_positive;
	my $col_negative;
	
	open IN,"<$infile_original_result" or die;
	open OUT,">$outfile_parsed_result" or die;
	
	while(defined(my $line=<IN>)){
		chomp($line);
		next if $line =~/^\s*$/;
	
		my @tmp = split /\s+/,$line;
		
		unless ( $Line_number ){	#the first line
	
			if ( $tmp[2] != 1 && $tmp[2] != -1 ){
				print STDERR "Error:label error\n";
				die;
			}	
			
			if ( $tmp[1] == 1 ){
				$col_positive = 1;
				$col_negative = 2;
			} elsif ( $tmp[2] == 1 ){
				$col_positive = 2;
				$col_negative = 1;
			}
			
			print OUT "Name\tpredictType\tnegative\tpositive\n";
			
		} else {
			
			my $name = $hash_name{ $Line_number };
			my $positive_prob = $tmp[$col_positive];
			my $negative_prob = $tmp[$col_negative];
			my $predictType = 1;
			if( $positive_prob < $negative_prob ){
				$predictType = -1;
			}
			
			print OUT "$name\t$predictType\t$negative_prob\t$positive_prob\n";
		}
		
		$Line_number++;
	}	
	close IN;
	close OUT;
	
}



####################################################


sub do_predict {
	
	my $libsvm_path 		= shift or die;
	my $infile_scaled_test 	= shift or die;
	my $infile_model 		= shift or die;
	my $outfile_test	 	= shift or die;


	my $current_dir = getcwd;
	chdir( $libsvm_path );

	my $os = $^O;
	if ( $os =~ /win/i ){
		die if system("svm-predict -b 1 $infile_scaled_test $infile_model $outfile_test");	
	} else {
		die if system("./svm-predict -b 1 $infile_scaled_test $infile_model $outfile_test");	
	}
	
	chdir( $current_dir ); 
	
}



####################################################


sub predict {
	
	my $libsvm_path = shift or die;
	my $infile_data = shift or die;
	my $svm_model	= shift or die;
	my $scaled_model= shift or die;
	my $outfile_out = shift or die;
	my $dest_dir 	= shift or die;
	
	
	# ****************************************************** #
	
	
	$libsvm_path  =~ s/\/$//g;
	$libsvm_path .= "/";
	

	# ****************************************************** #
	
	
	my $feature_file = "$dest_dir/predict.feature.txt";
	my $sample_size = &change_label_type($infile_data,$feature_file);
	

	# ****************************************************** #
	
	my $libsvm_file = "$dest_dir/libsvm.txt";
	&toLibsvm($feature_file,$libsvm_file);
	
	
	# ****************************************************** #
	
	my $scaled_file = "$dest_dir/scaled.txt";
	&scale_test_data($libsvm_path,$scaled_model,$libsvm_file,$scaled_file);
		
	
	# ****************************************************** #
	
	my $result_file = "$dest_dir/original_result.txt";
	&do_predict($libsvm_path,$scaled_file,$svm_model,$result_file);


	# ****************************************************** #

	&parse_result($infile_data,$result_file,$outfile_out);

}



####################################################
1;