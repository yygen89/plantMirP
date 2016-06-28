#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;
use Getopt::Long  qw(:config bundling);



use TOOLKIT;
use STRUCT;
use FEATURE;
use LIBSVM;



	
# parameters
my $dir;
my $help;	
my $test_fasta;
my $positive_fasta;
my $negative_fasta;


# model files
my $svm_model;
my $scaled_model;
my $potential_model;


#set this parameter to the path to the Libsvm binary
my $libsvm_path = "/home/yyg/libsvm/";


# parameter for feature extraction 
my $is_only_knowledge_based_energy_features = 0;


# parameters for statistic potential features
my $nbin = 20;
my $rbin = 1;
my $adjacent = 5;


GetOptions( 	
	"help|h!"			=> \$help,
	"dir|D=s"			=> \$dir,
	"testFasta|T=s"		=> \$test_fasta,
	"potential|P=s"		=> \$potential_model,
	"svmModel|S=s"	=> \$svm_model,
	"scaleModel|U=s"	=> \$scaled_model,
	"posFasta|K=s"		=> \$positive_fasta,
	"negFasta|G=s"		=> \$negative_fasta,
	"libsvmPath|L=s"	=> \$libsvm_path,
	"isEnergyF|J=i"		=> \$is_only_knowledge_based_energy_features,	
);






########################################
# start 
########################################


if ( $help ) {
	&usage(); exit;
}


if( defined $test_fasta ){
	$test_fasta = TOOLKIT::getAbsPath($test_fasta);
}else{
	&usage(); exit;
}


if( defined $potential_model ){
	$potential_model = TOOLKIT::getAbsPath($potential_model);
}else{
	&usage(); exit;
}


if( defined $svm_model ){
	$svm_model = TOOLKIT::getAbsPath($svm_model);
}else{
	&usage(); exit;
}


if( defined $scaled_model ){
	$scaled_model = TOOLKIT::getAbsPath($scaled_model);
}else{
	&usage(); exit;
}


my $current_dir = getcwd;
unless( defined $dir ){
	$dir = "$current_dir/predict";
} else {
	$dir = TOOLKIT::getAbsPath($dir);
}
mkdir($dir);




my $start_time = TOOLKIT::start();



#############################################################
# Load or calculate statistical potentials 
#############################################################


my %hash_potential = ();
my %hash_parameter = ();
my %hash_argument = ();
if(defined $potential_model ){
	
	##########################################################
	# Load statistical potentials 
	##########################################################

	&load_potential($potential_model,\%hash_potential,\%hash_parameter);


	# parameter for feature extraction
	$hash_parameter{is_only_knowledge_based_energy_features} = $is_only_knowledge_based_energy_features;


	# array & hash reference
	my @adjacents = ();
	my $adjacent = $hash_parameter{adjacent};
	for my $r ( 1 .. $adjacent ){
		push(@adjacents,$r);
	}
	$hash_argument{adjacents} = \@adjacents;
	$hash_argument{hash_potential} = \%hash_potential;


}elsif ( (defined $positive_fasta) and (defined $negative_fasta) ){
	
	
	#############################################################
	# Predict secondary structure for training dataset
	#############################################################


	my $positive_struct = "$dir/positive_fasta_struct";
	STRUCT::fold_linux_RNA($positive_fasta,$positive_struct);

	my $negative_struct = "$dir/negative_fasta_struct";
	STRUCT::fold_linux_RNA($negative_fasta,$negative_struct);

	
	
	#########################################################
	# Calculate statistical potentials 
	#########################################################


	# parameters for statistic potential features
	$hash_parameter{rbin} = $rbin;
	$hash_parameter{nbin} = $nbin;
	$hash_parameter{adjacent} = $adjacent;


	# parameter for feature extraction
	$hash_parameter{is_only_knowledge_based_energy_features} = $is_only_knowledge_based_energy_features;


	# array & hash reference
	my @adjacents = ( );
	for my $r ( 1 .. $adjacent ){
		push(@adjacents,$r);
	}
	my %hash_argument = ();
	$hash_argument{adjacents} = \@adjacents;


	for my $k ( @adjacents ){
		KnowledgeBasedEnergy::statistical_potentials($nbin,$rbin,
		$positive_struct,$negative_struct,
		\%hash_potential,'k_mer'=>$k);
	}


	my $potential_model = "$dir/potential.txt";
	open POTENTIAL,">$potential_model" or die;
	print POTENTIAL "#nbin\t$nbin\n";
	print POTENTIAL "#rbin\t$rbin\n";
	print POTENTIAL "#adjacent\t$adjacent\n";

	for my $base_ij ( sort keys %hash_potential ){
	for my $nbin ( sort {$a<=>$b} keys %{$hash_potential{$base_ij}} ){
		my $value = $hash_potential{$base_ij}{$nbin};
		print POTENTIAL "$base_ij\t$nbin\t$value\n";
	}}
	close POTENTIAL;

	$hash_argument{hash_potential} = \%hash_potential;
	
	
} else {
	
	&usage(); exit;
	
}



#############################################################
# Predict secondary structure 
#############################################################



my $test_struct = "$dir/test_fasta_struct";
STRUCT::fold_linux_RNA($test_fasta,$test_struct);




	

#############################################################
# Feature extraction
#############################################################



my $input_feature = "$dir/feature.txt";
open FEATURE,">$input_feature" or die;
FEATURE::feature_extraction($test_struct,\%hash_parameter,
\%hash_argument,\*FEATURE,'title'=>1,'label'=>"positive");
close FEATURE;





#############################################################
# Predict
#############################################################



my $outfile_out = "$dir/final_predict_result.txt";
LIBSVM::predict($libsvm_path,$input_feature,
$svm_model,$scaled_model,$outfile_out,$dir);

	
	





########################################
# End 
########################################


TOOLKIT::end($start_time);

print "\n\n============== End ============\n\n";






# =============================
# Load statistical potentials 
# =============================

sub load_potential {
	
	my $potential_file = shift or die;
	my $hash_potential = shift or die;
	my $hash_parameter = shift or die;
	
	open LOADIN,"<$potential_file" or die;
	while(my $line = <LOADIN>){
		chomp($line);	
		if($line =~/^#/){
			$line =~ s/^#//;
			my @tmp  = split/\s+/,$line;
			$hash_parameter->{$tmp[0]} = $tmp[1];
		}else{
			my @tmp  = split/\s+/,$line;
			my $base_ij = shift @tmp;
			my $bin     = shift @tmp;
			my $value   = shift @tmp;
			$hash_potential->{$base_ij}{$bin} = $value;
		}
	}
	close LOADIN;
	
}




# =============================
# Author: Yuangen Yao
# Date: 2013-2-20
# =============================

sub usage {
	
my $usage = << "USAGE";

Program: $0
Contact: Yao Yuangen <yygen89\@163.com>

Usage:

	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
		
	Options:
	
	--help		| -h : help
	--dir		| -D [dir] : the result directory
	--testFasta	| -T [file] : fasta file containing test sequences
	--potential	| -P [file] : file containing statistical potentials
	--svmModel	| -S [file] : file containing SVM model
	--scaleModel| -U [file] : file containing scaled model
	--posFasta	| -K [file]:  fasta file containing positive sequences 
	--negFasta	| -G [file]: fasta file containing negative sequences
	--libsvmPath| -L [path]: the absolute path to the Libsvm binary 
	--isEnergyF	| -J [0/1]:  If this option is set to 1, then using only knowledge-based energy features for prediction (default: 0). Note: the selected feature sets must match corresponding models.
	
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
	
	$0 [options] -T <test fasta> -S <SVM model> -U <scaled model> -P <potentials> -L <libsvm path>
	
	or
	
	$0 [options] -T <test fasta> -S <SVM model> -U <scaled model> -K <postive fasta> -G <negative fasta> -L <libsvm path>
		
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

USAGE
print $usage;

}
	

