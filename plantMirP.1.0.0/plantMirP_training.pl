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
use KnowledgeBasedEnergy;
use CROSS_VALIDATION;


	
# parameters
my $dir;
my $help;	
my $positive_fasta;
my $negative_fasta;
my $nbin = 20;
my $rbin = 1;
my $adjacent = 5;

#Cross validation
my $cv = 4;
my $round_num = 1;



GetOptions( 	
	"help|h!"			=> \$help,
	"dir|D=s"			=> \$dir,
	"posFasta|P=s"		=> \$positive_fasta,
	"negFasta|N=s"		=> \$negative_fasta,
	"nbin|B=i"			=> \$nbin,
	"rbin|R=i"			=> \$rbin,
	"adjacent|A=i"		=> \$adjacent,
	"nfold|F=s"			=> \$cv,
	"round|X=i"			=> \$round_num,
);



#set this parameter to the path to the Libsvm binary
my $libsvm_path = "/home/yyg/libsvm/";


# parameter for feature extraction 
my $is_only_knowledge_based_energy_features = 0;



########################################
# start 
########################################


if ( $help ) {
	&usage(); exit;
}


if( defined $positive_fasta ){
	$positive_fasta = TOOLKIT::getAbsPath($positive_fasta);
}else{
	&usage(); exit;
}


if( defined $negative_fasta ){
	$negative_fasta = TOOLKIT::getAbsPath($negative_fasta);	
}else{
	&usage(); exit;
}


my $current_dir = getcwd;
unless( defined $dir ){
	$dir = "$current_dir/training";
} else {
	$dir = TOOLKIT::getAbsPath($dir);
}
mkdir($dir);


my $start_time = TOOLKIT::start();

#############################################################
# Save all parameters
#############################################################


my %hash_parameter = ();

# parameters for CV
$hash_parameter{cv} = $cv;
$hash_parameter{round_num} = $round_num;


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



#############################################################
# Predict secondary structure 
#############################################################



my $positive_struct = "$dir/positive_fasta_struct";
STRUCT::fold_linux_RNA($positive_fasta,$positive_struct);


my $negative_struct = "$dir/negative_fasta_struct";
STRUCT::fold_linux_RNA($negative_fasta,$negative_struct);



	
#############################################################
# Calculate statistical potentials 
#############################################################


my %hash_potential = ();
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




#############################################################
# Feature extraction
#############################################################



my $input_feature = "$dir/feature.txt";
open FEATURE,">$input_feature" or die;
FEATURE::feature_extraction($positive_struct,
\%hash_parameter,\%hash_argument,
\*FEATURE,'title'=>1,'label'=>"positive");


FEATURE::feature_extraction($negative_struct,
\%hash_parameter,\%hash_argument,\*FEATURE,
'label'=>"negative");
close FEATURE;



#############################################################
# Calculate F score
#############################################################


my %hash_fscore = ();
FEATURE::F_score($input_feature,\%hash_fscore);

my $outfile_fscore = "$dir/F_score.txt";
open FSCORE,">$outfile_fscore" or die;
for my $feature ( sort {$hash_fscore{$b}<=>$hash_fscore{$a}} keys %hash_fscore ){
	print FSCORE "$feature\t$hash_fscore{$feature}\n";
}
close FSCORE;




#############################################################
# Training
#############################################################



my $svm_model = "$dir/svm_model.txt";
my $scaled_model = "$dir/scaled_model.txt";
LIBSVM::training($libsvm_path,$input_feature,$svm_model,$scaled_model,$dir);

	
	
	

#############################################################
# Cross-validation
#############################################################

if( $cv ){
	my $cv_dir = "$dir/cv";
	mkdir($cv_dir);
	
	my $cv_out = "$cv_dir/out.txt";
	my $cv_roc = "$cv_dir/roc.txt";
	CROSS_VALIDATION::cross_validation($positive_fasta,$negative_fasta,
	$cv_roc,$cv_out,$cv_dir,\%hash_parameter,\%hash_argument,$libsvm_path);
}




########################################
# End 
########################################

TOOLKIT::end($start_time);

print "\n\n============== End ============\n\n";






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
	--negFasta	| -N [file] : fasta file containing negative sequences
	--posFasta	| -P [file] : fasta file containing positive sequences
	--nbin		| -B [int] : the number of bin (default:20)
	--rbin		| -R [int] : the radius of bin (default:1)
	--adjacent	| -A [int] : the number of adjacent for statistics potentials (default:1 2 3 4 5)
	--nfold		| -F [int] : nfold cross validation (default:4)
	--round		| -X [int] : the number of round for cross-validation (default:1)
	
	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
	
	$0 [options] -P <postive fasta> -N <negative fasta>

	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

USAGE
print $usage;

}
	

