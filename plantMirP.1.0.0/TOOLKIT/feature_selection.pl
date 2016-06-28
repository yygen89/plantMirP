#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

use TOOLKIT;
use LIBSVM;
use FEATURE_SELECTION;


my $usage = "$0

< training_input_feature >
< predict_input_feature > 
< keep_feature_num > 
< max_similarity >

";

my $training_input_feature = shift or die $usage;
my $predict_input_feature = shift or die $usage;
my $keep_feature_num = shift or die $usage;
my $max_similarity = shift or die $usage;



#set this parameter to the path to the Libsvm binary
my $libsvm_path = "/home/yyg/libsvm/";


#############################################################
# Start
#############################################################


my $dir = getcwd;
my $start_time = TOOLKIT::start();



#############################################################
# Feature selection
#############################################################


my $relevance_measure_type = "MM";
my $keep_feature_name = "$dir/keep_feature_name.txt";
my $selected_feature_for_training = "$dir/selected_feature_for_training.txt";
my ($all_feature_num, $selected_feature_num) = FEATURE_SELECTION::feature_seletcion(
$training_input_feature, $selected_feature_for_training, $keep_feature_name,
$keep_feature_num, $relevance_measure_type, $max_similarity);
	

my $selected_feature_for_predict = "$dir/selected_feature_for_predict.txt";
FEATURE_SELECTION::load_feature($keep_feature_name,
$predict_input_feature,$selected_feature_for_predict);


print "\n\n********* All_feature_num: $all_feature_num *********\n";
print "********* Selected_feature_num: $selected_feature_num *********\n\n";




#############################################################
# Training
#############################################################


print "\n\n******** Training ********\n\n";
my $svm_model = "$dir/svm_model.txt";
my $scaled_model = "$dir/scaled_model.txt";
LIBSVM::training($libsvm_path,$selected_feature_for_training,
$svm_model,$scaled_model,$dir);

	


#############################################################
# Predict
#############################################################


print "\n\n******** Predict ********\n\n";
my $outfile_out = "$dir/final_predict_result.txt";
LIBSVM::predict($libsvm_path,$selected_feature_for_predict,
$svm_model,$scaled_model,$outfile_out,$dir);






########################################
# End 
########################################

TOOLKIT::end($start_time);

print "\n\n============== End ============\n\n";


