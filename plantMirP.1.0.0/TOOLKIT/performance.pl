#!/usr/bin/perl
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);


my $usage = "$0 

< total positive sample >
< total negative sample >
< ture positive >
< ture negative >

";

my $total_pos = shift;
my $total_neg = shift;
my $TP = shift;
my $TN = shift;


unless(defined $total_pos){
	print $usage;
	exit;
}

unless(defined $total_neg){
	print $usage;
	exit;
}

unless(defined $TP){
	print $usage;
	exit;
}

unless(defined $TN){
	print $usage;
	exit;
}



my %hash_perf = ();
my $str = &sn_sp_ac_mcc_1($total_pos,$total_neg,$TP,$TN,\%hash_perf);
print $str;





sub sn_sp_ac_mcc_1 {
	
	
	my ($total_pos,$total_neg,$TP,$TN,$hash_perf) = @_;
	
	
	my $FP = $total_neg - $TN;
	my $FN = $total_pos - $TP;
	
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
	
	if ( &is_numeric($sn) && &is_numeric($sp) && $sn * $sp >= 0 ){
		$Gm = sqrt( $sn * $sp );
		$Gm = &round($Gm,6);
	}
	
	if ( defined $hash_perf ) {		
		$hash_perf->{SN}  = $sn;
		$hash_perf->{SP}  = $sp;
		$hash_perf->{AC}  = $ac;
		$hash_perf->{PR}  = $pr;
		$hash_perf->{GM}  = $Gm;
		$hash_perf->{MCC} = $mcc;		
	}
	
	my $t = "";
	if( &is_numeric($sp) ){
		$t = &round(1-$sp,6);
	} else {
		$t = "NaN";
	}
	
	return "TOT_P:$total_pos\tTOT_N:$total_neg\tTP:$TP\tTN:$TN\tFP:$FP\tFN:$FN\nPR:$pr\tAC:$ac\tGM:$Gm\t1-SP:$t\tSN:$sn\tSP:$sp\tMCC:$mcc";
	
}





sub is_numeric {
	
    my $v = shift;

    if( looks_like_number( $v ) ){  	
        return 1;
    } else {
        return 0;
    }
    
}




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