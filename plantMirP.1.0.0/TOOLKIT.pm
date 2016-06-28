package TOOLKIT;
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;


####################################################
################     基本工具     ##################
####################################################


sub mean {
	
	my @data = @_;
	
	my $n = 0;
	foreach my $i ( @data ){
		unless( &is_numeric($i) ){
			print STDERR "Error:$i not a numeric";
			die;
		}
		$n++;
	}
	
	if( $n == 0 ) {
		print STDERR "Error:Empty array\n";
		die;
	}
	
	my $total = 0;	
	for( my $i = 0; $i  < scalar @data; $i++ ){
		$total += $data[$i];
	}
	
	return  $total / $n;
}

####################################################

sub pearson {
	
	my ( $dataA, $dataB ) = @_;
	
	# check dim.
	if( $#$dataA != $#$dataB ){
		print STDERR "Error:Unequal dimension\n";
		die;
	}
	
	# calculate means.
	my $meanA = &mean( @$dataA );
	my $meanB = &mean( @$dataB );
	
	# calculate numerator.
	my $sum = 0;
	for( my $i = 0; $i < scalar @$dataA; $i++ ){
		$sum += ( $dataA->[$i] - $meanA ) * ( $dataB->[$i] - $meanB );
	}
	my $up = $sum;
	
	# calculate denominator.
	$sum = 0;
	for( my $i = 0; $i < scalar @$dataA; $i++ ){
		$sum += ( $dataA->[$i] - $meanA ) * ( $dataA->[$i] - $meanA );
	}
	my $sqrtA = sqrt( $sum );	
	
	$sum = 0;
	for( my $i = 0; $i < scalar @$dataB; $i++ ){
		$sum += ( $dataB->[$i] - $meanB ) * ( $dataB->[$i] - $meanB );
	}
	my $sqrtB = sqrt( $sum );
	
	# check whether denominator is 0.
	if( $sqrtA <= 1.e-100 || $sqrtB <= 1.e-100 ){
		return 0;
	} else {
		return $up / ( $sqrtA * $sqrtB );
	}
	
}

####################################################

sub median {
	
	my $usage = "median < reference of numeric array >";
    my $rpole = shift or die $usage;
    
    my @pole = @$rpole;
    
    my $ii = 0;
    foreach my $i ( @pole ){
	    unless ( &is_numeric($i) ){
			print STDERR "Error:$i not a numeric\n";
			die;
		}
    	$ii++;
    }
    
    if( $ii == 0 ){
    	print STDERR "Error:Empty array";
    	die;
    }
    
   	@pole = sort( {$a<=>$b} @pole);
   	
   	my $ret;
    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }
    
    return $ret;
    
}

####################################################

sub shuffle {
	
	my $usage = " shuffle < reference of array >";
	my $array = shift or die $usage;
	
	for (my $j = @$array-1; $j > 0; --$j) {
  		my $i = int(rand($j));
   		my $t = $array->[$j]; $array->[$j] = $array->[$i]; $array->[$i] = $t;
	}
}



####################################################

sub getAbsPath {
	
	my $usage = "< directory or file >";
	my $file_or_dir = shift or die $usage;
	
	my $dir = dirname($file_or_dir);
	my $basename = basename($file_or_dir);
	
	my $current_dir = getcwd;
	my $child_dir = "$current_dir/$dir";
	
	if($dir eq "."){
		return "$current_dir/$basename";
	} elsif( $dir eq ".." ){
		chdir("../");
		my $parent_dir = getcwd;
		chdir($current_dir);		
		return "$parent_dir/$basename";
	} elsif ( -d $child_dir ){
		return "$child_dir/$basename";
	} else {
		return $file_or_dir;
	}
	
}

####################################################

sub toMatrix {
	
	my $infile 	 = shift or die;
	my $hash_row = shift or die;
	
	
	my $name;
	while (@_) {
		my $argument = shift @_;
		if ($argument=~/name/i) {$name=shift @_}
	}
	

	unless(-e $infile){
		print STDERR "Error:$infile not exists\n";
		die;
	}
	
	open IN,"<$infile" or die;
	
	my $flag = 1;
	my @columns = ();
	my %hash_col_tmp = ();
	
	while( defined( my $line = <IN> )){
		chomp($line);
		next if $line =~/^\s*$/;
		
		my @tmp = split/\t/,$line;
		my $row = shift @tmp;
		
		if ( $flag ) {
			if ( defined $name ){
				unshift(@tmp,$row);
			} else {	
				push(@columns,$row);
			}
			for my $i ( 0 .. $#tmp ){		
				unless( exists $hash_col_tmp{ $i } ){
					$hash_col_tmp{ $i } = $tmp[$i];
				}	
				push(@columns,$tmp[$i]);		
			}
			$flag = 0;
			next;
		} 
		
		for my $i ( 0 .. $#tmp ){
			my $col = $hash_col_tmp{ $i };
			unless( exists $hash_row->{ $row }{ $col } ){	
				$hash_row->{ $row }{ $col } = $tmp[$i];
			} else {
				print STDERR "Error:($row,$col) repeat\n";
				die;
			}
		}
	}
	close IN;
	
	return @columns;
	
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

sub logN { 
	
	my $n = shift;
	my $base = shift or die;
	
	if ( $n <= 0 ){
		print STDERR "Error:$n <= 0\n";
		die;
	} else {
		return log($n)/log($base);    
	}
}


####################################################

sub permutation {
	
	my ($in_ref_array,$out_ref_array,$n,$prefix) = @_;
	
	unless( defined $prefix ){
		$prefix = "";
	}
	
	if( $n == 0 ){
		for my $r ( @$in_ref_array ){
			push (@$out_ref_array,$prefix.$r);
		}
      	return;
   }
   
   for my $r ( @$in_ref_array ){
		&permutation($in_ref_array,$out_ref_array,$n-1,$prefix.$r);
	}
}


####################################################

sub factorial {
	
	my $n = shift;
	
	unless( &is_numeric($n) ){
		print STDERR "Error:$n not a numeric\n";
		die;
	}
	
	if ( $n <= 0 ){
		print STDERR "Error:$n <= 0\n";
		die;
	}
	
	my $nj = 1; 
	until( $n == 0 ) { 
		$nj *= $n; 
		$n--; 
	} 
	
	return $nj;
	
}


####################################################

sub matchPattern {
	
	my $nmismatch 	= shift;
	my $pattern 	= shift or die;
	my $subjectSeq 	= shift or die;
	my $hash_pos 	= shift or die;
	
	unless( defined $nmismatch ){
		print STDERR "Error:number of mismatch not defined\n";
		die;
	}
	
	my @pattern 	= split '', $pattern;
	my @subjectSeq 	= split '', $subjectSeq;

	my $lng = length( $pattern );
	my %err_count;
	for my $i ( 0 .. @subjectSeq - $lng ) {
		my $n_err = 0 ;
		for my $j ( 0 .. @pattern - 1 ) {
			$n_err++ if( $pattern[$j] ne $subjectSeq[$i+$j] );
		}
		next if ( $n_err > $nmismatch );
		$err_count{$i} = $n_err;
	}

	foreach( keys %err_count ) {
		if( $err_count{$_} == $nmismatch ){
			my $value = substr( $subjectSeq, $_, $lng );		
			$hash_pos->{ $_ } = $value;
		}
	}
}

####################################################

sub readfasta {
	
	my $usage = "< fasta file > < hash reference >";
	my $infile = shift or die $usage;
	my $hash_seq = shift or die $usage;
	
	unless(-e $infile){
		print STDERR "Error:$infile not exists";
		die;
	}
	
	open IN, $infile || die;
	
	my $c=0;
	my $seqId;
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}


####################################################

sub start {
	
	my $out = shift; 
	
    my($second, $minute, $hour, $dayOfMonth, $month, 
    $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $second = "0$second" if($second =~ /^\d$/);
    my $sTime = "$hour:$minute:$second";
    my $stime = time;
    my $start_time = localtime;
    
    print "\n\nstarted: $start_time\n\n";
    if( defined $out ){
    	print $out "started: $start_time\n";
    }
    
    return $stime;
}


####################################################

sub end {
	
	my $stime = shift or die;
	my $out = shift;
	
    my $etime = time - $stime;
    my ($second, $minute, $hour, $dayOfMonth, $month, 
    $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    $second = "0$second" if($second =~ /^\d$/);
    my $eTime = "$hour:$minute:$second";
    
    my $end_time = localtime;
    
    print "\n\nended: $end_time\n";
    print "total:", int($etime / 3600),"h:",
    int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
    
    if( defined $out ){
    	print $out "\n\nended: $end_time\n";
    	print $out "total:", int($etime / 3600),"h:",
    	int(($etime % 3600) / 60),"m:",int($etime % 60),"s\n\n";
    }
}


####################################################
1;
