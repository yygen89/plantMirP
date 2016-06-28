#!/usr/bin/perl
use strict;
use warnings;




my $class = "Viridiplantae";
my $miRBase_hairpin_fasta = "hairpin.fa";
my $organism_file = "organisms.txt";



my %hash_organism = ();
&toMatrix($organism_file,\%hash_organism);



my @species_list = ();
open OUT,">$class-organism.txt" or die;
for my $species ( sort keys %hash_organism ){
	my $tree = $hash_organism{$species}{"#tree"};
	my $Latin_name = $hash_organism{$species}{"#name"};
	if( $tree =~ /$class/i ){
		push(@species_list,$species);
		print OUT "$species\t$Latin_name\n";
	}
}
close OUT;
		
		
		
my $outfile_hairpin_fasta = "$class-hairpin.fa";
&get_miRNA_precursor(\@species_list,$miRBase_hairpin_fasta,$outfile_hairpin_fasta);



print "\n\n################# END #################\n\n";





########################################################

sub get_miRNA_precursor {
	
	my $usage = "$0 <specis name in miRBase (three letters) > ";
	
	my $input_species = shift or die $usage;
	my $infile_hairpin_fasta = shift or die;
	my $outfile_hairpin_fasta = shift or die;
	

	my %hash_id = ();
	my %hash_id2 = ();
	
	
	open HAIRPIN,">$outfile_hairpin_fasta" or die;
	for my $species ( @$input_species ){
		
		&trim($species);

		print "species: $species\n";


		my %hash_hairpin = ();
		&readfasta($infile_hairpin_fasta,\%hash_hairpin);

		
		for my $id ( keys %hash_hairpin ){
			
			if( $id =~/^$species-/i ){
				
				my $seq = $hash_hairpin{$id};
				$seq = uc($seq);
				$seq =~ s/T/U/ig;
								
				unless(not ($seq =~ /^([A|C|G|U]+)$/)){
				
					my @tmp = split/\s+/,$id;
					my $new_id = shift @tmp;
					unless(exists $hash_id2{$new_id}){
						$hash_id2{$new_id} = 1;
						print HAIRPIN ">$new_id\n$seq\n";
					}else{
						print STDERR "Error: $new_id repeat\n";
						die;
					}
						
				}
			}
		}
		
	}
	close HAIRPIN;
	
}




########################################################

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
	my %hash_id_check = ();
	while (defined (my $line=<IN>)) {
		chomp($line);
		next if $line=~/^\s*$/;
		if ($line =~/^>/) {
			$line =~s/\s*$//g;	
			$line =~s/^\s*//g;	
			$seqId = substr($line,1);
			$c++;
			
			unless(exists $hash_id_check{$seqId}){
				$hash_id_check{$seqId} = 1;
			}else{
				print STDERR "Error:$seqId repeat\n";
				die;
			}
				
		} else {
			$line =~s/\s*//g;
			$hash_seq->{$seqId}.=$line;
		}
	}
	close IN;
	
	return $c;
	
}


########################################################

sub trim {
	my $str = shift or die;
	$str =~ s/^\s*//g;
	$str =~ s/\s*$//g;
	return $str;
}


########################################################

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
	
	my %hash_check_col = ();
	my %hash_check_row = ();
	
	while( defined( my $line = <IN> )){
		chomp($line);
		next if $line =~/^\s*$/;
		
		my @tmp = split/\t/,$line;
		my $row = shift @tmp;
		
		
		unless(exists $hash_check_row{$row}){
			$hash_check_row{$row} = 1;
		}else{
			print STDERR "Error: $row repeat\n";
			die;
		}
			
		
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
			
			for my $col_t ( @columns ){
				unless(exists $hash_check_col{$col_t}){
					$hash_check_col{$col_t} = 1;
				}else{
					print STDERR "Error: $col_t repeat\n";
					die;
				}
			}
			
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

########################################################