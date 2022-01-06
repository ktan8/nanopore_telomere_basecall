#!/usr/bin/perl

use strict;
use warnings;

my $start = 50;
my $end = $start;
my $strand = "+";
#my $cigar = "40M60S";
my $cigar = "2S100M"; 

my ($start_new, $start_end) = convert_read_to_ref_coords($start, $end, $strand, $cigar);
print $start_new . "\n";


sub convert_read_to_ref_coords{
	my $start 	= $_[0];
	my $end		= $_[1];
	my $strand	= $_[2];
	my $cigar	= $_[3];

	if($strand eq "-"){
		$cigar = reverse_cigar($cigar);
	}
	
	my @cigar_int = ($cigar =~ /([0-9]+)[MIDNSHPX=]/g);
	my @cigar_char = ($cigar =~ /[0-9]+([MIDNSHPX=])/g);

	#my $ref_start = 0;
	#my $ref_end = 0;
	#my $query_start = 0;
	#my $query_end = 0;

	my $query_posn = 0;
	my $ref_posn = 0;
	my $start_reqr = -1;
	my $end_reqr = -1;	


	for(my $i=0; $i<scalar(@cigar_int); $i++){
		#my $ref_start_new 	= $ref_start;
		#my $ref_end_new 	= $ref_end;
		#my $query_start_new 	= $query_start;
		#my $query_end_new 	= $query_end;
		my $query_posn_new = 0;
		my $ref_posn_new = 0;


		# Consumes query
		if($cigar_char[$i] =~ /[MIS]/){
			#$query_start_new += $cigar_int[$i];
			#$query_end_new += $cigar_int[$i];
			$query_posn_new += $cigar_int[$i];
		}
		# Consumes reference
		if($cigar_char[$i] =~ /[MDN]/){
			$ref_posn_new += $cigar_int[$i];
		}

		# Check for start and ref position in reference coordinates
		if($start >= $query_posn && $start < $query_posn_new ){
			$start_reqr = $ref_posn + ($start - $query_posn);
		}
		if($end >= $query_posn && $end < $query_posn_new){
			$end_reqr = $ref_posn + ($end - $query_posn)
		}
	
		$query_posn = $query_posn_new;
		$ref_posn = $ref_posn_new;
	}
	# Deal with end of the query edge case
	if($query_posn == $end){
		$end_reqr = $ref_posn;
	}

	return $start_reqr, $end_reqr;
}


sub reverse_cigar{
	my $cigar = $_[0];
	my @cigar_int = ($cigar =~ /([0-9]+)[MIDNSHPX=]/g);
	my @cigar_char = ($cigar =~ /[0-9]+([MIDNSHPX=])/g);
	
	my $cigar_str = "";
	for(my $i=0; $i<scalar(@cigar_int); $i++){
		$cigar_str = $cigar_int[$i] . $cigar_char[$i] . $cigar_str;
	}
		
	return $cigar_str;	
}
