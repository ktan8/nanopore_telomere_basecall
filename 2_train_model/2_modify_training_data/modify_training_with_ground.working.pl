#!/ur/bin/perl


use strict;
use warnings;


# 1) chunk
# 2) raw
# 3) ground

my $chunk_bamfile = $ARGV[0];
my $raw_fasta	  = $ARGV[1];
my $ground_fasta  = $ARGV[2];

my %raw_fasta_hash = %{ parse_fasta_to_hash($raw_fasta) };
my %ground_fasta_hash = %{ parse_ground_to_hash($ground_fasta) };

open(my $BAM, "-|", "samtools view -F2304 $chunk_bamfile") || die $!;

my $mismapped_counter = 0;

while(my $line = <$BAM>){

	# Bam file from alignment of the chunk fasta to the
	# raw fasta
	chomp($line);
	my @lineArr 	= split(/\t/, $line);
	my $readname 	= $lineArr[0];
	my $chr 		= $lineArr[2];
	my $position 	= $lineArr[3];
	my $cigar 		= $lineArr[5];
	my $readsequence = $lineArr[9];


	# print "Read_position: " . $position . "\n";


	# Check if the read has been mapped to the correct 
	# fasta chromosome
	my ($readname_clean, $curr_chunk, $total_chunk) = split(/\:/, $readname);
	if($readname_clean ne $chr){	
		# Let's change it not to die, but rather keep a count of it
		#die "Readname and chrom name inconsistent";
		$mismapped_counter += 1;
	}	
	# print $position . "\n";
	# print $cigar . "\n";
	# print $ground_fasta_hash{$chr}->[1] . "\n";
	# print $ground_fasta_hash{$chr}->[2] . "\n";
	# print $ground_fasta_hash{$chr}->[3] . "\n";


	# Parse cigar tags for the mapping of the chunk fasta 
	# to the 
	my @cigar_int = ($cigar =~ /([0-9]+)[MIDNSHPX=]/g);
	my @cigar_char = ($cigar =~ /[0-9]+([MIDNSHPX=])/g);


	my ($ref_map_len, $left_softclip_len, $right_softclip_len) = cigar_len_info(\@cigar_int, \@cigar_char);

	my $ref_start = $position;
	my $ref_end = $position + $ref_map_len - 1;
	
	my ($ground_start, $ground_end) = 
	convert_read_to_ref_coords($ref_start, $ref_end, $ground_fasta_hash{$chr}->[1], $ground_fasta_hash{$chr}->[2]);

	# print "ground_coords" . "\n";
	# print $ref_start . "\n";
	# print $ref_end . "\n";
	# print $ground_fasta_hash{$chr}->[1] . "\n";
	# print $ground_fasta_hash{$chr}->[2] . "\n";
	# print $ground_start . "\n";
	# print $ground_end . "\n";


	my $ground_fasta_seq = $ground_fasta_hash{$chr}->[0];

	my $chunk_raw_fasta_seq	   = substr($raw_fasta_hash{$chr}, $ref_start - 1, $ref_map_len); # seq extracted from full fasta by mapping
	#my $chunk_ground_rough_fasta_seq = substr($ground_fasta_seq, $ref_start - 1, $ref_map_len);	
	my $chunk_ground_fasta_seq = substr($ground_fasta_seq, $ground_start - 1, $ground_end - $ground_start);
	#my $chunk_ground_fasta_seq = substr($ground_fasta_hash{$chr}->[0], $ref_start, $ground_end - $ground_start + 1);
	
	#my @result = ($readname, $readsequence, $chunk_raw_fasta_seq, $chunk_ground_rough_fasta_seq,  $chunk_ground_fasta_seq);
	my @result = ($readname, $readsequence, $chunk_raw_fasta_seq, $chunk_ground_fasta_seq);
	print join("\t", @result) . "\n";	
	#print $ground_fasta_seq . "\n";

}


print "Mismapped: " . $mismapped_counter . "\n";

close($BAM);

sub parse_fasta_to_hash{
	my $fasta_file = $_[0];
	my %fasta_hash;

	open(my $FASTA, $fasta_file) || die $!;
	while(my $line = <$FASTA>){
		my $name = $line;
		chomp($name);
		$name =~ s/^>//;
		my $sequence = <$FASTA>;
		chomp($sequence);
		$fasta_hash{$name} = $sequence;
	}
	close($FASTA);

	return \%fasta_hash;
}


sub cigar_len_info{
	my @cigar_int = @{ $_[0] };
	my @cigar_char = @{ $_[1] };
	my $ref_len = 0;
	my $softclip_left_len = 0;
	my $softclip_right_len = 0;

	for(my $i=0; $i < scalar(@cigar_int); $i++){	
		my $cigar_char_curr = $cigar_char[$i];
		my $cigar_int_curr = $cigar_int[$i];
		if($cigar_char_curr =~ /[MDN=X]/){
			$ref_len += $cigar_int_curr;
		}
	}

	if($cigar_char[0] eq "S"){
		$softclip_left_len = $cigar_int[0];
	}
	if($cigar_char[-1] eq "S"){
		$softclip_right_len = $cigar_int[-1];
	}

	return $ref_len, $softclip_left_len, $softclip_right_len;
}



sub parse_ground_to_hash{
	my $ground_file = $_[0];
	my %ground_hash;	

	open(my $ground, $ground_file) || die $!;
	while(my $line = <$ground>){
		chomp($line);
		my @lineArr = split(/\t/, $line);
		my $readname = $lineArr[0];
		my $strand = $lineArr[1];
		my $sequence = $lineArr[2];
		my $cigar = $lineArr[3];
		my $position = $lineArr[4];
		$ground_hash{$readname} = [$sequence, $strand, $cigar, $position];
		
	}
	close($ground);	
	
	return \%ground_hash;
}


sub convert_read_to_ref_coords{
	my $start 	= $_[0];
	my $end		= $_[1];
	my $strand	= $_[2];
	my $cigar	= $_[3];

	if($strand eq "-"){
		$cigar = reverse_cigar($cigar);
	}


	# Change softclip to match for cigar. This
	# follows the approach I applied in deeming the
	# softclip sequences as correct as they are telomeric.
	$cigar =~ s/S/M/g;

	
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
		#my $query_posn_new = 0;
		#my $ref_posn_new = 0;
		my $query_posn_new = $query_posn;
		my $ref_posn_new = $ref_posn;


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


		# print "========\n";
		# print $cigar_int[$i] . $cigar_char[$i] . "\n";
		# print $query_posn . "\n";
		# print $ref_posn . "\n";
		# print "=========\n";
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
