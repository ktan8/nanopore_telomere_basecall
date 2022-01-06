#!/usr/bin/perl


# Extract the ground truth sequence from
# the telomeric reads


use strict;
use warnings;

my $bam_file = $ARGV[0];
my $refgenome = $ARGV[1];

my $fai_file = $refgenome . ".fai";
my %fai_hash = %{ parse_fai($fai_file) };


open(my $BAM, "-|", "samtools view -F2304 $bam_file") || die $!;
while(my $line = <$BAM>){
	chomp($line);
	my @lineArr = split(/\t/, $line);
	my $readname 	= $lineArr[0];
	my $flag	= $lineArr[1];
	my $chr         = $lineArr[2];
	my $position    = $lineArr[3];
	my $cigar 	= $lineArr[5];
	my $sequence 	= $lineArr[9];
	my $seq_len 	= length($sequence);
	

	my @cigar_int = ($cigar =~ /([0-9]+)[MIDNSHPX=]/g);
	my @cigar_char = ($cigar =~ /[0-9]+([MIDNSHPX=])/g);

	my ($ref_len, $left_softclip_len, $right_softclip_len) = cigar_len_info(\@cigar_int, \@cigar_char);
	
	my $ref_start = $position - $left_softclip_len;
	my $ref_end = $position + $ref_len + $right_softclip_len - 1; # need to substract 1 as both start and end are inclusive values
	
	#print join("\t", ($cigar, $ref_len, $left_softclip_len, $right_softclip_len, $ref_start, $ref_end)) . "\n";
	

	my $ground_truth_seq = extract_fasta($chr, $ref_start, $ref_end);	
	
	my $strand = "+";
	if($flag & 16){
		$strand  = "-";
		$ground_truth_seq = reverse_complement($ground_truth_seq);
	}

	print join("\t", $readname, $strand, $ground_truth_seq, $cigar, $position) . "\n";
}
close($BAM);


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


sub parse_fai{
	my $fai_file = $_[0];
	
	open(my $FAI, $fai_file) || die ".fai file does not exist";
	my %fai_hash;
	while(my $line = <$FAI>){
		chomp($line);
		my @lineArr = split(/\t/, $line);
		my $chr = $lineArr[0];
		my $chr_size = $lineArr[1];
		$fai_hash{$chr} = $chr_size;	
	}
	close($FAI);

	return \%fai_hash;
}


sub extract_fasta{
	my $chr 	= $_[0];
	my $start 	= $_[1];
	my $end 	= $_[2];

	my $pad_left = 0;
	my $pad_right = 0;
	
	if($start < 1){
		$pad_left = 0 - $start + 1;
	}
	if($end > $fai_hash{$chr}){
		$pad_right = $end - $fai_hash{$chr};
	}
	
	my $region = $chr . ":" . $start . "-" . $end;
	my $sequence = `samtools faidx -n 1000000000 $refgenome $region | tail -n 1`;
	chomp($sequence);
	my $pad_left_seq = "N" x $pad_left;
	my $pad_right_seq = "N" x $pad_right;
	my $sequence_final = $pad_left_seq . $sequence . $pad_right_seq;
	
	return $sequence_final;
}


sub reverse_complement{
	my $seq = $_[0];
	my $revcomp = reverse $seq;
	$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;

	return $revcomp;
}
