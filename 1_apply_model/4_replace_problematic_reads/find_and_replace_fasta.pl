#!/usr/bin/perl

use strict;
use warnings;

my $original_fasta 	= $ARGV[0];
my $rebasecalled_fasta 	= $ARGV[1];
#my $fixed_fasta		= 


open(my $FIXED_READS, "|-", "zcat $rebasecalled_fasta") || die $!;
my %fixed_reads_hash;
while(my $readname = <$FIXED_READS>){
	chomp($readname);
	my $readname_clean = substr($readname, 1);
	my $sequence = <$FIXED_READS>;
	chomp($sequence);
	$fixed_reads_hash{$readname_clean} = $sequence;
}
close($FIXED_READS);



open(my $ORIGINAL_READS, "|-", "zcat $rebasecalled_fasta") || die $!;
while(my $readname = <$ORIGINAL_READS>){
	chomp($readname);
	my $readname_clean = substr($readname, 1);
	my $sequence = <$ORIGINAL_READS>;
	chomp($sequence);

	# Replace with fixed reads if name matches
	if(exists($fixed_reads_hash{$readname_clean})){
		my $fixed_read_curr_sequence = $fixed_reads_hash{$readname_clean};
		print ">" . $readname_clean . "\n";
		print $fixed_read_curr_sequence . "\n";
	}
	# Use original reads if does not match
	else{
		print ">" . $readname_clean . "\n";
		print $sequence . "\n";
	}
}
close($ORIGINAL_READS);


