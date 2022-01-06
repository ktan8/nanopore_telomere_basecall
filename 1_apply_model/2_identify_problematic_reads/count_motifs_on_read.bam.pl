#!/usr/bin/perl

use strict;
use warnings;

my $bamfile = $ARGV[0];

my @motifs = ("TTAGGG", "CCCTAA", "TTAAAA", "TTTTAA", "CTTCTT", "AAGAAG", "CCCTGG", "CCAGGG");
my $motif_count = 3;
my @motif_long;

for my $motif(@motifs){
	my $motif_long_str = $motif x $motif_count;
	push(@motif_long, $motif_long_str);
}

my $regex = join("|", @motif_long);


open(my $BAM, "-|", "samtools view -h -F256 $bamfile") || die $!;
while(my $line = <$BAM>){
	if($line =~ /^\@/){
		#print $line;
		next;
	}
	chomp($line);
	my @lineArr = split(/\t/, $line);
	my $sequence = $lineArr[9];
	my $motif_count = () = $sequence =~ /$regex/gi;
	print $lineArr[0] . "\t" . $motif_count . "\n";
	
		#print join("\t", @lineArr) . "\n";

}
close($BAM);

