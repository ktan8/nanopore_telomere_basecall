#!/usr/bin/perl

use strict;
use warnings;

my $fastafile = $ARGV[0];

my @motifs = ("TTAGGG", "CCCTAA", "TTAAAA", "TTTTAA", "CTTCTT", "AAGAAG", "CCCTGG", "CCAGGG");
my $motif_count = 3;
my @motif_long;

for my $motif(@motifs){
	my $motif_long_str = $motif x $motif_count;
	push(@motif_long, $motif_long_str);
}

my $regex = join("|", @motif_long);


open(my $FASTA, $fastafile) || die $!;
while(my $readname = <$FASTA>){
	chomp($readname);
	my $readname_clean = substr($readname, 1);
	my $fastaseq = <$FASTA>;
	chomp($fastaseq);
	my $motif_count = () = $fastaseq =~ /$regex/gi;
	print $readname_clean . "\t" . $motif_count . "\n";
}
close($FASTA);

