#!/usr/bin/perl

use strict;
use warnings;

my $file = $ARGV[0];
my $cutoff = 0.05;

open(my $FILE, $file) || die $!;

while(my $line = <$FILE>){
	chomp($line);
	my @lineArr = split(/\t/, $line);
	my $length1 = length($lineArr[1]);
	my $length2 = length($lineArr[2]);
	
	if(!defined($lineArr[3])){
		$lineArr[3] = "N";
	}
	my $length3 = length($lineArr[3]);

	my $diff = abs(($length3 / $length2) - 1);
	if($diff < $cutoff){
		#print join("\t", $length1, $length2, $length3, $diff) . "\n";
		print join("\t", @lineArr) . "\n";
	}
}

close($FILE);
