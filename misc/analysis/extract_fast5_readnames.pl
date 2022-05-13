#!/usr/bin/perl

use strict;
use warnings;

my $folder = $ARGV[0];

my $files = `ls $folder`;
my @files_arr = split(/\n/, $files);

foreach my $file (@files_arr){
	my $fullfile = $folder . "/" . $file;
	#my $readname = `h5dump -g \"Raw/Reads/\" $fullfile | grep -A11 H5T_STRING | grep '(0)'`;
	my $readname = `h5dump -g \"Raw/Reads/\" $fullfile | grep -A12 "read_id" | grep '(0)'`;

	my @readname_arr = split(/\"/, $readname);
	print $fullfile . "\t" . $readname_arr[1] . "\n";
}
