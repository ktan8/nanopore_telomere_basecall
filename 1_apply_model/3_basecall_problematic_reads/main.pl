#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $dirname     = dirname(__FILE__);
my $fast5_directory 		= $ARGV[0];
my $basecalling_model 		= $dirname . "/../1_bonito_basecalling_model/chm13_nanopore_trained_run225/";
my $required_readnames 		= $ARGV[1];
my $required_fast5_directory 	= $ARGV[2];
my $output_fastagz 		= $ARGV[3];

# Make directory for the fast5 files


# Extract the fast5 files required
system("fast5_subset -i $fast5_directory -s $required_fast5_directory -l $required_readnames --recursive");


# Basecall problematic reads with bonito
system("bonito basecaller --recursive $basecalling_model $required_fast5_directory | gzip > $output_fastagz");


