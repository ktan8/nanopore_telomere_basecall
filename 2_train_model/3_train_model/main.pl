#!/usr/bin/perl

use strict;
use warnings;

my $original_basecalling_model 	= $ARGV[0];
my $training_data_directory 	= $ARGV[1];
my $output_model_name 	= $ARGV[2];

# Train a new tuned model
system("bonito train -f --epochs 10 --lr 5e-4 --pretrained $original_basecalling_model --directory $training_data_directory $output_model_name");


