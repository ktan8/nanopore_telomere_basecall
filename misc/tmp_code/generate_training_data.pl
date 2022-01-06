#!/usr/bin/perl

use strict;
use warnings;

# 1) Basecall the fast5 data
# 2) Generate reference from the fasta generated
# 3) Use the 


my $fast5_dir 	= $ARGV[0];
my $model 	= $ARGV[1];
my $label 	= $ARGV[2];
my $fasta	= $label . ".fasta";
my $bam		= $label . ".bam";


system("bonito basecaller --recursive $model $fast5_dir > $fasta");
system("bonito basecaller --recursive --save-ctc --reference $fasta $model $fast5_dir | samtools view -bS - > $bam");


#bonito basecaller dna_r9.4.1 --save-ctc --reference reference.mmi /data/reads > /data/training/ctc-data/basecalls.sam
#bonito train --directory /data/training/ctc-data /data/training/model-dir



