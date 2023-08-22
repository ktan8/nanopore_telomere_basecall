#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $dirname     = dirname(__FILE__);
my $training_data_label = $ARGV[0];


# Get ground truth with right strand training data
system("samtools view -f0 -b $training_data_label.chunk.map_rawfasta.bam > $training_data_label.chunk.map_rawfasta.rightstrand.bam");

# Modify the training data with the new ground truth
system("perl $dirname/modify_training_with_ground.working.pl $training_data_label.chunk.map_rawfasta.rightstrand.bam $training_data_label.raw.fasta $training_data_label.raw.groundtruth.txt > $training_data_label.chunk.groundtruth.rightstrand.txt");

# Clean ground truth data
system("perl $dirname/clean.chunk_groundtruth.pl $training_data_label.chunk.groundtruth.rightstrand.txt > $training_data_label.chunk.groundtruth.rightstrand.clean.txt");

# Replace chunk data
system("python $dirname/replace_chunk_data.py $training_data_label.chunk.groundtruth.rightstrand.clean.txt tmp.chunks.npy tmp.references.npy ./tmp.reference_lengths.npy");
