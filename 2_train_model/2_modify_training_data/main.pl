#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $dirname     = dirname(__FILE__);

# Modify the training data with the new ground truth
system("perl $dirname/modify_training_with_ground.working.pl run225.chunk.map_rawfasta.bam run225.raw.fasta run225.raw.groundtruth.txt > run225.chunk.groundtruth.txt");

# Get ground truth with right strand training data
system("samtools view -f0 -b run225.chunk.map_rawfasta.bam > run225.chunk.map_rawfasta.rightstrand.bam");

# Clean ground truth data
system("perl $dirname/clean.chunk_groundtruth.pl run225.chunk.groundtruth.rightstrand.txt > run225.chunk.groundtruth.rightstrand.clean.txt");

# Replace chunk data
system("python $dirname/replace_chunk_data.py run225.chunk.groundtruth.rightstrand.clean.txt tmp.chunks.npy tmp.references.npy ./tmp.reference_lengths.npy");
