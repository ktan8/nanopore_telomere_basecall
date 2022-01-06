#!/usr/bin/perl

use strict;
use warnings;

my $fast5 = $ARGV[0];
my $bonito_model = $ARGV[1];
my $label = $ARGV[2];
my $ref_genome = $ARGV[3];




#~/kartong/software/bonito/bonito/models/dna_r9.4.1/


my $bonito = "bonito";
my $minimap2 = "minimap2";
my $extract_ground = "/homes6/kartong/code/telomere_basecalling/extract_ground_truth.pl";



my $fasta_raw = $label . ".raw.fasta";
my $fasta_chunk = 

my $fasta_raw_refgenome_bam = $label . ".raw.map_ref.bam";
my $groundtruth = $label . ".raw.groundtruth.txt";
my $chunk_map_rawfasta_bam = $label . ".chunk.map_rawfasta.bam";


# 1 - Basecall the fast5 files of interest
system("$bonito basecaller --recursive $bonito_model $fast5 > $fasta_raw");


# 2 - Map the full length raw fasta to the reference genome
system("$minimap2 -a -x map-ont -t 8 $ref_genome $fasta_raw | samtools view -bS - > ");


# 3 - Extract the ground truth sequence from the reference genome
system("perl $extract_ground $fasta_raw_refgenome_bam $ref_genome > $groundtruth");


# 4 - Generate the chunk files for training with bonito
#   - Gives the numpy files
#   - Also gives the mapping of the chunk to the raw fasta
system("bonito basecaller --recursive --save-ctc --reference $fasta_raw $bonito_model $fast5 | samtools view -bS - > $chunk_map_rawfasta_bam");

# 5 - 


