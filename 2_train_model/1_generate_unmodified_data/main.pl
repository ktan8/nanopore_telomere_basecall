#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $fast5 = $ARGV[0];
my $bonito_model = $ARGV[1];
my $label = $ARGV[2];
my $ref_genome = $ARGV[3];




#~/kartong/software/bonito/bonito/models/dna_r9.4.1/


my $bonito = "bonito";
my $minimap2 = "minimap2";
my $dirname     = dirname(__FILE__);
my $extract_ground = $dirname . "/" . "extract_ground_truth.pl";
#my $extract_ground = "/PHShome/ktt19/code/telomere_basecalling/extract_ground_truth.pl";



my $fasta_raw = $label . ".raw.fasta";
#my $fasta_chunk = 

my $fasta_raw_refgenome_bam = $label . ".raw.map_ref.bam";
my $groundtruth = $label . ".raw.groundtruth.txt";
my $chunk_map_rawfasta_bam = $label . ".chunk.map_rawfasta.bam";
my $chunk_map_rawfasta_sam = $label . ".chunk.map_rawfasta.sam";


# 1 - Basecall the fast5 files of interest
system("$bonito basecaller --recursive $bonito_model $fast5 > $fasta_raw");


# 2 - Map the full length raw fasta to the reference genome
system("$minimap2 -a -x map-ont -t 8 $ref_genome $fasta_raw | samtools view -bS - > $fasta_raw_refgenome_bam");


# 3 - Extract the ground truth sequence from the reference genome
system("perl $extract_ground $fasta_raw_refgenome_bam $ref_genome > $groundtruth");


# 4 - Generate the chunk files for training with bonito
#   - Gives the numpy files
#   - Also gives the mapping of the chunk to the raw fasta

# For some weird reason, it cannot be directly piped or no .npy files
# will be generated.
#system("bonito basecaller --recursive --save-ctc --reference $fasta_raw $bonito_model $fast5 | samtools view -bS - > $chunk_map_rawfasta_bam");

system("bonito basecaller --recursive --save-ctc --reference $fasta_raw $bonito_model $fast5 > $chunk_map_rawfasta_sam");
system("samtools view -bS $chunk_map_rawfasta_sam > $chunk_map_rawfasta_bam");


# 5 - Extracts the ground sequence for the chunk that we need 

my $orig_chunk 	= "chunks.npy";
my $orig_ref	= "references.npy";
my $orig_reflen	= "reference_lengths.npy";
my $tmp_chunk	= "tmp.chunks.npy";
my $tmp_ref	= "tmp.references.npy";
my $tmp_reflen	= "tmp.reference_lengths.npy";

# Move the original chunk and reference files to a 
# temporary location
system("mv $orig_chunk $tmp_chunk");
system("mv $orig_ref $tmp_ref");
system("mv $orig_reflen $tmp_reflen");








