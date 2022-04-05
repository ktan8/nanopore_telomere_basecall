#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my $dirname     = dirname(__FILE__);
my $input_fasta = $ARGV[0];
my $fast5_dir	= $ARGV[1];
my $label 	= $ARGV[2];
my $problematic_fasta = $label . ".problematic";
my $problematic_readnames = $label . ".problematic.repeatcounts.filtered.readname";
my $extracted_fast5_dir = $label . "/";
my $corrected_telomere_fasta = $label . ".correctedreads.fasta.gz";
my $final_output_fasta = $label . ".telomerefixed.fasta.gz";


system("bash -c 'perl $dirname/2_identify_problematic_reads/main.pl <(zcat -f $input_fasta) $problematic_fasta'");

system("perl $dirname/3_basecall_problematic_reads/main.pl $fast5_dir $problematic_readnames $extracted_fast5_dir $corrected_telomere_fasta");

system("perl $dirname/4_replace_problematic_reads/find_and_replace_fasta.pl $input_fasta $corrected_telomere_fasta | gzip > $final_output_fasta");

