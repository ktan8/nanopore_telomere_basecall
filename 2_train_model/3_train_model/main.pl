#!/usr/bin/perl

use strict;
use warnings;

# Train a new tuned model
system("bonito train -f --epochs 10 --lr 5e-4 --pretrained ~/kartong/software/bonito/bonito/models/dna_r9.4.1/ --directory ./ ./test3_final_run225_rightstrand");


