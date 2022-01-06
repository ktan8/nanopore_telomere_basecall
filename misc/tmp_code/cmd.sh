

 perl extract_ground_truth.pl ~/kartong/data/chm13_nanopore/chm13.nanopore.rel7.map_chm13.sorted.term_10kb.bam


# Test - Generate ground truth for chunk
perl ~/code/telomere_basecalling/modify_training_with_ground.pl ~/tmp/train.bam ~/kartong/data/chm13_nanopore/bonito_telomere/chm13.nanopore.rel7.map_chm13.sorted.term_10kb.bam.groundtruth ~/tmp/train.fasta



# Generate the groundtruth for each chunk data
perl modify_training_with_ground.working.pl chm13_prom_telomere_train.chunk.map_rawfasta.bam chm13_prom_telomere_train.raw.fasta chm13_prom_telomere_train.raw.groundtruth.txt > chm13_prom_telomere_train.chunk.groundtruth.txt


