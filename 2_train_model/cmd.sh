conda activate bonito
module load openmpi/cuda/64/3.1.4

# Generate the data for run225
perl wrapper.eris_dsclust.pl /homes6/kartong/kartong/data/chm13_nanopore/bonito_telomere/promethion_telomere/data/fast5/225/ ~/kartong/software/bonito/bonito/models/dna_r9.4.1/ ./run225 ~/kartong/genome/chm13.draft_v1.0.fasta


perl modify_training_with_ground.working.pl run225.chunk.map_rawfasta.bam run225.raw.fasta run225.raw.groundtruth.txt > run225.chunk.groundtruth.txt

# Get ground truth with right strand training data
samtools view -f0 -b run225.chunk.map_rawfasta.bam > run225.chunk.map_rawfasta.rightstrand.bam

# Clean ground truth data
perl clean.chunk_groundtruth.pl run225.chunk.groundtruth.rightstrand.txt > run225.chunk.groundtruth.rightstrand.clean.txt

# Replace chunk data
python replace_chunk_data.py run225.chunk.groundtruth.rightstrand.clean.txt tmp.chunks.npy tmp.references.npy ./tmp.reference_lengths.npy


nohup bonito train -f --epochs 10 --lr 5e-4 --pretrained ~/kartong/software/bonito/bonito/models/dna_r9.4.1/ --directory ./ ./test3_final_run225_rightstrand &
