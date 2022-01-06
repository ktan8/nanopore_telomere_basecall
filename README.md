
# Nanopore telomere basecalling

This package contains of two parts. The first part is the pipeline for retraining a model to basecall and rectify basecalling errors at telomeres. The second part consists of a pipeline that can be directly applied to basecall and fix basecallinge errors at telomeric regions.


## Appling tuned model to the telomeric nanopore reads

The 

1_bonito_basecalling_model
2_identify_problematic_reads
3_basecall_problematic_reads
4_replace_problematic_reads

### 1_bonito_basecalling_model
This directory contains the tuned basecalling model for bonito. The model can be downloaded from the following path and unzip into this folder.
https://zenodo.org/api/files/86cb9586-300f-493d-b9c4-0ab2f2848e3c/chm13_nanopore_trained_run225.zip

### 2_identify_problematic_reads
This directory contains a set of scripts used to identify the problematic telomeric reads for basecalling. Specifically, the scripts will identify long-reads with a high freqeuncy of telomeric repeats and telomeric repeat artefacts to redo the base calls. A list of readnames corresponding to the candidate reads will then be produced

To apply this step, run the following command:
```
perl main.pl <input_fasta> <output_label>
```

### 3_basecall_problematic_reads
This directory contains a set of scripts to extract the fast5 for the required reads. Subsequent to that, we will basecall these reads using the tuned bonito model.

To apply this step, run the following command:
```
perl main.pl <directory_of_all_fast5_files> <file_of_readnames_for_required_reads> <directory_to_output_extracted_fast5> <output_fastagz_file>
```

### 4_replace_problematic_reads
This step allows one to replace the reads with the basecalling repeat artefacts with the reads that have been fixed.

To apply this step, run the following command:
```
perl find_and_replace_fasta.pl <original_full_fasta_file> <rebasecalled_fasta_file>
```




## Tuning Nanopore bonito basecalling model for telomeric reads



