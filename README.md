
# Nanopore telomere basecalling

This package contains of two parts. The first part is the pipeline for retraining a model to basecall and rectify basecalling errors at telomeres. The second part consists of a pipeline that can be directly applied to basecall and fix basecallinge errors at telomeric regions.


## Appling tuned model to the telomeric nanopore reads

There are a series of four steps to apply the tuned basecalling model to the telomeric nanopore reads. These steps and the corresponding scripts can be found in the following directories.

### Dependencies
To apply the bonito basecalling model, you will need the following software in your environment.
1. Python
2. Perl
3. fast5_subset (included as part of the ont_fast5_api package at https://github.com/nanoporetech/ont_fast5_api)
4. Bonito basecaller (current pipeline was tested on Bonito v0.3.5 which can be obtained from https://github.com/nanoporetech/bonito/releases). Please also ensure that the Bonito basecaller is working before applying this pipeline as Bonito is dependent on a number of other packages (e.g. CuPy, CUDA, etc.). Please refer to the Bonito repository for detailed information on the required packages.
5. Samtools (http://www.htslib.org/)

### Applying full pipeline
To apply the full pipeline in a single step, one can use the command
```
perl 1_apply_model/fullpipeline.pl <input_fasta> <fast5_directory_of_nanopore_signal_data> <output_label>
```

<input_fasta> - Pre-called fasta files that you can generate using either Guppy or the Bonito basecaller. (This pipeline only re-basecalls the telomeric reads. So you will still need to generate fasta files from your raw fast5 files using either Guppy or the default Bonito caller)
<fast5_directory_of_nanopore_signal_data> - This specifies the folder where your fast5 files are. The required fast5 files to re-basecall are extracted from this directory
<output_label> - Any name or output path that you so desire.


Otherwise, the pipeline can also be applied by following each of the following steps.


### 1. bonito_basecalling_model
This directory contains the tuned basecalling model for bonito. The model can be downloaded from the following path (https://zenodo.org/record/5819148/files/chm13_nanopore_trained_run225.zip?download=1) and unzip into this folder.


### 2. identify_problematic_reads
This directory contains a set of scripts used to identify the problematic telomeric reads for basecalling. Specifically, the scripts will identify long-reads with a high freqeuncy of telomeric repeats and telomeric repeat artefacts to redo the base calls. A list of readnames corresponding to the candidate reads will then be produced

To apply this step, run the following command:
```
perl main.pl <input_fasta> <output_label>
```

### 3. basecall_problematic_reads
This directory contains a set of scripts to extract the fast5 for the required reads. Subsequent to that, we will basecall these reads using the tuned bonito model.

Note that fast5_subset from https://github.com/nanoporetech/ont_fast5_api is required to run this step. Also note that the extraction of the desired reads with fast5_subset can take some time.

To apply this step, run the following command:
```
perl main.pl <directory_of_all_fast5_files> <file_of_readnames_for_required_reads> <directory_to_output_extracted_fast5> <output_fastagz_file>
```

### 4. replace_problematic_reads
This step allows one to replace the reads with the basecalling repeat artefacts with the reads that have been fixed.

To apply this step, run the following command:
```
perl find_and_replace_fasta.pl <original_full_fasta_file> <rebasecalled_fasta_file>
```


## Tuning Nanopore bonito basecalling model for telomeric reads

There are a series of three steps if one would like to retrain the bonito basecalling model for telomeric regions.

### 1. generate_unmodified_data
This step generates the unmodified data needed for generating new training data. Also, the ground truth sequence is also extracted from the reference genome in this step.

To apply this step, run the following command:
```
perl main.pl <fast5_of_reads_from_telomeres> <original_bonito_model> <training_data_output_label> <reference_genome_fasta>
```

### 2. modify_training_data
This step modifies the training data (chunks.npy) using the ground truth data. This allows the modified training data to now be used by the bonito software to tune the original basecalling model.

To apply this step, run the following command:
```
perl main.pl <training_data_label>
```

### 3. train_model
This step retrains the bonito basecalling model with a low learning rate. A new basecalling model that is optimized for telomeric region will be generated after this step.

To apply this step, run the following command:
```
perl main.pl <original_basecalling_model> <training_data_directory> <tuned_model_name>
```


## FAQ
1. Why are strange non-telomeric sequences being generated after applying the pipeline?

Note that the tuned bonito basecalling model that we have published is compatible with Bonito v0.3.5, but not with higher versions of Bonito (e.g. v0.5 and above). Please check that your are using the correct version of Bonito. I will also be updating the tuned model for later versions of Bonito in the future after I have properly evaluated them.


## Contact
If you have further queries, please contact me at the following email address.

Kar-Tong Tan (ktan@broadinstitute.org)
