# Nextflow experimental pipeline for the analysis of the SARS-CoV-2 genome

This is experimental pipline used for the analysis SAR-CoV-2 genome, and exploration of the Nextflow framework.

# Dependencies

This pipeline requires the following software to be installed:

- python 3.9.x
- bwa
- samtools
- nproc

Also install packages from `requirements.txt`:

    pip install -r requirements.txt

# Installation

Download and install [Nextflow](https://www.nextflow.io/) with command: 

    wget -qO- https://get.nextflow.io | bash

put the resulting `nextflow` executable in your `$PATH`.

# Make a copy of run script and adjust run parameters

    cp run_nf_pipeline.sh.template run_nf_pipeline.sh

Edit flags in `run_nf_pipeline.sh` to adjust the run parameters.
- --reference_genome : path to the reference genome (fasta file)
- --reads : path to the reads (fastq file).
 
Reads MUST:
 - be a single path to two files.
 - be enclosed in single quotes (e.g. `'path/to/reads/*_{1,2}.fastq.gz'`).
 - be in the format `*_{1,2}.fastq.gz` (e.g. `sample1_1.fastq.gz` and `sample1_2.fastq.gz`).

# Run the pipeline locally 

    ./run_nf_pipeline.sh

# Results

Results will be in the `work` directory.
In addition, files `flowchart.png` and `report.html` will be generated in project root directory.
