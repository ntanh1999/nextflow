# Pan-genome analysis nextflow pipeline
A nextflow pipeline for Bacterial pan-genome analysis 
## Installation
0. Download and install the appropriate conda, such as anaconda from 
   https://repo.anaconda.com/archive/
1. Create a conda environment with all the necessary dependencies: From the repository directory run
```bash

conda create -y -c conda-forge -c defaults --name nextflow python=3.7 mamba

source activate nextflow

mamba install -y -c conda-forge -c bioconda -c anaconda -c defaults  --file requirements.txt

```
2. Install nextflow
See https://www.nextflow.io/docs/latest/getstarted.html#installation 
## Pipeline Description
### Input
Input files could be:
- Genome assembly (*.fna).
- Illumina paired-end reads (*.fastq.gz).
- Illumina adapters sequences.
### Parameters
#### --reads
- The path of the reads FASTQ files.
- Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string value by single quote characters.
#### --assembly
- The path of the assembly FASTA files.
- Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string value by single quote characters.
#### --adapter
- The path of Illumina adapters (Default 'adapters/TruSeq3-PE-2.fa').
- Use absolute path.
#### --result
- The folder where the results will be stored for user.

## Examples
Activate conda environment
```
conda activate nextflow
```
Download example dataset (4 Klebsiella pneumoniae strain - 2 reads, 2 assembly)
```
cd data
bash download_Kp4.sh
```
Run the pipeline
```
cd ..
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --assembly 'data/*.fna' --result results
```
