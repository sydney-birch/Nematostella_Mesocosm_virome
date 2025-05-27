#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=count_reads
#SBATCH --nodes=1
#SBATCH --mem=400GB
##SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules

## Run count reads script 
## mapped
./count_reads_in_fastqs.py -a mapped_fastq_files

## unmapped
./count_reads_in_fastqs.py -a unmapped_fastq_files

## filtered
./count_reads_in_fastqs.py -a filtered_unmapped_fastq_files

## nr ncbi filtered
./count_reads_in_fastqs.py -a ncbi_nr_filtered_unmapped_fastq_files
