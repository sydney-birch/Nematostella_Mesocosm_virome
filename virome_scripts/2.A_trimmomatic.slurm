#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00

##Load module
module purge
module load trimmomatic/0.39 

## Running slurm in 2_trimmomatic dir - run for each set of reads 
./2.A_trimmomatic.py -a ../raw_reads_2-19-24/21081FL-06-01-01_S1_L001_R1_001.fastq.gz+_ds.c38327b14c2c4b9e82c49d5a8328b70e/Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8 -b ../../../2_trimmomatic/
./2.A_trimmomatic.py -a ../raw_reads_2-19-24/21081FL-06-01-01_S1_L001_R1_001.fastq.gz+_ds.c38327b14c2c4b9e82c49d5a8328b70e/Fastq-021324-XP-fcB-22HGNTLT3-L002-I8I8 -b ../../../2_trimmomatic/
