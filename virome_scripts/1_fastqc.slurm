#!/bin/bash

#SBATCH --partition=Nebula
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00

##Run Fastqc
module load fastqc

## Running slurm in raw_reads dir

./1_fastqc.py -a ../../../1_fastqc/before_trim
