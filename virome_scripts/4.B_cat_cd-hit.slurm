#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=RNAdb
#SBATCH --nodes=1
#SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Run Fastqc
module load cd-hit/4.8.1

## Running slurm in raw_reads dir

./4.B_cat_cd-hit.sh RNA_db rRNAdb
