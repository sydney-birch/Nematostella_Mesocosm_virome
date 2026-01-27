#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=headers
#SBATCH --nodes=1
##SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

./3_get_full_headers.sh uk_transcriptome.fa
