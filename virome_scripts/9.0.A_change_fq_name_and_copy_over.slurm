#!/bin/bash

#SBATCH --partition=Nebula
#SBATCH --job-name=name_change
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00


## Running name_change_script
./0.A_change_fq_name_and_copy_over.py -a name_change.txt -b /projects/areitze2_research/03_10_25_16S_QK_SB/
