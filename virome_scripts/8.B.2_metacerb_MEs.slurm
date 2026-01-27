#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=metacerberus
#SBATCH --nodes=1
##SBATCH --mem=400GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00

##Load modules and activate mercat2
module load anaconda3
conda activate metacerberus

## source the slurm script to initialize the Ray worker nodes
##source ray-slurm-metacerberus.sh
	
##metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/NH_T0vFIELD_T14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_NH_T0_v_CTRL_DEGs
metacerberus.py --fraggenescan module_17_seqs.fa --hmm KOFam_eukaryote --dir_out metacerb_module_17
metacerberus.py --fraggenescan module_21_seqs.fa --hmm KOFam_eukaryote --dir_out metacerb_module_21
metacerberus.py --fraggenescan module_33_seqs.fa --hmm KOFam_eukaryote --dir_out metacerb_module_33
metacerberus.py --fraggenescan module_60_seqs.fa --hmm KOFam_eukaryote --dir_out metacerb_module_60
metacerberus.py --fraggenescan module_96_seqs.fa --hmm KOFam_eukaryote --dir_out metacerb_module_96
