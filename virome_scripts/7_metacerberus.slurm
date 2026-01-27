#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=metacerberus
#SBATCH --nodes=1
##SBATCH --mem=300GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules and activate mercat2
module load anaconda3
conda activate MetaCerberus-1.3.0
## source the slurm script to initialize the Ray worker nodes
source ray-slurm-metacerberus.sh
ray stop


## Run individual 14 assemblies:
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/NS-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out NS-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/ME-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out ME-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/NH-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out NH-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/FIELD-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out FIELD-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/MA-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out MA-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/SC-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out SC-T0_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/FL-T0_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out FL-T0_output

metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/NS-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out NS-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/ME-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out ME-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/NH-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out NH-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/FIELD-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out FIELD-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/MA-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out MA-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/SC-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out SC-T14_output
metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/FL-T14_assembly.fa --hmm "KOFam_all, COG, VOG, PHROG, CAZy" --dir_out FL-T14_output

