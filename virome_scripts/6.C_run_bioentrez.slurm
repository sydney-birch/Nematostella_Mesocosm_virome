#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=bioentrez
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules 
#module load blast
module load anaconda3/2022.10

## Run blast
##./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/SC-T0_ref_blastout_70_pi_hits -o SC-T0
##./6.B_get_taxid-prot_v2.py -i SC-T0_round-2 -o SC-T0-2

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/NS-T0_ref_blastout_70_pi_hits -o NS-T0
##./6.B_get_taxid-prot_v2.py -i NS-T0-Round_2 -o NS-T0-2
##./6.B_get_taxid-prot_v2.py -i NS-T0-Round_3 -o NS-T0-3

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/FL-T0_ref_blastout_70_pi_hits -o FL-T0

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/MA-T0_ref_blastout_70_pi_hits -o MA-T0

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/ME-T0_ref_blastout_70_pi_hits -o ME-T0

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/FIELD-T14_ref_blastout_70_pi_hits -o FIELD-T14

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/NH-T0_ref_blastout_70_pi_hits -o NH-T0

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/ME-T14_ref_blastout_70_pi_hits -o ME-T14

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/NH-T14_ref_blastout_70_pi_hits -o NH-T14

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/MA-T14_ref_blastout_70_pi_hits -o MA-T14
##./6.B_get_taxid-prot_v2.py -i MA-T14-Round_2 -o MA-T14-2

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/NS-T14_ref_blastout_70_pi_hits -o NS-T14

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/FL-T14_ref_blastout_70_pi_hits -o FL-T14

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/FIELD-T0_ref_blastout_70_pi_hits -o FIELD-T0

##./6.B_get_taxid-prot_v2.py -i ../70_pi_hit1_accessions/SC-T14_ref_blastout_70_pi_hits -o SC-T14
##./6.B_get_taxid-prot_v2.py -i SC-T14-Round_2 -o SC-T14-2
./6.B_get_taxid-prot_v2.py -i SC-T14-Round_3 -o SC-T14-3
