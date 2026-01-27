#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=blastRest
#SBATCH --nodes=1
#SBATCH --mem=64GB
#SBATCH --ntasks-per-node=12
#SBATCH --time=96:00:00

##Load modules 
module load blast

#Run Blast
##tblastn -query viral.1.protein.faa -db ./blastdb/FL-T0_assembly.fa_FIX.fa -out FL-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001		
##tblastn -query viral.1.protein.faa -db ./blastdb/FL-T14_assembly.fa_FIX.fa -out FL-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 0.00001 -num_threads 12

##tblastn -query viral.1.protein.faa -db ./blastdb/FL-T0_assembly.fa_FIX.fa -out FL-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1
##tblastn -query viral.1.protein.faa -db ./blastdb/FL-T14_assembly.fa_FIX.fa -out FL-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

##tblastn -query viral.1.protein.faa -db ./blastdb/MA-T0_assembly.fa_FIX.fa -out MA-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
##tblastn -query viral.1.protein.faa -db ./blastdb/MA-T14_assembly.fa_FIX.fa -out MA-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

##tblastn -query viral.1.protein.faa -db ./blastdb/ME-T0_assembly.fa_FIX.fa -out ME-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
##tblastn -query viral.1.protein.faa -db ./blastdb/ME-T14_assembly.fa_FIX.fa -out ME-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

tblastn -query viral.1.protein.faa -db ./blastdb/NH-T0_assembly.fa_FIX.fa -out NH-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
tblastn -query viral.1.protein.faa -db ./blastdb/NH-T14_assembly.fa_FIX.fa -out NH-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

tblastn -query viral.1.protein.faa -db ./blastdb/NS-T0_assembly.fa_FIX.fa -out NS-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
tblastn -query viral.1.protein.faa -db ./blastdb/NS-T14_assembly.fa_FIX.fa -out NS-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

tblastn -query viral.1.protein.faa -db ./blastdb/SC-T0_assembly.fa_FIX.fa -out SC-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
tblastn -query viral.1.protein.faa -db ./blastdb/SC-T14_assembly.fa_FIX.fa -out SC-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		

tblastn -query viral.1.protein.faa -db ./blastdb/FIELD-T0_assembly.fa_FIX.fa -out FIELD-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		
tblastn -query viral.1.protein.faa -db ./blastdb/FIELD-T14_assembly.fa_FIX.fa -out FIELD-T14_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1		




#move to blastout dir
##mv ./blastdb/*_blastout ./blastout/
