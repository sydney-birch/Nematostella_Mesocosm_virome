# Host analysis workflow

We are picking up on Step 8 of the overall analysis to analyze the host transcriptomics

## 8) Host (Nematostella) Analysis  
### Prep) Use Salmon to align reads and get counts

First, align the mapped (Nematostella) reads to the UK transcriptome to get counts and sequence IDs for the reads. The Host reads were identified in Step 3.B of the virome analysis (reads that mapped to the Nematostella transcriptome). 

   A) Copy over transcriptome: GCF_932526225.1/rna.fna (I renamed to uk_transcriptome.fa)

   B) submit Salmon slurm: sbatch 8_salmon-mapped.slurm
```
#create the index (only run 1 time)
salmon index -t uk_transcriptome.fa -i mapped_index -k 31

#run salmon script to align each replicate
./8_run_salmon.py -a ../4_RNA_filt/mapped_fastq_files -b ../../8_salmon -c aligned

#example line of code of salmon alignment line: 
salmon quant -p 12 --seqBias --gcBias -i mapped_index -l A -1 {0}/{1}_paired1.fastq.gz -2 {0}/{1}_paired2.fastq.gz -o {1}".format(fastq_dir,sample)
```

### A) Differential Gene Expression     
A.1) copy over Salmon output files to computer for R analyses - place in a folder called: mapping

`scp -r sbirch1@hpc.charlotte.edu:./../../scratch/sbirch1/Nematostella_transcriptomics/8_salmon/aligned_NH_T* ./`

A.2) Run edgeR script in R to get Differentially Expressed Genes between Timepoints (T0 vs T14) for each location - write out into text files 

    script: 8.A.2_edgeR_Host.R   

A.3) Input DEG lists back to the terminal and get header information for each DEG to be used with selectSeqs 

```
#first run this script to turn the edgeR input into a list of headers to retrieve: 
./8.A.3_run_get_headers_from_edgeR.sh --> This runs 8.A.3_run_get_headers_from_edgeR.py 

#Next, move all header lists into a new dir called mapped_selectseqs_header_lists: 
mv *-header_list ../mapped_selectseqs_header_lists/

#Now use script to get full headers from initial accession id to run with selectSeqs: 
./8.A.3b_get_full_headers.sh uk_transcriptome.fa --> this will run: script 8.A.3b_get_full_headers.py
```

A.4) Run selectSeqs to extract the sequences for each DEG from the transcriptome using the full headers:
```
#Run select seq bash script: 
./8.A.4_run_selectseqs.sh --> this will run ./selectSeqs.pl
```

A.5) Run DEGs on Metacerberus to get potential functional annotations 
```
#Field: 
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/FIELD_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_FIELD_T0vT14_DEGs

#FL
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/FL_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_FL_T0_v_T14_DEGs        

#MA
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/MA_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_MA_T0_v_T14_DEGs  

#ME
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/ME_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_ME_T0_v_T14_DEGs        

#NH
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/NH_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_NH_T0_v_T14_DEGs

#NS
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/NS_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_NS_T0_v_T14_DEGs

#SC
metacerberus.py --fraggenescan ../8_salmon/mapped_selectseqs_fasta_files/SC_T0vT14_DEGs_edgeR_output.txt-header_list_full_header-seqs.fa --hmm KOFam_eukaryote --dir_out mapped_SC_T0_v_T14_DEGs
````

A.6) copy output to computer - turn counts into spreadsheets and make figures in R using script: 8.A.6_Metacerberus_Host_DEGs.R

A.7) Make DEG heatmaps 
```
#run script to turn list of DEGs into a string to import to R: 
./8.A.7_prep_for_R_string_deg.py
```

   Make heatmaps in R using: 8.A.7_Host_DEG_Heatmaps.R    
   
A.8) Use the Host data to make a heatmap of expression of the Nematostella 56 Immune gene set:   
    8.A.8_Host_Immune_gene_set_heatmaps.R
 
 ### B) Host WGCNA (Weighted Gene Co-expression Network Analysis) 
B.1) For the WGCNA - Run the WGCNA R script (8.B.1_WGCNA)
   * For this you will need to access the Salmon Host data in the mapping directory on your computer.  
   * You will also need to generate a text file called library_to_stages.txt that details which samples belong to what sample group.   
   * Your output will be boxplots of Module expression data and a heatmap of particular modules along with spreadsheets of module gene information.  

B.2) Once you have exported a list of your genes from the modules of interest (I'm using top 5 most significant modules) - Run metacerberus to get functional information:  
   * Use excel to get the full header information with the UK transcriptome header excel file (use xlookup)  
   * Then, in the terminal, run selectSeqs.pl to pull out the sequences for each header in each module you are interested in.
     ```
    ./selectSeqs.pl -f module_17_headers ./uk_transcriptome.fa >> module_17_seqs.fa
	./selectSeqs.pl -f module_21_headers ./uk_transcriptome.fa >> module_21_seqs.fa
	./selectSeqs.pl -f module_33_headers ./uk_transcriptome.fa >> module_33_seqs.fa
	./selectSeqs.pl -f module_60_headers ./uk_transcriptome.fa >> module_60_seqs.fa
	./selectSeqs.pl -f module_96_headers ./uk_transcriptome.fa >> module_96_seqs.fa
	```
   * Once you have your fasta file(s), run metacerberus on those files.
     `Run Metacerberus slurm: sbatch 8.B.2_metacerb_MEs.slurm`
   * Export, metacerberus data, and make figures with R - use 8.B.2_Metacerberus_ME_stats_output.R script
     
