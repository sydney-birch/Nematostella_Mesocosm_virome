#  Workflow for Nematostella mesocosom analysis

## Prep work: 
A) Transfer over all data and adjust names by running the change_raw_fq_file_names.py with name_change.txt --> has the conversion of names from admera (we sequenced 8 samples at UNCC intially - all other samples were sent to admera for sequencing)    

The admera fastq name Format: 21081FL-06-02-11_S27_L002_R1_001.fastq.gz     
Changed to: NH_T14-SC_B3_S27_L002_R1_001.fastq.gz    

`./0_change_raw_fq_file_names.py -a 0_name_change.txt -b Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8`

   
B) Download genome 
 `wget https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_932526225.1/`

 
	
## 1) Run fastqc
Use 1_fastqc.py to loop thru and run fastqc for each sample --> run in raw reads dir (1_fastqc.slurm)  
  `./1_fastqc.py -a ../../../1_fastqc/before_trim`    
  
      The actual line of code for fastqc: 
      fastqc {fastq_file} -o ../../../1_fastqc/before_trim
*output in fastqc dir --> before_trim dir*

 
	
## 2) Trimm reads and re-run fastqc 
  A) Run trimmomatic to trim the adaptors (run in 2_trimmomatic dir) - the script loops thru and runs for each sample (2.A_trimmommatic.slurm)

`./2.A_trimmomatic.py -a ../raw_reads_2-19-24/admera_seq_run/Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8 -b ../../../2_trimmomatic/`   

	The actual line of code for trimmomatic used: 
 
 	java -jar /apps/pkg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 4 {path_to_R1} {path_to_R2} -baseout {Sample_name}_filtered.fq.gz ILLUMINACLIP:/apps/pkg/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa:1:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36    
  
	Output - you get 4 files in this format in the 2_trimmomatic dir for each sample:
	    NH_T0-SC_B5_filtered_1P.fq.gz
	    NH_T0-SC_B5_filtered_1U.fq.gz
	    NH_T0-SC_B5_filtered_2P.fq.gz
	    NH_T0-SC_B5_filtered_2U.fq.gz
		#We are interested in the _1P and _2P files (paired files)

  *Note: 2.A was run on the samples sequenced from Admera, 2.A.2 was run on the 8 samples that were sequenced at UNCC due to differences in naming conventions*  `./2.A.2_trimmomatic.py -a ../raw_reads_10-3-23 -b ../2_trimmomatic`


   B) Re-Run Fastqc in the fastqc dir (2.B_fastqc.slurm) 
`./2.B_fastqc.py -a ../1_fastqc/after_trim -b ../2_trimmomatic`
*output goes to fastqc dir - after_trim dir*  

 
	
## 3) Map reads to genome to get two pots of data --> mapped and unmapped    

   ### A) Get unmapped reads (non Nematostella reads aka viral/microbial reads)

	Map reads (HISAT2 and Bowtie) and use sam tools to get mapped (1st alignment) and unmapped reads
		To get unmapped reads (All non-nematostella reads), Re-align the unmapped multiple times (4x total) to the genome - run 2 alignments with HISAT2 then 2 alignments with bowtie using high sensitivity settings

  
	#1st alignment - HISAT2 
		#1st Ali - (run slurm in 1st_ali_2-20-24): 
			./3.A_initial_hisat2_v2.py -a ../../2_trimmomatic -b 1 -c 2nd_ali -d ../3_HISAT2/1st_ali_2-20-24/
                         
			 The actual hisat2 code run in script: 
			 hisat2 -q -p 12 -x ../Nematostella_genome -1 {path_to_R1}_filtered_1P.fq.gz -2 {path_to_R2}_filtered_2P.fq.gz -S Nematostella_genome_ali_{alignment_#}_{Sample_Name}.sam
				## Output: sam file with _ali_1 in a dir in next ali dir (2nd_ali dir) - it makes the second ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_1_{1}.sam)		

  
		# run sam tools - make fastq files to be used in second ali of all reads that did NOT map to the genome (Run slurm in 2nd_ali dir)
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 1

                         The actual samtools code run in script:  
			 	samtools view Nematostella_genome_ali_{0}_*.sam -f 0x4 -h -b -o unalign_reads.bam
				samtools collate -u -O unalign_reads.bam | samtools fastq -1 unaligned_{0}_{1}_ali_paired1.fq -2 unaligned_{0}_{1}_ali_paired2.fq -0 /dev/null -s /dev/null -n
	
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_1_ali_paired1.fq)
						


	#2nd alignment - HISAT2 
		#2nd Ali - Run slurm in 2nd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.C_hisat2_2_ali.slurm
			./3.C_hisat2_next_ali_v2.py -b 2 -c 3rd_ali_BT 
		        	## Output: sam file with _ali_2 in a subdir in next ali dir (3rd_ali dir) - it makes the third ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_2_{1}.sam)

	    
		# run sam tools - make fastq files to be used in third ali (Run slurm in 3rd_ali dir) 
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 2
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_2_ali_paired1.fq)     	
		        	
		        	
		        	
	#3rd alignment - bowtie
		#First - copy over genome file ending in .fna (in 1st_ali dir) - need to make a bowtie index
			
		#3rd Ali - Run slurm in 3rd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.D_bowtie_ali_3.slurm
			./3.D_bowtie_ali_initial.py -b 3 -c 4th_ali_BT

      			The actual code run in script: 
	 			bowtie2 -q -p 12 --very-sensitive-local -x ../Nematostella_genome -1 {sample_name}_paired1.fq -2 {sample_name}_paired2.fq -S Nematostella_genome_BT_ali_{0}_{1}.sam".format
    
				# Output: sam file with BT_ali_3 in a dir in next ali dir (4th_ali dir) - it makes the fourth ali dir and all sub dirs which contain sam files (Nematostella_genome_BT_ali_3_{1}.sam)

   
		# run sam tools - make fastq files to be used in fourth ali (Run slurm in 4th_ali dir)	
			sbatch 3.E_samtools_processing_BT.slurm
			./3.E_samtools_processing_BT.py -b 3
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_3_ali_paired1.fq)
		
	
	#4th alignment - bowtie (hopefully final alignment - if it is rename 5th_ali to final_sam_files)
		#4th Ali - Run slurm in 4th_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.F_bowtie_4_ali.slurm
			./3.F_bowtie_next_ali.py -b 4 -c 5th_ali
				## Output: sam file with BT_ali_4 in a dir in next ali dir (5th_ali dir) - it makes the fifth ali dir and all sub dirs which contain sam files (Nematostella_genome_BT_ali_4_{1}.sam)
    
			####  if it is rename 5th_ali to final_unmapped_sam_files #### 
							
					
		# run sam tools - make fastq files  (Run slurm in 5th_ali dir/ final_sam_files dir)	
			sbatch 3.E_samtools_processing_BT.slurm
			./3.E_samtools_processing_BT.py -b 4
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_4_ali_paired1.fq)   


 ### B) Get mapped reads (Nematostella reads) and move unmapped reads into next step dir (RNA filtration)    
First, make 3 dirs:    
     `mkdir /Nematostella_transcriptomics/4_RNA_filt`    
     `mkdir /Nematostella_transcriptomics/4_RNA_filt/mapped_fastq_files`     
     `mkdir /Nematostella_transcriptomics/4_RNA_filt/unmapped_fastq_files`


Next, Run samtools to make fastq files for each sample of the MAPPED reads to the genome from the first alignment (in 2nd_ali dir)     
        `sbatch 3.H_samtools_mapped_reads.slurm`    
	`./3.H_samtools_mapped_reads.py -b 1 `    
	    `#Output: 2 fastq files (aligned_{sample_name}_1_ali_paired1.fastq) placed in 4_RNA_filt/mapped_fastq_files`


Move Unmapped reads into RNA filtration main dir        	  
        `sbatch 3.G_mv_fastqs.slurm`      
	`./3.G_mv_fastqs_v2.py -b 4 -c 4_RNA_filt `     
		`#Output: all unmapped fastqs will be moved to the RNA_filt unmapped_fastq_files dir`

  ## 4) rRNA Filtration of Viral/Microbial Reads
  Ultimately, in this step we want to have presumably all viral reads, so we need to filter out any prokaryotic and eukaryotic reads.    

  
  A) Download files to create your rRNA database     
  `./4.A_RNAdb_download_files.sh`    

  
  B) Concatenate all the downloads and run cd-hit to remove duplicate sequences.     
  `sbatch 4.B_cat_cd-hit.slurm` this runs: `./4.B_cat_cd-hit.sh RNA_db rRNAdb`     

  
  C) Map reads to rRNA db using hisat2       
  `sbatch 4.C_filtered_hisat2.slurm` this runs: `./4.C_filtered_hisat2.py -a filtered_unmapped_fastq_files -b rRNAdb_ed1.fa -c unmapped_fastq_files`     

  
  D) Use samtools to make fastq files of the reads that did *not* map to rRNA db - these are the primarily viral reads      
  `sbatch 4.D_samtools_processing.slurm` this runs: `./4.D_samtools_processing.py -a filtered_unmapped_fastq_files`      

  
  E) As an additional filtration step, we also mapped the reads to the NCBI nr database with viruses removed     
  `sbatch 4.E_nr_filtered_hisat2.slurm` this runs: `./4.E_nr_filtered_hisat2.py -a ncbi_nr_filtered_unmapped_fastq_files -b ncbi_nr_DB.fa -c filtered_unmapped_fastq_files`      

  
  F) Use samtools to make fastq files of the reads that did *not* map to nr database       
  `sbatch 4.F_samtools_processing.slurm ` this runs: `./4.F_samtools_processing.py -a ncbi_nr_filtered_unmapped_fastq_files`     
  
  
 ## 5) Make Viral Assemblies 

 Here, we are going to collapse the replicates to make assemblies for each sample. So three replicates will go into the generation of each location/timepoint sample resulting in a total of 14 assemblies (i.e., 7 locations and 2 timepoints)

 Use spades --rnaviral to generate the viral assemblies, an example line of code:     
 ```spades.py --rnaviral --pe1-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B1_4_ali_paired1.fq.gz --pe2-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B2_4_ali_paired1.fq.gz --pe3-1 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B4_4_ali_paired1.fq.gz --pe1-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B1_4_ali_paired2.fq.gz --pe2-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B2_4_ali_paired2.fq.gz --pe3-2 ../../4_RNA_filt/filtered_unmapped_fastq_files/filtered_unaligned_NH_T0-NS_B4_4_ali_paired2.fq.gz -o NS_T0_output```     

 All assembly code can be found in `5.A_spades_indiv_filtered_unmaapped.slurm`     

 
 ## 6) Viral Taxonomy Analysis     
To investigate the viral taxa found in our analysis, we ultimately used NCBI taxon kit. But we first had to identify protiens and find their taxon IDs so we used BLAST and BioEntrez.     


A) Download the refseq viral database.    
`wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz `     
*This has 683,238 viral protiens*      


B) Make blast databases for each of the 14 viral assemblies     
`./6.0_blastdb.sh`     


C) Run a tBLASTn using the RefSeq viral database against each of the 14 viral asssemblies.
`sbatch 6.1_blast.slurm`   
This runs a blast search --> `tblastn -query viral.1.protein.faa -db ./blastdb/FIELD-T0_assembly.fa_FIX.fa -out FIELD-T0_ref_blastout -outfmt 6 -max_target_seqs 1 -evalue 1e-5 -num_threads 12 -best_hit_score_edge 0.25 -best_hit_overhang 0.1`        


D) Decided to trim the blastout to 70 percent identity - more managable and more stringent:    
`./6.1.B_run_trim_blastout.sh`


E) Get Accession IDs from the 70_pi_blastout table
`./6.2_get_accessions.sh`


F) Get the full headers from the accession IDs to run with selectSeqs   
`./6.3_get_full_headers.sh`   this runs --> 6.3.B_get_full_headers.py    

G) Use BioEntrez to get taxids    
   * Key variables to adjust:
      * dbfrom="protein", id=nuccoreid , linkname="protein_taxonomy"
      * Generate an API key on NCBI in your account settings
      * Run slurm  6.C_run_bioentrez.slurm in chunks if needed
```
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/FIELD-T0_ref_blastout_70_pi_hits -o FIELD-T0 	#done Round 1 = 1103  (total 1103)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/FIELD-T14_ref_blastout_70_pi_hits -o FIELD-T14 	#done Round 1 = 559 (total 559)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/FL-T0_ref_blastout_70_pi_hits -o FL-T0 			#done Round 1 = 461 (total 461/456)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/FL-T14_ref_blastout_70_pi_hits -o FL-T14		#done Round 1 = 1077  (total 1057/1077)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/MA-T0_ref_blastout_70_pi_hits -o MA-T0 			#done round 1 = 475 (total 475)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/MA-T14_ref_blastout_70_pi_hits -o MA-T14 		#done Round 1 = 809, Round 2 = 230  (total 1024/1039)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/ME-T0_ref_blastout_70_pi_hits -o ME-T0 			#done Round 1 = 654 (total 536/654)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/ME-T14_ref_blastout_70_pi_hits -o ME-T14 		#done Round 1 = 811 (total 784/811)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/NH-T0_ref_blastout_70_pi_hits -o NH-T0 			#done Round 1 = 744 (total 684/744)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/NH-T14_ref_blastout_70_pi_hits -o NH-T14 		#done Round 1 = 925 (total 884/925)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/NS-T0_ref_blastout_70_pi_hits -o NS-T0 			#done Round 1 = 155, Round 2 = 125, Round 3 = 100 (total 380)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/NS-T14_ref_blastout_70_pi_hits -o NS-T14 		#done Round 1 = 1112  (total 1051/1112)
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/SC-T0_ref_blastout_70_pi_hits -o SC-T0 			#done Round 1 = 254, Round 2 = 61 (total 315) 
./6.B_get_taxid-prot.py -i ../70_pi_hit1_accessions/SC-T14_ref_blastout_70_pi_hits -o SC-T14 		#done Round 1 = 742, Round 2 = 489, Round 3 = 111  (total 1312/1342)
```
   * Concatenate files that need it - if needed to run on small chuncks:
   ```
     #MA-T14 	
		cat MA-T14_accid-taxid.txt >> MA-T14_TOTAL_accid-taxid.txt
		cat MA-T14-2_accid-taxid.txt >> MA-T14_TOTAL_accid-taxid.txt #1039
		
		cat MA-T14_taxid.txt >> MA-T14_TOTAL_taxid.txt	
		cat MA-T14-2_taxid.txt >> MA-T14_TOTAL_taxid.txt #1039
		
	#NS-T0
		cat NS-T0_accid-taxid.txt >> NS-T0_TOTAL_accid-taxid.txt
		cat NS-T0-2_accid-taxid.txt >> NS-T0_TOTAL_accid-taxid.txt
		cat NS-T0-3_accid-taxid.txt >> NS-T0_TOTAL_accid-taxid.txt #380
		
		cat NS-T0_taxid.txt >> NS-T0_TOTAL_taxid.txt	
		cat NS-T0-2_taxid.txt >> NS-T0_TOTAL_taxid.txt
		cat NS-T0-3_taxid.txt >> NS-T0_TOTAL_taxid.txt #380
	
	#SC-T0
		cat SC-T0_accid-taxid.txt >> SC-T0_TOTAL_accid-taxid.txt
		cat SC-T0-2_accid-taxid.txt >> SC-T0_TOTAL_accid-taxid.txt #315
		
		cat SC-T0_taxid.txt >> SC-T0_TOTAL_taxid.txt	
		cat SC-T0-2_taxid.txt >> SC-T0_TOTAL_taxid.txt #315
			
	#SC-T14
		cat SC-T14_accid-taxid.txt >> SC-T14_TOTAL_accid-taxid.txt
		cat SC-T14-2_accid-taxid.txt  >> SC-T14_TOTAL_accid-taxid.txt
		cat SC-T14-3_accid-taxid.txt  >> SC-T14_TOTAL_accid-taxid.txt #1342
		
		cat SC-T14_taxid.txt >> SC-T14_TOTAL_taxid.txt	
		cat SC-T14-2_taxid.txt >> SC-T14_TOTAL_taxid.txt
		cat SC-T14-3_taxid.txt >> SC-T14_TOTAL_taxid.txt #1342
   ```

H) Run Taxon Kit from NCBI to get viral taxonomic information using the viral taxid 

```
cat final_taxids/FIELD-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FIELD-T0_linage.txt

cat final_taxids/FIELD-T14_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FIELD-T14_linage.txt
    
cat final_taxids/FL-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FL-T0_linage.txt

cat final_taxids/FL-T14_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> FL-T14_linage.txt

cat final_taxids/MA-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> MA-T0_linage.txt		

cat final_taxids/MA-T14_TOTAL_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> MA-T14_linage.txt

cat final_taxids/ME-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> ME-T0_linage.txt
    
cat final_taxids/ME-T14_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> ME-T14_linage.txt

cat final_taxids/NH-T0_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> NH-T0_linage.txt		

cat final_taxids/NH-T14_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> NH-T14_linage.txt

cat final_taxids/NS-T0_TOTAL_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> NS-T0_linage.txt
    
cat final_taxids/NS-T14_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> NS-T14_linage.txt		

cat final_taxids/SC-T0_TOTAL_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> SC-T0_linage.txt

cat final_taxids/SC-T14_TOTAL_taxid.txt \
    | ./taxonkit reformat --data-dir TAXONKIT_DB -I 1 -F -P -f "{k}\t{p}\t{c}\t{o}\t{f}\t{g}\t{s}\t{t}" >> SC-T14_linage.txt
```

I) Analyze Taxonomy data in R 
   * First you'll need to process the data in excel - make a total_Lineage file - import each lineage file and adjust data to columns
   	    * Each tab in this file is a location
   * Additionally, make a spreadsheet called genus.csv and copy each location into a column in this sheet - this will be used as the input in R script
   * Make a third spreadsheet of the genetic composition of the viruses present using the VMR spreadsheet using xlookup
        * Run two xlookups - one using the species column, the other using virus name column to conduct your search
        * Then run an if statement to merge the data: `=IF(V2=0,W2,V2)` where V2 = species col and W2 = virus name col
        * remove s__ by find and replace


Now run R scripts to look at taxanomic overlaps and run diversity statistics: genus_work.R


 ## 7) Viral Functional Analysis       
   * Run MetaCerberus on the 14 viral assemblies generated in step 5 (spades assemblies)   

      * sbatch 7_metacerberus.slurm    
      * example line: `metacerberus.py --prodigal ../5_assemblies/filtered_unmapped_assembly/NS-T0_assembly.fa --hmm All --dir_out updated_NS-T0_output`   

   * Copy over output files to your computer   
   * Generate spreadsheets to input into R for KEGG, VOG, PHROG, COG, FOAM (or whichever you are most interested in)
        * Column A should have the term name, each column after should be a location/time point filled in the the corresponding counts
        * Once you make this spreadsheet - transpose it and make a csv with three columns: Term Name (Function), Count, Location
   * Run R script to analyze MetaCerberus data: Metacerberus_stats_output.R

 

	
## 9) 16s rRNA Analysis (QIIME2)   
   * For the QIIME2 analysis - I'm following the Moving Pictures tutorial (https://docs.qiime2.org/2020.6/tutorials/moving-pictures/) and the Atacama soil tutorial (https://docs.qiime2.org/2024.10/tutorials/atacama-soils/)

Prep work: 
   * Current data is in a shared projects folder - copy over the files into your working dir and change names   
`./9.0.A_change_fq_name_and_copy_over.py -a name_change.txt -b /projects/areitze2_research/03_10_25_16S_QK_SB/`   

   * Remove underscores in name
`./9.0.B_remove_underscores_in_name.py -b ../raw_reads/`   
   * Make a metadata file and copy to terminal
     `scp Meso_2022_16s_metadata_v2.txt sbirch1@hpc.charlotte.edu:../../scratch/sbirch1/QIIME_Meso_2022`
     
Step 1: Import your data by creating a qiime artifact - follow Casava 1.8 paired-end demultiplexed fastq   

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_reads \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza
```
Visualize output (go to https://view.qiime2.org): 
```
qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux.qzv
```


Step 2: Sequence quality control and feature table construction - I'm using option 1: DADA2

```
#part a - quality filter and truncate(denoise) 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 300 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats denoising-stats-dada2.qza

#part b - Make visualizations and summaries 
qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv \
  --m-sample-metadata-file metadata-file-2022_v2.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats-dada2.qza \
  --o-visualization denoising-stats-dada2.qzv

 	#send to computer to visualize

	#rename files for next section  
		mv rep-seqs-dada2.qza rep-seqs.qza
		mv table-dada2.qza table.qza
```

Step 3: Generate a tree for phylogenetic diversity analysis   
```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza



## need to chose a sampling depth based on table.qzv
	#Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.
 		#use the table.qzv --> intaractive sample tab --> change to loc-time --> adjust sampling depth
 	
 		##Doing differnt sampling depths:
 				53024 --> Retained 3,764,704 (44.99%) features in 71 (98.61%) samples at the specifed sampling depth. (all samples retained - its one of the blanks missing)
 			*2  56322 --> Retained 3,942,540 (47.12%) features in 70 (97.22%) samples at the specifed sampling depth. (leaves 4 reps Field T14)
 			*1  61745 --> Retained 4,260,405 (50.92%) features in 69 (95.83%) samples at the specifed sampling depth. (leaves 3 reps Field T14)
 			*3  68009 --> Retained 4,556,603 (54.46%) features in 67 (93.06%) samples at the specifed sampling depth. (leaves 2 reps Field T14 and 4 reps NH T0)
			    75304 --> Retained 4,744,152 (56.70%) features in 63 (87.50%) samples at the specifed sampling depth.	(leaves 1 rep Field T14 and 4 reps NH T0, FL T0, MA T14, NS T14)
```

**Chose SAMPLING DEPTH = 61745**   

Step 4: Alpha and beta diversity analysis   

Choosing sampling depth --> 61745 --> Retained 4,260,405 (50.92%) features in 69 (95.83%) samples at the specifed sampling depth. (leaves 3 reps Field T14)   

```	
#core metrics phylogenetic analysis
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 61745 \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --output-dir core-metrics-results-61745 

#download visualizations


## Test for associations between categorical metadata columns and alpha diversity data 


	### Alpha diversity 

		#A: Faith phylogenetic Diveristy - measure of community richness (Measures of biodiversity that incorporates phylogenetic difference between species)
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-61745/faith_pd_vector.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --o-visualization core-metrics-results-61745/faith-pd-group-significance.qzv

		#B: Evenness metrics 
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-61745/evenness_vector.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --o-visualization core-metrics-results-61745/evenness-group-significance.qzv

		#C: Shannon index - Calculates richness and diversity using a natural logarithm (Accounts for both abundance and evenness of the taxa present)
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-61745/shannon_vector.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --o-visualization core-metrics-results-61745/shannon-group-significance.qzv  

#download visualizations 

 
 
 	### Beta Diversity 
 	
		#A: Testing beta diversity on specific columns using unweighted unifrac 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-61745/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --m-metadata-column loc-time \
  --o-visualization core-metrics-results-61745/unweighted-unifrac-loc-time-significance.qzv \
  --p-pairwise 
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-61745/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --m-metadata-column Location \
  --o-visualization core-metrics-results-61745/unweighted-unifrac-Location-significance.qzv \
  --p-pairwise 
  
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-61745/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --m-metadata-column Timepoint \
  --o-visualization core-metrics-results-61745/unweighted-unifrac-Timepoint-significance.qzv \
  --p-pairwise 
 
 
		#B:  Generating PCoA plots with emperor: unweighted unifrac 
qiime emperor plot \
  --i-pcoa core-metrics-results-61745/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --p-custom-axes loc-time \
  --o-visualization core-metrics-results-61745/unweighted-unifrac-emperor-loc-time.qzv

		#C: Generating PCoA plots with emperor: bray curtis 
qiime emperor plot \
  --i-pcoa core-metrics-results-61745/bray_curtis_pcoa_results.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --p-custom-axes loc-time \
  --o-visualization core-metrics-results-61745/bray-curtis-emperor-date.qzv

#download visualizations



## Alpha rarefraction 

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 90000 \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --o-visualization alpha-rarefaction.qzv

#download visualizations
```

Step 5: Taxonomic analysis

First download taxonomic file 
```
# Choosing: 
http://ftp.microbio.me/greengenes_release/2022.10/  
	# 2022.10.backbone.full-length.nb.qza
    Naive Bayes classifier trained on the V4 region, and separately, the full length 16S from the backbone sequences.

# use wget to download:
wget \
  -O "2022.10.backbone.full-length.nb.qza" \
  --no-check-certificate \
"https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.full-length.nb.qza"
```
			
Run analysis - assign taxonomy to the sequences 
```
qiime feature-classifier classify-sklearn \
  --i-classifier 2022.10.backbone.full-length.nb.qza\
  --i-reads rep-seqs.qza \
  --o-classification taxonomy_v2.qza

qiime metadata tabulate \
  --m-input-file taxonomy_v2.qza \
  --o-visualization taxonomy_v2.qzv

#download visualizations


#view the taxonomic composition of our samples with interactive bar plots
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy_v2.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --o-visualization taxa-bar-plots_v2.qzv
  
#download visualizations
```

Step 6: Differential abundance testing with ANCOM - run at specific taxonomic levels - only showing species run here 

```
### Level 7 - Species 
  Dont do feature table filter
  
qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy_v2.qza \
  --p-level 7 \
  --o-collapsed-table level-7-table.qza
  
qiime composition ancombc \
  --i-table level-7-table.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --p-formula 'loctime' \
  --o-differentials 7-ancom-loc-time.qza  
  
qiime composition da-barplot \
  --i-data 7-ancom-loc-time.qza  \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization level-7-da-barplot-loc-time.qzv 
  
#download visualizations 


#FIX COMPARISON: from FLT0 to FieldT0
qiime composition ancombc \
  --i-table level-7-table.qza \
  --m-metadata-file metadata-file-2022_v2.tsv \
  --p-formula 'loctime' \
  --p-reference-levels loctime::Field-T0 \
  --o-differentials 7-ancom-loc-time_fieldt0.qza  
  
qiime composition da-barplot \
  --i-data 7-ancom-loc-time_fieldt0.qza  \
  --p-significance-threshold 0.001 \
  --p-level-delimiter ';' \
  --o-visualization level-7-da-barplot-loc-time_fieldt0.qzv 
  
#download visualizations
```

Make Visualizations in R: qiime2_meso_2022.R


