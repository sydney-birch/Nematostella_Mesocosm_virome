## Overview of nematostella mesocosom transcriptome pipline

# Prep work: 
A) Transfer over all data and adjust names by running the change_raw_fq_file_names.py with name_change.txt --> has the conversion of names from admera
'./change_raw_fq_file_names.py -a name_change.txt -b Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8'
			
B) Download genome 
 'wget https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_932526225.1/'
	
	
### 1) Run fastqc
Use 1_fastqc.py to loop thru and run fastqc for each sample --> run in raw reads dir   
  './1_fastqc.py -a ../../../1_fastqc/before_trim'
*output in fastqc dir --> before_trim dir*
	
	
### 2) Trimm reads and re-run fastqc 
  A) Run trimmomatic slurm and .py to trim the adaptors (run in 2_trimmomatic dir) - the script loops thru and runs for each sample

'./2.A_trimmomatic.py -a ../raw_reads_2-19-24/admera_seq_run/Fastq-021324-XP-fcB-22HGNTLT3-L001-I8I8 -b ../../../2_trimmomatic/'

	*Output - you get 4 files in this format in the 2_trimmomatic dir for each sample:*
	NH_T0-SC_B5_filtered_1P.fq.gz
	NH_T0-SC_B5_filtered_1U.fq.gz
	NH_T0-SC_B5_filtered_2P.fq.gz
	NH_T0-SC_B5_filtered_2U.fq.gz
	*We are interested in the _1P and _2P files (paired files)*

   B) Re-Run Fastqc in the fastqc dir (2.B_fastqc.slurm) 
      './2.B_fastqc.py -a ../1_fastqc/after_trim -b ../2_trimmomatic'
       *output goes to fastqc dir - after_trim dir*  
	
	
### 3) Map reads to genome to get two pots of data --> mapped and unmapped
	 Run HISAT2 and use sam tools to get mapped (1st alignment) and unmapped reads
		 Re-align the unmapped - run 2 alignments with HISAT2 then 2 alignments with bowtie
		
		#3.A_intial_hisat2.slurm --> run in 1st_ali (takes about 1.5 days to run)
		#3.A_initial_hisat2_v2.py
		
		#3.B_samtools_processing.slurm --> run in 2nd_ali, 3rd_ali (takes about 5 hours to run)
		#3.B_samtools_processing.py
		
		#3.C_hisat2_#_ali.slurm --> run in 2nd_ali
		#3.C_hisat2_next_ali_v2.py
		
		#3.D_bowtie_intial.slurm --> run in 3rd_ali
		#3.D_bowtie_ali_inital.py
		
		#3.E_samtools_processing_BT.slurm --> run in 4th_ali, 5th_ali(rename to final_sam_files)
		#3.E_samtools_processing_BT.py
		
		#3.F_bowtie_#_ali.slurm --> 4th_ali
		#3.F_bowtie_next_ali.py
		
		#3.G_mv_fastqs.slurm
		#3.G_mv_fastqs.py
		
		#3.H_samtools_mapped_reads.slurm
		#3.H_samtools_mapped_reads.py
		
		
	
