# Overview of nematostella mesocosom transcriptome pipline

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

  *Note: 2.A was run on the samples sequenced from Admera, 2.A.2 was run on the 8 samples that were sequenced at UNCC due to differences in naming conventions*
  
  `./2.A.2_trimmomatic.py -a ../raw_reads_10-3-23 -b ../2_trimmomatic`



   B) Re-Run Fastqc in the fastqc dir (2.B_fastqc.slurm) 
   
`./2.B_fastqc.py -a ../1_fastqc/after_trim -b ../2_trimmomatic`
      
*output goes to fastqc dir - after_trim dir*  

 
	
## 3) Map reads to genome to get two pots of data --> mapped and unmapped
	Run HISAT2 and use sam tools to get mapped (1st alignment) and unmapped reads
		Re-align the unmapped - run 2 alignments with HISAT2 then 2 alignments with bowtie
		
	#1st alignment - HISAT2 
		#1st Ali - (run slurm in 1st_ali_2-20-24): 
			./3.A_initial_hisat2_v2.py -a ../../2_trimmomatic -b 1 -c 2nd_ali -d ../3_HISAT2/1st_ali_2-20-24/
				## Output: sam file with _ali_1 in a dir in next ali dir (2nd_ali dir) - it makes the second ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_1_{1}.sam)		
		
		# run sam tools - make fastq files to be used in second ali (Run slurm in 2nd_ali dir)
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 1
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_1_ali_paired1.fq)
						

	#2nd alignment - HISAT2 
		#2nd Ali - Run slurm in 2nd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.C_hisat2_2_ali.slurm
			./3.C_hisat2_next_ali_v2.py -b 2 -c 3rd_ali_BT 
		        	## Output: sam file with _ali_2 in a subdir in next ali dir (3rd_ali dir) - it makes the third ali dir and all sub dirs which contain sam files (Nematostella_genome_ali_2_{1}.sam)
		        		
		# run sam tools - make fastq files to be used in third ali (Run slurm in 3rd_ali dir) (started 9:30)
			sbatch 3.B_samtools_processing.slurm
			./3.B_samtools_processing.py -b 2
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_2_ali_paired1.fq)     	
		        	
		        	
		        	
	#3rd alignment - bowtie
		#First - copy over genome file ending in .fna (in 1st_ali dir) - need to make a bowtie index
			
		#3rd Ali - Run slurm in 3rd_ali dir (Have to hard code where index is - In 1st ali dir - adjust name if needed)
			sbatch 3.D_bowtie_ali_3.slurm
			./3.D_bowtie_ali_initial.py -b 3 -c 4th_ali_BT
				# Output: sam file with BT_ali_3 in a dir in next ali dir (4th_ali dir) - it makes the fourth ali dir and all sub dirs which contain sam files (Nematostella_genome_BT_ali_3_{1}.sam)
			
		# run sam tools - make fastq files to be used in fourth ali (Run slurm in 4th_ali dir)	
			sbatch 3.E_samtools_processing_BT.slurm
			./3.E_samtools_processing_BT.py -b 3
				## Output: 2 fastq files in each sample name dir (unaligned_{sample_name}_3_ali_paired1.fq)
		
	
		
		
	
