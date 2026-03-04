# 16S analysis workflow

This is picking up after the host transcriptome analysis. 

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


