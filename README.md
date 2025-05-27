# Nematostella Mesocosm Virome/Host Transcriptome Study - Overview
This repository contains the virome, host transcriptome, and 16s rRNA work for our 2022 Mesocosm study looking at how Nematostella from different locations interact with microbes and viruses from a foreign environment. 

### Study Info and Sample Structure: 
We have sequencing data on animals from 7 locations: Field (NH)	; 	NS	;	ME	;	NH	;	MA 	; SC ;	FL for two time points: T0 and T14 with 3 replicates for each sample. The 16s dataset contains 5 replicates for each sample. 

Where NS	;	ME	;	NH	;	MA 	; SC ;	FL animals were sampled at our intial time point (T0) before going into the mesocosm which was conducted in NH at the Jackson Estuarine Laboratory. The mesocosms were bins filled with water from the estuary that was bag filtered to remove larger debris. Water changes occured every other day, and each mesocosm bin had a pump to produce water flow for the animals. Animals were exposed to the NH estuarine water for 14 days when we collected our last time point (T14). In addition to the Mesocosom sampling, we also collected animals from the natural populations in NH (Field animals) to use as a control. 

Our Main study question is:   
**Does the origin of location impact the type and diversity of viruses and bacteria that assoiciate with Nematostella?** 

### Virome Pipeline overview:  

1) We extracted RNA and DNA using the Allprep RNA/DNA mini kit
2) We created libraries for total RNAseq using the Tecan Universal Plus total RNA-seq with NuQuant kit with an N. vectensis rRNA depletion step
3) We sequenced samples using NovaseqX
4) Prep work: Run Fastqc, trim reads, and run fastqc again
5) Map reads to the genome to get two pots of data --> mapped (Nematostella reads) and unmapped (virome/microbe reads)    
- For virome (unmapped reads) - map back to the genome multiple times until 0 reads align
- For Nematostella transcriptome (mapped reads) - use mapped reads from 1st alignment

Virome:
1) map the unmapped reads back to the genome until 0 reads align
2) Align reads to rRNA database comprised of SILVA LSU & SSU, GTDB SSU, and Rfam --> unmapped reads presumably viral
3) Assemble viral transcriptome using Spades rnaviral (created 14 assemblies, collapsing replicates)
4) Analyze Function
.... A) Run MetaCerberus on assemblies
.... B) Get functional terms: KEGG, COG, VOG, PHROG
.... C) Visualize in R 
5) Analyze Taxonomy
.... A) BLAST transcriptomes to RefSeq Viruses
.... B) Filter BLAST output (>= 70% identity)
.... C) BioEntrez and NCBI Taxon Kit to get taxon IDs and lineage info
.... D) visualize in R &  run diversity metrics 

