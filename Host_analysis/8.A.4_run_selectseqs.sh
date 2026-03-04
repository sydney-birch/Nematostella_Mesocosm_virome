#! /bin/bash

#file of headers one per line=$1
 
for i in ./total_mapped_headers/*_full_header
do
  echo "Getting sequences for $i ...."
  ./selectSeqs.pl -f $i ./uk_transcriptome.fa >> $i-seqs.fa
done	

mkdir mapped_selectseqs_fasta_files
mv ./total_mapped_headers/*-seqs.fa ./mapped_selectseqs_fasta_files
