#! /bin/bash

full_fasta_file=$1
 
for i in ./70_pi_hit1_accessions/*hits
do
  echo "Getting full headers for $i ...."
  ./3.B_get_full_headers.py -a $i -f $1 -o ${i%.}_total
done

mkdir 70_pi_total_header_hit1_accessions
mv ./70_pi_hit1_accessions/*_total ./70_pi_total_header_hit1_accessions
