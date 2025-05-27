#! /bin/bash

#load modules 
module load cd-hit
 
#Concatenate the files in given dir ($1) and name file $2
cat $1/*.fa > $2.fa

#Run cd-hit on the concatenated file 
cd-hit -i $2.fa -o $2_ed1.fa -c 1.00 -M 1100    #cdhit --> looks for and removes redundant seq - c is threshold 

#example of running code: ./4.B_cat_cd-hit.sh RNA_db rRNAdb
