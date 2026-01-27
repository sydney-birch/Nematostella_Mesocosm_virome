#! /bin/bash

for i in ./70_pi_blastout/*70_pi
do
  echo "getting hit accessions  $i ..."
  cut -f1 $i > ${i%.}_hits
done

mkdir 70_pi_hit1_accessions

mv ./70_pi_blastout/*hits ./70_pi_hit1_accessions
