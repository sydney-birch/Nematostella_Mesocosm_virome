#! /bin/bash

for i in ./blastout/*blastout
do
  echo "trimming blastout to 70%  $i ..."
  ./1.B_trim_blastout_70_pident.py -i $i

done

mkdir 70_pi_blastout

mv ./blastout/*70_pi ./70_pi_blastout
