#!/bin/bash
module load blast

for i in ./fastas/*.fa
do
    echo "formating for BLAST  $i ..."
    makeblastdb -in $i -parse_seqids -dbtype nucl
done

# move blastdbs to new dir
mkdir blastdb
cp ./fastas/* blastdb

# clean up fastas dir
rm ./fastas/*.fa.*

#db types: nucl or prot
