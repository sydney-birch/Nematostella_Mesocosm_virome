#! /bin/bash

#edgeR output file=$1
 
for i in ./mapped_DEG_files/*_DEGs_edgeR_output.txt
do
  echo "Getting headers for $i ...."
  ./get_headers_from_edgeR.py -a $i 
done	
