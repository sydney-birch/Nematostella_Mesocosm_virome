#! /usr/bin/env python3

#import modules
#module load anaconda3
from Bio import Entrez
import argparse
import time

#create an instance of Argument Parser and add positional argument
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with refseq accids - one id per line")
parser.add_argument("-o", help="Name of output file group with one taxid per line (for example: NH-T14)")

args = parser.parse_args()


#This script will open up the input file of refseq accids (one per line) and will iterate 
#through to get the taxid using bio.entrez elink function (nuccore the nuc db and the taxonomy db)
#It will write out 2 files - one with only the TaxIds using the -o file name and another file
#with the TaxId and corresponding refseq Id Tab seperated (the -o name with accid-taxid)

#The list that is outputed from elink for reference
#{'LinkSetDb': [{'Link': [{'Id': '445686'}], 'DbTo': 'taxonomy', 'LinkName': 'nuccore_taxonomy'}]


def get_taxids(input_file,output_file_name):
    #open input file with refseq ids - one per line
    with open(input_file, "r") as in_handle:
        #open output file for Taxa ids - one per line
        with open("{0}_taxid.txt".format(output_file_name), "w") as out_handle_1:
            #open output file for Taxa ids and accids - one per line
            with open("{0}_accid-taxid.txt".format(output_file_name), "w") as out_handle_2:
                try:
                    #time.sleep(5)
                    for ref_id in in_handle:
                        #time.sleep(3)
                        #print("ref id: ", ref_id)
                        Entrez.email = "sbirch1@uncc.edu"
                        NCBI_API_Key='282329c7fc748e86c370b04aa67dfdd6da09'
                        nuccoreid = ref_id.strip()
                        print("nuccoreid: ",nuccoreid)
                        handle = Entrez.elink(dbfrom="protein", id=nuccoreid, linkname="protein_taxonomy")
                        record = Entrez.read(handle)
                        handle.close()
                        #print(record[0]["LinkSetDb"][0]["Link"][0]["Id"])
                        taxid = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
                        print("taxid: ", taxid)
                        #Write out taxid to taxid only file:
                        out_handle_1.write("{0}\n".format(taxid))
                        #Write out taxid and refseq accid to combined file:
                        out_handle_2.write("{0}\t{1}\n".format(nuccoreid,taxid))
                        handle.close()
                        time.sleep(12)
                except RuntimeError as err:
                    print("Run time error, was last on: ", ref_id)         
#call funciton
result = get_taxids(args.i,args.o)        
        


#This works for doing entrez search - use refseq id to get taxid 
#Entrez.email = "sbirch1@uncc.edu"
#nuccoreid = "NC_027384.1"
#handle = Entrez.elink(dbfrom="nuccore", id=nuccoreid, linkname="nuccore_taxonomy")
#record = Entrez.read(handle)
#handle.close()
#print(record[0]["LinkSetDb"][0]["Link"][0]["Id"])
