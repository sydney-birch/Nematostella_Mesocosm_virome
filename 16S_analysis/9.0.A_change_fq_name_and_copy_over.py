#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="file that has the names of the admera names and my names to make dictionary")
parser.add_argument("-b", help="directory to navigate to with fastq files")

args = parser.parse_args()

# renames the fastq files from new sequencing from admera 

#make a dictionary of the term from admera to the name that I want {key:Value} {admera name:My name} 
#iterate through each file in the dir 
#split the name on _ --> search in dictionary for {0} --> join with dictionary value (my name), {spots: 1,2,3,4} with _

#  21081FL-06-02-11_S27_L002_R1_001.fastq.gz


#Step 1: make dictionary from input file of file names to change 

names_db = {}

with open(args.a, "r") as in_handle:
    for line in in_handle:
        line = line.rstrip()
        sp_line = line.split("\t")
        #print("my name: ", sp_line[0])
        #print("name to change: ", sp_line[1])
        
        names_db.setdefault(sp_line[1],sp_line[0])
        
    print("dictionary: ", names_db)

#Step 2: iterate through each file in the directory and change the names    

def name_change (dir_name, name_db):
    os.chdir(dir_name)
    for item in os.scandir():
        print("iterating on item: ", item) #for debug
        #if "Undetermined" in item.name:
         #   continue

        if item.is_dir():
            sample_name = item.name
            print("potential sample Name: ", sample_name)
            
            sp_name = sample_name.split("_")
            print ("split name: ", sp_name)
            
            if sp_name[0] in name_db:
                print("value of key - new name: ", name_db[sp_name[0]])
           
                new_name = name_db[sp_name[0]]
                print("new name: ", new_name)
            
                os.chdir(sample_name)
                for fq in os.scandir(): 
                    fq_name = fq.name
                    sp_name_2 = fq_name.split("_")
                
                    temp_list = []
                    temp_list.append(new_name)
                    temp_list.append(sp_name_2[1])
                    temp_list.append(sp_name_2[2])
                    temp_list.append(sp_name_2[3])
                    temp_list.append(sp_name_2[4])
                
                    fq_new_name = "_".join(temp_list)
                    print("Total fq new name: ", fq_new_name)
                
                    result = subprocess.run("cp {0} /scratch/sbirch1/QIIME_Meso_2022/raw_reads/{1}".format(fq_name, fq_new_name), shell=True)
                    print("name successfully changed: ", new_name)
                os.chdir("..")
            
#call run_fastqc function
result = name_change(args.b, names_db)
