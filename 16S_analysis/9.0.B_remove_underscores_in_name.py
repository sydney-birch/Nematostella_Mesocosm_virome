#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os
import gzip


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-b", help="directory to navigate to with fastq files")
args = parser.parse_args()

def name_change (dir_name):
    os.chdir(dir_name)
    for item in os.scandir():
        print("iterating on item: ", item) #for debug
   
        sample_name = item.name
        print("sample Name: ", sample_name)
            
        sp_name = sample_name.split("_")
        print ("split name: ", sp_name)
              
        temp_list = []
        temp_list.append(sp_name[0])
        temp_list.append(sp_name[1])
        temp_list.append(sp_name[2])
        
        new_name = "-".join(temp_list)
        print("new name: ", new_name)
        
        temp_list_2 = []
        temp_list_2.append(new_name)
        temp_list_2.append(sp_name[3])
        temp_list_2.append(sp_name[4])
        temp_list_2.append(sp_name[5])
        temp_list_2.append(sp_name[6])
        
        final_new_name = "_".join(temp_list_2)
        print("final new name: ", final_new_name)
                
        result = subprocess.run("mv {0} {1}".format(sample_name, final_new_name), shell=True)
        print("name successfully changed: ", final_new_name)
               
            
#call run_fastqc function
result = name_change(args.b)
