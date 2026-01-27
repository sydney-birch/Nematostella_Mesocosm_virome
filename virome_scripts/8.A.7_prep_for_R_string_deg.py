#! /usr/bin/env python3

#import modules
import argparse
import time

#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file: file continaing a singel column where each row is a gene you want to turn into a string")


args = parser.parse_args()

#This script will take in a file of a single column where each row you want to prep for R as a string 
# the final format will look like this: "header_1","header_2","header_3".....
#so you can import this string into R to make heatmaps 


	
tot_list = []
#part 1 populate list with headers
with open("{0}-R_string.txt".format(args.a), "w") as out_handle:
    with open(args.a, "r") as in_handle:
        for line in in_handle:
            sp_line = line.split("\t")
            #print("sp_line: ", sp_line)
            #print("sp_line[0]", sp_line[0])
            tot_list.append(sp_line[0])
            
        print("Full list populated: ", tot_list, len(tot_list))   

#part 2 join list in correct format for R and write to output file 
new_string = ",".join(tot_list)
print(new_string)

with open("{0}-R_string.txt".format(args.a), "w") as out_handle:
    out_handle.write("{0}".format(new_string))


## to double check it worked run wc -l on your input file and it should match the number in the tot list

