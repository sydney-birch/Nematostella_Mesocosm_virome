#! /usr/bin/env python3

#import modules 
import argparse
import subprocess
import os


#create an instance of Argument Parser and add positional argument 
parser = argparse.ArgumentParser()
parser.add_argument("-a", help="input file fron edgeR")
args = parser.parse_args()

header_list=[]

def get_headers(input_file):
    with open("{0}".format(input_file), "r") as in_handle: 
        for line in in_handle: 
            sp_line = line.split("\"")
            #print("first split: ", sp_line)
            #print("first item: ", sp_line[1])
            header = sp_line[1]
            #print("header: ", header)
            header_list.append(header)

        #print("complete list: ", header_list)
        header_list.remove("logFC")
        print("final list: ", header_list)

    with open("{0}-header_list".format(input_file),"w") as out_handle: 
        for header in header_list:
            out_handle.write("{0}\n".format(header))


#call funciton
result = get_headers(args.a) 


#write out file
#with open("{0}-header_list".format(input_file),"w") as out_handle: 
#    for header in header_list:
#        out_handle.write(header)
