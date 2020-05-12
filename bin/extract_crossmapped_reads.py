#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:37:04 2019

@author: bozena
"""

import sys
import argparse
import pysam


def check_reference_organisms(read_reference_name, host_reference_names,pathogen_reference_names):
    reference = ''
    if read_reference_name in host_reference_names:
        reference = 'host'
    elif read_reference_name in pathogen_reference_names:
        reference = 'pathogen'
    else:
        print(('There is no ' + read_reference_name + ' in the reference name set'))
    return reference


def add_read(multimapped_reads, read_name, reference_name):
    # if read is in the dict, append a new reference
    if read_name in multimapped_reads:
       if reference_name not in multimapped_reads[read_name]: 
           multimapped_reads[(read_name)].append((reference_name))
    # else create a new key
    else:
        multimapped_reads[read_name] = [reference_name]


def find_and_save_cross_mapped_reads(multimapped_reads_with_reference, output_file_name):
    #crossmapeed_reads = [{read_name:reference_name} for read_name, reference_name in multimapped_reads_with_reference.iteritems() if len(reference_name) > 1]
    crossmapeed_reads = [read_name for read_name, reference_name in list(multimapped_reads_with_reference.items()) if len(reference_name) > 1]
    with open(output_file_name, 'w') as f:
        for cross_mapped_read in crossmapeed_reads:
            f.write(str(cross_mapped_read) + '\n')



def read_reads_from_samfile(sam_file_name, host_reference_names, pathogen_reference_names, output_file_name): 
    samfile = pysam.AlignmentFile(sam_file_name, "rb")  
    multimapped_reads = dict()
    for read in samfile:
        reference_organisms = check_reference_organisms(read.reference_name,host_reference_names,pathogen_reference_names)
        add_read(multimapped_reads,read.query_name,reference_organisms)
    find_and_save_cross_mapped_reads(multimapped_reads,output_file_name)
        

    
parser = argparse.ArgumentParser()


parser.add_argument("-o", "--output", default="cross_mapped_reads.txt", help="output file name")
parser.add_argument("-h_ref", "--host_reference_name", metavar='<host_reference_name>', help="Path to the host reference gtf, multiple gtfs (read as a list)")
parser.add_argument("-p_ref", "--pathogen_reference_name", metavar='<pathogen_reference_name>', help="Path to the pathogen reference gtf ")



args = parser.parse_args()

host_reference_names = [line.rstrip() for line in open(args.host_reference_name)]
pathogen_reference_names = [line.rstrip() for line in open(args.pathogen_reference_name)]


read_reads_from_samfile(sys.stdin, host_reference_names, pathogen_reference_names, args.output)



