#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:37:04 2019

@author: B.Mika-Gospodorz

Input files: stdin with bam file with multi-mapped reads, reference_host_names.txt and reference_pathogen_names.txt files that contain references extracted with extract_reference_names_from_fasta_files.sh
Output: txt file with cross-mapped reads
Description: Used to identify and extract reads that mapped onto both host and pathogen genomes (cross-mapped reads). The script is executed by remove_crossmapped_reads_BAM.sh and remove_crossmapped_read_paires_BAM.sh scripts.
"""

import argparse
import sys

import pysam


# function to identify references given read mapped to
def check_reference_organisms(read_reference_name, host_reference_names, pathogen_reference_names):
    reference = ""
    if (
        read_reference_name in host_reference_names
    ):  # if reference of read is in the list of host references, set host as reference
        reference = "host"
    elif (
        read_reference_name in pathogen_reference_names
    ):  # if reference of read is in the list of pathogen references, set pathogen as reference
        reference = "pathogen"
    else:
        print(("There is no " + read_reference_name + " in the reference name set"))
    return reference


# function to add read and its reference to dictionary
def add_read(multimapped_reads, read_name, reference_name):
    if read_name in multimapped_reads:  # if read is in the dict, and reference is not defined, append the reference
        if reference_name not in multimapped_reads[read_name]:
            multimapped_reads[(read_name)].append((reference_name))
    else:  # else create new key (read) and set reference as value
        multimapped_reads[read_name] = [reference_name]


# function to find and save cross-mapped reads
def find_and_save_cross_mapped_reads(multimapped_reads_with_reference, output_file_name):
    # extract reads with more than 1 reference name defined (cross-mapped reads)
    crossmapeed_reads = [
        read_name
        for read_name, reference_name in list(multimapped_reads_with_reference.items())
        if len(reference_name) > 1
    ]
    # save cross-mapped reads
    with open(output_file_name, "w") as f:
        for cross_mapped_read in crossmapeed_reads:
            f.write(str(cross_mapped_read) + "\n")


# function to identify and save cross-mapped reads
def read_reads_from_samfile(sam_file_name, host_reference_names, pathogen_reference_names, output_file_name):
    # read bam file
    samfile = pysam.AlignmentFile(sam_file_name, "rb")
    # initialize dictionary of multi-mapped reads
    multimapped_reads = dict()
    for read in samfile:  # iterate over reads from bam file
        # find reference the read mapped to
        reference_organisms = check_reference_organisms(
            read.reference_name, host_reference_names, pathogen_reference_names
        )
        # add read and reference to multimapped_reads dict.
        add_read(multimapped_reads, read.query_name, reference_organisms)
    # find and save cross-mapped reads
    find_and_save_cross_mapped_reads(multimapped_reads, output_file_name)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--output", default="cross_mapped_reads.txt", metavar="output_file_name", help="output file name"
)
parser.add_argument(
    "-h_ref", "--host_reference_names", metavar="<host_reference_name>", help="Path to reference_host_names.txt file"
)
parser.add_argument(
    "-p_ref",
    "--pathogen_reference_names",
    metavar="<pathogen_reference_name>",
    help="Path to reference_pathogen_names.txt file",
)
args = parser.parse_args()

# create list of host and pathogen chromosome/plasmid names
host_reference_names = [line.rstrip() for line in open(args.host_reference_names)]
pathogen_reference_names = [line.rstrip() for line in open(args.pathogen_reference_names)]


# identify and extract cross-mapped reads
read_reads_from_samfile(sys.stdin, host_reference_names, pathogen_reference_names, args.output)
