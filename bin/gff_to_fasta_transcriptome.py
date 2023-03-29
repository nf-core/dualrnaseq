#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 11:26:46 2019

@author: B. Mika-Gospodorz and R. Hayward

Input files: Fasta and GFF
Output file: Fasta
Description: Used to create a transcriptome fasta from user specified fields within a GFF file.

-----
For example, from the host, you can extract genes and gene IDs:
-----
python gff_to_fasta_transcriptome.py -fasta chr1.fa -gff chr1.gff -f gene -a gene_id -o chr1_genes.fa

-----
Or, transcripts and associated IDs:
-----
python gff_to_fasta_transcriptome.py -fasta chr1.fa -gff chr1.gff -f transcript -a ID -o chr1_transcripts.fa

"""

import argparse

from Bio import SeqIO


class FileNotProvided(Exception):
    pass


class NoFastaFileProvided(FileNotProvided):
    pass


class NoGffFileProvided(FileNotProvided):
    pass


def create_transcriptome(
    fasta_records_dict: dict, gff_files: list, features: list, gene_attribute: str, output_file_name: str
) -> None:
    if not len(gff_files):
        raise NoGffFileProvided
    with open(output_file_name, "a") as out_name:  # Open output file
        for gff_file in gff_files:
            for line in open(gff_file):
                if len(line.split()) == 0:  # ignore blank lines
                    continue
                d = line.rstrip()  # remove '\n'
                if (d[0] != "#") and (d != "") and (d != " "):  # Ignore comments
                    d_list = d.split("\t")  # Split based on tabs
                    if d_list[2] in features:
                        reference_name = d_list[0]
                        split_8 = d_list[8].split(";")  # Further split based on ;
                        # find index of id of interest
                        index_feature = [split_8.index(el) for el in split_8 if gene_attribute in el]
                        description = [s.split("=")[1] for s in split_8]
                        # re-order ids from column 8
                        if not index_feature:
                            print("lack of " + gene_attribute + " attribute for record:" + split_8[0].split("=")[1])
                        else:
                            out_name.write(">" + description[index_feature[0]] + "\n")
                            # Determine sequence based on direction
                            start = int(d_list[3]) - 1
                            stop = int(d_list[4])
                            if d_list[6] == "+":
                                out_name.write(str(fasta_records_dict[reference_name].seq[start:stop]) + "\n")
                            elif d_list[6] == "-":
                                out_name.write(
                                    str(fasta_records_dict[reference_name].seq[start:stop].reverse_complement()) + "\n"
                                )


def make_argument_parser():
    # Script arguments and descriptions
    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", nargs="+", help="genome fasta file")
    parser.add_argument("-gff", nargs="+", help="gff file")
    parser.add_argument(
        "-f", "--gene_feature", nargs="+", help="gene feature defined in the 3rd column of the gff file"
    )
    parser.add_argument("-a", "--gene_attribute", help="gene attribute")
    parser.add_argument("-o", help="output file name")
    return parser


def get_fasta_records_dict(fasta: list) -> dict:
    fasta_records_dict = dict()
    if not len(fasta):
        raise NoFastaFileProvided
    for fasta_file in fasta:
        fasta_records_dict.update(SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")))
    return fasta_records_dict


# pass fasta, features and attributes to function
if __name__ == "__main__":
    parser = make_argument_parser()
    args = parser.parse_args()
    gene_features = [feature.replace("[", "").replace("]", "").replace(",", "") for feature in args.gene_feature]
    fasta_records_dict = get_fasta_records_dict(args.fasta)
    create_transcriptome(fasta_records_dict, args.gff, gene_features, args.gene_attribute, args.o)
