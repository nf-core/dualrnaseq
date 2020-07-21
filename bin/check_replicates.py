#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 14:55:44 2020

@author: bozena
"""

import argparse
parser = argparse.ArgumentParser(description="""Are there replicates?""")
parser.add_argument("-s", "--sample_names", metavar='<input_files>', nargs='+')
args = parser.parse_args()

sample_names = [feature.replace('[' , '').replace(']','').replace(',','') for feature in args.sample_names]
suffixes = [name.rsplit('_', 1)[1] for name in sample_names]

print("R1" in suffixes and "R2" in suffixes or "r1" in suffixes and "r2" in suffixes)

