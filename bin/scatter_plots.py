#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 22:56:24 2020

@author: B. Mika-Gospodorz

Input files: quantification tsv file 
Output file: pdf files
Description: Used to plot scatter plots of technical replicates
"""


import argparse
import pandas as pd
import matplotlib
#do not use Xwindows backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import itertools
from scipy import stats
import math


# function to plot scatter plots
def scatter_plot(TPM_log, TPMs, name_1, name_2, lim_plot, organism):
    g1 = sns.JointGrid('rep1',y='rep2', data=TPM_log, height = 6)
    g1 = g1.plot_joint(plt.scatter,edgecolor="black", linewidth = 0.5)
    # calculate pearson's correlation coefficient
    stat = lambda a,b: stats.pearsonr(TPMs.rep1,TPMs.rep2)
    # add label with pearson's correlation coefficient
    # g1 =g1.annotate(stat, template="{stat}: {val:.4f}", stat="Pearson's r", loc="upper left", fontsize=15)  # {stat}: {val:.2f} (p = {p:.3g}) with p-value
    # that 0.9.0 only, the following should ALSO work for v0.11.2
    pearr = stats.pearsonr(TPMs.rep1, TPMs.rep2)
    st=f"Pearson's r = {pearr[0]:.4f}"
    g1.ax_joint.text(0.05, 0.95, st, fontsize=15, transform=g1.ax_joint.transAxes, verticalalignment='top')
    # return to main
    g1.ax_marg_x.set_axis_off()
    g1.ax_marg_y.set_axis_off()
    plt.xlabel(r'$\mathrm{log_{10}TPM}$'+ '\n' + '\n' +  name_1.split('_TPM')[0] ,fontsize=15, labelpad= 5)
    plt.ylabel(name_2.split('_TPM')[0] + '\n' + '\n'+  r'$\mathrm{log_{10}TPM}$', fontsize=15, labelpad =5)
    plt.xlim(-0.15, lim_plot)
    plt.ylim(-0.15, lim_plot)
    plt.title(organism,fontsize=15,fontweight="bold")
    plt.tick_params(axis="both", labelsize=15)
    plt.subplots_adjust(bottom=0.1)
    # save plot
    plt.savefig('scatter_plot_' + name_1 + '_' + name_2 + '_' + organism + '.pdf', dpi = 300,bbox_inches='tight')
    plt.close(plt.gcf())



parser = argparse.ArgumentParser(description="""Plots scatter plots of replicates""")
parser.add_argument("-q", "--quantification_table", metavar='<quantification_table_host>', help="Path to the quantification table")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-org", "--organism", help = 'host or pathogen')


args = parser.parse_args()
quantification_table_path = args.quantification_table
gene_attribute = args.gene_attribute

# read quantification table as data frame
col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
types_dict = {gene_attribute: str}
types_dict.update({col: float for col in col_names if col not in types_dict})
quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)

# extract columns with 'TPM' suffix
TMP_column = [column for column in quantification_table.columns if 'TPM' in column]

#find axes' limits
TPM_table_plus_1 = quantification_table[TMP_column] + 1
TPM_table_log = TPM_table_plus_1.apply(np.log10, axis=1)
max_TPM = TPM_table_log.max().max()
lim_plot = math.ceil(max_TPM) + 0.5


# remove '_TPM' suffix from sample names
columns_conditions_no_TPM = [column[:-4] for column in TMP_column]

# extract conditions - remove indication of which replicate it is (_1, _2 or _3)
columns_conditions = [name.rsplit('_', 1)[0] for name in columns_conditions_no_TPM]
# find positions of technical replicates in quantification table
patterns = {}
for i in range(0,len(columns_conditions)): #iterate over number of samples
    prefix = columns_conditions[i] # extract condition
    d1 = {prefix:[i]} # dictionary with condition and its position in columns conditions, which is equivalent to positions in quantification table
    if prefix in patterns.keys(): # if condition already exist in the dictionary, then append the position
        patterns[prefix].insert(patterns[prefix][0],list(d1.values())[0][0])
    else:   # append condition and its position
        patterns.update(d1)        


# create list of replicates' position lists, e.g. [[1, 2, 0], [3, 4, 5], [6, 7, 8]]
condition_list = [patterns[condition] for condition in patterns.keys()]
# plot scatter plots of technical replicates
for cond in condition_list: # iterate over list of technical replicates' positions 
    if len(cond) > 1:
        sample_cond = [TMP_column[c] for c in cond] # extract column names for replicates in quantification table
        TPMs = quantification_table[sample_cond]
        combinations = list(itertools.combinations(TPMs, 2)) # create list of possible combinations of sample pairs
        for com in combinations:  # iterate over sample pair combinations
              TPMs = pd.concat([quantification_table[com[0]], quantification_table[com[1]]], axis=1) # extract TPM values of two replicates
              TPMs.columns = ['rep1','rep2']
              TPM_plus1 = TPMs + 1 # add 1 to deal with 0 values
              TPM_log = TPM_plus1.apply(np.log10, axis=1) # log-transformation
              #plot scatter plot
              scatter_plot(TPM_log,TPMs, com[0],com[1],lim_plot, args.organism)
