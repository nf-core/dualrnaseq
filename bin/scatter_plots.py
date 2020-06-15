#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 22:56:24 2020

@author: bozena
"""

import argparse
import pandas as pd
import matplotlib
#do not use Xwindows backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns
import itertools
from scipy import stats
import math



def scatter_plot(TPM_log, TPMs, name_1, name_2, lim_plot, organism):
    g1 = sns.JointGrid('rep1',y='rep2', data=TPM_log, height = 6)
    g1 = g1.plot_joint(plt.scatter,edgecolor="black", linewidth = 0.5)
    stat = lambda a,b: stats.pearsonr(TPMs.rep1,TPMs.rep2)
    g1 =g1.annotate(stat, template="{stat}: {val:.4f}", stat="Pearson's r", loc="upper left", fontsize=15)  # {stat}: {val:.2f} (p = {p:.3g}) with p-value
    g1.ax_marg_x.set_axis_off()
    g1.ax_marg_y.set_axis_off()
    #plt.rc('text', usetex=True)
    plt.xlabel(r'$\mathrm{log_{10}TPM}$'+ '\n' + '\n' +  name_1.split('_TPM')[0] ,fontsize=15, labelpad= 5)
    plt.ylabel(name_2.split('_TPM')[0] + '\n' + '\n'+  r'$\mathrm{log_{10}TPM}$', fontsize=15, labelpad =5)
    
    # https://matplotlib.org/3.1.1/tutorials/text/usetex.html
    #http://tug.org/texlive/quickinstall.html
    #  r'$\' {\fontsize{20pt}{3em}\selectfont{}[log_{10}(TPM+1)] .format(name_1.split('_TPM')[0]), fontsize=13)  {\fontsize{50pt}{3em}\selectfont{}a}
    # r'E_obs and E_syn @ t={0}, $Q_i^{{-1}}$={1}, $\ell$={2}'.format(time, q_intr, lpath) 
    # r'{\fontsize{30pt}{3em}\selectfont{}{Mean WRFv3.5 LHF\r} {\fontsize{18pt}{3em}\selectfont{}(September 16 - October 30, 2012)}'
    plt.xlim(-0.15, lim_plot)
    plt.ylim(-0.15, lim_plot)
    plt.title(organism,fontsize=15,fontweight="bold")
    plt.tick_params(axis="both", labelsize=15)
    plt.subplots_adjust(bottom=0.1)
    #plt.savefig('results/' + prefix + '_statter_1_2.svg', dpi = 300)
    plt.savefig('scatter_plot' + name_1 + '_' + name_2 + '_' + organism + '.pdf', dpi = 300,bbox_inches='tight')
    plt.close(plt.gcf())


parser = argparse.ArgumentParser(description="""RNA class statistics""")
parser.add_argument("-q", "--quantification_table", metavar='<quantification_table_host>', help="Path to the quantification table")
parser.add_argument("-a", "--gene_attribute", help="gene attribute")
parser.add_argument("-org", "--organism", help = 'host or pathogen')


args = parser.parse_args()

quantification_table_path = args.quantification_table
gene_attribute = args.gene_attribute

#samples_keys = [sample.replace('[' , '').replace(']','').replace(',','') for sample in args.conditions ]


col_names = pd.read_csv(quantification_table_path, sep = '\t', nrows=0).columns
types_dict = {gene_attribute: str}
types_dict.update({col: float for col in col_names if col not in types_dict})
quantification_table = pd.read_csv(quantification_table_path, sep = '\t',index_col=0,dtype=types_dict)

#extract TPM 
TMP_column = [column for column in quantification_table.columns if 'TPM' in column]

#find axes limits
TPM_table_plus_1 = quantification_table[TMP_column] + 1
TPM_table_log = TPM_table_plus_1.apply(np.log10, axis=1)
max_TPM = TPM_table_log.max().max()
lim_plot = math.ceil(max_TPM) + 0.5


#remove x_TPM
columns_conditions = [column[:-6] for column in TMP_column]


#find position of replicates
patterns = {}
for i in range(0,len(columns_conditions)):
    prefix = columns_conditions[i]
    d1 = {prefix:[i]}
    if prefix in patterns.keys():
        patterns[prefix].insert(patterns[prefix][0],list(d1.values())[0][0])
    else:
        patterns.update(d1)        

#list of replicates - positions
condition_list = [patterns[condition] for condition in patterns.keys()]

for cond in condition_list:
    sample_cond = [TMP_column[c] for c in cond]
    TPMs = quantification_table[sample_cond]
    combinations = list(itertools.combinations(TPMs, 2))
    for com in combinations: 
         TPMs = pd.concat([quantification_table[com[0]], quantification_table[com[1]]], axis=1)
         TPMs.columns = ['rep1','rep2']
         TPM_plus1 = TPMs + 1 # to deal with 0 values
         TPM_log = TPM_plus1.apply(np.log10, axis=1)
         #plot 
         scatter_plot(TPM_log,TPMs, com[0],com[1],lim_plot, args.organism)