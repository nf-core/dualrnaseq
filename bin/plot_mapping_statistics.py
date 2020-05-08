#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 12:57:05 2020

@author: bozena
"""
import argparse
import pandas as pd
import matplotlib
#do not use Xwindows backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_mapping_stats(df_comb,no_samples,profile, x_lab, xticks_np, percentage):
    my_cmap2= matplotlib.cm.get_cmap('tab20')
    color = [my_cmap2.colors[1], my_cmap2.colors[0],my_cmap2.colors[15],my_cmap2.colors[14]] 
    df_comb = df_comb.loc[reversed(df_comb.index)]
   # plt.yticks(yint)
    fig = df_comb.plot(kind='barh', stacked=True,figsize=(60, no_samples * 2 + 3),legend=True,color = color, width = 0.8, fontsize=40)
    fig.spines['top'].set_visible(False)
    fig.spines['right'].set_visible(False)
    fig.spines['bottom'].set_visible(True)
    fig.spines['left'].set_visible(False)
    plt.xlabel(x_lab, fontsize=40,labelpad=15)
    plt.xticks(xticks_np)
    if(percentage == False):
        fig.ticklabel_format(style='scientific', axis='x',scilimits=(0,0)) 
        fig.xaxis.major.formatter._useMathText = True
        tx = fig.xaxis.get_offset_text()
        tx.set_fontsize(40)
        plt.xticks(xticks_np,rotation = 90, fontsize=40)
        plt.tick_params(axis='x', which='major', pad=10)
    plt.ylabel('')
    fig.legend(loc = 'upper center',ncol=4,bbox_to_anchor=(0.5, 1),bbox_transform=plt.gcf().transFigure,frameon=False, prop={'size': 50})
    fig.set_yticklabels(df_comb.index)
    df_comb.to_csv("mapping_stats_" + profile + ".csv",sep='\t')  
    plt.savefig('mapping_stats_' + profile + '.pdf', dpi = 300, orientation = 'landscape',transparent=False,bbox_inches='tight')




parser = argparse.ArgumentParser(description="""plot mapping statistics""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.csv', help="Path to the outputfiles from mapping statistics ")
#parser.add_argument("-p", "--profile", metavar='<profile_of_mapping>', help="multi/uniquely -mapped")
args = parser.parse_args()
    

df_stats = pd.read_csv(args.input_files,sep="\t",index_col=0)

##percentage
pathogen_percent = (df_stats['pathogen']/df_stats['processed_reads']) * 100
host_percent = (df_stats['host']/df_stats['processed_reads']) * 100
unmapped_percent = (df_stats['unmapped_reads']/df_stats['processed_reads']) * 100
trimmed_percent = (df_stats['trimmed_reads']/df_stats['processed_reads']) * 100

df_comb = pd.concat([host_percent,pathogen_percent,unmapped_percent, trimmed_percent], axis=1)
df_comb.columns = ['host','pathogen','unmapped reads','trimmed reads']
no_samples = df_stats.shape[0]
plot_mapping_stats(df_comb,no_samples,'samples_percentage','[ % ]',np.arange(0, 105, step=5), percentage = True)
    
##total
df_comb = pd.concat([df_stats['host'],df_stats['pathogen'],df_stats['unmapped_reads'], df_stats['trimmed_reads']], axis=1)
df_comb.columns = ['host','pathogen','unmapped reads','trimmed reads']
no_samples = df_stats.shape[0]

step = int(df_stats['total_raw_reads'].max()/50)
step2 = int((step + 50) /100) * 100

plot_mapping_stats(df_comb,no_samples,'samples_total_reads','Number of reads',np.arange(0, df_stats['total_raw_reads'].max() + step2, step=step2), percentage = False)




