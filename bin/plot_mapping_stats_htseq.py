#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 18:23:21 2020

@author: B. Mika-Gospodorz

Input files: tsv file 
Output file: tsv and pdf files
Description: Used to plot mapping statistics for HTSeq
"""

import argparse
import pandas as pd
import matplotlib
#do not use Xwindows backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# function to plot mapping statistics
def plot_mapping_stats(df_comb,no_samples,profile, x_lab, xticks_np, percentage, n_col, m):
    # define color palette
    my_cmap2= matplotlib.cm.get_cmap('tab20')
    my_cmap4 = matplotlib.cm.get_cmap('tab20c')
    color = ['#af6f09','#840000','#f9bc08','#ef4026','#f7d560','#ff796c',my_cmap4.colors[6], my_cmap2.colors[15],my_cmap2.colors[14]]
     # reverse sample order
    df_comb = df_comb.loc[reversed(df_comb.index)]
    # make plot 
    fig = df_comb.plot(kind='barh', stacked=True,figsize=(60, no_samples + 3),legend=True,color = color, width = 0.8, fontsize=40)
    fig.spines['top'].set_visible(False)
    fig.spines['right'].set_visible(False)
    fig.spines['bottom'].set_visible(True)
    fig.spines['left'].set_visible(False)
    # set x label
    plt.xlabel(x_lab, fontsize=40,labelpad=15)
    # set x label ticks
    plt.xticks(xticks_np)
    if(percentage == False): # set scientific format of label ticks if there is no percentage values
        fig.ticklabel_format(style='scientific', axis='x',scilimits=(m,m)) 
        fig.xaxis.major.formatter._useMathText = True
        tx = fig.xaxis.get_offset_text()
        tx.set_fontsize(40)
        plt.xticks(xticks_np,rotation = 90, fontsize=40)
        plt.tick_params(axis='x', which='major', pad=10)
    # set legent position and format
    fig.legend(loc = 'upper center',ncol=n_col,bbox_to_anchor=(0.5,0),bbox_transform=plt.gcf().transFigure,frameon=False, prop={'size': 50})
    # set sample names as y labels
    fig.set_yticklabels(df_comb.index)
    # save used table
    df_comb.to_csv("mapping_stats_" + profile + ".csv",sep='\t')  
    # save plot
    plt.savefig('mapping_stats_star_' + profile + '.pdf', dpi = 300, orientation = 'landscape',transparent=False,bbox_inches='tight')


parser = argparse.ArgumentParser(description="""plot mapping statistics""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.tsv', help="Path to htseq mapping statistics tsv table ")
args = parser.parse_args()
    

# read mapping statistic table as data frame
df_stats = pd.read_csv(args.input_files,sep="\t",index_col=0)

# extract number of samples
no_samples = df_stats.shape[0]


# calculate percentage of mapping staistics
if 'trimmed_reads' in df_stats.columns: # if statistics contain information on number of trimmed reads
    pathogen_uniquely_mapped_percent = (df_stats['pathogen_assigned_reads']/df_stats['total_raw_reads']) * 100
    host_uniquely_mapped_percent = (df_stats['host_assigned_reads']/df_stats['total_raw_reads']) * 100
    unassigned_pathogen_reads = (df_stats['unassigned_pathogen_reads']/df_stats['total_raw_reads']) * 100
    unassigned_host_reads = (df_stats['unassigned_host_reads']/df_stats['total_raw_reads']) * 100
    pathogen_multi_mapped_percent = (df_stats['pathogen_multi_mapped_reads']/df_stats['total_raw_reads']) * 100
    host_multi_mapped_percent = (df_stats['host_multi_mapped_reads']/df_stats['total_raw_reads']) * 100
    cross_mapped_percent = (df_stats['cross_mapped_reads']/df_stats['total_raw_reads']) * 100
    unmapped_percent = (df_stats['unmapped_reads']/df_stats['total_raw_reads']) * 100
    trimmed_percent = (df_stats['trimmed_reads']/df_stats['total_raw_reads']) * 100
    
    df_comb = pd.concat([host_uniquely_mapped_percent, pathogen_uniquely_mapped_percent,unassigned_host_reads, unassigned_pathogen_reads,  host_multi_mapped_percent, pathogen_multi_mapped_percent, cross_mapped_percent, unmapped_percent, trimmed_percent], axis=1)
    df_comb.columns = ['host quantified reads','pathogen quantified reads','host unquantified reads', 'pathogen unquantified reads', 'host multi-mapped reads',  'pathogen multi-mapped reads', 'cross-mapped reads','unmapped reads','trimmed reads']
    # define no. of columns of plot legend 
    n_col = 5
else:
    pathogen_uniquely_mapped_percent = (df_stats['pathogen_assigned_reads']/df_stats['processed_reads']) * 100
    host_uniquely_mapped_percent = (df_stats['host_assigned_reads']/df_stats['processed_reads']) * 100
    unassigned_pathogen_reads = (df_stats['unassigned_pathogen_reads']/df_stats['processed_reads']) * 100
    unassigned_host_reads = (df_stats['unassigned_host_reads']/df_stats['processed_reads']) * 100
    pathogen_multi_mapped_percent = (df_stats['pathogen_multi_mapped_reads']/df_stats['processed_reads']) * 100
    host_multi_mapped_percent = (df_stats['host_multi_mapped_reads']/df_stats['processed_reads']) * 100
    cross_mapped_percent = (df_stats['cross_mapped_reads']/df_stats['processed_reads']) * 100
    unmapped_percent = (df_stats['unmapped_reads']/df_stats['processed_reads']) * 100
    df_comb = pd.concat([host_uniquely_mapped_percent, pathogen_uniquely_mapped_percent,unassigned_host_reads, unassigned_pathogen_reads,  host_multi_mapped_percent, pathogen_multi_mapped_percent, cross_mapped_percent, unmapped_percent], axis=1)
    df_comb.columns = ['host quantified reads','pathogen quantified reads', 'host unquantified reads', 'pathogen unquantified reads', 'host multi-mapped reads', 'pathogen multi-mapped reads','cross-mapped reads','unmapped reads']
    # define no. of columns of plot legend 
    n_col = 4


# plot mapping statistics expressed as percentages
plot_mapping_stats(df_comb,no_samples,'samples_percentage','[ % ]',np.arange(0, 105, step=5), percentage = True, n_col = n_col, m=0)
 



# plot mapping statistics 
if 'trimmed_reads' in df_stats.columns:  # if no. of trimmed reads is defined
    df_comb = pd.concat([df_stats['host_assigned_reads'],df_stats['pathogen_assigned_reads'],df_stats['unassigned_host_reads'],df_stats['unassigned_pathogen_reads'],df_stats['host_multi_mapped_reads'] ,df_stats['pathogen_multi_mapped_reads'], df_stats['cross_mapped_reads'],df_stats['unmapped_reads'], df_stats['trimmed_reads']], axis=1)
    df_comb.columns = ['host quantified reads','pathogen quantified reads', 'host unquantified reads', 'pathogen unquantified reads', 'host multi-mapped reads', 'pathogen multi-mapped reads', 'cross-mapped reads','unmapped reads','trimmed reads']
    # calculate max total no. of reads and define x label values, round to 7 digits
    max_limit = round(df_stats['total_raw_reads'].max(),-6)
    if max_limit > 0 :       # if total number of reads is higher than 10^6 specify label ticks adjusted to this magnitude
	    # divide max_limit by 20 (number of label ticks) to get step
	    step2 = int(max_limit / 20) 
	    # define label ticks
	    array_labels = np.arange(0, max_limit + step2, step=step2)
	    array_labels2= np.around(array_labels)
	    # set m value to specify format of scientific format of x label ticks
	    m = 6
    else:                   # define label ticks if number of reads is smaller than 10^6
	    step = int(df_stats['total_raw_reads'].max()/20)
	    array_labels = np.arange(0, df_stats['total_raw_reads'].max() + step, step=step)
	    array_labels2= np.around(array_labels)
	    # set m value to specify format of scientific format of x label ticks
	    m = 0
    # define no. of columns of plot legend 
    n_col = 5
else:
    df_comb = pd.concat([df_stats['host_assigned_reads'],df_stats['pathogen_assigned_reads'],df_stats['unassigned_host_reads'],df_stats['unassigned_pathogen_reads'], df_stats['host_multi_mapped_reads'], df_stats['pathogen_multi_mapped_reads'], df_stats['cross_mapped_reads'],df_stats['unmapped_reads']], axis=1)
    df_comb.columns = ['host quantified reads', 'pathogen quantified reads','host unquantified reads',  'pathogen unquantified reads', 'host multi-mapped reads','pathogen multi-mapped reads','cross-mapped reads','unmapped reads']
    # calculate max total no. of processed reads and define x label values, round to 7 digits
    max_limit = round(df_stats['processed_reads'].max(),-6)
    if max_limit > 0 :       # if total number of reads is higher than 10^6 specify label ticks adjusted to this magnitude
	    # divide max_limit by 20 (number of label ticks) to get step
	    step2 = int(max_limit / 20) 
	    # define label ticks
	    array_labels = np.arange(0, max_limit + step2, step=step2)
	    array_labels2= np.around(array_labels)
	    # set m value to specify format of scientific format of x label ticks
	    m = 6
    else:                   # define label ticks if number of reads is smaller than 10^6
	    step = int(df_stats['processed_reads'].max()/20)
	    array_labels = np.arange(0, df_stats['processed_reads'].max() + step, step=step)
	    array_labels2= np.around(array_labels)
	    # set m value to specify format of scientific format of x label ticks
	    m = 0
    # define no. of columns of plot legend 
    n_col = 4


# plot mapping statistics
plot_mapping_stats(df_comb,no_samples,'samples_total_reads','Number of reads', array_labels2, percentage = False, n_col = n_col,m = m)

