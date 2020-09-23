#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 12:18:58 2020

@author: B. Mika-Gospodorz

Input files: tsv file with RNA class statistics for all samples
Output file: pdf file
Description: Used to plot RNA class mapping statistics for all samples
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


# function to assign color to each RNA class
def create_color_pallet(RNA_classes, colors):
    # define dictionary with RNA classes and colors
    color_map_collection = {}
    i = 0 
    for RNA_class in RNA_classes: # iterate over RNA classes and assign colors to them
        color_map_collection.update( {RNA_class: colors[i]} )
        i = i + 1
    return color_map_collection


# function to plot RNA class statistics 
def plot_mapping_stats(df_comb,no_samples, color_map_collection, profile):
    df_comb = df_comb.T
    # sort RNA classes by values of the first sample into descending order
    df_comb = df_comb.sort_values(by= [(df_comb.columns[0])], ascending=False)
    df_comb = df_comb.T
    # reverse the sample order
    df_comb = df_comb.loc[reversed(df_comb.index)]
    # define no. of RNA classes
    n_col = len(df_comb.columns)
    # define no. of columns of plot legend depending on no. of RNA classes
    if n_col <= 8 and n_col % 8 == 0:
       ncol_leg = n_col
    else:
        ncol_leg = 6
    # define color map based on RNA class order in RNA class statistic data frame
    color_map = []           
    for key in df_comb.columns:
        color_map.append(color_map_collection[key])
    # plot statistics
    fig = df_comb.plot(kind='barh', stacked=True,figsize=(60, no_samples + 3),legend=True,color = color_map, width = 0.8, fontsize=40)
    fig.spines['top'].set_visible(False)
    fig.spines['right'].set_visible(False)
    fig.spines['bottom'].set_visible(True)
    fig.spines['left'].set_visible(False)
    plt.xlabel('[ % ]', fontsize=40)
    plt.xticks(np.arange(0, 105, step=5),rotation = 0)
    fig.legend(loc = 'lower center',ncol=ncol_leg,bbox_to_anchor=(0.5,1),bbox_transform=plt.gcf().transFigure,frameon=False, prop={'size': 50})
    fig.set_yticklabels(df_comb.index)
    plt.savefig('RNA_class_stats_combined_' + str(profile) + '.pdf', dpi = 300, orientation = 'landscape',transparent=False,bbox_inches='tight')



##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
parser = argparse.ArgumentParser(description="""plot RNA_class_stats""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.tsv', help="Path to RNA class mapping statistics tsv table")
parser.add_argument("-org", "--organism", metavar='<organism>', help = 'host or pathogen')
args = parser.parse_args()



# load colormaps
my_cmap = matplotlib.cm.get_cmap('tab20c')
my_cmap2 = matplotlib.cm.get_cmap('Set2')
my_cmap3 = matplotlib.cm.get_cmap('tab20b')
my_cmap4 = matplotlib.cm.get_cmap('tab20')

# define list of colors
colors = [my_cmap.colors[0],my_cmap.colors[3],my_cmap.colors[1],my_cmap2.colors[7],my_cmap.colors[2],my_cmap2.colors[6],
          my_cmap.colors[4],my_cmap2.colors[5],my_cmap2.colors[0],my_cmap2.colors[2],my_cmap.colors[5],my_cmap4.colors[6],
          my_cmap4.colors[19],my_cmap4.colors[16],my_cmap4.colors[18],my_cmap4.colors[10],my_cmap3.colors[13],my_cmap3.colors[14],my_cmap3.colors[16],
          my_cmap3.colors[12],my_cmap3.colors[8],my_cmap3.colors[0], my_cmap3.colors[4],my_cmap2.colors[3],my_cmap2.colors[4]]


# read RNA class mapping statistics 
df_stats = pd.read_csv(args.input_files,sep="\t",index_col=0)
# remove _NumReads suffix from sample name
new_index = [name.split('_NumReads')[0] for name in df_stats.index ]
df_stats.index = new_index

# remove mock samples for pathogen
if args.organism == 'pathogen':
    mock_samples = [sample_name for sample_name in df_stats.index if 'mock' in sample_name]
    df_stats = df_stats.drop(mock_samples, axis=0)

# create dictionary with RNA classes and colors assigned to them 
colour_map = create_color_pallet(df_stats.columns,colors)

# create data frame with RNA class statistics
df_comb = df_stats

# extract no. of samples
no_samples = df_stats.shape[0]

# plot results
plot_mapping_stats(df_comb,no_samples,colour_map, args.organism)




