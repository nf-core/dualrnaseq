#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 17:21:53 2020

@author: B. Mika-Gospodorz

Input files: tsv file with RNA class statistics
Output files: pdf files
Description: Used to plot RNA class mapping statistics for each sample separately
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


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
def plot_RNA_class_stats_each(sample,color_map_collection): 
        sample_name = sample.index[0]
        title_plot = sample_name.split('_NumReads')[0]
        sample = sample.T
        # sort RNA classes by values of the sample into descending order
        sample = sample.sort_values(by=sample_name, ascending=0)
        sample = sample.T
        # define color map based on RNA class order in RNA class statistic data frame
        color_map = []    
        for key in sample.columns:
            color_map.append(color_map_collection[key]) 
        # plot statistics
        fig = sample.plot(kind='barh', stacked=True,color=color_map,legend=True,figsize=(14,2),width=0.3)
        fig.annotate(title_plot, xy=(28, 0.18), xytext=(35, 0.19),fontsize=19)
        fig.legend(loc='lower center',ncol=5,bbox_to_anchor=(0.5,-0.2),bbox_transform=plt.gcf().transFigure,frameon=False, prop={'size': 15})
        fig.axis('off')
        plt.savefig(title_plot.replace(" ", "") + '.pdf', dpi = 300,size=(2000,7000),bbox_inches='tight')
        plt.close('all')



parser = argparse.ArgumentParser(description="""plot RNA_class_stats""")
parser.add_argument("-i", "--input_files", metavar='<input_files>', default='*.csv', help="Path to RNA class mapping statistics tsv table")
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

# create dictionary with RNA classes and colors assigned to them 
colour_map = create_color_pallet(df_stats.columns,colors)

# plot RNA class statistics for each sample
df_stats = df_stats.T
for row in df_stats.columns:
    sample = pd.DataFrame(df_stats[row]).T
    plot_RNA_class_stats_each(sample,colour_map)

