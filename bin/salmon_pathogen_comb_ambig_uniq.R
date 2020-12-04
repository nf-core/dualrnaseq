#!/usr/bin/env Rscript

#-------------------------
#
# Description: 
# Combine Pathogen quantification summaries (including ambig and unique metrics) into one DF
#
# Created by R. Hayward
# Date: 1st September 2020
#
# Input args:   [1] = salmon/*/aux_info/*_pathogen_quant_ambig_uniq.sf
# Output file: "pathogen_quant_combined_ambig_uniq.tsv"
#
#-------------------------

args<-commandArgs(TRUE)

# Get file list
file_list <- list.files(path="salmon", pattern = "*pathogen_quant_ambig_uniq.sf", recursive = T, include.dirs = T, full.names = T)
#Get sample names
sample_names <- list.files(path="salmon", pattern = "*pathogen_quant_ambig_uniq.sf", recursive = T, include.dirs = T, full.names = T)
#just keep simple sample names
sample_namesv2 <- substring(sample_names, 8)
sample_namesv3 <- gsub('.{29}$','',sample_namesv2)

# Read all csv files in the folder and create a list of dataframes
ldf <- lapply(file_list , read.table, sep = '\t', header=T)

# Combine each dataframe in the list into a single dataframe
df.final <- do.call("cbind", ldf)

#Count no. of samples
sample_no <- length(sample_names)

#Loop through each sample within df and add modified sample names
n=1
i=7 #The number of cols from each sample
for (j in 1:sample_no) {
  colnames(df.final)[n:i] <- paste0(sample_namesv3[j],"_",colnames(df.final)[n:i])
  #increase counters
  n=n+7
  i=i+7
}

#Change the first three cols back to original names
colnames(df.final)[1:3] <- c("Name","Length","EffectiveLength")

#Rename cols with sample name_xxx
df.final_v2 <- df.final[,!(names(df.final)[1:length(colnames(df.final))] %in% grep("*_Name", names(df.final), value = TRUE))]
df.final_v2 <- df.final_v2[,!(names(df.final_v2)[1:length(colnames(df.final_v2))] %in% grep("*_Length", names(df.final_v2), value = TRUE))]
df.final_v2 <- df.final_v2[,!(names(df.final_v2)[1:length(colnames(df.final_v2))] %in% grep("*_EffectiveLength", names(df.final_v2), value = TRUE))]

#Save
write.table(df.final_v2, "pathogen_quant_combined_ambig_uniq.tsv", row.names = F, sep = "\t", quote = F)

