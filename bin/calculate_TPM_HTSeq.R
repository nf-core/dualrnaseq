#!/usr/bin/env Rscript

#-------------------------
#
# Description: 
# Script to count TPM values for HTSeq quantification results
#
# Created by B. Mika-Gospodorz
# Date: 29th May 2020
#
# Input args:   [1] = quantification_results_uniquely_mapped.tsv - HTSeq quantification results
#			    [2] = host gene attribute used for quantification, e.g. 'gene_id'
#				[3] =  pathogen gff file with gene features (from 3rd col) replaced with 'quant'
#				[4] =  host gff file with gene features (from 3rd col) replaced with 'quant'
# Output file: quantification_results_uniquely_mapped_NumReads_TPM.tsv
#
#-------------------------

library(plyr)
library(rtracklayer)


args = commandArgs(trailingOnly=TRUE)
# read quantification table
table_htseq <- read.table(args[1], sep="\t", stringsAsFactors=F, header = T,row.names = gene_attribute) 
gene_attribute <- args[2]
pathogen_gff <- import(args[3])
host_gff <- import(args[4])


# function to extract gene length from gff file
extract_transcript_length <- function(x,host_gff,gene_attribute){
  #find annotations that contain host transcripts
  h_tr <- host_gff[mcols(host_gff)[gene_attribute][[1]]==x]
  # Filter quant features
  quants <- h_tr[h_tr$type == "quant",]
  # merge overlapping features, so that the same bases are covered by reduced collection of features
  t_length <- sum(width(reduce(quants)))
  return(t_length)
}

######################### extract gene length #######################################

### pathogen
names_gff_pathogen <- mcols(pathogen_gff)[gene_attribute][[1]]  # extract gene names from gff file
common_pathogen = intersect(rownames(table_htseq), names_gff_pathogen) # find positions of pathogen genes in quantification table
pathogen_table_htseq <- table_htseq[common_pathogen,] # extract pathogen quantification results
pathogen_gff_match <- match(common_pathogen,names_gff_pathogen) # find positions of corresponding genes in gff file
pathogen_table_htseq <- cbind(length = pathogen_gff@ranges@width[pathogen_gff_match],pathogen_table_htseq) # extract gene length from gff file and combine it with quantification table


### host
names_gff_host <- mcols(host_gff)[gene_attribute][[1]] # extract gene names from gff file
lack_of_attribute <- which(is.na(names_gff_host)) #find positions without gene_attribute (eg.genes that don't have transcript_id attribute)
if (length(lack_of_attribute)> 0) {
  host_gff <- host_gff[-lack_of_attribute] #remove positions without gene_attribute 
  names_gff_host <- names_gff_host[-lack_of_attribute] # extract gene names that contain gene_attribute
}
common_host = intersect(rownames(table_htseq),names_gff_host) # find positions of host genes in quantification table
host_table_htseq <- table_htseq[common_host,] #extract quantification results of host genes 
transcript_length <- sapply(common_host, extract_transcript_length, host_gff=host_gff,gene_attribute = gene_attribute,simplify = TRUE) #extract gene length
host_table_htseq <- cbind(length = transcript_length,host_table_htseq) # add gene length into quantification table


# combine host and pathogen tables
htseq_table_with_gene_length = rbind(pathogen_table_htseq,host_table_htseq)


######################### calculate TPM values ####################################### 

# function to calculate TPMs
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  return(rate / sum(rate) * 1e6)
}

# function to add _TPM suffix
rename_add_TPM <- function(x) {
  paste(x,"_TPM",sep='')
}

# function to add _NumReads suffix
rename_add_NumReads <- function(x) {
  paste(x,"_NumReads",sep='')
}


# calculate TPMs for each gene in each sample
TPMs <- apply(htseq_table_with_gene_length[2:dim(htseq_table_with_gene_length)[2]], 2, function(x) tpm(x, htseq_table_with_gene_length$length))
colnames(TPMs) <- colnames(htseq_table_with_gene_length[,-1])
# add TPM suffix to column names with TPM values
colnames(TPMs) <-sapply(colnames(TPMs),function(x) rename_add_TPM(x))
# add NumReads suffix to column names with NumReads values
colnames_NR <- sapply(colnames(htseq_table_with_gene_length[,-1]),function(x) rename_add_NumReads(x))
# add 'length' name for column which stores gene length
colnames(htseq_table_with_gene_length) <- (c('length',colnames_NR))

# combine results
quant_results <- cbind(name = rownames(TPMs), htseq_table_with_gene_length,TPMs)
names(quant_results)[1] <- gene_attribute

# save results
write.table(quant_results,file = "quantification_results_uniquely_mapped_NumReads_TPM.tsv",sep = "\t", row.names = F, quote = F)

