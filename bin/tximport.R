#!/usr/bin/env Rscript

#-------------------------
#
# Description: 
# Script to calculate gene-level estimates of host quantification results of Salmon using tximport
#
# Created by B. Mika-Gospodorz
# Date: 29th May 2020
#
# Input args:   [1] =  directory where salmon results are located
#			    [2] = tsv annotation file created with extract_annotations_from_gff.py script
#				[3] = sample name
# Output files: "*_host_quant_gene_level.sf"
#
#-------------------------

library(tximport)

args = commandArgs(trailingOnly=TRUE)
salmon_path = args[1]
annotations <- read.table(args[2], header = TRUE)


# extract 'transcript_id' and 'gene_id' column from annotation table
tx2gene <-annotations[,c('transcript_id','gene_id')]

# list of quantification files
files = list.files(salmon_path, pattern = "host_quant.sf", recursive = T, full.names = T)
names = basename(dirname(files))
names(files) = names

# run tximport
txi = tximport(files, type = "salmon",tx2gene = tx2gene, dropInfReps=TRUE)

# extract TPM, length and counts 
TPMs <- txi$abundance
length <- txi$length
counts <- txi$counts

# rename colnames 
rename_add_TPM <- function(x) {
  paste(x,"_TPM",sep='')
}

rename_add_Length <- function(x) {
  paste(x,"_Length",sep='')
}

rename_add_NumReads <- function(x) {
  paste(x,"_NumReads",sep='')
}
colnames(TPMs) <-sapply(colnames(TPMs),function(x) rename_add_TPM(x))
colnames(length) <-sapply(colnames(length),function(x) rename_add_Length(x))
colnames(counts) <-sapply(colnames(counts),function(x) rename_add_NumReads(x))

# combine results and save the table
gene_results <- cbind('Name'=rownames(txi$abundance), TPMs, length, counts)
write.table(gene_results,file = paste(args[3],"_host_quant_gene_level.sf",sep=''),sep = "\t", row.names = F, quote = F)
