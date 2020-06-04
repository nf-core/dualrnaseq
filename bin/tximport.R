#!/usr/bin/env Rscript
library(tximport)

args = commandArgs(trailingOnly=TRUE)

salmon_path = args[1]
annotations <- read.table(args[2], header = TRUE)


tx2gene <-annotations[,c('transcript_id','gene_id')]

files = list.files(salmon_path, pattern = "host_quant.sf", recursive = T, full.names = T)
names = basename(dirname(files))
names(files) = names

txi = tximport(files, type = "salmon",tx2gene = tx2gene, dropInfReps=TRUE)

TPM <- cbind(rownames(txi$abundance),txi$abundance)
write.table(TPM,file = "host_gene_TPM.csv",sep = "\t", row.names = F)

length <- cbind(rownames(txi$length),txi$length)
write.table(length,file = "host_gene_length.csv",sep = "\t", row.names = F)

counts <- cbind(rownames(txi$counts),txi$counts)
write.table(counts,file = "host_gene_counts.csv",sep = "\t", row.names = F)
