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

TPMs <- txi$abundance
length <- txi$length
counts <- txi$counts

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

gene_results <- cbind(rownames(txi$abundance), TPMs, length, counts)

write.table(gene_results,file = "host_quantification_gene_level.csv",sep = "\t", row.names = F, quote = F)
