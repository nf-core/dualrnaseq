#!/usr/bin/env Rscript

library(plyr)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
gene_attribute <- args[2]
table_htseq <- read.table(args[1], sep="\t", stringsAsFactors=F, header = T,row.names = gene_attribute) 

pathogen_gff <- import(args[3])

host_gff <- import(args[4])


if(gene_attribute == 'gene_id'){
  common_pathogen = intersect(rownames(table_htseq),pathogen_gff$gene_id)
  pathogen_table_htseq <- table_htseq[common_pathogen,]
  pathogen_gff_match <- match(common_pathogen,pathogen_gff$gene_id)
  pathogen_table_htseq <- cbind(length = pathogen_gff@ranges@width[pathogen_gff_match],pathogen_table_htseq) 
  
  common_host = intersect(rownames(table_htseq),host_gff$gene_id)
  host_table_htseq <- table_htseq[common_host,]
  host_gff_match <- match(common_host,host_gff$gene_id)
  host_table_htseq <- cbind(length = host_gff@ranges@width[host_gff_match],host_table_htseq) 
}

htseq_table_with_gene_length = rbind(pathogen_table_htseq,host_table_htseq)

#tpm calculation
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  return(rate / sum(rate) * 1e6)
}

TPMs <- apply(htseq_table_with_gene_length[2:dim(htseq_table_with_gene_length)[2]], 2, function(x) tpm(x, htseq_table_with_gene_length$length))
colnames(TPMs) <- colnames(htseq_table_with_gene_length[,-1])

rename_add_TPM <- function(x) {
  paste(x,"_TPM",sep='')
}
colnames(TPMs) <-sapply(colnames(TPMs),function(x) rename_add_TPM(x))

rename_add_NumReads <- function(x) {
  paste(x,"_NumReads",sep='')
}

colnames_NR <- sapply(colnames(htseq_table_with_gene_length[,-1]),function(x) rename_add_NumReads(x))

colnames(htseq_table_with_gene_length) <- (c('length',colnames_NR))

quant_results <- cbind(name = rownames(TPMs), htseq_table_with_gene_length,TPMs)

rename(quant_results, c("name" = gene_attribute))
 
write.table(quant_results,file = "quantification_results_uniquely_mapped_NumReads_TPM.csv",sep = "\t", row.names = F, quote = F)

