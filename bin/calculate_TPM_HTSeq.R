#!/usr/bin/env Rscript

library(plyr)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
gene_attribute <- args[2]
table_htseq <- read.table(args[1], sep="\t", stringsAsFactors=F, header = T,row.names = gene_attribute) 
pathogen_gff <- import(args[3])

host_gff <- import(args[4])
pathogen_gff <- import(args[3])


extract_gene_length <- function(x,host_gff){
  #find annotations that contain a host gene_id (hosT_gene)
  h_gene <- host_gff[host_gff$gene_id == x]
  # Filter quant features
  quants <- h_gene[h_gene$type == "quant",]
  ## reduce - merge overlapping features, so that the same bases are covered by a reduced collection of features.
  g_length <- sum(width(reduce(quants)))
  return(g_length)
}


if(gene_attribute == 'gene_id'){
  #find pathogen genes
  common_pathogen = intersect(rownames(table_htseq),pathogen_gff$gene_id)
  pathogen_table_htseq <- table_htseq[common_pathogen,]
  pathogen_gff_match <- match(common_pathogen,pathogen_gff$gene_id)
  pathogen_table_htseq <- cbind(length = pathogen_gff@ranges@width[pathogen_gff_match],pathogen_table_htseq) 
  
  #find host genes
  common_host = intersect(rownames(table_htseq),host_gff$gene_id)
  #extract quantification results of host genes
  host_table_htseq <- table_htseq[common_host,]
  #extract gene length
  gene_length <- sapply(common_host, extract_gene_length, host_gff=host_gff,simplify = TRUE)
  host_table_htseq <- cbind(length = gene_length,host_table_htseq)
}



# splitting up isoforms as preparation for the next step
#tmp <- split(quants,as.character(quants$Parent))
# for each isoform, calculate the sum of all reduced exons
#Gene_length <- sum(width(reduce(tmp)))

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

names(quant_results)[1] <- gene_attribute
 
write.table(quant_results,file = "quantification_results_uniquely_mapped_NumReads_TPM.csv",sep = "\t", row.names = F, quote = F)

