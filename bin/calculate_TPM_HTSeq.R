#!/usr/bin/env Rscript

library(plyr)
library(rtracklayer)

args = commandArgs(trailingOnly=TRUE)
gene_attribute <- args[2]
table_htseq <- read.table(args[1], sep="\t", stringsAsFactors=F, header = T,row.names = gene_attribute) 
pathogen_gff <- import(args[3])

host_gff <- import(args[4])
pathogen_gff <- import(args[3])

#pathogen_gff <- import('~/Desktop/PhD_studies/Dual_RNA_seq/results/results_Hela_Salmonella_simulations_SingleEnd_75_ID/references/Salmonella_combined_BMG3_quant_feature_new_attribute.gff3')
#table_htseq <- read.table('~/Desktop/PhD_studies/Dual_RNA_seq/results/results_Hela_Salmonella_simulations_SingleEnd_75_ID/HTSeq/uniquely_mapped/quantification_results_uniquely_mapped.csv', sep="\t", stringsAsFactors=F, header = T,row.names = 'ID')
#host_gff <- import('~/Desktop/PhD_studies/Dual_RNA_seq/results/results_Hela_Salmonella_simulations_SingleEnd_75_ID/references/gencode.v33.chr_patch_hapl_scaff.annotation_with_tRNA_quant_feature.gff3')

extract_transcript_length <- function(x,host_gff,gene_attribute){
  #find annotations that contain host transcripts
  h_tr <- host_gff[mcols(host_gff)[gene_attribute][[1]]==x]
  # Filter quant features
  quants <- h_tr[h_tr$type == "quant",]
  ## reduce - merge overlapping features, so that the same bases are covered by a reduced collection of features.
  t_length <- sum(width(reduce(quants)))
  return(t_length)
}



#pathogen
#find pathogen genes
names_gff_pathogen <- mcols(pathogen_gff)[gene_attribute][[1]]
common_pathogen = intersect(rownames(table_htseq), names_gff_pathogen)
# `$`(pathogen_gff , 'ID')
pathogen_table_htseq <- table_htseq[common_pathogen,]
pathogen_gff_match <- match(common_pathogen,names_gff_pathogen)
pathogen_table_htseq <- cbind(length = pathogen_gff@ranges@width[pathogen_gff_match],pathogen_table_htseq) 
  
#host
#remove positions without the gene_attribute (eg.genes that don't have transcript_id attribute)
names_gff_host <- mcols(host_gff)[gene_attribute][[1]]
lack_of_attribute <- which(is.na(names_gff_host))
if (length(lack_of_attribute)> 0) {
  host_gff <- host_gff[-lack_of_attribute]
  names_gff_host <- names_gff_host[-lack_of_attribute]
}

#find host transcripts
common_host = intersect(rownames(table_htseq),names_gff_host)
#extract quantification results of host genes
host_table_htseq <- table_htseq[common_host,]

#extract gene length
transcript_length <- sapply(common_host, extract_transcript_length, host_gff=host_gff,gene_attribute = gene_attribute,simplify = TRUE)
host_table_htseq <- cbind(length = transcript_length,host_table_htseq)




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

