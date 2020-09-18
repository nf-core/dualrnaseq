#!/usr/bin/env Rscript

#-------------------------
#
# Description: 
# Script to combine quant.sf and ambig_info.tsv after running Salmon
# Also generates gene level summaries, which are traditionally performed with tximport,
# but cannot currently do this when examining the ambig_info file
#
# Created by R. Hayward
# Date: 1st September 2020
#
# Input args:   [1] = "salmon/*/quant.sf"
#			    [2] = "salmon/*/aux_info/ambig_info.tsv"
#				[3] = "$sample_name"
#				[4] = "$annotations"
# Output files: "sample_ID_host_quant_ambig_uniq.sf"
#				"sample_ID_pathogen_quant_ambig_uniq.sf"
#				"sample_ID_host_quant_ambig_uniq_gene_level.sf"
#
#-------------------------


args<-commandArgs(TRUE)


#Load in the annotations file
annotations <- read.table(args[4], header = TRUE)

tx2gene <-annotations[,c('transcript_id','gene_id')]

#Read in quant.sf
try(quant_df <- read.delim(args[1],header = T, stringsAsFactors = F))
#Read in ambig_info.tsv
try(ambig_df <- read.table(args[2],header = T, stringsAsFactors = F))

#Join tables together
input_df <- cbind(quant_df, ambig_df)

#Extract out Host and Pathogen annotations
host_df <- input_df[input_df$Name %in% tx2gene$transcript_id,]
pathogen_df <- input_df[!input_df$Name %in% tx2gene$transcript_id,]

#Save host and pathogen quant files
write.table(host_df, paste0(args[3],"_host_quant_ambig_uniq.sf"), sep = "\t",col.names = T, row.names = F, quote = F)
write.table(pathogen_df, paste0(args[3],"_pathogen_quant_ambig_uniq.sf"), sep = "\t",col.names = T, row.names = F, quote = F)

#Separate out each column from host data
#TPM
abundanceMatTx <- matrix(data=host_df$TPM, ncol = 1, dimnames = list(host_df$Name,"TPM"))
#No. of reads
countsMatTx <- matrix(data=host_df$NumReads, ncol = 1, dimnames = list(host_df$Name,"NumReads"))
#Effecive length
lengthMatTx <- matrix(data=host_df$EffectiveLength, ncol = 1, dimnames = list(host_df$Name,"Length"))
#Ambigious counts
ambigMatTx <- matrix(data=host_df$AmbigCount, ncol = 1, dimnames = list(host_df$Name,"AmbigCount"))
#Unique counts
uniqMatTx <- matrix(data=host_df$UniqueCount, ncol = 1, dimnames = list(host_df$Name,"UniqCount"))
#Transcript names
txId <- host_df$Name

#Just keep transcripts identified in sample
tx2gene <- tx2gene[tx2gene$transcript_id %in% txId, ]
#tx2gene$gene_id <- droplevels(tx2gene$gene_id)

#Match gene ID
geneId <- tx2gene$gene_id[match(txId, tx2gene$transcript_id)]

#Summarise each of the transcripts back to their genes
abundanceMat <- rowsum(abundanceMatTx, geneId)
countsMat <- rowsum(countsMatTx, geneId)
lengthMat <- rowsum(lengthMatTx, geneId)
uniqMat <- rowsum(uniqMatTx, geneId)
ambigMat <- rowsum(ambigMatTx, geneId)

#Join gene df together
host_gene_mat <- cbind(abundanceMat,lengthMat,countsMat,ambigMat,uniqMat)

#Convert to df
host_gene_df <- as.data.frame(host_gene_mat)

#Add sample name in front of each col heading
colnames(host_gene_df) <- paste(args[3], colnames(host_gene_df), sep = "_")

#add Name col and re-arrange
host_gene_df$Name <- row.names(host_gene_df)
host_gene_df <- host_gene_df[,c(6,1:5)]

#Save down host data
write.table(host_gene_df, file = paste(args[3],"_host_quant_ambig_uniq_gene_level.sf",sep=''),sep = "\t", row.names = F, quote = F)
