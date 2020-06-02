library(plyr)

args = commandArgs(trailingOnly=TRUE)
gene_attribute <- args[2]
table_htseq <- read.table(args[1], sep="\t", stringsAsFactors=F, header = T,row.names = gene_attribute) 

pathogen_gff <- rtracklayer::import(args[3]'../results/references/Salmonella_combined_BMG3_quant_feature_new_attribute.gff3')

host_gff <- rtracklayer::import(args[4]'../results/references/gencode.v33.chr_patch_hapl_scaff.annotation_with_tRNA_quant_feature.gff3')


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
colnames(TPMs) <- colnames(htseq_table_with_gene_length)

write.table(TPMs,file = "HTSeq_TPM.csv",sep = "\t", row.names = F)

write.table(htseq_table_with_gene_length,file = "HTSeq_quantification_with_gene_length.csv",sep = "\t", row.names = F)
