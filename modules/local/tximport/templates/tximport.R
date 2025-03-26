#!/usr/bin/env Rscript

#-------------------------
#
# Description: 
# Script to calculate gene-level estimates of host quantification results of Salmon using tximport
#
# Created by B. Mika-Gospodorz
#
# Input args:  
#			    annotations_file: tsv annotation file created with extract_annotations_from_gff.py script
#				  prefix: sample name
# Output files: "*_host_quant_gene_level.sf"
#
#-------------------------



################################################
################################################
## Functions                                  ##
################################################
################################################

#' Parse out options from a string without recourse to optparse


parse_args <- function(x){
  args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
  args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))
  
  # Ensure the option vectors are length 2 (key/ value) to catch empty ones
  args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})
  
  parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
  parsed_args[! is.na(parsed_args)]
}


################################################
################################################
###   Load libraries                          ##
################################################
################################################

library(tximport)

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
#####################################################
# Set defaults 

opt <- list( 
  annotations_file = '$annotations',
  prefix = '$prefix',
  host_quant = '$host_quant'
)

opt_types <- lapply(opt, class)


# Apply parameter overrides
## to overwrite extra parames
args_opt <- parse_args('$task.ext.args')
for ( ao in names(args_opt)){
  if (! ao %in% names(opt)){
    stop(paste("Invalid option:", ao))
  }else{
    
    # Preserve classes from defaults where possible
    if (! is.null(opt[[ao]])){
      args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])
    }
    opt[[ao]] <- args_opt[[ao]]
  }
}

################################################
################################################
##      Tximport                              ##
################################################
################################################


annotations <- read.table(opt\$annotations_file, header = TRUE)
# extract 'transcript_id' and 'gene_id' column from annotation table
tx2gene <-annotations[,c('transcript_id','gene_id')]

# run tximport
txi = tximport(opt\$host_quant, type = "salmon",tx2gene = tx2gene, dropInfReps=TRUE)


# extract TPM, length and counts 
TPMs <- txi\$abundance
length <- txi\$length
counts <- txi\$counts

colnames(TPMs) <- c("TPM")
colnames(length) <- c("EffectiveLength")
colnames(counts) <- c("NumReads")


# combine results and save the table
gene_results <- cbind('Name'=rownames(txi\$abundance), TPMs, length, counts)
write.table(gene_results,file = paste(opt\$prefix,"_host_quant_gene_level.sf",sep=''),sep = "\t", row.names = F, quote = F)


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
tximport.version <- as.character(packageVersion('tximport'))
writeLines(
  c(
    '"${task.process}":',
    paste('    r-base:', r.version),
    paste('    tximport:', tximport.version)
  ),
  'versions.yml')