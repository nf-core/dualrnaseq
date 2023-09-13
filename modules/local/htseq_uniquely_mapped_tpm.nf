process HTSEQ_UNIQUELY_MAPPED_TPM {
	label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

	input: 
		tuple val(meta), path(input_quantification)
		val(host_attribute)
		path(gff_host)
		path(gff_pathogen)

	output:
		path "quantification_results_uniquely_mapped_NumReads_TPM.tsv", emit: split_table_htseq_host
		path "quantification_results_uniquely_mapped_NumReads_TPM.tsv", emit: split_table_htseq_pathogen

	script:
	"""
	Rscript $workflow.projectDir/bin/calculate_TPM_HTSeq.R $input_quantification $host_attribute $gff_pathogen $gff_host
	"""
	}