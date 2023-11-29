process SPLIT_QUATIFICATION_TABLES_HTSEQ {
	label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

	input: 
		path(quant_table)
		path(host_annotations)
		path(pathogen_annotations)

	output:
		path "host_quantification_uniquely_mapped_htseq.tsv", emit: host_quantification_stats_htseq
		path "pathogen_quantification_uniquely_mapped_htseq.tsv", emit: pathogen_quantification_stats_htseq
		env pathonen_tab, emit: pathonen_tab
		env host_tab, emit: host_tab

	script:
	"""
	bash $workflow.projectDir/bin/split_quant_tables.sh $quant_table $host_annotations $pathogen_annotations quantification_uniquely_mapped_htseq.tsv
	pathonen_tab=\$(if [ \$(cat pathogen_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
	host_tab=\$(if [ \$(cat host_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
	"""
	}