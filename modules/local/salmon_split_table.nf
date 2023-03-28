process SALMON_SPLIT_TABLE {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'nfcore/dualrnaseq:dev' }"

    input:
	    path input_quantification
        val organism
    output:
	    path "combined_${organism}.tsv", emit: combined_quant_data
    script:
    """
    $workflow.projectDir/bin/split_quant_tables_salmon.sh \
    $transcriptome_pathogen $transcriptome_host \
    salmon/*/quant.sf "quant.sf"
    """
}
