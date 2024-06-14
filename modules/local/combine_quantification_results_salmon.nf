process COMBINE_QUANTIFICATION_RESULTS_SALMON {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker.io/nfcore/dualrnaseq:dev' :
        'docker.io/nfcore/dualrnaseq:dev' }"

    input:
	    path input_quantification
        val organism
    output:
	    path "combined_${organism}.tsv", emit: combined_quant_data

    script:
    def args = task.ext.args ?: ''
    """
    python $workflow.projectDir/bin/collect_quantification_data_salmon.py \
        -i $input_quantification \
        -a $params.gene_attribute_gff_to_create_transcriptome_host \
        -org $organism
    """
}
