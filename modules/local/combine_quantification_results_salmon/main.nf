process COMBINE_QUANTIFICATION_RESULTS_SALMON {
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
    python $workflow.projectDir/bin/collect_quantification_data_salmon.py \
        -i $input_quantification \
        -a $params.gene_attribute_gff_to_create_transcriptome_host \
        -org $organism
    """
}