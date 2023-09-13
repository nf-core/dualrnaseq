process COMBINE_QUANTIFICATION_RESULTS_HTS {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

    input: 
	    path input_quantification
    output:
	    path "quantification_results_htseq.tsv", emit: combined_quant_data
    script:
    """
    python $workflow.projectDir/bin/collect_quantification_data_hts.py \
        -i $input_quantification \
        -a $params.gene_attribute_gff_to_create_transcriptome_host        
    """
}

