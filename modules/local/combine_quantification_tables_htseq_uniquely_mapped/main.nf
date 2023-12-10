process COMBINE_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED {
    tag "comb_quants_htseq_uniq_mapped"
    label 'process_high'
    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
    storeDir "${params.outdir}/HTSeq"
    

    input: 
    path input_quantification
    val(host_attribute)

    output:
    path "quantification_results_*.tsv", emit: tsv
    path "versions.yml"           , emit: versions

    script:
    """
    python $workflow.projectDir/bin/collect_quantification_data.py \\
    -i $input_quantification \\
    -q htseq \\
    -a $host_attribute 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}