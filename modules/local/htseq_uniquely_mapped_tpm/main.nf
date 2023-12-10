process HTSEQ_UNIQUELY_MAPPED_TPM {
    tag "htseq_uniquely_mapped_TPM"
    label 'process_high'

    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
    storeDir "${params.outdir}/HTSeq"
    
    input: 
    path input_quantification 
    val(host_attribute)
    path gff_host
    path gff_pathogen

    output:
    path "quantification_results_uniquely_mapped_NumReads_TPM.tsv", emit: tsv
    path "versions.yml"           , emit: versions

    script:
    """
    $workflow.projectDir/bin/calculate_TPM_HTSeq.R \\
    $input_quantification \\
    $host_attribute \\
    $gff_pathogen \\
    $gff_host

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}