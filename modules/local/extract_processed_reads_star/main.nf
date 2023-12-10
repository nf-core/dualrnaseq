process EXTRACT_PROCESSED_READS_STAR {
    publishDir "${params.outdir}/mapping_statistics/STAR/processed_reads", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR/processed_reads"
    tag "extract_processed_reads_STAR"

    label 'process_high'
    
    input: 
    tuple val(sample_name), path (Log_final_out)

    output:
    path "${sample_name}.txt", emit: txt
    path "versions.yml"           , emit: versions

    script:
    """
    $workflow.projectDir/bin/extract_processed_reads.sh $Log_final_out $sample_name star
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}