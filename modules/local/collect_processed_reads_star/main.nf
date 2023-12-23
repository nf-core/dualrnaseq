process COLLECT_PROCESSED_READS_STAR {
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR"
    tag "collect_processed_reads_STAR"

    label 'process_high' 
    
    input: 
    path(process_reads)

    output:
    path("processed_reads_star.tsv")
    path "versions.yml"           , emit: versions

    script:
    """
    cat $process_reads > processed_reads_star.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}