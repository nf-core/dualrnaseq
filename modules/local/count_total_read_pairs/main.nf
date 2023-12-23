process COUNT_TOTAL_READS_PAIRS {
    tag "count_total_reads"
    publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics"

    label 'process_high'

    input:
    path(tsv)
    output:
    path "total_raw_read_pairs_fastq.tsv", emit:tsv

    script:
    """
    $workflow.projectDir/bin/collect_total_raw_read_pairs.py -i $tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}