process COUNT_TOTAL_READS_PAIRS {
    tag "count_total_reads"
    publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics"

    label 'process_high'

    input:
    path(tsv)
    output:
    path "total_raw_read_pairs_fastq.tsv"

    script:
    """
    $workflow.projectDir/bin/collect_total_raw_read_pairs.py -i $tsv
    """
}