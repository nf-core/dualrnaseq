process COUNT_TOTAL_READS_PAIRS {
    tag "count_total_reads"
    publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics"

    label 'process_high'

    input:
    file(tsv) from to_collect_total_reads.collect()
    output:
    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_salmon
    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_salmon_alignment
    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_star
    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_star_for_salmon

    script:
    """
    $workflow.projectDir/bin/collect_total_raw_read_pairs.py -i $tsv
    """
}