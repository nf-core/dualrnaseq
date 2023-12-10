process STAR_MAPPING_STATS {
    storeDir "${params.outdir}/mapping_statistics/STAR"
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    tag "star_mapping_stats"

    label 'process_high' 

    input:
    path total_raw_reads
    path total_processed_reads
    path uniquely_mapped_reads
    path multi_mapped_reads 
    path cross_mapped_reads

    output:
    file ('star_mapping_stats.tsv'), emit: tsv

    script:
    """
    python $workflow.projectDir/bin/mapping_stats.py -total_raw $total_raw_reads -total_processed $total_processed_reads -m_u $uniquely_mapped_reads -m_m $multi_mapped_reads -c_m $cross_mapped_reads -t star -o star_mapping_stats.tsv
    """
}
