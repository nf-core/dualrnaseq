process HTSEQ_QUANTIFICATION_STATS_UNIQUELY_MAPPED {
    storeDir "${params.outdir}/mapping_statistics/HTSeq"
    publishDir "${params.outdir}/mapping_statistics/HTSeq", mode: params.publish_dir_mode
    tag "quantification_stats_htseq"

    label 'process_high'

    input:
    path quant_table_host
    path quant_table_pathogen
    val attribute
    path star_stats
    

    output:
    file ('htseq_uniquely_mapped_reads_stats.tsv'), emit: htseq_mapped_stats_to_plot
    path "versions.yml"           , emit: versions

    script:
    """
    python $workflow.projectDir/bin/mapping_stats.py \\
     -q_p $quant_table_pathogen \\
     -q_h $quant_table_host \\
     -a $attribute \\
     -star $star_stats \\
     -t htseq \\
     -o htseq_uniquely_mapped_reads_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS

    """
}