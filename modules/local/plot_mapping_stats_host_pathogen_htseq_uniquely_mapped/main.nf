process PLOT_MAPPING_STATS_HOST_PATHOGEN_HTSEQ_UNIQUELY_MAPPED{
    tag "$name2"
    publishDir "${params.outdir}/mapping_statistics/HTSeq", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/HTSeq"

    label 'process_high'

    input:
    path(stats)

    output:
    path "*.tsv"
    path "*.pdf"
    path "versions.yml"           , emit: versions

    script:
    """
    python $workflow.projectDir/bin/plot_mapping_stats_htseq.py -i $stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}