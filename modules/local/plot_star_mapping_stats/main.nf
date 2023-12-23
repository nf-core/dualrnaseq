process PLOT_STAR_MAPPING_STATS {
    tag "plot_star_mapping_stats"
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR"

    label 'process_high'

    input:
    path(stats)

    output:
    path "*.tsv"
    path "*.pdf"
    path "versions.yml"          , emit: versions

    script:
    """
    python $workflow.projectDir/bin/plot_mapping_stats_star.py -i $stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}