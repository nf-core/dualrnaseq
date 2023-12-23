process plot_star_mapping_stats {
    tag "plot_star_mapping_stats"
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR"

    label 'process_high'

    input:
    file(stats) from star_mapped_stats_to_plot

    output:
    file "*.tsv"
    file "*.pdf"

    script:
    """
    python $workflow.projectDir/bin/plot_mapping_stats_star.py -i $stats
    """
}