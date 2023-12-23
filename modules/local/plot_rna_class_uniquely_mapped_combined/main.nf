process PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED {
    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host"
    tag "plt_rna_stats_htseq_host_all"

    label 'process_high'

    input:
    path stats_table
    val (organism)

    output:
    path "RNA_class_stats_combined_host.pdf"
    path "versions.yml" , emit: versions

    script:
    """
    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py \\
    -i $stats_table \\
    -org $organism

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}
