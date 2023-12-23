process PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH{
    publishDir "${params.outdir}/mapping_statistics/${tool}/RNA_classes_pathogen", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/${tool}/RNA_classes_pathogen"
    tag "plot_rna_stats_htseq_pathogen"

    label 'process_high'

    input:
    file stats_table
    val(tool)

    output:
    path "*.pdf"
    path "versions.yml" , emit: versions

    script:
    """
    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py \\
    -i $stats_table
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}

		