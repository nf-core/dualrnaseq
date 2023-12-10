process COLLECT_STATS_STAR_UNIQUELY_MAPPED {
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR"
    tag "collect_uniq_mapped_reads_STAR"

    label 'process_high' 
    
    input: 
    path stats

    output:
    path "uniquely_mapped_reads_star.tsv", emit: tsv
    path "versions.yml"           , emit: versions

    script:
    """
    python $workflow.projectDir/bin/combine_tables.py \\
    -i $stats \\
    -o uniquely_mapped_reads_star.tsv \\
    -s uniquely_mapped_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    """
}