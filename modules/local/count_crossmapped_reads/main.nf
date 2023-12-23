process COUNT_CROSSMAPPED_READS {
    tag "count_crossmapped_reads"
    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR"

    label 'process_high'

    input:
    path(cross_mapped_reads)

    output:
    path "cross_mapped_reads_sum.txt", emit: txt
    path "versions.yml"           , emit: versions
    
    script:
    """
    $workflow.projectDir/bin/count_cross_mapped_reads.sh $cross_mapped_reads
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
}