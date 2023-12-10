process REMOVE_CROSSMAPPED_READS {
    tag "$sample_name"
    label 'process_high'
    publishDir "${params.outdir}/STAR/multimapped_reads", mode: params.publish_dir_mode
    storeDir "${params.outdir}/STAR/multimapped_reads"


    input:
    tuple val(sample_name), path(alignment)
    path(host_reference)
    path(pathogen_reference)


    output:
    tuple val(sample_name), path("${bam_file_without_crossmapped}"), emit: alignment_multi_mapping_stats
    path("${cross_mapped_reads}"), emit: count_crossmapped_reads
    path "versions.yml"           , emit: versions

    script:
    bam_file_without_crossmapped = sample_name + "_no_crossmapped.bam"
    cross_mapped_reads = sample_name + "_cross_mapped_reads.txt"
    
    if (params.single_end){
        """
        $workflow.projectDir/bin/remove_crossmapped_reads_BAM.sh $alignment $workflow.projectDir/bin $host_reference $pathogen_reference $cross_mapped_reads $bam_file_without_crossmapped

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            version
        END_VERSIONS
        """
    } else {
        """
        $workflow.projectDir/bin/remove_crossmapped_read_pairs_BAM.sh $alignment $workflow.projectDir/bin $host_reference $pathogen_reference $cross_mapped_reads $bam_file_without_crossmapped

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            version
        END_VERSIONS
        """
    }
}