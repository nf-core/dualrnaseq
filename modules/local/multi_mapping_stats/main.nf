process MULTI_MAPPING_STATS {
    tag "$sample_name"
    publishDir "${params.outdir}/mapping_statistics/STAR/multi_mapped", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR/multi_mapped"

    label 'process_high'

    input:
    tuple val(sample_name), path(alignment)
    path(host_reference_names)
    path(pathogen_reference_names)

    output:
    path("*_multi_mapped.txt"), emit: txt
    path "versions.yml"           , emit: versions

    script: 
    def name = "${sample_name.id}_multi_mapped.txt"
    if (params.single_end){
    """
    $workflow.projectDir/bin/count_multi_mapped_reads.sh ${alignment} ${host_reference_names} ${pathogen_reference_names} ${sample_name.id} ${name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
    } else {
    """
    $workflow.projectDir/bin/count_multi_mapped_read_pairs.sh ${alignment} ${host_reference_names} ${pathogen_reference_names} ${sample_name.id} ${name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
    """
    }
}