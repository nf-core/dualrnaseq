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
    path("${name}"), emit: txt
    path "versions.yml"           , emit: versions

    shell: 
    name = sample_name + '_multi_mapped.txt'
    if (params.single_end){
    '''
    !{workflow.projectDir}/bin/count_multi_mapped_reads.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    '''
    } else {
    '''
    !{workflow.projectDir}/bin/count_multi_mapped_read_pairs.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    '''
    }
}