process UNIQUE_MAPPING_STATS_STAR {

    label 'process_high'
    tag "$sample_name"
    publishDir "${params.outdir}/mapping_statistics/STAR/uniquely_mapped", mode: params.publish_dir_mode
    storeDir "${params.outdir}/mapping_statistics/STAR/uniquely_mapped"

    input:
    tuple val(sample_name), file(alignment)
    path(host_reference_names)
    path(pathogen_reference_names)

    output:
    path("${name}"), emit: txt
    path "versions.yml"           , emit: versions
    
    shell: 
    name = sample_name + '_uniquely_mapped.txt'
    if (params.single_end){
    '''
    !{workflow.projectDir}/bin/count_uniquely_mapped_reads.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    '''
    } else {
    '''
    !{workflow.projectDir}/bin/count_uniquely_mapped_read_pairs.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version
    END_VERSIONS
    '''
    }
}