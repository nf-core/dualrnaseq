process COUNT_READS {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.txt'), emit: read_counts

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    count_1=\$(cat ${reads[0]} | echo \$((\$(wc -l) / 4)))
    count_2=\$(cat ${reads[1]} | echo \$((\$(wc -l) / 4)))
    echo "${reads[0].baseName}\t\$count_1" > ${reads[0].baseName}_counts.txt
    echo "${reads[1].baseName}\t\$count_2" > ${reads[1].baseName}_counts.txt
    """
}
