process HTSEQ_COUNT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::htseq=2.0.2-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py38h7a2e8c7_0' :
        'quay.io/biocontainers/htseq:2.0.2--py38h7a2e8c7_0' }"

    // Note:
    // creating separate module here as the nf-core one asks fo a bam index, which isnt required.
    // also the naming convensions of input and outputs could be clearer.

    input:
    tuple val(meta), path(bam)
    path(gff)

    output:
    tuple val(meta), path("*_counts.txt"), emit: counts
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def output_file = meta.id + "_counts.txt"
    """
	htseq-count \\
        ${args} \\
        ${bam} \\
        ${gff} \\
        > ${output_file}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq-count: \$( htseq-count --help | grep -i version | tail -n 1 | cut -d' ' -f2 )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq-count: \$( htseq-count --help | grep -i version | tail -n 1 | cut -d' ' -f2 )
    END_VERSIONS
    """
}
