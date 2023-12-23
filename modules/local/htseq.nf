process HTSEQ_COUNT {
    tag "$sample_name.id"
    label 'process_high'

    conda "bioconda::htseq=2.0.2-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py38h7a2e8c7_0' :
        'quay.io/biocontainers/htseq:2.0.2--py38h7a2e8c7_0' }"

    input:
    tuple val(sample_name), path(st)
    tuple path(gff)
    val(quantifier)

    output:
	tuple val(sample_name), path("*_count.txt"), emit: txt
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
	def output_file = "${sample_name.id}_count.txt"
    """
	htseq-count  \\
        -n ${task.cpus}  \\
        -t $quantifier  \\
        -f bam  \\
        -r pos $st $gff  \\
        -i $params.host_gff_attribute  \\
        -s $params.stranded  \\
        --max-reads-in-buffer=${params.max_reads_in_buffer}  \\
        -a ${params.minaqual}  \\
        ${params.htseq_params}  \\
        $args  \\
        > $output_file

	sed -i '1{h;s/.*/'"${sample_name.id}"'/;G}' "$output_file"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htseq-count: \$( htseq-count --help | grep -i version | awk '{ print \$10 }' )
    END_VERSIONS
    """
}
