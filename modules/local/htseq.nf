process HTSEQ {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::htseq=2.0.2-0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htseq:2.0.2--py38h7a2e8c7_0' :
        'quay.io/biocontainers/htseq:2.0.2--py38h7a2e8c7_0' }"

    input:
	tuple val(meta), path(st)
    path(gff)
	val(host_attribute)
	val(quantifier)
    val(stranded)
    
    output:
	tuple val(meta), path("*_count.txt"), emit: results
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
	def output_file = meta.id + "_count.txt"
    """
	htseq-count  \\
        -n ${task.cpus}  \\
        -t $quantifier  \\
        -f bam  \\
        -r pos $st $gff  \\
        -i $host_attribute  \\
        -s $stranded  \\
        --max-reads-in-buffer=${params.max_reads_in_buffer}  \\
        -a ${params.minaqual}  \\
        $args  \\
        > $output_file

	sed -i '1{h;s/.*/'"$meta.id"'/;G}' "$output_file"

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    htseq-count: \$( htseq-count --help | grep -i version | tail -n 1 | cut -d' ' -f2 )
END_VERSIONS
    """
}
