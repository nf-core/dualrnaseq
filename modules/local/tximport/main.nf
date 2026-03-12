process TXIMPORT {
    tag "${meta.id}"
    label 'process_single'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'nfcore/dualrnaseq:dev'
        : 'nfcore/dualrnaseq:dev'}"

    input:
    tuple val(meta), file(host_quant)
    file annotations

    output:
    tuple val(meta), file("${meta.id}_host_quant_gene_level.sf"), emit: salmon_gene_level
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    template('tximport.R')

    stub:
    """
    touch ${meta.id}_host_quant_gene_level.sf

    Rscript -e "r.version <- strsplit(version[['version.string']], ' ')[[1]][3];
                tximport.version <- as.character(packageVersion('tximport'));
                writeLines(
                    c(
                        '"${task.process}":',
                        paste('    r-base:', r.version),
                        paste('    tximport:', tximport.version)
                    ),
                'versions.yml')"

    """
}
