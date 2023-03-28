process SALMON_SPLIT_TABLE {
    label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'nfcore/dualrnaseq:dev' }"

    input:
    tuple val(meta), path(quant)
    path transcript_fasta_pathogen
    path transcript_fasta_host

    output:
	tuple val(meta), path("host_quant.sf"),        emit: host
    tuple val(meta), path("pathogen_quant.sf"),    emit: pathogen

    script:
    """
    grep ">" ${transcript_fasta_pathogen} \
    | awk -F ">" '{ print \$2 }' \
    | awk 'NR==FNR{a[\$0]=\$0}NR>FNR{if(\$1==a[\$1])print \$0}' - ${quant} \
    > pathogen_quant

    awk 'NR==1 {print; exit}' ${quant} \
    | cat - pathogen_quant \
    > pathogen_quant.sf

    grep ">" ${transcript_fasta_host} \
    | awk -F ">" '{ print \$2 }' \
    | awk -F ' ' '{print \$1}' \
    | awk 'NR==FNR{a[\$0]=\$0}NR>FNR{if(\$1==a[\$1])print \$0}' - ${quant} \
    > host_quant

    awk 'NR==1 {print; exit}' ${quant} \
    | cat - host_quant \
    > host_quant.sf
    """
}
