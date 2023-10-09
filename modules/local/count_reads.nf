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

    if (reads instanceof Path || reads.size() == 1) {
        // Code to process a single file
        """
        count_1=\$((\$(zcat ${reads[0]} | wc -l) / 4))
        echo "${reads[0].baseName}\t\$count_1" > ${reads[0].baseName}_counts.txt
        """
    } else if (reads instanceof Path || reads.size() == 2) {
        // Code to process two files
        """
        count_1=\$((\$(zcat ${reads[0]} | wc -l) / 4))
        count_2=\$((\$(zcat ${reads[1]} | wc -l) / 4))
        echo "${reads[0].baseName}\t\$count_1" > ${reads[0].baseName}_counts.txt
        echo "${reads[1].baseName}\t\$count_2" > ${reads[1].baseName}_counts.txt
        """
    } else {
        """
        echo "Error: Input 'reads' should contain either one or two files."
        exit 1
        """
    }
}
