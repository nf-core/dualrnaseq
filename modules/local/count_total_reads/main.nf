process COUNT_TOTAL_READS {
	tag "count_total_reads"
	publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
	storeDir "${params.outdir}/mapping_statistics"

	label 'process_high'

	input:
	path(fastq)
	output:
	path "total_raw_reads_fastq.tsv", emit:tsv

	script:
	"""
	$workflow.projectDir/bin/count_total_reads.sh $fastq >> total_raw_reads_fastq.tsv
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        version: 1.0.0
    END_VERSIONS
	"""
}