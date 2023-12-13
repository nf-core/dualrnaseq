process COUNT_TOTAL_READS {
	tag "count_total_reads"
	publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
	storeDir "${params.outdir}/mapping_statistics"

	label 'process_high'

	input:
	path(fastq)
	output:
	path "total_raw_reads_fastq.tsv"

	script:
	"""
	$workflow.projectDir/bin/count_total_reads.sh $fastq >> total_raw_reads_fastq.tsv
	"""
}