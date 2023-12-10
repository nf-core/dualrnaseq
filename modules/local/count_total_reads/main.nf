process COUNT_TOTAL_READS {
	tag "count_total_reads"
	publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
	storeDir "${params.outdir}/mapping_statistics"

	label 'process_high'

	input:
	file(fastq) from raw_read_count_file.collect()
	output:
	file "total_raw_reads_fastq.tsv" into to_collect_total_reads

	script:
	"""
	$workflow.projectDir/bin/count_total_reads.sh $fastq >> total_raw_reads_fastq.tsv
	"""
}