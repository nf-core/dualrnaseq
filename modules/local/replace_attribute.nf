// ONE FOR HOST_GENOME
// ONE FOR HOST_TRNA
process REPLACE_ATTRIBUTE_GFF_STAR_SALMON {
	tag "repl_attribute_host_tRNA_gff"

	label 'process_high'

	input:
		path(gff) 
		val(attribute_in)
		val(attribute_out)

	output:
		path "${outfile_name}"

	script:
		outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_${attribute_out}_attribute.gff3")
		"""
		$workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} $attribute_out $attribute_in
		"""
}
