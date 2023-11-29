process COMBINE_ANNOTATIONS_QUANT {
	label 'process_high'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

	input: 
		path(quantification_table)
		path(annotation_table)
		val(attribute)

	output:
		path "pathogen_combined_quant_annotations.tsv"

	script:
	"""
	python3 $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org pathogen
	"""
	}