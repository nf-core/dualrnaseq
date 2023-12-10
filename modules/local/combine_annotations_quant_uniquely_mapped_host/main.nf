process COMBINE_ANNOTATIONS_QUANT_UNIQUELY_MAPPED_HOST {
    tag "comb_annots_quant_pathogen"
    label 'process_high'

    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
    storeDir "${params.outdir}/HTSeq"

    
    input: 
    path quantification_table
    path annotation_table
    val attribute
    val organism

    output:
    path "pathogen_combined_quant_annotations.tsv", emit: tsv
    path "versions.yml"           , emit: versions

    script:
    """
    $workflow.projectDir/bin/combine_quant_annotations.py \\
    -q $quantification_table \\
    -annotations $annotation_table \\
    -a $attribute -org $organism
    """
}