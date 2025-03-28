	process REPLACE_GENE_FEATURE_GFF {
	    tag "repl_gene_feature_gff"
	    label 'process_high'

	    input:
            path(gff)
            val(features)
	    output:
    	    path "${outfile_name}"

	    script:
            outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature.gff3")
            """
            $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features
            """
	}