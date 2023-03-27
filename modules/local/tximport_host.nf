	process TXIMPORT_HOST {
        
        tag "tximport_host"
	    publishDir "${params.outdir}/salmon/${prefix}", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon/${prefix}"
	    
   	    label 'process_high'

        conda "bioconda::r-tximport"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
            'https://depot.galaxyproject.org/singularity/bioconductor-tximport:1.0.3--r3.2.2_0' :
            'davelabhub/tximport:v20201109' }"

	    input: 
	        tuple val(prefix), file("salmon/${prefix}/*") 
	        file (annotations) 

	    output:
	        file "${prefix}_host_quant_gene_level.sf", emit: salmon_files_to_combine_gene_level

	    script:
	    
	    template 'tximport.R'
	   
	}