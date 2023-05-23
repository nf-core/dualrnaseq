	process TXIMPORT_HOST {
        
        tag "$meta.id"
	    
   	    label 'process_high'

        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'nfcore/dualrnaseq:dev' :
        'nfcore/dualrnaseq:dev' }"

	    input: 
	        tuple val(meta), file("host_quant.sf") 
	        file (annotations) 

	    output:
	        tuple val(meta), file ("${meta.id}_host_quant_gene_level.sf"), emit: salmon_files_to_combine_gene_level

	    script:
	    """
	    tximport.R $annotations $meta.id
	    """
	}