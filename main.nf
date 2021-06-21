#!/usr/bin/env nextflow
/*

==========================================

nf-core/dualrnaseq

==========================================

------------------------------------------
Authors: B.Mika-Gospodorz and R.Hayward

Homepage: https://github.com/nf-core/dualrnaseq
------------------------------------------


Brief description:
------------------------------------------
Extract host and pathogen expression using three methods, with options for various statistical outputs.
1) STAR + HTSeq
2) Salmon selective alignment
3) STAR + Salmon - alignment-based mode

*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/dualrnaseq --input '*_R{1,2}.fastq.gz' -profile docker"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}


////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

//----------
// Check if genomes exists in the config file
//----------

if (params.genomes && params.genome_host && !params.genomes.containsKey(params.genome_host)) {
    exit 1, "The provided genome '${params.genome_host}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

if (params.genomes && params.genome_pathogen && !params.genomes.containsKey(params.genome_pathogen)) {
    exit 1, "The provided genome '${params.genome_pathogen}' is not available in the genomes.config file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}


//----------
// Reference genomes, annotation files and transcriptomes
//----------
params.fasta_host = params.genome_host ? params.genomes[ params.genome_host ].fasta_host ?: false : false
if (params.fasta_host) { ch_fasta_host = file(params.fasta_host, checkIfExists: true) }

params.fasta_pathogen = params.genome_pathogen ? params.genomes[ params.genome_pathogen ].fasta_pathogen ?: false : false
if (params.fasta_pathogen) { ch_fasta_pathogen = file(params.fasta_pathogen, checkIfExists: true) }

params.gff_host_tRNA = params.genome_host ? params.genomes[ params.genome_host ].gff_host_tRNA ?: false : false
if (params.gff_host_tRNA) { ch_gff_host_tRNA = file(params.gff_host_tRNA, checkIfExists: true) }

params.gff_host_genome = params.genome_host ? params.genomes[ params.genome_host ].gff_host ?: false : false
if (params.gff_host_genome) { ch_gff_host_genome = file(params.gff_host_genome, checkIfExists: true) }

params.gff_pathogen = params.genome_pathogen ? params.genomes[ params.genome_pathogen ].gff_pathogen ?: false : false
if (params.gff_pathogen) { ch_gff_pathogen = file(params.gff_pathogen, checkIfExists: true) }

if(params.read_transcriptome_fasta_host_from_file){
    params.transcriptome_host = params.genome_host ? params.genomes[ params.genome_host ].transcriptome_host ?: false : false
    if (params.transcriptome_host) { ch_transcriptome_host = file(params.transcriptome_host, checkIfExists: true) }
}

if(params.read_transcriptome_fasta_pathogen_from_file){
    params.transcriptome_pathogen = params.genome_pathogen ? params.genomes[ params.genome_pathogen ].transcriptome_pathogen ?: false : false
    if (params.transcriptome_pathogen) { ch_transcriptome_pathogen = file(params.transcriptome_pathogen, checkIfExists: true) }
}


//----------
// Trimming - check if only one of the trimming tools is set to "true"
//----------
if (params.run_cutadapt & params.run_bbduk) {
	exit 1, "Trimming: both --run_cutadapt and --run_bbduk are set to true, please use only one of the adapter trimming tools"
}

//----------
// BBDuk - fasta file of adapters
//----------
if(params.run_bbduk) {
	if (params.adapters) { adapter_database = file(params.adapters, checkIfExists: true) }
}


//----------
// Salmon library type
//----------
if (params.run_salmon_selective_alignment | params.run_salmon_alignment_based_mode){
	if (!params.libtype){
    exit 1, "Salmon: Please specify --libtype"
} else if (params.libtype == 'A'){
  // continue
} else if (params.single_end & (params.libtype != 'U' && params.libtype != 'SR' && params.libtype != 'SF')) {
	    exit 1, "Salmon: Invalid library type --libtype ${params.libtype}! Library types available for single-end reads are:'U', 'SR', 'SF'."
} else if (!params.single_end & (params.libtype != 'IU' && params.libtype != 'ISR' && params.libtype != 'ISF' && params.libtype != 'MU' && params.libtype != 'MSR' && params.libtype != 'MSF' && params.libtype != 'OU' && params.libtype != 'OSR' && params.libtype != 'OSF' )) {
	    exit 1, "Salmon: Invalid library type --libtype ${params.libtype}! Library types available for paired-end reads are: 'IU', 'ISR', 'ISF', 'MU', 'MSR', 'MSF','OU', 'OSR', 'OSF'."
}
}


//----------
// HTSeq - check if there is an alignment file
//----------
if (params.run_htseq_uniquely_mapped){
	if (!params.run_star){
    exit 1, "HTSeq: There is no alignment file. Please specify --run_star"
  }
}


//----------
// Salmon alignment based mode - check if there is an alignment file
//----------
if (params.run_salmon_alignment_based_mode){
	if (!params.run_star){
    exit 1, "Salmon: There is no alignment file. Please specify --run_star"
  }
}

//----------
// Mapping stats
//----------
if(params.mapping_statistics) {
	if (params.rna_classes_to_replace_host) { ch_RNA_classes = file(params.rna_classes_to_replace_host, checkIfExists: true) }
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}


//----------
// Stage config files
ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)

/*
 * Create a channel for input read files
 */
if (params.input_paths) {
    if (params.single_end) {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
            .into { ch_read_files_fastqc; trimming_reads; raw_read_count; scatter_plots_set  }
    } else {
        Channel
            .from(params.input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
            .into { ch_read_files_fastqc; trimming_reads; raw_read_count;scatter_plots_set }
    }
} else {
    Channel
        .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
        .into { ch_read_files_fastqc; trimming_reads; raw_read_count; scatter_plots_set}
}



//----------
// Channel for host and pathogen fasta files
//----------
Channel
    .value( ch_fasta_pathogen)
    .collect()
    .set { genome_fasta_pathogen_to_unzip}

Channel
    .value( ch_fasta_host )
    .set { genome_fasta_host_to_unzip}


//----------
// Channel for host GFF and/or tRNA-based files
//----------
if(params.gff_host_tRNA){
	Channel
	    .value(ch_gff_host_tRNA)
	    .set {host_gff_trna_file_to_unzip}
	Channel
	    .value(ch_gff_host_genome)
	    .set {host_gff_trna_to_unzip}
}else{
	Channel
	    .value(ch_gff_host_genome)
	    .set {host_gff_to_unzip}
}


//----------
// Channel for pathogen GFF files
//----------
Channel
    .value(ch_gff_pathogen)
    .set {pathogen_gff_to_unzip}


//----------
// Channel for host and pathogen transcriptomes
//----------
if(params.read_transcriptome_fasta_host_from_file){
Channel
    .value(ch_transcriptome_host)
    .into {host_transcriptome_to_combine; transcriptome_host_to_split_q_table_salmon; transcriptome_host_to_split_table_salmon; transcriptome_host_to_split_q_table_salmon_alignment_based; transcriptome_host_to_split_table_salmon_alignment; transcriptome_fasta_host_ref_names}
}

if(params.read_transcriptome_fasta_pathogen_from_file){
Channel
    .value(ch_transcriptome_pathogen)
    .into {pathogen_transcriptome_to_combine; transcriptome_pathogen_to_split_table_salmon; transcriptome_pathogen_to_split_table_salmon_alignment; transcriptome_pathogen_to_split_q_table_salmon; transcriptome_pathogen_to_split_q_table_salmon_alignment_based;transcriptome_fasta_pathogen_ref_names}
}


//----------
// Channel to capture Cutadapt-based params
//----------
if (params.run_cutadapt){
	if(params.single_end){
		Channel
		    .value( params.a )
		    .set { adapter_sequence_3 }
	} else {
		Channel
		    .from (params.a, params.A)
		    .collect()
		    .set { adapter_sequence_3 }
	}
	Channel
	 .value(params.quality_cutoff)
	 .set { quality_cutoff}
}


//----------
// Channel to capture BBDuk-based params
//----------
if (params.run_bbduk){
	Channel
	    .value( adapter_database )
	    .set { adapter_database }
}


//----------
// Channel to capture Salmon-based params
//----------
if (params.run_salmon_selective_alignment){
	Channel
	    .value(params.kmer_length)
	    .set {kmer_length_salmon_index}

}


if (params.run_salmon_selective_alignment | params.run_salmon_alignment_based_mode) {
	Channel
	    .value(params.gene_attribute_gff_to_create_transcriptome_host)
	    .into {host_gff_attribute_salmon_alignment_tRNA; gene_attribute_gff_to_create_transcriptome_host_salmon; host_atr_collect_data_salmon; combine_annot_quant_pathogen; combine_annot_quant_host; atr_scatter_plot_pathogen; atr_scatter_plot_host; attribute_quant_stats_salmon; host_annotations_RNA_class_stats_pathogen; attribute_host_RNA_class_stats; host_atr_collect_data_salmon_alignment_mode; combine_annot_quant_pathogen_salmon_alignment_based; combine_annot_quant_host_salmon_alignment_based; atr_scatter_plot_pathogen_alignment; atr_scatter_plot_host_alignment; attribute_quant_stats_salmon_alignment;host_annotations_RNA_class_stats_pathogen_alignment; attribute_host_RNA_class_stats_alignment}

	// Convert csv into a Groovy list, then make a channel
	gene_feature_gff_to_create_transcriptome_host_list = params.gene_feature_gff_to_create_transcriptome_host ? params.gene_feature_gff_to_create_transcriptome_host.split(',').collect{it.trim()} : []
	Channel
	    .value(gene_feature_gff_to_create_transcriptome_host_list)
	    .collect()
	    .into { gene_feature_gff_host_salmon_alignment; gene_feature_gff_to_create_transcriptome_host_salmon}

	Channel
	    .value(params.gene_attribute_gff_to_create_transcriptome_pathogen)
	    .into {pathogen_gff_attribute_salmon_alignment; gene_attribute_gff_to_create_transcriptome_pathogen_salmon}

	// Convert csv into a Groovy list, then make a channel
	gene_feature_gff_to_create_transcriptome_pathogen_list = params.gene_feature_gff_to_create_transcriptome_pathogen ? params.gene_feature_gff_to_create_transcriptome_pathogen.split(',').collect{it.trim()} : []
	Channel
	    .value(gene_feature_gff_to_create_transcriptome_pathogen_list)
	    .collect()
	    .into {gene_feature_to_quantify_pathogen_salmon_alignment; gene_feature_to_extract_annotations_pathogen; gene_feature_gff_to_create_transcriptome_pathogen_salmon}

	Channel
	    .value(params.libtype)
	    .into {libtype_salmon; libtype_salmon_alignment_mode}
}


//----------
// Channel to capture htseq and star-based params
//----------
if(params.run_htseq_uniquely_mapped | params.run_star){

	// Convert csv into a Groovy list, then make a channel
	gene_feature_gff_to_quantify_host_list = params.gene_feature_gff_to_quantify_host ? params.gene_feature_gff_to_quantify_host.split(',').collect{it.trim()} : []
	Channel
	    .value(gene_feature_gff_to_quantify_host_list)
	    .collect()
	    .into {gene_feature_to_quantify_host; gene_feature_to_extract_annotations_host_htseq}

	// Convert csv into a Groovy list, then make a channel
	gene_feature_gff_to_quantify_pathogen_list = params.gene_feature_gff_to_quantify_pathogen ? params.gene_feature_gff_to_quantify_pathogen.split(',').collect{it.trim()} : []
	Channel
	    .value(gene_feature_gff_to_quantify_pathogen_list)
	    .collect()
	    .into {gene_feature_to_quantify_pathogen; gene_feature_to_extract_annotations_pathongen_htseq}

	Channel
	    .value(params.pathogen_gff_attribute)
	    .into { pathogen_gff_attribute; pathogen_gff_attribute_to_extract_annotations_htseq}

	Channel
	    .value(params.stranded)
	    .set { stranded_htseq_unique}

	Channel
	    .value(params.host_gff_attribute)
	    .into { host_gff_attribute_to_pathogen; host_gff_attribute_htseq; host_gff_attribute_htseq_combine; host_gff_attribute_to_extract_annotations_htseq; host_gff_attribute_mapping_stats_htseq; host_gff_attribute_RNA_class_pathogen_htseq; host_gff_attribute_RNA_class_host_htseq; combine_annot_quant_pathogen_host_gff_attribute; host_gff_attribute_htseq_TPM; atr_scatter_plot_pathogen_htseq_u_m; atr_scatter_plot_host_htseq_u_m}
}


//----------
// Channel to capture mapping statistics-based  params
//----------
if(params.mapping_statistics) {

Channel
    .value(ch_RNA_classes)
    .into { RNA_classes_to_replace; RNA_classes_to_replace_alignment; RNA_classes_to_replace_htseq_uniquely_mapped}
}




/*
--------------------------------------------------------------------

HEADER LOG INFO

--------------------------------------------------------------------
*/

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = workflow.runName
summary['Input']            = params.input
summary['Host Fasta Ref']          = params.fasta_host
summary['Pathogen Fasta Ref']      = params.fasta_pathogen
summary['Host tRNA gff Ref']       = params.gff_host_tRNA
summary['Host genome gff Ref']     = params.gff_host_genome
summary['Pathogen genome gff Ref'] = params.gff_pathogen
summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName

if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}

summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-dualrnaseq-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/dualrnaseq Workflow Summary'
    section_href: 'https://github.com/nf-core/dualrnaseq'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }




/*
--------------------------------------------------------------------

PARSE SOFTWARE VERSION NUMBERS

--------------------------------------------------------------------
*/


process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    python --version > v_python.txt
    R --version > v_r.txt
    cutadapt --version > v_cutadapt.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    STAR --version > v_star.txt
    htseq-count . . --version > v_htseq.txt
    samtools --version > v_samtools.txt
    gffread --version > v_gffread.txt
    salmon --version > v_salmon.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
--------------------------------------------------------------------
Workflow - Processes
--------------------------------------------------------------------
*/


/*
 * STEP 1 - Check if there are technical replicates - to plot scatter plots of technical replicates
 */

if(params.mapping_statistics) {

	scatter_plots_set
			.map { tag, file -> tag }
			.set {scatter_plots}

	process check_replicates {
	    tag "check_replicates"

	    label 'process_high'

	    input:
	    val(sample_name) from scatter_plots.collect()

	    output:
	    stdout repl_scatter_plots_salmon_pathogen
	    stdout repl_scatter_plots_salmon_host
	    stdout repl_scatter_plots_salmon_alignment_host
	    stdout repl_scatter_plots_salmon_alignment_pathogen
	    stdout repl_scatter_plots_htseq_pathogen
	    stdout repl_scatter_plots_htseq_host

	    shell:
	    '''
	    python !{workflow.projectDir}/bin/check_replicates.py -s !{sample_name} 2>&1
	    '''
	}
}



/*
 * STEP 2 - Prepare reference files
 */

//Uncompress Pathogen genome
process uncompress_pathogen_fasta_genome {
    tag "uncompress_pathogen_genome"
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode

    label 'process_high'

    input:
    each file(f_ext) from genome_fasta_pathogen_to_unzip.collect()

    output:
    file "${base_name_file}.fasta" into genome_fasta_pathogen_to_combine
    file "${base_name_file}.fasta" into genome_fasta_pathogen_ref_names
    file "${base_name_file}.fasta" into genome_fasta_pathogen_to_transcriptome

    shell:
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()
    if (ext_file == "fasta" | ext_file == "fa"){
	    '''
	    cp -n !{f_ext} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "zip"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
      gunzip -f -S .zip !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else if(ext_file == "gz"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
	    gunzip -f !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else {
      '''
      echo "Your pathogen genome files appear to have the wrong extension. \n Currently, the pipeline only supports .fasta or .fa, or compressed files with .zip or .gz extensions."
      '''
    }
}


//Uncompress Host genome
process uncompress_host_fasta_genome {
    tag "uncompress_host_genome"
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode

    label 'process_high'

    input:
    file(f_ext) from genome_fasta_host_to_unzip

    output:
    file "${base_name_file}.fasta" into genome_fasta_host_to_combine
    file "${base_name_file}.fasta" into genome_fasta_host_to_decoys
    file "${base_name_file}.fasta" into genome_fasta_host_ref_names
    file "${base_name_file}.fasta" into genome_fasta_host_to_transcriptome
    file "${base_name_file}.fasta" into genome_fasta_host_to_transcriptome_tRNA

    shell:
    //Tests to see if the input files are compressed.
    //At this stage, only accepts fasta or .fa files, or .gz or .zip genome files.
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()
    if (ext_file == "fasta" | ext_file == "fa"){
    	'''
    	cp -n !{f_ext} !{base_name_file}.fasta
    	'''
    }else if(ext_file == "zip"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
  	  '''
  	  gunzip -f -S .zip !{f_ext}
  	  cp -n !{old_base_name_file} !{base_name_file}.fasta
  	  '''
    }else if(ext_file == "gz"){
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.fasta|.fa/,"")
	    '''
	    gunzip -f !{f_ext}
	    cp -n !{old_base_name_file} !{base_name_file}.fasta
	    '''
    }else {
      '''
      echo "Your host genome files appear to have the wrong extension. \n Currently, the pipeline only supports .fasta or .fa, or compressed files with .zip or .gz extensions."
      '''
    }
}



//Uncompress Pathogen GFF
process uncompress_pathogen_gff {
    tag "uncompress_pathogen_GFF"
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode

    label 'process_high'

    input:
    file(f_ext) from pathogen_gff_to_unzip

    output:
    file "${base_name_file}.gff3" into gff_feature_quant_pathogen_salmon_alignment
    file "${base_name_file}.gff3" into gff_pathogen_create_transcriptome
    file "${base_name_file}.gff3" into gff_feature_quant_pathogen_htseq
    file "${base_name_file}.gff3" into extract_annotations_pathogen_gff_htseq

    shell:
    //Tests to see if the input files are compressed.
    //At this stage, only accepts .gff or gff3, or .gz or .zip annotation files.
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()

    if (ext_file == "gff" | ext_file == "gff3"){
      '''
      cp -n !{f_ext} !{base_name_file}.gff3
      '''
    }else if(ext_file == "zip"){
      '''
      gunzip -f -S .zip !{f_ext}
      cp -n !{base_name_file} !{base_name_file}.gff3
      '''
    }else if(ext_file == "gz"){
      //if gff or gff3, need to save as .gff3
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.gff|.gff3/,"")
      '''
      gunzip -f !{f_ext}
      cp -n !{old_base_name_file} !{base_name_file}.gff3
      '''
    }else {
      '''
      echo "Your pathogen GFF file appears to be in the wrong format or has the wrong extension. \n Currently, the pipeline only supports .gff or .gff3, or compressed files with .zip or .gz extensions."
      '''
    }
}




//Uncompress host GFF (no tRNA)
process uncompress_host_gff {
    tag "uncompress_host_GFF"
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode

    label 'process_high'

    input:
    file(f_ext) from host_gff_to_unzip

    output:
    file "${base_name_file}.gff3" into gff_host_genome_star_salmon_change_atr
    file "${base_name_file}.gff3" into gff_host_create_transcriptome
    file "${base_name_file}.gff3" into gff_host_genome_htseq
    file "${base_name_file}.gff3" into extract_annotations_host_gff_htseq
    file "${base_name_file}.gff3" into gff_host_star_alignment_gff
    file "${base_name_file}.gff3" into gff_host_star_htseq_alignment_gff
    file "${base_name_file}.gff3" into genome_gff_star_index

    shell:
    //Tests to see if the input files are compressed.
    //At this stage, only accepts .gff or gff3, or .gz or .zip annotation files.
    ext_file = f_ext.getExtension()
    base_name_file = f_ext.getBaseName()

    if (ext_file == "gff" | ext_file == "gff3"){
      '''
      cp -n !{f_ext} !{base_name_file}.gff3
      '''
    }else if(ext_file == "zip"){
      '''
      gunzip -f -S .zip !{f_ext}
      cp -n !{base_name_file} !{base_name_file}.gff3
      '''
    }else if(ext_file == "gz"){
      //if gff or gff3, need to save as .gff3
      old_base_name_file = base_name_file
      base_name_file = old_base_name_file.replaceAll(/.gff|.gff3/,"")
      '''
      gunzip -f !{f_ext}
      cp -n !{old_base_name_file} !{base_name_file}.gff3
      '''
    }else {
      '''
      echo "Your host GFF file appears to be in the wrong format or has the wrong extension. \n Currently, the pipeline only supports .gff or .gff3, or compressed files with .zip or .gz extensions."
      '''
    }
}




if(params.gff_host_tRNA){
  //Uncompress host GFF (with tRNA)
  process uncompress_host_gff_trna {
      tag "uncompress_host_GFF_trna"
      publishDir "${params.outdir}/references", mode: params.publish_dir_mode

      label 'process_high'

      input:
      file(f_ext) from host_gff_trna_to_unzip

      output:
      file "${base_name_file}.gff3" into gff_host_genome_star_salmon_change_atr
      file "${base_name_file}.gff3" into gff_host_create_transcriptome
      file "${base_name_file}.gff3" into combine_gff_host_genome_htseq
      file "${base_name_file}.gff3" into gff_host_star_alignment_gff
      file "${base_name_file}.gff3" into gff_host_star_htseq_alignment_gff
      file "${base_name_file}.gff3" into genome_gff_star_index

      shell:
      //Tests to see if the input files are compressed.
      //At this stage, only accepts .gff or gff3, or .gz or .zip annotation files.
      ext_file = f_ext.getExtension()
      base_name_file = f_ext.getBaseName()

      if (ext_file == "gff" | ext_file == "gff3"){
        '''
        cp -n !{f_ext} !{base_name_file}.gff3
        '''
      }else if(ext_file == "zip"){
        '''
        gunzip -f -S .zip !{f_ext}
        cp -n !{base_name_file} !{base_name_file}.gff3
        '''
      }else if(ext_file == "gz"){
        //if gff or gff3, need to save as .gff3
        old_base_name_file = base_name_file
        base_name_file = old_base_name_file.replaceAll(/.gff|.gff3/,"")
        '''
        gunzip -f !{f_ext}
        cp -n !{old_base_name_file} !{base_name_file}.gff3
        '''
      }else {
        '''
        echo "Your host GFF file appears to be in the wrong format or has the wrong extension. \n Currently, the pipeline only supports .gff or .gff3, or compressed files with .zip or .gz extensions."
        '''
      }
  }


    //Uncompress host GFF (with tRNA)
  process uncompress_host_gff_trna_file {
      tag "uncompress_host_GFF_trna_file"
      publishDir "${params.outdir}/references", mode: params.publish_dir_mode

      label 'process_high'

      input:
      file(f_ext) from host_gff_trna_file_to_unzip

      output:
      file "${base_name_file}.gff3" into change_attrubute_gff_host_tRNA_salmon_alignment
      file "${base_name_file}.gff3" into gff_host_create_transcriptome_tRNA
      file "${base_name_file}.gff3" into combine_gff_host_tRNA_htseq

      shell:
      //Tests to see if the input files are compressed.
      //At this stage, only accepts .gff or gff3, or .gz or .zip annotation files.
      ext_file = f_ext.getExtension()
      base_name_file = f_ext.getBaseName()

      if (ext_file == "gff" | ext_file == "gff3"){
        '''
        cp -n !{f_ext} !{base_name_file}.gff3
        '''
      }else if(ext_file == "zip"){
        '''
        gunzip -f -S .zip !{f_ext}
        cp -n !{base_name_file} !{base_name_file}.gff3
        '''
      }else if(ext_file == "gz"){
        //if gff or gff3, need to save as .gff3
        old_base_name_file = base_name_file
        base_name_file = old_base_name_file.replaceAll(/.gff|.gff3/,"")
        '''
        gunzip -f !{f_ext}
        cp -n !{old_base_name_file} !{base_name_file}.gff3
        '''
      }else {
        '''
        echo "Your host GFF tRNA file appears to be in the wrong format or has the wrong extension. \n Currently, the pipeline only supports .gff or .gff3, or compressed files with .zip or .gz extensions."
        '''
      }
  }

} // end of "if(params.gff_host_tRNA){"




/*
 * combine pathogen and host genome fasta files
 */

process combine_pathogen_host_fasta_genome {
    tag "combine_genome_fa_files"
    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
    storeDir "${params.outdir}/references"

    label 'process_high'

    input:
    file(host_fa) from genome_fasta_host_to_combine
    file(pathogen_fa) from genome_fasta_pathogen_to_combine.collect()

    output:
    file "host_pathogen.fasta" into host_pathogen_fasta_index
    file "host_pathogen.fasta" into host_pathogen_fasta_star_index
    file "host_pathogen.fasta" into genome_fasta_file_host_pathogen_to_decoy_transcriptome

    script:
    """
    cat $pathogen_fa $host_fa > host_pathogen.fasta
    """
}


if(params.run_htseq_uniquely_mapped) {

	/*
	 * combine host genome and tRNA gff files
	 */

	if(params.gff_host_tRNA){
		process combine_host_genome_tRNA_gff_files_htseq {
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"
		    tag "comb_host_genome_tRNA_gff"

		    label 'process_high'

		    input:
		    file(host_gff_genome) from combine_gff_host_genome_htseq
		    file(host_gff_tRNA) from combine_gff_host_tRNA_htseq

		    output:
		    file "${outfile_name}" into gff_host_genome_htseq
		    file "${outfile_name}" into extract_annotations_host_gff_htseq

		    script:
		    outfile_name = host_gff_genome[0].toString().replaceAll(/.gff3|.gff/,"_with_tRNA.gff3")
		    """
		    cat $host_gff_genome $host_gff_tRNA > ${outfile_name}
		    """
		}
	}


	/*
	 * Replace gene features (3rd column of gff) with 'quant' in host gff
	 */

	process replace_gene_feature_gff_host_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "repl_gene_feature_gff_host"

	    label 'process_high'

	    input:
	    file(gff) from gff_host_genome_htseq
	    val(features) from gene_feature_to_quantify_host

	    output:
	    file "${outfile_name}" into combine_gff_host
	    file "${outfile_name}" into gff_host_to_TPM

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature.gff3")
	    """
	    $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features
	    """
	}


	/*
	 * Replace gene features (3rd column of gff) with 'quant' in pathogen gff
	 */

	process replace_gene_feature_gff_pathogen_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "repl_gene_feature_gff_pathogen"

	    label 'process_high'

	    input:
	    file(gff) from gff_feature_quant_pathogen_htseq
	    val(features) from gene_feature_to_quantify_pathogen

	    output:
	    file "${outfile_name}" into to_replace_gff_attribute

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature.gff3")
	    """
	    $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features
	    """
	}


	/*
	 * Replace pathogen gene attribute in 9th column of pathogen gff with host attribute
	 */

	process replace_attribute_pathogen_gff_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "repl_attribute_pathogen_gff"

	    label 'process_high'

	    input:
	    file(gff) from to_replace_gff_attribute
	    val(host_attribute) from host_gff_attribute_to_pathogen
	    val(pathogen_attribute) from pathogen_gff_attribute

	    output:
	    file "${outfile_name}" into combine_gff_pathogen
	    file "${outfile_name}" into gff_pathogen_to_TPM

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_new_attribute.gff3")
	    """
	    $workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} $host_attribute $pathogen_attribute
	    """
	}


	/*
	 * combine pathogen and host genome gff files
	 */

	process combine_pathogen_host_gff_files_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "comb_pathogen_host_gff_files"

	    label 'process_high'

	    input:
	    file(host_gff) from combine_gff_host
	    file(pathogen_gff_genome) from combine_gff_pathogen

	    output:
	    file "host_pathogen_htseq.gff" into quantification_gff_u_m

	    script:
	    """
	    cat $pathogen_gff_genome $host_gff > host_pathogen_htseq.gff
	    """
	}


	/*
	* extract annotations from pathogen gff files - for RNA class statistics and to merge with quantification results
	*/

	process extract_annotations_pathogen_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "extract_annotations_pathogen"

	    label 'process_high'

	    input:
	    file gff from extract_annotations_pathogen_gff_htseq
	    val(features) from gene_feature_to_extract_annotations_pathongen_htseq
	    val(pathogen_attribute) from pathogen_gff_attribute_to_extract_annotations_htseq

	    output:
	    file "${outfile_name}*_htseq.tsv" into pathogen_annotations_RNA_class_stats_htseq
	    file "${outfile_name}*_htseq.tsv" into annotation_pathogen_combine_quant_htseq_u_m
	    file "${outfile_name}*_htseq.tsv" into annotation_pathogen_split_quant_htseq

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"")
	    """
	    python $workflow.projectDir/bin/extract_annotations_from_gff.py -gff $gff -f $features -a $pathogen_attribute -org pathogen -q_tool htseq -o ${outfile_name}
	    """
	}


	/*
	* extract annotations from host gff files - for RNA class statistics and to merge with quantification results
	*/

	process extract_annotations_host_htseq {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "extract_annotations_host"

	    label 'process_high'

	    input:
	    file gff from extract_annotations_host_gff_htseq
	    val(features) from gene_feature_to_extract_annotations_host_htseq
	    val(host_attribute) from host_gff_attribute_to_extract_annotations_htseq

	    output:
	    file "${outfile_name}*_htseq.tsv" into host_annotations_RNA_class_stats_htseq
	    file "${outfile_name}*_htseq.tsv" into annotation_host_combine_quant_htseq
	    file "${outfile_name}*_htseq.tsv" into annotation_host_split_quant_htseq

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"")
	    """
	    python $workflow.projectDir/bin/extract_annotations_from_gff.py -gff $gff -f $features -a $host_attribute -org host -q_tool htseq -o ${outfile_name}
	    """
	}
}



if(params.run_star | params.run_salmon_alignment_based_mode) {

	if(params.mapping_statistics) {

		/*
		* extract reference names from host genome fasta files - for mapping stats
		*/

		process extract_reference_names_host_star {
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"
		    tag "extract_ref_names_host_star"

		    label 'process_high'

		    input:
		    file(host_fa) from genome_fasta_host_ref_names

		    output:
		    file "reference_host_names.txt" into reference_host_names_uniquelymapped
		    file "reference_host_names.txt" into reference_host_names_crossmapped_find
		    file "reference_host_names.txt" into reference_host_names_multimapped

		    script:
		    """
		    $workflow.projectDir/bin/extract_reference_names_from_fasta_files.sh reference_host_names.txt $host_fa
		    """
		}



		/*
		* extract reference names from pathogen genome fasta files - for mapping stats
		*/

		process extract_reference_names_pathogen_star {
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"
		    tag "extract_ref_names_pathgn_star"

		    label 'process_high'

		    input:
		    file(pathogen_fa) from genome_fasta_pathogen_ref_names.collect()

		    output:
		    file "reference_pathogen_names.txt" into reference_pathogen_names_uniquelymapped
		    file "reference_pathogen_names.txt" into reference_pathogen_crossmapped_find
		    file "reference_pathogen_names.txt" into reference_pathogen_names_multimapped

		    script:
		    """
		    $workflow.projectDir/bin/extract_reference_names_from_fasta_files.sh reference_pathogen_names.txt $pathogen_fa
		    """
		}
	}
}



if(params.run_salmon_selective_alignment | params.run_salmon_alignment_based_mode) {


	/*
	 * Replace 'Parent' gene attribute in 9th column of host gff with 'parent'
	 */

	process replace_attribute_host_genome_gff_star_salmon {
	    tag "repl_attribute_host_genome_gff"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(gff) from gff_host_genome_star_salmon_change_atr

	    output:
	    file "${outfile_name}" into combine_gff_host_genome_star_salmon

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_parent_attribute.gff3")

	    """
	    $workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} parent Parent
	    """
	}




	if(params.gff_host_tRNA){

		/*
		 * Replace gene attribute in 9th column of host tRNA gff with 'parent'
		 */

		process replace_attribute_host_tRNA_gff_star_salmon {
		    tag "repl_attribute_host_tRNA_gff"
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"

		    label 'process_high'

		    input:
		    file(gff) from change_attrubute_gff_host_tRNA_salmon_alignment
		    val(host_attribute) from host_gff_attribute_salmon_alignment_tRNA

		    output:
		    file "${outfile_name}" into combine_host_gffs

		    script:
		    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_parent_attribute.gff3")
		    """
		    $workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} parent $host_attribute
		    """
		}

		/*
		 * Combine tRNA gff file with genome host gff
	 	 */

		process combine_host_genome_tRNA_gff_star_salmon {
		    tag "comb_host_genome_tRNA_gff"
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"

		    label 'process_high'

		    input:
		    file(host_gff_genome) from combine_gff_host_genome_star_salmon
		    file(host_gff_tRNA) from combine_host_gffs

		    output:
		    file "${outfile_name}" into gff_host_tRNA_genome_salmon_alignment

		    script:
		    outfile_name = host_gff_genome[0].toString().replaceAll(/.gff3|.gff/,"_with_tRNA_STAR_salmon.gff3")
		    """
		    cat $host_gff_genome $host_gff_tRNA > ${outfile_name}
		    """
		}
	}

	gff_host_genome_salmon_alignment = params.gff_host_tRNA ? gff_host_tRNA_genome_salmon_alignment : combine_gff_host_genome_star_salmon


	/*
	 * Replace gene features (3rd column of gff) with 'quant' in host gff
	 */

	process replace_gene_feature_gff_host_salmon {
	    tag "repl_gene_feature_gff_host"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(gff) from gff_host_genome_salmon_alignment
	    val(features) from gene_feature_gff_host_salmon_alignment

	    output:
	    file "${outfile_name}" into combine_gff_host_salmon_alignment
	    file "${outfile_name}" into extract_annotations_host_gff_salmon

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature_salmon_alignment.gff3")
	    """
	    $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features
	    """
	}


	/*
	 * Replace gene attribute in 9th column of pathogen gff with 'parent'
	 */

	process replace_attribute_pathogen_gff_star_salmon {
	    tag "repl_attribute_pathogen_gff"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(gff) from gff_feature_quant_pathogen_salmon_alignment
	    val(pathogen_attribute) from pathogen_gff_attribute_salmon_alignment

	    output:
	    file "${outfile_name}" into to_replace_gff_feature_salmon_alignment
	    file "${outfile_name}" into extract_annotations_pathogen_gff_salmon


	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_parent_attribute.gff3")
	    """
	    $workflow.projectDir/bin/replace_attribute_gff.sh $gff ${outfile_name} parent $pathogen_attribute
	    """
	}


	/*
	 * Extract annotations from pathogen gff file - for RNA class statistics and to merge with quantification results
	 */

	process extract_annotations_pathogen_salmon {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "extract_gff_annots_pathogen"

	    label 'process_high'

	    input:
	    file gff from extract_annotations_pathogen_gff_salmon
	    val(features) from gene_feature_to_extract_annotations_pathogen

	    output:
	    file "${outfile_name}*_salmon.tsv" into pathogen_annotations_RNA_class_stats
	    file "${outfile_name}*_salmon.tsv" into pathogen_annotations_RNA_class_stats_salmon_alignment
	    file "${outfile_name}*_salmon.tsv" into annotation_pathogen_combine_quant
	    file "${outfile_name}*_salmon.tsv" into annotation_pathogen_combine_quant_salmon_alignment_based

	    script:
	    outfile_name = "pathogen_gff_annotations"
	    """
	    python $workflow.projectDir/bin/extract_annotations_from_gff.py -gff $gff -f $features -a parent -org pathogen -q_tool salmon -o ${outfile_name}
	    """
	}


	/*
	 * Extract annotations from host gff file - for RNA class statistics and to merge with quantification results
	 */

	process extract_annotations_host_salmon {
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"
	    tag "extract_gff_annots_host"

	    label 'process_high'

	    input:
	    file gff from extract_annotations_host_gff_salmon

	    output:
	    file "${outfile_name}*_salmon.tsv" into host_annotations_RNA_class_stats
	    file "${outfile_name}*_salmon.tsv" into host_annotations_RNA_class_stats_salmon_alignment
	    file "${outfile_name}*_salmon.tsv" into tximport_annotations
	    file "${outfile_name}*_salmon.tsv" into host_annotations_uniq_ambig
            file "${outfile_name}*_salmon.tsv" into tximport_annotations_salmon_alignment
	    file "${outfile_name}*_salmon.tsv" into host_annotations_uniq_ambig_AB
	    file "${outfile_name}*_salmon.tsv" into annotation_host_combine_quant
	    file "${outfile_name}*_salmon.tsv" into annotation_host_combine_quant_salmon_alignment_based
	    file "${outfile_name}*_salmon.tsv" into annotation_host_combine_quant_gene_level_salmon
	    file "${outfile_name}*_salmon.tsv" into annotation_host_combine_quant_gene_level_salmon_alignment

	    script:
	    outfile_name = "host_gff_annotations"
	    """
	    python $workflow.projectDir/bin/extract_annotations_from_gff.py -gff $gff -f quant -a parent -org host -q_tool salmon -o ${outfile_name}
	    """
	}



	if(!params.read_transcriptome_fasta_host_from_file){

		/*
		 * create host transcriptome fasta file
		 */

		process create_transcriptome_fasta_host {
		    tag "create_transcripts_host"
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"

		    label 'process_high'

		    input:
		    file(gff) from gff_host_create_transcriptome
		    file(host_fa) from genome_fasta_host_to_transcriptome

		    output:
		    file "${outfile_name}" into host_transcriptome
		    file "${outfile_name}" into transcriptome_host_to_split_table_salmon_without_tRNA
		    file "${outfile_name}" into transcriptome_host_to_split_table_salmon_alignment_without_tRNA
		    file "${outfile_name}" into transcriptome_host_to_split_q_table_salmon_without_tRNA
		    file "${outfile_name}" into transcriptome_host_to_split_q_table_salmon_alignment_based_without_tRNA
		    file "${outfile_name}" into host_transcriptome_to_combine_host
		    file "${outfile_name}" into transcriptome_fasta_host_ref_names_without_tRNA

		    script:
		    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_transcriptome.fasta")
		    """
		    gffread -w $outfile_name -g $host_fa $gff
		    """
	}


		if(params.gff_host_tRNA){

			/*
			 * Create host tRNA transcriptome fasta file
			 */

			process create_transcriptome_fasta_host_tRNA {
			    tag "create_transcripts_tRNA_host"
			    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
			    storeDir "${params.outdir}/references"

			    label 'process_high'

			    input:
			    file(gff) from gff_host_create_transcriptome_tRNA
			    file(host_fa) from genome_fasta_host_to_transcriptome_tRNA
			    val(features) from gene_feature_gff_to_create_transcriptome_host_salmon
			    val(attribute) from gene_attribute_gff_to_create_transcriptome_host_salmon

			    output:
			    file "${outfile_name}" into host_transcriptome_to_combine_tRNA

			    script:
			    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_transcriptome.fasta")
			    """
			    python $workflow.projectDir/bin/gff_to_fasta_transcriptome.py -fasta $host_fa -gff $gff  -f $features -a $attribute -o $outfile_name
			    """
			}

			/*
			 * Combine host genome transcriptome with tRNA transcriptome
			 */

			process combine_host_fasta_transcriptomes {
			    tag "comb_host_fa_tRNA_transcripts"
			    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
			    storeDir "${params.outdir}/references"

			    label 'process_high'

			    input:
			    file(host_tr_fa) from host_transcriptome_to_combine_host
			    file(host_tRNA_tr_fa) from host_transcriptome_to_combine_tRNA

			    output:
			    file "host_transcriptome.fasta" into host_transcriptome_genome_tRNA
                            file "host_transcriptome.fasta" into transcriptome_host_to_split_table_salmon_with_tRNA
                            file "host_transcriptome.fasta" into transcriptome_host_to_split_table_salmon_alignment_with_tRNA
                            file "host_transcriptome.fasta" into transcriptome_host_to_split_q_table_salmon_with_tRNA
                            file "host_transcriptome.fasta" into transcriptome_host_to_split_q_table_salmon_alignment_based_with_tRNA
                            file "host_transcriptome.fasta" into transcriptome_fasta_host_ref_names_with_tRNA

 			    script:
			    """
			    cat $host_tr_fa $host_tRNA_tr_fa > host_transcriptome.fasta
			    """
			}
		}

	host_transcriptome_to_combine = params.gff_host_tRNA ? host_transcriptome_genome_tRNA : host_transcriptome
	transcriptome_host_to_split_q_table_salmon = params.gff_host_tRNA ? transcriptome_host_to_split_q_table_salmon_with_tRNA : transcriptome_host_to_split_q_table_salmon_without_tRNA
	transcriptome_host_to_split_table_salmon = params.gff_host_tRNA ? transcriptome_host_to_split_table_salmon_with_tRNA : transcriptome_host_to_split_table_salmon_without_tRNA
	transcriptome_host_to_split_q_table_salmon_alignment_based = params.gff_host_tRNA ? transcriptome_host_to_split_q_table_salmon_alignment_based_with_tRNA : transcriptome_host_to_split_q_table_salmon_alignment_based_without_tRNA
	transcriptome_host_to_split_table_salmon_alignment = params.gff_host_tRNA ? transcriptome_host_to_split_table_salmon_alignment_with_tRNA : transcriptome_host_to_split_table_salmon_alignment_without_tRNA
	transcriptome_fasta_host_ref_names = params.gff_host_tRNA ? transcriptome_fasta_host_ref_names_with_tRNA : transcriptome_fasta_host_ref_names_without_tRNA
}



	if(!params.read_transcriptome_fasta_pathogen_from_file){

		/*
		 * create pathogen transcriptome fasta file
		 */

		process create_transcriptome_fasta_pathogen {
		    tag "create_transcripts_fa_pathogen"
		    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/references"

		    label 'main'
		    label 'process_high'

		    input:
		    file(gff) from gff_pathogen_create_transcriptome
		    file(pathogen_fa) from genome_fasta_pathogen_to_transcriptome.collect()
		    val(features) from gene_feature_gff_to_create_transcriptome_pathogen_salmon
		    val(attribute) from gene_attribute_gff_to_create_transcriptome_pathogen_salmon

		    output:
		    file "${outfile_name}" into pathogen_transcriptome_to_combine
		    file "${outfile_name}" into transcriptome_pathogen_to_split_table_salmon
		    file "${outfile_name}" into transcriptome_pathogen_to_split_table_salmon_alignment
		    file "${outfile_name}" into transcriptome_pathogen_to_split_q_table_salmon
		    file "${outfile_name}" into transcriptome_pathogen_to_split_q_table_salmon_alignment_based
		    file "${outfile_name}" into transcriptome_fasta_pathogen_ref_names

		    script:
		    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_transcriptome.fasta")
		    """
		    python $workflow.projectDir/bin/gff_to_fasta_transcriptome.py -fasta $pathogen_fa -gff $gff -f $features -a $attribute  -o $outfile_name
		    """
		}
	}


	/*
	 * combine pathogen and host transcriptome fasta files
	 */


	process combine_pathogen_host_fasta_transcriptome {
	    tag "comb_pathogen_host_transcripts"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(host_tr_fa) from host_transcriptome_to_combine
	    file(pathogen_tr_fa) from pathogen_transcriptome_to_combine

	    output:
	    file "host_pathogen_transcriptome.fasta" into transcriptome_fasta_file_host_pathogen_to_decoy_transcriptome
	    file "host_pathogen_transcriptome.fasta" into transcriptome_salmon_alignment_based_mode

	    script:
	    """
	    cat $pathogen_tr_fa $host_tr_fa > host_pathogen_transcriptome.fasta
	    """
	}

}


if(params.run_salmon_alignment_based_mode){

	/*
	 * Replace gene features (3rd column of gff) with 'quant' in pathogen gff
	 */

	process replace_gene_feature_gff_pathogen_salmon {
	    tag "repl_gene_feature_gff_pathogen"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(gff) from to_replace_gff_feature_salmon_alignment
	    val(features) from gene_feature_to_quantify_pathogen_salmon_alignment

	    output:
	    file "${outfile_name}" into combine_gff_pathogen_salmon_alignment

	    script:
	    outfile_name = gff[0].toString().replaceAll(/.gff3|.gff/,"_quant_feature_salmon_alignment.gff3")
	    """
	    $workflow.projectDir/bin/replace_feature_gff.sh $gff ${outfile_name} $features
	    """
	}


	/*
	 * Combine pathogen gff with host gff
	 */


	process combine_pathogen_host_gff_files_salmon {
	    tag "combine_pathogen_host_gff"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(host_gff) from combine_gff_host_salmon_alignment
	    file(pathogen_gff_genome) from combine_gff_pathogen_salmon_alignment

	    output:
	    file "host_pathogen_star_alignment_mode.gff" into gff_host_pathogen_star_salmon_alignment_gff

	    script:
	    """
	    cat $pathogen_gff_genome $host_gff > host_pathogen_star_alignment_mode.gff
	    """
	}
}


/*
 * STEP 3 - FastQC
 */

if (!params.skip_fastqc) {
	process fastqc {
	    tag "$name"
	    label 'process_medium'

	    publishDir "${params.outdir}/fastqc", mode: params.publish_dir_mode
      //saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
	    storeDir "${params.outdir}/fastqc"

	    input:
	    set val(name), file(reads) from ch_read_files_fastqc

	    output:
	    file "*_fastqc.{zip,html}" into ch_fastqc_results

	    script:
	    fastqc_params = params.fastqc_params
	    """
	    fastqc --quiet --threads $task.cpus --noextract $reads $fastqc_params
	    """
	}
}else{
   Channel.empty()
      .set {ch_fastqc_results}
}


/*
 *  STEP 4 - Trimming
 */


if (params.run_cutadapt | params.run_bbduk) {

	/*
	 *  run Cutadapt
	 */

	if (params.run_cutadapt) {
		process trimming {
		    tag "$name_reads"
		    publishDir "${params.outdir}/trimming_cutadapt", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/trimming_cutadapt"

		    label 'process_high'

		    input:
		    set val(name), file(reads) from trimming_reads
		    val adapter_seq_3 from adapter_sequence_3
		    val q_value from quality_cutoff

		    output:
		    set val(name_sample), file("${name_sample}{_1,_2,}_trimmed.fastq.gz") into trimming_results_star_htseq, trimming_results_to_salmon, trimming_results_to_qc, trimming_results_star_salmon

		    script:
		    cutadapt_params = params.cutadapt_params
			if (params.single_end){
		    		name_reads = reads.toString().replaceAll(/:/,"_")
		    		name_reads = name_reads.replaceAll(/-/,"_")
		    		name_out = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
		    		name_sample = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"")

		    		"""
		    		cutadapt -j ${task.cpus} -q $q_value -a $adapter_seq_3 -m 1 -o ${name_out} $reads $cutadapt_params
		    		"""
			} else{
				name_reads = name.replaceAll(/:/,"_")
				name_reads = name_reads.replaceAll(/-/,"_")
				name_sample = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"")
				name_1 = reads[0][0].toString().replaceAll(/:/,"_")
				name_1 = name_1.replaceAll(/-/,"_")
				name_1 = name_1.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
				name_2 = reads[1][0].toString().replaceAll(/:/,"_")
				name_2 = name_2.replaceAll(/-/,"_")
				name_2 = name_2.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
				"""
				cutadapt -j ${task.cpus} -q $q_value -a ${adapter_seq_3[0]} -A ${adapter_seq_3[1]} -o ${name_1} -p ${name_2} -m 1 ${reads[0]} ${reads[1]} $cutadapt_params
				"""
			}
		}

	}

	/*
	 *  run BBDuk
	 */

	if (params.run_bbduk) {
		process bbduk {
		    tag "$name_reads"
		    publishDir "${params.outdir}/trimming_bbduk", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/trimming_bbduk"

		    label 'process_high'

		    input:
		    set val(name), file(reads) from trimming_reads
		    file adapters from adapter_database

		    output:
		    set val(name_sample), file("${name_sample}{_1,_2,}_trimmed.fastq.gz") into trimming_results_star_htseq, trimming_results_to_salmon, trimming_results_to_qc, trimming_results_star_salmon
		    file "*.log"

		    script:
		    minlen = params.minlen
		    qtrim = params.qtrim
		    trimq = params.trimq
		    ktrim = params.ktrim
		    k = params.k
		    mink = params.mink
		    hdist = params.hdist
		    bbduk_params = params.bbduk_params
		    if (params.single_end){
		    		name_reads = reads.toString().replaceAll(/:/,"_")
		    		name_reads = name_reads.replaceAll(/-/,"_")
		    		name_out = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
		    		name_sample = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"")
				fileoutput = name_sample + '.log'
		    		"""
		    		bbduk.sh -Xmx1g in=$reads out=${name_out} ref=$adapters minlen=$minlen qtrim=$qtrim trimq=$trimq ktrim=$ktrim k=$k mink=$mink hdist=$hdist &> $fileoutput $bbduk_params
		    		"""
		    } else{
				name_reads = name.replaceAll(/:/,"_")
				name_reads = name_reads.replaceAll(/-/,"_")
				name_sample = name_reads.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"")
				name_1 = reads[0][0].toString().replaceAll(/:/,"_")
				name_1 = name_1.replaceAll(/-/,"_")
				name_1 = name_1.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
				name_2 = reads[1][0].toString().replaceAll(/:/,"_")
				name_2 = name_2.replaceAll(/-/,"_")
				name_2 = name_2.replaceAll(/.fastq.gz|.fq.gz|.fastq|.fq/,"_trimmed.fastq.gz")
				fileoutput = name_sample + '.log'
				"""
				bbduk.sh -Xmx1g in1=${reads[0]} in2=${reads[1]} out1=${name_1} out2=${name_2} ref=$adapters minlen=$minlen qtrim=$qtrim trimq=$trimq ktrim=$ktrim k=$k mink=$mink hdist=$hdist $bbduk_params tpe tbo &> $fileoutput
				"""
			}
		}

	}
}else{
   trimming_reads
      .into {trimming_results_to_qc; trimming_results_star_htseq; trimming_results_to_salmon; trimming_results_star_salmon}
}



/*
 * STEP 5 -FastQC after trimming
 */



if (params.run_cutadapt | params.run_bbduk & !params.skip_fastqc) {
	process fastqc_after_trimming {
	    tag "$sample_name"
	    publishDir "${params.outdir}/fastqc_after_trimming", mode: params.publish_dir_mode
      //saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
	    storeDir "${params.outdir}/fastqc_after_trimming"

	    label 'process_medium'

	    input:
	    set val(name),file(reads) from trimming_results_to_qc

	    output:
	    file "*_trimmed_fastqc.{zip,html}" into ch_fastqc_trimmed_results

	    script:
	    sample_name = name.replaceFirst(/.fastq.gz|.fq.gz|.fastq|.fq/, "")
	    fastqc_params = params.fastqc_params
	    """
	    fastqc --threads ${task.cpus} --quiet --noextract $reads $fastqc_params
	    """
	}
}else{
   Channel.empty()
      .set {ch_fastqc_trimmed_results}
}



if(params.run_cutadapt | params.run_bbduk & params.mapping_statistics) {

	raw_read_count
	    .map { tag, file -> file }
	    .set {raw_read_count_file}
	 /*
	 * count total number of raw single-end reads
	 */
	 process count_total_reads {
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

	if (!params.single_end){

		/*
		 * count total number of paired-end reads
		 */
		process count_total_read_pairs {
		    tag "count_total_reads"
		    publishDir "${params.outdir}/mapping_statistics", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics"

		    label 'process_high'

		    input:
		    file(tsv) from to_collect_total_reads.collect()
		    output:
		    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_salmon
		    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_salmon_alignment
		    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_star
		    file "total_raw_read_pairs_fastq.tsv" into collect_total_reads_raw_star_for_salmon

		    script:
		    """
		    $workflow.projectDir/bin/collect_total_raw_read_pairs.py -i $tsv
		    """
		}
	}else{
	   to_collect_total_reads
		  .into {collect_total_reads_raw_salmon; collect_total_reads_raw_salmon_alignment; collect_total_reads_raw_star; collect_total_reads_raw_star_for_salmon}
	}
}else{
   Channel.empty()
  .into {collect_total_reads_raw_salmon; collect_total_reads_raw_salmon_alignment; collect_total_reads_raw_star; collect_total_reads_raw_star_for_salmon}
}



/*
 * STEP 6 - Salmon Selective Alignment
 */


if(params.run_salmon_selective_alignment) {

	/*
	 * create genome as decoy sequence
	 */

	process create_decoy_transcriptome_file {
	    tag "create_decoy_transcripts_file"
	    publishDir "${params.outdir}/references", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/references"

	    label 'process_high'

	    input:
	    file(host_fa) from genome_fasta_host_to_decoys
      	    //file(host_tr_fa) from host_transcriptome_to_combine
            file(host_pathogen_genome_fasta) from genome_fasta_file_host_pathogen_to_decoy_transcriptome
	    file(host_pathogen_transcriptome_fasta) from transcriptome_fasta_file_host_pathogen_to_decoy_transcriptome

	    output:
	    file "gentrome.fasta" into salmon_index_gentrome
	    file "decoys.txt" into salmon_index_decoys

	    shell:
	    '''
	    grep ">" !{host_fa} | cut -d " " -f 1 > decoys.txt
	    sed -i -e 's/>//g' decoys.txt
	    cat !{host_pathogen_transcriptome_fasta} !{host_fa} > gentrome.fasta
	    '''
	}


	/*
	 * Salmon - build an index
	 */

	process salmon_index {
	    tag "salmon_index"
            storeDir "${params.outdir}/salmon"

            label 'process_high'

	    input:
	    file(gentrome) from salmon_index_gentrome
	    file(decoys) from salmon_index_decoys
	    val(kmer_length) from kmer_length_salmon_index

	    output:
	    file "transcripts_index" into salmon_index

	    script:
	    salmon_sa_params_index = params.salmon_sa_params_index
	    keepDuplicates = params.keepDuplicates ? "--keepDuplicates" : ''
	    """
	    salmon index -t $gentrome -i transcripts_index --decoys $decoys -k $kmer_length -p ${task.cpus} $keepDuplicates $salmon_sa_params_index
	    """
	}


	/*
	 * Salmon - quantification
	 */


	process salmon_quantification {
	    storeDir "${params.outdir}/salmon"
	    tag "${sample}"

            label 'process_high'

	    input:
	    file(index) from salmon_index.collect()
	    set val(sample), file(reads) from trimming_results_to_salmon
	    val(libtype) from libtype_salmon

	    output:
	    set val(sample_name), file("${sample_name}") into split_table
	    set val(sample_name), file("${sample_name}") into split_table_uniq_ambig
	    file("${sample_name}") into salmon_files_to_combine
	    file("${sample_name}") into multiqc_salmon_quant
	    set val(sample_name), file("${sample_name}") into collect_processed_read_counts

	    script:
	    UnmappedNames = params.writeUnmappedNames ? "--writeUnmappedNames" : ''
	    softclip = params.softclipOverhangs ? "--softclipOverhangs" : ''
	    incompatPrior = params.incompatPrior
	    dumpEq = params.dumpEq ? "--dumpEq" : ''
	    salmon_sa_params_mapping = params.salmon_sa_params_mapping
	    if (params.single_end){
	    	sample_name = sample.replaceFirst(/.fastq.gz|.fq.gz|.fastq|.fq/, "")
		writeMappings = params.writeMappings ? "--writeMappings=$sample_name/mapping.sam" : ''
		"""
		salmon quant -p ${task.cpus} -i $index -l $libtype -r $reads $softclip --incompatPrior $incompatPrior $UnmappedNames --validateMappings $dumpEq $writeMappings -o $sample_name $salmon_sa_params_mapping
		"""
	    } else{
		sample_name = sample.replaceFirst(/.fastq.gz|.fq.gz|.fastq|.fq/, "")
		writeMappings = params.writeMappings ? "--writeMappings=$sample_name/mapping.sam" : ''
		"""
		salmon quant -p ${task.cpus} -i $index -l $libtype -1 ${reads[0]} -2 ${reads[1]} $softclip --incompatPrior $incompatPrior $UnmappedNames --validateMappings $dumpEq $writeMappings -o $sample_name $salmon_sa_params_mapping
 		"""
	    }
	}



	/*
	 * Salmon - split quantification tables into host and pathogen results
	 */


	process split_table_salmon_each {
            publishDir "${params.outdir}/salmon/${sample_name}", mode: params.publish_dir_mode
            storeDir "${params.outdir}/salmon/${sample_name}"
            tag "split_quantification ${sample_name}"

            label 'process_high'

            input:
            set val(sample_name), file ("salmon/*") from split_table
            file transcriptome_pathogen from transcriptome_pathogen_to_split_q_table_salmon
            file transcriptome_host from transcriptome_host_to_split_q_table_salmon

            output:
            set val(sample_name), file("host_quant.sf") into salmon_host_tximport

            script:
            """
            $workflow.projectDir/bin/split_quant_tables_salmon.sh $transcriptome_pathogen $transcriptome_host  salmon/*/quant.sf "quant.sf"
            """
        }


	if(params.generate_salmon_uniq_ambig) {

	    /*
	     * Extract and combine the ambig and unique counts
	     */

		process extract_ambig_uniq_transcripts_genes {
		    publishDir "${params.outdir}/salmon/${sample_name}/aux_info", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon/${sample_name}/aux_info"
		    tag "extract_ambig_uniq_transcripts_genes ${sample_name}"

		    label 'process_high'

		    input:
		    set val(sample_name), file("salmon/*") from split_table_uniq_ambig
		    file (annotations) from host_annotations_uniq_ambig


		    output:
		    file "${sample_name}_host_quant_ambig_uniq.sf"
		    file "${sample_name}_pathogen_quant_ambig_uniq.sf"
		    file "${sample_name}_host_quant_ambig_uniq_gene_level.sf"
		    set val(sample_name), file("${sample_name}_host_quant_ambig_uniq.sf") into host_files_comb
		    set val(sample_name), file("${sample_name}_pathogen_quant_ambig_uniq.sf") into path_files_comb

		    script:
		    """
		    $workflow.projectDir/bin/salmon_extract_ambig_uniq_transcripts_genes.R salmon/*/quant.sf salmon/*/aux_info/ambig_info.tsv $sample_name $annotations
		    """
		}

	    /*
	     * Combine the host ambig and unique counts
	     */
		process host_comb_ambig_uniq {
		    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon"
		    tag "host_comb_ambig_uniq"

		    label 'process_high'

		    input:
		    file("salmon/*") from host_files_comb.collect()

		    output:
		    file "host_quant_combined_ambig_uniq.tsv"

		    script:
		    """
		    $workflow.projectDir/bin/salmon_host_comb_ambig_uniq.R salmon/*/aux_info/*_host_quant_ambig_uniq.sf
		    """
		}

	    /*
	     * Combine the pathogen ambig and unique counts
	     */
		process pathogen_comb_ambig_uniq {
		    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon"
		    tag "pathogen_comb_ambig_uniq"

		    label 'process_high'

		    input:
		    file("salmon/*") from path_files_comb.collect()

		    output:
		    file "pathogen_quant_combined_ambig_uniq.tsv"

		    script:
		    """
		    $workflow.projectDir/bin/salmon_pathogen_comb_ambig_uniq.R salmon/*/aux_info/*_pathogen_quant_ambig_uniq.sf
		    """
		}
	}



	/*
	 * Salmon - combine quantification results
	 */


	process combine_quantification_tables_salmon {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon"
	    tag "combine_quantification_salmon"
	    label 'process_high'

	    input:
	    file input_quantification from salmon_files_to_combine.collect()
	    val gene_attribute from host_atr_collect_data_salmon

	    output:
	    file "combined_quant.tsv" into split_table_salmon

	    script:
	    """
	    python $workflow.projectDir/bin/collect_quantification_data.py -i $input_quantification -q salmon -a $gene_attribute -org both
	    """
	}



	/*
	 * Salmon - split quantification tables into host and pathogen results
	 */


	process split_quantification_tables_salmon {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
            tag "split_quant_table_salmon"

            label 'process_high'

            input:
            file quant_table from split_table_salmon
            file transcriptome_pathogen from transcriptome_pathogen_to_split_table_salmon
            file transcriptome_host from transcriptome_host_to_split_table_salmon

            output:
            file 'host_quant_salmon.tsv' into host_quantification_mapping_stats_salmon
            file 'pathogen_quant_salmon.tsv' into pathogen_quantification_mapping_stats_salmon
            file 'host_quant_salmon.tsv' into host_quantification_RNA_stats_salmon
            file 'pathogen_quant_salmon.tsv' into pathogen_quantification_RNA_stats_salmon
            file 'host_quant_salmon.tsv' into quant_host_add_annotations
            file 'pathogen_quant_salmon.tsv' into quant_pathogen_add_annotations
            file 'host_quant_salmon.tsv' into quant_scatter_plot_host
            file 'pathogen_quant_salmon.tsv' into quant_scatter_plot_pathogen
	    env pathonen_tab into pathonen_tab into scatterplots_pathogen_salmon
	    env host_tab into host_tab into scatterplots_host_salmon

            script:
            """
            $workflow.projectDir/bin/split_quant_tables_salmon.sh $transcriptome_pathogen $transcriptome_host $quant_table "quant_salmon.tsv"
            pathonen_tab=\$(if [ \$(cat pathogen_quant_salmon.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
            host_tab=\$(if [ \$(cat host_quant_salmon.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
            """
        }


	/*
	 * Salmon - combine pathogen annotations extracted from gff with quantification tables
	 */

	process combine_annotations_quant_pathogen {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon"
	    tag "comb_annots_quant_pathgn_salmon"

   	    label 'process_high'

	    input:
	    file quantification_table from quant_pathogen_add_annotations
	    file annotation_table from annotation_pathogen_combine_quant
	    val attribute from combine_annot_quant_pathogen

	    output:
	    file "pathogen_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org pathogen
	    """
	}


	/*
	 * Salmon - combine host annotations extracted from gff with quantification tables
	 */

	process combine_annotations_quant_host_salmon {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon"
	    tag "comb_annots_quant_host_salmon"

   	    label 'process_high'

	    input:
	    file quantification_table from quant_host_add_annotations
	    file annotation_table from annotation_host_combine_quant
	    val attribute from combine_annot_quant_host

	    output:
	    file "host_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org host
	    """
	}



	/*
	 * Tximport - host
	 */

	process tximport_host {
	    publishDir "${params.outdir}/salmon/${sample_name}", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon/${sample_name}"
	    tag "tximport_host"

   	    label 'process_high'

	    input:
	    set val(sample_name), file("salmon/${sample_name}/*") from salmon_host_tximport
	    file (annotations) from tximport_annotations

	    output:
	    file "${sample_name}_host_quant_gene_level.sf" into salmon_files_to_combine_gene_level

	    script:
	    """
	    $workflow.projectDir/bin/tximport.R salmon $annotations $sample_name
	    """
	}


	/*
	 * Salmon - Combine host gene level quantification results estimated with Tximport
	 */

	process combine_host_quant_gene_level_salmon {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon"
	    tag "comb_host_quant_genes_salmon"

	    label 'process_high'

	    input:
	    file input_quantification from salmon_files_to_combine_gene_level.collect()

	    output:
	    file "host_combined_gene_level.tsv" into quant_gene_level_host_add_annotations_salmon

	    script:
	    """
	    python $workflow.projectDir/bin/collect_quantification_data.py -i $input_quantification -q salmon -a gene_id -org host_gene_level
	    """
	}


	/*
	 * Salmon - combine host annotations extracted from gff with gene-level quantification estimates
	 */

	process combine_annotations_quant_gene_level_salmon {
	    publishDir "${params.outdir}/salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon"
	    tag "comb_annots_gene_host_salmon"

	    label 'process_high'

	    input:
	    file quantification_table from quant_gene_level_host_add_annotations_salmon
	    file annotation_table from annotation_host_combine_quant_gene_level_salmon

	    output:
	    file "host_combined_quant_gene_level_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_annotations_salmon_gene_level.py -q $quantification_table -annotations $annotation_table -a gene_id -org host
	    """
	}



	if(params.mapping_statistics) {

		/*
		 * Salmon - plot scatter plots of technical replicates for pathogen
		 */

		process scatter_plot_pathogen_salmon {
		    publishDir "${params.outdir}/mapping_statistics/salmon/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/scatter_plots"
		    tag "scatter_plot_salmon_pathogen"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_pathogen
		    val attribute from atr_scatter_plot_pathogen
		    val replicates from repl_scatter_plots_salmon_pathogen
		    val pathogen_table_non_empty from scatterplots_pathogen_salmon

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    pathogen_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org pathogen
		    """
		}


		/*
		 * Salmon - plot scatter plots of technical replicates for host
		 */

		process scatter_plot_host_salmon {
		    publishDir "${params.outdir}/mapping_statistics/salmon/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/scatter_plots"
		    tag "scatter_plot_salmon_host"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_host
		    val attribute from atr_scatter_plot_host
		    val replicates from repl_scatter_plots_salmon_host
		    val host_table_non_empty from scatterplots_host_salmon

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    host_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org host
		    """
		}


		/*
		 * Salmon - Extract processed reads from log file
		 */

		process extract_processed_reads {
		    publishDir "${params.outdir}/mapping_statistics/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon"
		    tag "extract_processed_reads"

		    label 'process_high'

		    input:
		    set val(sample_name), file ("salmon/*") from collect_processed_read_counts

		    output:
		    file "${sample_name}.txt" into collect_results

		    script:
		    """
		    $workflow.projectDir/bin/extract_processed_reads.sh salmon/*/aux_info/meta_info.json $sample_name salmon
		    """
		}


		/*
		 * Salmon - Collect processed reads from all samples
		 */

		process collect_processed_reads {
		    publishDir "${params.outdir}/mapping_statistics/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon"
		    tag "collect_processed_reads"

		    label 'process_high'

		    input:
		    file process_reads from collect_results.collect()

		    output:
		    file "processed_reads_salmon.tsv" into mapping_stats_total_reads

		    script:
		    """
		    cat $process_reads > processed_reads_salmon.tsv
		    """
		}


		/*
		 * Salmon - Collect mapping statistics
		 */

		process salmon_quantification_stats {
		    publishDir "${params.outdir}/mapping_statistics/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon"
		    tag "quantification_stats_salmon"

		    label 'process_high'

		    input:
		    file quant_table_host from host_quantification_mapping_stats_salmon
		    file quant_table_pathogen from pathogen_quantification_mapping_stats_salmon
		    val attribute from attribute_quant_stats_salmon
		    file total_processed_reads from mapping_stats_total_reads
		    file total_raw_reads from collect_total_reads_raw_salmon.ifEmpty('.')

		    output:
		    file ('salmon_host_pathogen_total_reads.tsv') into salmon_mapped_stats_to_plot

		    script:
		    """
		    python $workflow.projectDir/bin/mapping_stats.py -q_p $quant_table_pathogen -q_h $quant_table_host -total_processed $total_processed_reads -total_raw $total_raw_reads -a $attribute -t salmon -o salmon_host_pathogen_total_reads.tsv
		    """
		}

		/*
		 * Salmon - Plot mapping statistics
		 */

		process plot_salmon_mapping_stats_host_pathogen {
		    tag "$name2"
		    publishDir "${params.outdir}/mapping_statistics/salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon"

		    label 'process_high'

		    input:
		    file(stats) from salmon_mapped_stats_to_plot

		    output:
		    file "*.tsv"
		    file "*.pdf"

		    script:
		    """
		    python $workflow.projectDir/bin/plot_mapping_statistics_salmon.py -i $stats
		    """
		}



		/*
		 * Salmon - Collect pathogen RNA class statistics
		 */


		process RNA_class_statistics_salmon_pathogen {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_pathogen", mode: params.publish_dir_mode
		    tag "rna_class_stats_pathgn_salmon"

		    label 'process_high'

		    input:
		    file quant_table from pathogen_quantification_RNA_stats_salmon
		    val attribute from host_annotations_RNA_class_stats_pathogen
		    file gene_annotations from pathogen_annotations_RNA_class_stats

		    output:
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen_combined
		    file "pathogen_RNA_classes_sum_counts_*.tsv"
		    stdout plot_RNA_stats_pathogen_boolean
		    stdout plot_RNA_stats_pathogen_combined_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -q_tool salmon -org pathogen 2>&1
		    '''
		}


		/*
		 * Salmon - Collect host RNA class statistics
		 */

		process RNA_class_statistics_salmon_host {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_host", mode: params.publish_dir_mode
		    tag "rna_class_stats_host_salmon"

		    label 'process_high'

		    input:
		    file quant_table from host_quantification_RNA_stats_salmon
		    val attribute from attribute_host_RNA_class_stats
		    file gene_annotations from host_annotations_RNA_class_stats
		    file rna_classes_to_replace from RNA_classes_to_replace

		    output:
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_combined
		    file "host_RNA_classes_sum_counts_*.tsv"
		    file "host_gene_types_groups_*"
		    stdout plot_RNA_stats_host_combined_boolean
		    stdout plot_RNA_stats_host_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -rna !{rna_classes_to_replace} -q_tool salmon -org host 2>&1
		    '''
		}


		/*
		 * Salmon - Plot pathogen RNA class statistics for each sample separately
		 */

		process plot_RNA_class_salmon_pathogen_each {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_pathogen"
		    tag "plot_RNA_stats_pathogen_salmon"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen
		    val plot_rna from plot_RNA_stats_pathogen_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}


		/*
		 * Salmon - Plot pathogen RNA class statistics for all samples
		 */

		process plot_RNA_class_salmon_pathogen_combined {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_pathogen"
		    tag "plot_RNA_stats_comb_pathogen"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen_combined
		    val plot_rna from plot_RNA_stats_pathogen_combined_boolean

		    output:
		    file "RNA_class_stats_combined_pathogen.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org pathogen
		    """
		}

		/*
		 * Salmon - Plot host RNA class statistics for each sample separately
		 */

		process plot_RNA_class_salmon_host_each {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_host"
		    tag "plot_RNA_stats_host_salmon"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host
		    val plot_rna from plot_RNA_stats_host_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}

		/*
		 * Salmon - Plot host RNA class statistics for all samples
		 */

		process plot_RNA_class_salmon_host_combined {
		    publishDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon/RNA_classes_host"
		    tag "plot_RNA_stats_comb_host"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host_combined
		    val plot_rna from plot_RNA_stats_host_combined_boolean

		    output:
		    file "RNA_class_stats_combined_host.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org host
		    """
		}
	}
}else{
   Channel.empty()
      .set {multiqc_salmon_quant}
}


/*
 * STEP 7 - Salmon alignment_based_mode
 */


if (params.run_salmon_alignment_based_mode){

	/*
	 * STAR - index
	 */

	process STARindex_salmon_alignment {
		publishDir "${params.outdir}/STAR_for_salmon", mode: params.publish_dir_mode
		storeDir "${params.outdir}/STAR_for_salmon"
		tag "STAR_index"

         	label 'process_high'

		input:
		file(fasta) from host_pathogen_fasta_index
		file(gff) from genome_gff_star_index

		output:
		file "index/*" into star_index_transcriptome_alignment

		script:
		sjdbOverhang = params.sjdbOverhang
		star_salmon_index_params = params.star_salmon_index_params
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir index/ --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $sjdbOverhang $star_salmon_index_params
		"""
	}


	/*
	 * STAR - alignment
	 */

	process ALIGNMENTstar_for_salmon {
	    tag "$sample_name"
	    publishDir "${params.outdir}/STAR_for_salmon", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/STAR_for_salmon"

	    label 'process_medium'

	    input:
	    set val(sample_name),file(reads) from  trimming_results_star_salmon
	    file(gff) from gff_host_pathogen_star_salmon_alignment_gff.collect()
	    file(index) from star_index_transcriptome_alignment.collect()

	    output:
	    file "${sample_name}/*" into multiqc_star_for_salmon_alignment
	    set val(sample_name), file("${sample_name}/${sample_name}Aligned.toTranscriptome.out.bam") into salmon_quantify_alignment_based_mode
	    file "${sample_name}/*"
	    set val(sample_name), file("${sample_name}/${sample_name}Log.final.out") into collect_processed_read_counts_STAR_for_salmon

	    script:
	    outSAMunmapped = params.outSAMunmapped
	    outSAMattributes = params.outSAMattributes
	    quantTranscriptomeBan = params.quantTranscriptomeBan
	    outFilterMultimapNmax = params.outFilterMultimapNmax
	    outFilterType = params.outFilterType
	    alignSJoverhangMin = params.alignSJoverhangMin
	    alignSJDBoverhangMin = params.alignSJDBoverhangMin
	    outFilterMismatchNmax = params.outFilterMismatchNmax
	    outFilterMismatchNoverReadLmax = params.outFilterMismatchNoverReadLmax
	    alignIntronMin = params.alignIntronMin
	    alignIntronMax = params.alignIntronMax
	    alignMatesGapMax = params.alignMatesGapMax
	    limitBAMsortRAM = params.limitBAMsortRAM
	    readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
	    winAnchorMultimapNmax = params.winAnchorMultimapNmax
	    star_salmon_alignment_params = params.star_salmon_alignment_params

	    if (params.single_end){
	    	"""
	    	mkdir $sample_name
	    	STAR --runThreadN ${task.cpus} --genomeDir . --sjdbGTFfile $gff $readFilesCommand --readFilesIn $reads --outSAMtype BAM Unsorted --outSAMunmapped $outSAMunmapped --outSAMattributes $outSAMattributes --outFileNamePrefix $sample_name/$sample_name --sjdbGTFfeatureExon quant --sjdbGTFtagExonParentTranscript parent --quantMode TranscriptomeSAM --quantTranscriptomeBan $quantTranscriptomeBan --outFilterMultimapNmax $outFilterMultimapNmax --outFilterType $outFilterType --limitBAMsortRAM $limitBAMsortRAM --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --winAnchorMultimapNmax $winAnchorMultimapNmax $star_salmon_alignment_params
	    	"""
	    } else {
	    	"""
	    	mkdir $sample_name
	    	STAR --runThreadN ${task.cpus} --genomeDir . --sjdbGTFfile $gff $readFilesCommand --readFilesIn ${reads[0]} ${reads[1]} --outSAMtype BAM Unsorted --outSAMunmapped $outSAMunmapped --outSAMattributes $outSAMattributes --outFileNamePrefix $sample_name/$sample_name --sjdbGTFfeatureExon quant --sjdbGTFtagExonParentTranscript parent --quantMode TranscriptomeSAM --quantTranscriptomeBan $quantTranscriptomeBan --outFilterMultimapNmax $outFilterMultimapNmax --outFilterType $outFilterType --limitBAMsortRAM $limitBAMsortRAM --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --winAnchorMultimapNmax $winAnchorMultimapNmax $star_salmon_alignment_params
	    	"""
	    }
	}


	/*
	 * Salmon alignment-based mode - quantification
	 */

	process salmon_quantification_alignment_based_mode {
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "${sample}"

	    label 'process_high'

	    input:
	    file(transcriptome) from transcriptome_salmon_alignment_based_mode.collect()
	    set val(sample), file(bam_file) from salmon_quantify_alignment_based_mode
	    val(libtype) from libtype_salmon_alignment_mode

	    output:
	    set val(sample_name), file("${sample_name}") into split_table_alignment_based
	    set val(sample_name), file("${sample_name}") into split_table_uniq_ambig_ab
	    file("${sample_name}") into salmon_files_to_combine_alignment_mode
	    file("${sample_name}") into multiqc_salmon_alignment_quant
	    set val(sample_name), file("${sample_name}") into collect_processed_read_counts_alignment_based

	    script:
	    incompatPrior = params.incompatPrior
	    salmon_alignment_based_params = params.salmon_alignment_based_params
	    sample_name = sample.replaceFirst(/.fastq.gz|.fq.gz|.fastq|.fq/, "")
	    """
	    salmon quant -p ${task.cpus} -t $transcriptome -l $libtype -a $bam_file --incompatPrior $incompatPrior -o $sample_name $salmon_alignment_based_params
	    """
	}


	/*
	 * Salmon alignment-based mode - split quantification tables into host and pathogen results
	 */

	process split_table_salmon_each_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode/${sample_name}", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode/${sample_name}"
	    tag "split_quant_tbl_for_sal_modes"

	    label 'process_high'

	    input:
	    set val(sample_name), file ("salmon/*") from split_table_alignment_based
	    //set val(sample_name), file ("salmon/*") from split_table_alignment_based_uniq_ambig
	    file transcriptome_pathogen from transcriptome_pathogen_to_split_q_table_salmon_alignment_based
	    file transcriptome_host from transcriptome_host_to_split_q_table_salmon_alignment_based

	    output:
	    set val(sample_name), file("host_quant.sf") into salmon_alignment_host_tximport
	    set val(sample_name), file("pathogen_quant.sf")

            script:
            """
            $workflow.projectDir/bin/split_quant_tables_salmon.sh $transcriptome_pathogen $transcriptome_host salmon/*/quant.sf "quant.sf"
            """
	}




     if(params.generate_salmon_uniq_ambig) {

	    /*
	     * Extract and combine the ambig and unique counts
	     */
	   	 process extract_ambig_uniq_transcripts_genes_alignment_based {
		    publishDir "${params.outdir}/salmon_alignment_mode/${sample_name}/aux_info", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon_alignment_mode/${sample_name}/aux_info"
		    tag "extract_ambig_uniq_transcripts_genes_AB ${sample_name}"

		    label 'process_high'

		    input:
		    set val(sample_name), file("salmon/*") from split_table_uniq_ambig_ab
		    file (annotations) from host_annotations_uniq_ambig_AB


		    output:
		    file "${sample_name}_host_quant_ambig_uniq.sf"
		    file "${sample_name}_pathogen_quant_ambig_uniq.sf"
		    file "${sample_name}_host_quant_ambig_uniq_gene_level.sf"
		    set val(sample_name), file("${sample_name}_host_quant_ambig_uniq.sf") into host_files_comb_uniq_ambig_AB
		    set val(sample_name), file("${sample_name}_pathogen_quant_ambig_uniq.sf") into path_files_comb_uniq_ambig_AB

		    script:
		    """
		    $workflow.projectDir/bin/salmon_extract_ambig_uniq_transcripts_genes.R salmon/*/quant.sf salmon/*/aux_info/ambig_info.tsv $sample_name $annotations
		    """
		}

	    /*
	     * Combine the host ambig and unique counts
	     */
		process host_comb_ambig_uniq_alignment_based {
		    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon_alignment_mode"
		    tag "host_comb_ambig_uniq_AB"

		    label 'process_high'

		    input:
		    file("salmon/*") from host_files_comb_uniq_ambig_AB.collect()

		    output:
		    file "host_quant_combined_ambig_uniq.tsv"

		    script:
		    """
		    $workflow.projectDir/bin/salmon_host_comb_ambig_uniq.R salmon/*/aux_info/*_host_quant_ambig_uniq.sf
		    """
		}

	    /*
	     * Combine the pathogen ambig and unique counts
	     */
		process pathogen_comb_ambig_uniq_alignment_based {
		    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/salmon_alignment_mode"
		    tag "pathogen_comb_ambig_uniq_AB"

		    label 'process_high'

		    input:
		    file("salmon/*") from path_files_comb_uniq_ambig_AB.collect()

		    output:
		    file "pathogen_quant_combined_ambig_uniq.tsv"

		    script:
		    """
		    $workflow.projectDir/bin/salmon_pathogen_comb_ambig_uniq.R salmon/*/aux_info/*_pathogen_quant_ambig_uniq.sf
		    """
		}
	    }


	/*
	 * Tximport - host
	 */

	process tximport_host_salmon_alignment {
	    publishDir "${params.outdir}/salmon_alignment_mode/${sample_name}", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode/${sample_name}"
	    tag "tximport_host"

   	    label 'process_high'

	    input:
	    set val(sample_name), file("salmon/${sample_name}/*") from salmon_alignment_host_tximport
	    file (annotations) from tximport_annotations_salmon_alignment

	    output:
	    file "${sample_name}_host_quant_gene_level.sf" into salmon_files_to_combine_gene_level_alignment

	    script:
	    """
	    $workflow.projectDir/bin/tximport.R salmon $annotations $sample_name
	    """
	}


	/*
	 * Salmon alignment-based mode - Combine host gene level quantification results estimated with Tximport
	 */

	process combine_host_quant_gene_level_salmon_alignment {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "comb_host_quant_genes_sal_alig"

	    label 'process_high'

	    input:
	    file input_quantification from salmon_files_to_combine_gene_level_alignment.collect()

	    output:
	    file "host_combined_gene_level.tsv" into quant_gene_level_host_add_annotations_salmon_alignment

	    script:
	    """
	    python $workflow.projectDir/bin/collect_quantification_data.py -i $input_quantification -q salmon -a gene_id -org host_gene_level
	    """
	}


	/*
	 * Salmon alignment-based mode - Combine quantification results
	 */

	process combine_quantification_tables_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "combine_quantification_salmon"

	    label 'process_high'

	    input:
	    file input_quantification from salmon_files_to_combine_alignment_mode.collect()
	    val gene_attribute from host_atr_collect_data_salmon_alignment_mode

	    output:
	    file "combined_quant.tsv" into split_table_salmon_salmon_alignment

	    script:
	    """
	    python $workflow.projectDir/bin/collect_quantification_data.py -i $input_quantification -q salmon -a $gene_attribute -org both
	    """
	}


	/*
	 * Salmon alignment-based mode - Split quantification tables into host and pathogen results
	 */


	process split_quantification_tables_salmon_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    tag "split_quantification"

	    label 'process_high'

	    input:
	    file quant_table from split_table_salmon_salmon_alignment
	    file transcriptome_pathogen from transcriptome_pathogen_to_split_table_salmon_alignment
	    file transcriptome_host from transcriptome_host_to_split_table_salmon_alignment

	    output:
	    file 'host_quant_salmon.tsv' into host_quantification_mapping_stats_salmon_alignment_based
	    file 'pathogen_quant_salmon.tsv' into pathogen_quantification_mapping_stats_salmon_alignment_based
	    file 'host_quant_salmon.tsv' into host_quantification_RNA_stats_salmon_alignment_based
	    file 'pathogen_quant_salmon.tsv' into pathogen_quantification_RNA_stats_salmon_alignment_based
	    file 'host_quant_salmon.tsv' into quant_host_add_annotations_salmon_alignment_based
	    file 'pathogen_quant_salmon.tsv' into quant_pathogen_add_annotations_alignment_based
	    file 'host_quant_salmon.tsv' into quant_scatter_plot_host_salmon_alignment_based
	    file 'pathogen_quant_salmon.tsv' into quant_scatter_plot_pathogen_salmon_alignment_based
	    env pathonen_tab into pathonen_tab into scatterplots_pathogen_salmon_alignment_based
	    env host_tab into host_tab into scatterplots_host_salmon_alignment_based

            script:
            """
            $workflow.projectDir/bin/split_quant_tables_salmon.sh $transcriptome_pathogen $transcriptome_host $quant_table "quant_salmon.tsv"
            pathonen_tab=\$(if [ \$(cat pathogen_quant_salmon.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
            host_tab=\$(if [ \$(cat host_quant_salmon.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
            """
	}


	/*
	 * Salmon alignment-based mode - Combine pathogen annotations with quantification results
	 */

	process combine_annotations_quant_pathogen_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "comb_annots_quant_pathgn_salmn"

	    label 'process_high'

	    input:
	    file quantification_table from quant_pathogen_add_annotations_alignment_based
	    file annotation_table from annotation_pathogen_combine_quant_salmon_alignment_based
	    val attribute from combine_annot_quant_pathogen_salmon_alignment_based

	    output:
	    file "pathogen_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org pathogen
	    """
	}

	/*
	 * Salmon alignment-based mode - Combine host annotations with quantification results
	 */

	process combine_annotations_quant_host_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "comb_annots_quant_host_salmn"

	    label 'process_high'

	    input:
	    file quantification_table from quant_host_add_annotations_salmon_alignment_based
	    file annotation_table from annotation_host_combine_quant_salmon_alignment_based
	    val attribute from combine_annot_quant_host_salmon_alignment_based

	    output:
	    file "host_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org host
	    """
	}

	/*
	 * Salmon alignment-based mode - Combine host annotations with gene-level estimates
	 */

	process combine_annotations_quant_gene_level_salmon_alignment_mode {
	    publishDir "${params.outdir}/salmon_alignment_mode", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/salmon_alignment_mode"
	    tag "comb_annots_gene_host_salmn"

	    label 'process_high'

	    input:
	    file quantification_table from quant_gene_level_host_add_annotations_salmon_alignment
	    file annotation_table from annotation_host_combine_quant_gene_level_salmon_alignment

	    output:
	    file "host_combined_quant_gene_level_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_annotations_salmon_gene_level.py -q $quantification_table -annotations $annotation_table -a gene_id -org host
	    """
	}


	if(params.mapping_statistics) {

		/*
		 * Extract processed reads from STAR log file
		 */

		process extract_processed_reads_STAR_for_salmon {
		    publishDir "${params.outdir}/mapping_statistics/STAR_for_salmon/processed_reads", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR_for_salmon/processed_reads"
		    tag "extract_processed_reads_STAR"

		    label 'process_high'

		    input:
		    set val(sample_name), file (Log_final_out) from collect_processed_read_counts_STAR_for_salmon

		    output:
		    file "${sample_name}.txt" into collect_results_star_for_salmon

		    script:
		    """
		    $workflow.projectDir/bin/extract_processed_reads.sh $Log_final_out $sample_name star
		    """
		}

		/*
		 * Collect STAR processed reads
		 */

		process collect_processed_reads_STAR_for_salmon {
		    publishDir "${params.outdir}/mapping_statistics/STAR_for_salmon", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR_for_salmon"
		    tag "collect_processed_reads_STAR"

		    label 'process_high'

		    input:
		    file process_reads from collect_results_star_for_salmon.collect()

		    output:
		    file "processed_reads_star.tsv" into mapping_stats_total_processed_reads_alignment_for_salmon

		    script:
		    """
		    cat $process_reads > processed_reads_star.tsv
		    """
		}


		/*
		 * Salmon alignment-based mode  - plot scatter plots of technical replicates for pathogen
		 */

		process scatter_plot_pathogen_salmon_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/scatter_plots"
		    tag "scatter_plots_salmon_pathogen"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_pathogen_salmon_alignment_based
		    val attribute from atr_scatter_plot_pathogen_alignment
		    val replicates from repl_scatter_plots_salmon_alignment_pathogen
		    val pathogen_table_non_empty from scatterplots_pathogen_salmon_alignment_based

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    pathogen_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org pathogen
		    """
		}


		/*
		 * Salmon alignment-based mode  - plot scatter plots of technical replicates for host
		 */

		process scatter_plot_host_salmon_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/scatter_plots"
		    tag "scatter_plots_salmon_host"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_host_salmon_alignment_based
		    val attribute from atr_scatter_plot_host_alignment
		    val replicates from repl_scatter_plots_salmon_alignment_host
		    val host_table_non_empty from scatterplots_host_salmon_alignment_based

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    host_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org host
		    """
		}


		/*
		 * Salmon alignment-based mode  - Extract processed reads from log file
		 */

		process extract_processed_reads_salmon_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based"
		    tag "extract_processed_reads"

		    label 'process_high'

		    input:
		    set val(sample_name), file ("salmon_alignment_mode/*") from collect_processed_read_counts_alignment_based

		    output:
		    file "${sample_name}.txt" into collect_results_alignment_based

		    script:
		    """
		    $workflow.projectDir/bin/extract_processed_reads.sh salmon_alignment_mode/*/aux_info/meta_info.json $sample_name salmon_alignment
		    """
		}


		/*
		 * Salmon alignment-based mode  - Collect processed reads
		 */

		process collect_processed_reads_salmon_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based"
		    tag "collect_processed_reads"

		    label 'process_high'

		    input:
		    file process_reads from collect_results_alignment_based.collect()

		    output:
		    file "processed_reads_salmon_alignment.tsv" into mapping_stats_total_reads_alignment
		    script:
		    """
		    cat $process_reads > processed_reads_salmon_alignment.tsv
		    """
		}


		/*
		 * Salmon alignment-based mode  - Collect mapping stats
		 */

		process salmon_quantification_stats_salmon_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based"
		    tag "quantification_stats_salmon"

		    label 'process_high'

		    input:
		    file quant_table_host from host_quantification_mapping_stats_salmon_alignment_based
		    file quant_table_pathogen from pathogen_quantification_mapping_stats_salmon_alignment_based
		    val attribute from attribute_quant_stats_salmon_alignment
		    file total_processed_reads from mapping_stats_total_reads_alignment
		    file total_processed_reads_star from mapping_stats_total_processed_reads_alignment_for_salmon
		    file total_raw_reads from collect_total_reads_raw_salmon_alignment.ifEmpty('.')

		    output:
		    file ('salmon_alignment_host_pathogen_total_reads.tsv') into salmon_mapped_stats_to_plot_alignment

		    script:
		    """
		    python $workflow.projectDir/bin/mapping_stats.py -q_p $quant_table_pathogen -q_h $quant_table_host -total_processed $total_processed_reads -total_raw $total_raw_reads -a $attribute --star_processed $total_processed_reads_star -t salmon_alignment -o salmon_alignment_host_pathogen_total_reads.tsv
		    """
		}


		/*
		 * Salmon alignment-based mode  - Plot mapping stats
		 */


		process plot_salmon_mapping_stats_host_pathogen_salmon_alignment_based {
		    tag "$name2"
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based"

		    label 'process_high'

		    input:
		    file(stats) from salmon_mapped_stats_to_plot_alignment

		    output:
		    file "*.tsv"
		    file "*.pdf"

		    script:
		    """
		    python $workflow.projectDir/bin/plot_mapping_statistics_salmon_alignment.py -i $stats
		    """
		}


		/*
		 * Salmon alignment-based mode  - Collect pathogen RNA class mapping stats
		 */

		process RNA_class_statistics_salmon_pathogen_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_pathogen", mode: params.publish_dir_mode
		    tag "rna_class_stats_pathogen"

		    label 'process_high'

		    input:
		    file quant_table from pathogen_quantification_RNA_stats_salmon_alignment_based
		    val attribute from host_annotations_RNA_class_stats_pathogen_alignment
		    file gene_annotations from pathogen_annotations_RNA_class_stats_salmon_alignment

		    output:
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen_alignment
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen_combined_alignment
		    file "pathogen_RNA_classes_sum_counts_*.tsv"
		    stdout plot_RNA_stats_pathogen_alignment_boolean
		    stdout plot_RNA_stats_pathogen_combined_alignment_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -q_tool salmon -org pathogen 2>&1
		    '''
		}

		/*
		 * Salmon alignment-based mode  - Collect host RNA class mapping stats
		 */

		process RNA_class_statistics_salmon_host_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_host", mode: params.publish_dir_mode
		    tag "rna_class_stats_host"

		    label 'process_high'

		    input:
		    file quant_table from host_quantification_RNA_stats_salmon_alignment_based
		    val attribute from attribute_host_RNA_class_stats_alignment
		    file gene_annotations from host_annotations_RNA_class_stats_salmon_alignment
		    file rna_classes_to_replace from RNA_classes_to_replace_alignment

		    output:
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_alignment
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_combined_alignment
		    file "host_RNA_classes_sum_counts_*.tsv"
		    file "host_gene_types_groups_*"
		    stdout plot_RNA_stats_host_alignment_boolean
		    stdout plot_RNA_stats_host_combined_alignment_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -rna !{rna_classes_to_replace} -q_tool salmon -org host 2>&1
		    '''
		}


		/*
		 * Salmon alignment-based mode  - Plot pathogen RNA class mapping stats for each sample separately
		 */

		process plot_RNA_class_salmon_pathogen_each_alignment_based{
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_pathogen"
		    tag "plot_rna_class_stats_path_sal"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen_alignment
		    val plot_rna from plot_RNA_stats_pathogen_alignment_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}

		/*
		 * Salmon alignment-based mode  - Plot host RNA class mapping stats for each sample separately
		 */

		process plot_RNA_class_salmon_host_each_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_host"
		    tag "plot_rna_class_stats_host_sal"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host_alignment
		    val plot_rna from plot_RNA_stats_host_alignment_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}


		/*
		 * Salmon alignment-based mode  - Plot pathogen RNA class mapping stats for all samples
		 */

		process plot_RNA_class_salmon_pathogen_combined_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_pathogen"
		    tag "plot_rna_class_stats_path_all"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen_combined_alignment
		    val plot_rna from plot_RNA_stats_pathogen_combined_alignment_boolean

		    output:
		    file "RNA_class_stats_combined_pathogen.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org pathogen
		    """
		}

		/*
		 * Salmon alignment-based mode  - Plot host RNA class mapping stats for all samples
		 */

		process plot_RNA_class_salmon_host_combined_alignment_based {
		    publishDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/salmon_alignment_based/RNA_classes_host"
		    tag "plot_rna_class_stats_host_all"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host_combined_alignment
		    val plot_rna from plot_RNA_stats_host_combined_alignment_boolean

		    output:
		    file "RNA_class_stats_combined_host.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org host
		    """
		}
	}
}else{
   Channel.empty()
     .into {multiqc_star_for_salmon_alignment; multiqc_salmon_alignment_quant}
}



/*
 * STEP 8 - STAR
 */


if(params.run_star) {

	/*
	 * STAR - build an index
	 */

	process STARindex {
      publishDir "${params.outdir}/STAR", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/STAR"
	    tag "build_star_index"


	    label 'process_high'

	    input:
	    file(fasta) from host_pathogen_fasta_star_index
	    file(gff) from gff_host_star_alignment_gff

	    output:
      file "index/*" into star_index_genome_alignment

	    script:
	    sjdbOverhang = params.sjdbOverhang
	    star_index_params = params.star_index_params
	    """
	    mkdir index
	    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir index/ --genomeFastaFiles $fasta --sjdbGTFfile $gff --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang $sjdbOverhang $star_index_params
	    """
	}

	/*
	 * STAR - alignment
	 */

	process ALIGNMENTstar {
	    tag "${sample_name}"
      publishDir "${params.outdir}/STAR", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/STAR"

      label 'process_high'

	    input:
	    set val(sample_name),file(reads) from  trimming_results_star_htseq
	    file(gff) from gff_host_star_htseq_alignment_gff.collect()
	    file(index) from star_index_genome_alignment.collect()

	    output:
	    set val(sample_name), file("${sample_name}/${sample_name}Aligned*.out.bam") into star_aligned_u_m
	    set val(sample_name), file("${sample_name}/${sample_name}Aligned*.out.bam") into alignment_unique_mapping_stats
	    set val(sample_name), file("${sample_name}/${sample_name}Aligned*.out.bam") into alignment_crossmapped_extract
	    file "${sample_name}/*" into multiqc_star_alignment
	    set val(sample_name), file("${sample_name}/${sample_name}Log.final.out") into collect_processed_read_counts_STAR

	    script:
	    outSAMunmapped = params.outSAMunmapped
	    outSAMattributes = params.outSAMattributes
	    outFilterMultimapNmax = params.outFilterMultimapNmax
	    outFilterType = params.outFilterType
	    alignSJoverhangMin = params.alignSJoverhangMin
	    alignSJDBoverhangMin = params.alignSJDBoverhangMin
	    outFilterMismatchNmax = params.outFilterMismatchNmax
	    outFilterMismatchNoverReadLmax = params.outFilterMismatchNoverReadLmax
	    alignIntronMin = params.alignIntronMin
	    alignIntronMax = params.alignIntronMax
	    alignMatesGapMax = params.alignMatesGapMax
	    outWigType = params.outWigType
	    outWigStrand = params.outWigStrand
	    limitBAMsortRAM = params.limitBAMsortRAM
	    readFilesCommand = reads[0].toString().endsWith('.gz') ? "--readFilesCommand zcat" : ''
	    winAnchorMultimapNmax = params.winAnchorMultimapNmax
	    star_alignment_params = params.star_alignment_params

	    if (params.single_end){
	    """
	    	mkdir $sample_name
	    	STAR --runThreadN ${task.cpus} --genomeDir . --sjdbGTFfile $gff $readFilesCommand --readFilesIn $reads --outSAMtype BAM SortedByCoordinate --outSAMunmapped $outSAMunmapped --outSAMattributes $outSAMattributes --outWigType $outWigType --outWigStrand $outWigStrand --outFileNamePrefix $sample_name/$sample_name --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFilterMultimapNmax $outFilterMultimapNmax --outFilterType $outFilterType --limitBAMsortRAM $limitBAMsortRAM --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --winAnchorMultimapNmax $winAnchorMultimapNmax $star_alignment_params
	    """
	    } else {
	    """
	    mkdir $sample_name
	    STAR --runThreadN ${task.cpus} --genomeDir . --sjdbGTFfile $gff $readFilesCommand --readFilesIn ${reads[0]} ${reads[1]} --outSAMtype BAM SortedByCoordinate --outSAMunmapped $outSAMunmapped --outSAMattributes $outSAMattributes --outWigType $outWigType --outWigStrand $outWigStrand --outFileNamePrefix $sample_name/$sample_name --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript Parent --outFilterMultimapNmax $outFilterMultimapNmax --outFilterType $outFilterType --limitBAMsortRAM $limitBAMsortRAM --alignSJoverhangMin $alignSJoverhangMin --alignSJDBoverhangMin $alignSJDBoverhangMin --outFilterMismatchNmax $outFilterMismatchNmax --outFilterMismatchNoverReadLmax $outFilterMismatchNoverReadLmax --alignIntronMin $alignIntronMin --alignIntronMax $alignIntronMax --alignMatesGapMax $alignMatesGapMax --winAnchorMultimapNmax $winAnchorMultimapNmax $star_alignment_params
	    """
	    }
	}



	if(params.mapping_statistics) {

		/*
		 * STAR - Remove cross mapped reads from bam file
		 */

		process remove_crossmapped_reads {
		    tag "$sample_name"
		    publishDir "${params.outdir}/STAR/multimapped_reads", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/STAR/multimapped_reads"

                    label 'process_high'

		    input:
		    set val(sample_name), file(alignment) from alignment_crossmapped_extract
		    file(host_reference) from reference_host_names_crossmapped_find.collect()
		    file(pathogen_reference) from reference_pathogen_crossmapped_find.collect()


		    output:
		    set val(sample_name), file("${bam_file_without_crossmapped}") into alignment_multi_mapping_stats
		    file "${cross_mapped_reads}" into count_crossmapped_reads

		    script:
		    bam_file_without_crossmapped = sample_name + "_no_crossmapped.bam"
		    cross_mapped_reads = sample_name + "_cross_mapped_reads.txt"

		    if (params.single_end){
		    	"""
		    	$workflow.projectDir/bin/remove_crossmapped_reads_BAM.sh $alignment $workflow.projectDir/bin $host_reference $pathogen_reference $cross_mapped_reads $bam_file_without_crossmapped
		    	"""
		    } else {
		    	"""
		    	$workflow.projectDir/bin/remove_crossmapped_read_pairs_BAM.sh $alignment $workflow.projectDir/bin $host_reference $pathogen_reference $cross_mapped_reads $bam_file_without_crossmapped
		    	"""
		    }
		}


		/*
		 * STAR - Extract processed reads from log file
		 */

		process extract_processed_reads_STAR {
		    publishDir "${params.outdir}/mapping_statistics/STAR/processed_reads", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR/processed_reads"
		    tag "extract_processed_reads_STAR"

	   	    label 'process_high'

		    input:
		    set val(sample_name), file (Log_final_out) from collect_processed_read_counts_STAR

		    output:
		    file "${sample_name}.txt" into collect_results_star

		    script:
		    """
		    $workflow.projectDir/bin/extract_processed_reads.sh $Log_final_out $sample_name star
		    """
		}


		/*
		 * STAR - Collect STAR processed reads
		 */

		process collect_processed_reads_STAR {
		    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR"
		    tag "collect_processed_reads_STAR"

		    label 'process_high'

		    input:
		    file process_reads from collect_results_star.collect()

		    output:
		    file "processed_reads_star.tsv" into mapping_stats_total_processed_reads_alignment
		    file "processed_reads_star.tsv" into mapping_stats_htseq_total_processed_reads_alignment

		    script:
		    """
		    cat $process_reads > processed_reads_star.tsv
		    """
		}


		/*
		 * STAR - Count uniquely mapped reads
		 */

		process unique_mapping_stats_STAR {
		    tag "$sample_name"
		    publishDir "${params.outdir}/mapping_statistics/STAR/uniquely_mapped", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR/uniquely_mapped"

		    label 'process_high'

		    input:
		    set val(sample_name), file(alignment) from alignment_unique_mapping_stats
		    file(host_reference_names) from reference_host_names_uniquelymapped.collect()
		    file(pathogen_reference_names) from reference_pathogen_names_uniquelymapped.collect()

		    output:
		    file("${name}") into STAR_mapping_stats_unique

		    shell:
		    name = sample_name + '_uniquely_mapped.txt'
		    if (params.single_end){
		    '''
		    !{workflow.projectDir}/bin/count_uniquely_mapped_reads.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
		    '''
		    } else {
		    '''
		    !{workflow.projectDir}/bin/count_uniquely_mapped_read_pairs.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
		    '''
		    }
		}

		/*
		 * STAR - Collect uniquely mapped reads
		 */

		process collect_stats_STAR_uniquely_mapped {
		    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR"
		    tag "collect_uniq_mapped_reads_STAR"

		    label 'process_high'

		    input:
		    file stats from STAR_mapping_stats_unique.collect()

		    output:
		    file "uniquely_mapped_reads_star.tsv" into mapping_stats_uniquely_mapped_star

		    script:
		    """
		    python $workflow.projectDir/bin/combine_tables.py -i $stats -o uniquely_mapped_reads_star.tsv -s uniquely_mapped_reads
		    """
		}


		/*
		 * STAR - Count cross mapped reads
		 */

		process count_crossmapped_reads {
		    tag "count_crossmapped_reads"
		    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR"

		    label 'process_high'

		    input:
		    file(cross_mapped_reads) from count_crossmapped_reads.collect()

		    output:
                    file "cross_mapped_reads_sum.txt" into STAR_mapping_stats_cross_mapped

		    script:
		    """
		    $workflow.projectDir/bin/count_cross_mapped_reads.sh $cross_mapped_reads
		    """
		}


		/*
		 * STAR - Count multi-mapped reads (multi_mapped reads without cross_mapped reads )
		 */

		process multi_mapping_stats {
		    tag "$sample_name"
		    publishDir "${params.outdir}/mapping_statistics/STAR/multi_mapped", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR/multi_mapped"

		    label 'process_high'

		    input:
		    set val(sample_name),file(alignment) from alignment_multi_mapping_stats
		    file(host_reference_names) from reference_host_names_multimapped.collect()
		    file(pathogen_reference_names) from reference_pathogen_names_multimapped.collect()

		    output:
		    file("${name}") into STAR_mapping_stats_multi

		    shell:
		    name = sample_name + '_multi_mapped.txt'
		    if (params.single_end){
		    '''
		    !{workflow.projectDir}/bin/count_multi_mapped_reads.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
		    '''
		    } else {
		    '''
		    !{workflow.projectDir}/bin/count_multi_mapped_read_pairs.sh !{alignment} !{host_reference_names} !{pathogen_reference_names} !{sample_name} !{name}
		    '''
		    }
		}


		/*
		 * STAR - Collect multi-mapped reads
		 */

		process collect_stats_STAR_multi_mapped {
			    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
			    storeDir "${params.outdir}/mapping_statistics/STAR"
			    tag "collect_multi_mapped_reads_STAR"

			    label 'process_high'

			    input:
			    file stats from STAR_mapping_stats_multi.collect()

			    output:
			    file "multi_mapped_reads_star.tsv" into mapping_stats_multi_mapped_star

			    script:
			    """
			    python $workflow.projectDir/bin/combine_tables.py -i $stats -o multi_mapped_reads_star.tsv -s multi_mapped_reads
			    """
			}


		/*
		 * STAR - Collect mapping stats
		 */

		process star_mapping_stats {
		    storeDir "${params.outdir}/mapping_statistics/STAR"
		    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
		    tag "star_mapping_stats"

		    label 'process_high'

		    input:
		    file total_raw_reads from collect_total_reads_raw_star.ifEmpty('.')
		    file total_processed_reads from mapping_stats_total_processed_reads_alignment
		    file uniquely_mapped_reads from mapping_stats_uniquely_mapped_star
		    file multi_mapped_reads from mapping_stats_multi_mapped_star
		    file cross_mapped_reads from STAR_mapping_stats_cross_mapped

		    output:
		    file ('star_mapping_stats.tsv') into star_mapped_stats_to_plot
		    file ('star_mapping_stats.tsv') into mapping_stats_star_htseq_stats

		    script:
		    """
		    python $workflow.projectDir/bin/mapping_stats.py -total_raw $total_raw_reads -total_processed $total_processed_reads -m_u $uniquely_mapped_reads -m_m $multi_mapped_reads -c_m $cross_mapped_reads -t star -o star_mapping_stats.tsv
		    """
		}

		/*
		 * STAR - Plot mapping stats
		 */

		process plot_star_mapping_stats {
		    tag "plot_star_mapping_stats"
		    publishDir "${params.outdir}/mapping_statistics/STAR", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/STAR"

		    label 'process_high'

		    input:
		    file(stats) from star_mapped_stats_to_plot

		    output:
		    file "*.tsv"
		    file "*.pdf"

		    script:
		    """
		    python $workflow.projectDir/bin/plot_mapping_stats_star.py -i $stats
		    """
		}
	}
}else{
   Channel.empty()
     .set {multiqc_star_alignment}
}



/*
 * STEP 9 - HTSeq
 */


if(params.run_htseq_uniquely_mapped){

	/*
	 * HTSeq - quantification
	 */

	process HTseq_unique_mapping {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/HTSeq"
	    tag "$sample_name"

            label 'process_high'

	    input:
	    file(gff) from quantification_gff_u_m.collect()
	    set val(sample_name), file(st) from star_aligned_u_m
	    val(host_attribute) from host_gff_attribute_htseq
	    val(stranded) from stranded_htseq_unique


	    output:
	    file ("$name_file2") into htseq_files_to_combine
	    file ("$name_file2") into multiqc_htseq_unique

	    script:
	    name_file2 = sample_name + "_count_u_m.txt"
	    host_attr = host_attribute
	    max_reads_in_buffer = params.max_reads_in_buffer
	    minaqual = params.minaqual
	    htseq_params = params.htseq_params
	    """
	    htseq-count -n $task.cpus -t quant -f bam -r pos $st $gff -i $host_attr -s $stranded --max-reads-in-buffer=$max_reads_in_buffer -a $minaqual $htseq_params > $name_file2
	    sed -i '1{h;s/.*/'"$sample_name"'/;G}' "$name_file2"
	    """
	}


	/*
	 * HTSeq - combine quantification results
	 */

	htseq_files_to_combine
	    .collect()
	    .set { htseq_result_quantification }

	process combine_quantification_tables_htseq_uniquely_mapped {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/HTSeq"
	    tag "comb_quants_htseq_uniq_mapped"

            label 'process_high'

	    input:
	    file input_quantification from htseq_result_quantification
	    val(host_attribute) from host_gff_attribute_htseq_combine

	    output:
	    file "quantification_results_*.tsv" into htseq_result_quantification_TPM

	    script:
	    """
	    python $workflow.projectDir/bin/collect_quantification_data.py -i $input_quantification -q htseq -a $host_attribute
	    """
	}


	/*
	 * HTSeq - Calculate TPM values
	 */

	process htseq_uniquely_mapped_TPM {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/HTSeq"
	    tag "htseq_uniquely_mapped_TPM"

            label 'process_high'

	    input:
	    file input_quantification from htseq_result_quantification_TPM
	    val(host_attribute) from host_gff_attribute_htseq_TPM
	    file gff_host from gff_host_to_TPM
	    file gff_pathogen from gff_pathogen_to_TPM

	    output:
	    file "quantification_results_uniquely_mapped_NumReads_TPM.tsv" into split_table_htseq_host
	    file "quantification_results_uniquely_mapped_NumReads_TPM.tsv" into split_table_htseq_pathogen

	    script:
	    """
	    $workflow.projectDir/bin/calculate_TPM_HTSeq.R $input_quantification $host_attribute $gff_pathogen $gff_host
	    """
	}


	/*
	 * HTSeq - Split quantification tables into host and pathogen results
	 */


	process split_quantification_tables_htseq_uniquely_mapped {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    tag "split_quants_uniq_mapped_host"

            label 'process_high'

	    input:
	    file quant_table from split_table_htseq_host
	    file host_annotations from annotation_host_split_quant_htseq
	    file pathogen_annotations from annotation_pathogen_split_quant_htseq

	    output:
	    file 'host_quantification_uniquely_mapped_htseq.tsv' into host_quantification_stats_htseq
	    file 'host_quantification_uniquely_mapped_htseq.tsv' into host_quantification_stats_htseq_total
	    file 'host_quantification_uniquely_mapped_htseq.tsv' into host_htseq_quantification_RNA_stats
	    file 'host_quantification_uniquely_mapped_htseq.tsv' into quant_host_add_annotations_htseq_u_m
	    file 'host_quantification_uniquely_mapped_htseq.tsv' into quant_scatter_plot_host_htseq_u_m
	    file 'pathogen_quantification_uniquely_mapped_htseq.tsv' into pathogen_quantification_stats_htseq
	    file 'pathogen_quantification_uniquely_mapped_htseq.tsv' into pathogen_quantification_stats_htseq_total
	    file 'pathogen_quantification_uniquely_mapped_htseq.tsv' into pathogen_htseq_quantification_RNA_stats
	    file 'pathogen_quantification_uniquely_mapped_htseq.tsv' into quant_pathogen_add_annotations_htseq_u_m
	    file 'pathogen_quantification_uniquely_mapped_htseq.tsv' into quant_scatter_plot_pathogen_htseq_u_m
	    env pathonen_tab into pathonen_tab into scatterplots_pathogen_htseq
	    env host_tab into host_tab into scatterplots_host_htseq

	    script:
	    """
	    $workflow.projectDir/bin/split_quant_tables.sh $quant_table $host_annotations $pathogen_annotations quantification_uniquely_mapped_htseq.tsv
            pathonen_tab=\$(if [ \$(cat pathogen_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
            host_tab=\$(if [ \$(cat host_quantification_uniquely_mapped_htseq.tsv | wc -l) -gt 1  ]; then echo "true"; else echo "false"; fi)
	    """
	}


	/*
	 * HTSeq - Combine pathogen annotations with quantification results
	 */

	process combine_annotations_quant_pathogen_uniquely_mapped_host {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/HTSeq"
	    tag "comb_annots_quant_pathogen"

	    label 'process_high'

	    input:
	    file quantification_table from quant_pathogen_add_annotations_htseq_u_m
	    file annotation_table from annotation_pathogen_combine_quant_htseq_u_m
	    val attribute from combine_annot_quant_pathogen_host_gff_attribute

	    output:
	    file "pathogen_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org pathogen
	    """
	}


	/*
	 * HTSeq - Combine host annotations with quantification results
	 */

	process combine_annotations_quant_host_uniquely_mapped_host {
	    publishDir "${params.outdir}/HTSeq", mode: params.publish_dir_mode
	    storeDir "${params.outdir}/HTSeq"
	    tag "comb_annots_quant_host"

	    label 'process_high'

	    input:
	    file quantification_table from quant_host_add_annotations_htseq_u_m
	    file annotation_table from annotation_host_combine_quant_htseq
	    val attribute from combine_annot_quant_pathogen_host_gff_attribute

	    output:
	    file "host_combined_quant_annotations.tsv"

	    script:
	    """
	    $workflow.projectDir/bin/combine_quant_annotations.py -q $quantification_table -annotations $annotation_table -a $attribute -org host
	    """
	}


	if(params.mapping_statistics) {

		/*
	 	 * HTSeq - Plot scatter plots of technical replicates for pathogen
	 	 */

		process scatter_plot_pathogen_htseq {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots"
		    tag "scatter_plot_pathogen_htseq"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_pathogen_htseq_u_m
		    val attribute from atr_scatter_plot_pathogen_htseq_u_m
		    val replicates from repl_scatter_plots_htseq_pathogen
		    val pathogen_table_non_empty from scatterplots_pathogen_htseq

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    pathogen_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org pathogen
		    """
		}


		/*
	 	 * HTSeq - Plot scatter plots of technical replicates for host
	 	 */


		process scatter_plot_host_htseq {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/scatter_plots"
		    tag "scatter_plot_host_htseq"

		    label 'process_high'

		    input:
		    file quant_table from quant_scatter_plot_host_htseq_u_m
		    val attribute from atr_scatter_plot_host_htseq_u_m
		    val replicates from repl_scatter_plots_htseq_host
		    val host_table_non_empty from scatterplots_host_htseq

		    output:
		    file ('*.pdf')

		    when:
		    replicates.toBoolean()
		    host_table_non_empty.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/scatter_plots.py -q $quant_table -a $attribute -org host
		    """
		}


		/*
	 	 * HTSeq - Collect mapping and quantification stats
	 	 */


		process htseq_quantification_stats_uniquely_mapped {
		    storeDir "${params.outdir}/mapping_statistics/HTSeq"
		    publishDir "${params.outdir}/mapping_statistics/HTSeq", mode: params.publish_dir_mode
		    tag "quantification_stats_htseq"

		    label 'process_high'

		    input:
		    file quant_table_host from host_quantification_stats_htseq_total
		    file quant_table_pathogen from pathogen_quantification_stats_htseq_total
		    val attribute from host_gff_attribute_mapping_stats_htseq
		    file star_stats from mapping_stats_star_htseq_stats

		    output:
		    file ('htseq_uniquely_mapped_reads_stats.tsv') into htseq_mapped_stats_to_plot

		    script:
		    """
		    python $workflow.projectDir/bin/mapping_stats.py -q_p $quant_table_pathogen -q_h $quant_table_host -a $attribute  -star $star_stats -t htseq -o htseq_uniquely_mapped_reads_stats.tsv
		    """
		}


		/*
	 	 * HTSeq - Plot mapping and stats
	 	 */

		process plot_mapping_stats_host_pathogen_htseq_uniquely_mapped{
		    tag "$name2"
		    publishDir "${params.outdir}/mapping_statistics/HTSeq", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq"

		    label 'process_high'

		    input:
		    file(stats) from htseq_mapped_stats_to_plot

		    output:
		    file "*.tsv"
		    file "*.pdf"

		    script:
		    """
		    python $workflow.projectDir/bin/plot_mapping_stats_htseq.py -i $stats
		    """
		}


		/*
	 	 * HTSeq - Collect pathogen RNA class mapping stats
	 	 */

		process RNA_class_statistics_htseq_uniquely_mapped_pathogen {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_pathogen", mode: params.publish_dir_mode
		    tag "rna_class_stats_htseq_pathogen"

		    label 'process_high'

		    input:
		    file quant_table from pathogen_htseq_quantification_RNA_stats
		    val attribute from host_gff_attribute_RNA_class_pathogen_htseq
		    file gene_annotations from pathogen_annotations_RNA_class_stats_htseq


		    output:
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen_htseq_u_m
		    file "pathogen_RNA_classes_percentage_*.tsv" into plot_RNA_stats_pathogen_combined_htseq_u_m
		    file "pathogen_RNA_classes_sum_counts_*.tsv"
		    stdout plot_RNA_stats_pathogen_htseq_u_m_boolean
		    stdout plot_RNA_stats_pathogen_combined_htseq_u_m_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -q_tool htseq -org pathogen 2>&1
		    '''
		}


		/*
	 	 * HTSeq - Collect host RNA class mapping stats
	 	 */

		process RNA_class_statistics_htseq_uniquely_mapped_host {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host", mode: params.publish_dir_mode
		    tag "rna_class_stats_htseq_host"

		    label 'process_high'

		    input:
		    file quant_table from host_htseq_quantification_RNA_stats
		    val attribute from host_gff_attribute_RNA_class_host_htseq
		    file gene_annotations from host_annotations_RNA_class_stats_htseq
		    file rna_classes_to_replace from RNA_classes_to_replace_htseq_uniquely_mapped

		    output:
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_htseq_u_m
		    file "host_RNA_classes_percentage_*.tsv" into plot_RNA_stats_host_combined_htseq_u_m
		    file "host_RNA_classes_sum_counts_*.tsv"
		    file "host_gene_types_groups_*"
		    stdout plot_RNA_stats_host_htseq_u_m_boolean
		    stdout plot_RNA_stats_host_combined_htseq_u_m_boolean

		    shell:
		    '''
		    python !{workflow.projectDir}/bin/RNA_class_content.py -q !{quant_table} -a !{attribute} -annotations !{gene_annotations} -rna !{rna_classes_to_replace} -q_tool htseq -org host 2>&1
		    '''
		}


		/*
	 	 * HTSeq - Plot pathogen RNA class mapping stats for each sample separately
	 	 */

		process plot_RNA_class_htseq_uniquely_mapped_pathogen_each{
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_pathogen"
		    tag "plot_rna_stats_htseq_pathogen"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen_htseq_u_m
		    val plot_rna from plot_RNA_stats_pathogen_htseq_u_m_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}


		/*
	 	 * HTSeq - Plot host RNA class mapping stats for each sample separately
	 	 */

		process plot_RNA_class_htseq_uniquely_mapped_host_each {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host"
		    tag "plot_rna_stats_htseq_host"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host_htseq_u_m
		    val plot_rna from plot_RNA_stats_host_htseq_u_m_boolean

		    output:
		    file "*.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_each.py -i $stats_table
		    """
		}


		/*
	 	 * HTSeq - Plot pathogen RNA class mapping stats for all samples
	 	 */

		process plot_RNA_class_htseq_uniquely_mapped_pathogen_combined {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_pathogen", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_pathogen"
		    tag "plt_rna_stats_htseq_pathgn_all"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_pathogen_combined_htseq_u_m
		    val plot_rna from plot_RNA_stats_pathogen_combined_htseq_u_m_boolean

		    output:
		    file "RNA_class_stats_combined_pathogen.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org pathogen
		    """
		}

		/*
	 	 * HTSeq - Plot host RNA class mapping stats for all samples
	 	 */

		process plot_RNA_class_htseq_uniquely_host_combined {
		    publishDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host", mode: params.publish_dir_mode
		    storeDir "${params.outdir}/mapping_statistics/HTSeq/RNA_classes_host"
		    tag "plt_rna_stats_htseq_host_all"

		    label 'process_high'

		    input:
		    file stats_table from plot_RNA_stats_host_combined_htseq_u_m
		    val plot_rna from plot_RNA_stats_host_combined_htseq_u_m_boolean

		    output:
		    file "RNA_class_stats_combined_host.pdf"

		    when:
		    plot_rna.toBoolean()

		    script:
		    """
		    python $workflow.projectDir/bin/plot_RNA_class_stats_combined.py -i $stats_table -org host
		    """
		}

}
}else{
   Channel.empty()
     .set {multiqc_htseq_unique}
}




process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: params.publish_dir_mode

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from ch_fastqc_results.last().collect().ifEmpty([])
    file ('fastqc_after_trimming/*') from ch_fastqc_trimmed_results.last().collect().ifEmpty([])
    file ('salmon/*') from multiqc_salmon_quant.collect().ifEmpty([])
    file ('salmon_alignment_mode/*') from multiqc_salmon_alignment_quant.collect().ifEmpty([])
    file ('STAR/*') from multiqc_star_alignment.collect().ifEmpty([])
    file ('STAR_for_salmon/*') from multiqc_star_for_salmon_alignment.collect().ifEmpty([])
    file ('uniquely_mapped/*') from multiqc_htseq_unique.collect().ifEmpty([])
    file ('software_versions/*') from ch_software_versions_yaml.collect()
    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -d --export -f $rtitle $rfilename $custom_config_file .
    """
}


/*
 * STEP 11 - Output Description HTML
 */


process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs

    output:
    file "results_description.html"
//    file images from ch_output_docs_images

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}






/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/dualrnaseq] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/dualrnaseq] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/dualrnaseq] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/dualrnaseq] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/dualrnaseq] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html
            log.info "[nf-core/dualrnaseq] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/dualrnaseq]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/dualrnaseq]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "${c_red}====================================================${c_reset}\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "${c_red}====================================================${c_reset}\n"
                }
            }
        }
    }
}
