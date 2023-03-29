/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDualrnaseq.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta_host ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { SALMON_SELECTIVE_ALIGNMENT } from '../subworkflows/local/salmon_selective_alignment'
include { SALMON_ALIGNMENT_BASE } from '../subworkflows/local/salmon_alignment_base'
include { EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_SALMON;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_HOST_HTSEQ;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_SALMON;
    EXTRACT_ANNOTATIONS as EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ
    } from '../modules/local/extract_annotations/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                            } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_AFTER_TRIMMING   } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                          } from '../modules/nf-core/cutadapt/main'
include { MULTIQC                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DUALRNASEQ {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    ch_gff_host = Channel.fromPath(params.gff_host, checkIfExists: true)
    EXTRACT_ANNOTATIONS_HOST_SALMON (
        ch_gff_host,
        params.extract_annotations_host_salmon_feature,
        params.extract_annotations_host_salmon_attribute,
        params.extract_annotations_host_salmon_organism,
        'salmon'
    )

	ch_gene_feature_gff_to_quantify_host = Channel
	    .value(params.gene_feature_gff_to_quantify_host)
	    .collect()
    EXTRACT_ANNOTATIONS_HOST_HTSEQ (
        ch_gff_host,
        ch_gene_feature_gff_to_quantify_host,
        params.host_gff_attribute,
        params.extract_annotations_host_htseq_organism,
        'htseq'
    )

    ch_gff_pathogen = Channel.fromPath(params.gff_pathogen, checkIfExists: true)
	ch_gene_feature_gff_to_quantify_pathogen = Channel
	    .value(params.gene_feature_gff_to_quantify_pathogen)
	    .collect()
    EXTRACT_ANNOTATIONS_PATHOGEN_HTSEQ (
        ch_gff_pathogen,
        ch_gene_feature_gff_to_quantify_pathogen,
        params.pathogen_gff_attribute,
        params.extract_annotations_pathogen_htseq_organism,
        'htseq'
    )

	ch_gene_feature_gff_to_create_transcriptome_pathogen = Channel
	    .value(params.gene_feature_gff_to_create_transcriptome_pathogen)
	    .collect()
    EXTRACT_ANNOTATIONS_PATHOGEN_SALMON (
        ch_gff_pathogen,
        ch_gene_feature_gff_to_create_transcriptome_pathogen,
        params.pathogen_gff_attribute,
        params.extract_annotations_pathogen_salmon_organism,
        'salmon'
    )

    if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(INPUT_CHECK.out.reads)
            ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    if (!(params.skip_tools && params.skip_tools.split(',').contains('cutadapt'))) {
            CUTADAPT(INPUT_CHECK.out.reads)
            ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
            FASTQC_AFTER_TRIMMING(CUTADAPT.out.reads)
            ch_versions = ch_versions.mix(FASTQC_AFTER_TRIMMING.out.versions.first())
    }


    //
    // SUBWORKFLOW: Create salmon index and run the quantification
    //
    // for testing purposes use only host transcript_fasta; chimeric transcript fasta should be an input
    params.transcript_fasta = params.transcript_fasta_host

    ch_genome_fasta                     = Channel.fromPath(params.fasta_host, checkIfExists: true)
    ch_transcript_fasta                 = Channel.fromPath(params.transcript_fasta, checkIfExists: true)
    ch_transcript_fasta_pathogen        = Channel.fromPath(params.transcript_fasta_pathogen, checkIfExists: true)
    ch_transcript_fasta_host            = Channel.fromPath(params.transcript_fasta_host, checkIfExists: true)
    // TODO change to gff in the future
    ch_gtf                              = Channel.fromPath(params.gff_host, checkIfExists: true)

    if ( params.run_salmon_selective_alignment ) {
        SALMON_SELECTIVE_ALIGNMENT (
            INPUT_CHECK.out.reads,
            ch_genome_fasta,
            ch_transcript_fasta,
            ch_gtf,
            ch_transcript_fasta_pathogen,
            ch_transcript_fasta_host
        )
        ch_versions = ch_versions.mix(SALMON_SELECTIVE_ALIGNMENT.out.versions)
    }

    if ( params.run_salmon_alignment_based_mode ) {
        SALMON_ALIGNMENT_BASE (
            INPUT_CHECK.out.reads,
            ch_genome_fasta,
            ch_transcript_fasta,
            ch_gtf
        )
        ch_versions = ch_versions.mix(SALMON_ALIGNMENT_BASE.out.versions)
    }



    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDualrnaseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDualrnaseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
