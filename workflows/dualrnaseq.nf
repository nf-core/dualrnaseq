/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

//subworkflow and module inclusion
include { INPUT_CHECK                     } from '../subworkflows/local/input_check'
include { PREPARE_REFERENCE_FILES         } from '../subworkflows/local/prepare_reference_files'
include { SALMON_SELECTIVE_ALIGNMENT      } from '../subworkflows/local/salmon_selective_alignment'
include { SALMON_ALIGNMENT_BASED          } from '../subworkflows/local/salmon_alignment_based'
include { STAR_HTSEQ as STAR_ALIGNMENT    } from '../subworkflows/local/star_htseq'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                          } from '../modules/nf-core/fastqc'
include { FASTQC as FASTQC_AFTER_TRIMMING } from '../modules/nf-core/fastqc'
include { CUTADAPT                        } from '../modules/nf-core/cutadapt'
include { MULTIQC                         } from '../modules/nf-core/multiqc'
include { CUSTOM_DUMPSOFTWAREVERSIONS     } from '../modules/nf-core/custom/dumpsoftwareversions'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DUALRNASEQ {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //salmon_sa_out = Channel.empty()
    salmon_ab_out = Channel.empty()

    // Initialize required channels
    ch_workflow_summary = Channel.empty()
    ch_methods_description = Channel.empty()
    ch_multiqc_config = Channel.fromPath("${projectDir}/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo) : Channel.empty()




    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK(ch_samplesheet)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // Store input reads
    ch_reads = INPUT_CHECK.out.reads.map { meta, reads -> tuple(meta, reads) }



    // if skip_tools passed, but not contain fastqc
    if (!(params.skip_tools && params.skip_tools.split(',').contains('fastqc'))) {
        FASTQC(INPUT_CHECK.out.reads)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    // if skip_tools passed, but not contain cutadapt
    if (!(params.skip_tools && params.skip_tools.split(',').contains('cutadapt'))) {
        CUTADAPT(ch_reads)
        ch_reads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())
    }

    // if skip_tools passed, but not contain fastqc and cutadapt - so should run fastqc after trimming
    if (!(params.skip_tools && (params.skip_tools.split(',').contains('fastqc') || params.skip_tools.split(',').contains('cutadapt')))) {
        FASTQC_AFTER_TRIMMING(ch_reads)
        ch_versions = ch_versions.mix(FASTQC_AFTER_TRIMMING.out.versions.first())
    }


    // ---------------
    // Prepare reference files
    // ---------------
    // uncompress files, merge host and pathogen files together,
    // update features in reference files to be compatable with software
    PREPARE_REFERENCE_FILES(
        params.host_fasta_genome,
        params.host_gff,
        params.pathogen_fasta_genome,
        params.pathogen_gff,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] })
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    // Run Salmon selective alighment
    if (params.run_salmon_SA) {
        SALMON_SELECTIVE_ALIGNMENT(
            ch_reads,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_genome,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_pathogen_transcripts_gff,
            PREPARE_REFERENCE_FILES.out.pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.annotations_host_salmon,
        )
        ch_versions = ch_versions.mix(SALMON_SELECTIVE_ALIGNMENT.out.versions)
        //salmon_sa_out = SALMON_SELECTIVE_ALIGNMENT.out
    }

    // Run Salmon alignment based
    if (params.run_salmon_AB) {
        SALMON_ALIGNMENT_BASED(
            ch_reads,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_genome,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_pathogen_transcripts_gff,
            PREPARE_REFERENCE_FILES.out.pathogen_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.host_fasta_transcripts,
            PREPARE_REFERENCE_FILES.out.annotations_host_salmon,
        )
        ch_versions = ch_versions.mix(SALMON_ALIGNMENT_BASED.out.versions)
        salmon_ab_out = SALMON_ALIGNMENT_BASED.out
    }

    // Run if STAR genome alignment
    if (params.run_star) {
        STAR_ALIGNMENT(
            ch_reads,
            PREPARE_REFERENCE_FILES.out.host_pathogen_fasta_genome,
            PREPARE_REFERENCE_FILES.out.host_pathogen_genes_gff,
        )
    }

    //Capture software versions
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
