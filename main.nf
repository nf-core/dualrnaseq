#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/dualrnaseq
========================================================================================
    Github : https://github.com/nf-core/dualrnaseq
    Website: https://nf-co.re/dualrnaseq
    Slack  : https://nfcore.slack.com/channels/dualrnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { DUALRNASEQ } from './workflows/dualrnaseq'

//
// WORKFLOW: Run main nf-core/dualrnaseq analysis pipeline
//
workflow NFCORE_DUALRNASEQ {
    DUALRNASEQ ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_DUALRNASEQ ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
