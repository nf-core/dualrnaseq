process CHECK_REPLICATES {
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