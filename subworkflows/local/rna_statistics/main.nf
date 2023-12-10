// include {
//     RNA_CLASS_STATISTICS_UNIQUELY_MAPPED as RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST;
//     RNA_CLASS_STATISTICS_UNIQUELY_MAPPED as RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_PATOGEN;
// } from '../../../modules/local/RNA_class_statistics_uniquely_mapped/main.nf'

// include {
//     PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH as PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_HOST;
//     PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH as PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_PATOGEN;
// } from '../../../modules/local/plot_RNA_class_uniquely_mapped_each/main.nf'


// include {
//     PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED as PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_HOST;
//     PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED as PLOT_RNA_CLASS_UNIQUELY_MAPPED_COMBINED_PATHOGEN;
// } from '../../../modules/local/plot_RNA_class_uniquely_mapped_combined/main.nf'


// include {
//     SCATTER_PLOT as SCATTER_PLOT_HOST;
//     SCATTER_PLOT as SCATTER_PLOT_PATHOGEN;
// } from '../../modules/local/scatter_plot/main.nf'

/*
    Common statistics used for Star HTSEQ and SALOMON SELECTIVE ALIGNMENT
*/
workflow RNA_STATISTICS {
    // take:
    //     host_quantification
    //     pathogen_quantification
    //     attribute
    //     tool // salmon, htseq

    // main:
    //     SCATTER_PLOT_HOST(
    //         .host_quantification,
    //         attribute,
    //         Channel.value("replicates"), // output from CHECK_REPLICATES
    //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.scatterplots_host_htseq,
    //         "host"
    //     )
    //     SCATTER_PLOT_PATHOGEN(
    //         pathogen_quantification,
    //         attribute,
    //         Channel.value("replicates"), // output from CHECK_REPLICATES
    //         SPLIT_QUANTIFICATION_TABLES_HTSEQ_UNIQUELY_MAPPED.out.scatterplots_pathogen_htseq,
    //         "pathogen"
    //     )

        // // statistics 
        // RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_HOST(
        //     host_quantification,
        //     attribute,
        //     annotations_host_htseq
        //     "host"
        //     tool
        // )
        // RNA_CLASS_STATISTICS_UNIQUELY_MAPPED_PATOGEN(
        //     pathogen_quantification,
        //     attribute,
        //     annotations_pathogen_htseq
        //     "pathogen"
        //     tool
        // )

        // // each
        // PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_HOST(
        //     RNA_CLASS_STATISTICS_HTSEQ_UNIQUELY_MAPPED_HOST.out
        // )
        // PLOT_RNA_CLASS_UNIQUELY_MAPPED_EACH_PATOGEN(
        //     RNA_CLASS_STATISTICS_HTSEQ_UNIQUELY_MAPPED_PATHOGEN.out
        // )

        // // combined
        // PLOT_RNA_CLASS_UNIQUELY_COMBINED_HOST(
        //     RNA_CLASS_STATISTICS_HTSEQ_UNIQUELY_MAPPED_HOST.OUT
        // )
        // PLOT_RNA_CLASS_UNIQUELY_COMBINED_PATHOGEN(
        //     RNA_CLASS_STATISTICS_HTSEQ_UNIQUELY_MAPPED_PATHOGEN.out
        // )
}