# nf-core/dualrnaseq: Output

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [FastQC](#fastqc)
* [Trimming reads](#trimming-reads)
* [Mapping/Quantification](#mappingquantification)
* [MultiQC](#multiqc)
* [Pipeline info](#pipeline-info)
* [Mapping statistics](#mapping-statistics)

## FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads.
It provides information about the quality score distribution across your reads and the per base sequence content (%T/A/G/C).
Information about adapter contamination and other overrepresented sequences is also displayed.

**Output directory for raw reads:** `results/fastqc`

**Output directory for trimmed reads:** `results/fastqc_after_trimming`

Contents:

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the appropriate directory based on the trimming tool used.

## Trimming reads

Two software options are provided to remove low quality reads and adapters: [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

All samples are processed in the output directory and appended with `_trimmed.fastq.gz`

**Output directory:** `results/trimming`

## Mapping/quantification

Depending on which mapping/quantification mode was selected, will depend on which folders and files appear here.

### A) Salmon selective alignment

**Output directory:** `results/salmon`

**Description:** All files and folders produced by Salmon. Also contains separated quantification results from the host and pathogen,
plus gene and transcript quantification results.

### B) STAR + Salmon - alignment based

**Output directory:** `results/salmon_alignment_mode`

**Description:** All files and folders produced by Salmon. Also contains separated quantification results from the host and pathogen,
plus gene and transcript quantification results.

**Output directory:** `results/STAR_for_salmon`

**Description:** All files and folders produced by STAR that are required for the alignment-based mode of Salmon. This folder also contains the STAR index.

### C) STAR + HTSeq

**Output directory:** `results/STAR`

**Description:** All files and folders produced by STAR when aligning to a genome. This folder also contains the STAR index.

**Output directory:** `results/HTSeq`

**Description:** Contains counted features and TPMs for the host and pathogen; as a combined file, and as separate count files.

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project.

Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory:** `results/MultiQC`

* `ProjectName_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `ProjectName_multiqc_report_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline
* `multiqc_plots/`
  * Directory containing plots from within the .html report, saved as `PDF`, `PNG` and `SVG`

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## Pipeline info

Directory containing all of the pipeline-specific reports, timeline and descriptions.

Detailed descriptions of each file can be found on the [Nextflow website](https://www.nextflow.io/docs/latest/tracing.html)

**Output directory:** `results/pipeline_info`

* `execution_report.html`
* `execution_timeline.html`
* `execution_trace.txt`
* `pipeline_dag.svg`
* `pipeline_report.html`
* `pipeline_report.txt`
* `results_description.html`
* `software_versions.csv`

## Mapping statistics

Depending on which mapping/quantification mode was selected, will depend on which folders and files appear here.
In general, files and images within these folders show the number of reads, mapping statistics for both the host and pathogen,
scatter plots showing correlations between replicates within the same conditions, RNA-class statistics, as well as gene and transcript-based metrics.

**Output directory:** `results/mapping_statistics`

Contents:

**HTSeq:** `results/HTSeq`

* `uniquely_mapped`
* `RNA_classes_host`
* `RNA_classes_pathogen`
* `scatter_plots`

**Salmon:** `results/salmon`

* `RNA_classes_host`
* `RNA_classes_pathogen`
* `scatter_plots`

 **Salmon alignment based:** `results/salmon_alignment_based`

* `RNA_classes_host`
* `RNA_classes_pathogen`
* `scatter_plots`

 **STAR:** `results/STAR`

* `multi_mapped`
* `processed_reads`
* `uniquely_mapped`

**STAR for Salmon:** `results/STAR_for_salmon`

* `processed_reads`
