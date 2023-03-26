# nf-core/dualrnaseq: Running the pipeline

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/dualrnaseq/usage](https://nf-co.re/dualrnaseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Table of contents

1. [Running the pipeline](#1-running-the-pipeline)
   * [Quick start](#11-quick-start)
   * [Basic run](#12-basic-run)
   * [Updating the pipeline](#13-updating-the-pipeline)
   * [Reproducibility](#14-reproducibility)
2. [Configuration profile](#2-configuration-profile)
3. [Input sequence reads](#3-input-sequence-reads)
4. [Reference genomes and annotation](#4-reference-genomes-and-annotation)
   * [Genomes](#41-genomes)
   * [Annotation](#42-annotation)
5. [Trimming reads and adapters](#5-trimming-reads-and-adapters)
   * [Cutadapt](#51-cutadapt)
   * [BBDuk](#52-bbduk)
6. [Read mapping and quantification](#6-read-mapping-and-quantification)
   * [Salmon - Selective alignment](#61-salmon---selective-alignment)
   * [Salmon - quantification in alignment-based mode](#62-salmon---quantification-in-alignment-based-mode)
   * [STAR - alignment-based genome mapping + quantification with HTSeq](#63-star---alignment-based-genome-mapping-+-feature-counting-with-HTSeq)
7. [Mapping statistics](#7-mapping-statistics)
8. [Example usage](#8-example-usage)
9. [Output files](#9-output-files)
10. [Job resources](#10-job-resources-and-submission)

## 1. Running the pipeline

### 1.1 Quick Start

i. Install [`nextflow`](https://nf-co.re/usage/installation)

ii. Install either [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility (please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))

iii. Download the pipeline and test it on a minimal dataset with a single command

```bash
nextflow run nf-core/dualrnaseq -profile test,<docker/singularity/conda/institute>
```

> Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

iv. Start running your own analysis!

```bash
nextflow run nf-core/dualrnaseq -profile <docker/singularity/conda/institute> --input '*_R{1,2}.fastq.gz' --genome_host GRCh38 --genome_pathogen SL1344
```

To see all of the available parameters when running the pipeline, please see [the parameter documentation](https://nf-co.re/dualrnaseq/parameters).

You can also use the [interactive launch tool](https://nf-co.re/launch?pipeline=dualrnaseq), which has embedded help and descriptions for each parameter.
There is a web-based interface or a purely command-line interface. Launch via the web link above or on the command line:

```bash
nf-core launch dualrnaseq
```

### 1.2 Basic run

Once ready, a basic command for running the pipeline would be the following:

```bash
nextflow run nf-core/dualrnaseq
    -profile docker \
    --input "/folder_to_reads/*_R{1,2}.fq.gz" \
    --fasta_host host.fa \
    --fasta_pathogen pathogen.fa \
    --gff_host host.gff \
    --gff_pathogen pathogen.gff \
    --run_star --outdir results
```

This will launch the pipeline with the `docker` configuration profile. Click [here](#https://nf-co.re/usage/configuration), or see [below](#2-configuration-profile) for more information about profiles.
It takes all compressed .fq files in the specified folder and will map the reads (using STAR) to the supplied host and pathogen genomes.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and logs.
```

### 1.3 Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. Subsequent runs will use the cached version if available - even if the pipeline has since been updated. To make sure that you're running the latest version, make sure that you regularly run this command so the cached version is the most recent:

```bash
nextflow pull nf-core/dualrnaseq
```

### 1.4 Reproducibility

It's a good idea to specify a pipeline version when running on your data. This will ensure results can be reproduced more easily, and can be easily referenced or referred to, particularly when multiple versions are available. If a version number was't strictly defined in a previous run, it can be found in the various reports in the `results` directory.

To find the latest version number, go to the [nf-core/dualrnaseq releases page](https://github.com/nf-core/dualrnaseq/releases) - which will be a numeric value (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `nextflow run -r 1.3.1 nf-core/dualrnaseq/main.nf`.

## 2. Configuration profile

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/dualrnaseq`](https://hub.docker.com/r/nfcore/dualrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/dualrnaseq`](https://hub.docker.com/r/nfcore/dualrnaseq/)
* `podman`
  * A generic configuration profile to be used with [Podman](https://podman.io/)
  * Pulls software from Docker Hub: [`nfcore/dualrnaseq`](https://hub.docker.com/r/nfcore/dualrnaseq/)
* `shifter`
  * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
  * Pulls software from Docker Hub: [`nfcore/dualrnaseq`](https://hub.docker.com/r/nfcore/dualrnaseq/)
* `charliecloud`
  * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
  * Pulls software from Docker Hub: [`nfcore/dualrnaseq`](https://hub.docker.com/r/nfcore/dualrnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

#### Custom resource requests

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

Whilst these default requirements will hopefully work for most people with most data, you may find that you want to customise the compute resources that the pipeline requests. You can do this by creating a custom config file. For example, to give the workflow process `star` 32GB of memory, you could use the following config:

```nextflow
process {
  withName: star {
    memory = 32.GB
  }
}
```

To find the exact name of a process you wish to modify the compute resources, check the live-status of a nextflow run displayed on your terminal or check the nextflow error for a line like so: `Error executing process > 'bwa'`. In this case the name to specify in the custom config file is `bwa`.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

### Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

#### Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## 3. Input sequence reads

`--input`

Input files can be read as either uncompressed or compressed (gzip) fasta or fastq files. They should be named descriptively without spaces and special characters (such as : and @), with the corresponding replicate (if any) denoted with a capital `R` or lower case `r`, and read number `1` or `2` appended at the end. The best practise for this pipeline is to use underscores to separate different experimental conditions, for example:

| Paired-end | Single-end |
| --- | --- |
| Host_pathogen_mock_10h_R1_1.fq | Host_pathogen_mock_10h_r1.fq |
| Host_pathogen_mock_10h_R1_2.fq | Host_pathogen_mock_10h_r2.fq |
| Host_pathogen_mock_10h_R2_1.fq | Host_pathogen_mock_10h_r3.fq |
| Host_pathogen_mock_10h_R2_2.fq | Host_pathogen_MOI_1_10h_R1.fq |
| Host_pathogen_mock_10h_R3_1.fq | Host_pathogen_MOI_1_10h_R2.fq |
| Host_pathogen_mock_10h_R3_2.fq | Host_pathogen_MOI_1_10h_R3.fq |
| Host_pathogen_MOI_1_10h_r1_1.fq |  |
| Host_pathogen_MOI_1_10h_r1_2.fq |  |
| Host_pathogen_MOI_1_10h_r2_1.fq |  |
| Host_pathogen_MOI_1_10h_r2_2.fq |  |
| Host_pathogen_MOI_1_10h_r3_1.fq |  |
| Host_pathogen_MOI_1_10h_r3_2.fq |  |

Once correctly named, instead of typing all files, the folder containing these files can be specified, such as `--input "folder_to_files/*.fq.gz"`. Or, for paired-end `--input "folder_to_files/*{1,2}.fastq.gz"`

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### Additional parameters

The following two parameters are associated with the sequencing library type. Their defaults are shown below, but should be changed to the experiment-specific values.

`--single_end = false`

`--stranded = "yes"`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when launched.

> Note: it is not possible to run a mixture of single-end and paired-end files in one run.

## 4. Reference genomes and annotation

### 4.1 Genomes

The main goal of Dual RNA-seq is simultaneous profiling of host and pathogen gene expression. Thus, the pipeline requires references for each of the organisms.

These parameters can be used in two ways:

#### A) Using configuration files

Most nf-core pipelines, including this one, come with a pre-built set of reference genome indices that work out of the box.
They are hosted in the cloud (see <https://ewels.github.io/AWS-iGenomes/>) and will be downloaded on demand.
For this pipeline however, it is very likely that you will need to specify additional custom genomes to work with.

> See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and [Reference genomes](https://nf-co.re/usage/reference_genomes) for instructions.

In short, the best way to do this is to extend the genomes config that comes with the pipeline, using the same configuration structure.
You can then save this file in a number of different locations, according to your preference.
See <https://nf-co.re/usage/configuration#custom-configuration-files> for details on this, but in short the common places are:

* In the Nextflow home directory as a file called `config` (no file extensions). Typically this is in your user's home directory: `~/.nextflow/config`
* As a file called `nextflow.config` in the current working directory when you run the pipeline.
* Any custom path, which you will then need to specify on the command line with `-c`. eg: `-c path/to/config` (multiple files can be given like this)

Depending on which location you use, you may want to disable the bundled igenomes to avoid conflicts.
You can do this by setting `igenomes_ignore = true` in the config, or `--igenomes_ignore` on the command line.

The syntax for a custom reference file follows this syntax:

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta_host  = '<path to the genome fasta file>'
      gff_host = '<path to the genome annotation file>'
      gff_host_tRNA = '<path to the tRNA annotation file>' // Optional
      transcriptome_host = '<path to the transcriptome fasta file>' // Optional
    }
    'SL1344' {
      fasta_pathogen  = '<path to the genome fasta file>'
      gff_pathogen = '<path to the genome gff annotation file>'
      transcriptome_pathogen = '<path to the transcriptome fasta file>' // Optional
    }
  }
  // Default genomes (optional).
  // Ignored if --genome_host 'OTHER-GENOME' and --genome_pathogen 'OTHER-GENOME' specified on command line
  genome_host = 'GRCh38'
  genome_pathogen = 'SL1344'
}
```

Defining default genomes in your configuration file is optional. You can specify the references using the following flags on command line:

* `--genome_pathogen SL1344`
* `--genome_host GRCh38`

Any number of additional genome references can be added to this file and specified using the given genome key, through either `--genome_host` or `--genome_pathogen`.

Please note that:

* If `--transcriptome_host` or `--transcriptome_pathogen` is not given, the transcriptome fasta will be created by the pipeline using the provided genome and annotation files.
* If `gff_host_tRNA` file is provided, the pipeline combines the files from `gff_host` and `gff_host_tRNA` to create a single host gff file.

#### B) Using pipeline-specific parameters

If you preferr, you can specify each parameter on the command line when you run the pipeline.
Reference and annotation files (`fasta` and `GFF3`) can be compressed (`.gz` or `.zip`) or uncompressed.

Host:

* `--fasta_host /path/to/file`
* `--gff_host /path/to/file`
* `--gff_host_tRNA /path/to/file`
* `--transcriptome_host /path/to/file`

Pathogen:

* `--fasta_pathogen /path/to/file`
* `--gff_pathogen /path/to/file`
* `--transcriptome_pathogen /path/to/file`

These parameters can also be set in a nextflow config, or supplied in a JSON / YAML file with `-params-file`.

For example, with a file `my-config.yml`:

```yaml
fasta_host: /path/to/file
gff_host: /path/to/file
gff_host_tRNA: /path/to/file
transcriptome_host: /path/to/file
```

You can run the pipeline with:

```bash
nextflow run nf-core/dualrnaseq -params-file my-config.yml
```

##### Host tRNA

We have specified this parameter for users familiar with the [Gencode gene annotations](https://www.gencodegenes.org/). Their annotative files include lncRNAs, snoRNAs, rRNAs and other non-coding genes except tRNAs. tRNAs are available in another gff file (predicted tRNA genes). The tRNA gff file looks a little different than the main annotation file, so we don't recommend adding different GFF files in its place.

##### Warning! The nf-core/dualrnaseq pipeline does not support iGenomes

Many nf-core pipelines provide an option to use iGenomes [Reference genomes](https://nf-co.re/usage/reference_genomes). However, the nf-core/dualrnaseq pipeline does not support this functionality. Thie pipeline requires GFF files and not GTF files.

### 4.2 Annotation

Host-based annotations are generally more uniform in design than pathogen annotations, where different terms for the same feature are often used. For example, in Human annotation files, main features are defined as genes, transcripts and exons, and associated with gene_id, transcript_id and other uniform identifiers. When extracting features, this uniform naming convention makes this straight forward. However, bacterial naming conventions are less uniform. Names for features include genes, CDS, ID, Name, locus_tag amongst others.

We have defined four parameters for both the host and pathogen to account for these differences in conventions - and should be modified when needed.
Default values are shown below:

Host:

`--gene_attribute_gff_to_create_transcriptome_host "transcript_id"`

`--gene_feature_gff_to_create_transcriptome_host "[exon, tRNA]"`

`--gene_feature_gff_to_quantify_host "[exon, tRNA]"`

`--host_gff_attribute "gene_id"`

Pathogen:

`--gene_attribute_gff_to_create_transcriptome_pathogen "locus_tag"`

`--gene_feature_gff_to_create_transcriptome_pathogen "[gene, sRNA, tRNA, rRNA]"`

`--gene_feature_gff_to_quantify_pathogen "[gene, sRNA, tRNA, rRNA]"`

`--pathogen_gff_attribute "locus_tag"`

#### Features within pathogen annotation files

As previously mentioned, bacterial annotative files can be challenging to work with due to non-uniformity. We find that summarising the contents of the GFF can help to identify which features you may want to examine. The following bash script can do this easily:

```bash
awk -F '\t' '{print $3}' file.gff3 | sort | uniq -c
```

## 5. Trimming reads and adapters

By default, trimming and adapter removal is not run. To run either tool, you will need to specify either `run_cutadapt` or `run_bbduk` during run time.

For convenience, two software options are provided to remove low quality reads and adapters:

### 5.1 Cutadapt

Cutadapt is best suited when the library preparation steps and adapter types are known. By default, parameters are set up to remove TruSeq adaptors and with a quality cutoff of 10 (from the 3' end).
For the full list of parameters, click [here](https://nf-co.re/dualrnaseq/parameters#cutadapt).

### 5.2 BBDuk

The advantage of using BBDuk is that no prior knowledge of library preparation steps are needed. There is a file stored within the pipeline (`$baseDir/assets/adapters.fa`) containing common adapter types, from which BBDuk will search through and remove. This is extremely useful when analysing data from public repositories when minimal background has been given. For the full list of parameters, click [here](https://nf-co.re/dualrnaseq/parameters#bbduk).

## 6. Read mapping and quantification

The nf-core/dualrnaseq pipeline provides three strategies to map and quantify your dual RNA-seq data.

1) The first approach **Salmon with Selective alignment**, performs mapping using the Selective alignment algorithm which quantifies the reads. `--run_salmon_selective_alignment`

2) The second approach **Salmon with alignment-based mode**, utilizes aligned reads from STAR to quantify input reads. `--run_salmon_alignment_based_mode`

3) The third strategy utilises alignment-based mapping executed with **STAR**, counting uniquely aligned reads with **HTSeq**. `--run_star` and `--run_htseq_uniquely_mapped`

### 6.1 Salmon - selective alignment

Salmon is a transcriptome-based mapping tool that performs both mapping and quantification. In the first phase it performs indexing of reference transcripts (pathogen transcripts are defined as gene or CDS features), where a chimeric transcriptome of host and pathogen files is created. During this step, coordinates of gene features are also extracted from the host and pathogen annotation files.
To avoid spurious mapping of reads that originate from unannotated locus to sequences similar to annotated transcripts, a decoy-aware transcriptome is created and incorporated into the index. In the pipeline, decoy sequences are created from only the host genome (please see our paper for more information *paper coming soon*).

> Note: **Selective alignment** is an improvement to the original alignment-free approach (also called quasi-mapping). In Selective-alignment, the best transcript for a read from a set of mappings is selected based on the alignment-based score instead of the longest exact match - which increases accuracy. See [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/salmon.html) for more information.

To summarize transcript-level estimates obtained with Salmon into gene-level abundance estimates, the nf-core/dualrnaseq pipeline uses [`Tximport.`](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)

### 6.2 Salmon - quantification in alignment-based mode

In this [mode](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), Salmon performs quantification utilising an aligned BAM file. In this pipeline, the alignment file is generated with STAR. The first step involves creating an index of a chimeric genome (created from the host and pathogen genome files). Next, STAR performs an alignment, but for the purpose of Salmon (it generates alignments translated into transcript coordinates). To learn more on this behavior, please see `Output in transcript coordinates` from the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

> Note: there are numerous STAR-based flags that can be modified within the pipeline - which can be viewed [here](https://nf-co.re/dualrnaseq/parameters#star---general).
> Salmon performs quantification based on a reference transcriptome. It is recommended to allow the pipeline to create a transcriptome using the provided genome (fasta) and annotative (gff) files.
> When quantifying alignments, the parameters `--libtype` and `--incompatPrior` should be adjusted as required.

Gene-level estimates are obtained using [`Tximport.`](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html)

### 6.3 STAR - alignment-based genome mapping + feature counting with HTSeq

STAR is a splice-aware alignment tool which aligns reads to a reference genome. In this pipeline, STAR generates a chimeric genome index (combining host and pathogen genomes), then identifies and maps spliced alignments across splice junctions.
If using this method, the following three parameters (or links in a genome configuration file) are required: `--genome_host`, `--genome_pathogen` and `--gff_host`.

Users have the option to only run STAR if desired. If also running HTSeq-count, you will need to pass `--run_htseq_uniquely_mapped`.

To quantify uniquely mapped host and pathogen reads, the pipeline uses HTSeq-count. In addition to the host gff, other parameters must be specified including `--gff_pathogen`, `--gene_feature_gff_to_quantify_host`, `--host_gff_atribute`, `--gene_feature_gff_to_quantify_pathogen` and `--pathogen_gff_atribute`.

## 7. Mapping statistics

To summarise the mapping statistics including total mapped reads, unmapped reads, host-specific and pathogen-specific mapped reads, add the following parameter to the command line: `--mapping_statistics`.

This will create the following:

* Count the total number of reads before and after trimming
* Scatterplots comparing all replicates (separate for both host and pathogen genes)
* Plots of the % of mapped/quantified reads
* Plots of RNA-class statistics

To see examples of this output, [click here](output.md#mapping-statistics).

## 8. Example usage

As discussed above in the [Read mapping and quantification section](#62-salmon---quantification-in-alignment-based-mode) above, there are different ways to run the pipeline:

### Example 1

* Using Docker

* Single end reads

* Host (Human) and pathogen (*E. coli*) genomes defined through `genomes.conf`

* Salmon - Selective alignment

```bash
nextflow run dualrnaseq/main.nf -profile docker,cluster \
--genome_host "GRCh38" --genome_pathogen "Escherichia_coli_K_12_DH10B" \
--input "folder_to_reads/*.fq.gz" --single_end --genomes_ignore = true \
--outdir "/outdir_folder/" \
--run_salmon_selective_alignment \
```

### Example 2

* Using Docker

* Paired-end reads (unstranded)

* Host (Mouse) and pathogen (*C. trachomatis*) genomes defined through `genomes.conf`

* Salmon - quantification in alignment-based mode

* Custom kmer length

 ```bash
qsub -q all.q nextflow run dualrnaseq/main.nf -profile docker \
--genome_host "GRCm38" --genome_pathogen "C_trachomatis_strain_d" \
--input "folder_to_reads/*{1,2}.fastq.gz" --genomes_ignore = true \
--outdir "/outdir_folder/" \
--run_salmon_alignment_based_mode --libtype "IU" --incompatPrior 0.0 --kmer_length 19 \
```

### Example 3

* Using Singularity

* Single end reads

* Host (Human) and pathogen (*M. pneumoniae*) genomes defined through `genomes.conf`

* STAR - alignment-based genome mapping

* HTSeq - quantification of uniquely mapped reads

 ```bash
nextflow run dualrnaseq/main.nf -profile singularity \
--genome_host "GRCH38" --genome_pathogen "Mycoplasma_pneumoniae" \
--input "folder_to_reads/*.fq.gz" --single_end --genomes_ignore = true \
--outdir /outdir_folder/ \
--run_star --run_htseq_uniquely_mapped \
```

### Example 4

* Run on a SGE cluster

* Using Singularity

* Single end reads

* Host (Human) and pathogen (*S. typhimurium*) genomes defined through `genomes.conf`

* All three modes with custom gene attributes

```bash
qsub -q all.q nextflow run dualrnaseq/main.nf -profile singularity,cluster \
--genome_host "GRCH38" --genome_pathogen "Salmonella_typhimurium" \
--input "folder_to_reads/*.fq.gz" --single_end --genomes_ignore = true \
--outdir "/outdir_folder/" \
--run_salmon_alignment_based_mode --libtype "SF" \
--run_salmon_selective_alignment \
--run_star --run_htseq_uniquely_mapped \
--gene_feature_gff_to_create_transcriptome_pathogen "[ID, sRNA, tRNA, rRNA]" \
--gene_feature_gff_to_quantify_pathogen "[ID, sRNA, tRNA, rRNA]"
```

## 9. Output files

Click [here](output.md) for a description on output files.

## 10. Job resources and submission

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus, Nextflow processes must run until the pipeline is finished. To achieve this, we recommend running in the background through `screen` / `tmux` or a similar tool. Alternatively, you can run dualrnaseq submitted by your job scheduler on a cluster.

It is also recommended to limit the Nextflow Java virtual machine memory. This can be achieved by adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

### Automatic resubmission

Each step in the pipeline has a default set of requirements for the number of CPUs, memory and time. For most steps within the pipeline, if the job exits with an error code of `143` (exceeded requested resources), it will automatically resubmit with higher requirements (2 x original, then 3 x original). If these requests continue to fail, then the pipeline will stop.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly, it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition above). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

### AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

`--awsqueue`  The JobQueue that you intend to use on AWS Batch.

`--awsregion` The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

`--awscli`  The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't
