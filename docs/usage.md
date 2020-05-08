# nf-core/dualrnaseq: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--reads`](#--reads)
  * [`--single_end`](#--single_end)
* [Reference genomes](#reference-genomes)
  * [`--genome_host` (using iGenomes)](#--genome_host-using-igenomes)
  * [`--genome_pathogen` (using iGenomes)](#--genome_pathogen-using-igenomes)
  * [`--fasta_host`](#--fasta_host)
  * [`--fasta_pathogen`](#--fasta_pathogen)
  * [`--gff_host`](#--gff_host)
  * [`--gff_host_tRNA`](#--gff_host_tRNA)
  * [`--gff_pathogen`](#--gff_pathogen)
  * [`--transcriptome_host`](#--transcriptome_host)
  * [`--transcriptome_pathogen`](#--transcriptome_pathogen)
  * [`--igenomes_ignore`](#--igenomes_ignore)
* [Adapter trimming](#Adapter-trimming)
    * [`--a`](#--a)
    * [`--A`](#--A)
    * [`--quality-cutoff`](#--quality-cutoff)
    * [`--skipTrimming`](#--skipTrimming)
* [Read mapping and quantification](#Read-mapping-and-quantification)
  * [STAR - alignment-based genome mapping](#STAR---alignment-based-genome-mapping)
    * [`--a`](#--a)
  * [HTSeq - quantification of uniquely-mapped reads](#HTSeq---quantification-of-uniquely-mapped-reads)
  * [HTSeq - quantification of multi-mapped reads](#HTSeq---quantification-of-multi-mapped-reads)
  * [Salmon - Selective-Alignment and quantification](#Salmon---Selective-Alignment-and-quantification)
  * [Salmon - quantification in alignment-based mode](#Salmon---quantification-in-alignment-based-mode)
* [Maping statistics](#Maping-statistics)

* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)


## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

<!-- TODO nf-core: Document required command line parameters to run the pipeline-->

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/dualrnaseq --reads '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```
More information about running the pipeline can be found in the `docs/` directory:

- [Output and how to interpret the results](output.md) - change !!!!!!!!!!!!!!!!!!!!!!!!!!
- [Extra Documentation on annotation](annotation.md) - change !!!!!!!!!!!!!!!!!!!!

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/dualrnaseq
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/dualrnaseq releases page](https://github.com/nf-core/dualrnaseq/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/dualrnaseq`](http://hub.docker.com/r/nfcore/dualrnaseq/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/dualrnaseq`](http://hub.docker.com/r/nfcore/dualrnaseq/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

<!-- TODO nf-core: Document required command line parameters -->

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--single_end`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--single_end --reads '*.fastq'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

## Reference genomes
The main goal of Dual RNA-seq is simultaneous profiling of host and pathogen gene expression. Thus, the pipeline requires to provide references for each of the organisms ([genome_host](#--genome_host-(using-iGenomes)) and [genome_pathogen](#--genome_pathogen-(using-iGenomes))).

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.
  
### `--genome_host` (using iGenomes)
There are 3 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome_host` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common host genomes that are supported are:

* Human
  * `--genome_host GRCh38`
* Mouse
  * `--genome_host GRCm38`

<!--  > There are numerous others - check the config file for more.-->

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

<!-- TODO nf-core: Update reference genome example according to what is needed -->

```nextflow
params {
  genomes {
    'GRCh38' {
      fasta_host  = '<path to the genome fasta file>'
      gff_host = '<path to the genome gff annotation file>'
      gff_host_tRNA = '<path to the tRNA gff annotation fasta file>'
      transcriptome_host = '<path to the tRNA gff annotation fasta file>'
    }
  }
// Any number of additional genomes, key is used with --genome
}

```
### `--genome_pathogen` (using iGenomes)
To run the pipeline with pathogen references defined in iGenomes, you must specify which iGenomes keys to use with the `--genome_pathogen` flag.

You can fing the iGenomes keys in [iGenomes config file](../conf/igenomes.config). Common pahogen genome that is supported is:

* _Salmonella_ Typhimurium
  * `--genome_pathogen SL1344`

If your genome of interest is not provided with iGenomes you can create your own configuration file and save it as conf/genomes.config.

The syntax for this reference configuration is as follows:

<!-- TODO nf-core: Update reference genome example according to what is needed -->

```nextflow
params {
  genomes {
    'SL1344' {
      fasta_pathogen  = '<path to the genome fasta file>'
      gff_pathogen = '<path to the genome gff annotation file>'
      transcriptome_pathogen = '<path to the genome gff annotation file>'
    }
  }
// Any number of additional genomes, key is used with --genome
}

```
<!-- TODO nf-core: Describe reference path flags -->

### `--fasta_host`

If you prefer, you can specify the full path to your host genome fasta file when you run the pipeline:

```bash
--fasta_host '[path to genome fasta reference of host]'
```

### `--fasta_pathogen`

If you prefer, you can specify the full path to your pathogen genome fasta file when you run the pipeline:

```bash
--fasta_pathogen '[path to genome fasta reference of pathogen]'
```

### `--gff_host`

If you prefer, you can specify the full path to genome annotation file of your host in the gff3 format:

```bash
--gff_host '[path to host genome gff file]'
```
### `--gff_host_tRNA`

If you wish you can specify tRNA host annotations. 
If you prefer, you can specify the full path to tRNA annotations of your host in the gff3 format::

```bash
--gff_host_tRNA '[path to host tRNA gff file]'
```

### `--gff_pathogen`

If you prefer, you can specify the full path to genome annotations of your pathogen in the gff3 format:

```bash
--gff_pathogen '[path to pathogen genome gff file]'
```

### `--read_transcriptome_fasta_host_from_file`

### `--read_transcriptome_fasta_pathogen_from_file`


### `--transcriptome_host`

If you prefer, you can specify the full path to your host transcriptome fasta file:
`--read_transcriptome_fasta_host_from_file`, def. false

```bash
--transcriptome_host '[path to host ranscriptome fasta file]'
```

### `--transcriptome_pathogen`

`--read_transcriptome_fasta_pathogen_from_file`, def. false
If you prefer, you can specify the full path to your pathogen transcriptome fasta file:

```bash
--transcriptome_pathogen '[path to transcriptome fasta file of pathogen]'
```

### `--igenomes_ignore`

Do not load `igenomes.config` when running the pipeline. You may choose this option if you observe clashes between custom parameters and those supplied in `igenomes.config`.


### Chimeric transcriptome

### Chimeric gff file

### `--gene_attribute_gff_to_create_transcriptome_host`
default- "transcript_id" 

### `--gene_feature_gff_to_create_transcriptome_host` 

def ["exon", "tRNA"] in chimeric gff file into quant

### `gene_attribute_gff_to_create_transcriptome_pathogen`
def "locus_tag" into parent

### `gene_feature_gff_to_create_transcriptome_pathogen` 

= ["gene","sRNA","tRNA","rRNA"] into quant

## Adapter trimming

To remove adapter sequences that were introduced during the library preparation the pipeline utilizes cutadapt.
To learn more on cutadapt and its parameters check the [`cutadapt documentation.`](https://cutadapt.readthedocs.io/en/stable/guide.html) 

By default, the pipeline trims Illumina TruSeq adapters. See [`Illumina TruSeq.`](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq) 

### `--a`
For single-end reads as well as the first reads of paired-end data, adapter sequence can be specified with `--a` flag. See [`adapter-types.`](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types) 

```bash
--a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
```
### `--A`
For paired-end data, the adapter sequence for the second reads can be defined by `--A`. See [`trimming paired-end reads.`](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads)

```bash
--A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```
### `--quality-cutoff`
Cutadapt can also remove low-quality ends of the reads. By default, the 3’ end of each read is trimmed using cutoff 10. You can change this threshold defining the `--quality-cutoff` flag. See [`quality trimming.`](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming)

```bash
--quality-cutoff 10
```
If you specify two comma-separated cutoffs, the first value represents 5’ cutoff, and the second one 3’ cutoff.
```bash
--quality-cutoff 10,15
```
### `--skipTrimming`
If you don't want to trim your reads, please specify the following flag:
```bash
--skipTrimming
```
or set the parameter to true in your config file. 

## Read mapping and quantification

### STAR - alignment-based genome mapping

### `--run_star`

To run STAR you must specify [`--fasta_host`](#--fasta_host) and [`--fasta_pathogen`](#--fasta_pathogen) which provide paths to host and pathogen genome fasta files. 


### HTSeq - quantification of uniquely-mapped reads

### HTSeq - quantification of multi-mapped reads

### Salmon - Selective-Alignment and quantification

### `--run_salmon_selective_alignment`

If host and pathogen transcriptomes are not provided with [`--transcriptome_host`](#--transcriptome_host) and [`--transcriptome_pathogen`](#--transcriptome_pathogen) flags, the pipeline can create transcriptomes based on genomic fasta files and gff annotation files. 

To create chimeric transriptome, it changes gene attributes to `--gene_attribute_gff_to_create_transcriptome_host` in gff tRNA and bacterial gff file. 
into quant

default- "transcript_id", changes into Parent attribute 

`--gene_feature_gff_to_create_transcriptome_host` = ["exon", "tRNA"] in chimeric gff file
into quant


`--gene_attribute_gff_to_create_transcriptome_pathogen` def "locus_tag" into parent
`--gene_feature_gff_to_create_transcriptome_pathogen`= ["gene","sRNA","tRNA","rRNA"] into quant


`--read_transcriptome_fasta_host_from_file` - `--transcriptome_host` 
`--gff_host_tRNA`:
gene_feature_gff_to_create_transcriptome_host
gene_attribute_gff_to_create_transcriptome_host
gff_host_genome
fasta_host

`--read_transcriptome_fasta_pathogen_from_file` -  `--transcriptome_pathogen`
gene_attribute_gff_to_create_transcriptome_pathogen#
gene_feature_gff_to_create_transcriptome_pathogen


### `--kmer_length`

def 31

### `--libtype`

### `-- writeUnmappedNames`

def true

### `--softclipOverhangs`
def = true

### `--incompatPrior_value`
def 0.0

### `--dumpEq`
def true

### `--writeMappings`

def false

### Salmon - quantification in alignment-based mode

### `--run_salmon_alignment_based_mode`

The pipeline utilizes STAR to create transcriptome alignment. .....
To run STAR you must specify [`--fasta_host`](#--fasta_host) and [`--fasta_pathogen`](#--fasta_pathogen) which provide paths to host and pathogen genome fasta files. 

If host and pathogen transcriptomes are not provided with [`--transcriptome_host`](#--transcriptome_host) and [`--transcriptome_pathogen`](#--transcriptome_pathogen) flags, the pipeline can create transcriptomes based on genomic fasta files and gff annotation files. 
For this, gff files and genomic fasta files are required. 

To create chimeric transriptome, it changes gene attributes to `--gene_attribute_gff_to_create_transcriptome_host` in gff tRNA and bacterial gff file. 
into quant

default- "transcript_id", changes into Parent attribute 

`--gene_feature_gff_to_create_transcriptome_host` = ["exon", "tRNA"] in chimeric gff file
into quant


`--gene_attribute_gff_to_create_transcriptome_pathogen` def "locus_tag" into parent

`--gene_feature_gff_to_create_transcriptome_pathogen`= ["gene","sRNA","tRNA","rRNA"] into quant


`--read_transcriptome_fasta_host_from_file` - `--transcriptome_host` 
`--gff_host_tRNA`:
gene_feature_gff_to_create_transcriptome_host
gene_attribute_gff_to_create_transcriptome_host
gff_host_genome
fasta_host

`--read_transcriptome_fasta_pathogen_from_file` -  `--transcriptome_pathogen`
gff_pathogen
fasta_pathogen
gene_attribute_gff_to_create_transcriptome_pathogen#
gene_feature_gff_to_create_transcriptome_pathogen

In the nf-core/dualrnaseq pipeline you can specify the following cutadapt parameters: 

## Maping statistics
### `--mapping_statistics`  
def. true

calc. count_total_reads , count_total_trimmed_reads, scatterplots

salmon - extract_processed_reads , unmapped

plot_salmon_mapping_stats_host_pathogen,
RNA class statistics - host, RNA_classes_to_replace.csv , def "data/RNA_classes_to_replace.csv" host
combined, each sample



## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

<!-- TODO nf-core: Describe any other command line flags here -->

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
