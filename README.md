# ![nf-core/dualrnaseq](docs/images/nf-core-dualrnaseq_logo.png)

**Dual RNA-seq pipeline**.

[![GitHub Actions CI Status](https://github.com/nf-core/dualrnaseq/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/dualrnaseq/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/dualrnaseq/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/dualrnaseq/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/309982089.svg)](https://zenodo.org/badge/latestdoi/309982089)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)
[![Docker](https://img.shields.io/docker/automated/nfcore/dualrnaseq.svg)](https://hub.docker.com/r/nfcore/dualrnaseq)
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23dualrnaseq-4A154B?logo=slack)](https://nfcore.slack.com/channels/dualrnaseq)

## Introduction

**nf-core/dualrnaseq** is specifically used for the analysis of Dual RNA-seq data, interrogating host-pathogen interactions through simultaneous RNA-seq.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

This pipeline has been initially tested with eukaryotic host's including Human and Mouse, and pathogens including *Salmonella enterica*, *Orientia tsutsugamushi*, *Streptococcus penumoniae*, *Escherichia coli* and *Mycobacterium leprae*. The workflow should work with any eukaryotic and bacterial organisms with an available reference genome and annotation.

### Method

The workflow merges host and pathogen genome annotations taking into account differences in annotation conventions, then processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)),   quantifies gene expression ([STAR](https://github.com/alexdobin/STAR) and [HTSeq](https://htseq.readthedocs.io/en/master/); [STAR](https://github.com/alexdobin/STAR), [Salmon](https://combine-lab.github.io/salmon/) and [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html); or [Salmon](https://combine-lab.github.io/salmon/) in quasimapping mode and [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)), and summarises the results ([MultiQC](http://multiqc.info/)), as well as generating a number of custom summary plots and separate results tables for the pathogen and host. See the [output documentation](docs/output.md) for more details.

### Workflow

The workflow diagram below gives a simplified visual overview of how dualrnaseq has been designed.

![nf-core/dualrnaseq](docs/images/Workflow_diagram_dualrnaseq.png)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(please only use [`Conda`](https://conda.io/miniconda.html) as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_

3. Download the pipeline and test it on a minimal dataset with a single command:

    ```bash
    nextflow run nf-core/dualrnaseq -profile test,<docker/singularity/podman/shifter/charliecloud/conda/institute>
    ```

    > Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.

4. Start running your own analysis!

    <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

    ```bash
    nextflow run nf-core/dualrnaseq -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input '*_R{1,2}.fastq.gz'
    ```

See [usage docs](https://nf-co.re/dualrnaseq/usage) for all of the available options when running the pipeline.

## Documentation

The nf-core/dualrnaseq pipeline comes with documentation about the pipeline: [usage](https://nf-co.re/dualrnaseq/usage) and [output](https://nf-co.re/dualrnaseq/output).

## Credits

nf-core/dualrnaseq was originally written by Bozena Mika-Gospodorz and Regan Hayward.

We thank the following people for their extensive assistance in the development
of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#dualrnaseq` channel](https://nfcore.slack.com/channels/dualrnaseq) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
