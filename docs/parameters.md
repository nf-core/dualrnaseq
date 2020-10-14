# nf-core/dualrnaseq: Parameters

## Table of contents

1. [Nextflow](#1-nextflow)
2. [Genome references and annotation](#2-genome-references-and-annotation)
3. [Input sequence reads](#3-input-sequence-reads)
4. [FastQC and adapter trimming](#4-fastqc-and-adapter-trimming)
5. [Salmon](#5-salmon)
6. [STAR](#6-star)
7. [HTSeq](#7-htseq)
8. [RNA mapping statistics](#8-rna-mapping-statistics)
9. [Other](#9-other)

### General comment

> All of the parameters listed here can be found in either the main configuration file `nextflow.config` or `base.config`. Alternatively, each parameter can be specified by the user when they require adjustments to the default settings. The format for parameters is either a flag telling the pipeline to run something, such as `--run_STAR`, or to specify a particular value `--max_cpus 16`, string `--outWigStrand "Stranded"` or file `--outdir "/path_to_file/file"`.
Although many of the parameters listed below are set as `False` in the configuration files - their usage on the command line will generally not require setting them to either True or False. Instead, by passing a parameter it becomes set to True.

## 1. Nextflow

### Core parameters

> Note: most of the core Nextflow parameters only require a single hyphen

#### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
This is used in the MultiQC report and in the summary HTML / e-mail.

#### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results where the inputs are the same, continuing from where it ran to previously.
Alternatively, you can supply a run name to resume a specific run: `-resume [run-name]`.
  
If unsure of the run name, you can use the `nextflow log` command to show previous run names.

#### `-c`

Specify the path to a specific config file.

> Note: you can use this to override pipeline defaults.

### Pipeline resources

#### `--max_memory 128.GB`

Use to set the maximum memory for each process.

#### `--max_time 240.h`

Use to set the max run time for each process.

#### `--max_cpus 16`

Use to set the max number of CPUs for each process.

### Reslts directory

#### `--outdir`

Where results will be saved (should be enclosed by quotation marks `"..."`).

### Custom configuration

#### `--custom_config_version`

Provide git commit id for custom institutional configs hosted at `nf-core/configs`.

This was implemented for reproducibility purposes. Default: `master`

```bash
## Download and use config file with following git commit id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

#### `--custom_config_base`

If you're running offline, nextflow will be unable to fetch the institutional config files from the internet. If required, files should be downloaded when an active internet connection is available, then pointing Nextflow to those files with this option.

For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note: to make this process easier, the nf-core/tools helper package has a `download` command to download all required pipeline files + singularity containers + institutional configs in one go.

## 2. Genome references and annotation

### Genomes

#### `--fasta_host`

```bash
--fasta_host "<path to host genome fasta file>"
```

#### `--fasta_pathogen`

### References

#### `--gff_host`

#### `--gff_host_tRNA`

#### `--gff_pathogen`

### Transcriptome

The first two parameters are set to `False`. If supplying custom transcriptome files, add the appropriate flags below.

#### `--read_transcriptome_fasta_host_from_file`

#### `--read_transcriptome_fasta_pathogen_from_file`

#### `--transcriptome_host`

#### `--transcriptome_pathogen`

> The above nine parameters are all set as `False`. If specified, the folder/file should be enclosed by quotations `"..."`.

## 3. Input sequence reads

### Details

Input files can be read as either .fastq or .fastq.gz. They should be named descriptively without spaces and special characters (such as : and @), with the corresponding replicate (if any) appended at the end. The best practise for this pipeline is to use underscores to separate different experimental conditions.

#### `--reads`

Please note the following requirements:

* The path must be enclosed in quotes
* The path must have at least one `*` wildcard character
* When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.
* If left unspecified, a default pattern is used: `"data/*{1,2}.fastq.gz"`
* It is not possible to run a mixture of single-end and paired-end files in one run

> Note: by default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--single_end` on the command line when launched. For example: `--single_end --reads '*.fastq'`

## 4. FastQC and Adapter trimming

### FastQC

By default, the pipeline utilizes FastQC tool for quality control of raw sequencing reads. To learn more on FastQC, please check [`FastQC website.`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

#### `--skipFastqc`

An option to not run FastQC. (Default: False)

This is set to False within the configuration files, but only needs to be passed on the command line to become True.

#### `--fastqc_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional fastqc parameters you wish to use, except --quiet --threads --noextract flags which are already specified in the dualrnaseq pipeline.

### Adapter trimming

Trimming is performed by either Cutadapt or BBDuk with the following related options:  

### Cutadapt

Cutadapt requires prior knowledge of the adaptors used during library preparation.

By default, the pipeline trims Illumina TruSeq adapters. See [`Illumina TruSeq.`](https://cutadapt.readthedocs.io/en/stable/guide.html#illumina-truseq)

To learn more on Cutadapt and its parameters visit the [`cutadapt documentation.`](https://cutadapt.readthedocs.io/en/stable/guide.html)

#### `--run_cutadapt`

To run Cutadapt (Default: False)

#### `--a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"`

For single-end reads as well as the first reads of paired-end data, adapter sequence can be specified with `--a` flag. For more information, see [`adapter-types.`](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types)

#### `--A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"`

For paired-end data, the adapter sequence for the second reads can be defined here. For more information, see [`trimming paired-end reads.`](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads).

#### `--quality_cutoff 10` or `--quality-cutoff 10,15`

Cutadapt can also remove low-quality read ends. By default, the 3’ end of each read is trimmed using a cutoff of 10. For more information on cutoff values, see [`quality trimming.`](https://cutadapt.readthedocs.io/en/stable/guide.html#quality-trimming)

If you specify two comma-separated cutoffs, the first value represents the 5’ cutoff, and the second one the 3’ cutoff.

#### `--cutadapt_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional cutadapt parameters you wish to use, except -m and -j which are already specified in the dualrnaseq pipeline.

### BBDuk

BBDuk does not require any prior knowledge about adapter types, searching for common adapter sequences from the file `$baseDir/data/adapters.fa`.

To learn more about BBDuk and its parameters visit the [`BBDuk website.`](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/)

#### `--run_bbduk`

To run BBDuk (Default: False)

#### `--minlen 18`

Reads shorter than this after trimming will be discarded (Pairs will be discarded if both are shorter).

#### `--qtrim "r"`

To trim read ends to remove bases with quality below trimq.

Possible options:`rl` (trim both ends), `f` (neither end), `r` (right end only) `l` (left end only), `w` (sliding window).

#### `--trimq 10`

Cutoff to trim regions with average quality BELOW given value.

Option is avaiable if qtrim is set to something other than f. Reads shorter than this after trimming will be discarded (Pairs will be discarded if both are shorter).

#### `--ktrim "r"`

To trim reads to remove bases matching reference kmers. Avaiable options: `f` (don't trim), `r` (trim to the right - 3' adapters), `l` (trim to the left - 5' adapters).

#### `--k 17`

Kmer length used for finding contaminants (adapters). Contaminants shorter than k will not be found. k must be at least 1.

#### `--mink 11`

Look for shorter kmers at read tips down to this length, when k-trimming or masking. 0 means disabled. Enabling this will disable maskmiddle

#### `--hdist 1`

Maximum Hamming distance for ref kmers (subs only).

#### `--adapters`

Fasta file with adapter sequences (Default: `$baseDir/data/adapters.fa`).

#### `--BBDuk_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional BBDuk parameters you wish to use, except -Xmx1g which is already specified in the dualrnaseq pipeline.

## 5. Salmon

### General parameters

These parameters are available for Salmon in both Selective Alignment and alignment-based mode.  

#### `--libtype SF`

To define the sequencing library of your data.

To learn more on library types available in Salmon, please read [_`What’s this LIBTYPE?`_](https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype)

#### `--incompatPrior 0.0`

By default, this is set to `0.0`, to ensure that only mappings or alignments that are compatible with the specified library type are considered by Salmon. You can find more information on this parameter in the [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#incompatprior)

#### `--generate_salmon_uniq_ambig`

Option to extract all of the unique and ambigious reads after quantification. Works for both Selective alignment and alignment-based modes (Default: False).

#### `--gene_feature_gff_to_create_transcriptome_host "[exon, tRNA]"`

The pipeline uses gene features from the 3rd column of the host annotative (gff) file, to extract the coordinates of transcripts to be quantified.

By default, the pipeline uses `exon` from the `--gff_host` file and `tRNA` from the `--gff_host_tRNA` file.

#### `--gene_feature_gff_to_create_transcriptome_pathogen "[gene, sRNA, tRNA, rRNA]"`

The pipeline uses gene features from the 3rd column of the pathogen annotative (gff) file, to extract the coordinates of transcripts to be quantified.

By default, the pipeline uses features as `gene`, `sRNA`, `tRNA` and `rRNA` from the `--gff_pathogen` file.

#### `--gene_attribute_gff_to_create_transcriptome_host "transcript_id"`

This flag defines the gene attribute from the 9th column of the host annotative (gff) file, where the transcript names are extracted.

By default, the pipeline extracts `transcript_id` from the `--gff_host` file.

#### `--gene_attribute_gff_to_create_transcriptome_pathogen "locus_tag"`

This flag defines the gene attribute from the 9th column of the pathogen annotative (gff) file, where transcripts, genes or CDS regions are extracted.

By default, the pipeline extracts `locus_tag` from the `--gff_pathogen` file.

### Salmon Selective Alignment

Parameters listed below are available only for Salmon with Selective Alignment.

#### `--run_salmon_selective_alignment`

To run Salmon with selective alignment (Default: False).

#### `--kmer_length 21`

To define the k-mer length (`-k` parameter in Salmon, see [`preparing transcriptome indices`](https://salmon.readthedocs.io/en/latest/salmon.html?highlight=index#preparing-transcriptome-indices-mapping-based-mode)). By default, this parameter is set to 21.

#### `--writeUnmappedNames`

By default the pipeline does not save names of unmapped reads. You can learn more about this option in [`Salmon documentation`](https://salmon.readthedocs.io/en/latest/salmon.html#writeunmappednames). If you don't want to keep this option, set the `--writeUnmappedNames` flag to False.
(Default: False).

#### `--softclipOverhangs`

By default, the pipeline does not allow soft-clipping of reads (Default: False).

_"Soft-clipping allows reads that overhang the beginning or ends of the transcript. In this case, the overhanging section of the read will simply be unaligned, and will not contribute or detract from the alignment score"_.
If it is set to `False`, the end-to-end alignment of the entire read is forced, so that the occurance of any overhangings may affect the alignment score.

#### `--dumpEq`

To save the equivalence classes and their counts, change this option to `True`. See [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#dumpeq) for more information
(Default: False).

#### `--writeMappings`

If set to `True`, the pipeline will create a `mapping.sam` file containing mapping information. To learn more on this option, please view the [`Salmon documentation.`](https://salmon.readthedocs.io/en/latest/salmon.html#writemappings)
(Default: False).

#### `--keepDuplicates`

Option to remove/collapse identical transcripts during the indexing stage (Default: False).


#### `--salmon_SA_params_index "--param_a 4 --param_b 5 -param_x"`

Define set of additional salmon index parameters you wish to use in selective alignment mode.

#### `--salmon_SA_params_mapping "--param_a 4 --param_b 5 -param_x"`

Define set of additional salmon quant parameters you wish to use.

### Salmon alignment based mode

#### `--run_salmon_alignment_based_mode`

Option to run Salmon in alignment-based mode (Default: False).

#### `--salmon_alignment_based_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional salmon quant parameters you wish to use in salmon alignment-based mode.

## 6. STAR

### General parameters

These parameters are available for STAR in both quantification modes, using HTSeq and Salmon in alignment-based mode.

#### `--run_star`

Option to run STAR (Default: False).

#### `--outSAMunmapped "Within"`

By default, the pipeline saves unmapped reads within the main BAM file. If you want to switch off this option, set the `--outSAMunmapped` flag to `None`.  See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more details.

For paired-end reads, the `KeepPairs` parameter will record the unmapped mates for each alignment, and will keep it adjacent to its mapped read (only affects multi-mapping reads).

#### `--outSAMattributes "Standard"`

To specify the attributes of the output BAM file. The default value is `Standard`, but there are a range of options if needed. Please see the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for the full list.

By default, the pipeline uses `Standard` option to keep NH HI AS nM SAM attributes.

#### `--outFilterMultimapNmax 999`

To specify the maximum number of loci a read is allowed to map to.
By default, this  option is set to 999 in the pipeline. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterType "BySJout"`

By default, the pipeline keeps reads containing junctions that passed filtering into the file `SJ.out.tab`. This option reduces the number of ”spurious” junctions. (ENCODE standard options for long RNA-seq pipeline). You can read more about the flag and its options in the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf)

#### `--alignSJoverhangMin 8`

The number of minimum overhang for unannotated junctions can be changed here. By default, the pipeline uses 8. (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignSJDBoverhangMin 1`

The number of minimum overhang for annotated junctions can be changed here. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterMismatchNmax 999`

To define a threshold for the number of mismatches to be allowed. By default, the pipeline uses a large number `999` to switch this filter off. (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--outFilterMismatchNoverReadLmax 1`

Here, you can define a threshold for a ratio of mismatches to *read* length. The alignment will be considered if the ratio is less than or equal to this value. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignIntronMin 20`

By default, the nf-core dualrnaseq pipeline uses `20` as a minimum intron length. If the genomic gap is smaller than this value, it is considered as a deletion.
(ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignIntronMax 1000000`

The maximum intron length is set to 1,000,000 (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--alignMatesGapMax 1000000`

The maximum genomic distance between mates is 1,000,000 (ENCODE standard options for long RNA-seq pipeline). See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--limitBAMsortRAM 0`

Option to limit RAM when sorting BAM file. If `0`, will be set to the genome index size, which can be quite large when running on a desktop or laptop.

#### `--winAnchorMultimapNmax 999`

The maximum number of loci anchors that are allowed to map. By default, the pipeline uses a large number `999` to switch this filter off.

#### `--sjdbOverhang 100`

Option to specify the length of the donor/acceptor sequence on each side of the junctions used in constructing the splice junctions database. By default the option is set to `100`. However, we recommend to set a value depending on the read length: read/mate length - 1.

### STAR for HTSeq

Parameters available for STAR - HTSeq quantification method.

#### `--outWigType "None"`

Used to generate signal outputs, such as "wiggle" and "bedGraph". To view all available signal types, please see the [`STAR documentation`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf).

By default, the pipeline does not generate any of these files.

#### `--outWigStrand "Stranded"`

Options are `Stranded` or `Unstranded` when defining the strandedness of wiggle/bedGraph output. See [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--STAR_index_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional star parameters to create index.

#### `--STAR_alignment_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional star alignment parameters.

### STAR for Salmon alignment-based mode

#### `--quantTranscriptomeBan "Singleend"`

The nf-core/dualrnaseq pipeline runs STAR to generate a transcriptomic alignments. By default, it allows for insertions, deletions and soft-clips (`Singleend` option). To prohibit this behaviour, please specify `IndelSoftclipSingleend`. See the [`STAR documentation.`](https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf) for more information.

#### `--STAR_salmon_alignment_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional alignment parameters for STAR in salmon alignment-based mode.

## 7. HTSeq

### Parameters

#### `--run_htseq_uniquely_mapped`

Used to run HTSeq-count and extract uniquely mapped reads from both the host and pathogen (Default: False).

#### `--stranded "yes"`

A parameter for the library type. Options include `"yes"` or `"no"`.

#### `--max_reads_in_buffer 30000000`

Option to define the number of maximum reads allowed to stay in memory until the mates are found. Has an effect for paired-end reads (Default: 30000000).

#### `--minaqual 10`

To specify a threshold for a minimal MAPQ alignment quality.
By default, this parameter is set to 10.

### Gene features and attributes

The four parameters below are used to extract gene features from both the host and pathogen. These values may need to be changed, especially for the pathogen, as many different names exist, such as `ID`, `Gene`, `Name`, `locus_tag` etc.

A good idea is to view the accompanying annotative file and examine the fields within.

> Note: If a `tRNA.gff` file is included, it is assumed that it has the same gene atribute as the annotative (gff) file, i.e. `gene_id`

#### Host

#### `--gene_feature_gff_to_quantify_host "[exon, tRNA]"`

#### `--host_gff_atribute "gene_id"`

#### Pathogen

#### `--gene_feature_gff_to_quantify_pathogen "[gene, sRNA, tRNA, rRNA]"`

#### `--pathogen_gff_atribute "locus_tag"`

#### `--htseq_params "--param_a 4 --param_b 5 -param_x"`

Define set of additional htseq parameters you wish to use in the pipeline.

## 8. RNA mapping statistics

### Parameters and files

#### `--mapping_statistics`

Option to generate mapping statistics (Default: False).

This will create the following:

* Count the total number of reads before and after trimming
* Scatterplots comparing all replicates (separate for both host and pathogen reads)
* Plots of the % of mapped/quantified reads
* Plots of RNA-class statistics (as many types can be identified, the parameter below `--RNA_classes_to_replace_host` can help to summarise these)

#### `--RNA_classes_to_replace_host "$baseDir/data/RNA_classes_to_replace.csv"`

Located within the `data/` folder of dualrnaseq, this tab delimited file contains headers which groups similar types of RNA classes together. This helps to keep the RNA-class names simplified for plotting purposes.

Initially, the user can run the pipeline without the 'others' class (remove the 'others' column) to identify the concentration of all RNA types,including e.g. scRNAs). Depending on the requirements, the user can decide which types should be included/excluded or grouped together.  

## 9. Other

### Reports

#### `--email`

Set this parameter with your e-mail address to get a summary e-mail with details of the run when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

#### `--email_on_fail`

This works the same as `--email`, except emails are only sent if the workflow is not successful.

#### `--max_multiqc_email_size 25.MB`

Threshold size for MultiQC report to be attached in the notification email. If file generated by pipeline exceeds the threshold, it will not be attached.

#### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

#### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

#### `--multiqc_config`

Specify path to a custom MultiQC configuration file.

### AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

#### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

#### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

#### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.
