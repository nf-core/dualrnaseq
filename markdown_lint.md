
#### `nf-core lint` overall result: Passed :white_check_mark:

Posted for pipeline commit 93ea44e

```diff
+| ✅ 107 tests passed       |+
!| ❗ 18 tests had warnings |!
-| ❌  0 tests failed       |-
```

<details>

### :heavy_exclamation_mark: Test warnings:

* [Test #4](https://nf-co.re/errors#4) - Config `process.container` looks wrong. Should be `nfcore/dualrnaseq:dev` but is `dualrnaseq/dualrnaseq:dev`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS full test should test full datasets: `./.github/workflows/awsfulltest.yml`
* [Test #6](https://nf-co.re/errors#6) - Found a bioconda environment.yml file but no badge in the README
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `bioconda::samtools=1.9`, `1.11` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `bioconda::STAR=2.7.3a`, `2.7.6a` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `bioconda::multiqc=1.8`, `1.9` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `conda-forge::python=3.7.6`, `3.9.0` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `conda-forge::markdown=3.1.1`, `3.3.3` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `conda-forge::pymdown-extensions=6.0`, `8.0.1` available
* [Test #8](https://nf-co.re/errors#8) - Conda package is not latest available: `pip=20.0.2`, `20.2.4` available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.76, 1.78 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.18.2, 1.19.3 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 1.0.2, 1.1.3 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 0.15.4, 0.16.0.1 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 2.8, 2.10 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 3.1.1, 3.3.2 available
* [Test #8](https://nf-co.re/errors#8) - PyPi package is not latest available: 0.10.0, 0.11.0 available
* [Test #10](https://nf-co.re/errors#10) - TODO string found in `original_ci_yml_file_docker`: _You can customise CI pipeline run tests as required_

### :white_check_mark: Tests passed:

* [Test #1](https://nf-co.re/errors#1) - File found: `nextflow.config`
* [Test #1](https://nf-co.re/errors#1) - File found: `nextflow_schema.json`
* [Test #1](https://nf-co.re/errors#1) - File found: `LICENSE` or `LICENSE.md` or `LICENCE` or `LICENCE.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `README.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `CHANGELOG.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/README.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/output.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `docs/usage.md`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/branch.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/ci.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/linting.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `main.nf`
* [Test #1](https://nf-co.re/errors#1) - File found: `environment.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `Dockerfile`
* [Test #1](https://nf-co.re/errors#1) - File found: `conf/base.config`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/awstest.yml`
* [Test #1](https://nf-co.re/errors#1) - File found: `.github/workflows/awsfulltest.yml`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `Singularity`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `parameters.settings.json`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `bin/markdown_to_html.r`
* [Test #1](https://nf-co.re/errors#1) - File not found check: `.travis.yml`
* [Test #3](https://nf-co.re/errors#3) - Licence check passed
* [Test #2](https://nf-co.re/errors#2) - Dockerfile check passed
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.name`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.nextflowVersion`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.description`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.version`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.homePage`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `timeline.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `trace.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `report.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `dag.enabled`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.cpus`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.memory`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.time`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `params.outdir`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `params.input`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `manifest.mainScript`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `timeline.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `trace.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `report.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `dag.file`
* [Test #4](https://nf-co.re/errors#4) - Config variable found: `process.container`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.version`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.nf_required_version`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.container`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.singleEnd`
* [Test #4](https://nf-co.re/errors#4) - Config variable (correctly) not found: `params.igenomesIgnore`
* [Test #4](https://nf-co.re/errors#4) - Config `timeline.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `report.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `trace.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `dag.enabled` had correct value: `true`
* [Test #4](https://nf-co.re/errors#4) - Config `manifest.name` began with `nf-core/`
* [Test #4](https://nf-co.re/errors#4) - Config variable `manifest.homePage` began with https://github.com/nf-core/
* [Test #4](https://nf-co.re/errors#4) - Config `dag.file` ended with `.svg`
* [Test #4](https://nf-co.re/errors#4) - Config variable `manifest.nextflowVersion` started with >= or !>=
* [Test #4](https://nf-co.re/errors#4) - Config `manifest.version` ends in `dev`: `'1.0dev'`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions 'branch' workflow is triggered for PRs to master: `./.github/workflows/branch.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions 'branch' workflow looks good: `./.github/workflows/branch.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions CI is triggered on expected events: `./.github/workflows/ci.yml`
* [Test #5](https://nf-co.re/errors#5) - CI is building the correct docker image: `docker build --no-cache . -t dualrnaseq/dualrnaseq:dev`
* [Test #5](https://nf-co.re/errors#5) - CI is pulling the correct docker image: docker pull dualrnaseq/dualrnaseq:dev
* [Test #5](https://nf-co.re/errors#5) - CI is tagging docker image correctly: docker tag dualrnaseq/dualrnaseq:dev dualrnaseq/dualrnaseq:dev
* [Test #5](https://nf-co.re/errors#5) - Continuous integration checks minimum NF version: `./.github/workflows/ci.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions linting workflow is triggered on PR and push: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - Continuous integration runs Markdown lint Tests: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - Continuous integration runs nf-core lint Tests: `./.github/workflows/linting.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS test is triggered on workflow_dispatch: `./.github/workflows/awstest.yml`
* [Test #5](https://nf-co.re/errors#5) - GitHub Actions AWS full test is triggered only on published release and workflow_dispatch: `./.github/workflows/awsfulltest.yml`
* [Test #6](https://nf-co.re/errors#6) - README Nextflow minimum version badge matched config. Badge: `19.10.0`, Config: `19.10.0`
* [Test #8](https://nf-co.re/errors#8) - Conda environment name was correct (nf-core-dualrnaseq-1.0dev)
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::fastqc=0.11.9`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::fastqc=0.11.9`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::samtools=1.9`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::salmon=1.3.0`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::salmon=1.3.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::STAR=2.7.3a`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::gffread=0.12.1`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::gffread=0.12.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::multiqc=1.8`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::htseq=0.12.4`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::htseq=0.12.4`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::bioconductor-rtracklayer=1.48.0`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::bioconductor-rtracklayer=1.48.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `bioconda::bioconductor-tximport=1.16.0`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `bioconda::bioconductor-tximport=1.16.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `conda-forge::python=3.7.6`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `conda-forge::r-plyr=1.8.6`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `conda-forge::r-plyr=1.8.6`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `conda-forge::markdown=3.1.1`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `conda-forge::pymdown-extensions=6.0`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `agbiome::bbtools=37.62`
* [Test #8](https://nf-co.re/errors#8) - Conda package is latest available: `agbiome::bbtools=37.62`
* [Test #8](https://nf-co.re/errors#8) - Conda dependency had pinned version number: `pip=20.0.2`
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: biopython==1.76
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: numpy==1.18.2
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: pandas==1.0.2
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: pysam==0.15.4
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: cutadapt==2.8
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: matplotlib==3.1.1
* [Test #8](https://nf-co.re/errors#8) - Pip dependency had pinned version number: seaborn==0.10.0
* [Test #9](https://nf-co.re/errors#9) - Found all expected strings in Dockerfile file
* [Test #12](https://nf-co.re/errors#12) - Name adheres to nf-core convention
* [Test #13](https://nf-co.re/errors#13) - Did not find any cookiecutter template strings (128 files)
* [Test #14](https://nf-co.re/errors#14) - Schema lint passed
* [Test #14](https://nf-co.re/errors#14) - Schema title + description lint passed
* [Test #15](https://nf-co.re/errors#15) - Schema matched params returned from nextflow config

### Run details:

* nf-core/tools version 1.11
* Run at `2020-10-29 20:59:36`

</details>
