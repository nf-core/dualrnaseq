## `nf-core pipelines lint` overall result: Failed :x:

Posted for pipeline commit 8ffabab

```diff
+| ✅ 174 tests passed       |+
!| ❗  25 tests had warnings |!
-| ❌  44 tests failed       |-
```

<details>

### :x: Test failures:

* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File must be removed: `lib/NfcoreTemplate.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File must be removed: `lib/Utils.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File must be removed: `lib/WorkflowMain.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File must be removed: `lib/WorkflowDualrnaseq.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File must be removed: `lib/nfcore_external_java_deps.jar`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_cpus`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_memory`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (incorrectly) found: `params.max_time`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value incorrect: `params.run_salmon_SA` is set as `true` in `nextflow_schema.json` but is `false` in `nextflow.config`.
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value incorrect: `params.run_salmon_AB` is set as `false` in `nextflow_schema.json` but is `true` in `nextflow.config`.
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value incorrect: `params.tracedir` is set as `${params.outdir}/pipeline_info` in `nextflow_schema.json` but is `null/pipeline_info` in `nextflow.config`.
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `CODE_OF_CONDUCT.md` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/CONTRIBUTING.md` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/ISSUE_TEMPLATE/bug_report.yml` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/PULL_REQUEST_TEMPLATE.md` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/workflows/branch.yml` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/workflows/linting_comment.yml` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/workflows/linting.yml` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `assets/email_template.html` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `assets/email_template.txt` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `assets/nf-core-dualrnaseq_logo_light.png` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `docs/images/nf-core-dualrnaseq_logo_light.png` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `docs/images/nf-core-dualrnaseq_logo_dark.png` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.gitignore` does not match the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.prettierignore` does not match the template
* [actions_awsfulltest](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_awsfulltest) - `.github/workflows/awsfulltest.yml` is not triggered correctly
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `host_gff_attribute` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `gene_feature_gff_to_quantify_host` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `extract_annotations_host_salmon_feature` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `extract_annotations_host_salmon_attribute` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `pathogen_gff_attribute` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `gene_feature_gff_to_quantify_pathogen` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `htseq_quantifier` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `outSAMunmapped` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `sjdbGTFfeatureExon` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `sjdbGTFtagExonParentTranscript` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `quantMode` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `quantTranscriptomeBan` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `limitBAMsortRAM` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `schema_ignore_params` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Param `fasta` from `nextflow config` not found in nextflow_schema.json
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Default value for param `run_salmon_SA` invalid: Schema default (`True`) does not match the config default (`false`)
* [schema_params](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_params) - Default value for param `run_salmon_AB` invalid: Schema default (`False`) does not match the config default (`true`)
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` does not contain a matching 'report_comment'.  
The expected comment is:  
```This report has been generated by the <a href="https://github.com/nf-core/dualrnaseq/tree/dev" target="_blank">nf-core/dualrnaseq</a> analysis pipeline. For information about how to interpret these results, please see the <a href="https://nf-co.re/dualrnaseq/dev/docs/output" target="_blank">documentation</a>.```  
The current comment is:  
```This report has been generated by the <a href="https://github.com/nf-core/dualrnaseq" target="_blank">nf-core/dualrnaseq</a> analysis pipeline. For information about how to interpret these results, please see the <a href="https://nf-co.re/dualrnaseq" target="_blank">documentation</a>.```

### :heavy_exclamation_mark: Test warnings:

* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found: `conf/igenomes_ignored.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found: `ro-crate-metadata.json`
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `LICENSE` does not match the template
* [readme](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/readme) - README contains the placeholder `zenodo.XXXXXXX`. This should be replaced with the zenodo doi (after the first release).
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Write a 1-2 sentence summary of what data the pipeline is for and what it does_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Add full-sized test dataset and amend the paragraph below if applicable_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Fill in short bullet-pointed list of the default steps in the pipeline_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Update the example "typical command" below used to run the pipeline_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _If applicable, make list of people who have also contributed_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file._
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `README.md`: _Add bibliography of tools and data used in your pipeline_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `methods_description_template.yml`: _#Update the HTML below to your prefered methods description, e.g. add publication citation for this pipeline_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `base.config`: _Check the defaults for all processes_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `base.config`: _Customise requirements for specific processes._
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `ci.yml`: _You can customise CI pipeline run tests as required_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `awsfulltest.yml`: _You can customise AWS full pipeline tests as required_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `WorkflowMain.groovy`: _Add Zenodo DOI for pipeline after first release_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `usage.md`: _Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website._
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `output.md`: _Write this documentation describing your workflow's output_
* [pipeline_todos](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_todos) - TODO string in `dualrnaseq.nf`: _Add all file path parameters for the pipeline to the list below_
* [system_exit](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/system_exit) - `System.exit` in WorkflowDualrnaseq.groovy: _//     System.exit(1)_  [line 18]
* [system_exit](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/system_exit) - `System.exit` in WorkflowDualrnaseq.groovy: _System.exit(1)_  [line 74]
* [system_exit](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/system_exit) - `System.exit` in WorkflowMain.groovy: _System.exit(1)_  [line 85]
* [system_exit](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/system_exit) - `System.exit` in NfcoreSchema.groovy: _System.exit(1)_  [line 180]
* [nfcore_yml](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nfcore_yml) - nf-core version not set in `.nf-core.yml`

### :white_check_mark: Tests passed:

* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.gitattributes`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.gitignore`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.nf-core.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.editorconfig`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.prettierignore`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.prettierrc.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `CHANGELOG.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `CITATIONS.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `CODE_OF_CONDUCT.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `LICENSE` or `LICENSE.md` or `LICENCE` or `LICENCE.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `nextflow_schema.json`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `nextflow.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `README.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/.dockstore.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/CONTRIBUTING.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/bug_report.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/config.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/ISSUE_TEMPLATE/feature_request.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/PULL_REQUEST_TEMPLATE.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/branch.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/ci.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting_comment.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/linting.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.html`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `assets/email_template.txt`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `assets/sendmail_template.txt`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `assets/nf-core-dualrnaseq_logo_light.png`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `conf/modules.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `conf/test.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `conf/test_full.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/images/nf-core-dualrnaseq_logo_light.png`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/images/nf-core-dualrnaseq_logo_dark.png`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/output.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/README.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `docs/usage.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `main.nf`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `assets/multiqc_config.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `conf/base.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `conf/igenomes.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/awstest.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `.github/workflows/awsfulltest.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File found: `modules.json`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/bug_report.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.github/ISSUE_TEMPLATE/feature_request.md`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.github/workflows/push_dockerhub.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.markdownlint.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.nf-core.yaml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.yamllint.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `bin/markdown_to_html.r`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `conf/aws.config`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `docs/images/nf-core-dualrnaseq_logo.png`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Checks.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Completion.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `lib/Workflow.groovy`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `parameters.settings.json`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `pipeline_template.yml`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `Singularity`
* [files_exist](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_exist) - File not found check: `.travis.yml`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.name`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.nextflowVersion`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.description`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.version`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.homePage`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.enabled`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.cpus`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.memory`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `process.time`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.outdir`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `params.input`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `manifest.mainScript`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `timeline.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `trace.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `report.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable found: `dag.file`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.nf_required_version`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.container`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.singleEnd`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.igenomesIgnore`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.name`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable (correctly) not found: `params.enable_conda`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``timeline.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``report.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``trace.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``dag.enabled`` had correct value: ``true``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``manifest.name`` began with ``nf-core/``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable ``manifest.homePage`` began with https://github.com/nf-core/
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``dag.file`` ended with ``.html``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config variable ``manifest.nextflowVersion`` started with >= or !>=
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config ``manifest.version`` ends in ``dev``: ``2.0dev``
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config `params.custom_config_version` is set to `master`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config `params.custom_config_base` is set to `https://raw.githubusercontent.com/nf-core/configs/master`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Lines for loading custom profiles found
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - nextflow.config contains configuration profile `test`
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.igenomes_base= s3://ngi-igenomes/igenomes
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.gene_feature_gff_to_create_transcriptome_host= ['exon', 'tRNA']
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.gene_feature_gff_to_create_transcriptome_pathogen= ['gene', 'sRNA', 'tRNA', 'rRNA']
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.gene_attribute_gff_to_create_transcriptome_host= transcript_id
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.gene_attribute_gff_to_create_transcriptome_pathogen= locus_tag
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.salmon_sa_index_args= -k 21
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.salmon_sa_args= --softclipOverhangs
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.run_star= false
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.run_htseq= false
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_version= master
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.custom_config_base= https://raw.githubusercontent.com/nf-core/configs/master
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_cpus= 16
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_memory= 128.GB
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_time= 240.h
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.publish_dir_mode= copy
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.max_multiqc_email_size= 25.MB
* [nextflow_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nextflow_config) - Config default value correct: params.validate_params= true
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.gitattributes` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.prettierrc.yml` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/.dockstore.yml` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/ISSUE_TEMPLATE/config.yml` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `.github/ISSUE_TEMPLATE/feature_request.yml` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `assets/sendmail_template.txt` matches the template
* [files_unchanged](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/files_unchanged) - `docs/README.md` matches the template
* [actions_ci](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_ci) - '.github/workflows/ci.yml' is triggered on expected events
* [actions_ci](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_ci) - '.github/workflows/ci.yml' checks minimum NF version
* [actions_awstest](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_awstest) - '.github/workflows/awstest.yml' is triggered correctly
* [actions_awsfulltest](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_awsfulltest) - `.github/workflows/awsfulltest.yml` does not use `-profile test`
* [readme](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/readme) - README Nextflow minimum version badge matched config. Badge: `22.10.1`, Config: `22.10.1`
* [plugin_includes](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/plugin_includes) - No wrong validation plugin imports have been found
* [pipeline_name_conventions](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/pipeline_name_conventions) - Name adheres to nf-core convention
* [template_strings](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/template_strings) - Did not find any Jinja template strings (0 files)
* [schema_lint](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_lint) - Schema lint passed
* [schema_lint](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_lint) - Schema title + description lint passed
* [schema_lint](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/schema_lint) - Input mimetype lint passed: 'text/csv'
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: ci.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: awstest.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting_comment.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: linting.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: fix-linting.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: awsfulltest.yml
* [actions_schema_validation](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/actions_schema_validation) - Workflow validation passed: branch.yml
* [merge_markers](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/merge_markers) - No merge markers found in pipeline files
* [modules_json](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_json) - Only installed modules found in `modules.json`
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` found and not ignored.
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_section_order`
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `export_plots`
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains `report_comment`
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` follows the ordering scheme of the minimally required plugins.
* [multiqc_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/multiqc_config) - `assets/multiqc_config.yml` contains 'export_plots: true'.
* [modules_structure](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_structure) - modules directory structure is correct 'modules/nf-core/TOOL/SUBTOOL'
* [base_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/base_config) - `conf/base.config` found and not ignored.
* [base_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/base_config) - `CUSTOM_DUMPSOFTWAREVERSIONS` found in `conf/base.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `conf/modules.config` found and not ignored.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `SAMPLESHEET_CHECK` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `COMBINE_FILES` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `REPLACE_ATTRIBUTE_GFF_STAR_SALMON_HOST` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `UNCOMPRESS_HOST_FASTA_GENOME` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `REPLACE_GENE_FEATURE_GFF_HOST_SALMON` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `NFCORE_DUALRNASEQ` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `FASTQC` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `FASTQC_AFTER_TRIMMING` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `CUTADAPT` found in `conf/modules.config` and Nextflow scripts.
* [modules_config](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/modules_config) - `CUSTOM_DUMPSOFTWAREVERSIONS` found in `conf/modules.config` and Nextflow scripts.
* [nfcore_yml](https://nf-co.re/tools/docs/3.1.2/pipeline_lint_tests/nfcore_yml) - Repository type in `.nf-core.yml` is valid: `pipeline`

### Run details

* nf-core/tools version 3.1.2
* Run at `2025-01-23 14:43:24`

</details>
