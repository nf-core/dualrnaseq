process SAMPLESHEET_CHECK {
    tag "${samplesheet}"
    label 'process_single'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    conda "python=3.8.3"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'nfcore/dualrnaseq:dev'
        : 'nfcore/dualrnaseq:dev'}"

    input:
    path samplesheet

    output:
    path '*.valid.csv', emit: csv
    path "versions.yml", emit: versions

    script:
    """
    # Ensure script exists and is executable
    if [ ! -f "${workflow.projectDir}/bin/check_samplesheet.py" ]; then
        echo "Error: check_samplesheet.py not found"
        exit 1
    fi

    # Copy and run script
    cp "${workflow.projectDir}/bin/check_samplesheet.py" .
    chmod +x check_samplesheet.py

    ./check_samplesheet.py \\
        ${samplesheet} \\
        samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
