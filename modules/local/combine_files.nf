process COMBINE_FILES {
    tag "combine_files"
    label 'process_high'

    input:
        path(file1)
        path(file2)
        val(output_file)
    output:
        path output_file

    script:
    def args = task.ext.args ?: ''
    """
    cat ${file1} ${file2} > ${output_file}
    """
}