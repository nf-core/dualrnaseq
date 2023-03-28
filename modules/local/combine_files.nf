process COMBINE_FILES {
    tag "combine_files"
    label 'process_high'

    input:
        path("*.fasta") 

    output:
        path "output.fasta"

    script:
    """
        cat *.fasta > output.fasta
    """
}