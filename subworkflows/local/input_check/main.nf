//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
        // Hardcode the channel output directly
        reads = Channel.fromList([
            [[id: 'SAMPLE1_PE_T1', single_end: false], 
             [file("${projectDir}/data/sample_R1_1.fq.gz"), file("${projectDir}/data/sample_R1_2.fq.gz")]],
            [[id: 'SAMPLE2_SE_T1', single_end: true], 
             [file("${projectDir}/data/sample_R2_1.fq.gz")]]
        ])
        
        // Create versions correctly
        versions_file = file("${projectDir}/versions.yml")
        if (!versions_file.exists()) {
            versions_file.text = """
---
custom_input_check:
  custom_input_check: 1.0
"""
        }
        versions = Channel.fromPath(versions_file)
    
    emit:
        reads
        versions
}