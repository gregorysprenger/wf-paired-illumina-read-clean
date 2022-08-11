process INFILE_HANDLING {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    input:
        tuple val(basename), path(input)

    output:
        path input, emit: input
        val basename, emit: base
        path ".command.out"
        path ".command.err"
        
    shell:
        '''

        source bash_functions.sh
        
        echo "INFO: R1 = !{input[0]}"
        echo "INFO: R2 = !{input[1]}"

        verify_file_minimum_size !{input[0]} 'fastq' '10M'
        verify_file_minimum_size !{input[1]} 'fastq' '10M'

        '''
}