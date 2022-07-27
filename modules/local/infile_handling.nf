process INFILE_HANDLING {
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}" }

    input:
        path input

    output:
        path "*R1*", emit: R1
        path "*R2*", emit: R2
        path ".command.out"
        path ".command.err"
        
    shell:
        '''
        source bash_functions.sh
        
        shopt -s nullglob
        R1=( "!{input}"/*_R1*{fq,fastq}.gz )
        R2=( "!{input}"/*_R2*{fq,fastq}.gz )
        echo "INFO: ${R1} found"
        echo "INFO: ${R2} found"
        shopt -u nullglob

        verify_file_minimum_size ${R1} 'fastq' '10M'
        verify_file_minimum_size ${R2} 'fastq' '10M'
    
        cp ${R1} .
        cp ${R2} .
        
        '''

}