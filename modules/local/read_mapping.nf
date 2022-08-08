process READ_MAPPING {
    publishDir "${params.outpath}",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/bwa@sha256:60fc4074ac71ee4ab4a14320ea3bcca50d315c6ef550a2be82791683a15add44"

    input:
        path single_gz
        path base_fna

    output:
        path "*.single.bam", emit: single_bam

    shell:
    '''
    source bash_functions.sh

    # Get basename of input file
    base=$(basename "!{single_gz}" | cut -d _ -f 1 | sed 's/[-.,]//g')

    # Single read mapping if available
    if [[ $(find -L !{single_gz} -type f -size +5k 2> /dev/null) ]]; then
        bwa index !{base_fna}

        NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
        echo "INFO: Number of threads found: ${NSLOTS}"

        #bwa mem -t $NSLOTS -x intractg -v 2 !{base_fna}\
        !{single_gz} 
        #|\samtools sort -@ $NSLOTS --reference !{base_fna} -l 9\
        #-o ${base}.single.bam

        #rm -f !{base_fna}.{ann,amb,bwt,pac,sa}
        #verify_file_minimum_size "${base}.single.bam" 'binary sequence alignment map' '1k'
        #samtools index ${base}.single.bam

    fi

    '''
}