process MLST {
    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/mlst@sha256:27f290753760c44204d6e04b6ead7935d03b48d5f0a5ccce068def9ce33babe6"

    input:
        path base_fna

    output:
        path "Summary.MLST.tab"

    shell:
    '''
    # MLST for each assembly

    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
    echo "INFO: Number of threads found: ${NSLOTS}"

    if [ -s !{base_fna} ]; then
        mlst --threads $NSLOTS "!{base_fna}" \
        >> Summary.MLST.tab
    fi


    '''
}