process ANNOTATE {
    
    publishDir "${params.outpath}/annot",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    label "process_medium"

    container "snads/prokka@sha256:ef7ee0835819dbb35cf69d1a2c41c5060691e71f9138288dd79d4922fa6d0050"

    input:
        path base_fna
        
    output:
        path "*.gbk", emit: annotation
        path ".command.out"
        path ".command.err"

    shell:
    '''

    source bash_functions.sh
    
    # Annotate cleaned and corrected assembly

    # Get basename of input file
    base=$(basename "!{base_fna}" | cut -d . -f 1 | sed 's/[-.,]//g')

    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
    echo "INFO: Number of threads found: !{task.cpus}"

    prokka --outdir prokka --prefix "${base}"\
    --force --addgenes --locustag "${base}" --mincontiglen 1\
    --evalue 1e-08 --cpus !{task.cpus} !{base_fna}

    for ext in gb gbf gbff gbk ; do
    if [ -s "prokka/${base}.${ext}" ]; then
        #verify_file_minimum_size "${base}.${ext}" 'annotated assembly' '3M'
        mv -f prokka/${base}.${ext} ${base}.gbk
        rm -rf ${base}
        break
    fi
    done

    '''
}