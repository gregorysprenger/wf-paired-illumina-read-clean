process KRAKEN_ONE {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    label "process_medium"

    container "staphb/kraken@sha256:d372099288c3a7c0cc90ea7e516c643e7096c90a551b45d531bd26b4e7f46255"

    input:
        path R1_paired_gz
        path R2_paired_gz
        path single_gz
        val base

    output:
        path "*taxonomy-reads.tab"
        path "*kraken.tab.gz"
        path ".command.out"
        path ".command.err"

    shell:
    '''

    source summarize_kraken.sh

    # Investigate taxonomic identity of cleaned reads
    if [ ! -s !{base}_taxonomy-reads.tab ]; then

        NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
        echo "INFO: Number of threads found: !{task.cpus}"

        echo "INFO: Starting Kraken1"
        kraken --db /kraken-database/minikraken_20171013_4GB --threads !{task.cpus} --fastq-input --gzip-compressed \
        !{R1_paired_gz} !{R2_paired_gz} !{single_gz} > !{base}_kraken.output

        echo "INFO: Run kraken-report"
        kraken-report --db /kraken-database/minikraken_20171013_4GB !{base}_kraken.output > kraken.tab 2>&1 | tr '^M' '\n' 1>&2
        echo "INFO: Kraken1 finished"

        echo "INFO: Summarize Kraken1"
        summarize_kraken 'kraken.tab' > !{base}_taxonomy-reads.tab

        echo "INFO: gzip kraken1"
        mv kraken.tab !{base}_kraken.tab
        gzip !{base}_kraken.tab
    fi

    '''
}

process KRAKEN_TWO {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    label "process_medium"

    container "staphb/kraken2@sha256:5b107d0141d6042a6b0ac6a5852990dc541fbff556a85eb0c321a7771200ba56"

    input:
        path R1_paired_gz
        path R2_paired_gz
        path single_gz
        val base

    output:
        path "*taxonomy2-reads.tab"
        path "*kraken2.tab.gz"

    shell:
    '''

    source summarize_kraken.sh

    if [ ! -s !{base}_taxonomy2-reads.tab ]; then

        NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
        echo "INFO: Number of threads found: !{task.cpus}"

        echo "INFO: Starting Kraken2"
        kraken2 --db /kraken2-db --threads !{task.cpus} --gzip-compressed --output /dev/null \
        --use-names --report kraken2.tab \
        !{R1_paired_gz} !{R2_paired_gz} !{single_gz}
        echo "INFO: Kraken2 finished"

        echo "INFO: Summarize kraken2"
        summarize_kraken 'kraken2.tab' > !{base}_taxonomy2-reads.tab

        echo "INFO: gzip kraken2"
        mv kraken2.tab !{base}_kraken2.tab
        gzip !{base}_kraken2.tab

    fi
    
    '''
}