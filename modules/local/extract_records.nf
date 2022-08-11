process EXTRACT_RECORDS {

    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/biopython@sha256:bb041f55fd45d0fb577656e2d1f1a9f477d3ba80878b3b42218adff3322ae06e"

    input:
        path extract_record
        path annotation
        val base

    output:
        path "16S.*.fa", emit: extracted_rna
        path ".command.out"
        path ".command.err"

    shell:
    '''

    # 16S extraction
    if [[ -s "!{annotation}" ]]; then
        python3 !{extract_record} -i "!{annotation}" \
        -u product -o "16S.!{base}.fa" -q '16S ribosomal RNA' \
        --search-type any_q_is_rec -f rRNA
    fi

    '''
}