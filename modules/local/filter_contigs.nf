process FILTER_CONTIGS {
    
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/biopython@sha256:bb041f55fd45d0fb577656e2d1f1a9f477d3ba80878b3b42218adff3322ae06e"

    input:
        path filter_contigs
        path contigs
        path R1_paired_gz
        path outpath

    output:
        path "*.uncorrected.fna", emit: uncorrected_contigs
        path ".command.out"
        path ".command.err"

    shell:
    '''

    # Get basename of input file
    base=$(basename "!{R1_paired_gz}" | cut -d _ -f 1 | sed 's/[-.,]//g')

    python3 !{filter_contigs} \
    -i !{contigs}\
    -b "${base}" -l 1000 -o ${base}.uncorrected.fna

    '''
}