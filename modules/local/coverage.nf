process COVERAGE {
    publishDir "${params.outpath}/qa",
        mode: "${params.publish_dir_mode}",
        pattern: "*"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "snads/bedtools@sha256:9b80fb5c5ef1b6f4a4a211d8739fa3fe107da34d1fb6609d6b70ddc7afdce12c"

    input:
        path single_bam
        path paired_bam

    output:
        path "*"

    shell:
    '''
    # Calculate coverage

    # Get basename of input file
    base=$(basename "!{single_bam}" | cut -d . -f 1 | sed 's/[-.,]//g')

    single_cov='0 bp TooFewToMap Singleton Reads (0.0x)\t'
    if [ -s !{single_bam} ]; then
        single_cov=$(bedtools genomecov -d -split -ibam !{single_bam} |\
        awk '{sum+=$3} END{print sum " bp Singleton Reads Mapped (" sum/NR "x)\t"}')
    fi

    cov_nfo=$(bedtools genomecov -d -split -ibam ${base}.paired.bam |\
    awk -v SEcov="${single_cov}" 'BEGIN{sum=0} {sum+=$3} END{
    print sum " bp Paired Reads Mapped (" sum/NR "x)\t" SEcov NR " bp Genome"}')

    #rm -f ${base}.{paired,single}.bam{,.bai}
    echo -e "${base}\t${cov_nfo}" >> \
    Summary.Illumina.CleanedReads-AlnStats.tab
    
    '''
}