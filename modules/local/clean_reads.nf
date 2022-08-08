process CLEAN_READS {
    publishDir "${params.outpath}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*.txt"
    publishDir "${params.outpath}/asm",
        mode: "${params.publish_dir_mode}",
        pattern: "*.fna"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}

    container "gregorysprenger/bwa-samtools-pilon@sha256:209ac13b381188b4a72fe746d3ff93d1765044cbf73c3957e4e2f843886ca57f"
    
    input:
        path uncorrected_contigs
        path R1_paired_gz
        path R2_paired_gz
        path single_gz
        path outpath

    output:
        path "*.InDels-corrected.cnt.txt"
        path "*.SNPs-corrected.cnt.txt"
        path "*.fna", emit: base_fna
        path "*.paired.bam", emit: paired_bam
        path "*.single.bam", emit: single_bam

    shell:
    '''
    source bash_functions.sh

    # Get basename of input file
    base=$(basename "!{R1_paired_gz}" | cut -d _ -f 1 | sed 's/[-.,]//g')

    # Correct cleaned SPAdes contigs with cleaned PE reads
    #verify_file_minimum_size "!{uncorrected_contigs}" 'filtered SPAdes assembly' '1k' #1M

    echo -n '' > ${base}.InDels-corrected.cnt.txt
    echo -n '' > ${base}.SNPs-corrected.cnt.txt

    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
    echo "INFO: Number of threads found: ${NSLOTS}"

    for _ in {1..3}; do
        bwa index !{uncorrected_contigs}
        bwa mem -t ${NSLOTS} -x intractg -v 2 !{uncorrected_contigs}\
        !{R1_paired_gz} !{R2_paired_gz} |\
        samtools sort -@ ${NSLOTS} --reference !{uncorrected_contigs} -l 9\
        -o ${base}.paired.bam

        #rm -f !{uncorrected_contigs}.{ann,amb,bwt,pac,sa}

        #verify_file_minimum_size "${base}.paired.bam" 'binary sequence alignment map' '1M' #25M

        samtools index ${base}.paired.bam

        pilon --genome !{uncorrected_contigs} --frags ${base}.paired.bam\
        --output "${base}" --changes \
        --fix snps,indels --mindepth 0.50 --threads ${NSLOTS} >&2

        #verify_file_minimum_size "!{uncorrected_contigs}" 'polished assembly' '1K' #1M

        echo $(grep -c '-' ${base}.changes >> ${base}.InDels-corrected.cnt.txt)
        echo $(grep -vc '-' ${base}.changes >> ${base}.SNPs-corrected.cnt.txt)

        rm -f ${base}.{changes,uncorrected.fna}
        rm -f "${base}"Pilon.bed
        mv -f ${base}.fasta ${base}.uncorrected.fna

        sed -i 's/_pilon//1' ${base}.uncorrected.fna

    done

    mv -f ${base}.uncorrected.fna ${base}.fna

    #verify_file_minimum_size "${base}.fna" 'corrected SPAdes assembly' '1k' #1M

    # Single read mapping if available
    if [[ !{single_gz} ]]; then
        bwa index ${base}.fna

        bwa mem -t ${NSLOTS} -x intractg -v 2 ${base}.fna\
        !{single_gz} |\
        samtools sort -@ ${NSLOTS} --reference ${base}.fna -l 9\
        -o ${base}.single.bam

        #verify_file_minimum_size "${base}.single.bam" 'binary sequence alignment map' '1k'
        samtools index ${base}.single.bam

    fi

    '''
}