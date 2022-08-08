process REMOVE_PHIX {

    publishDir "${params.outpath}/trim_reads",
        mode: "${params.publish_dir_mode}",
        pattern: "*.tsv"
    publishDir "${params.process_log_dir}",
        mode: "${params.publish_dir_mode}",
        pattern: ".command.*",
        saveAs: { filename -> "${task.process}${filename}"}
    
    container "snads/bbtools@sha256:9f2a9b08563839cec87d856f0fc7607c235f464296fd71e15906ea1d15254695"
    
    input:
        path PHIX
        path R1
        path R2

    output:
        path "*noPhiX-R1.fsq", emit: noPhiX_R1
        path "*noPhiX-R2.fsq", emit: noPhiX_R2
        path "*raw.tsv"
        path "*phix.tsv"
        path ".command.out"
        path ".command.err"

    shell:
    '''
    source bash_functions.sh

    # Get basename of input file
    base=$(basename "!{R1}" | cut -d _ -f 1 | sed 's/[-.,]//g')
    
    # Remove PhiX
    if ! verify_file_minimum_size !{PHIX} 'PhiX genome' '5k'; then
        echo "ERROR: PhiX file is too small" >&2
        exit 1
    fi

    echo "INFO: Starting bbduck"

    NSLOTS=$(cat /sys/devices/system/cpu/present | cut -d '-' -f2)
    echo "INFO: Number of threads found: ${NSLOTS}"

    bbduk.sh threads=4 k=31 hdist=1\
    ref="!{PHIX}" in="!{R1}" in2="!{R2}"\
    out=${base}_noPhiX-R1.fsq out2=${base}_noPhiX-R2.fsq\
    qin=auto qout=33 overwrite=t

    echo "INFO: bbduck finished"

    for suff in R1.fsq R2.fsq ; do
        verify_file_minimum_size "${base}_noPhiX-${suff}" 'PhiX cleaned read' '25M'
    done

    TOT_READS=$(grep '^Input: ' .command.err \
    | awk '{print $2}')
    TOT_BASES=$(grep '^Input: ' .command.err \
    | awk '{print $4}')

    if [[ -z "${TOT_READS}" || -z "${TOT_BASES}" ]]; then
        echo 'ERROR: unable to parse input counts from bbduk log' >&2
        exit 1
    fi

    PHIX_READS=$(grep '^Contaminants: ' .command.err \
    | awk '{print $2}' | sed 's/,//g')
    PHIX_BASES=$(grep '^Contaminants: ' .command.err \
    | awk '{print $5}' | sed 's/,//g')

    echo "INFO: $TOT_BASES bp and $TOT_READS reads provided as raw input" >&2
    echo "INFO: ${PHIX_BASES:-0} bp of PhiX were detected and removed in ${PHIX_READS:-0} reads" >&2

    echo -e "${base}\t${TOT_BASES} bp Raw\t${TOT_READS} reads Raw" \
    > ${base}_raw.tsv
    echo -e "${base}\t${PHIX_BASES:-0} bp PhiX\t${PHIX_READS:-0} reads PhiX" \
    > ${base}_phix.tsv
    '''

}